module npy_class
  use precision, only: WP
  use string, only: str_medium, str_long
  use pgrid_class, only: pgrid
  use partmesh_class, only: partmesh
  implicit none
  private

  ! limitations:
  !  -  doesn't work on big endian systems (shouldn't be too hard to fix, but
  !     I am lazy)
  !  -  this is *not* a performant way to write large amounts of data; it is
  !     an easy way to do analysis on small amounts of data through python
  !     without all the trouble of using an ensight reader
  !  -  an ideal future direction would be to implement an HDF5 based output;
  !     support across languages is very good for HDF5

  public :: npy

  ! file related constants
  character(len=*), parameter :: SLASH = "/"
  character(len=*), parameter :: DEFAULT_PATH = "./npy/"

  ! constants needed for the npy format
  integer, parameter :: FIELDLEN = str_medium
  character(len=*), parameter :: NPYDIR = "./npy/"
  character(len=*), parameter :: NPYMETA = "meta.json"
  character(len=*), parameter :: NPYMAGIC = transfer(-109_1, 'x') // 'NUMPY'
  character(len=*), parameter :: NPYVERSION = transfer(1_1, 'x') // transfer(0_1, 'x')
  character(len=*), parameter :: NPYDESCR4 = "'descr': '=f'"
  character(len=*), parameter :: NPYDESCR8 = "'descr': '=d'"
  character(len=*), parameter :: NPYORDER = "'fortran_order': True"
  character(len=*), parameter :: NPYSHAPE1 = "'shape': ("
  character(len=*), parameter :: NPYSHAPE2 = ")"

  ! default units for writing files
  integer, parameter :: NPY_DEFAULT_UNIT = 93
  integer, parameter :: JSON_DEFAULT_UNIT = 94

  ! since it's fortran, I guess we have to write another linked list real
  ! quick...
  type :: llitem
    type(llitem), pointer :: next => null()
    character(len=FIELDLEN) :: name
    real(WP), dimension(:),       pointer :: rdata1d   => null()
    real(WP), dimension(:,:),     pointer :: rdata2d   => null()
    real(WP), dimension(:,:,:),   pointer :: rdata3d   => null()
    real(WP), dimension(:,:,:),   pointer :: rdata4dc1 => null()
    real(WP), dimension(:,:,:),   pointer :: rdata4dc2 => null()
    real(WP), dimension(:,:,:),   pointer :: rdata4dc3 => null()
  end type llitem
  type :: pmitem
    type(pmitem), pointer :: next => null()
    character(len=FIELDLEN) :: name
    type(partmesh), pointer :: pm => null()
  end type pmitem

  ! the npy writer class, written to roughly mirror the ensight writer class
  type :: npy

    class(pgrid), pointer :: pg => null()

    character(len=str_medium) :: meta_fn
    character(len=str_medium) :: dir
    logical :: mkdirs

    integer, dimension(:,:), allocatable :: blkbds

    integer :: Nt
    real(WP), dimension(:), allocatable :: ts

    type(llitem), pointer :: fieldscalars => null()   ! real 3d
    type(llitem), pointer :: fieldvectors => null()   ! real 4d

    type(pmitem), pointer :: pms => null()            ! particle data

  contains

    procedure :: write_meta
    procedure :: write_data
    procedure :: add_scalar
    procedure :: add_vector
    procedure :: add_particle
    procedure, private :: get_filename

  end type npy

  interface npy; procedure npy_from_args; end interface

contains

  function npy_from_args(pg, folder, path, metaname, skipmkdirs) result(this)
    use mpi_f08, only: MPI_INTEGER, mpi_gather, mpi_barrier
    use messager, only: die
    implicit none
    type(npy) :: this
    class(pgrid), intent(in), target :: pg
    character(len=*), intent(in) :: folder
    character(len=str_medium), intent(in), optional :: path, metaname
    logical, intent(in), optional :: skipmkdirs
    character(len=str_medium) :: path_actual, meta_actual
    integer :: ierr

    ! store pg pointer
    this%pg => pg

    ! set directories
    this%mkdirs = .true.
    if (present(skipmkdirs)) this%mkdirs = .not. skipmkdirs
    path_actual = DEFAULT_PATH
    if (present(path)) path_actual = trim(path)
    this%dir = trim(path_actual) // SLASH // trim(folder) // SLASH
    meta_actual = NPYMETA
    if (present(metaname)) meta_actual = trim(metaname)
    this%meta_fn = trim(this%dir) // trim(meta_actual)

    ! make directories if they don't exist
    if (this%pg%amroot .and. this%mkdirs) call execute_command_line('mkdir -p ' // trim(this%dir))
    call mpi_barrier(this%pg%comm, ierr)

    ! allocate initial array of times
    this%Nt = 0
    allocate(this%ts(20))

    ! collect block sizes
    allocate(this%blkbds(4, this%pg%npx * this%pg%npy))
    call mpi_gather((/ this%pg%imin_, this%pg%imax_, this%pg%jmin_,           &
      this%pg%jmax_ /), 4, MPI_INTEGER, this%blkbds, 4, MPI_INTEGER, 0,       &
      this%pg%xycomm, ierr)

    ! make sure the main root is an xy root
    if (this%pg%rank .eq. 0 .and. this%pg%xyrank .ne. 0)                      &
      call die('[npy] comm root must be an xycomm root')

  end function npy_from_args

  subroutine add_scalar(this, name, d)
    implicit none
    class(npy), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), target, intent(in) :: d
    type(llitem), pointer :: newitem

    allocate(newitem)
    newitem%name = trim(name)
    newitem%rdata3d => d(this%pg%imin_:this%pg%imax_,                         &
      this%pg%jmin_:this%pg%jmax_, this%pg%kmin_:this%pg%kmax_)
    newitem%next => this%fieldscalars
    this%fieldscalars => newitem

  end subroutine add_scalar

  subroutine add_vector(this, name, d1, d2, d3)
    use messager, only: die
    implicit none
    class(npy), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), target, intent(in) :: d1, d2, d3
    type(llitem), pointer :: newitem

    if (any(shape(d1) .ne. shape(d2)) .or. any(shape(d1) .ne. shape(d3)))     &
      call die('[npy] vector dimensions do not match')

    allocate(newitem)
    newitem%name = name
    newitem%rdata4dc1 => d1(this%pg%imin_:this%pg%imax_,                      &
      this%pg%jmin_:this%pg%jmax_, this%pg%kmin_:this%pg%kmax_)
    newitem%rdata4dc2 => d2(this%pg%imin_:this%pg%imax_,                      &
      this%pg%jmin_:this%pg%jmax_, this%pg%kmin_:this%pg%kmax_)
    newitem%rdata4dc3 => d3(this%pg%imin_:this%pg%imax_,                      &
      this%pg%jmin_:this%pg%jmax_, this%pg%kmin_:this%pg%kmax_)
    newitem%next => this%fieldvectors
    this%fieldvectors => newitem

  end subroutine add_vector

  subroutine add_particle(this, name, pm)
    implicit none
    class(npy), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(partmesh), intent(in), target :: pm
    type(pmitem), pointer :: newitem

    allocate(newitem)
    newitem%name = trim(name)
    newitem%pm => pm
    newitem%next => this%pms
    this%pms => newitem

  end subroutine add_particle

  subroutine write_data(this, t)
    use mpi_f08, only: MPI_REQUEST, MPI_STATUS, mpi_irecv, mpi_isend,         &
      mpi_wait, MPI_MODE_WRONLY, MPI_MODE_CREATE,                          &
      mpi_file_open, mpi_file_close, mpi_file_write_at, MPI_FILE,             &
      MPI_STATUS_SIZE, MPI_OFFSET_KIND, MPI_INTEGER, mpi_allgather
    use parallel, only: MPI_REAL_WP, info_mpiio
    use messager, only: die
    implicit none
    class(npy), intent(inout) :: this
    real(WP), intent(in) :: t
    type(MPI_REQUEST) :: sreq
    type(MPI_STATUS) :: sstat, rstat, wstat
    type(MPI_FILE) :: fh
    integer(MPI_OFFSET_KIND) :: offset
    character(len=str_medium) :: fn
    type(llitem), pointer :: item
    type(pmitem), pointer :: pitem
    integer :: i, li, hi, lj, hj, ntot, bufsize, headsize, ierr
    integer, dimension(this%pg%nproc) :: partcounts
    integer, dimension(:), allocatable :: dims
    real(WP), dimension(:),   allocatable :: wbuf1d
    real(WP), dimension(:,:), allocatable :: wbuf2d
    real(WP), dimension(:,:,:),   allocatable :: sbuf3d, rbuf3d, slice3d
    real(WP), dimension(:,:,:,:), allocatable :: sbuf4d, rbuf4d, slice4d

    ! add the time
    this%Nt = this%Nt + 1
    if (size(this%ts) .lt. this%Nt) then
      realloctime : block
        real(WP), dimension(:), allocatable :: newts
        allocate(newts(1:size(this%ts,1)))
        newts(:) = this%ts(:)
        deallocate(this%ts); allocate(this%ts(2*this%Nt));
        this%ts(1:this%Nt-1) = this%ts(1:this%Nt-1)
        deallocate(newts)
      end block realloctime
    end if
    this%ts(this%Nt) = t

    ! write the field scalars
    allocate(dims(3))
    dims(:) = (/ this%pg%nx, this%pg%ny, this%pg%nz /)
    allocate(sbuf3d(this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_,  &
      this%pg%kmin_:this%pg%kmax_))
    if (this%pg%xyrank .eq. 0)                                                &
      allocate(slice3d(this%pg%imin:this%pg%imax,this%pg%jmin:this%pg%jmax,   &
        this%pg%kmin_:this%pg%kmax_))
    item => this%fieldscalars
    do while (associated(item))
      call this%get_filename(item%name, this%Nt, fn)
      if (this%pg%rank .eq. 0) call write_array_header(fn, dims, headsize)
      call mpi_bcast(headsize, 1, MPI_INTEGER, 0, this%pg%comm, ierr)
      bufsize = this%pg%nx_ * this%pg%ny_ * this%pg%nz_
      sbuf3d(:,:,:) = item%rdata3d(:,:,:)
      call mpi_isend(sbuf3d, bufsize, MPI_REAL_WP, 0, this%pg%xyrank,   &
        this%pg%xycomm, sreq, ierr)
      if (this%pg%xyrank .eq. 0) then
        do i = 1, this%pg%npx * this%pg%npy
          li = this%blkbds(1,i); hi = this%blkbds(2,i);
          lj = this%blkbds(3,i); hj = this%blkbds(4,i);
          allocate(rbuf3d(li:hi,lj:hj,this%pg%kmin_:this%pg%kmax_))
          call mpi_recv(                                                      &
            rbuf3d, (hi-li+1)*(hj-lj+1)*this%pg%nz_,                          &
            MPI_REAL_WP, i-1, i-1, this%pg%xycomm, rstat, ierr)
          slice3d(li:hi,lj:hj,:) = rbuf3d(:,:,:)
          deallocate(rbuf3d)
        end do
        bufsize = this%pg%nx * this%pg%ny * this%pg%nz_
        !TODO this only works if each processor owns some of the data; maybe this is a problem?
        offset = int(headsize + WP * (this%pg%kproc - 1) * bufsize, MPI_OFFSET_KIND)
        call mpi_barrier(this%pg%zcomm, ierr)
        call mpi_file_open(this%pg%zcomm, trim(fn), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE),  &
          info_mpiio, fh, ierr)
        if (ierr .ne. 0) call die('[npy write data] could not open file ' // trim(fn))
        call mpi_file_write_at(fh, offset, slice3d, bufsize, MPI_REAL_WP, wstat, ierr)
        call mpi_file_close(fh, ierr)
      end if
      item => item%next
    end do
    call mpi_wait(sreq, sstat, ierr)
    deallocate(dims, sbuf3d)
    if (this%pg%xyrank .eq. 0) deallocate(slice3d)

    ! write the field vectors
    allocate(dims(4))
    dims(:) = (/ this%pg%nx, this%pg%ny, this%pg%nz, 3 /)
    allocate(sbuf4d(3,this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_,&
      this%pg%kmin_:this%pg%kmax_))
    if (this%pg%xyrank .eq. 0)                                                &
      allocate(slice4d(3,this%pg%imin:this%pg%imax,this%pg%jmin:this%pg%jmax, &
        this%pg%kmin_:this%pg%kmax_))
    item => this%fieldvectors
    do while (associated(item))
      call this%get_filename(item%name, this%Nt, fn)
      if (this%pg%rank .eq. 0) call write_array_header(fn, dims, headsize)
      call mpi_bcast(headsize, 1, MPI_INTEGER, 0, this%pg%comm, ierr)
      bufsize = this%pg%nx_ * this%pg%ny_ * this%pg%nz_ * 3
      sbuf4d(1,:,:,:) = item%rdata4dc1(:,:,:)
      sbuf4d(2,:,:,:) = item%rdata4dc2(:,:,:)
      sbuf4d(3,:,:,:) = item%rdata4dc3(:,:,:)
      call mpi_isend(sbuf4d, bufsize, MPI_REAL_WP, 0, this%pg%xyrank,   &
        this%pg%xycomm, sreq, ierr)
      if (this%pg%xyrank .eq. 0) then
        do i = 1, this%pg%npx * this%pg%npy
          li = this%blkbds(1,i); hi = this%blkbds(2,i);
          lj = this%blkbds(3,i); hj = this%blkbds(4,i);
          allocate(rbuf4d(3,li:hi,lj:hj,this%pg%kmin_:this%pg%kmax_))
          call mpi_recv(                                                      &
            rbuf4d, 3*(hi-li+1)*(hj-lj+1)*this%pg%nz_,                        &
            MPI_REAL_WP, i-1, i-1, this%pg%xycomm, rstat, ierr)
          slice4d(:,li:hi,lj:hj,:) = rbuf4d(:,:,:,:)
          deallocate(rbuf4d)
        end do
        bufsize = this%pg%nx * this%pg%ny * this%pg%nz * 3
        !TODO this only works if each processor owns some of the data; maybe this is a problem?
        offset = int(headsize + WP * (this%pg%kproc - 1) * bufsize, MPI_OFFSET_KIND)
        call mpi_barrier(this%pg%zcomm, ierr)
        call mpi_file_open(this%pg%zcomm, trim(fn), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE),  &
          info_mpiio, fh, ierr)
        if (ierr .ne. 0) call die('[npy write data] could not open file ' // trim(fn))
        call mpi_file_write_at(fh, offset, slice4d, bufsize, MPI_REAL_WP, wstat, ierr)
        call mpi_file_close(fh, ierr)
      end if
      item => item%next
    end do
    call mpi_wait(sreq, sstat, ierr)
    deallocate(dims, sbuf4d)
    if (this%pg%xyrank .eq. 0) deallocate(slice4d)

    pitem => this%pms
    do while (associated(pitem))

      ! share particle counts
      call mpi_allgather(pitem%pm%n, 1, MPI_INTEGER, partcounts, 1, MPI_INTEGER, this%pg%comm, ierr)
      ntot = sum(partcounts)

      ! write particle positions
      call this%get_filename('pos', this%Nt, fn)
      if (this%pg%rank .eq. 0) call write_array_header(fn, (/ 3, ntot /), headsize)
      call mpi_bcast(headsize, 1, MPI_INTEGER, 0, this%pg%xycomm, ierr)
      bufsize = 3 * pitem%pm%n
      offset = int(headsize + WP * 3 * sum(partcounts(1:this%pg%rank)), MPI_OFFSET_KIND)
      call mpi_file_open(this%pg%comm, trim(fn), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE),  &
        info_mpiio, fh, ierr)
      if (ierr .ne. 0) call die('[npy write data] could not open file ' // trim(fn))
      call mpi_file_write_at(fh, offset, pitem%pm%pos, bufsize, MPI_REAL_WP, wstat, ierr)
      call mpi_file_close(fh, ierr)

      ! write the particle scalars
      allocate(wbuf1d(pitem%pm%n))
      do i = 1, pitem%pm%nvar
        call this%get_filename(trim(pitem%pm%varname(i)), this%Nt, fn)
        if (this%pg%rank .eq. 0) call write_array_header(fn, (/ ntot /), headsize)
        call mpi_bcast(headsize, 1, MPI_INTEGER, 0, this%pg%xycomm, ierr)
        bufsize = pitem%pm%n
        offset = int(headsize + WP * sum(partcounts(1:this%pg%rank)), MPI_OFFSET_KIND)
        call mpi_file_open(this%pg%comm, trim(fn), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE),  &
          info_mpiio, fh, ierr)
        if (ierr .ne. 0) call die('[npy write data] could not open file ' // trim(fn))
        wbuf1d(:) = pitem%pm%var(i,:)
        call mpi_file_write_at(fh, offset, wbuf1d, bufsize, MPI_REAL_WP, wstat, ierr)
        call mpi_file_close(fh, ierr)
      end do
      deallocate(wbuf1d)

      ! write particle vectors
      allocate(wbuf2d(3,pitem%pm%n))
      do i = 1, pitem%pm%nvec
        call this%get_filename(trim(pitem%pm%vecname(i)), this%Nt, fn)
        if (this%pg%rank .eq. 0) call write_array_header(fn, (/ 3, ntot /), headsize)
        call mpi_bcast(headsize, 1, MPI_INTEGER, 0, this%pg%xycomm, ierr)
        bufsize = 3 * pitem%pm%n
        offset = int(headsize + WP * 3 * sum(partcounts(1:this%pg%rank)), MPI_OFFSET_KIND)
        call mpi_file_open(this%pg%comm, trim(fn), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE),  &
          info_mpiio, fh, ierr)
        if (ierr .ne. 0) call die('[npy write data] could not open file ' // trim(fn))
        wbuf2d(:,:) = pitem%pm%vec(:,i,:)
        call mpi_file_write_at(fh, offset, wbuf2d, bufsize, MPI_REAL_WP, wstat, ierr)
        call mpi_file_close(fh, ierr)
      end do
      deallocate(wbuf2d)

      pitem => pitem%next

    end do

    ! rewrite metadata
    call this%write_meta()

  end subroutine write_data

  subroutine list_size(list, N)
    implicit none
    type(llitem), pointer, intent(in) :: list
    integer, intent(out) :: N
    type(llitem), pointer :: item

    N = 0; item => list;
    do while (associated(item)); N = N + 1; item => item%next; end do;

  end subroutine list_size

  subroutine write_meta(this)
    implicit none
    class(npy), intent(inout) :: this
    character(len=FIELDLEN), dimension(:), allocatable ::                     &
      fieldscalars, fieldvectors
    integer :: i, Nfs, Nfv
    real(WP), dimension(:), allocatable :: xs, ys, zs, ts
    type(llitem), pointer :: item

    call list_size(this%fieldscalars, Nfs)
    call list_size(this%fieldvectors, Nfv)

    allocate(fieldscalars(Nfs), fieldvectors(Nfv))

    item => this%fieldscalars;
    do i = 1, Nfs
      fieldscalars(i) = item%name
      item => item%next
    end do

    item => this%fieldvectors;
    do i = 1, Nfv
      fieldvectors(i) = item%name
      item => item%next
    end do

    !TODO shift grid points for different positions within cell here
    allocate(xs(this%pg%nx), ys(this%pg%ny), zs(this%pg%nz), ts(this%Nt))
    xs(:) = this%pg%xm(this%pg%imin:this%pg%imax)
    ys(:) = this%pg%ym(this%pg%jmin:this%pg%jmax)
    zs(:) = this%pg%zm(this%pg%kmin:this%pg%kmax)
    ts(:) = this%ts(1:this%Nt)

    if (this%pg%rank .eq. 0) call write_meta_sub(this%meta_fn, NPYMETA, &
      fieldscalars, fieldvectors, this%pms,    &
      xs, ys, zs, ts)

    deallocate(fieldscalars, fieldvectors)
    deallocate(xs, ys, zs, ts)

  end subroutine write_meta

  subroutine get_filename(this, field, n, fn)
    implicit none
    class(npy), intent(in) :: this
    character(len=*), intent(in) :: field
    integer, intent(in) :: n
    character(len=str_medium), intent(out) :: fn

    write(fn, '(a,a,a,a,i0.8,a)') trim(this%dir), SLASH, trim(field), "_", n, ".npy"

  end subroutine get_filename

  !TODO combine this with write_meta; there's not a good reason to split this
  subroutine write_meta_sub(fn, name, fieldscalars, fieldvectors,             &
      plist, xs, ys, zs, ts, unit)
    implicit none
    character(len=*), intent(in) :: fn, name
    character(len=FIELDLEN), dimension(:), intent(in) ::                      &
      fieldscalars, fieldvectors
    type(pmitem), pointer, intent(in) :: plist
    real(WP), dimension(1:), intent(in) :: xs, ys, zs, ts
    integer, intent(in), optional :: unit
    integer :: unit_actual
    type(pmitem), pointer :: pitem

    unit_actual = NPY_DEFAULT_UNIT
    if (present(unit)) unit_actual = unit

    open(unit=unit_actual, file=fn, access='stream', form='formatted',        &
      action='write', status='replace')

    write(unit_actual, '(a)') "{"
    write(unit_actual, '(3a)') """name"" : """, name, ""","
    call write_meta_string_list(unit_actual, "fieldscalars", fieldscalars)
    call write_meta_string_list(unit_actual, "fieldvectors", fieldvectors)
    write(unit_actual, '(a)') """partgroups"" : {"
    pitem => plist
    do while (associated(pitem))
      call write_meta_part_desc(unit_actual, pitem%name,  pitem%pm,           &
        associated(pitem%next))
      pitem => plist%next
    end do
    write(unit_actual, '(a)') "},"
    !TODO add field locations off of cell center
    write(unit_actual, '(a)') """fieldlocations"" : ""c"","
    call write_meta_number_list(unit_actual, "x", xs)
    call write_meta_number_list(unit_actual, "y", ys)
    call write_meta_number_list(unit_actual, "z", zs)
    call write_meta_number_list(unit_actual, "t", ts)
    write(unit_actual, '(a)') """fieldfilenames"" : ""{field:}_{n:08d}.npy"","
    write(unit_actual, '(a)') """particlefilenames"" : ""{partgroup:}_{prop:}_{n:08d}.npy"""
    write(unit_actual, '(a)') "}"

    close(unit=unit_actual)

  end subroutine write_meta_sub

  subroutine write_meta_part_desc(unit, label, pm, inclcomma)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: label
    type(partmesh), intent(in) :: pm
    logical, intent(in) :: inclcomma
    integer :: i

    write(unit, '(3a)') "    """, trim(label), """ : {"
    write(unit, '(a)', advance='no') "        ""scalars"" : ["
    do i = 1, pm%nvar
      write(unit, '(3a)', advance='no') """", trim(pm%varname(i)), """"
      if (i .ne. pm%nvar) write(unit, '(a)', advance='no') ", "
    end do
    write(unit, '(a)') "],"
    write(unit, '(a)', advance='no') "        ""vectors"" : ["
    do i = 1, pm%nvec
      write(unit, '(3a)', advance='no') """", trim(pm%vecname(i)), """"
      if (i .ne. pm%nvar) write(unit, '(a)', advance='no') ", "
    end do
    write(unit, '(a)') "]"
    if (inclcomma) then
      write(unit, '(a)') "    },"
    else
      write(unit, '(a)') "    }"
    end if

  end subroutine write_meta_part_desc

  subroutine write_meta_number_list(unit, label, list)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: label
    real(WP), dimension(:), intent(in) :: list
    integer :: i

    write(unit, '(3a)', advance='no') """", trim(label), """ : ["
    do i = 1, size(list)
      write(unit, '(f8.3)', advance='no') list(i)
      if (i .ne. size(list)) write(unit, '(a)', advance='no') ", "
    end do
    write(unit, '(a)') "],"

  end subroutine write_meta_number_list

  subroutine write_meta_string_list(unit, label, list)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: label
    character(len=FIELDLEN), dimension(:), intent(in) :: list
    integer :: i

    write(unit, '(3a)', advance='no') """", trim(label), """ : ["
    do i = 1, size(list, 1)
      write(unit, '(3a)', advance='no') """", trim(list(i)), """"
      if (i .ne. size(list, 1)) write(unit, '(a)', advance='no') ", "
    end do
    write(unit, '(a)') "],"

  end subroutine write_meta_string_list

  subroutine write_array_header(fn, dims, headersize, unit)
    use messager, only: die
    implicit none
    character(len=*), intent(in) :: fn
    integer, dimension(:), intent(in) :: dims
    integer, intent(out) :: headersize
    integer, intent(in), optional :: unit
    character(len=str_long) :: header, dimstring
    character(len=len(NPYDESCR8)) :: descr
    integer :: unit_actual, i

    unit_actual = NPY_DEFAULT_UNIT
    if (present(unit)) unit_actual = unit

    ! assemble the header
    if (WP .eq. 8) then
      descr = NPYDESCR8
    elseif (WP .eq. 4) then
      descr = NPYDESCR4
    else
      call die('unsupported precision')
    end if
    if (size(dims) .eq. 1) then
      write(dimstring, '(i16)') dims(1)
    elseif (size(dims) .eq. 2) then
      write(dimstring, '(i16,a,i16)') dims(1), ",", dims(2)
    elseif (size(dims) .eq. 3) then
      write(dimstring, '(i16,a,i16,a,i16)') dims(1), ",", dims(2), ",", dims(3)
    elseif (size(dims) .eq. 4) then
      write(dimstring, '(i16,a,i16,a,i16,a,i16)') dims(1), ",", dims(2), ",", &
        dims(3), ",", dims(4)
    else
      call die('unsupported number of dimensions')
    end if
    write(header, '(9a)') "{", descr, ",", NPYORDER, ",", NPYSHAPE1,          &
      trim(dimstring), NPYSHAPE2, "}"

    ! compute length so that the header is a multiple of 64 bytes
    i = len_trim(header)
    headersize = 64 * ((i + 10) / 64 + 1) - 10

    open(unit=unit_actual, file=fn, access='stream', form='unformatted',      &
      action='write', status='replace')
    write(unit_actual) NPYMAGIC
    write(unit_actual) NPYVERSION
    write(unit_actual) int(headersize, 2)
    write(unit_actual) header(1:headersize)
    close(unit=unit_actual)

    headersize = headersize + 10

  end subroutine write_array_header

end module npy_class

