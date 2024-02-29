!> Adds particle coupling to an existing coupler
module lptcoupler_class
  use precision,      only: WP
  use mpi_f08,        only: MPI_COMM, MPI_GROUP, MPI_UNDEFINED,               &
    mpi_group_translate_ranks, MPI_INTEGER, MPI_LOGICAL
  use parallel,       only: MPI_REAL_WP
  use string,         only: str_medium
  use sgrid_class,    only: sgrid
  use pgrid_class,    only: pgrid
  use config_class,   only: config
  use lpt_class,      only: lpt, part
  implicit none
  private

  ! Expose type/constructor/methods
  public :: lptcoupler

  ! buffer sizing behavior
  real(WP), parameter :: MEM_ADJ_UP = 1.3_WP          !< Particle array size increase factor
  !TODO add buffer shrinking
  real(WP), parameter :: MEM_ADJ_DN = 0.7_WP          !< Particle array size decrease factor

  !> Coupler object definition
  type :: lptcoupler

    ! This is the name of the coupler
    character(len=str_medium) :: name='UNNAMED_LPTCPL'!< Coupler name (default=UNNAMED_CPL)

    ! Random integer to identify this coupler
    integer :: chk

    ! source and destination objects
    ! we allow two operation modes, one in which the destination is another lpt and one in which
    ! the destination is an array; the mode is determine automatically from which set_dst is
    ! dispatched
    logical :: array_mode
    class(pgrid), pointer :: dst_pg
    type(lpt), pointer :: src, dst

    ! destination sgrid (needed on all processes for particle localization)
    class(sgrid), pointer :: dsg
    integer, dimension(:,:,:), allocatable :: urankmap !< map from destination cell to ugrp rank
    integer, dimension(:,:,:), allocatable :: drankmap !< map from destination cell to dgrp rank

    ! bounding corners of overlap region
    real(WP), dimension(3) :: olapmin, olapmax

    ! Logicals to help us know if we have received a src or dst grid
    logical :: initialized = .false.  ! default value here skips finalizer
    logical :: have_src, have_dst

    ! This is our communication information
    type(MPI_COMM) :: comm                            !< Intracommunicator over (at least) the union of both groups
    type(MPI_GROUP) :: sgrp, dgrp, ugrp               !< Source and destination groups and their union
    integer :: unp, snp, dnp, urank, srank, drank     !< number of processors and ranks on each group
    integer :: uroot, sroot, droot                    !< union rank of the root of each group

    ! Memory adaptation parameters and send/receive buffers
    !   there are gaps in the send buffer        (offsets at multiples of dnp)
    !   there are not gaps in the recieve buffer (offsets at cumsum(recvcounts))
    integer :: sendbufsize, recvbufsize               !< Current allocated buffer sizes
    integer, dimension(:), allocatable :: sendcounts  !< Send counts (grp size)
    integer, dimension(:), allocatable :: recvcounts  !< Recieve counts (grp size)
    type(part), dimension(:), pointer :: sendbuf      !< Send buffer (sendbufsize times dgrp size)
    type(part), dimension(:), pointer :: recvbuf      !< Receive buffer (recvbufsize)

    ! Particle array used for array mode
    integer :: pulledparticlecount
    type(part), dimension(:), pointer :: pulledparticles

    ! Miscellaneous behavior flags for debugging or specialized use cases
    !   - to skip removal of particles from overlap domain in dst, set
    !     dstflag = 0 but leave dontsync = .false.
    !   - if for some reason you only want one global copy of each particle
    !     that *moves* from src to dst, set srcflag = 1 and
    !     dstflag = 0
    integer :: srcflag = 0                            !< what to set the flag to when removing a particle
    integer :: dstflag = 1                            !< what to set the flag to when removing a particle
    logical :: dontsync = .false.                     !< skip syncing/recycling in dst
    integer, dimension(:), allocatable :: pushflagfilter  !< flags to skip when pushing

  contains

    procedure :: set_src                              !< Sets the source
    generic :: set_dst => set_dst_lpt, set_dst_arr    !< Sets the destination
    procedure, private :: set_dst_lpt
    procedure, private :: set_dst_arr
    procedure :: initialize                           !< Allocate buffers
    procedure :: in_overlap                           !< Check if a point is in the overlap region
    procedure :: push                                 !< Src routine that pushes a field into our send data storage
    procedure :: pull                                 !< Dst routine that pulls a field from our received data storage
    procedure :: transfer                             !< Routine that performs the src->dst data transfer
    final :: finalize                                 !< Deallocate buffers

  end type lptcoupler

  !> Declare constructors
  !TODO interface lptcoupler; procedure construct_from_coupler; end interface lptcoupler;
  interface lptcoupler; procedure construct_from_two_groups; end interface lptcoupler;

contains

  !> Coupler constructor from two groups
  function construct_from_two_groups(src_grp, dst_grp, name) result(self)
    use messager, only: die
    use parallel, only: comm
    use mpi_f08, only: mpi_group_union, mpi_comm_create_group,                &
      mpi_comm_group, mpi_group_rank, mpi_group_size, mpi_comm_rank,          &
      mpi_bcast
    implicit none
    intrinsic random_number
    type(lptcoupler) :: self
    type(MPI_GROUP), intent(in) :: src_grp, dst_grp
    character(len=*), intent(in) :: name
    integer, dimension(2) :: ranks
    double precision :: rf
    integer :: ierr

    ! Set name for the coupler
    self%name=trim(adjustl(name))

    ! Build group union
    self%sgrp = src_grp; self%dgrp = dst_grp;
    call mpi_group_union(self%sgrp, self%dgrp, self%ugrp, ierr)

    ! Set ranks and number of processors
    call mpi_group_size(self%sgrp, self%snp, ierr)
    call mpi_group_size(self%dgrp, self%dnp, ierr)
    call mpi_group_rank(self%sgrp, self%srank, ierr)
    call mpi_group_rank(self%dgrp, self%drank, ierr)

    ! Create intracommunicator for the new group
    call mpi_comm_create_group(comm, self%ugrp, 0, self%comm, ierr)
    !call mpi_comm_group(self%comm, self%ugrp, ierr)

    ! get ranks and number of processors in this new group
    call mpi_group_size(self%ugrp, self%unp, ierr)
    call mpi_group_rank(self%ugrp, self%urank, ierr)
    call mpi_group_translate_ranks(self%sgrp, 1, (/ 0 /), self%ugrp, ranks(1:1))
    call mpi_group_translate_ranks(self%dgrp, 1, (/ 0 /), self%ugrp, ranks(2:2))
    self%uroot = 0; self%sroot = ranks(1); self%droot = ranks(2);

    ! check ranks in the union group and the union communicator are the same
    call check_ranks_equal(self%comm, self%ugrp)

    ! Default to no src or dst
    self%have_src = .false.; self%have_dst = .false.;
    self%initialized = .false.

    ! Set check
    if (self%urank .eq. 0) then
      call random_number(rf)
      self%chk = floor(2**15 * rf)
    end if
    call mpi_bcast(self%chk, 1, MPI_INTEGER, 0, self%comm)

  end function construct_from_two_groups

  !> Set the source lpt
  subroutine set_src(this, ps)
    use messager, only: warn
    implicit none
    class(lptcoupler), intent(inout) :: this
    class(lpt), target, intent(in) :: ps

    if (this%have_src) call warn('[lptcoupler] source grid has already been set')

    this%src => ps; this%have_src = .true.; this%initialized = .false.;

  end subroutine set_src

  !> Set the destination lpt
  subroutine set_dst_lpt(this, ps)
    use messager, only: warn
    implicit none
    class(lptcoupler), intent(inout) :: this
    class(lpt), target, intent(in) :: ps

    if (this%have_dst) call warn('[lptcoupler] destination grid has already been set')

    this%array_mode = .false.; this%have_dst = .true.; this%initialized = .false.;

    this%dst => ps; this%dst_pg => ps%cfg;

    ! we check that the this%dst_pg%comm rank and the
    ! this%dst_pg%group ranks are equal
    call check_ranks_equal(this%dst_pg%comm, this%dst_pg%group)

  end subroutine set_dst_lpt
  subroutine set_dst_arr(this, pg)
    use messager, only: warn
    implicit none
    class(lptcoupler), intent(inout) :: this
    class(pgrid), target, intent(in) :: pg

    if (this%have_dst) call warn('[lptcoupler] destination grid has already been set')

    this%array_mode = .true.; this%have_dst = .true.; this%initialized = .false.;

    this%dst => null(); this%dst_pg => pg;

    ! we check that the this%dst_pg%comm rank and the
    ! this%dst_pg%group ranks are equal
    call check_ranks_equal(this%dst_pg%comm, this%dst_pg%group)

  end subroutine set_dst_arr

  !> Function to check ranks are equal between a group and a communicator
  subroutine check_ranks_equal(comm, group)
    use mpi_f08, only: mpi_comm_rank, mpi_group_rank
    use messager, only: die
    implicit none
    type(MPI_COMM), intent(in) :: comm
    type(MPI_GROUP), intent(in) :: group
    integer :: rank1, rank2, ierr

    call mpi_comm_rank(comm, rank1, ierr)
    call mpi_group_rank(group, rank2, ierr)
    if (rank1 .ne. rank2) call die('[lptcoupler] rank mismatch between cgrp and comm')

  end subroutine check_ranks_equal

  !> Test if a position is in the overlap region
  function in_overlap(this, x) result(I)
    implicit none
    class(lptcoupler), intent(in) :: this
    real(WP), dimension(3) :: x
    logical :: I

    !TODO match equality here with what's done in lpt
    I = all(this%olapmin .le. x) .and. all(x .le. this%olapmax)

  end function in_overlap

  !> Allocate buffers
  subroutine initialize(this)
    use messager, only: warn, die
    use mpi_f08,  only: mpi_allreduce, MPI_MIN, MPI_MAX, mpi_bcast
    implicit none
    class(lptcoupler), intent(inout) :: this
    real(WP), dimension(3) :: bcorner, tcorner
    integer :: i, j, k, ierr

    if (this%srank .ne. MPI_UNDEFINED .and. .not. this%have_src)              &
      call die('[lptcoupler] source grid has not been set')
    if (this%drank .ne. MPI_UNDEFINED .and. .not. this%have_dst)              &
      call die('[lptcoupler] destination grid has not been set')

    if (this%initialized) then
      call warn('[lptcoupler] already initialized')
      deallocate(this%sendcounts, this%sendbuf, this%recvbuf, this%urankmap,  &
        this%drankmap)
    end if

    ! find overlap
    bcorner(:) = -huge(bcorner); tcorner(:) = huge(tcorner);
    if (this%srank .eq. 0) then
      bcorner(1) = max(this%src%cfg%x(this%src%cfg%imin  ), bcorner(1))
      bcorner(2) = max(this%src%cfg%y(this%src%cfg%jmin  ), bcorner(2))
      bcorner(3) = max(this%src%cfg%z(this%src%cfg%kmin  ), bcorner(3))
      tcorner(1) = min(this%src%cfg%x(this%src%cfg%imax+1), tcorner(1))
      tcorner(2) = min(this%src%cfg%y(this%src%cfg%jmax+1), tcorner(2))
      tcorner(3) = min(this%src%cfg%z(this%src%cfg%kmax+1), tcorner(3))
    end if
    if (this%drank .eq. 0) then
      bcorner(1) = max(this%dst_pg%x(this%dst_pg%imin  ), bcorner(1))
      bcorner(2) = max(this%dst_pg%y(this%dst_pg%jmin  ), bcorner(2))
      bcorner(3) = max(this%dst_pg%z(this%dst_pg%kmin  ), bcorner(3))
      tcorner(1) = min(this%dst_pg%x(this%dst_pg%imax+1), tcorner(1))
      tcorner(2) = min(this%dst_pg%y(this%dst_pg%jmax+1), tcorner(2))
      tcorner(3) = min(this%dst_pg%z(this%dst_pg%kmax+1), tcorner(3))
    end if
    call mpi_allreduce(bcorner, this%olapmin, 3, MPI_REAL_WP, MPI_MAX,        &
      this%comm, ierr)
    call mpi_allreduce(tcorner, this%olapmax, 3, MPI_REAL_WP, MPI_MIN,        &
      this%comm, ierr)
    if (any(this%olapmin .ge. this%olapmax))                                  &
      call die('[lptcoupler] no overlap between source and destination grids')

    ! set destination grid on all processes
    bcast_dsg : block
      integer, dimension(5) :: intparams
      real(WP), dimension(:), pointer :: x, y, z
      logical, dimension(3) :: pers
      type(sgrid), pointer :: sgptr
      if (this%urank .eq. this%droot) then
        intparams = (/ this%dst_pg%coordsys, this%dst_pg%no,                  &
          this%dst_pg%nx, this%dst_pg%ny, this%dst_pg%nz /)
        pers = (/ this%dst_pg%xper, this%dst_pg%yper, this%dst_pg%zper /)
        x => this%dst_pg%x(this%dst_pg%imin:this%dst_pg%imax+1)
        y => this%dst_pg%y(this%dst_pg%jmin:this%dst_pg%jmax+1)
        z => this%dst_pg%z(this%dst_pg%kmin:this%dst_pg%kmax+1)
      end if
      call mpi_bcast(intparams, 5, MPI_INTEGER, this%droot, this%comm, ierr)
      call mpi_bcast(pers, 3, MPI_LOGICAL, this%droot, this%comm, ierr)
      if (this%urank .ne. this%droot)                                         &
        allocate(x(intparams(3)+1), y(intparams(4)+1), z(intparams(5)+1))
      call mpi_bcast(x, intparams(3)+1, MPI_REAL_WP, this%droot, this%comm,   &
        ierr)
      call mpi_bcast(y, intparams(4)+1, MPI_REAL_WP, this%droot, this%comm,   &
        ierr)
      call mpi_bcast(z, intparams(5)+1, MPI_REAL_WP, this%droot, this%comm,   &
        ierr)
      if (this%drank .ne. MPI_UNDEFINED) then
        this%dsg => this%dst_pg
      else
        allocate(sgptr)
        sgptr = sgrid(intparams(1), intparams(2), x, y, z, pers(1), pers(2),  &
          pers(3), 'dstgcopy')
        this%dsg => sgptr
      end if
      if (this%urank .ne. this%droot) deallocate(x, y, z)
    end block bcast_dsg

    ! set rank map from destination index to global rank
    allocate(this%urankmap(this%dsg%imino:this%dsg%imaxo,                     &
      this%dsg%jmino:this%dsg%jmaxo,this%dsg%kmino:this%dsg%kmaxo))
    allocate(this%drankmap(this%dsg%imino:this%dsg%imaxo,                     &
      this%dsg%jmino:this%dsg%jmaxo,this%dsg%kmino:this%dsg%kmaxo))
    if (this%drank .eq. 0) then
      do k = this%dst_pg%kmino, this%dst_pg%kmaxo
        do j = this%dst_pg%jmino, this%dst_pg%jmaxo
          do i = this%dst_pg%imino, this%dst_pg%imaxo
            this%drankmap(i,j,k) = this%dst_pg%get_rank((/ i, j, k /))
          end do
          call mpi_group_translate_ranks(this%dgrp, this%dst_pg%nxo,          &
            this%drankmap(:,j,k), this%ugrp, this%urankmap(:,j,k))
        end do
      end do
    end if
    call mpi_bcast(this%urankmap, this%dsg%nxo*this%dsg%nyo*this%dsg%nzo,     &
      MPI_INTEGER, this%droot, this%comm, ierr)
    call mpi_bcast(this%drankmap, this%dsg%nxo*this%dsg%nyo*this%dsg%nzo,     &
      MPI_INTEGER, this%droot, this%comm, ierr)

    ! allocate memory
    allocate(this%sendcounts(this%unp), this%recvcounts(this%unp))
    call this%push(only_count=.true.)
    this%sendbufsize = max(ceiling(MEM_ADJ_UP * maxval(this%sendcounts)), 1)
    allocate(this%sendbuf(this%sendbufsize*this%dnp))
    this%recvbufsize = 0; allocate(this%recvbuf(this%recvbufsize));

    this%initialized = .true.

  end subroutine initialize

  !> Deallocate buffers
  subroutine finalize(this)
    implicit none
    type(lptcoupler), intent(inout) :: this

    if (.not. this%initialized) return

    deallocate(this%urankmap, this%drankmap, this%sendcounts,                 &
      this%recvcounts, this%sendbuf, this%recvbuf)

    if (this%drank .eq. MPI_UNDEFINED) deallocate(this%dsg)

    this%initialized = .false.

  end subroutine finalize

  !> Routine that puts overlapping particles into the send buffer and updates send counts
  !> As a design decision this doesn't take an lpt object as we have a pointer
  ! from initialization, but the functionality is the same as in the coupler
  ! case; lpt is not modified, the send buffer is updated, and no communication
  ! is conducted.
  subroutine push(this, only_count)
    implicit none
    class(lptcoupler), intent(inout) :: this
    logical, intent(in), optional :: only_count
    logical :: only_count_actual
    integer :: i, j, k, urank, drank, newbufsize
    integer, dimension(this%dnp) :: drank_to_urank
    integer, dimension(3) :: dstind
    type(part), dimension(:), pointer :: newbuf

    only_count_actual = .false.
    if (present(only_count)) only_count_actual = only_count

    this%sendcounts(:) = 0

    if (this%srank .eq. MPI_UNDEFINED) return

    call mpi_group_translate_ranks(this%dgrp, this%dnp,                       &
      (/ (i, i = 1, this%dnp) /) - 1, this%ugrp, drank_to_urank)

    do i = 1, this%src%np_
      if (.not. this%in_overlap(this%src%p(i)%pos)) cycle
      if (allocated(this%pushflagfilter)) then
        if (any(this%src%p(i)%flag .eq. this%pushflagfilter)) cycle
      end if
      dstind = this%dsg%get_ijk_global(this%src%p(i)%pos)
      urank = this%urankmap(dstind(1), dstind(2), dstind(3))
      drank = this%drankmap(dstind(1), dstind(2), dstind(3))
      if (.not. only_count_actual) then
        if (this%sendcounts(urank+1) + 1 .gt. this%sendbufsize) then
          newbufsize = ceiling((maxval(this%sendcounts) + 1) * MEM_ADJ_UP)
          allocate(newbuf(this%dnp * newbufsize))
          do j = 0, this%dnp - 1
            do k = 1, this%sendcounts(drank_to_urank(j+1)+1)
              newbuf(newbufsize * j + k) = this%sendbuf(this%sendbufsize * j + k)
            end do
          end do
          deallocate(this%sendbuf); this%sendbuf => newbuf;
          this%sendbufsize = newbufsize;
          call checksendbuffer(this, 'resize')
        end if
        j = this%sendbufsize * drank + this%sendcounts(urank+1) + 1
        this%sendbuf(j) = this%src%p(i)
        this%sendbuf(j)%ind = dstind
        this%src%p(i)%flag = this%srcflag
      end if
      this%sendcounts(urank+1) = this%sendcounts(urank+1) + 1
    end do

    if (.not. only_count_actual) call checksendbuffer(this, 'endpush')

  end subroutine push

  !> Check send buffer contains particles meant for communication
  subroutine checksendbuffer(this, event)
    use messager, only: die
    implicit none
    class(lptcoupler), intent(in) :: this
    character(len=*), intent(in) :: event
    integer :: j, k
    real(WP), dimension(3) :: x
    integer, dimension(3) :: ind
    integer, dimension(this%dnp) :: drank_to_urank

    call mpi_group_translate_ranks(this%dgrp, this%dnp,                       &
      (/ (j, j = 1, this%dnp) /) - 1, this%ugrp, drank_to_urank)

    do j = 0, this%dnp - 1
      do k = 1, this%sendcounts(drank_to_urank(j+1)+1)
        x(:) = this%sendbuf(this%sendbufsize * j + k)%pos
        if (.not. this%in_overlap(x))                                         &
          call die("Particle in send buffer not in overlap. (" // event // ")")
        ind = this%dsg%get_ijk_global(x)
        if (this%drankmap(ind(1), ind(2), ind(3)) .ne. j)                     &
          call die("Particle has incorrect drank. (" // event // ")")
      end do
    end do

  end subroutine checksendbuffer

  !> Routine that replaces particles in the overlap region with those from the coupler
  subroutine pull(this)
    implicit none
    class(lptcoupler), intent(inout) :: this
    integer :: i, recv_tot

    if (this%drank .eq. MPI_UNDEFINED) return

    recv_tot = sum(this%recvcounts)

    ! TODO debug
    do i = 1, recv_tot
      this%recvbuf(i)%ind = this%dst_pg%get_ijk_global(this%recvbuf(i)%pos)
      if (this%dst_pg%rank .ne. this%dst_pg%get_rank(this%recvbuf(i)%ind)) then
        write(*,*) 'ERROR: particle sent to wrong process'
        if (this%in_overlap(this%recvbuf(i)%pos)) then
          write(*,*) '       ', 'is in overlap region'
        else
          write(*,*) '       ', 'is not in overlap region'
        end if
        write(*,*) '       ', this%recvbuf(i)%ind, this%dst_pg%get_ijk_global(this%recvbuf(i)%pos)
        write(*,*) '       ', this%recvbuf(i)%pos, this%olapmin, this%olapmax
        write(*,*) '       ', this%dst_pg%rank, this%dst_pg%get_rank(this%recvbuf(i)%ind)
      end if
    end do

    if (this%array_mode) then

      this%pulledparticlecount = recv_tot
      this%pulledparticles => this%recvbuf(1:recv_tot)

    else

      !TODO do we want this dstflag conditional?
      if (this%dstflag .ne. 0) then
        do i = 1, this%dst%np_
          if (this%in_overlap(this%dst%p(i)%pos)) then
            this%dst%p(i)%flag = this%dstflag
          end if
        end do
      end if

      if (.not. this%dontsync) call this%dst%recycle()

      call this%dst%resize(this%dst%np_+recv_tot, dontshrink=.true.)
      this%dst%p((this%dst%np_+1):(this%dst%np_+recv_tot)) = this%recvbuf(1:recv_tot)
      this%dst%np_ = this%dst%np_ + recv_tot

      if (.not. this%dontsync) call this%dst%sync()

    end if

  end subroutine pull

  !> Routine that transfers the data from src to dst - both src_group and
  !> dst_group processors need to call
  subroutine transfer(this)
    use mpi_f08,   only: mpi_alltoall, mpi_alltoallv, mpi_bcast
    use lpt_class, only: MPI_PART
    use messager, only: die
    implicit none
    class(lptcoupler), intent(inout) :: this
    integer, dimension(1:this%unp) :: senddisps, recvdisps, dranks
    integer :: i, chk, ierr

    !TODO DEBUG first check send buffer
    call checksendbuffer(this, 'pretransfer')

    !TODO DEBUG send check
    chk = -1
    if (this%urank .eq. 0) chk = this%chk
    call mpi_bcast(chk, 1, MPI_INTEGER, 0, this%comm)
    if (chk .ne. this%chk) call die("[lptcoupler] check does not match")

    ! send sizes
    call mpi_alltoall(this%sendcounts, 1, MPI_INTEGER, this%recvcounts, 1,    &
      MPI_INTEGER, this%comm, ierr)

    !TODO DEBUG check processes that shouldn't be getting particles aren't and vice versa
    if (this%drank .eq. MPI_UNDEFINED .and. this%recvcounts(this%urank+1) .ne. 0) call die("processor not in dst group being sent particles")
    if (this%srank .eq. MPI_UNDEFINED .and. this%sendcounts(this%urank+1) .ne. 0) call die("processor not in src group sending particles")

    ! resize recieve buffer if needed
    if (sum(this%recvcounts) .gt. this%recvbufsize) then
      this%recvbufsize = ceiling(MEM_ADJ_UP * sum(this%recvcounts))
      deallocate(this%recvbuf); allocate(this%recvbuf(this%recvbufsize));
    end if

    ! compute send displacements
    if (this%srank .eq. MPI_UNDEFINED) then
      senddisps(:) = 0
    else
      call mpi_group_translate_ranks(this%ugrp, this%unp,                     &
        (/ (i, i = 1, this%unp) /) - 1, this%dgrp, dranks, ierr)
      senddisps(1) = 0
      do i = 2, this%unp
        senddisps(i) = senddisps(i-1)
        if (dranks(i-1) .ne. MPI_UNDEFINED)                                   &
          senddisps(i) = senddisps(i) + this%sendbufsize
      end do
    end if

    ! compute recieve displacements
    recvdisps(:) = (/ (sum(this%recvcounts(1:i)), i = 0, this%unp - 1) /)

    !TODO DEBUG
    !write(*,"(I3,A,I8,A,I8,A,I8)") this%urank, " - ranks - ", this%srank, ", ", this%drank, ", ", MPI_UNDEFINED
    !write(*,"(I3,A,I8,A,I8)") this%urank, " - bufsizes - ", this%sendbufsize, ", ", this%recvbufsize
    !write(*,"(I3,A,*(I8))") this%urank, " - sendcounts - ", this%sendcounts
    !write(*,"(I3,A,*(I8))") this%urank, " - senddisps - ",       senddisps
    !write(*,"(I3,A,*(I8))") this%urank, " - recvcounts - ", this%recvcounts
    !write(*,"(I3,A,*(I8))") this%urank, " - recvdisps - ",       recvdisps

    ! send particles
    call mpi_alltoallv(this%sendbuf, this%sendcounts, senddisps, MPI_PART,    &
                       this%recvbuf, this%recvcounts, recvdisps, MPI_PART,    &
                       this%comm, ierr)

  end subroutine transfer

end module lptcoupler_class

