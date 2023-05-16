!> Adds particle coupling to an existing coupler
module lptcoupler_class
  use precision,      only: WP
  use mpi_f08,        only: MPI_Comm, MPI_Group
  use string,         only: str_medium
  use sgrid_class,    only: sgrid
  use config_class,   only: config
  use lpt_class,      only: lpt, part
  implicit none
  private

  ! Expose type/constructor/methods
  public :: lptcoupler
  
  ! buffer sizing behavior
  real(WP), parameter :: MEM_ADJ_UP = 1.3_WP          !< Particle array size increase factor
  real(WP), parameter :: MEM_ADJ_DN = 0.7_WP          !< Particle array size decrease factor

  !> Coupler object definition
  type :: lptcoupler

    ! This is the name of the coupler
    character(len=str_medium) :: name='UNNAMED_LPTCPL'!< Coupler name (default=UNNAMED_CPL)

    ! source and destination lpt objects
    type(lpt), pointer :: src, dst

    ! destination pgrid (needed on all processes for particle localization)
    class(sgrid), pointer :: dsg
    integer, dimension(:,:,:), allocatable :: rankmap, drankmap

    ! bounding corners of overlap region
    real(WP), dimension(3) :: olapmin, olapmax

    ! Logicals to help us know if we have received a src or dst grid
    logical :: have_src, have_dst, initialized

    ! This is our communication information
    type(MPI_Comm) :: comm                            !< Intracommunicator over (at least) the union of both groups
    type(MPI_Group) :: cgrp                           !< Group for whole communicator
    type(MPI_Group) :: sgrp, dgrp, grp                !< Source and destination groups and their union
    integer :: np, snp, dnp, rank, srank, drank       !< number of processors and ranks on each group
    integer :: sroot, droot                           !< union rank of source and destination roots

    ! Memory adaptation parameters and send/receive buffers
    !   there are gaps in the send buffer        (offsets at multiples of dnp)
    !   there are not gaps in the recieve buffer (offsets at cumsum(recvcounts))
    integer :: sendbufsize, recvbufsize               !< Current allocated buffer sizes
    integer, dimension(:), allocatable :: sendcounts  !< Send counts (grp size)
    integer, dimension(:), allocatable :: recvcounts  !< Recieve counts (grp size)
    type(part), dimension(:), pointer :: sendbuf      !< Send buffer (sendbufsize times dgrp size)
    type(part), dimension(:), pointer :: recvbuf      !< Receive buffer (recvbufsize)

    ! Miscellaneous behavior flags for debugging or specialized use cases
    !   - to skip removal of particles from overlap domain in dst, set
    !     dstflag = 0 but leave dontsync = .false.
    !   - if for some reason you only want one global copy of each particle
    !     that *moves* from src to dst, set srcflag = 1 and
    !     dstflag = 0
    integer :: srcflag = 0                            !< what to set the flag to when removing a particle
    integer :: dstflag = 1                            !< what to set the flag to when removing a particle
    logical :: dontsync = .false.                     !< skip syncing/recycling in dst

  contains

    procedure :: set_src                              !< Sets the source
    procedure :: set_dst                              !< Sets the destination
    procedure :: initialize                           !< Allocate buffers
    procedure :: in_overlap                           !< Check if in overlap region
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
    use mpi_f08, only: MPI_GROUP_UNION, MPI_COMM_CREATE_GROUP, MPI_COMM_GROUP,&
      MPI_GROUP_RANK, MPI_GROUP_SIZE, MPI_COMM_RANK, MPI_UNDEFINED,           &
      MPI_GROUP_TRANSLATE_RANKS
    implicit none
    type(lptcoupler) :: self
    type(MPI_Group), intent(in) :: src_grp, dst_grp
    character(len=*), intent(in) :: name
    integer, dimension(2) :: ranks
    integer :: ierr

    ! Set name for the coupler
    self%name=trim(adjustl(name))

    ! Build group union
    self%sgrp = src_grp; self%dgrp = dst_grp;
    call MPI_GROUP_UNION(self%sgrp, self%dgrp, self%grp, ierr)

    ! Set ranks and number of processors
    call MPI_GROUP_SIZE(self%grp, self%np, ierr)
    call MPI_GROUP_SIZE(self%sgrp, self%snp, ierr)
    call MPI_GROUP_SIZE(self%dgrp, self%dnp, ierr)
    call MPI_GROUP_RANK(self%grp, self%rank, ierr)
    call MPI_GROUP_RANK(self%sgrp, self%srank, ierr)
    call MPI_GROUP_RANK(self%dgrp, self%drank, ierr)
    call MPI_GROUP_TRANSLATE_RANKS(self%sgrp, 1, (/ 0 /), self%grp, ranks(1:1))
    call MPI_GROUP_TRANSLATE_RANKS(self%dgrp, 1, (/ 0 /), self%grp, ranks(2:2))
    self%sroot = ranks(1); self%droot = ranks(2);

    ! Create intracommunicator for the new group
    call MPI_COMM_CREATE_GROUP(comm, self%grp, 0, self%comm, ierr)
    call MPI_COMM_GROUP(self%comm, self%cgrp, ierr)

    ! Default to no src or dst
    self%have_src = .false.; self%have_dst = .false.;
    self%initialized = .false.

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
  subroutine set_dst(this, ps)
    use mpi_f08, only: MPI_COMM_RANK
    use messager, only: warn
    implicit none
    class(lptcoupler), intent(inout) :: this
    class(lpt), target, intent(in) :: ps

    if (this%have_dst) call warn('[lptcoupler] destination grid has already been set')

    this%dst => ps; this%have_dst = .true.; this%initialized = .false.;

    ! we check that the this%dst%cfg%comm rank and the
    ! this%dst%cfg%group ranks are equal
    call check_ranks_equal(this%dst%cfg%comm, this%dst%cfg%group)

  end subroutine set_dst

  !> Function to check ranks are equal between a group and a communicator
  subroutine check_ranks_equal(comm, group)
    use mpi_f08, only: MPI_COMM_RANK, MPI_GROUP_RANK
    use messager, only: die
    implicit none
    type(MPI_Comm), intent(in) :: comm
    type(MPI_Group), intent(in) :: group
    integer :: rank1, rank2, ierr

    call MPI_COMM_RANK(comm, rank1, ierr)
    call MPI_GROUP_RANK(group, rank2, ierr)
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
    use mpi_f08,  only: MPI_UNDEFINED, mpi_allreduce, MPI_MIN, MPI_MAX,       &
      MPI_INTEGER, MPI_LOGICAL, mpi_bcast
    use parallel, only: MPI_REAL_WP
    implicit none
    class(lptcoupler), intent(inout) :: this
    real(WP), dimension(3) :: bcorner, tcorner
    type(sgrid), pointer :: sgptr
    integer :: i, j, k, ierr

    if (this%srank .ne. MPI_UNDEFINED .and. .not. this%have_src) call die('[lptcoupler] source grid has not been set')
    if (this%drank .ne. MPI_UNDEFINED .and. .not. this%have_dst) call die('[lptcoupler] destination grid has not been set')

    if (this%initialized) then
      call warn('[lptcoupler] already initialized')
      deallocate(this%sendcounts, this%sendbuf, this%recvbuf)
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
      bcorner(1) = max(this%dst%cfg%x(this%dst%cfg%imin  ), bcorner(1))
      bcorner(2) = max(this%dst%cfg%y(this%dst%cfg%jmin  ), bcorner(2))
      bcorner(3) = max(this%dst%cfg%z(this%dst%cfg%kmin  ), bcorner(3))
      tcorner(1) = min(this%dst%cfg%x(this%dst%cfg%imax+1), tcorner(1))
      tcorner(2) = min(this%dst%cfg%y(this%dst%cfg%jmax+1), tcorner(2))
      tcorner(3) = min(this%dst%cfg%z(this%dst%cfg%kmax+1), tcorner(3))
    end if
    call MPI_ALLREDUCE(bcorner, this%olapmin, 3, MPI_REAL_WP, MPI_MAX, this%comm, ierr)
    call MPI_ALLREDUCE(tcorner, this%olapmax, 3, MPI_REAL_WP, MPI_MIN, this%comm, ierr)
    if (any(this%olapmin .gt. this%olapmax)) then
      call die('[lptcoupler] no overlap between source and destination grids')
    end if

    ! set destination grid on all processes
    bcast_dsg : block
      integer, dimension(5) :: intparams
      real(WP), dimension(:), pointer :: x, y, z
      logical, dimension(3) :: pers
      if (this%rank .eq. this%droot) then
        intparams(1) = this%dst%cfg%coordsys
        intparams(2) = this%dst%cfg%no
        intparams(3) = this%dst%cfg%nx
        intparams(4) = this%dst%cfg%ny
        intparams(5) = this%dst%cfg%nz
        pers(1) = this%dst%cfg%xper
        pers(2) = this%dst%cfg%yper
        pers(3) = this%dst%cfg%zper
        x => this%dst%cfg%x(this%dst%cfg%imin:this%dst%cfg%imax+1)
        y => this%dst%cfg%y(this%dst%cfg%jmin:this%dst%cfg%jmax+1)
        z => this%dst%cfg%z(this%dst%cfg%kmin:this%dst%cfg%kmax+1)
      end if
      call mpi_bcast(intparams, 5, MPI_INTEGER, this%droot, this%comm, ierr)
      call mpi_bcast(pers, 3, MPI_LOGICAL, this%droot, this%comm, ierr)
      if (this%rank .ne. this%droot)                                          &
        allocate(x(intparams(3)+1), y(intparams(4)+1), z(intparams(5)+1))
      call mpi_bcast(x, intparams(3)+1, MPI_REAL_WP, this%droot, this%comm, ierr)
      call mpi_bcast(y, intparams(4)+1, MPI_REAL_WP, this%droot, this%comm, ierr)
      call mpi_bcast(z, intparams(5)+1, MPI_REAL_WP, this%droot, this%comm, ierr)
      if (this%drank .ne. MPI_UNDEFINED) then
        this%dsg => this%dst%cfg
      else
        allocate(sgptr)
        sgptr = sgrid(intparams(1), intparams(2), x, y, z, pers(1), pers(2),  &
          pers(3), 'dstgcopy')
        this%dsg => sgptr
      end if
      if (this%rank .ne. this%droot) deallocate(x, y, z)
    end block bcast_dsg

    ! set rank map from destination index to global and destination ranks
    allocate(                                                                 &
      this%rankmap(this%dsg%imin:this%dsg%imax,this%dsg%jmin:this%dsg%jmax,   &
      this%dsg%kmin:this%dsg%kmax), this%drankmap(this%dsg%imin:this%dsg%imax,&
      this%dsg%jmin:this%dsg%jmax,this%dsg%kmin:this%dsg%kmax)                &
      )
    if (this%drank .eq. 0) then
      do k = this%dst%cfg%kmin, this%dst%cfg%kmax
        do j = this%dst%cfg%jmin, this%dst%cfg%jmax
          do i = this%dst%cfg%imin, this%dst%cfg%imax
            this%rankmap(i,j,k) = this%dst%cfg%get_rank((/ i, j, k /))
          end do
          call mpi_group_translate_ranks(this%grp, this%dst%cfg%nx + 1,       &
            this%rankmap(:,j,k), this%dgrp, this%drankmap(:,j,k), ierr)
        end do
      end do
    end if
    call mpi_bcast(this%rankmap, this%dsg%nx*this%dsg%ny*this%dsg%nz,       &
      MPI_INTEGER, this%droot, this%comm, ierr)

    ! allocate memory
    allocate(this%sendcounts(this%np), this%recvcounts(this%np))
    call this%push(only_count=.true.)
    this%sendbufsize = ceiling(MEM_ADJ_UP * maxval(this%sendcounts))
    ! as an initial guess, assume the number of particles per processor are
    ! similar in the two domains; if not this will be fixed at the first transfer
    this%recvbufsize = this%sendbufsize
    allocate(this%sendbuf(this%sendbufsize*this%dnp), this%recvbuf(this%recvbufsize))

    this%initialized = .true.

    if (this%rank .eq. 0) write(*,*) '[lptcoupler] initialized'

  end subroutine initialize

  !> Deallocate buffers
  subroutine finalize(this)
    implicit none
    type(lptcoupler), intent(inout) :: this

    if (.not. this%initialized) return

    deallocate(this%sendcounts, this%recvcounts, this%sendbuf, this%recvbuf)

  end subroutine finalize

  !> Routine that puts overlapping particles into the send buffer and updates send counts
  !> As a design decision this doesn't take an lpt object as we have a pointer
  ! from initialization, but the functionality is the same as in the coupler
  ! case; lpt is not modified, the send buffer is updated, and no communication
  ! is conducted.
  subroutine push(this, only_count)
    use mpi_f08, only: MPI_UNDEFINED
    implicit none
    class(lptcoupler), intent(inout) :: this
    logical, intent(in), optional :: only_count
    logical :: only_count_actual
    integer :: i, j, rank, drank, os, ns, oe, ne, oldbufsize
    integer, dimension(3) :: dstind
    type(part), dimension(:), pointer :: oldbuf

    only_count_actual = .false.
    if (present(only_count)) only_count_actual = only_count

    this%sendcounts(:) = 0

    if (this%srank .eq. MPI_UNDEFINED) return

    do i = 1, this%src%np_
      if (.not. this%in_overlap(this%src%p(i)%pos)) cycle
      dstind = this%dsg%get_ijk_global(this%src%p(i)%pos)
      rank = this%rankmap(dstind(1), dstind(2), dstind(3))
      drank = this%drankmap(dstind(1), dstind(2), dstind(3))
      this%sendcounts(rank+1) = this%sendcounts(rank+1) + 1
      if (.not. only_count_actual) then
        if (maxval(this%sendcounts) .gt. this%sendbufsize) then
          oldbufsize = this%sendbufsize; oldbuf => this%sendbuf;
          this%sendbufsize = ceiling(maxval(this%sendcounts) * this%dnp * MEM_ADJ_UP)
          nullify(this%sendbuf); allocate(this%sendbuf(this%sendbufsize));
          do j = 1, this%dnp
            os = 1 + oldbufsize * (j - 1); oe = os + oldbufsize - 1;
            ns = 1 + this%sendbufsize * (j - 1); ne = os + oldbufsize - 1;
            this%sendbuf(ns:ne) = oldbuf(os:oe)
          end do
          deallocate(oldbuf)
        end if
        j = this%sendbufsize * drank + this%sendcounts(rank+1)
        this%sendbuf(j) = this%src%p(i)
        this%sendbuf(j)%ind = dstind
        this%src%p(i)%flag = this%srcflag
      end if
    end do

  end subroutine push

  !> Routine that replaces particles in the overlap region with those from the coupler
  subroutine pull(this)
    implicit none
    class(lptcoupler), intent(inout) :: this
    integer :: i, recv_tot

    if (this%dstflag .ne. 0) then
      do i = 1, this%dst%np_
        if (this%in_overlap(this%dst%p(i)%pos)) then
          this%dst%p(i)%flag = this%dstflag
        end if
      end do
    end if

    if (.not. this%dontsync) call this%dst%recycle()

    recv_tot = sum(this%recvcounts)
    call this%dst%resize(this%dst%np_+recv_tot, dontshrink=.true.)
    this%dst%p((this%dst%np_+1):(this%dst%np_+recv_tot)) = this%recvbuf(1:recv_tot)
    this%dst%np_ = this%dst%np_ + recv_tot

    if (.not. this%dontsync) call this%dst%sync()

  end subroutine pull

  !> Routine that transfers the data from src to dst - both src_group and dst_group processors need to call
  subroutine transfer(this)
    use mpi_f08,   only: MPI_UNDEFINED, MPI_INTEGER, MPI_ALLTOALL, MPI_ALLTOALLv
    use lpt_class, only: MPI_PART
    implicit none
    class(lptcoupler), intent(inout) :: this
    integer, dimension(this%np) :: senddisps, recvdisps
    integer :: i, ierr

    ! send sizes
    call MPI_ALLTOALL(this%sendcounts,1,MPI_INTEGER,this%recvcounts,1,MPI_INTEGER,this%comm,ierr)

    ! resize recieve buffer if needed
    if (sum(this%recvcounts) .gt. this%recvbufsize) then
      this%recvbufsize = ceiling(MEM_ADJ_UP * sum(this%recvcounts))
      deallocate(this%recvbuf); allocate(this%recvbuf(this%recvbufsize));
    end if

    ! compute send displacements
    senddisps(1) = 0
    do i = 2, this%np
      senddisps(i) = senddisps(i-1)
      if (this%srank .ne. MPI_UNDEFINED) senddisps(i) = senddisps(i) + this%sendbufsize
    end do

    ! compute recieve displacements
    recvdisps(1) = 0
    recvdisps(2:this%np) = (/ (sum(this%recvcounts(1:i)), i = 1, this%np - 1) /)

    ! send particles    
    call MPI_ALLTOALLv(this%sendbuf, this%sendcounts, senddisps, MPI_PART,    &
                       this%recvbuf, this%recvcounts, recvdisps, MPI_PART,    &
                       this%comm, ierr)

  end subroutine transfer

end module lptcoupler_class
