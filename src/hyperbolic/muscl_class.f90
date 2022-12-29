!>  MUSCL-type solver class
!>  Provides support for various BC, generic hyperbolic structures, and source
!>  terms projected onto eigenstructures
!>
!>  Originally written by John P Wakefield, December 2022
!>
!>  this file includes:
!>
!>    - module block
!>        - function interfaces
!>        - bc type
!>        - solver class def
!>        - system class defs
!>    - contains (with headers containing exactly these phrases)
!>        - standard interface class functions
!>        - muscl-specific class functions
!>        - limiter definitions
!>
!>  When possible, helper functions immediately follow the class function that
!>  calls them.
!>
module muscl_class
  use precision,      only: WP
  use string,         only: str_medium
  use config_class,   only: config
  use iterator_class, only: iterator
  implicit none
  private

  ! expose type/constructor/methods
  public :: muscl, constructor, bcond

  ! expose interfaces
  public :: eigenvals_ftype, rsolver_ftype, limiter_ftype

  ! limiter names
  integer(1), parameter, public :: UPWIND   = 0_1
  integer(1), parameter, public :: LAXWEND  = 1_1
  integer(1), parameter, public :: BEAMWARM = 2_1
  integer(1), parameter, public :: FROMM    = 3_1
  integer(1), parameter, public :: MINMOD   = 4_1
  integer(1), parameter, public :: SUPERBEE = 5_1
  integer(1), parameter, public :: MC       = 6_1
  integer(1), parameter, public :: VANLEER  = 7_1

  ! List of known available bcond types for this solver
  ! here 'right' means bcond direction +1
  ! bits 1 and 2 - 00 nothing, 10 interpolate, 11 reflect interpolate
  ! bits 3 and 4 - 00 nothing, 01 zero normal velocity
  ! bit 5 - right moving waves, 0 leave, 1 zero
  ! bit 6 - left moving waves, 0 leave, 1 zero
  integer(1), parameter, public :: wall_reflect = 7_1   !< reflecting wall
  integer(1), parameter, public :: wall_absorb = 39_1   !< absorbing wall
  integer(1), parameter, public :: no_waves = 34_1      !< outflow/free condition condition

  ! riemann solutions for an N dimensional system are stored in N by N+4
  ! dimensional arrays, in which
  ! rs(:, 1)       left-going eigenvalues
  ! rs(:, 2)       right-going eigenvalues
  ! rs(:, 3)       wave strengths
  ! rs(:, 4)       (projected) source strengths
  ! rs(:, 5:(N+4)) (linearized) eigenvectors

  interface
    pure subroutine eigenvals_ftype(N, params, u, evals)
      use precision, only: WP
      implicit none
      integer, intent(in) :: N
      real(WP), dimension(:), intent(in) :: params
      real(WP), dimension(N), intent(in) :: u
      real(WP), dimension(N), intent(out) :: evals
    end subroutine
    pure subroutine rsolver_ftype(N, pl, ul, pr, ur, rs)
      use precision, only: WP
      implicit none
      integer, intent(in) :: N
      real(WP), dimension(:), intent(in) :: pl, pr
      real(WP), dimension(N), intent(in) :: ul, ur
      real(WP), dimension(:,:), intent(out) :: rs
    end subroutine
    pure function limiter_ftype(r) result(phi)
      use precision, only: WP
      implicit none
      real(WP), intent(in) :: r
      real(WP) :: phi
    end function
  end interface

  !> boundary conditions for the hyperbolic solver
  type :: bcond
    type(bcond), pointer :: next                              !< linked list of bconds
    character(len=str_medium) :: name = 'UNNAMED_BCOND'       !< bcond name (default UNNAMED_BCOND)
    integer(1) :: type                                        !< bcond type
    type(iterator) :: itr                                     !< this is the iterator for the bcond - this identifies the (i,j,k)
    character(len=1) :: face                                  !< bcond face (x/y/z)
    integer(1) :: dir                                         !< bcond direction (+1,-1,0 for interior)
  end type bcond

  !> solver object definition
  type :: muscl

    ! config / pgrid object
    class(config), pointer :: cfg                             !< config / pgrid information

    ! name
    character(len=str_medium) :: name = 'UNNAMED_MUSCL'       !< solver name (default UNNAMED_MUSCL)

    ! system information
    integer :: N                                              !< system dimension
    integer :: P                                              !< number of parameters
    procedure(eigenvals_ftype), pointer, nopass :: evals_x, evals_y, evals_z
    procedure(rsolver_ftype), pointer, nopass :: rsolv_x, rsolv_y, rsolv_z
    !TODO all the reflecting bc stuff
    logical, dimension(:), allocatable :: reflect_mask_x, reflect_mask_y, reflect_mask_z

    ! limiting
    procedure(limiter_ftype), pointer, nopass :: limiter      !< limiter function
    real(WP) :: upratio_divzero_eps                           !< softening epsilon for upwind ratio

    ! boundary condition list
    integer :: nbc                                            !< number of bcond for our solver
    real(WP), dimension(:), allocatable :: mfr                !< mFR through each bcond
    real(WP), dimension(:), allocatable :: area               !< area for each bcond
    real(WP) :: correctable_area                              !< area of bcond that can be corrected
    type(bcond), pointer :: first_bc                          !< list of bcond for our solver

    ! flow variables, parameter arrays
    real(WP), dimension(:,:,:,:), pointer :: Uc, dU       !< state variables
    real(WP), dimension(:,:,:,:), pointer :: params       !< params for evals/rsolver

    !TODO get this working
    ! transsonic flags
    !integer(1), dimension(:,:,:), allocatable :: trans        !< bitwise entropy violation check xyz
    !integer :: trans_total
    !logical :: have_trans_flags

    ! wave bcs
    ! first 6 bits used, same convention as type, xxyyzz
    integer(1), dimension(:,:,:), allocatable :: wavebcs

    ! CFL numbers
    real(WP) :: CFL_x, CFL_y, CFL_z                       !< global CFL numbers (over dt)
    real(WP) :: CFL_safety                                !< CFL modifier when using cached value
    logical :: have_CFL_x_estim, have_CFL_y_estim, have_CFL_z_estim
    logical :: have_CFL_x_exact, have_CFL_y_exact, have_CFL_z_exact

    ! monitoring quantities
    real(WP), dimension(:), pointer :: Umin, Umax         !< state variable range
    logical :: have_Urange

  contains

    ! standard interface
    procedure :: print => muscl_print                         !< output solver to the screen
    procedure :: add_bcond                                    !< add a boundary condition
    procedure :: get_bcond                                    !< get a boundary condition
    procedure :: apply_bcond                                  !< apply all boundary conditions
    procedure :: get_cfl                                      !< get maximum CFL
    procedure :: get_range                                    !< calculate min/max field values
    procedure :: get_min => get_range                         !< compatibility
    procedure :: get_max => get_range                         !< compatibility
    !procedure :: get_mfr                                      !< mfr at ea bcond in last step
    procedure :: compute_dU                                   !< take step
    procedure :: compute_dU_x, compute_dU_y, compute_dU_z     !< take step

    ! muscl specific
    !TODO get this working
    !procedure :: check_transonic                              !< check for transonic waves
    procedure :: recalc_cfl                                   !< calculate maximum CFL

  end type muscl

  !> declare constructors
  interface muscl; procedure constructor; end interface

contains

  !! standard interface class functions

  !> default constructor for MUSCL-type flow solver
  function constructor(cfg, name, N, P, evals_x, evals_y, evals_z, rsolv_x,   &
      & rsolv_y, rsolv_z, lim, upratio_divzero_eps, CFL_safety) result(this)
    use messager, only: die
    implicit none
    type(muscl) :: this
    class(config), target, intent(in) :: cfg
    character(len=*), intent(in), optional :: name
    integer, intent(in) :: N, P
    procedure(eigenvals_ftype), pointer, intent(in) :: evals_x, evals_y, evals_z
    procedure(rsolver_ftype), pointer, intent(in) :: rsolv_x, rsolv_y, rsolv_z
    integer(1), intent(in) :: lim
    real(WP), intent(in) :: upratio_divzero_eps, CFL_safety
    integer :: imino, imaxo, jmino, jmaxo, kmino, kmaxo, minmino, maxmaxo

    ! check number of ghost cells
    if (cfg%no .lt. 2) call die("must have at least two ghost cells")

    ! point to pgrid object
    this%cfg => cfg

    ! set the name for the solver
    if (present(name)) this%name = trim(adjustl(name))

    ! system information
    this%N = N; this%P = P;
    this%evals_x => evals_x; this%evals_y => evals_y; this%evals_z => evals_z;
    this%rsolv_x => rsolv_x; this%rsolv_y => rsolv_y; this%rsolv_z => rsolv_z;

    ! limiting
    select case (lim)
    case (UPWIND)
      this%limiter => limiter_upwind
    case (LAXWEND)
      this%limiter => limiter_laxwend
    case (BEAMWARM)
      this%limiter => limiter_beamwarm
    case (FROMM)
      this%limiter => limiter_fromm
    case (MINMOD)
      this%limiter => limiter_minmod
    case (SUPERBEE)
      this%limiter => limiter_superbee
    case (MC)
      this%limiter => limiter_mc
    case (VANLEER)
      this%limiter => limiter_vanleer
    case default
      call die("could not find limiter")
    end select
    this%upratio_divzero_eps = upratio_divzero_eps

    ! boundary condition list
    this%nbc = 0
    this%first_bc => NULL()

    ! get array sizes
    imino = this%cfg%imino_; imaxo = this%cfg%imaxo_;
    jmino = this%cfg%jmino_; jmaxo = this%cfg%jmaxo_;
    kmino = this%cfg%kmino_; kmaxo = this%cfg%kmaxo_;
    minmino = min(imino, jmino, kmino); maxmaxo = max(imaxo, jmaxo, kmaxo)

    ! allocate flow variables and parameter arrays
    allocate(this%params(P,imino:imaxo,jmino:jmaxo,kmino:kmaxo))
    allocate(this%Uc(N,imino:imaxo,jmino:jmaxo,kmino:kmaxo))
    allocate(this%dU(N,imino:imaxo,jmino:jmaxo,kmino:kmaxo))

    ! allocate transsonic flags, initialize values
    !allocate(this%trans(imino:imaxo,jmino:jmaxo,kmino:kmaxo))
    !this%have_trans_flags = .false.

    ! wave bcs flags
    allocate(this%wavebcs(imino:imaxo,jmino:jmaxo,kmino:kmaxo))
    this%wavebcs(:,:,:) = 0_1

    ! initialize CFL computations
    this%CFL_safety = CFL_safety
    this%have_CFL_x_estim = .false.; this%have_CFL_x_exact = .false.;
    this%have_CFL_y_estim = .false.; this%have_CFL_y_exact = .false.;
    this%have_CFL_z_estim = .false.; this%have_CFL_z_exact = .false.;

    ! initialize monitoring
    allocate(this%Umin(N), this%Umax(N))
    this%have_Urange = .false.

  end function

  !> solver print function
  subroutine muscl_print(this)
    implicit none
    class(muscl), intent(in) :: this

    !TODO



  end subroutine

  !> Add a boundary condition
  subroutine add_bcond(this, name, type, locator, face, dir)
    use string,   only: lowercase
    use messager, only: die
    implicit none
    class(muscl), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(1), intent(in) :: type
    external :: locator
    interface
      logical function locator(pargrid,ind1,ind2,ind3)
        use pgrid_class, only: pgrid
        class(pgrid), intent(in) :: pargrid
        integer, intent(in) :: ind1,ind2,ind3
      end function locator
    end interface
    character(len=1), intent(in) :: face
    integer(1), intent(in) :: dir
    type(bcond), pointer :: new_bc
    integer :: i, j, k, n
    integer(1) :: fi, wbc

    ! check input
    select case (lowercase(face))
      case ('x'); case ('y'); case ('z');
      case default; call die('[muscl add_bcond] Unknown bcond face - expecting x, y, or z')
    end select
    if (abs(dir) .ne. 1 .and. dir .ne. 0) then
      call die('[muscl add_bcond] Unknown bcond dir - expecting -1, +1, or 0')
    end if

    ! prepare new bcond
    allocate(new_bc)
    new_bc%name = trim(adjustl(name))
    new_bc%type = type
    new_bc%face = lowercase(face)
    fi = int(ichar(new_bc%face) - ichar('x'), 1)
    new_bc%itr = iterator(pg=this%cfg, name=new_bc%name, locator=locator,     &
      & face=new_bc%face)
    new_bc%dir = dir

    ! insert it in front
    new_bc%next => this%first_bc
    this%first_bc => new_bc
    this%nbc = this%nbc + 1

    ! update wave zeroing array
    wbc = ishft(ior(type, 96_1), -4_1)
    if (dir .eq. -1_1) then
      select case (wbc)
      case (1)
        wbc = 2_1
      case (2)
        wbc = 1_1
      case default
      end select
    end if
    wbc = ishft(wbc, 2 * fi)
    do n = 1,new_bc%itr%n_
      i = new_bc%itr%map(1,n); j = new_bc%itr%map(2,n); k = new_bc%itr%map(3,n)
      this%wavebcs(i,j,k) = this%wavebcs(i,j,k) + wbc
    end do

  end subroutine add_bcond

  !> get a boundary condition
  subroutine get_bcond(this,name,my_bc)
    use messager, only: die
    implicit none
    class(muscl), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(bcond), pointer, intent(out) :: my_bc

    my_bc => this%first_bc
    do while (associated(my_bc))
      if (trim(my_bc%name).eq.trim(name)) return
      my_bc => my_bc%next
    end do

    call die('[muscl get_bcond] Boundary condition was not found')

  end subroutine get_bcond

  !> enforce boundary condition
  !> note that this interpolates (ensuring ghost cells have the right value),
  !> but all wave-related constraints are imposed during the step
  !TODO
  subroutine apply_bcond(this, t, dt)
    use messager, only: die
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), intent(in) :: t, dt
    integer(1) :: masked_type
    type(bcond), pointer :: my_bc

    ! Traverse bcond list
    my_bc => this%first_bc
    do while (associated(my_bc))

      ! Only processes inside the bcond work here
      if (my_bc%itr%amIn) then

        ! Select appropriate action based on the bcond type
        masked_type = ior(my_bc%type, 3_1)
        select case (masked_type)

          case (0)                      !< do nothing

          case (2)                      !< interpolate



          case (3)                      !< reflect interpolate



          case (1)
          call die('[muscl apply_bcond] Unknown bcond type')
          case default 
          call die('[muscl apply_bcond] Unknown bcond type')
        end select

        !<<<<<<<<< TODO example of how this might be done
        !! Implement based on bcond direction
        !select case (my_bc%face)
        !case ('x')
        !  stag=min(my_bc%dir,0)
        !  do n=1,my_bc%itr%n_
        !    i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
        !    this%U(i     ,j    ,k    )=this%U(i-my_bc%dir     ,j    ,k    )
        !    this%V(i+stag,j:j+1,k    )=this%V(i-my_bc%dir+stag,j:j+1,k    )
        !    this%W(i+stag,j    ,k:k+1)=this%W(i-my_bc%dir+stag,j    ,k:k+1)
        !  end do
        !case ('y')
        !  stag=min(my_bc%dir,0)
        !  do n=1,my_bc%itr%n_
        !    i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
        !    this%U(i:i+1,j+stag,k    )=this%U(i:i+1,j-my_bc%dir+stag,k    )
        !    this%V(i    ,j     ,k    )=this%V(i    ,j-my_bc%dir     ,k    )
        !    this%W(i    ,j+stag,k:k+1)=this%W(i    ,j-my_bc%dir+stag,k:k+1)
        !  end do
        !case ('z')
        !  stag=min(my_bc%dir,0)
        !  do n=1,my_bc%itr%n_
        !    i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
        !    this%U(i:i+1,j    ,k+stag)=this%U(i:i+1,j    ,k-my_bc%dir+stag)
        !    this%V(i    ,j:j+1,k+stag)=this%V(i    ,j:j+1,k-my_bc%dir+stag)
        !    this%W(i    ,j    ,k     )=this%W(i    ,j    ,k-my_bc%dir     )
        !  end do
        !end select
        !! If needed, clip
        !if (my_bc%type.eq.clipped_neumann) then
        !  select case (my_bc%face)
        !  case ('x')
        !    do n=1,my_bc%itr%n_
        !      i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
        !      if (this%U(i,j,k)*my_bc%rdir.lt.0.0_WP) this%U(i,j,k)=0.0_WP
        !    end do
        !  case ('y')
        !    do n=1,my_bc%itr%n_
        !      i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
        !      if (this%V(i,j,k)*my_bc%rdir.lt.0.0_WP) this%V(i,j,k)=0.0_WP
        !    end do
        !  case ('z')
        !    do n=1,my_bc%itr%n_
        !      i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
        !      if (this%W(i,j,k)*my_bc%rdir.lt.0.0_WP) this%W(i,j,k)=0.0_WP
        !    end do
        !  end select
        !end if

        !>>>>>>>>>>>>>> end example

      end if

      ! Sync full fields after each bcond - this should be optimized
      call this%cfg%sync(this%Uc)
      call this%cfg%sync(this%dU)

      ! Move on to the next bcond
      my_bc => my_bc%next

    end do

  end subroutine apply_bcond

  !> get CFL number
  subroutine get_cfl(this, dt, cfl)
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP), intent(out) :: cfl
    logical :: have_CFL_estim, have_CFL_exact

    have_CFL_estim = this%have_CFL_x_estim .and. this%have_CFL_y_estim .and.  &
      & this%have_CFL_z_estim

    have_CFL_exact = this%have_CFL_x_exact .and. this%have_CFL_y_exact .and.  &
      this%have_CFL_z_exact

    if (.not. have_CFL_estim) call this%recalc_cfl()

    cfl = dt * max(this%CFL_x, this%CFL_y, this%CFL_z)

    if (.not. have_CFL_exact) cfl = this%CFL_safety * cfl

  end subroutine get_cfl

  !> get min/max field values
  subroutine get_range(this)
    use messager, only: die
    use mpi_f08,  only: MPI_ALLREDUCE, MPI_MIN, MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), dimension(this%N) :: task_Umin, task_Umax
    integer :: i, j, k, ierr

    if (.not. this%have_Urange) then

      task_Umin = + huge(1.0_WP)
      task_Umax = - huge(1.0_WP)

      do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
            task_Umin = min(task_Umin, this%Uc(:,i,j,k))
            task_Umax = max(task_Umax, this%Uc(:,i,j,k))
          end do
        end do
      end do

      call MPI_ALLREDUCE(task_Umin, this%Umin, this%N, MPI_REAL_WP, MPI_MIN,  &
        & this%cfg%comm, ierr)
      if (ierr .ne. 0) call die("error reducing Umin")
      call MPI_ALLREDUCE(task_Umax, this%Umax, this%N, MPI_REAL_WP, MPI_MAX,  &
        & this%cfg%comm, ierr)
      if (ierr .ne. 0) call die("error reducing Umax")

      this%have_Urange = .true.

    end if

  end subroutine get_range

  !> get mass flow rate through boundaries at previous step
  !TODO

  !> compute dU
  ! this is an unsplit step---its use is generally a bad idea.  A better
  ! approach is to use the directional subroutines below in operator splitting;
  ! e.g. Strang.  Strang stepping should utilize Jack's timestepping class.
  subroutine compute_dU(this, dt)
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), intent(in) :: dt

    this%dU(:, :, :, :) = 0.0_WP

    call this%compute_dU_x(dt)

    call this%compute_dU_y(dt)

    call this%compute_dU_z(dt)

  end subroutine compute_dU

  !> compute dU in x direction
  subroutine compute_dU_x(this, dt)
    use messager, only: die
    use mpi_f08,  only: MPI_ALLREDUCE, MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP) :: cfl
    integer :: N, P, M, j, k, ierr

    call this%cfg%sync(this%Uc)

    N = this%N; P = this%P; M = size(this%Uc, 2); cfl = 0.0_WP;

    do k = this%cfg%kmin_,this%cfg%kmax_
      do j = this%cfg%jmin_,this%cfg%jmax_
        call wavestep_1d(N, P, M, this%limiter, this%upratio_divzero_eps,     &
          & this%rsolv_x, this%cfg%dx, this%wavebcs(:,j,k), dt,               &
          & this%params(:,:,j,k), this%Uc(:,:,j,k), this%dU(:,:,j,k), cfl)
      end do
    end do

    call this%cfg%sync(this%dU)

    call MPI_ALLREDUCE(cfl, this%CFL_x, 1, MPI_REAL_WP, MPI_MAX,              &
      & this%cfg%comm, ierr)
    if (ierr .ne. 0) call die("error reducing CFL")

    this%have_CFL_x_estim = .true.; this%have_CFL_x_exact = .false.;

  end subroutine compute_dU_x

  !> compute dU in y direction
  subroutine compute_dU_y(this, dt)
    use messager, only: die
    use mpi_f08,  only: MPI_ALLREDUCE, MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP) :: cfl
    integer :: N, P, M, i, k, jmino_, jmaxo_, ierr
    real(WP), dimension(:,:), allocatable :: U1d, dU1d, p1d
    integer(1), dimension(:), allocatable :: bc1d

    call this%cfg%sync(this%Uc)

    N = this%N; P = this%P; cfl = 0.0_WP;
    jmino_ = this%cfg%jmino_; jmaxo_ = this%cfg%jmaxo_;

    allocate(U1d(N,jmino_:jmaxo_), p1d(P,jmino_:jmaxo_))
    allocate(dU1d(N,jmino_:jmaxo_), bc1d(jmino_:jmaxo_))

    M = size(U1d, 2);
    do k = this%cfg%kmin_, this%cfg%kmax_
      do i = this%cfg%imin_, this%cfg%imax_
        U1d = this%Uc(:,i,jmino_:jmaxo_,k)
        dU1d = this%dU(:,i,jmino_:jmaxo_,k)
        p1d = this%params(:,i,jmino_:jmaxo_,k)
        bc1d = this%wavebcs(i,jmino_:jmaxo_,k)
        call wavestep_1d(N, P, M, this%limiter, this%upratio_divzero_eps,     &
          & this%rsolv_y, this%cfg%dy, bc1d, dt, p1d, U1d, dU1d, cfl)
        this%dU(:,i,jmino_:jmaxo_,k) = dU1d
      end do
    end do
    deallocate(U1d, dU1d, p1d, bc1d)

    call this%cfg%sync(this%dU)

    call MPI_ALLREDUCE(cfl, this%CFL_y, 1, MPI_REAL_WP, MPI_MAX,              &
      & this%cfg%comm, ierr)
    if (ierr .ne. 0) call die("error reducing CFL")

    this%have_CFL_y_estim = .true.; this%have_CFL_y_exact = .false.;

  end subroutine compute_dU_y

  !> compute dU in z direction
  subroutine compute_dU_z(this, dt)
    use messager, only: die
    use mpi_f08,  only: MPI_ALLREDUCE, MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP) :: cfl
    integer :: N, P, M, i, j, kmino_, kmaxo_, ierr
    real(WP), dimension(:,:), allocatable :: U1d, dU1d, p1d
    integer(1), dimension(:), allocatable :: bc1d

    call this%cfg%sync(this%Uc)

    N = this%N; P = this%P; cfl = 0.0_WP;
    kmino_ = this%cfg%kmino_; kmaxo_ = this%cfg%kmaxo_;

    allocate(U1d(N,kmino_:kmaxo_), p1d(P,kmino_:kmaxo_))
    allocate(dU1d(N,kmino_:kmaxo_), bc1d(kmino_:kmaxo_))

    M = size(U1d, 2);
    do j = this%cfg%jmin_, this%cfg%jmax_
      do i = this%cfg%imin_, this%cfg%imax_
        U1d = this%Uc(:,i,j,kmino_:kmaxo_)
        dU1d = this%dU(:,i,j,kmino_:kmaxo_)
        p1d = this%params(:,i,j,kmino_:kmaxo_)
        bc1d = this%wavebcs(i,j,kmino_:kmaxo_)
        call wavestep_1d(N, P, M, this%limiter, this%upratio_divzero_eps,     &
          & this%rsolv_z, this%cfg%dz, bc1d, dt, p1d, U1d, dU1d, cfl)
        this%dU(:,i,j,kmino_:kmaxo_) = dU1d
      end do
    end do
    deallocate(U1d, dU1d, p1d, bc1d)

    call this%cfg%sync(this%dU)

    call MPI_ALLREDUCE(cfl, this%CFL_z, 1, MPI_REAL_WP, MPI_MAX,              &
      & this%cfg%comm, ierr)
    if (ierr .ne. 0) call die("error reducing CFL")

    this%have_CFL_z_estim = .true.; this%have_CFL_z_exact = .false.;

  end subroutine compute_dU_z

  ! requires two ghost cells even if using the `upwind' limiter
  ! left as an independent function, but called is by class
  pure subroutine wavestep_1d(N, P, M, limfun, eps, rsolver, dxs, wbcs, dt,   &
      & params, U, dU, CFLmax)
    implicit none
    integer, intent(in) :: N, P, M
    procedure(limiter_ftype), pointer, intent(in) :: limfun
    real(WP), intent(in) :: eps
    procedure(rsolver_ftype), pointer, intent(in) :: rsolver
    real(WP), dimension(M), intent(in) :: dxs                   ! size M
    integer(1), dimension(:), intent(in) :: wbcs                ! size MI
    real(WP), intent(in) :: dt
    real(WP), dimension(P,M), intent(in) :: params              ! size (P, M)
    real(WP), dimension(N,M), intent(in) :: U                   ! size (N, M)
    real(WP), dimension(N,M), intent(inout) :: dU               ! size (N, M)
    real(WP), intent(inout) :: CFLmax
    real(WP), dimension(N,N+4) :: ll, lc, rc, rr                ! size (N, N+4)
    real(WP), dimension(N) :: phil, phir
    integer :: j

    call rsolver(N, params(:,1), U(:,1), params(:,2), U(:,2), ll)
    call rsolver(N, params(:,2), U(:,2), params(:,3), U(:,3), lc)
    call rsolver(N, params(:,3), U(:,3), params(:,4), U(:,4), rc)

    !TODO this is broken
    !call handle_wavebc(wbcs(1), ll)
    !call handle_wavebc(wbcs(2), lc)
    !call handle_wavebc(wbcs(3), rc)

    CFLmax = max(CFLmax, maxval(abs(ll(:,1:2))) / min(dxs(1), dxs(2)))
    CFLmax = max(CFLmax, maxval(abs(lc(:,1:2))) / min(dxs(2), dxs(3)))
    CFLmax = max(CFLmax, maxval(abs(rc(:,1:2))) / min(dxs(3), dxs(4)))

    call compute_limval(N, limfun, eps, ll, lc, rc, phil)

    do j = 3, M-2
      call rsolver(N, params(:,j+1), U(:,j+1), params(:,j+2), U(:,j+2), rr)
      !TODO this is broken
      !call handle_wavebc(wbcs(j+1), rr)
      call compute_limval(N, limfun, eps, lc, rc, rr, phir)
      CFLmax = max(CFLmax, maxval(abs(rr(:,1:2))) / min(dxs(j+1), dxs(j+2)))
      call wavestep_firstorder_single(N, dxs(j), dt, lc, rc, dU(:,j))
      call wavestep_limitedcorrection_single(N, dxs(j), dt, lc, rc, &
        & phil, phir, dU(:, j))
      ll(:, :) = lc(:, :)
      lc(:, :) = rc(:, :)
      rc(:, :) = rr(:, :)
      phil = phir
    end do

  end subroutine wavestep_1d

  pure subroutine handle_wavebc(bc, rs)
    implicit none
    integer(1), intent(in) :: bc
    real(WP), dimension(:,:), intent(inout) :: rs
    integer :: N

    !TODO this is broken

    N = size(rs, 1)

    ! zero negative waves?
    if (ior(bc, 1_1) .eq. 1_1) rs(:, 3:4) = max(rs(:, 3:4), 0.0_WP)

    ! zero positive waves?
    if (ior(bc, 2_1) .eq. 2_1) rs(:, 3:4) = min(rs(:, 3:4), 0.0_WP)

  end subroutine handle_wavebc

  ! this is called by the class, but is left as an independent function
  pure subroutine compute_limval(N, limfun, eps, l_rs, c_rs, r_rs, limvals)
    implicit none
    procedure(limiter_ftype), pointer, intent(in) :: limfun
    real(WP), intent(in) :: eps
    integer, intent(in) :: N
    real(WP), dimension(:,:), intent(in) :: l_rs, c_rs, r_rs
    real(WP), dimension(:), intent(out) :: limvals
    real(WP) :: r
    integer :: p

    do p = 1,N
        !TODO replace this with eigenvector projections
        r = c_rs(p, 1) + c_rs(p, 2)
        if (r .gt. 0.0_WP) then
          r = l_rs(p, 3) * c_rs(p, 3) / (c_rs(p, 3) * c_rs(p, 3) + eps)
        else
          r = r_rs(p, 3) * c_rs(p, 3) / (c_rs(p, 3) * c_rs(p, 3) + eps)
        end if
        limvals(p) = limfun(r)
    end do

  end subroutine compute_limval

  ! standard Roe/MUSCL/Leveque style wave step with limiters
  ! the step has been split up to improve performance
  ! this function has been made painful to look at in the interest of speed in
  ! the julia implementation, then it stayed that way because I didn't want to
  ! rewrite it
  ! this is called by the class, but is left as an independent function
  pure subroutine wavestep_firstorder_single(N, dx, k, lc_rs, rc_rs, dU)
    implicit none
    integer, intent(in) :: N
    real(WP), intent(in) :: dx, k
    real(WP), dimension(N,N+4), intent(in) :: lc_rs, rc_rs
    real(WP), dimension(N), intent(inout) :: dU
    real(WP), dimension(N) :: tmp1, tmp2

    ! waves from left
    tmp1 = lc_rs(:, 2) * lc_rs(:, 3)
    !l_dxbeta : block
    !  integer :: j
    !  tmp2 = lc_rs(:, 1) + lc_rs(:, 2)
    !  ! this probably isn't necessary for anything physical
    !  do j = 1, N
    !    if (tmp2(j) .eq. 0.0_WP) then
    !      tmp2(j) = 0.5_WP
    !    else if (tmp2(j) .gt. 0.0_WP) then
    !      tmp2(j) = 1.0_WP
    !    else
    !      tmp2(j) = 0.0_WP
    !    end if
    !  end do
    !  tmp2 = tmp2 * lc_rs(:, 4)
    !  tmp2 = dx * tmp2
    !end block l_dxbeta
    !tmp1 = tmp1 - tmp2
    tmp2 = matmul(lc_rs(:, 5:(N+4)), tmp1)
    tmp2 = k / dx * tmp2
    dU = dU - tmp2

    ! waves from right
    tmp1 = rc_rs(:, 1) * rc_rs(:, 3)
    !r_dxbeta : block
    !  integer :: j
    !  tmp2 = rc_rs(:, 1) + rc_rs(:, 2)
    !  do j = 1, N
    !    if (tmp2(j) .eq. 0.0_WP) then
    !      tmp2(j) = 0.5_WP
    !    else if (tmp2(j) .lt. 0.0_WP) then
    !      tmp2(j) = 1.0_WP
    !    else
    !      tmp2(j) = 0.0_WP
    !    end if
    !  end do
    !  tmp2 = tmp2 * rc_rs(:, 4)
    !  tmp2 = dx * tmp2
    !end block r_dxbeta
    !tmp1 = tmp1 - tmp2
    tmp2 = matmul(rc_rs(:, 5:(N+4)), tmp1)
    tmp2 = k / dx * tmp2
    dU = dU - tmp2

  end subroutine wavestep_firstorder_single

  ! this is called by the class, but is left as an independent function
  pure subroutine wavestep_limitedcorrection_single(N, dx, k, lc_rs,   &
      & rc_rs, phil, phir, dU)
    implicit none
    integer, intent(in) :: N
    real(WP), intent(in) :: dx, k
    real(WP), dimension(:,:), intent(in) :: lc_rs, rc_rs
    real(WP), dimension(:), intent(in) :: phil, phir
    real(WP), dimension(:), intent(inout) :: dU
    real(WP), dimension(N) :: tmp1, tmp2

    tmp1 = (1.0_WP - k / dx * (rc_rs(:,2) - rc_rs(:,1)))
    tmp1 = tmp1 * (rc_rs(:,2) - rc_rs(:,1)) * rc_rs(:, 3) * phir(:)
    tmp2 = matmul(rc_rs(:, 5:(N+4)), tmp1)
    tmp2 = - 0.5_WP * k / dx * tmp2
    dU = dU + tmp2

    tmp1 = (1.0_WP - k / dx * (lc_rs(:,2) - lc_rs(:,1)))
    tmp1 = tmp1 * (lc_rs(:,2) - lc_rs(:,1)) * lc_rs(:, 3) * phil(:)
    tmp2 = matmul(lc_rs(:, 5:(N+4)), tmp1)
    tmp2 = 0.5_WP * k / dx * tmp2
    dU = dU + tmp2

  end subroutine

  !! muscl-specific class functions

  !> check for transsonic rarefaction waves
  !TODO
! pure subroutine check_transonic(this)
!   implicit none
!   class(muscl), intent(inout) :: this
!   real(WP) :: dim_CFL
!   real(WP), dimension(N) :: evals_l, evals_r
!   integer :: imin, imax, jmin, jmax, kmin, kmax
!   integer :: imino, imaxo, jmino, jmaxo, kmino, kmaxo

!   this%trans_flags(:,:,:) = 0
!   this%trans_total = 0

!   imino = this%cfg%imino_; imaxo = this%cfg%imaxo_;
!   jmino = this%cfg%jmino_; jmaxo = this%cfg%jmaxo_;
!   kmino = this%cfg%kmino_; kmaxo = this%cfg%kmaxo_;
!   imin = this%cfg%imin_; imax = this%cfg%imax_;
!   jmin = this%cfg%jmin_; jmax = this%cfg%jmax_;
!   kmin = this%cfg%kmin_; kmax = this%cfg%kmax_;

!   dim_CFL = 0.0_WP
!   do k=this%cfg%kmin_,this%cfg%kmax_
!     do j=this%cfg%jmin_,this%cfg%jmax_
!       call check_transsonic_1d(this%evals_x, 1, this%params(:,imino:imaxo,j,k), &
!         &this%Uc(:,imino:imaxo,j,k), flags(imino:imaxo-1,j,k), this%trans_total, &
!         &dim_CFL)
!     end do
!   end do

!   this%loc_CFL_x = dim_CFL

!   dim_CFL = 0.0_WP
!   do k = this%cfg%kmin_, this%cfg%kmax_
!     do i = this%cfg%imin_, this%cfg%imax_
!       this%Ucache1(:,jmino:jmaxo) = this%Uc(:,i,jmino:jmaxo,:)
!       this%pcache(:,jmino:jmaxo) = this%params(:,i,jmino:jmaxo,:)
!       this%tcache(:) = 0
!       call check_transsonic_1d(this%evals_y, 2, this%pcache(:,jmino:jmaxo),     &
!         &this%Ucache1(:,jmino:jmaxo), this%tcache(jmino:jmaxo-1),               &
!         &this%trans_total, dim_CFL)
!       flags(i,jmin:jmax-1) = flags(i,jmin:jmax-1) + this%tcache(jmin:jmax-1)
!     end do
!   end do

!   this%loc_CFL_y = dim_CFL

!   dim_CFL = 0.0_WP
!   do j = this%cfg%jmin_, this%cfg%jmax_
!     do i = this%cfg%imin_, this%cfg%imax_
!       this%Ucache1(:,kmino:kmaxo) = this%Uc(:,i,:,kmino:kmaxo)
!       this%pcache(:,kmino:kmaxo) = this%params(:,i,:,kmino:kmaxo)
!       this%tcache(:) = 0
!       call check_transsonic_1d(this%evals_z, 4, this%pcache(:,kmino:kmaxo),     &
!         &this%Ucache1(:,kmino:kmaxo), this%tcache(kmino:kmaxo-1),               &
!         &this%trans_total, dim_CFL)
!       flags(i,kmin:kmax-1) = flags(i,kmin:kmax-1) + this%tcache(kmin:kmax-1)
!     end do
!   end do

!   this%loc_CFL_z = dim_CFL

!   this%have_trans_flags = .true.
!   this%have_CFL_estim = .true.
!   this%have_CFL_exact = .true.
!   this%CFL_reduced = .false.

! end subroutine check_transsonic

  ! called by class, but kept as independent function
! pure subroutine check_transsonic_1d(evalf, flaginc, dxs, params, U, flags, total, CFL)
!   procedure(eigenvals_ftype), intent(in) :: evalf
!   integer(1), intent(in) :: flaginc
!   integer, intent(in) :: P, N, M, MI
!   integer, intent(in) :: M2, M3
!   real(WP), dimension(M), intent(in) :: dxs
!   real(WP), dimension(P,M2), intent(in) :: params
!   real(WP), dimension(N,M3), intent(in) :: U
!   integer(1), dimension(MI), intent(inout) :: flags
!   integer, intent(inout) :: total
!   real(WP), intent(inout) :: CFL
!   real(WP), dimension(N) :: evals
!   logical, dimension(N) :: gtzl, gtzr
!   integer :: i

!   call evalf(U(:,1), params(:,1), evals)
!   do i = 1, N; gtzr(i) = evals(i) .gt. 0.0_WP; end do

!   do i = 2, M
!     CFL = max(CFL, maxval(abs(dxs(i) * evals)))
!     gtzl(:) = gtzr(:)
!     call evalf(U(:,1), params(:,1), evals)
!     gtzr(:) = evals(:) .gt. 0.0_WP
!     if (any((.not. gtzl) .and. gtzr)) then
!       flags(i-1) = flags(i-1) + flaginc
!       total = total + 1
!     end if
!   end do

! end subroutine check_transsonic_1d

  !> calculate CFL from cell values
  subroutine recalc_cfl(this)
    use messager, only: die
    use mpi_f08,  only: MPI_ALLREDUCE, MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(muscl), intent(inout) :: this
    real(WP), dimension(this%N) :: evals
    real(WP) :: task_CFL_x, task_CFL_y, task_CFL_z
    integer :: i, j, k, ierrx, ierry, ierrz

    call this%cfg%sync(this%Uc)

    task_CFL_x = 0.0_WP; task_CFL_y = 0.0_WP; task_CFL_z = 0.0_WP;
    do k = this%cfg%kmin_, this%cfg%kmax_
      do j = this%cfg%jmin_, this%cfg%jmax_
        do i = this%cfg%imin_, this%cfg%imax_
          call this%evals_x(this%N, this%params(:,i,j,k), this%Uc(:,i,j,k), evals)
          task_CFL_x = max(task_CFL_x, this%cfg%dxi(i) * maxval(abs(evals)))
          call this%evals_y(this%N, this%params(:,i,j,k), this%Uc(:,i,j,k), evals)
          task_CFL_y = max(task_CFL_y, this%cfg%dyi(i) * maxval(abs(evals)))
          call this%evals_z(this%N, this%params(:,i,j,k), this%Uc(:,i,j,k), evals)
          task_CFL_z = max(task_CFL_z, this%cfg%dzi(i) * maxval(abs(evals)))
        end do
      end do
    end do

    call MPI_ALLREDUCE(task_CFL_x, this%CFL_x, 1, MPI_REAL_WP, MPI_MAX,       &
      & this%cfg%comm, ierrx)
    call MPI_ALLREDUCE(task_CFL_y, this%CFL_y, 1, MPI_REAL_WP, MPI_MAX,       &
      & this%cfg%comm, ierry)
    call MPI_ALLREDUCE(task_CFL_z, this%CFL_z, 1, MPI_REAL_WP, MPI_MAX,       &
      & this%cfg%comm, ierrz)

    if (ierrx .ne. 0) call die("error reducing CFL")
    if (ierry .ne. 0) call die("error reducing CFL")
    if (ierrz .ne. 0) call die("error reducing CFL")

    this%have_CFL_x_estim = .true.; this%have_CFL_x_exact = .true.;
    this%have_CFL_y_estim = .true.; this%have_CFL_y_exact = .true.;
    this%have_CFL_z_estim = .true.; this%have_CFL_z_exact = .true.;

  end subroutine recalc_cfl

  !! limiter definitions

  ! limiter functions
  pure function limiter_upwind(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = r
    phi = 0.0
  end function

  pure function limiter_laxwend(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = r
    phi = 1.0
  end function

  pure function limiter_beamwarm(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = r
  end function

  pure function limiter_fromm(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = 0.5_WP * (r + 1.0_WP)
  end function

  pure function limiter_minmod(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = max(0.0_WP, min(1.0_WP, r))
  end function

  pure function limiter_superbee(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = max(0.0_WP, min(1.0_WP, 2.0_WP * r), min(2.0_WP, r))
  end function

  pure function limiter_mc(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = max(0.0_WP, min((1.0_WP + r) / 2.0_WP, 2.0_WP, 2.0_WP * r))
  end function

  pure function limiter_vanleer(r) result(phi)
    implicit none
    real(WP), intent(in) :: r
    real(WP) :: phi
    phi = (r + abs(r)) / (1.0_WP + abs(r))
  end function

end module

