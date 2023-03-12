!>
!> Module for estimation of closures in Shankar's two-fluid models.
!>
!> This module serves two purposes:
!>    1) estimation of statistics and io
!>    2) management of parameters and integral length scales
!>
!> Originally written by John P Wakefield in Feb 2023
!>
module estimclosures_class
  use precision,      only: WP
  use mathtools,      only: pi
  use fft3d_class,    only: fft3d
  use monitor_class,  only: monitor
  use string,         only: str_medium, str_long
  use sgrid_class,    only: sgrid
  use pgrid_class,    only: pgrid
  use config_class,   only: config
  use lpt_class,      only: lpt
  implicit none
  private

  public :: estimclosures, estimclosures_mesh, FLT_NUM_PARAMS, FLT_GAUSSIAN,  &
    ETAODX, RHOF, FORCE_TIMESCALE

  !> give all filters six parameters, not all of which may be used
  integer, parameter :: FLT_NUM_PARAMS = 6

  !> filter types
  integer, parameter :: FLT_GAUSSIAN = 1

  !> mesh spacing to parameter conversion
  real(WP), parameter :: ETAODX = 0.5_WP

  !> parameter primitive quantities
  real(WP), parameter :: RHOF = 1.0_WP
  real(WP), parameter :: FORCE_TIMESCALE = 1.0_WP

  !> sizes of reported statistics arrays
  integer, parameter :: NUM_S = 14
  integer, parameter :: NUM_M = 11
  integer, parameter :: NUM_B = 0
  integer, parameter :: NUM_H = 0

  !> structure containing statistics to store
  type :: statpoint
    real(WP) :: time
    integer :: step
    real(WP), dimension(NUM_S) :: S       ! spatial
    real(WP), dimension(NUM_M) :: M       ! microscopic
    real(WP), dimension(NUM_B) :: B       ! B eq verification
    real(WP), dimension(NUM_H) :: H       ! Vie/Houssem equation
  end type statpoint

  !> filter information and corresponding monitor file
  type :: filter
    character(len=str_medium) :: out_fname
    type(monitor) :: mon
    type(statpoint) :: d
    integer :: type
    real(WP), dimension(FLT_NUM_PARAMS) :: params
    complex(WP), dimension(:,:,:), allocatable :: flt_f
  contains
    procedure :: setup => filter_setup
    final :: filter_destruct
  end type filter

  !> helper structure to read in filter info from files
  type :: filter_info_row
    character(len=str_medium) :: out_fname
    character(len=str_medium) :: typename
    real(WP), dimension(FLT_NUM_PARAMS) :: params
  end type

  !> parameters are stored in arrays, in the following order:
  !> params_nondim - 4 items - Relambda, Stk, phiinf, Wovk
  !> params_primit - 7 items - rhof, rhop, ktarget, epstarget, nu, dp, g

  !> object managing computation of closures for each filter and navigation
  !> of paramter space
  type, abstract :: estimclosures
    ! configs for simulation
    class(pgrid), pointer :: sim_pg
    ! fft object
    type(fft3d), pointer :: fft
    type(sgrid), pointer :: fft_sg
    type(pgrid), pointer :: fft_pg
    ! work arrays
    complex(WP), dimension(:,:,:), allocatable :: phi_r
    complex(WP), dimension(:,:,:), allocatable :: vel_x_r, vel_y_r, vel_z_r,  &
      work_x_r, work_y_r, work_z_r
    complex(WP), dimension(:,:,:), allocatable :: phi_f,                      &
      vel_x_f, vel_y_f, vel_z_f, work_x_f, work_y_f, work_z_f
    ! filter list
    integer :: num_filters
    type(filter), dimension(:), allocatable :: filters
    ! number of timescales before writing output
    real(WP) :: interval_tinfs
  contains
    procedure :: init => ec_init    ! parent constructor
    ! commented because I'm too lazy to write an interface block
    !procedure, deferred :: get_next_params
    !procedure, deferred :: get_interval
    procedure :: compute_statistics
    procedure :: destruct => ec_destruct
  end type estimclosures

  type, extends(estimclosures) :: estimclosures_mesh
    ! nondim mesh info
    real(WP), dimension(4) :: pmin, pspacing
    integer, dimension(4) :: Icurr, Imax, Idir
    ! primitive quantities that are constant (but not parameters)
    integer :: Np
    real(WP) :: nu
  contains
    procedure :: get_next_params
    procedure :: get_interval
  end type estimclosures_mesh

  interface estimclosures_mesh; procedure ecmesh_from_args; end interface

contains

  subroutine ec_init(ec, sim_pg)
    use param,       only: param_read
    use mpi_f08,     only: mpi_bcast, MPI_REAL, MPI_INTEGER, MPI_CHARACTER, MPI_COMM_WORLD
    use parallel,    only: MPI_REAL_WP, group
    use sgrid_class, only: cartesian
    use parallel,    only: group
    use messager,    only: die
    implicit none
    class(estimclosures), intent(inout) :: ec
    class(pgrid), target, intent(in) :: sim_pg
    real(WP), dimension(3) :: dx
    real(WP), dimension(:), allocatable :: fftxs, fftys, fftzs
    integer, dimension(3) :: FFTN
    character(len=str_medium) :: filterfile
    type(filter_info_row) :: f_info_raw
    integer :: li, hi, lj, hj, lk, hk, i, fh, ierr

    ! store pointer to simulation config
    ec%sim_pg => sim_pg

    ! setup ffts (they need their own, much finer, pgrid)
    call param_read('EC FFT mesh', FFTN)
    dx(:) = (/ sim_pg%xL, sim_pg%yL, sim_pg%zL /) / FFTN
    allocate(fftxs(FFTN(1)+1), fftys(FFTN(2)+1), fftzs(FFTN(3)+1))
    do i = 1, FFTN(1) + 1; fftxs(i) = real(i-1,WP) * dx(1); end do;
    do i = 1, FFTN(2) + 1; fftys(i) = real(i-1,WP) * dx(2); end do;
    do i = 1, FFTN(3) + 1; fftzs(i) = real(i-1,WP) * dx(3); end do;
    allocate(ec%fft_sg, ec%fft_pg, ec%fft)
    ec%fft_sg = sgrid(coord=cartesian, no=0, x=fftxs, y=fftys, z=fftzs,       &
      xper=.true., yper=.true., zper=.true., name='EC_FFT_G')
    ec%fft_pg = pgrid(ec%fft_sg, group, (/ sim_pg%npx, sim_pg%npy,            &
      sim_pg%npz /))
    ec%fft = fft3d(ec%fft_pg)

    ! allocate arrays
    li = ec%fft_pg%imin_; lj = ec%fft_pg%jmin_; lk = ec%fft_pg%kmin_;
    hi = ec%fft_pg%imax_; hj = ec%fft_pg%jmax_; hk = ec%fft_pg%kmax_;
    allocate(ec%phi_r(li:hi, lj:hj, lk:hk))
    allocate(ec%work_x_r(li:hi, lj:hj, lk:hk))
    allocate(ec%work_y_r(li:hi, lj:hj, lk:hk))
    allocate(ec%work_z_r(li:hi, lj:hj, lk:hk))
    allocate(ec%phi_f(li:hi, lj:hj, lk:hk))
    allocate(ec%vel_x_f(li:hi, lj:hj, lk:hk))
    allocate(ec%vel_y_f(li:hi, lj:hj, lk:hk))
    allocate(ec%vel_z_f(li:hi, lj:hj, lk:hk))
    allocate(ec%work_x_f(li:hi, lj:hj, lk:hk))
    allocate(ec%work_y_f(li:hi, lj:hj, lk:hk))
    allocate(ec%work_z_f(li:hi, lj:hj, lk:hk))

    ! read filter info (as root)
    if (sim_pg%rank .eq. 0) then
      call param_read('EC filter list', filterfile)
      open(newunit=fh, file=filterfile, action='READ', access='SEQUENTIAL')
      ! first just count them
      ec%num_filters = 0
      do
        read(fh,*,iostat=ierr) f_info_raw
        if (ierr .ne. 0) exit
        ec%num_filters = ec%num_filters + 1
      end do
      allocate(ec%filters(ec%num_filters))
      rewind(fh)
      do i = 1, ec%num_filters
        read(fh,*,iostat=ierr) f_info_raw
        if (ierr .ne. 0) call die("[EC] error reading filterfile row")
        select case (f_info_raw%typename)
        case ('gaussian', 'GAUSSIAN', 'g', 'G')
          ec%filters(i)%type = FLT_GAUSSIAN
        case default
          call die("[EC] unknown filter type '"//f_info_raw%typename//"'")
        end select
        ec%filters(i)%out_fname = f_info_raw%out_fname
        ec%filters(i)%params(:) = f_info_raw%params
      end do
      close(fh)
    end if

    ! broadcast filter info and setup filters
    call mpi_bcast(ec%num_filters, 1, MPI_INTEGER, 0, sim_pg%comm, ierr)
    if (sim_pg%rank .ne. 0) allocate(ec%filters(ec%num_filters))
    do i = 1, ec%num_filters
      call mpi_bcast(ec%filters(i)%out_fname, str_medium, MPI_CHARACTER, 0,   &
        sim_pg%comm, ierr)
      call mpi_bcast(ec%filters(i)%type, 1, MPI_INTEGER, 0, sim_pg%comm, ierr)
      call mpi_bcast(ec%filters(i)%params, FLT_NUM_PARAMS, MPI_REAL_WP, 0,    &
        sim_pg%comm, ierr)
      call ec%filters(i)%setup(ec%fft)
    end do

  end subroutine ec_init

  function ecmesh_from_args(cfg) result(ec)
    use param, only: param_read
    implicit none
    type(estimclosures_mesh) :: ec
    class(config), intent(in) :: cfg
    real(WP), dimension(4) :: pmax
    real(WP) :: etamin
    integer :: n

    ! init parent
    call ec%init(cfg)

    ! read params
    call param_read('EC min Relambda', ec%pmin(1))
    call param_read('EC max Relambda',    pmax(1))
    call param_read('EC num Relambda', ec%Imax(1))
    call param_read('EC min Stk',      ec%pmin(2))
    call param_read('EC max Stk',         pmax(2))
    call param_read('EC num Stk',      ec%Imax(2))
    call param_read('EC min Wovk',     ec%pmin(4))
    call param_read('EC max Wovk',        pmax(4))
    call param_read('EC num Wovk',     ec%Imax(4))
    call param_read('EC min vf',       ec%pmin(3))
    call param_read('EC max vf',          pmax(3))
    call param_read('EC num vf',       ec%Imax(3))
    call param_read('EC integral timescales', ec%interval_tinfs)

    ! init mesh
    ec%Icurr(:) = 1; ec%Icurr(1) = 0; ec%Idir(:) = +1;
    do n = 1, 4
      if (ec%Imax(n) .gt. 1) then
        ec%pspacing(n) = (pmax(n) - ec%pmin(n)) / (ec%Imax(n) - 1)
      else
        ec%pspacing(n) = 0.0_WP
      end if
    end do

    ! set up parameter computation (mostly computing viscosity)
    etamin = ETAODX * cfg%min_meshsize**(1.0_WP / 3)
    ec%nu = etamin**2 * pmax(1) / sqrt(15.0_WP)

  end function ecmesh_from_args 

  ! assumes out_fname, type, num_params, params are set
  ! assumes fft has been setup
  subroutine filter_setup(flt, fft)
    use messager, only: die
    use mpi_f08, only: MPI_SUM, mpi_allreduce
    use parallel, only: MPI_REAL_WP
    implicit none
    class(filter), intent(inout) :: flt
    class(fft3d), intent(inout) :: fft
    real(WP) :: r2_y, r2_z, a, b, flt_int_l, flt_int_g
    integer :: j, k

    ! setup monitor file
    flt%mon = monitor(fft%pg%amRoot, flt%out_fname)
    call flt%mon%add_column(flt%d%S(1),  'PP')
    call flt%mon%add_column(flt%d%S(2),  'PPUBx')
    call flt%mon%add_column(flt%d%S(3),  'PPUBy')
    call flt%mon%add_column(flt%d%S(4),  'PPUBz')
    call flt%mon%add_column(flt%d%S(5),  'PC2')
    call flt%mon%add_column(flt%d%S(6),  'PCUCx')
    call flt%mon%add_column(flt%d%S(7),  'PCUCy')
    call flt%mon%add_column(flt%d%S(8),  'PCUCz')
    call flt%mon%add_column(flt%d%S(9),  'PC2UCx')
    call flt%mon%add_column(flt%d%S(10), 'PC2UCy')
    call flt%mon%add_column(flt%d%S(11), 'PC2UCz')
    call flt%mon%add_column(flt%d%S(12), 'UBx')
    call flt%mon%add_column(flt%d%S(13), 'UBy')
    call flt%mon%add_column(flt%d%S(14), 'UBz')
    call flt%mon%add_column(flt%d%M(1),  'DP')
    call flt%mon%add_column(flt%d%M(2),  'NP')
    call flt%mon%add_column(flt%d%M(3),  'RP')
    call flt%mon%add_column(flt%d%M(4),  'VF')
    call flt%mon%add_column(flt%d%M(5),  'VVxx')
    call flt%mon%add_column(flt%d%M(6),  'VVxy')
    call flt%mon%add_column(flt%d%M(7),  'VVxz')
    call flt%mon%add_column(flt%d%M(8),  'VVyy')
    call flt%mon%add_column(flt%d%M(9),  'VVyz')
    call flt%mon%add_column(flt%d%M(10), 'VVzz')
    call flt%mon%add_column(flt%d%M(11), 'DRG')

    ! allocate Fourier space filter
    allocate(flt%flt_f(fft%pg%imin_:fft%pg%imax_,fft%pg%jmin_:fft%pg%jmax_,   &
      fft%pg%kmin_:fft%pg%kmax_))

    ! fill real filter; we build fft%pg such that the corner is at origin
    select case (flt%type)
    case (FLT_GAUSSIAN)
      a = (6.0_WP / (pi * flt%params(1)**2))**1.5_WP
      b = -6.0_WP / flt%params(1)**2
      do k = fft%pg%kmin_, fft%pg%kmax_
        r2_z = fft%pg%z(k)**2
        do j = fft%pg%jmin_, fft%pg%jmax_
          r2_y = fft%pg%y(j)**2
          flt%flt_f(:,j,k) = a * exp(b * (r2_z + r2_y +                       &
            fft%pg%x(fft%pg%imin_:fft%pg%imax_)**2))
        end do
      end do
    case default
      call die("[estimclosures] invalid filter type")
    end select

    ! check filter is normalized
    flt_int_l = sum(realpart(flt%flt_f))
    call mpi_allreduce(flt_int_l, flt_int_g, 1, MPI_REAL_WP, MPI_SUM,         &
      fft%pg%comm)
    flt%flt_f = flt%flt_f / flt_int_g

    ! transform filter
    call fft%forward_transform(flt%flt_f)
    if (fft%oddball) then
      flt%flt_f(fft%pg%imin_,fft%pg%jmin_,fft%pg%kmin_) = (1.0_WP, 0.0_WP)
    end if

  end subroutine filter_setup

  ! params_nondim - 4 items - Relambda, Stk, phiinf, Wovk
  ! params_primit - 7 items - rhof, rhop, ktarget, epstarget, nu, dp, g
  subroutine get_next_params(ec, params, done)
    use mathtools, only: pi
    implicit none
    class(estimclosures_mesh), intent(inout) :: ec
    real(WP), dimension(7), intent(out)      :: params
    real(WP), dimension(4)                   :: nondim
    logical, intent(out)                     :: done
    integer :: n

    ! check if we are done
    done = all(ec%Icurr .eq. ec%Imax)

    ! increment first coord
    ec%Icurr(1) = ec%Icurr(1) + ec%Idir(1)

    ! if this pushes it past its bounds, increment the next coordinate instead
    ! and change directions of previous coordinate
    do n = 1, 3
      if (ec%Icurr(n) .lt. 1 .or. ec%Icurr(n) .gt. ec%Imax(n)) then
        ec%Icurr(n) = ec%Icurr(n) - ec%Idir(n)
        ec%Icurr(n+1) = ec%Icurr(n+1) + ec%Idir(n+1)
        ec%Idir(n) = -1 * ec%Idir(n)
      end if
    end do

    ! get nondimensional values
    nondim(:) = ec%pmin(:) + (ec%Icurr(:) - 1) * ec%pspacing(:)

    ! fluid parameters
    ! using the current approach viscosity is fixed throughout; it is computed
    ! at initialization
    params(1) = RHOF
    params(5) = ec%nu
    params(3) = nondim(1)**2 * ec%nu / (10 * FORCE_TIMESCALE)
    params(4) = 2.0_WP / 3 * params(3) / FORCE_TIMESCALE

    ! particle and gravity parameters
    params(6) = (6 * ec%sim_pg%vol_total * nondim(3) / (pi * ec%Np))**(1.0_WP / 3)
    params(2) = params(1) * 18 * nondim(2) * sqrt(ec%nu**3 / params(4)) / params(6)**2
    params(7) = params(2) * params(6)**2 / (18 * params(1) * ec%nu)

  end subroutine get_next_params

  subroutine get_interval(ec, interval)
    implicit none
    class(estimclosures_mesh), intent(in) :: ec
    real(WP), intent(out) :: interval

    interval = ec%interval_tinfs * FORCE_TIMESCALE

  end subroutine get_interval

  subroutine filter_destruct(flt)
    implicit none
    type(filter), intent(inout) :: flt

    deallocate(flt%flt_f)

  end subroutine filter_destruct

  subroutine ec_destruct(ec, dealloc_sim_pg)
    implicit none
    class(estimclosures), intent(inout) :: ec
    logical, intent(in), optional :: dealloc_sim_pg

    if (present(dealloc_sim_pg) .and. dealloc_sim_pg) deallocate(ec%sim_pg)

    deallocate(ec%fft_sg, ec%fft_pg, ec%fft)

    deallocate(ec%phi_r, ec%work_x_r, ec%work_y_r, ec%work_z_r, ec%phi_f,     &
      ec%vel_x_f, ec%vel_y_f, ec%vel_z_f, ec%work_x_f, ec%work_y_f,           &
      ec%work_z_f)

    deallocate(ec%filters)

  end subroutine ec_destruct

  subroutine compute_statistics(ec, time, step, ps)
    use mpi_f08,    only:  MPI_CHARACTER, MPI_INTEGER, mpi_reduce, MPI_SUM
    use parallel,   only: MPI_REAL_WP
    use lpt_class,  only: part
    use messager,   only: die
    implicit none
    class(estimclosures), intent(inout) :: ec
    real(WP), intent(in) :: time
    integer, intent(in) :: step
    class(lpt), intent(in) :: ps
    integer, dimension(3) :: ind, indmin, indmax
    real(WP) :: fft_rescale, pvol, dp_loc, vf_loc, drg_loc, dp, vf, drg
    real(WP), dimension(6) :: meanVV_loc, meanVV
    real(WP), dimension(3) :: tmpvec, dxdydz
    integer :: n, i, j, k, N3

    ! rescaling values
    N3 = ec%fft%pg%nx * ec%fft%pg%ny * ec%fft%pg%nz
    fft_rescale = 1.0_WP
    indmin(:) = (/ ec%fft%pg%imin, ec%fft%pg%jmin, ec%fft%pg%kmin /)
    indmax(:) = (/ ec%fft%pg%imax, ec%fft%pg%jmax, ec%fft%pg%kmax /)
    dxdydz(:) = (/ ec%fft%pg%dx(indmin(1)), ec%fft%pg%dy(indmin(2)), ec%fft%pg%dz(indmin(3)) /)

    ! project particle volumes and velocities to mesh, compute unfiltered
    ! quantities
    !TODO drag
    dp_loc = 0.0_WP; vf_loc = 0.0_WP; drg_loc = 0.0_WP; meanVV_loc(:) = 0.0_WP;
    ec%phi_f(:,:,:) = 0.0_WP;
    ec%vel_x_f = 0.0_WP; ec%vel_y_f = 0.0_WP; ec%vel_z_f = 0.0_WP;
    do n = 1, ps%np_
      ind = min(max(floor(ps%p(n)%pos / dxdydz), indmin), indmax)
      ind(:) = ec%fft%pg%get_ijk_global(ps%p(n)%pos, ind)
      i = ind(1); j = ind(2); k = ind(3);
      dp_loc = dp_loc + ps%p(n)%d
      pvol = pi * ps%p(n)%d**3 / 6.0_WP
      vf_loc = vf_loc + pvol
      meanVV_loc(1) = meanVV_loc(1) + ps%p(n)%vel(1) * ps%p(n)%vel(1)
      meanVV_loc(2) = meanVV_loc(2) + ps%p(n)%vel(1) * ps%p(n)%vel(2)
      meanVV_loc(3) = meanVV_loc(3) + ps%p(n)%vel(1) * ps%p(n)%vel(3)
      meanVV_loc(4) = meanVV_loc(4) + ps%p(n)%vel(2) * ps%p(n)%vel(2)
      meanVV_loc(5) = meanVV_loc(5) + ps%p(n)%vel(2) * ps%p(n)%vel(3)
      meanVV_loc(6) = meanVV_loc(6) + ps%p(n)%vel(3) * ps%p(n)%vel(3)
      ec%phi_f(i,j,k) = ec%phi_f(i,j,k) + pvol
      ec%vel_x_f(i,j,k) = ec%vel_x_f(i,j,k) + pvol * ps%p(n)%vel(1)
      ec%vel_y_f(i,j,k) = ec%vel_y_f(i,j,k) + pvol * ps%p(n)%vel(2)
      ec%vel_z_f(i,j,k) = ec%vel_z_f(i,j,k) + pvol * ps%p(n)%vel(3)
    end do
    ec%phi_f = ec%phi_f * N3 / ec%sim_pg%vol_total
    ec%vel_x_f = ec%vel_x_f * N3 / ec%sim_pg%vol_total
    ec%vel_y_f = ec%vel_y_f * N3 / ec%sim_pg%vol_total
    ec%vel_z_f = ec%vel_z_f * N3 / ec%sim_pg%vol_total
    call mpi_reduce(ps%np_, n, 1, MPI_INTEGER, MPI_SUM, 0, ec%fft%pg%comm)
    call mpi_reduce(dp_loc, dp, 1, MPI_REAL_WP, MPI_SUM, 0, ec%fft%pg%comm)
    call mpi_reduce(vf_loc, vf, 1, MPI_REAL_WP, MPI_SUM, 0, ec%fft%pg%comm)
    call mpi_reduce(drg_loc, drg, 1, MPI_REAL_WP, MPI_SUM, 0, ec%fft%pg%comm)
    call mpi_reduce(meanVV_loc, meanVV, 6, MPI_REAL_WP, MPI_SUM, 0,           &
      ec%fft%pg%comm)
    if (ec%fft%pg%rank .eq. 0) then
      dp = dp / n; vf = vf / ec%sim_pg%vol_total; meanVV = meanVV / n;
    end if

    ! do forward transforms to get variables in fourier space
    call ec%fft%forward_transform(ec%phi_f)
    call ec%fft%forward_transform(ec%vel_x_f)
    call ec%fft%forward_transform(ec%vel_y_f)
    call ec%fft%forward_transform(ec%vel_z_f)

    ! for each filter, convolve and compute means, write result
    do i = 1, ec%num_filters

      ! set time and step in statpoint
      ec%filters(i)%d%time = time; ec%filters(i)%d%step = step;

      ! compute filtered phi
      ec%phi_r = ec%phi_f * ec%filters(i)%flt_f
      call ec%fft%backward_transform(ec%phi_r)
      ec%phi_r = fft_rescale * ec%phi_r

      ! compute filtered velocity
      ec%vel_x_r = ec%vel_x_f * ec%filters(i)%flt_f
      ec%vel_y_r = ec%vel_y_f * ec%filters(i)%flt_f
      ec%vel_z_r = ec%vel_z_f * ec%filters(i)%flt_f
      call ec%fft%backward_transform(ec%vel_x_r)
      call ec%fft%backward_transform(ec%vel_y_r)
      call ec%fft%backward_transform(ec%vel_z_r)
      ec%vel_x_r = fft_rescale * ec%vel_x_r
      ec%vel_y_r = fft_rescale * ec%vel_y_r
      ec%vel_z_r = fft_rescale * ec%vel_z_r

      ! set S1, S12, S13, S14 (on all tasks)
      tmpvec(1) = sum(realpart(ec%phi_r))
      call mpi_allreduce(tmpvec(1), ec%filters(i)%d%S(1), 1, MPI_REAL_WP,     &
        MPI_SUM, 0, ec%fft%pg%comm)
      ec%filters(i)%d%S(1) = ec%filters(i)%d%S(1) / N3
      tmpvec(:) = (/ sum(realpart(ec%vel_x_r)), sum(realpart(ec%vel_y_r)),    &
        sum(realpart(ec%vel_z_r)) /)
      call mpi_allreduce(tmpvec, ec%filters(i)%d%S(12:14), 3, MPI_REAL_WP,    &
        MPI_SUM, 0, ec%fft%pg%comm)
      ec%filters(i)%d%S(12:14) = ec%filters(i)%d%S(12:14) / N3

      ! set S2, S3, S4 (filtered density velocity product)
      ec%work_x_r = ec%phi_r * ec%vel_x_r
      ec%work_y_r = ec%phi_r * ec%vel_y_r
      ec%work_z_r = ec%phi_r * ec%vel_z_r
      tmpvec(:) = (/ sum(realpart(ec%work_x_r)), sum(realpart(ec%work_y_r)),  &
        sum(realpart(ec%work_z_r)) /)
      call mpi_reduce(tmpvec, ec%filters(i)%d%S(2:4), 3, MPI_REAL_WP,         &
        MPI_SUM, 0, ec%fft%pg%comm)
      ec%filters(i)%d%S(2:4) = ec%filters(i)%d%S(2:4) / N3

      ! done with uncentered quantities; center phi and vel
      ec%phi_r = ec%phi_r - ec%filters(i)%d%S(1)
      ec%vel_x_r = ec%vel_x_r - ec%filters(i)%d%S(12)
      ec%vel_y_r = ec%vel_y_r - ec%filters(i)%d%S(13)
      ec%vel_z_r = ec%vel_z_r - ec%filters(i)%d%S(14)

      ! set S6 through S11
      ec%work_x_r = ec%phi_r * ec%vel_x_r
      ec%work_y_r = ec%phi_r * ec%vel_y_r
      ec%work_z_r = ec%phi_r * ec%vel_z_r
      tmpvec(:) = (/ sum(realpart(ec%work_x_r)), sum(realpart(ec%work_y_r)),  &
        sum(realpart(ec%work_z_r)) /)
      call mpi_reduce(tmpvec, ec%filters(i)%d%S(6:8), 3, MPI_REAL_WP,         &
        MPI_SUM, 0, ec%fft%pg%comm)
      ec%filters(i)%d%S(6:8) = ec%filters(i)%d%S(6:8) / N3
      ec%work_x_r = ec%phi_r * ec%work_x_r
      ec%work_y_r = ec%phi_r * ec%work_y_r
      ec%work_z_r = ec%phi_r * ec%work_z_r
      tmpvec(:) = (/ sum(realpart(ec%work_x_r)), sum(realpart(ec%work_y_r)),  &
        sum(realpart(ec%work_z_r)) /)
      call mpi_reduce(tmpvec, ec%filters(i)%d%S(9:11), 3, MPI_REAL_WP,        &
        MPI_SUM, 0, ec%fft%pg%comm)
      ec%filters(i)%d%S(9:11) = ec%filters(i)%d%S(9:11) / N3

      ! set S5
      ec%phi_r = ec%phi_r * ec%phi_r  ! in place should be faster than work_r
      tmpvec(1) = sum(realpart(ec%phi_r))
      call mpi_reduce(tmpvec(1), ec%filters(i)%d%S(5), 1, MPI_REAL_WP,        &
        MPI_SUM, 0, ec%fft%pg%comm)
      ec%filters(i)%d%S(5) = ec%filters(i)%d%S(5) / N3

      ! copy microscopic statistics
      ec%filters(i)%d%M(1) = dp;           ec%filters(i)%d%M(2) = n;
      ec%filters(i)%d%M(3) = ps%rho;       ec%filters(i)%d%M(4) = vf;
      ec%filters(i)%d%M(5:10) = meanVV(:); ec%filters(i)%d%M(11) = drg;

      !todo B

      !todo H

      ! save stats
      if (ec%fft%pg%rank .eq. 0) call ec%filters(i)%mon%write()

    end do

  end subroutine compute_statistics

end module estimclosures_class

