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
  use hitstats,       only: filterset
  use monitor_class,  only: monitor
  use string,         only: str_medium, str_long
  use sgrid_class,    only: sgrid
  use pgrid_class,    only: pgrid
  use config_class,   only: config
  use lpt_class,      only: lpt
  implicit none
  private

  public :: estimclosures, estimclosures_mesh, ETAODX, RHOF, FORCE_TIMESCALE

  !> mesh spacing to parameter conversion
  real(WP), parameter :: ETAODX = 0.4_WP

  !> parameter primitive quantities
  real(WP), parameter :: RHOF = 1.0_WP
  real(WP), parameter :: LEN_SCALE_PROPORTION = 0.2_WP
  real(WP), parameter :: FORCE_TIMESCALE = 1.0_WP

  !> parameters are stored in arrays, in the following order:
  !> params_nondim - 4 items - Relambda, Stk, phiinf, Wovk
  !> params_primit - 7 items - rhof, rhop, ktarget, epstarget, nu, dp, g

  !> object managing computation of closures for each filter and navigation
  !> of paramter space
  type, abstract :: estimclosures
    ! configs for simulation
    class(pgrid), pointer :: sim_pg
    ! current state
    integer :: step, interval, sweepnum
    real(WP) :: time
    real(WP), dimension(4) :: nondm_actual
    real(WP), dimension(2) :: dimturb    ! urms, eta
    character(len=str_medium) :: out_fname
    type(monitor) :: mon
    ! filterset object
    type(filterset), pointer :: hs
    ! integral length scales
    real(WP) :: linf, etamin, eta
    ! number of timescales before writing output
    real(WP) :: interval_tinfs, sim_burnin_mult, param_burnin_mult
    ! store current parameter set
    real(WP), dimension(7) :: param_target
    real(WP), dimension(4) :: nondm_target
  contains
    procedure :: construct => ec_construct  ! constructors of abstract classes in fortran are wonky
    procedure :: init => ec_init            ! parent initializer
    ! commented because I'm too lazy to write an interface block
    !procedure, deferred :: get_next_params
    !procedure, deferred :: get_interval
    procedure :: compute_statistics
    procedure :: monitor_setup => ec_monitor_setup
    procedure :: io_setup => ec_io_setup
    procedure :: io_write => ec_io_write
    procedure :: destruct => ec_destruct
  end type estimclosures

  ! first entry is datapoint index, 2-5 are nondim params, 6 is sweep
  type, extends(estimclosures) :: estimclosures_mesh
    ! nondim mesh info
    real(WP), dimension(2:5) :: pmin, pmax, pspacing
    logical, dimension(2:5)  :: plog
    integer, dimension(6)    :: Icurr, Imax, Idir
    integer :: interval_number, interval_total
    ! primitive quantities that are constant (but not parameters)
    integer :: Np
  contains
    procedure :: get_next_params
    procedure :: get_interval
    procedure :: params_are_new
  end type estimclosures_mesh
  interface estimclosures_mesh; procedure ecmesh_from_args; end interface;

contains

  subroutine ec_construct(ec, sim_pg)
    use param, only: param_read
    implicit none
    class(estimclosures), intent(inout) :: ec
    class(pgrid), target, intent(in) :: sim_pg

    ! store pointer to simulation config
    ec%sim_pg => sim_pg

    ! set integral length scales
    ec%linf = LEN_SCALE_PROPORTION * ec%sim_pg%vol_total**(1.0_WP / 3)
    !TODO move min and max meshsizes to sgrid or pgrid
    ec%etamin = ETAODX * min(                                                 &
      minval(sim_pg%dx(sim_pg%imin_:sim_pg%imax_)),                           &
      minval(sim_pg%dy(sim_pg%jmin_:sim_pg%jmax_)),                           &
      minval(sim_pg%dz(sim_pg%kmin_:sim_pg%kmax_))                            &
    )

    ! read params
    call param_read('EC integral timescales', ec%interval_tinfs)
    call param_read('EC new sim multiplier', ec%sim_burnin_mult)
    call param_read('EC new params multiplier', ec%param_burnin_mult)

    ! allocate (but don't initialize) hitstats object
    allocate(ec%hs)

  end subroutine ec_construct

  subroutine ec_init(ec, ps, rhof_array, visc, U, V, W, sx, sy, sz)
    use param,       only: param_read
    implicit none
    class(estimclosures), intent(inout) :: ec
    type(lpt), target, intent(inout) :: ps
    real(WP), dimension(:,:,:), target, intent(in) :: rhof_array, visc, U, V, W, sx, sy, sz
    integer, dimension(3) :: FFTN
    character(len=str_medium) :: filterfile

    call param_read('EC filter list', filterfile)
    call param_read('EC FFT mesh', FFTN)
    call ec%hs%init(ec%sim_pg, filterfile, FFTN, rhof_array, visc, U, V, W, sx, sy, sz, ps)

  end subroutine ec_init

  function ecmesh_from_args(cfg) result(ec)
    use messager, only: die
    use param, only: param_read
    implicit none
    type(estimclosures_mesh) :: ec
    class(config), intent(in) :: cfg
    integer :: n
    real(WP) :: Relammax

    ! init parent
    call ec%construct(cfg)

    ! read params
    call param_read('EC min Relambda',       ec%pmin(2))
    call param_read('EC max Relambda',       ec%pmax(2))
    call param_read('EC num Relambda',       ec%Imax(2))
    call param_read('EC log Relambda',       ec%plog(2))
    call param_read('EC min Stk',            ec%pmin(3))
    call param_read('EC max Stk',            ec%pmax(3))
    call param_read('EC num Stk',            ec%Imax(3))
    call param_read('EC log Stk',            ec%plog(3))
    call param_read('EC min vf',             ec%pmin(4))
    call param_read('EC max vf',             ec%pmax(4))
    call param_read('EC num vf',             ec%Imax(4))
    call param_read('EC log vf',             ec%plog(4))
    call param_read('EC min Wovk',           ec%pmin(5))
    call param_read('EC max Wovk',           ec%pmax(5))
    call param_read('EC num Wovk',           ec%Imax(5))
    call param_read('EC log Wovk',           ec%plog(5))
    call param_read('EC data per statpoint', ec%Imax(1))
    call param_read('EC num sweeps',         ec%Imax(6))

    ! make sure the largest requirest Relambda is representable on the mesh
    Relammax = sqrt(15.0_WP) * (ec%linf / ec%etamin)**(2.0_WP / 3)
    if (ec%pmax(2) .gt. Relammax) call die("[EC] requested Relambda is greater&
      & than what is possible to resolve on mesh")

    ! init mesh
    ec%Icurr(:) = 1; ec%Idir(:) = +1;
    do n = 2, 5
      if (ec%plog(n)) then
        if (ec%Imax(n) .gt. 1) then
          ec%pspacing(n) = exp(log(ec%pmax(n) / ec%pmin(n)) / (ec%Imax(n) - 1))
        else
          ec%pspacing(n) = 1.0_WP
        end if
      else
        if (ec%Imax(n) .gt. 1) then
          ec%pspacing(n) = (ec%pmax(n) - ec%pmin(n)) / (ec%Imax(n) - 1)
        else
          ec%pspacing(n) = 0.0_WP
        end if
      end if
    end do
    ec%Icurr(1) = 0
    ec%interval_number = 0
    ec%interval_total = product(ec%Imax)

  end function ecmesh_from_args

  ! assumes out_fname, type, num_params, params are set
  subroutine ec_monitor_setup(ec)
    implicit none
    class(estimclosures), intent(inout) :: ec

    ! setup monitor file
    ec%out_fname = trim(adjustl('ecstate'))
    ec%mon = monitor(ec%sim_pg%amRoot, ec%out_fname)
    call ec%mon%add_column(ec%step,            'step')
    call ec%mon%add_column(ec%time,            'time')
    call ec%mon%add_column(ec%interval,        'interval')
    call ec%mon%add_column(ec%sweepnum,        'sweep')
    call ec%mon%add_column(ec%nondm_target(1), 'tgt_Relam' )
    call ec%mon%add_column(ec%nondm_target(2), 'tgt_Stk'   )
    call ec%mon%add_column(ec%nondm_target(3), 'tgt_phiinf')
    call ec%mon%add_column(ec%nondm_target(4), 'tgt_Wovk'  )
    call ec%mon%add_column(ec%nondm_actual(1), 'act_Relam' )
    call ec%mon%add_column(ec%nondm_actual(2), 'act_Stk'   )
    call ec%mon%add_column(ec%nondm_actual(3), 'act_phiinf')
    call ec%mon%add_column(ec%nondm_actual(4), 'act_Wovk'  )
    call ec%mon%add_column(ec%dimturb(1),      'urms')
    call ec%mon%add_column(ec%dimturb(2),      'eta' )
    call ec%mon%add_column(ec%param_target(5), 'nu')
    call ec%mon%add_column(ec%param_target(7), 'g' )

  end subroutine ec_monitor_setup

  ! params_nondim - 4 items - Relambda, Stk, phiinf, Wovk
  ! params_primit - 7 items - rhof, rhop, ktarget, epstarget, nu, dp, g
  subroutine get_next_params(ec, params, done)
    use mathtools, only: pi
    use messager,  only: die, log
    implicit none
    class(estimclosures_mesh), intent(inout) :: ec
    real(WP), dimension(7), intent(out)      :: params
    logical, intent(out)                     :: done
    integer :: n
    character(len=str_long) :: message

    ! check if we are done
    done = ec%interval_number .eq. ec%interval_total

    ! increment number of total intervals
    ec%interval_number = ec%interval_number + 1

    ! increment first coord
    ec%Icurr(1) = ec%Icurr(1) + ec%Idir(1)

    ! if this pushes it past its bounds, increment the next coordinate instead
    ! and change directions of previous coordinate
    do n = 1, 5
      if (ec%Icurr(n) .lt. 1 .or. ec%Icurr(n) .gt. ec%Imax(n)) then
        ec%Icurr(n) = ec%Icurr(n) - ec%Idir(n)
        ec%Icurr(n+1) = ec%Icurr(n+1) + ec%Idir(n+1)
        ec%Idir(n) = -1 * ec%Idir(n)
      end if
    end do

    ! get nondimensional values
    do n = 2, 5
      if (ec%plog(n)) then
        ec%nondm_target(n-1) = ec%pmin(n) * ec%pspacing(n)**(ec%Icurr(n) - 1)
      else
        ec%nondm_target(n-1) = ec%pmin(n) + (ec%Icurr(n) - 1) * ec%pspacing(n)
      end if
    end do
    ec%sweepnum = ec%Icurr(6)

    ! check to make sure we didn't make a mistake
    if (.not. done) then
      if (any(ec%nondm_target - 1e3_WP * epsilon(1.0_WP) .gt. ec%pmax))             &
        call die('[EC] big whoops')
      if (any(ec%nondm_target + 1e3_WP * epsilon(1.0_WP) .lt. ec%pmin))             &
        call die('[EC] little whoops')
    end if

    ! fluid parameters
    ! using the current approach viscosity is fixed throughout; it is computed
    ! at initialization
    ec%param_target(1) = RHOF
    ec%eta = ec%linf * (sqrt(15.0_WP) / ec%nondm_target(1))**1.5_WP
    if (ec%eta .lt. ec%etamin) call die("[EC] computed eta less than etamin")
    ec%param_target(5) = ec%eta**2 * ec%nondm_target(1) / (FORCE_TIMESCALE * sqrt(15.0_WP))
    ec%param_target(4) = ec%param_target(5)**3 / ec%eta**4
    ec%param_target(3) = 1.5_WP * ec%param_target(4) * FORCE_TIMESCALE

    ! particle and gravity parameters
    ec%param_target(6) = (6 * ec%sim_pg%vol_total * ec%nondm_target(3) / (pi * ec%Np))**(1.0_WP / 3)
    ec%param_target(2) = ec%param_target(1) * 18 * ec%nondm_target(2) * sqrt(ec%param_target(5)**3 / ec%param_target(4)) / ec%param_target(6)**2
    ec%param_target(7) = 18 * ec%param_target(1) / ec%param_target(2) * (ec%param_target(5)**5 * ec%param_target(4) / ec%param_target(6)**8)**0.25_WP * ec%nondm_target(4)

    ! print values if not done
    if (ec%sim_pg%amroot .and. .not. done) then
      write(*,*) "[EC] index: ", ec%Icurr
      write(*,*) "[EC] nondim params: ", ec%nondm_target
      write(*,*) "[EC] dim params: ", ec%param_target
    end if

    ! log change if we moved to new parameters
    if (ec%params_are_new() .and. ec%sim_pg%amRoot) then
      write(message,'("[EC] changed to new param array (dimensional): ",e12.5,&
        &", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5)') &
        ec%param_target
      call log(message)
      write(message,'("[EC] changed to new param array (nondimensional): ",e12&
        &.5,", ",e12.5,", ",e12.5,", ",e12.5)') ec%nondm_target
      call log(message)
    end if

    ! update parent class state
    ec%interval = ec%Icurr(1);
    ec%sweepnum = ec%Icurr(6);

    ! copy to output
    params(:) = ec%param_target(:)

  end subroutine get_next_params

  function params_are_new(ec) result(new)
    implicit none
    class(estimclosures_mesh), intent(in) :: ec
    logical :: new

    new = all((ec%Icurr .eq. 1       .and. ec%Idir .eq. +1) .or.              &
              (ec%Icurr .eq. ec%Imax .and. ec%Idir .eq. -1)      )

  end function params_are_new

  subroutine get_interval(ec, interval)
    implicit none
    class(estimclosures_mesh), intent(in) :: ec
    real(WP), intent(out) :: interval

    interval = ec%interval_tinfs * FORCE_TIMESCALE

    if (ec%params_are_new()) then
      if (all(ec%Icurr(:) .eq. 1)) then             ! new sim
        if (ec%sim_pg%amroot) write(*,*) "[EC] using new sim burnin"
        interval = interval * ec%sim_burnin_mult
      else                                          ! new params
        if (ec%sim_pg%amroot) write(*,*) "[EC] using new params burnin"
        interval = interval * ec%param_burnin_mult
      end if
    end if

  end subroutine get_interval

  subroutine ec_destruct(ec, dealloc_sim_pg)
    implicit none
    class(estimclosures), intent(inout) :: ec
    logical, intent(in), optional :: dealloc_sim_pg

    if (present(dealloc_sim_pg) .and. dealloc_sim_pg) deallocate(ec%sim_pg)

  end subroutine ec_destruct

  subroutine ec_io_setup(ec)
    use npy_class, only: npy
    implicit none
    class(estimclosures), intent(inout) :: ec

    call ec%hs%setup_sliceio()

  end subroutine ec_io_setup

  subroutine ec_io_write(ec, t)
    use npy_class, only: npy
    implicit none
    class(estimclosures), intent(inout) :: ec
    real(WP), intent(in) :: t

    call ec%hs%write_sliceio(t)

  end subroutine ec_io_write

  subroutine compute_statistics(ec, Re_lambda, Stk, phiinf, Wovk, urms, eta, time, step)
    implicit none
    class(estimclosures), intent(inout) :: ec
    real(WP), intent(in)  :: Re_lambda, Stk, phiinf, Wovk, urms, eta, time
    integer, intent(in) :: step

    ec%step = step; ec%time = time;
    ec%nondm_actual(:) = (/ Re_lambda, Stk, phiinf, Wovk /)
    ec%nondm_target(:) = ec%nondm_target(:)
    ec%dimturb(:) = (/ urms, eta /)
    call ec%mon%write()

    call ec%hs%compute_stats(step)

  end subroutine compute_statistics

end module estimclosures_class
