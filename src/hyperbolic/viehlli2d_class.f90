!> Method in Vie 2015 with HLLI Riemann solver
module viehlli2d_class
  use precision,      only: WP
  use iterator_class, only: iterator
  use string,         only: str_medium
  use config_class,   only: config
  implicit none
  private

  ! expose type/constructor/methods
  public :: viehlli2d, viehlli2d_bcond, vh2_presc_ftype,                      &
    viehlli2d_cons2prim, viehlli2d_prim2cons, viehlli2d_isvalid

  interface
    pure subroutine vh2_presc_ftype(idx, pos, U)
      use precision, only: WP
      implicit none
      integer, dimension(3), intent(in) :: idx
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(6), intent(out) :: U
    end subroutine vh2_presc_ftype
  end interface

  ! list of known available bcond types for this solver
  ! here 'right' means bcond direction +1
  ! all bcs are enforced using ghost cells
  ! bits 1 and 2 - 00 nothing, 10 interpolate, 11 reflect
  ! bits 3 and 4 - 00 nothing, 10 zero normal velocity, 01 zero tangential velocity
  ! bit  5       - 1 if a function is provided to prescribe values
  integer(1), parameter, public :: NOTHING = 0_1
  integer(1), parameter, public :: INTERP  = 2_1
  integer(1), parameter, public :: REFLECT = 3_1
  integer(1), parameter, public :: ZERONORMVEL = 4_1
  integer(1), parameter, public :: ZEROTANGVEL = 8_1
  integer(1), parameter, public :: PRESC = 16_1
  ! example bcs
  integer(1), parameter, public :: SLIPWALL = REFLECT
  integer(1), parameter, public :: NOSLIPWALL = IAND(REFLECT, ZEROTANGVEL)
  integer(1), parameter, public :: OUTFLOW = INTERP
  integer(1), parameter, public :: INFLOW = PRESC

  !> boundary conditions
  type :: viehlli2d_bcond
    type(viehlli2d_bcond), pointer :: next                !< linked list of bconds
    character(len=str_medium) :: name = 'UNNAMED_BCOND'   !< bcond name (default UNNAMED_BCOND)
    integer(1) :: bc_type                                 !< bcond type
    type(iterator) :: itr                                 !< this is the iterator for the bcond
    integer(1) :: dir                                     !< bcond direction 0-6
    procedure(vh2_presc_ftype), nopass, pointer :: presc  !< function to get prescribed values
  end type viehlli2d_bcond

  type :: viehlli2d

    ! flow variables and slopes
    real(WP), dimension(:,:,:,:), pointer :: Uc, dU

    ! config
    class(config), pointer :: cfg                         !< config

    ! flag to use (or not use) contact resolution
    logical :: resolve_contacts

    ! boundary condition list
    integer :: nbc                                        !< number of bcs
    type(viehlli2d_bcond), pointer :: first_bc            !< bc list

    ! cfl numbers
    real(WP) :: cfl_x, cfl_y                              !< cfl numbers

    ! Monitoring quantities
    real(WP), dimension(6) :: Umin, Umax, Uint

  contains

    procedure :: print => viehlli2d_print                 !< print solver info
    procedure :: add_bcond => viehlli2d_add_bcond         !< add a boundary condition
    procedure :: get_bcond => viehlli2d_get_bcond         !< get a boundary condition
    procedure :: apply_bcond => viehlli2d_apply_bcs       !< apply all boundary conditions
    procedure :: get_cfl => viehlli2d_get_cfl             !< calculate maximum cfl
    procedure :: get_max => viehlli2d_get_range           !< calculate field value range
    procedure :: compute_dU_x => viehlli2d_comp_dU_x      !< take a step in x
    procedure :: compute_dU_y => viehlli2d_comp_dU_y      !< take a step in y
    procedure, nopass :: rhs => viehlli2d_rhs             !< get rhs
    procedure, nopass :: backeuler => viehlli2d_backeuler !< get back euler difference

  end type viehlli2d

  interface viehlli2d; procedure viehlli2d_from_args; end interface;

contains

  ! pure
  subroutine viehlli2d_cons2prim(cons, prim)
    implicit none
    real(WP), intent(in), dimension(6) :: cons
    real(WP), intent(out), dimension(6) :: prim

    prim(1) = cons(1)
    prim(2) = cons(2) / cons(1)
    prim(3) = cons(3) / cons(1)
    prim(4) = cons(4) / cons(1) - prim(2)**2
    prim(5) = cons(5) / cons(1) - prim(2)*prim(3)
    prim(6) = cons(6) / cons(1) - prim(3)**2

  end subroutine viehlli2d_cons2prim

  pure subroutine viehlli2d_prim2cons(prim, cons)
    implicit none
    real(WP), intent(in), dimension(6) :: prim
    real(WP), intent(out), dimension(6) :: cons

    cons(1) = prim(1)
    cons(2) = prim(1) * prim(2)
    cons(3) = prim(1) * prim(3)
    cons(4) = prim(1) * (prim(4) + prim(2)**2)
    cons(5) = prim(1) * (prim(5) + prim(2)*prim(3))
    cons(6) = prim(1) * (prim(6) + prim(3)**2)

  end subroutine viehlli2d_prim2cons

  ! input conservative
  !pure
  subroutine viehlli2d_maxwavespeeds(U, sx, sy)
    implicit none
    real(WP), intent(in), dimension(6) :: U
    real(WP), intent(out) :: sx, sy
    real(WP) :: vx, vy

    vx = U(2) / U(1); vy = U(3) / U(1);
    if (U(4) / U(1) - vx**2 .lt. 0.0_WP) then
      write(*,*) 'U(4) / U(1) - vx**2 = ', U(4) / U(1) - vx**2
    end if
    if (U(6) / U(1) - vy**2 .lt. 0.0_WP) then
      write(*,*) 'U(6) / U(1) - vy**2 = ', U(6) / U(1) - vy**2
    end if
    sx = abs(vx) + sqrt(3 * abs(U(4) / U(1) - vx**2))
    sy = abs(vy) + sqrt(3 * abs(U(6) / U(1) - vy**2))

  end subroutine viehlli2d_maxwavespeeds

  ! input primitive, output conserved
  pure subroutine viehlli2d_flux(Up, F)
    implicit none
    real(WP), intent(in), dimension(6) :: Up
    real(WP), intent(out), dimension(6) :: F

    F(1) = Up(1) * Up(2)
    F(2) = Up(1) * (Up(2)**2      + Up(4))
    F(3) = Up(1) * (Up(2) * Up(3) + Up(5))
    F(4) = Up(1) * (Up(2)**3 + 3 * Up(2) * Up(4))
    F(5) = Up(1) * (Up(2)**2 * Up(3) + 2 * Up(2) * Up(5) +     Up(3) * Up(4))
    F(6) = Up(1) * (Up(2) * Up(3)**2 +     Up(2) * Up(6) + 2 * Up(3) * Up(5))

  end subroutine viehlli2d_flux

  ! Ubp in primitive variables, Udiff in conservative variables, contacts in
  ! conservative variables
  pure subroutine viehlli2d_contactwaves(Ubp, Udiff, sL, sR, contacts)
    implicit none
    real(WP), dimension(6), intent(in) :: Ubp, Udiff
    real(WP), intent(in) :: sL, sR
    real(WP), dimension(6), intent(out) :: contacts
    real(WP), dimension(6,6) :: RdL
    real(WP) :: v1, v2, Pp11, Pp12, Pp22, d1, d2, d3, d4, c, phi
    real(WP), dimension(4) :: d, peigs, neigs
    real(WP), parameter :: PHI_EPS = 1.0e-12_WP

    v1 = Ubp(2); v2 = Ubp(3); Pp11 = Ubp(4); Pp12 = Ubp(5); Pp22 = Ubp(6);
    c = sqrt(Ubp(4))

    peigs = (/ v1 - c, v1, v1, v1 + c /)
    neigs = min(peigs, 0.0_WP)
    peigs = peigs - neigs
    d = 1.0_WP - neigs / sL - peigs / sR
    d1 = d(1); d2 = d(2); d3 = d(3); d4 = d(4);

    RdL(1, :) = (/ d2 - d2*v1**2/(3*c**2), 2*d2*v1/(3*c**2), 0.0_WP, -d2/(3*c**2), 0.0_WP, 0.0_WP /)
    RdL(2, :) = (/ d2*v1 - d2*v1**3/(3*c**2), 2*d2*v1**2/(3*c**2), 0.0_WP, -d2*v1/(3*c**2), 0.0_WP, 0.0_WP /)
    RdL(3, 1) = (Pp12*c*d1*v1/2 + Pp12*c*d4*v1/2 + Pp12*d1*v1**2/2 - Pp12*d4*v1**2/2 - c**3*d1*v2/2 + c**3*d2*v2 - c**3*d4*v2/2 - c**2*d1*v1*v2/2 + c**2*d4*v1*v2/2 - c*d1*v1**2*v2 - c*d2*v1**2*v2/3 - c*d4*v1**2*v2 - d1*v1**3*v2 + d4*v1**3*v2)/c**3
    RdL(3, 2) = (-Pp12*c*d1/2 - Pp12*c*d4/2 - Pp12*d1*v1 + Pp12*d4*v1 + c**2*d1*v2/2 - c**2*d4*v2/2 + c*d1*v1*v2 + 2*c*d2*v1*v2/3 + c*d4*v1*v2 + 2*d1*v1**2*v2 - 2*d4*v1**2*v2)/c**3
    RdL(3, 3) = (c*d1 + c*d4 + d1*v1 - d4*v1)/(2*c)
    RdL(3, 4) = (Pp12*d1/2 - Pp12*d4/2 - c*d2*v2/3 - d1*v1*v2 + d4*v1*v2)/c**3
    RdL(3, 5) = (-d1 + d4)/(2*c)
    RdL(3, 6) = 0.0_WP
    RdL(4, :) = (/ d2*v1**2 - d2*v1**4/(3*c**2), 2*d2*v1**3/(3*c**2), 0.0_WP, -d2*v1**2/(3*c**2), 0.0_WP, 0.0_WP /)
    RdL(5, 1) = -Pp12*d1*v1/(2*c) + Pp12*d4*v1/(2*c) + Pp12*d1*v1**3/(2*c**3) - Pp12*d4*v1**3/(2*c**3) + c*d1*v2/2 - c*d4*v2/2 + d2*v1*v2 + d1*v1**2*v2/(2*c) - d4*v1**2*v2/(2*c) - d2*v1**3*v2/(3*c**2) - d1*v1**4*v2/c**3 + d4*v1**4*v2/c**3
    RdL(5, 2) = (Pp12*c**2*d1/2 - Pp12*c**2*d4/2 + Pp12*c*d1*v1/2 + Pp12*c*d4*v1/2 - Pp12*d1*v1**2 + Pp12*d4*v1**2 - c**3*d1*v2/2 - c**3*d4*v2/2 - c**2*d1*v1*v2/2 + c**2*d4*v1*v2/2 - c*d1*v1**2*v2 + 2*c*d2*v1**2*v2/3 - c*d4*v1**2*v2 + 2*d1*v1**3*v2 - 2*d4*v1**3*v2)/c**3
    RdL(5, 3) = (-c**2*d1 + c**2*d4 + d1*v1**2 - d4*v1**2)/(2*c)
    RdL(5, 4) = (-Pp12*c*d1/2 - Pp12*c*d4/2 + Pp12*d1*v1/2 - Pp12*d4*v1/2 + c*d1*v1*v2 - c*d2*v1*v2/3 + c*d4*v1*v2 - d1*v1**2*v2 + d4*v1**2*v2)/c**3
    RdL(5, 5) = (c*d1 + c*d4 - d1*v1 + d4*v1)/(2*c)
    RdL(5, 6) = 0.0_WP
    RdL(6, 1) = (-Pp12**2*c*d1*v1 + Pp12**2*c*d4*v1 - Pp12**2*d1*v1**2 + 4*Pp12**2*d3*v1**2/3 - Pp12**2*d4*v1**2 + Pp12*c**3*d1*v2 - Pp12*c**3*d4*v2 + 2*Pp12*c**2*d1*v1*v2 - 2*Pp12*c**2*d3*v1*v2 + 2*Pp12*c**2*d4*v1*v2 + 5*Pp12*c*d1*v1**2*v2 - 5*Pp12*c*d4*v1**2*v2 + 4*Pp12*d1*v1**3*v2 - 16*Pp12*d3*v1**3*v2/3 + 4*Pp12*d4*v1**3*v2 - Pp22*c**2*d3*v1**2/3 - c**4*d1*v2**2 + 2*c**4*d3*v2**2 - c**4*d4*v2**2 - 3*c**3*d1*v1*v2**2 + 3*c**3*d4*v1*v2**2 - 4*c**2*d1*v1**2*v2**2 + 13*c**2*d3*v1**2*v2**2/3 - 4*c**2*d4*v1**2*v2**2 - 6*c*d1*v1**3*v2**2 + 6*c*d4*v1**3*v2**2 - 4*d1*v1**4*v2**2 + 16*d3*v1**4*v2**2/3 - 4*d4*v1**4*v2**2)/c**4
    RdL(6, 2) = (3*Pp12**2*c*d1 - 3*Pp12**2*c*d4 + 6*Pp12**2*d1*v1 - 8*Pp12**2*d3*v1 + 6*Pp12**2*d4*v1 - 6*Pp12*c**2*d1*v2 + 6*Pp12*c**2*d3*v2 - 6*Pp12*c**2*d4*v2 - 18*Pp12*c*d1*v1*v2 + 18*Pp12*c*d4*v1*v2 - 24*Pp12*d1*v1**2*v2 + 32*Pp12*d3*v1**2*v2 - 24*Pp12*d4*v1**2*v2 + 2*Pp22*c**2*d3*v1 + 3*c**3*d1*v2**2 - 3*c**3*d4*v2**2 + 12*c**2*d1*v1*v2**2 - 14*c**2*d3*v1*v2**2 + 12*c**2*d4*v1*v2**2 + 24*c*d1*v1**2*v2**2 - 24*c*d4*v1**2*v2**2 + 24*d1*v1**3*v2**2 - 32*d3*v1**3*v2**2 + 24*d4*v1**3*v2**2)/(3*c**4)
    RdL(6, 3) = (-Pp12*c*d1 + Pp12*c*d4 - Pp12*d1*v1 + 2*Pp12*d3*v1 - Pp12*d4*v1 + c**2*d1*v2 - 2*c**2*d3*v2 + c**2*d4*v2 + 3*c*d1*v1*v2 - 3*c*d4*v1*v2 + 2*d1*v1**2*v2 - 4*d3*v1**2*v2 + 2*d4*v1**2*v2)/c**2
    RdL(6, 4) = (-3*Pp12**2*d1 + 4*Pp12**2*d3 - 3*Pp12**2*d4 + 3*Pp12*c*d1*v2 - 3*Pp12*c*d4*v2 + 12*Pp12*d1*v1*v2 - 16*Pp12*d3*v1*v2 + 12*Pp12*d4*v1*v2 - Pp22*c**2*d3 + c**2*d3*v2**2 - 6*c*d1*v1*v2**2 + 6*c*d4*v1*v2**2 - 12*d1*v1**2*v2**2 + 16*d3*v1**2*v2**2 - 12*d4*v1**2*v2**2)/(3*c**4)
    RdL(6, 5) = (Pp12*d1 - 2*Pp12*d3 + Pp12*d4 - c*d1*v2 + c*d4*v2 - 2*d1*v1*v2 + 4*d3*v1*v2 - 2*d4*v1*v2)/c**2
    RdL(6, 6) = d3

    !TODO phi
    phi = 1.0_WP

    contacts(:) = matmul(RdL, Udiff)
    contacts(:) = (-phi * sR * sL / (sR - sL)) * contacts(:)

  end subroutine viehlli2d_contactwaves

  ! input in primitive variables, flux in conservative variables
  !pure
  subroutine viehlli2d_cellflux(resolve_contacts, Ubp, Ulp, Urp, F)
    implicit none
    logical, intent(in) :: resolve_contacts
    real(WP), intent(in), dimension(6) :: Ubp, Ulp, Urp
    real(WP), intent(out), dimension(6) :: F
    real(WP), dimension(6) :: Ul, Ur, Ub, Udiff, Fl, Fr
    real(WP) :: sL, sR

    sL = min(Ubp(2) - sqrt(3 * Ubp(4)), Ulp(2) - sqrt(3 * Ulp(4)), 0.0_WP)
    sR = max(Ubp(2) + sqrt(3 * Ubp(4)), Urp(2) + sqrt(3 * Urp(4)), 0.0_WP)

    call viehlli2d_flux(Ulp, Fl)
    call viehlli2d_flux(Urp, Fr)
    call viehlli2d_prim2cons(Ulp, Ul)
    call viehlli2d_prim2cons(Urp, Ur)
    Udiff(:) = Ur(:) - Ul(:)

    F(:) = (sR * Fl(:) - sL * Fr(:)) + sL * sR * Udiff(:)
    F(:) = (sR - sL)**(-1) * F(:)
    
    !TODO can we break this correction out of this loop?
    if (resolve_contacts) then
      call viehlli2d_prim2cons(Ubp, Ub)
      call viehlli2d_contactwaves(Ub, Udiff, sL, sR, Fl)
      F = F + Fl
    end if

  end subroutine viehlli2d_cellflux

  subroutine viehlli2d_print(this)
    implicit none
    class(viehlli2d), intent(in) :: this

    write(*,*) 'Two-Fluid solver in the Vie 2015 paper with HLLI Riemann solver in 2d'

  end subroutine viehlli2d_print

  function viehlli2d_from_args(cfg, resolve_contacts) result(this)
    implicit none
    type(viehlli2d) :: this
    class(config), target, intent(in) :: cfg
    logical, intent(in) :: resolve_contacts

    !< config
    this%cfg => cfg

    !< whether to resolve contacts
    this%resolve_contacts = resolve_contacts

    !< allocate arrays
    allocate(this%Uc(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    allocate(this%dU(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

    !< nullify boundary condition list
    this%nbc = 0; nullify(this%first_bc);

  end function viehlli2d_from_args

  !> check if state is valid, useful for debugging, input conserved coords
  !TODO pure
  function viehlli2d_isvalid(U) result(v)
    implicit none
    real(WP), dimension(6), intent(in) :: U
    real(WP), dimension(6) :: prim
    logical :: v

    call viehlli2d_cons2prim(U, prim)

    v = prim(1) .gt. 0.0_WP .and. prim(4) .ge. 0.0_WP .and.                   &
      prim(6) .ge. 0.0_WP .and. prim(6) * prim(4) - prim(5)**2 .ge. 0.0_WP

  end function viehlli2d_isvalid

  !> source terms
  pure subroutine viehlli2d_rhs(gvec, taup, fl_vel, U, rhs)
    implicit none
    real(WP), intent(in) :: taup
    real(WP), dimension(2), intent(in) :: gvec, fl_vel
    real(WP), dimension(6), intent(in) :: U
    real(WP), dimension(6), intent(out) :: rhs
    real(WP), dimension(2) :: pt_vel

    pt_vel(:) = U(2:3) / U(1)

    rhs(1) = 0.0_WP

    rhs(2:3) = U(1) * gvec + U(1) / taup * (fl_vel - pt_vel)

    rhs(4) = pt_vel(1) * gvec(1) + (pt_vel(1) * fl_vel(1) - U(4) / U(1)) / taup
    rhs(4) = 2 * U(1) * rhs(4)

    rhs(5) = gvec(1) * pt_vel(2) + gvec(2) * pt_vel(1)
    rhs(5) = rhs(5) + (pt_vel(1) * fl_vel(2) + pt_vel(2) * fl_vel(1) -        &
      & 2 * U(5) / U(1)) / taup
    rhs(5) = U(1) * rhs(5)

    rhs(6) = gvec(2) * pt_vel(2) + (pt_vel(2) * fl_vel(2) - U(6) / U(1)) / taup
    rhs(6) = 2 * U(1) * rhs(6)

  end subroutine viehlli2d_rhs

  pure subroutine viehlli2d_backeuler(gvec, taup, fl_vel, dt, U, dU)
    implicit none
    real(WP), intent(in) :: taup, dt
    real(WP), dimension(2), intent(in) :: gvec, fl_vel
    real(WP), dimension(6), intent(in) :: U
    real(WP), dimension(6), intent(out) :: dU

    dU(1) = 0.0_WP
    dU(2) = (gvec(1)*U(1)*taup + U(1)*fl_vel(1) - U(2))/(dt + taup)
    dU(3) = (gvec(2)*U(1)*taup + U(1)*fl_vel(2) - U(3))/(dt + taup)
    dU(4) = 2*(-dt*U(4) + dt*gvec(1)**2*U(1)*taup**2 + 2*dt*gvec(1)*U(1)*taup*fl_vel(1) + dt*U(1)*fl_vel(1)**2 - U(4)*taup + gvec(1)*U(2)*taup**2 + U(2)*taup*fl_vel(1))
    dU(5) = -2*dt*U(5) + 2*dt*gvec(1)*gvec(2)*U(1)*taup**2 + 2*dt*gvec(1)*U(1)*taup*fl_vel(2) + 2*dt*gvec(2)*U(1)*taup*fl_vel(1) + 2*dt*U(1)*fl_vel(1)*fl_vel(2) - 2*U(5)*taup + gvec(1)*U(3)*taup**2 + gvec(2)*U(2)*taup**2 + U(2)*taup*fl_vel(2) + U(3)*taup*fl_vel(1)
    dU(6) = 2*(-dt*U(6) + dt*gvec(2)**2*U(1)*taup**2 + 2*dt*gvec(2)*U(1)*taup*fl_vel(2) + dt*U(1)*fl_vel(2)**2 - U(6)*taup + gvec(2)*U(3)*taup**2 + U(3)*taup*fl_vel(2))
    dU(4:6) = dU(4:6) / (2*dt**2 + 3*dt*taup + taup**2)

  end subroutine viehlli2d_backeuler

  ! input and output in primitive variables
  pure subroutine viehlli2d_applyslopes(Ump, Ds, dx, d, Us)
    implicit none
    real(WP), intent(in), dimension(6) :: Ump
    real(WP), intent(in), dimension(6) :: Ds
    real(WP), intent(in) :: dx
    integer, intent(in) :: d
    real(WP), intent(out), dimension(6) :: Us
    real(WP) :: alpha

    alpha = 1.0_WP + Ds(1)**2 * dx**2 / (12 * Ump(1)**2)

    Us(1) = Ump(1)
    Us(2) = Ump(2) + d * Ds(1) * Ds(2) / Ump(1) * (12.0_WP * dx**2)
    Us(3) = Ump(3) + d * Ds(1) * Ds(3) / Ump(1) * (12.0_WP * dx**2)
    Us(4) = Ump(4) + d * 12.0_WP * dx**2 * (alpha * Ds(2)**2 + Ds(1) / Ump(1) * Ds(4))
    Us(5) = Ump(5) + d * 12.0_WP * dx**2 * (alpha * Ds(2) * Ds(3) + Ds(1) / Ump(1) * Ds(5))
    Us(6) = Ump(6) + d * 12.0_WP * dx**2 * (alpha * Ds(3)**2 + Ds(1) / Ump(1) * Ds(6))

  end subroutine viehlli2d_applyslopes

  ! input primitive form, slopes in primitive form
  !pure
  subroutine viehlli2d_compute_slopes(dxi, Ulp, Ump, Urp, Ds)
    implicit none
    real(WP), intent(in) :: dxi
    real(WP), intent(in), dimension(6) :: Ulp, Ump, Urp
    real(WP), intent(out), dimension(6) :: Ds
    real(WP), dimension(6) :: Usp
    real(WP) :: Dmax, xt, ytl, ytr, A, B, alpha, H, corr, dxiadj, Dstar, disc
    integer :: d
    real(WP), parameter :: B_EPS = 1.0E-12_WP

    !TODO DEBUG
    Ds(:) = 0.0_WP
    return

    Ds(1) = 0.5_WP * (sign(1.0_WP, Urp(1) - Ump(1)) + sign(1.0_WP, Ump(1) - Ulp(1)))
    Ds(1) = Ds(1) * min(Urp(1) - Ump(1), Ump(1) - Ulp(1), 2 * Ump(1)) * dxi
    dxiadj = dxi / (1.0_WP - Ds(1) / (6 * Ump(1) * dxi))
    alpha = 1.0_WP + Ds(1)**2 / (12 * Ump(1)**2 * dxi**2)
    Dmax = sqrt(max(12 * Ump(4) * dxi**2 / alpha, 0.0_WP))
    Ds(2) = 0.5_WP * (sign(1.0_WP, Urp(2) - Ump(2)) + sign(1.0_WP, Ump(2) - Ulp(2)))
    Ds(2) = Ds(2) * min(abs(Urp(2) - Ump(2)), abs(Ump(2) - Ulp(2)), Dmax) * dxiadj
    Dmax = sqrt(max(12 * Ump(6) * dxi**2 / alpha, 0.0_WP))
    Ds(3) = 0.5_WP * (sign(1.0_WP, Urp(3) - Ump(3)) + sign(1.0_WP, Ump(3) - Ulp(3)))
    Ds(3) = Ds(3) * min(abs(Urp(3) - Ump(3)), abs(Ump(3) - Ulp(3)), Dmax) * dxiadj
    H = Ump(4) * Ds(3)**2 + Ump(6) * Ds(2)**2 - 2 * Ump(5) * Ds(2) * Ds(3)
    if (H .le. 0.0_WP) then
      corr = 1.0_WP
    else
      if ((Ump(4) * Ump(6) - Ump(5)**2) * (12 * dxi**2) / (alpha * H) .lt. 0.0_WP) write(*,*) (Ump(4) * Ump(6) - Ump(5)**2) * (12 * dxi**2) / (alpha * H)
      corr = min(1.0_WP, sqrt((Ump(4) * Ump(6) - Ump(5)**2) * (12 * dxi**2) / (alpha * H)))
    end if
    Ds(2:3) = corr * Ds(2:3)
    Ds(4:6) = 0.5_WP * (sign(1.0_WP, Urp(4:6) - Ump(4:6)) + sign(1.0_WP, Ump(4:6) - Ulp(4:6)))
    Ds(4:6) = Ds(4:6) * min(abs(Urp(4:6) - Ump(4:6)), abs(Ump(4:6) - Ulp(4:6))) * dxiadj
    ! set values const across two cell edges and initial corr value
    B = Ds(4) * Ds(6) - Ds(5)**2
    corr = 1.0_WP
    ! check left edge
    call viehlli2d_applyslopes(Ump, Ds, 1.0_WP / dxi, -1, Usp)
    A = Usp(4) * Ds(6) + Usp(6) * Ds(4) - 2 * Usp(5) * Ds(5)
    xt = d * 0.5_WP / dxi - Ds(1) / (12 * dxi**2 * Ump(1))
    Dstar = Usp(4) * Usp(6) - Usp(5)**2
    disc = A**2 - 4 * B * Dstar
    if (abs(B) .gt. B_EPS) then                           ! not linear
      if (disc .gt. 0) then                               ! poly has roots
        ytr = sqrt(disc)
        ytl = (- A - ytr) / (2 * B * xt)
        ytr = (- A + ytr) / (2 * B * xt)
        if (B * xt**2 .gt. 0.0_WP) then
          if (ytl .le. corr .and. corr .le. ytr) then
            corr = ytl
          end if
        else
          corr = min(corr, ytr)
        end if
      end if
    end if
    ! check right edge
    call viehlli2d_applyslopes(Ump, Ds, 1.0_WP / dxi, +1, Usp)
    A = Usp(4) * Ds(6) + Usp(6) * Ds(4) - 2 * Usp(5) * Ds(5)
    xt = d * 0.5_WP / dxi + Ds(1) / (12 * dxi**2 * Ump(1))
    Dstar = Usp(4) * Usp(6) - Usp(5)**2
    disc = A**2 - 4 * B * Dstar
    if (abs(B) .gt. B_EPS) then                           ! not linear
      if (disc .gt. 0) then                               ! poly has roots
        ytr = sqrt(A**2 - 4 * B * Dstar)
        ytl = (- A - ytr) / (2 * B * xt)
        ytr = (- A + ytr) / (2 * B * xt)
        if (B * xt**2 .gt. 0.0_WP) then
          if (ytl .le. corr .and. corr .le. ytr) then
            corr = ytl
          end if
        else
          corr = min(corr, ytr)
        end if
      end if
    end if
    ! adjust by correction factor
    Ds(4:6) = corr * Ds(4:6)

  end subroutine viehlli2d_compute_slopes

  ! within this function we use Vie's notation
  !pure
  subroutine viehlli2d_step_core(imino_, imaxo_, jmino_, jmaxo_, kmino_, &
      kmaxo_, dxi, wU, wdU, dt, resolve_contacts)
    implicit none
    integer, intent(in) :: imino_, imaxo_, jmino_, jmaxo_, kmino_, kmaxo_
    real(WP), intent(in), dimension(imino_:imaxo_) :: dxi
    real(WP), dimension(6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(in) :: wU
    real(WP), dimension(6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(out) :: wdU
    real(WP), intent(in) :: dt
    logical, intent(in) :: resolve_contacts
    real(WP), dimension(6) :: Ulp, Ump, Urp, Uls, Urs, Ds, F
    integer :: i, j, k
    real(WP), parameter :: DET_EPS = 1e-12_WP

    !< compute fluxes
    do k = kmino_, kmaxo_
      do j = jmino_, jmaxo_
        call viehlli2d_cons2prim(wU(:,imino_  ,j,k), Ulp)
        call viehlli2d_cons2prim(wU(:,imino_+1,j,k), Ump)
        call viehlli2d_cons2prim(wU(:,imino_+2,j,k), Urp)
        call viehlli2d_compute_slopes(dxi(imino_+1), Ulp, Ump, Urp, Ds)
        call viehlli2d_applyslopes(Ump, Ds, 1.0_WP / dxi(imino_+1), +1, Uls)
        Ulp(:) = Ump(:); Ump(:) = Urp(:);
        do i = imino_+2, imaxo_-1
          call viehlli2d_cons2prim(wU(:,i+1,j,k), Urp)
          call viehlli2d_compute_slopes(dxi(i), Ulp, Ump, Urp, Ds)
          call viehlli2d_applyslopes(Ump, Ds, 1.0_WP / dxi(i), -1, Urs)
          ! flux from left boundary
          !call viehlli2d_cellflux(resolve_contacts, 0.5_WP * (Uls + Urs), Uls, Urs, F)
          call viehlli2d_cellflux(resolve_contacts, 0.5_WP * (Ulp + Ump), Ulp, Ump, F)
          wdU(:,i-1,j,k) = wdU(:,i-1,j,k) - F
          wdU(:,i  ,j,k) = F
          ! prep for next iteration
          !call viehlli2d_applyslopes(Ump, Ds, 1.0_WP / dxi(i), +1, Uls)
          Ulp(:) = Ump(:); Ump(:) = Urp(:);
        end do
      end do
    end do

    !< rescale fluxes into step
    wdU(:,imino_:imino_+1,:,:) = 0.0_WP; wdU(:,imaxo_-1:imaxo_,:,:) = 0.0_WP;
    do i = imino_+2, imaxo_-2
      wdU(:,i,:,:) = dt * dxi(i) * wdU(:,i,:,:)
    end do

  end subroutine viehlli2d_step_core

  subroutine viehlli2d_add_bcond(this, name, bc_type, locator, dir, prescfun)
    use iterator_class, only: locator_ftype
    use string,         only: lowercase
    use messager,       only: die
    implicit none
    class(viehlli2d), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(1), intent(in) :: bc_type
    procedure(locator_ftype) :: locator
    character(len=2), intent(in) :: dir
    type(viehlli2d_bcond), pointer :: new_bc
    procedure(vh2_presc_ftype), pointer, intent(in), optional :: prescfun

    ! prepare new bcond
    allocate(new_bc)
    new_bc%name = trim(adjustl(name))
    new_bc%bc_type = bc_type
    new_bc%itr = iterator(pg=this%cfg, name=new_bc%name, locator=locator)
    if (present(prescfun)) then
      new_bc%presc => prescfun
    else
      nullify(new_bc%presc)
    end if
    select case (lowercase(dir))
    case ('c');              new_bc%dir = 0_1
    case ('-x', 'x-', 'xl'); new_bc%dir = 2_1
    case ('+x', 'x+', 'xr'); new_bc%dir = 1_1
    case ('-y', 'y-', 'yl'); new_bc%dir = 4_1
    case ('+y', 'y+', 'yr'); new_bc%dir = 3_1
    case ('-z', 'z-', 'zl'); new_bc%dir = 6_1
    case ('+z', 'z+', 'zr'); new_bc%dir = 5_1
    case default; call die('[viehlli2d add_bcond] unknown bcond direction')
    end select

    ! insert it in front
    new_bc%next => this%first_bc
    this%first_bc => new_bc
    this%nbc = this%nbc + 1

  end subroutine viehlli2d_add_bcond

  !> get a boundary condition
  !  instead of calling die, we just return a null pointer if bc is not found
  subroutine viehlli2d_get_bcond(this, name, my_bc)
    implicit none
    class(viehlli2d), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(viehlli2d_bcond), pointer, intent(out) :: my_bc

    my_bc => this%first_bc
    do while (associated(my_bc))
      if (trim(my_bc%name).eq.trim(name)) return
      my_bc => my_bc%next
    end do

  end subroutine viehlli2d_get_bcond

  !> apply boundary conditions
  subroutine viehlli2d_apply_bcs(this)
    implicit none
    class(viehlli2d), intent(inout) :: this

    !TODO

  end subroutine viehlli2d_apply_bcs

  !> get the cfl number
  subroutine viehlli2d_get_cfl(this, dt, cfl)
    use mpi_f08, only: MPI_MAX, mpi_allreduce
    use parallel, only: MPI_REAL_WP
    implicit none
    class(viehlli2d), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP), intent(out) :: cfl
    real(WP) :: cfl_x_, cfl_y_, sx, sy
    integer :: i, j, k, ierr

    cfl_x_ = 0.0_WP; cfl_y_ = 0.0_WP;
    do k = this%cfg%kmin_, this%cfg%kmax_
      do j = this%cfg%jmin_, this%cfg%jmax_
        do i = this%cfg%imin_, this%cfg%imax_
          call viehlli2d_maxwavespeeds(this%Uc(:,i,j,k), sx, sy)
          cfl_x_ = max(cfl_x_, sx * this%cfg%dxmi(i))
          cfl_y_ = max(cfl_y_, sy * this%cfg%dymi(j))
        end do
      end do
    end do
    call mpi_allreduce(cfl_x_, this%cfl_x, 1, MPI_REAL_WP, MPI_MAX, this%cfg%comm, ierr)
    call mpi_allreduce(cfl_y_, this%cfl_y, 1, MPI_REAL_WP, MPI_MAX, this%cfg%comm, ierr)

    cfl = dt * max(this%cfl_x, this%cfl_y)

  end subroutine viehlli2d_get_cfl

  !> get the range of the solution
  subroutine viehlli2d_get_range(this)
    use mpi_f08, only: MPI_MIN, MPI_MAX, MPI_SUM, mpi_allreduce
    use parallel, only: MPI_REAL_WP
    implicit none
    class(viehlli2d), intent(inout) :: this
    integer :: i, j, k, n, ierr
    real(WP), dimension(6) :: Umin_, Umax_, Uint_

    Umin_(:) = +huge(this%Umin)
    Umax_(:) = -huge(this%Umax)
    Uint_(:) = 0.0_WP
    do k = this%cfg%kmin_, this%cfg%kmax_
      do j = this%cfg%jmin_, this%cfg%jmax_
        do i = this%cfg%imin_, this%cfg%imax_
          do n = 1, 6
            Umin_(n) = min(Umin_(n), this%Uc(n,i,j,k))
            Umax_(n) = max(Umax_(n), this%Uc(n,i,j,k))
            Uint_(n) = Uint_(n) + this%Uc(n,i,j,k)
          end do
        end do
      end do
    end do
    call mpi_allreduce(Umin_, this%Umin, 6, MPI_REAL_WP, MPI_MIN, this%cfg%comm, ierr)
    call mpi_allreduce(Umax_, this%Umax, 6, MPI_REAL_WP, MPI_MAX, this%cfg%comm, ierr)
    call mpi_allreduce(Uint_, this%Uint, 6, MPI_REAL_WP, MPI_SUM, this%cfg%comm, ierr)

  end subroutine viehlli2d_get_range

  !> step in x direction
  subroutine viehlli2d_comp_dU_x(this, dt)
    implicit none
    class(viehlli2d), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP), dimension(:,:,:,:), allocatable :: wdU

    allocate(wdU(6,this%cfg%imino_:this%cfg%imaxo_,                      &
      this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

    call viehlli2d_step_core(                                                 &
      this%cfg%imino_, this%cfg%imaxo_, this%cfg%jmino_, this%cfg%jmaxo_,     &
      this%cfg%kmino_, this%cfg%kmaxo_, this%cfg%dxmi, this%Uc, wdU, dt,      &
      this%resolve_contacts                                                   &
    )

    this%dU(:,:,:,:) = this%dU(:,:,:,:) + wdU(:,:,:,:)

    deallocate(wdU)

  end subroutine viehlli2d_comp_dU_x

  !> step in y direction
  subroutine viehlli2d_comp_dU_y(this, dt)
    implicit none
    class(viehlli2d), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP), dimension(:,:,:,:), allocatable :: wU, wdU
    integer :: i, j, k

    allocate(wU(6,this%cfg%jmino_:this%cfg%jmaxo_,                       &
      this%cfg%imino_:this%cfg%imaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(wdU(6,this%cfg%jmino_:this%cfg%jmaxo_,                      &
      this%cfg%imino_:this%cfg%imaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

    do k = this%cfg%kmino_, this%cfg%kmaxo_
      do j = this%cfg%jmino_, this%cfg%jmaxo_
        do i = this%cfg%imino_, this%cfg%imaxo_
          wU(1,j,i,k) = this%Uc(1,i,j,k)
          wU(2,j,i,k) = this%Uc(3,i,j,k)
          wU(3,j,i,k) = this%Uc(2,i,j,k)
          wU(4,j,i,k) = this%Uc(6,i,j,k)
          wU(5,j,i,k) = this%Uc(5,i,j,k)
          wU(6,j,i,k) = this%Uc(4,i,j,k)
        end do
      end do
    end do

    call viehlli2d_step_core(                                                 &
      this%cfg%jmino_, this%cfg%jmaxo_, this%cfg%imino_, this%cfg%imaxo_,     &
      this%cfg%kmino_, this%cfg%kmaxo_, this%cfg%dymi, wU, wdU, dt,           &
      this%resolve_contacts                                                   &
    )

    do k = this%cfg%kmino_, this%cfg%kmaxo_
      do j = this%cfg%jmino_, this%cfg%jmaxo_
        do i = this%cfg%imino_, this%cfg%imaxo_
          this%dU(1,i,j,k) = this%dU(1,i,j,k) + wdU(1,j,i,k)
          this%dU(2,i,j,k) = this%dU(2,i,j,k) + wdU(3,j,i,k)
          this%dU(3,i,j,k) = this%dU(3,i,j,k) + wdU(2,j,i,k)
          this%dU(4,i,j,k) = this%dU(4,i,j,k) + wdU(6,j,i,k)
          this%dU(5,i,j,k) = this%dU(5,i,j,k) + wdU(5,j,i,k)
          this%dU(6,i,j,k) = this%dU(6,i,j,k) + wdU(4,j,i,k)
        end do
      end do
    end do

    deallocate(wU); deallocate(wdU);

  end subroutine viehlli2d_comp_dU_y

end module viehlli2d_class

