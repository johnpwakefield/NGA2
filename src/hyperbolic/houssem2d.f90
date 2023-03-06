
!> contains necessary functions for using the muscl class (or other general
!> hyperbolic solvers) to solve the Houssem 2019's two-fluid equations in 2d

!> Written by John Wakefield in January 2023

module hyperbolic_houssem2d
  use precision,    only: WP
  use string,       only: str_medium
  use config_class, only: config
  use muscl_class,  only: muscl, VANLEER, eigenvals_ftype, rsolver_ftype,     &
    & limiter_ftype
  implicit none

  real(WP), parameter :: CFL_SAFETY = 0.92_WP
  real(WP), parameter :: DIVZERO_EPS = 1e-9_WP
  character(len=str_medium), parameter :: houssem2d_muscl_name =              &
    & 'MUSCL_HOUSSEM2D'

contains

  !> factory
  function make_houssem2d_muscl(cfg, limiter) result(solver)
    implicit none
    type(muscl) :: solver
    class(config), target, intent(in) :: cfg
    integer(1), optional, intent(in) :: limiter
    integer(1) :: limiter_actual
    character(len=str_medium) :: name_actual
    procedure(eigenvals_ftype), pointer :: evals_x_ptr, evals_y_ptr, evals_z_ptr
    procedure(rsolver_ftype), pointer :: rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr

    if (present(limiter)) then
      limiter_actual = limiter
    else
      limiter_actual = VANLEER
    end if

    name_actual = houssem2d_muscl_name

    evals_x_ptr => houssem2d_evals_x
    evals_y_ptr => houssem2d_evals_y
    evals_z_ptr => houssem2d_evals_z
    rsolv_x_ptr => houssem2d_rsolv_x
    rsolv_y_ptr => houssem2d_rsolv_y
    rsolv_z_ptr => houssem2d_rsolv_z

    ! build solver
    solver = muscl(cfg, name_actual, 6, 2, evals_x_ptr, evals_y_ptr,          &
      & evals_z_ptr, rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr, limiter_actual,   &
      & DIVZERO_EPS, CFL_SAFETY)

    ! set velocity mask
    solver%vel_mask_x(:) = (/ .false., .true., .false., .false., .false., .false. /)
    solver%vel_mask_y(:) = (/ .false., .false., .true., .false., .false., .false. /)
    solver%vel_mask_z(:) = .false.

  end function make_houssem2d_muscl

  !> source terms
  pure subroutine houssem2d_rhs(eqPp, gvec, taup, fl_vel, U, rhs)
    implicit none
    real(WP), dimension(3), intent(in) :: eqPp
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

    rhs(4:6) = rhs(4:6) + (2 * U(1) / taup) * eqPp

  end subroutine houssem2d_rhs

  pure subroutine houssem2d_backeuler(eqPp, gvec, taup, fl_vel, dt, U, dU)
    implicit none
    real(WP), dimension(3), intent(in) :: eqPp
    real(WP), intent(in) :: taup, dt
    real(WP), dimension(2), intent(in) :: gvec, fl_vel
    real(WP), dimension(6), intent(in) :: U
    real(WP), dimension(6), intent(out) :: dU
    real(WP) :: a, b
    real(WP), dimension(5) :: c

    a = 1.0_WP / (dt / taup + 1.0_WP)
    b = 1.0_WP / (2 * dt**2 + 3 * dt * taup + taup**2)

    c(1:2) = U(1) * fl_vel(1:2)
    c(3:5) = U(1) * eqPp
    c(:) = (dt / taup) * c(:)
    c(:) = c(:) + U(2:6)

    dU(1) = U(1)
    dU(2:3) = a * c(1:2)
    dU(4) = 2 * dt * taup * fl_vel(1) * b * c(1)
    dU(5) = dt * taup * b * (fl_vel(2) * c(1) + fl_vel(1) * c(2))
    dU(6) = 2 * dt * taup * fl_vel(2) * b * c(2)
    dU(4:6) = dU(4:6) + taup / (2 * dt + taup) * c(3:5)

    dU(:) = dU(:) - U(:)

  end subroutine houssem2d_backeuler

  pure subroutine houssem2d_evals_x(N, params, U, evals)
    implicit none
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: evals

    call houssem2d_evals_1d(U, evals)

  end subroutine houssem2d_evals_x

  pure subroutine houssem2d_evals_y(N, params, U, evals)
    implicit none
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: evals
    real(WP), dimension(6) :: u_permute

    u_permute(:) = u(:)
    u_permute(2) = u(3); u_permute(3) = u(2);
    u_permute(4) = u(6); u_permute(6) = u(4);

    call houssem2d_evals_1d(u_permute, evals)

  end subroutine houssem2d_evals_y

  pure subroutine houssem2d_evals_z(N, params, U, evals)
    implicit none
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: evals
    real(WP), dimension(N) :: u_permute

    evals(:) = 0.0_WP

  end subroutine houssem2d_evals_z

  !> roe average
  pure subroutine houssem2d_roeavg(Ul, Ur, Ua)
    implicit none
    real(WP), dimension(6), intent(in) :: Ul, Ur
    real(WP), dimension(6), intent(out) :: Ua
    real(WP) :: wl, wr, v1, v2, v1l, v1r, v2l, v2r
    real(WP) :: e11l, e11r, e12l, e12r, e22l, e22r
    real(WP), parameter :: TT = 2.0_WP / 3.0_WP

    v1l = Ul(2) / Ul(1); v1r = Ur(2) / Ur(1);
    v2l = Ul(3) / Ul(1); v2r = Ur(3) / Ur(1);
    e11l = Ul(4) / Ul(1); e11r = Ur(4) / Ur(1);
    e12l = Ul(5) / Ul(1); e12r = Ur(5) / Ur(1);
    e22l = Ul(6) / Ul(1); e22r = Ur(6) / Ur(1);

    wl = sqrt(Ul(1)); wr = sqrt(Ur(1)); Ua(1) = wl * wr;
    wl = wl / (wl + wr); wr = 1.0_WP - wl;

    v1 = wl * v1l + wr * v1r; v2 = wl * v2l + wr * v2r;

    Ua(2) = v1 * Ua(1); Ua(3) = v2 * Ua(1);

    Ua(4) = (e11l*(Ua(1) + Ul(1)) + e11r*(Ua(1) + Ur(1)) - TT*Ua(1)*(v1l -    &
      v1r)**2)
    Ua(5) = (e12l*(Ua(1) + Ul(1)) + e12r*(Ua(1) + Ur(1)) + TT*Ua(1)*(-v1l*v2l &
      + v1l*v2r + v1r*v2l - v1r*v2r))
    Ua(6) = (e22l*(Ua(1) + Ul(1)) + e22r*(Ua(1) + Ur(1)) - TT*Ua(1)*(v2l -    &
      v2r)**2)

    Ua(4:6) = Ua(1) / (2*Ua(1) + Ul(1) + Ur(1)) * Ua(4:6)

  end subroutine houssem2d_roeavg

  !> eigenvalues
  pure subroutine houssem2d_evals_1d(U, lambda)
    implicit none
    real(WP), dimension(6), intent(in) :: U
    real(WP), dimension(6), intent(out) :: lambda
    real(WP) :: c, v1
    real(WP), parameter :: RT3 = sqrt(3.0_WP)

    v1 = U(2) / U(1); c = sqrt(U(4) / U(1) - v1**2);

    lambda = (/ v1 - RT3 * c, v1 - c, v1, v1, v1 + c, v1 + RT3 * c /)

  end subroutine houssem2d_evals_1d

  pure subroutine houssem2d_rsolv_x(N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call houssem2d_rsolv_roe_1d(pl, pr, Ul, Ur, rs)

  end subroutine houssem2d_rsolv_x

  pure subroutine houssem2d_rsolv_y(N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs
    real(WP), dimension(N) :: Uln, Urn
    real(WP), dimension(6) :: rs_row

    Uln(:) = Ul(:); Uln(2) = Ul(3); Uln(3) = Ul(2); Uln(4) = Ul(6); Uln(6) = Ul(4);
    Urn(:) = Ur(:); Urn(2) = Ur(3); Urn(3) = Ur(2); Urn(4) = Ur(6); Urn(6) = Ur(4);

    call houssem2d_rsolv_roe_1d(pl, pr, Uln, Urn, rs)

    rs_row(:) = rs(2,5:10); rs(2,5:10) = rs(3,5:10); rs(3,5:10) = rs_row(:);
    rs_row(:) = rs(4,5:10); rs(4,5:10) = rs(6,5:10); rs(6,5:10) = rs_row(:);

  end subroutine houssem2d_rsolv_y

  pure subroutine houssem2d_rsolv_z(N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs
    real(WP), dimension(N) :: Uln, Urn
    real(WP), dimension(6) :: rs_row

    rs(:,:) = 0.0_WP

  end subroutine houssem2d_rsolv_z

  !> Roe Solver
  pure subroutine houssem2d_rsolv_roe_1d(pl, pr, Ul, Ur, rs)
    implicit none
    real(WP), dimension(2), intent(in) :: pl, pr
    real(WP), dimension(6), intent(in) :: Ul, Ur
    real(WP), dimension(6,10), intent(out) :: rs
    real(WP), dimension(6) :: Ua, a, b
    real(WP) :: v1, v2, c, k, koc, cok, e11, e12, e22
    real(WP) :: drhopn, drhopnv1, drhopnv2, de11, de12, de22
    real(WP), parameter :: RT3 = sqrt(3.0_WP)
    real(WP), parameter :: EPS = 1e-8_WP

    ! Roe averages
    call houssem2d_roeavg(Ul, Ur, Ua)
    v1 = Ua(2) / Ua(1); v2 = Ua(3) / Ua(1);

    ! c and k
    c = sqrt(Ua(4) / Ua(1) - v1**2 + EPS)
    k = Ua(5) / Ua(1) - v1 * v2
    koc = k / c; cok = c / k;

    ! eigenvalues
    rs(:,2) = (/ v1 - RT3 * c, v1 - c, v1, v1, v1 + c, v1 + RT3 * c /)
    rs(:,1) = min(rs(:,2), 0.0_WP)
    rs(:,2) = rs(:,2) - rs(:,1)

    ! αs
    drhopn = Ur(1) - Ul(1);
    drhopnv1 = Ur(2) - Ul(2); drhopnv2 = Ur(3) - Ul(3);
    e11 = Ua(4) / Ua(1); e12 = Ua(5) / Ua(1); e22 = Ua(6) / Ua(1);
    de11 = Ur(4) / Ur(1) - Ul(4) / Ul(1)
    de12 = Ur(5) / Ur(1) - Ul(5) / Ul(1)
    de22 = Ur(6) / Ur(1) - Ul(6) / Ul(1)
    rs(1,3) = k**2*(-RT3*c**3*drhopn*e22*v1 + RT3*c**3*drhopnv1*e22 - c**2*de11*e22 - c**2*drhopn*e22*v1**2 + 6*c**2*drhopn*k*v1*v2 + 2*c**2*drhopnv1*e22*v1 - 6*c**2*drhopnv1*k*v2 + 2*RT3*c*de11*k*v2 - 2*RT3*c*drhopn*k**2*v1 + 2*RT3*c*drhopn*k*v1**2*v2 + 2*RT3*c*drhopnv1*k**2 - 4*RT3*c*drhopnv1*k*v1*v2 - 2*de11*k**2 - 2*drhopn*k**2*v1**2 + 4*drhopnv1*k**2*v1)/(6*c**2*(-c**4*e22**2 - 4*c**2*e22*k**2 + 12*c**2*k**2*v2**2 - 4*k**4))
    rs(2,3) = (-c**3*drhopn*v2 + c**3*drhopnv2 - c**2*de12 - c**2*drhopn*v1*v2 + c**2*drhopnv1*v2 + c**2*drhopnv2*v1 + c*drhopn*k*v1 - c*drhopnv1*k + de11*k + drhopn*k*v1**2 - 2*drhopnv1*k*v1)/(2*c**4)
    rs(3,3) = (3*c**2*drhopn - de11 - drhopn*v1**2 + 2*drhopnv1*v1)/(3*c**2)
    rs(4,3) = (3*c**4*de22 + 6*c**4*drhopn*v2**2 - 6*c**4*drhopnv2*v2 - c**2*de11*e22 - 6*c**2*de12*k - c**2*drhopn*e22*v1**2 - 6*c**2*drhopn*k*v1*v2 + 2*c**2*drhopnv1*e22*v1 + 6*c**2*drhopnv1*k*v2 + 6*c**2*drhopnv2*k*v1 + 4*de11*k**2 + 4*drhopn*k**2*v1**2 - 8*drhopnv1*k**2*v1)/(3*c**4*v1**2)
    rs(5,3) = (-c**3*drhopn*v2 + c**3*drhopnv2 + c**2*de12 + c**2*drhopn*v1*v2 - c**2*drhopnv1*v2 - c**2*drhopnv2*v1 + c*drhopn*k*v1 - c*drhopnv1*k - de11*k - drhopn*k*v1**2 + 2*drhopnv1*k*v1)/(2*c**4)
    rs(6,3) = k**2*(RT3*c**3*drhopn*e22*v1 - RT3*c**3*drhopnv1*e22 - c**2*de11*e22 - c**2*drhopn*e22*v1**2 + 6*c**2*drhopn*k*v1*v2 + 2*c**2*drhopnv1*e22*v1 - 6*c**2*drhopnv1*k*v2 - 2*RT3*c*de11*k*v2 + 2*RT3*c*drhopn*k**2*v1 - 2*RT3*c*drhopn*k*v1**2*v2 - 2*RT3*c*drhopnv1*k**2 + 4*RT3*c*drhopnv1*k*v1*v2 - 2*de11*k**2 - 2*drhopn*k**2*v1**2 + 4*drhopnv1*k**2*v1)/(6*c**2*(-c**4*e22**2 - 4*c**2*e22*k**2 + 12*c**2*k**2*v2**2 - 4*k**4))

    ! βs
    rs(:,4) = 0.0_WP

    ! eigenvectors
    rs(:,6) = (/ 0.0_WP, 0.0_WP, c, 0.0_WP, c*(-c + v1), 2*(c*v2 - k) /)
    rs(:,7) = (/ 1.0_WP, v1, v2, v1**2, v1*v2, 0.0_WP /)
    rs(:,8) = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, v1**2 /)
    rs(:,9) = (/ 0.0_WP, 0.0_WP, c, 0.0_WP, c*(c + v1), 2*(c*v2 + k) /)
    a(1) = cok**2*e22 + 2.0_WP
    a(2) = cok**2*e22*v1 - 6*c*cok*v2 + 2*v1
    a(3) = cok**2*e22*v2 - 4*v2
    a(4) = 3*c**2*cok**2*e22 + cok**2*e22*v1**2 + 6*c**2 - 12*c*cok*v1*v2 + 2*v1**2
    a(5) = 3*c*cok*e22 + cok**2*e22*v1*v2 - 6*c*cok*v2**2 + 6*k - 4*v1*v2
    a(6) = cok**2*e22**2 + 4*e22 - 12*v2**2 + 4*koc**2
    b(1) = -2*cok*v2
    b(2) = c*(cok**2*e22 + 2.0_WP - 2/k*v1*v2)
    b(3) = cok*(e22 - 2*v2**2) + 2*koc
    b(4) = 2*c*(cok**2*e22*v1 - 3*c*cok*v2 + 2*v1 - v1**2*v2/k)
    b(5) = c*cok**2*e22*v2 + cok*e22*v1 - 4*c*v2 - 2*cok*v1*v2**2 + 2*koc*v1
    b(6) = 0.0_WP
    rs(:, 5) = a - RT3 * b
    rs(:,10) = a + RT3 * b

  end subroutine houssem2d_rsolv_roe_1d

end module hyperbolic_houssem2d

