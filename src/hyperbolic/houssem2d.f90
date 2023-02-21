
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

  real(WP), parameter :: RT3 = sqrt(3.0_WP)

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
  pure subroutine houssem2d_rhs(gvec, taup, fl_vel, U, rhs)
    implicit none
    real(WP), intent(in) :: taup
    real(WP), dimension(2), intent(in) :: gvec, fl_vel
    real(WP), dimension(6), intent(in) :: U
    real(WP), dimension(6), intent(out) :: rhs
    real(WP), dimension(2) :: pt_vel

    pt_vel(:) = U(2:3) / U(1)

    rhs(1) = 0.0_WP

    rhs(2:3) = U(1) * gvec
    rhs(2:3) = rhs(2:3) + U(1) / taup * (fl_vel - pt_vel)

    rhs(4) = 2 * pt_vel(1) * gvec(1)
    rhs(4) = rhs(4) + 2 * (pt_vel(1) * fl_vel(1) - U(4) / U(1)) / taup
    rhs(4) = rhs(4) * U(1)

    rhs(5) = gvec(1) * pt_vel(2) + gvec(2) * pt_vel(1)
    rhs(5) = rhs(5) + (pt_vel(1) * fl_vel(2) + pt_vel(2) * fl_vel(1) -        &
      & 2 * U(5) / U(1)) / taup
    rhs(5) = rhs(5) * U(1)

    rhs(6) = 2 * gvec(2) * pt_vel(2)
    rhs(6) = rhs(6) + 2 * (pt_vel(2) * fl_vel(2) - U(6) / U(1)) / taup
    rhs(6) = rhs(6) * U(1)

  end subroutine houssem2d_rhs

  !pure subroutine houssem2d_rhs_backward(gvec, taup, fl_vel, U, dt, rhs)
  !  implicit none
  !  real(WP), intent(in) :: taup, dt
  !  real(WP), dimension(2), intent(in) :: gvec, fl_vel
  !  real(WP), dimension(6), intent(in) :: U
  !  real(WP), dimension(6), intent(out) :: rhs
  !  real(WP), dimension(2) :: pt_vel
  !  real(WP) :: a, b

  !  a = 1.0_WP / taup + 1.0_WP / dt
  !  b = 1.0_WP / (1.0_WP + taup / dt)

  !  pt_vel(:) = U(2:3) / U(1)

  !  rhs(1) = 0.0_WP

  !  rhs(2:3) = U(1) * (a * fl_vel - b * pt_vel)
  !  rhs(2:3) = rhs(2:3) + U(1) * gvec

  !  rhs(4) = 


  !end subroutine houssem2d_rhs

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

    u_permute(:) = u(:); u_permute(2) = u(3); u_permute(3) = u(2);

    call houssem2d_evals_1d(u_permute, evals)

  end subroutine houssem2d_evals_y

  pure subroutine houssem2d_evals_z(N, params, U, evals)
    implicit none
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: evals
    real(WP), dimension(N) :: u_permute

    ! do nothing

  end subroutine houssem2d_evals_z

  !> roe average
  pure subroutine houssem2d_roeavg(Ul, Ur, Ua)
    implicit none
    real(WP), dimension(6), intent(in) :: Ul, Ur
    real(WP), dimension(6), intent(out) :: Ua
    real(WP) :: wl, wr, v1, v2, v1l, v1r, v2l, v2r
    real(WP) :: e11l, e11r, e12l, e12r, e22l, e22r

    v1l = Ul(2) / Ul(1); v1r = Ur(2) / Ur(1);
    v2l = Ul(3) / Ul(1); v2r = Ur(3) / Ur(1);
    e11l = Ul(4) / Ul(1); e11r = Ur(4) / Ur(1);
    e12l = Ul(5) / Ul(1); e12r = Ur(5) / Ur(1);
    e22l = Ul(6) / Ul(1); e22r = Ur(6) / Ur(1);

    wl = sqrt(Ul(1)); wr = sqrt(Ur(1)); Ua(1) = wl * wr;
    wl = wl / (wl + wr); wr = 1.0_WP - wl;

    v1 = wl * v1l + wr * v1r; v2 = wl * v2l + wr * v2r;

    Ua(2) = v1 * Ua(1); Ua(3) = v2 * Ua(1);

    Ua(4) = (3*e11l*Ua(1) + 3*e11l*Ul(1) + 3*e11r*Ua(1) + 3*e11r*Ur(1) - 2*Ua(1)*v1l**2 + 4*Ua(1)*v1l*v1r - 2*Ua(1)*v1r**2) / (3*(2*Ua(1) + Ul(1) + Ur(1)))
    Ua(5) = (3*e12l*Ua(1) + 3*e12l*Ul(1) + 3*e12r*Ua(1) + 3*e12r*Ur(1) - 2*Ua(1)*v1l*v2l + 2*Ua(1)*v1l*v2r + 2*Ua(1)*v1r*v2l - 2*Ua(1)*v1r*v2r) / (3*(2*Ua(1) + Ul(1) + Ur(1)))
    Ua(6) = (3*e22l*Ua(1) + 3*e22l*Ul(1) + 3*e22r*Ua(1) + 3*e22r*Ur(1) - 2*Ua(1)*v2l**2 + 4*Ua(1)*v2l*v2r - 2*Ua(1)*v2r**2) / (3*(2*Ua(1) + Ul(1) + Ur(1)))

    Ua(4:6) = Ua(1) * Ua(4:6)

  end subroutine houssem2d_roeavg

  !> eigenvalues
  pure subroutine houssem2d_evals_1d(U, lambda)
    implicit none
    real(WP), dimension(6), intent(in) :: U
    real(WP), dimension(6), intent(out) :: lambda
    real(WP) :: c, v1

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

    ! do nothing

  end subroutine houssem2d_rsolv_z

  !> Roe Solver
  pure subroutine houssem2d_rsolv_roe_1d(pl, pr, Ul, Ur, rs)
    implicit none
    real(WP), dimension(2), intent(in) :: pl, pr
    real(WP), dimension(6), intent(in) :: Ul, Ur
    real(WP), dimension(6,10), intent(out) :: rs
    real(WP), dimension(6) :: Ua, a, b
    real(WP) :: v1, v2, c, k, koc, cok
    real(WP) :: drhopn, drhopnv1, drhopnv2, de11, de12, de22

    ! Roe averages
    call houssem2d_roeavg(Ul, Ur, Ua)
    v1 = Ua(2) / Ua(1); v2 = Ua(3) / Ua(1);

    ! c and k
    !TODO decide whether to keep epsilon
    c = sqrt(Ua(4) / Ua(1) - v1**2 + 1e-13_WP)
    k = sqrt(Ua(5) / Ua(1) - v1 * v2 + 1e-13_WP)
    koc = k / c; cok = c / k;

    ! eigenvalues
    rs(:,2) = (/ v1 - RT3 * c, v1 - c, v1, v1, v1 + c, v1 + RT3 * c /)
    rs(:,1) = min(rs(:,2), 0.0_WP)
    rs(:,2) = rs(:,2) - rs(:,1)

    ! αs
    drhopn = Ur(1) - Ul(1);
    drhopnv1 = Ur(2) - Ul(2); drhopnv2 = Ur(3) - Ul(3);
    de11 = Ur(4) - Ul(4); de12 = Ur(5) - Ul(5); de22 = Ur(6) - Ul(6);
    rs(1,3) = k**2*(-RT3*c**3*drhopn*Ua(6)*v1 + RT3*c**3*drhopnv1*Ua(6) - c**2*de11*Ua(6) - c**2*drhopn*Ua(6)*v1**2 + 6*c**2*drhopn*k*v1*v2 + 2*c**2*drhopnv1*Ua(6)*v1 - 6*c**2*drhopnv1*k*v2 + 2*RT3*c*de11*k*v2 - 2*RT3*c*drhopn*k**2*v1 + 2*RT3*c*drhopn*k*v1**2*v2 + 2*RT3*c*drhopnv1*k**2 - 4*RT3*c*drhopnv1*k*v1*v2 - 2*de11*k**2 - 2*drhopn*k**2*v1**2 + 4*drhopnv1*k**2*v1) / (6*c**2*(-c**4*Ua(6)**2 - 4*c**2*Ua(6)*k**2 + 12*c**2*k**2*v2**2 - 4*k**4))
    rs(2,3) =  (-c**3*drhopn*v2 + c**3*drhopnv2 - c**2*de12 - c**2*drhopn*v1*v2 + c**2*drhopnv1*v2 + c**2*drhopnv2*v1 + c*drhopn*k*v1 - c*drhopnv1*k + de11*k + drhopn*k*v1**2 - 2*drhopnv1*k*v1) / (2*c**4)
    rs(3,3) = (3*c**2*drhopn - de11 - drhopn*v1**2 + 2*drhopnv1*v1)/(3*c**2)
    rs(4,3) = (3*c**4*de22 + 6*c**4*drhopn*v2**2 - 6*c**4*drhopnv2*v2 - c**2*de11*Ua(6) - 6*c**2*de12*k - c**2*drhopn*Ua(6)*v1**2 - 6*c**2*drhopn*k*v1*v2 + 2*c**2*drhopnv1*Ua(6)*v1 + 6*c**2*drhopnv1*k*v2 + 6*c**2*drhopnv2*k*v1 + 4*de11*k**2 + 4*drhopn*k**2*v1**2 - 8*drhopnv1*k**2*v1) / (3*c**4*v1**2+DIVZERO_EPS)
    rs(5,3) = (-c**3*drhopn*v2 + c**3*drhopnv2 + c**2*de12 + c**2*drhopn*v1*v2 - c**2*drhopnv1*v2 - c**2*drhopnv2*v1 + c*drhopn*k*v1 - c*drhopnv1*k - de11*k - drhopn*k*v1**2 + 2*drhopnv1*k*v1) / (2*c**4)
    rs(6,3) = k**2*(RT3*c**3*drhopn*Ua(6)*v1 - RT3*c**3*drhopnv1*Ua(6) - c**2*de11*Ua(6) - c**2*drhopn*Ua(6)*v1**2 + 6*c**2*drhopn*k*v1*v2 + 2*c**2*drhopnv1*Ua(6)*v1 - 6*c**2*drhopnv1*k*v2 - 2*RT3*c*de11*k*v2 + 2*RT3*c*drhopn*k**2*v1 - 2*RT3*c*drhopn*k*v1**2*v2 - 2*RT3*c*drhopnv1*k**2 + 4*RT3*c*drhopnv1*k*v1*v2 - 2*de11*k**2 - 2*drhopn*k**2*v1**2 + 4*drhopnv1*k**2*v1) / (6*c**2*(-c**4*Ua(6)**2 - 4*c**2*Ua(6)*k**2 + 12*c**2*k**2*v2**2 - 4*k**4))

    ! βs
    !TODO source terms
    rs(:,4) = 0.0_WP

    ! eigenvectors
    rs(:,6) = (/ 0.0_WP, 0.0_WP, c, 0.0_WP, c*(-c + v1), 2*(c*v2 - k) /)
    rs(:,7) = (/ 1.0_WP, v1, v2, v1**2, v1*v2, 0.0_WP /)
    rs(:,8) = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, v1**2 /)
    rs(:,9) = (/ 0.0_WP, 0.0_WP, c, 0.0_WP, c*(c + v1), 2*(c*v2 + k) /)
    rs(1,5) = cok**2*Ua(6) + 2*RT3*cok*v2 + 2.0_WP
    rs(2,5) = -RT3*c*cok**2*Ua(6) + cok**2*Ua(6)*v1 - 6*c*cok*v2 - 2*RT3*c +  &
      & 2*RT3*cok*v1*v2 + 2*v1
    rs(3,5) = cok**2*Ua(6)*v2 - RT3*cok*Ua(6) + 2*RT3*cok*v2**2 - 4*v2 -      &
      & 2*RT3*koc
    rs(4,5) = 3*c**2*cok**2*Ua(6) - 2*RT3*c*cok**2*Ua(6)*v1 +                 &
      & 6*RT3*c**2*cok*v2 + cok**2*Ua(6)*v1**2 + 6*c**2 - 12*c*cok*v1*v2 -    &
      & 4*RT3*c*v1 + 2*RT3*cok*v1**2*v2 + 2*v1**2
    rs(5,5) = -RT3*c*cok**2*Ua(6)*v2 + 3*c*cok*Ua(6) + cok**2*Ua(6)*v1*v2 -   &
      & 6*c*cok*v2**2 - RT3*cok*Ua(6)*v1 + 4*RT3*c*v2 + 2*RT3*cok*v1*v2**2 +  &
      & 6*k - 4*v1*v2 - 2*RT3*koc*v1
    a(1) = cok**2*Ua(6) + 2.0_WP
    a(2) = cok**2*Ua(6)*v1 - 6*c*cok*v2 + 2*v1
    a(3) = cok**2*Ua(6)*v2 - 4*v2
    a(4) = cok**2*Ua(6)*(3*c**2 + v1**2) - 12*cok**2*k*v1*v2 + 2*(3*c**2 + v1**2)
    a(5) = cok**2*Ua(6)*v1*v2 + 3*cok**2*k*(Ua(6) - 2*v2**2) + 6*k - 4*v1*v2
    a(6) = cok**2*Ua(6)**2 + 4*Ua(6) - 12*v2**2 + 4*koc**2
    b(1) = -2*cok*v2
    b(2) = cok*(c*cok*Ua(6) + 2*k - 2*v1*v2)
    b(3) = cok*(Ua(6) - 2*v2**2) + 2*koc
    b(4) = 2*cok*(c*cok*Ua(6)*v1 + 2*k*v1 - v2*(3*c**2 + v1**2))
    b(5) = c*cok**2*Ua(6)*v2 - 4*c*v2 + cok*v1*(Ua(6) - 2*v2**2) + 2*koc*v1
    b(6) = 0.0_WP
    rs(:, 5) = a + RT3 * b
    rs(:,10) = a - RT3 * b

  end subroutine houssem2d_rsolv_roe_1d

end module hyperbolic_houssem2d

