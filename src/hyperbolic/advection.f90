
!> contains necessary functions for using the muscl class (or other general hyperbolic
!> solvers) to solve basic advection problems

!> written by John P Wakefield in December 2022

module hyperbolic_advection
  use precision,    only: WP
  use string,       only: str_medium
  use config_class, only: config
  use muscl_class,  only: muscl, constructor, eigenvals_ftype, rsolver_ftype, limiter_ftype
  implicit none

  real(WP), parameter :: advec_muscl_cflsafety = 0.99_WP
  real(WP), parameter :: advec_muscl_divzero_eps = 1e-12_WP
  character(len=str_medium), parameter :: advec_muscl_name = 'MUSCL_CONST_ADVEC'

contains

  !> factory
  function make_advec_muscl(cfg, N, limiter, velocity) result(solver)
    implicit none
    type(muscl) :: solver
    class(config), target, intent(in) :: cfg
    integer, intent(in) :: N
    integer(1), intent(in) :: limiter
    real(WP), dimension(3), intent(in) :: velocity
    character(len=str_medium) :: name_actual
    procedure(eigenvals_ftype), pointer :: evals_x_ptr, evals_y_ptr, evals_z_ptr
    procedure(rsolver_ftype), pointer :: rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr
    integer :: i

    name_actual = advec_muscl_name

    evals_x_ptr => advec_evals_x
    evals_y_ptr => advec_evals_y
    evals_z_ptr => advec_evals_z
    rsolv_x_ptr => advec_rsolv_x
    rsolv_y_ptr => advec_rsolv_y
    rsolv_z_ptr => advec_rsolv_z

    ! build solver
    solver = muscl(cfg, name_actual, N, 3, evals_x_ptr, evals_y_ptr,      &
      & evals_z_ptr, rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr, limiter,  &
      & advec_muscl_divzero_eps, advec_muscl_cflsafety)

    ! set velocity
    do i = 1, 3
      solver%params(i,:,:,:) = velocity(i)
    end do

  end function make_advec_muscl

  pure subroutine advec_evals_x(N, params, U, evals)
    implicit none
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: params
    real(WP), dimension(N), intent(in) :: U
    real(WP), dimension(N), intent(out) :: evals

    ! silence unused variable warning
    evals(1) = U(1)

    evals(:) = params(1)

  end subroutine advec_evals_x

  pure subroutine advec_evals_y(N, params, U, evals)
    implicit none
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: params
    real(WP), dimension(N), intent(in) :: U
    real(WP), dimension(N), intent(out) :: evals

    ! silence unused variable warning
    evals(1) = U(1)

    evals(:) = params(2)

  end subroutine advec_evals_y

  pure subroutine advec_evals_z(N, params, U, evals)
    implicit none
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: params
    real(WP), dimension(N), intent(in) :: U
    real(WP), dimension(N), intent(out) :: evals

    ! silence unused variable warning
    evals(1) = U(1)

    evals(:) = params(3)

  end subroutine advec_evals_z

  pure subroutine advec_rsolv_simple(v, Ul, Ur, rs)
    real(WP), intent(in) :: v
    real(WP), dimension(:), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs
    integer :: i

    rs(:,:) = 0.0_WP
    rs(:,1) = min(v, 0.0_WP)
    rs(:,2) = max(v, 0.0_WP)
    rs(:,3) = Ur(:) - Ul(:)

    do i = 1, size(Ul, 1)
      rs(i,i+4) = 1.0_WP
    end do

  end subroutine

  pure subroutine advec_rsolv_x(N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call advec_rsolv_simple(0.5_WP * (pl(1) + pr(1)), Ul, Ur, rs)

  end subroutine advec_rsolv_x

  pure subroutine advec_rsolv_y(N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call advec_rsolv_simple(0.5_WP * (pl(2) + pr(2)), Ul, Ur, rs)

  end subroutine advec_rsolv_y

  pure subroutine advec_rsolv_z(N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: N
    real(WP), dimension(:), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call advec_rsolv_simple(0.5_WP * (pl(3) + pr(3)), Ul, Ur, rs)

  end subroutine advec_rsolv_z

end module hyperbolic_advection

