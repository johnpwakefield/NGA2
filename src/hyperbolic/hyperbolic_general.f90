
!> module containing definitions relevant to a wide range of hyperbolic
!> solvers, e.g. function type definitions and limiter definitions

!> Originally written by John P Wakefield in Feb 2023.

module hyperbolic
  use precision, only: WP
  implicit none
  public

  !> limiter names
  integer(1), parameter, public :: UPWIND   = 0_1
  integer(1), parameter, public :: LAXWEND  = 1_1
  integer(1), parameter, public :: BEAMWARM = 2_1
  integer(1), parameter, public :: FROMM    = 3_1
  integer(1), parameter, public :: MINMOD   = 4_1
  integer(1), parameter, public :: SUPERBEE = 5_1
  integer(1), parameter, public :: MC       = 6_1
  integer(1), parameter, public :: VANLEER  = 7_1

  !> function types
  interface

    pure subroutine eigenvals_ftype(P, N, params, u, evals)
      use precision, only: WP
      implicit none
      integer, intent(in) :: P, N
      real(WP), dimension(P), intent(in) :: params
      real(WP), dimension(N), intent(in) :: u
      real(WP), dimension(N), intent(out) :: evals
    end subroutine

    pure subroutine rsolver_ftype(P, N, pl, ul, pr, ur, rs)
      use precision, only: WP
      implicit none
      integer, intent(in) :: P, N
      real(WP), dimension(P), intent(in) :: pl, pr
      real(WP), dimension(N), intent(in) :: ul, ur
      real(WP), dimension(:,:), intent(out) :: rs
    end subroutine

    pure subroutine flux_ftype(P, N, params, u, f)
      use precision, only: WP
      integer, intent(in) :: P, N
      real(WP), dimension(P), intent(in) :: params
      real(WP), dimension(N), intent(in) :: u
      real(WP), dimension(N), intent(out) :: f
    end subroutine

    pure function limiter_ftype(r) result(phi)
      use precision, only: WP
      implicit none
      real(WP), intent(in) :: r
      real(WP) :: phi
    end function

  end interface

contains

  !! limiter lookup function
  function get_limiter(i) result(lim)
    use messager, only: die
    implicit none
    integer(1), intent(in) :: i
    procedure(limiter_ftype), pointer :: lim

    lim => NULL()

    select case (i)
    case (UPWIND)
      lim => limiter_upwind
    case (LAXWEND)
      lim => limiter_laxwend
    case (BEAMWARM)
      lim => limiter_beamwarm
    case (FROMM)
      lim => limiter_fromm
    case (MINMOD)
      lim => limiter_minmod
    case (SUPERBEE)
      lim => limiter_superbee
    case (MC)
      lim => limiter_mc
    case (VANLEER)
      lim => limiter_vanleer
    case default
      call die("could not find limiter")
    end select

  end function get_limiter

  !! limiter definitions
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

end module hyperbolic

