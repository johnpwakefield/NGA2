



module filterstats
  implicit none
  private

  public :: FLT_GAUSSIAN, FLT_BOX, FLT_NUM_PARAMS
  public :: filter_info_row

  !> imaginary constant
  complex(WP), parameter :: im = (0.0_WP, 1.0_WP)

  !> give all filters six parameters, not all of which may be used
  integer, parameter :: FLT_NUM_PARAMS = 6

  !> filter types
  integer, parameter :: FLT_GAUSSIAN = 1
  integer, parameter :: FLT_BOX      = 2

  !> helper structure to read in filter info from files
  type :: filter_info_row
    character(len=str_medium) :: out_fname
    character(len=str_medium) :: typename
    real(WP), dimension(FLT_NUM_PARAMS) :: params
  end type

  !> microscopic statistics and corresponding monitor file
  type :: microstats
    character(len=str_medium) :: out_fname
    type(monitor) :: mon
    integer :: step, Np
    real(WP) :: rhop, dp, VF, taup
    real(WP), dimension(3) :: VB, slp, drg  ! mean velocity, slip, and drag x, y, z
    real(WP), dimension(6) :: VV            ! xx, xy, xz, yy, yz, zz
  contains
    procedure :: init => microstats_init_mon
    final :: microstats_destroy_mon
  end type microstats

contains


  !> microscopic statistics

  subroutine microstats_init_mon(this, amroot)
    implicit none
    class(microstats), intent(inout) :: this
    logical, intent(in) :: amroot

    this%mon = monitor(amroot, this%out_fname)
    call this%mon%add_column(this%step,   'step')
    call this%mon%add_column(this%dp,     'dp')
    call this%mon%add_column(this%Np,     'Np')
    call this%mon%add_column(this%rhop,   'rhop')
    call this%mon%add_column(this%VF,     'VF')
    call this%mon%add_column(this%taup,   'taup')
    call this%mon%add_column(this%VB(1),  'VBx')
    call this%mon%add_column(this%VB(2),  'VBy')
    call this%mon%add_column(this%VB(3),  'VBz')
    call this%mon%add_column(this%slp(1), 'slpx')
    call this%mon%add_column(this%slp(2), 'slpy')
    call this%mon%add_column(this%slp(3), 'slpz')
    call this%mon%add_column(this%drg(1), 'drgx')
    call this%mon%add_column(this%drg(2), 'drgy')
    call this%mon%add_column(this%drg(3), 'drgz')
    call this%mon%add_column(this%VV(1),  'VVxx')
    call this%mon%add_column(this%VV(2),  'VVxy')
    call this%mon%add_column(this%VV(3),  'VVxz')
    call this%mon%add_column(this%VV(4),  'VVyy')
    call this%mon%add_column(this%VV(5),  'VVyz')
    call this%mon%add_column(this%VV(6),  'VVzz')

  end subroutine microstats_init_mon

  subroutine microstats_destroy_mon(this)
    implicit none
    type(microstats), intent(inout) :: this

    ! nothing to do

  end subroutine microstats_destroy_mon


end module filterstats

