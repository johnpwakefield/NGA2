!> Various definitions and tools for initializing NGA2 config
module geometry
  use config_class, only: config
  use precision,    only: WP
  implicit none
  private

  !> Single config
  type(config), public :: cfg

  public :: geometry_init

contains


  !> Initialization of problem geometry
  subroutine geometry_init
    use sgrid_class, only: sgrid
    use param,       only: param_read
    implicit none
    type(sgrid) :: grid

    ! Create a grid from input params
    create_grid: block
      use sgrid_class, only: cartesian
      real(WP), dimension(3) :: L
      integer, dimension(3) :: N
      real(WP), dimension(:,:), allocatable :: x
      integer :: d, i

      call param_read('L',L)
      call param_read('N',N)

      allocate(x(maxval(N)+1,3))

      do d = 1, 3
        do i = 1, N(d)+1
          x(i,d) = real(i - 1, WP) / real(N(d), WP) * L(d)
        end do
      end do

      grid = sgrid(coord=cartesian, no=1,                                     &
        x=x(1:(N(1)+1),1), y=x(1:(N(2)+1),2), z=x(1:(N(3)+1),3),              &
        xper=.true., yper=.true., zper=.true., name='NPYTEST')

    end block create_grid

    ! Create a config from that grid on our entire group
    create_cfg: block
      use parallel, only: group
      integer, dimension(3) :: partition

      ! Read in partition
      call param_read('Partition',partition,short='p')

      ! Create partitioned grid
      cfg = config(grp=group, decomp=partition, grid=grid)

    end block create_cfg

    ! Create masks for this config
    create_walls: block
      cfg%VF=1.0_WP
    end block create_walls

  end subroutine geometry_init


end module geometry
