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
         integer :: i, j, k, nx, ny
         real(WP), dimension(:), allocatable :: x, y, z
         real(WP) :: pi

         pi = 4 * atan(1.0_WP)

         ! read in grid definition
         call param_read('nx', nx)
         call param_read('ny', ny)

         ! allocate
         allocate(x(nx+1), y(ny+1), z(1+1))

         ! create simple rectilinear grid
         do i = 1, nx+1
            x(i) = real(i-1,WP) / real(nx, WP)
         end do
         do j = 1, ny+1
            y(j) = real(j-1,WP) / real(ny, WP)
         end do
         do k = 1, 1+1
            z(k) = real(k-1,WP) / real(1, WP)
         end do

         ! generate serial grid object
         grid=sgrid(coord=cartesian, no=2, x=x, y=x, z=x, xper=.true.,        &
           & yper=.true., zper=.true.)

         deallocate(x, y, z)

      end block create_grid

      ! create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition

         ! read in partition
         call param_read('Partition', partition, short='p')

         ! create partitioned grid
         cfg = config(grp=group, decomp=partition, grid=grid)

      end block create_cfg

      ! Create masks for this config
      create_walls: block
         cfg%VF = 1.0_WP
      end block create_walls

   end subroutine geometry_init

end module geometry

