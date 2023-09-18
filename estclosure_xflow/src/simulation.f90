!> We use the whole process pool for everything and take steps in each domain
!> sequentially rather than trying to load balance.

!>  - xflow is in the x direction
!>  - gravity is in the y direction
!>  - hitstats collects data in the xy and xz planes

!> TODO put author statement

!> TODO make xflowstats module that does hitstats but for the xflow window

!> TODO functionality in EC
!> will be moved elsewhere (maybe here, maybe another class).  In a perfect
!> world I would take the time to carefully organize this into an MVC
!> structure, but we don't live in a perfect world.

module simulation
  use precision,            only: WP
  use parallel,             only: rank, amroot
  use fft3d_class,          only: fft3d
  use fft2d_class,          only: fft2d
  use incomp_class,         only: incomp
  use timetracker_class,    only: timetracker
  use sgrid_class,          only: sgrid, cartesian
  use config_class,         only: config
  use param,                only: param_read
  use ensight_class,        only: ensight
  use partmesh_class,       only: partmesh
  use npy_class,            only: npy
  use event_class,          only: periodic_event, threshold_event
  use monitor_class,        only: monitor
  use datafile_class,       only: datafile
  !use cubestats_class,      only: cubestats
  !use xflwstats_class,      only: xflwstats
  use eccontroller_class,   only: eccontroller_mesh
  use lpt_class,            only: lpt
  use string,               only: str_medium
  use coupler_class,        only: coupler
  use lptcoupler_class,     only: lptcoupler
  implicit none
  private


  real(WP), parameter :: FILTER_MESH_RATIO = 3.5_WP
  real(WP), parameter :: GRID_OVERLAP = 0.02_WP


  !> The easiest way to get data from the cube is probably to just evolve the
  !> cube with the same cross velocity (in the y direction) as the crossflow
  !> part; that way the coupler doesn't have to move and we don't have to move
  !> or reinterpolate data manually.  In addition, it's easy enough to
  !> dynamically adjust this velocity to match some multiple of the cube's rms
  !> velocity.
  real(WP) :: xflowmult, xflowvel


  !> we use 1 for a particle that is going to be cleaned up (as is enforced
  !> by lpt), 0 for a regular particle, and 2 for a particle that has been
  !> copied and should not be copied again
  real(WP), dimension(2) :: part_reset_bdrys


  !> Domain configuration

  type :: gendomain
    character(len=str_medium) :: desc
    type(config) :: cfg
    type(incomp) :: fs
    real(WP) :: cfl, meanU, meanV, meanW
    real(WP), dimension(:,:,:), allocatable :: rho, resU, resV, resW, Ui, Vi, Wi
    real(WP), dimension(:,:,:,:), allocatable :: SR
    real(WP), dimension(:,:,:,:,:), allocatable :: gradu
    type(lpt)      :: lp
    type(monitor)  :: mfile, cflfile, lptfile
    type(partmesh) :: pmesh
    type(ensight)  :: ens_out
  end type gendomain

  type, extends(gendomain) :: cubedomain
    type(fft3d) :: ps
    real(WP) :: urms, TKE, EPS, eta, Re_lambda,                               &
      EPS_ratio, TKE_ratio, ell_ratio, dx_eta, forcingtimescale, Stk, phiinf, &
      Wovk, EPSp
    !TODO re-enable
    !type(cubestats) :: stats
    type(monitor) :: hitfile
  end type cubedomain

  type, extends(gendomain) :: xflowdomain
    type(fft2d) :: ps
    !TODO re-enable
    !type(xflwstats) :: stats
  end type xflowdomain

  !TODO move this to xflowstats
  !type :: xflowslice
  !  type(config)   :: cfg
  !  type(partmesh) :: pmesh
  !  type(npy)      :: npy_out
  !end type xflowslice

  type(cubedomain)  :: cube
  type(xflowdomain) :: xflw
  !TODO move this to xflowstats
  !type(xflowslice)  :: xflwslice


  !> Time-related
  type(timetracker), public :: time
  type(periodic_event)      :: ens_evt
  logical                   :: ens_at_ints


  !TODO update EC for this case, also update cube to not compute covs
  !> Closure Estimation
  type(eccontroller_mesh)  :: ec
  type(threshold_event)    :: ec_evt
  logical                  :: ec_done


  !> Turbulent and forcing parameters in cube
  !> For simplicity, we assume particle diameters, fluid viscosity, and
  !> gravity are equal in both domains.
  real(WP), dimension(7) :: ec_params   ! rhof, rhop, ktarget, epstarget, nu, dp, g
  real(WP) :: EPS_target, TKE_target, nu, dp, G     ! G is the forcing constant


  !> Couplers
  type(coupler)             :: dom_cpl
  type(lptcoupler)          :: dom_lptcpl


  public :: simulation_init, simulation_run, simulation_final


  !> Wallclock time for monitoring
  type :: timer; real(WP) :: time_in, time, percent; end type timer;
  type(timer) :: wt_total, wt_cube_vel, wt_cube_pres, wt_cube_lpt,            &
    wt_xflw_vel, wt_xflw_pres, wt_xflw_lpt, wt_stat, wt_force, wt_cpl, wt_rest
  type(monitor) :: tfile


contains


  !> Turbulent statistics and parameters
  !TODO should this be in hitstats?

  subroutine compute_cube_stats()
    use mathtools, only: pi
    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
    use parallel,  only: MPI_REAL_WP
    real(WP) :: myTKE, myEPS, myEPSp
    integer :: i, j, k, ierr

    ! Compute mean velocities
    call cube%fs%cfg%integrate(A=cube%fs%U, integral=cube%meanU)
    cube%meanU = cube%meanU / cube%fs%cfg%vol_total
    call cube%fs%cfg%integrate(A=cube%fs%V, integral=cube%meanV)
    cube%meanV = cube%meanV / cube%fs%cfg%vol_total
    call cube%fs%cfg%integrate(A=cube%fs%W, integral=cube%meanW)
    cube%meanW = cube%meanW / cube%fs%cfg%vol_total

    ! Compute strainrate and grad(u)
    call cube%fs%get_strainrate(SR=cube%SR)
    call cube%fs%get_gradu(gradu=cube%gradu)

    ! compute TKE, EPS, EPSp, forcing_offset
    myTKE = 0.0_WP; myEPS = 0.0_WP; myEPSp = 0.0_WP;
    do k = cube%fs%cfg%kmin_, cube%fs%cfg%kmax_
      do j = cube%fs%cfg%jmin_, cube%fs%cfg%jmax_
        do i = cube%fs%cfg%imin_, cube%fs%cfg%imax_
          myTKE = myTKE + 0.5_WP*((cube%fs%U(i,j,k)-cube%meanU)**2 +              &
            (cube%fs%V(i,j,k)-cube%meanV)**2 + (cube%fs%W(i,j,k)-cube%meanW)**2   &
            )*cube%fs%cfg%vol(i,j,k)
          myEPS = myEPS + 2.0_WP*cube%fs%visc(i,j,k)*cube%fs%cfg%vol(i,j,k)*(     &
            cube%SR(1,i,j,k)**2+cube%SR(2,i,j,k)**2+cube%SR(3,i,j,k)**2+2.0_WP*(  &
            cube%SR(4,i,j,k)**2+cube%SR(5,i,j,k)**2+cube%SR(6,i,j,k)**2)          &
          ) / cube%fs%rho
          myEPSp = myEPSp + cube%fs%cfg%vol(i,j,k)*cube%fs%visc(i,j,k)*(                &
            cube%gradu(1,1,i,j,k)**2 + cube%gradu(1,2,i,j,k)**2 + cube%gradu(1,3,i,j,k)**2 + &
            cube%gradu(2,1,i,j,k)**2 + cube%gradu(2,2,i,j,k)**2 + cube%gradu(2,3,i,j,k)**2 + &
            cube%gradu(3,1,i,j,k)**2 + cube%gradu(3,2,i,j,k)**2 + cube%gradu(3,3,i,j,k)**2)
        end do
      end do
    end do
    call MPI_ALLREDUCE(myTKE,cube%TKE,1,MPI_REAL_WP,MPI_SUM,cube%fs%cfg%comm,ierr)
    cube%TKE = cube%TKE / cube%fs%cfg%vol_total
    call MPI_ALLREDUCE(myEPS,cube%EPS,1,MPI_REAL_WP,MPI_SUM,cube%fs%cfg%comm,ierr)
    cube%EPS = cube%EPS / cube%fs%cfg%vol_total
    call MPI_ALLREDUCE(myEPSp,cube%EPSp,1,MPI_REAL_WP,MPI_SUM,cube%fs%cfg%comm,ierr)
    cube%EPSp = cube%EPSp / (cube%fs%cfg%vol_total * cube%fs%rho)

    ! urms
    cube%urms = sqrt(2.0_WP / 3 * cube%TKE)

    ! additional monitoring quantities
    cube%eta = (nu**3 / cube%EPS)**0.25_WP
    cube%Re_lambda = sqrt(20.0_WP / 3) * cube%TKE * (cube%eta / nu)**2
    cube%EPS_ratio = cube%EPS / EPS_target
    cube%TKE_ratio = cube%TKE / TKE_target
    cube%ell_ratio = (2.0_WP / 3 * cube%TKE)**1.5_WP / cube%EPS / (cube%cfg%vol_total)**(1.0_WP/3)
    cube%dx_eta = cube%fs%cfg%max_meshsize / cube%eta
    cube%forcingtimescale = 3 * EPS_target / (2 * TKE_target)

    ! nondimensional particle quantities
    cube%Stk = ec_params(2) / ec_params(1) * ec_params(6)**2                  &
      / (18 * cube%eta**2)
    cube%phiinf = cube%lp%np * pi * ec_params(6)**3 /                         &
      (6 * cube%fs%cfg%vol_total)
    cube%Wovk = ec_params(7) * ec_params(2) / ec_params(1) * ec_params(6)**2  &
      * cube%eta / (18 * nu**2)

  end subroutine compute_cube_stats

  ! params_primit - 7 items - rhof, rhop, ktarget, epstarget, nu, dp, g
  subroutine update_parameters()
    implicit none
    integer :: i

    ! fluid parameters
    cube%fs%rho = ec_params(1)
    cube%rho(:,:,:) = cube%fs%rho
    xflw%fs%rho = ec_params(1)
    xflw%rho(:,:,:) = xflw%fs%rho
    cube%fs%visc(:,:,:) = ec_params(1) * ec_params(5)
    xflw%fs%visc(:,:,:) = ec_params(1) * ec_params(5)

    ! particle parameters
    cube%lp%rho = ec_params(2); xflw%lp%rho = ec_params(2);
    do i = 1, cube%lp%np_; cube%lp%p(i)%d = ec_params(6); end do;
    do i = 1, xflw%lp%np_; xflw%lp%p(i)%d = ec_params(6); end do;
    call cube%lp%sync(); call xflw%lp%sync();

    ! gravity
    cube%lp%gravity(:) = (/ 0.0_WP, -ec_params(7), 0.0_WP /)
    xflw%lp%gravity(:) = (/ 0.0_WP, -ec_params(7), 0.0_WP /)

    ! set xflow velocity
    call cube%fs%cfg%integrate(A=cube%fs%U, integral=cube%meanU)
    cube%meanU = cube%meanU / cube%fs%cfg%vol_total
    xflowvel = xflowmult * sqrt(2.0_WP / 3 * ec_params(3))
    fix_mean_vel: block
      call cube%fs%cfg%integrate(A=cube%fs%U, integral=cube%meanU)
      cube%meanU = cube%meanU / cube%fs%cfg%vol_total
      call xflw%fs%cfg%integrate(A=xflw%fs%U, integral=xflw%meanU)
      xflw%meanU = xflw%meanU / xflw%fs%cfg%vol_total
      cube%fs%U(:,:,:) = cube%fs%U(:,:,:) + (xflowvel - cube%meanU)
      xflw%fs%U(:,:,:) = xflw%fs%U(:,:,:) + (xflowvel - xflw%meanU)
    end block fix_mean_vel

    ! store parameters that are difficult to access for reference elsewhere
    TKE_target = ec_params(3); EPS_target = ec_params(4);
    nu = ec_params(5); dp = ec_params(6);

    ! print current parameters
    print_statistics: block
      use string, only: str_long
      use messager, only: log
      character(len=str_long) :: message
      if (amroot) then
        write(message,'("At t = ",e16.8," updated turbulent parameters to:")') time%t; call log(message);
        write(message,'("    Fluid density: ",e12.6)') cube%fs%rho; call log(message);
        write(message,'("    Target TKE: ",e12.6)') TKE_target; call log(message);
        write(message,'("    Target EPS: ",e12.6)') EPS_target; call log(message);
        write(message,'("    Kinematic Viscosity: ",e12.6)') nu; call log(message);
      end if
    end block print_statistics

  end subroutine update_parameters

  subroutine init_vel(fs, ktarget, resU, resV, resW, mask)
    use random,    only: random_normal
    use mpi_f08,   only: mpi_allreduce, mpi_sum
    use mathtools, only: pi
    class(incomp), intent(inout) :: fs
    real(WP), intent(in) :: ktarget
    real(WP), dimension(fs%cfg%imino_:,fs%cfg%jmino_:,fs%cfg%kmino_:), intent(out) :: resU, resV, resW
    real(WP), dimension(fs%cfg%imino_:,fs%cfg%jmino_:,fs%cfg%kmino_:), intent(in), optional :: mask
    real(WP) :: urms0
    real(WP), dimension(3) :: ubar
    integer :: i,j,k

    ! Gaussian initial field
    urms0 = sqrt(2.0_WP / 3 * ktarget)
    do k = fs%cfg%kmin_, fs%cfg%kmax_
      do j = fs%cfg%jmin_, fs%cfg%jmax_
        do i = fs%cfg%imin_, fs%cfg%imax_
          if (present(mask)) then
            if (mask(i,j,k) .gt. 0.5_WP) cycle
          end if
          fs%U(i,j,k) = random_normal(m=0.0_WP, sd=urms0)
          fs%V(i,j,k) = random_normal(m=0.0_WP, sd=urms0)
          fs%W(i,j,k) = random_normal(m=0.0_WP, sd=urms0)
        end do
      end do
    end do
    call fs%cfg%sync(fs%U); call fs%cfg%sync(fs%V); call fs%cfg%sync(fs%W);

    ! Compute mean and remove it from the velocity field
    call fs%cfg%integrate(A=fs%U, integral=ubar(1))
    call fs%cfg%integrate(A=fs%V, integral=ubar(2))
    call fs%cfg%integrate(A=fs%W, integral=ubar(3))
    ubar(:) = ubar(:) / fs%cfg%vol_total
    fs%U = fs%U - ubar(1); fs%V = fs%V - ubar(2); fs%W = fs%W - ubar(3);

    ! Project to divergence-free space
    call fs%get_div()
    fs%psolv%rhs = -fs%cfg%vol * fs%div * fs%rho / time%dt
    fs%psolv%sol = 0.0_WP
    call fs%psolv%solve()
    call fs%shift_p(fs%psolv%sol)
    call fs%get_pgrad(fs%psolv%sol, resU, resV, resW)
    fs%P = fs%P+fs%psolv%sol
    fs%U = fs%U - time%dt * resU / fs%rho
    fs%V = fs%V - time%dt * resV / fs%rho
    fs%W = fs%W - time%dt * resW / fs%rho

  end subroutine init_vel


  !> Ensight output for cube and xflow

  subroutine setup_ens(d)
    implicit none
    class(gendomain), intent(inout) :: d

    d%pmesh = partmesh(nvar=3, nvec=2, name=trim(d%desc)//'_lpt')
    d%pmesh%varname(1) = "id"
    d%pmesh%varname(2) = "dp"
    d%pmesh%varname(3) = "flag"
    d%pmesh%vecname(1) = "vel"
    d%pmesh%vecname(2) = "fld"

    d%ens_out = ensight(cfg=d%cfg, name=trim(d%desc))
    call d%ens_out%add_vector('velocity', d%Ui, d%Vi, d%Wi)
    call d%ens_out%add_scalar('pressure', d%fs%P)
    call d%ens_out%add_particle('particles', d%pmesh)

  end subroutine setup_ens

  subroutine write_ens(d)
    implicit none
    class(gendomain), intent(inout) :: d
    integer :: i

    call d%fs%interp_vel(d%Ui, d%Vi, d%Wi)
    call d%pmesh%reset()
    call d%pmesh%set_size(d%lp%np_)
    do i = 1, d%lp%np_
      d%pmesh%pos(:,i) = d%lp%p(i)%pos
      d%pmesh%var(1,i) = d%lp%p(i)%id
      d%pmesh%var(2,i) = d%lp%p(i)%d
      d%pmesh%var(3,i) = d%lp%p(i)%flag
      d%pmesh%vec(:,1,i) = d%lp%p(i)%vel
      d%lp%p(i)%ind = d%cfg%get_ijk_global(d%lp%p(i)%pos, d%lp%p(i)%ind)
      d%pmesh%vec(:,2,i) = d%lp%cfg%get_velocity(pos=d%lp%p(i)%pos,         &
        i0=d%lp%p(i)%ind(1), j0=d%lp%p(i)%ind(2), k0=d%lp%p(i)%ind(3),      &
        U=d%fs%U, V=d%fs%V, W=d%fs%W)
    end do
    call d%ens_out%write_data(time%t)

  end subroutine write_ens


  !> Monitor output shared between the two domains

  subroutine setup_gendom_monitors(d)
    implicit none
    class(gendomain), intent(inout) :: d

    ! Create simulation monitor
    d%mfile = monitor(d%fs%cfg%amRoot, trim(d%desc) // '_simulation')
    call d%mfile%add_column(time%n,'Timestep number')
    call d%mfile%add_column(time%t,'Time')
    call d%mfile%add_column(time%dt,'Timestep size')
    call d%mfile%add_column(time%cfl,'Maximum CFL')
    call d%mfile%add_column(d%fs%Umax,'Umax')
    call d%mfile%add_column(d%meanU,'Umean')
    call d%mfile%add_column(d%fs%Vmax,'Vmax')
    call d%mfile%add_column(d%meanV,'Vmean')
    call d%mfile%add_column(d%fs%Wmax,'Wmax')
    call d%mfile%add_column(d%meanW,'Wmean')
    call d%mfile%add_column(d%fs%Pmax,'Pmax')
    call d%mfile%add_column(d%fs%divmax,'Maximum divergence')
    call d%mfile%add_column(d%fs%psolv%it,'Pressure iteration')
    call d%mfile%add_column(d%fs%psolv%rerr,'Pressure error')
    call d%mfile%write()

    ! Create CFL monitor
    d%cflfile = monitor(d%fs%cfg%amRoot, trim(d%desc) // '_cfl')
    call d%cflfile%add_column(time%n,'Timestep number')
    call d%cflfile%add_column(time%t,'Time')
    call d%cflfile%add_column(d%fs%CFLc_x,'Convective xCFL')
    call d%cflfile%add_column(d%fs%CFLc_y,'Convective yCFL')
    call d%cflfile%add_column(d%fs%CFLc_z,'Convective zCFL')
    call d%cflfile%add_column(d%fs%CFLv_x,'Viscous xCFL')
    call d%cflfile%add_column(d%fs%CFLv_y,'Viscous yCFL')
    call d%cflfile%add_column(d%fs%CFLv_z,'Viscous zCFL')
    call d%cflfile%write()

    ! Create LPT monitor
    call d%lp%get_max()
    d%lptfile = monitor(d%lp%cfg%amRoot, trim(d%desc) // '_lpt')
    call d%lptfile%add_column(time%n,'Timestep number')
    call d%lptfile%add_column(time%t,'Time')
    call d%lptfile%add_column(d%lp%np,'Particle number')
    call d%lptfile%add_column(d%lp%Umin,'Particle Umin')
    call d%lptfile%add_column(d%lp%Umax,'Particle Umax')
    call d%lptfile%add_column(d%lp%Vmin,'Particle Vmin')
    call d%lptfile%add_column(d%lp%Vmax,'Particle Vmax')
    call d%lptfile%add_column(d%lp%Wmin,'Particle Wmin')
    call d%lptfile%add_column(d%lp%Wmax,'Particle Wmax')
    call d%lptfile%add_column(d%lp%dmin,'Particle dmin')
    call d%lptfile%add_column(d%lp%dmax,'Particle dmax')
    call d%lptfile%write()

  end subroutine setup_gendom_monitors

  subroutine setup_hitcube_monitors(d)
    implicit none
    class(cubedomain), intent(inout) :: d

    ! Create hit monitor
    d%hitfile = monitor(d%fs%cfg%amRoot, trim(d%desc) // '_hit')
    call d%hitfile%add_column(time%n, 'Timestep number')
    call d%hitfile%add_column(time%t, 'Time')
    call d%hitfile%add_column(d%Re_lambda, 'Re_lambda')
    call d%hitfile%add_column(d%TKE_ratio, 'TKEratio')
    call d%hitfile%add_column(d%EPS_ratio, 'EPSratio')
    call d%hitfile%add_column(d%ell_ratio, 'ell_ratio')
    call d%hitfile%add_column(d%dx_eta, 'dx_eta')
    call d%hitfile%add_column(d%eta, 'eta')
    call d%hitfile%add_column(d%TKE, 'TKE')
    call d%hitfile%add_column(d%EPS, 'EPS')
    call d%hitfile%add_column(d%urms, 'Urms')
    call d%hitfile%write()

  end subroutine setup_hitcube_monitors

  subroutine write_gendom_monitors(d)
    implicit none
    class(gendomain), intent(inout) :: d

    call d%fs%get_max()
    call d%lp%get_max()
    call d%mfile%write()
    call d%cflfile%write()
    call d%lptfile%write()

  end subroutine write_gendom_monitors

  subroutine write_hitcube_monitors(d)
    implicit none
    class(cubedomain), intent(inout) :: d

    call d%hitfile%write()

  end subroutine write_hitcube_monitors


  !> Initialization of each domain
  !TODO fix geometry dims, overlap, etc

  subroutine geometry_cube_init()
    use parallel, only: group
    implicit none
    type(sgrid) :: grid
    integer :: i, Nx
    real(WP) :: Lx
    integer, dimension(3) :: partition
    real(WP), dimension(:), allocatable :: x, y

    call param_read('Cube Lx', Lx)
    call param_read('Cube Nx', Nx)
    allocate(x(Nx+1), y(Nx+1))
    y(:) = (/ (real(i, WP) / Nx * Lx, i = 0, Nx) /)
    x(:) = y(:) - y(Nx+1) + Lx * GRID_OVERLAP
    grid = sgrid(coord=cartesian, no=2, x=x, y=y, z=y, xper=.true.,           &
      yper=.true., zper=.true., name='cube')
    deallocate(x, y)

    call param_read('Partition', partition, short='pc')
    cube%cfg = config(grp=group, decomp=partition, grid=grid)
    cube%cfg%VF = 1.0_WP

    cube%desc = 'cube'

    call param_read('Particle reset boundaries', part_reset_bdrys)
    part_reset_bdrys(:) = part_reset_bdrys(:) + cube%cfg%x(cube%cfg%imin)
    if (amroot) then
      write(*,*) "DEBUG Cube boundaries: ", cube%cfg%x(cube%cfg%imin), ", ", cube%cfg%x(cube%cfg%imax+1)
      write(*,*) "DEBUG Reset threshholds: ", part_reset_bdrys
    end if

  end subroutine geometry_cube_init

  subroutine geometry_xflow_init()
    use parallel, only: group
    implicit none
    type(sgrid) :: grid
    integer :: i, Nx, Ny
    integer, dimension(3) :: partition
    real(WP) :: Lx, Ly
    real(WP), dimension(:), allocatable :: x, y

    call param_read('Xflow Ll', Lx)
    call param_read('Xflow Nl', Nx)
    call param_read('Cube Lx',  Ly)
    call param_read('Xflow Nw', Ny)
    allocate(x(Nx+1), y(Ny+1))
    x(:) = (/ (real(i, WP) / Nx * Lx, i = 0, Nx) /)
    y(:) = (/ (real(i, WP) / Ny * Ly, i = 0, Ny) /)
    grid = sgrid(coord=cartesian, no=2, x=x, y=y, z=y, xper=.false.,          &
      yper=.true., zper=.true., name='xflow')
    deallocate(x, y)

    call param_read('Partition', partition, short='pp')
    xflw%cfg = config(grp=group, decomp=partition, grid=grid)
    xflw%cfg%VF = 1.0_WP

    xflw%desc = 'xflw'

  end subroutine geometry_xflow_init

  function xflow_inflow_locator(pg, i, j, k) result(flag)
    use pgrid_class, only: pgrid
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i, j, k
    logical :: flag

    flag = i .eq. pg%imin

  end function xflow_inflow_locator

  function xflow_outflow_locator(pg, i, j, k) result(flag)
    use pgrid_class, only: pgrid
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i, j, k
    logical :: flag

    flag = i .eq. pg%imax_+1

  end function xflow_outflow_locator

  !TODO move this to xflowstats
  !subroutine geometry_xflowslice_init()
  !  use parallel, only: group
  !  implicit none
  !  type(sgrid) :: grid
  !  integer :: i, Nx, Ny
  !  integer, dimension(3) :: partition
  !  real(WP) :: obs_start, obs_end
  !  real(WP), dimension(:), allocatable :: x, y, z

  !  call param_read('Obs start', obs_start)
  !  call param_read('Obs end',   obs_end  )
  !  call param_read('Obs Nl', Nx)
  !  call param_read('Obs Nw', Ny)
  !  allocate(x(Nx+1), y(Ny+1), z(Ny+1))
  !  x(:) = (/ (real(i, WP) / Nx, i = 0, Nx) /)
  !  x(:) = x(:) * (obs_end - obs_start) + obs_start
  !  y(:) = (/ (real(i, WP) / Ny, i = 0, Ny) /)
  !  z(:) = (/ (real(i, WP) / Ny, i = 0, Ny) /)
  !  y(:) = y(:) * (xflw%cfg%y(xflw%cfg%jmax+1) - xflw%cfg%y(xflw%cfg%jmin))   &
  !    + xflw%cfg%y(xflw%cfg%jmin)
  !  z(:) = z(:) * (xflw%cfg%z(xflw%cfg%kmax+1) - xflw%cfg%z(xflw%cfg%kmin))   &
  !    + xflw%cfg%z(xflw%cfg%kmin)
  !  grid = sgrid(coord=cartesian, no=1, x=x, y=y, z=z, xper=.false.,          &
  !    yper=.true., zper=.true., name='HIT')
  !  deallocate(x,y,z)

  !  call param_read('Partition', partition)
  !  xflwslice%cfg = config(grp=group, decomp=partition, grid=grid)
  !  xflwslice%cfg%VF = 1.0_WP

  !end subroutine geometry_xflowslice_init

  subroutine gendomain_allocate(d)
    implicit none
    class(gendomain), intent(inout) :: d
    integer :: imino_, imaxo_, jmino_, jmaxo_, kmino_, kmaxo_

    imino_ = d%cfg%imino_; imaxo_ = d%cfg%imaxo_;
    jmino_ = d%cfg%jmino_; jmaxo_ = d%cfg%jmaxo_;
    kmino_ = d%cfg%kmino_; kmaxo_ = d%cfg%kmaxo_;

    allocate(d%resU(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%resV(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%resW(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%Ui(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%Vi(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%Wi(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%SR(1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%gradu(1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
    allocate(d%rho(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  end subroutine gendomain_allocate


  !> Standard simulation functions

  subroutine simulation_init
    use param, only: param_read
    use messager, only: log
    implicit none

    ! Build geometries
    call geometry_cube_init()
    call geometry_xflow_init()
    !call geometry_xflowslice_init()
    if (amroot) call log("Initialized geometries.")

    ! Before allocating anything, create directories needed for ensight and hitstats output
    !TODO see if this needs to be done / can be done in a better way; update directories either way
    !TODO update these directory paths
    !create_directories: block
    !  use parallel, only: comm
    !  use mpi_f08, only: mpi_barrier
    !  if (amroot) then
    !    call execute_command_line("mkdir -p ./npy")
    !    call execute_command_line("mkdir -p ./npy/hitslicexz")
    !    call execute_command_line("mkdir -p ./npy/hitslicexy")
    !    call execute_command_line("mkdir -p ./ensight")
    !    call execute_command_line("mkdir -p ./ensight/HIT")
    !    call execute_command_line("mkdir -p ./ensight/HIT/velocity")
    !    call execute_command_line("mkdir -p ./ensight/HIT/pressure")
    !    call execute_command_line("mkdir -p ./ensight/HIT/particles")
    !  end if
    !  call mpi_barrier(comm)
    !end block create_directories

    ! Allocate domain memory
    call gendomain_allocate(cube)
    call gendomain_allocate(xflw)
    if (amroot) call log("Allocated domain memory.")

    ! Initialize time tracker with 2 subiterations
    initialize_timetracker: block
      time = timetracker(amRoot=amroot)
      call param_read('Max timestep size', time%dtmax)
      call param_read('Max cfl number', time%cflmax)
      time%dtmin = 1e3 * epsilon(1.0_WP)
      time%dt = time%dtmax
      time%itmax = 2
    end block initialize_timetracker
    if (amroot) call log("Initialized timetracker.")

    ! Initialize EC and get parameters
    initialize_ec: block
      use messager, only: die
      real(WP) :: interval

      ec = eccontroller_mesh(cube%cfg)
      call param_read('Number of particles', ec%Np)
      call param_read('EC integral timescales',ec%interval_tinfs)

      ec_evt = threshold_event(time=time, name='EC output')
      call ec%get_next_params(ec_params, ec_done)
      if (ec_done) call die("[EC] It was over before it started.")

      call ec%get_interval(interval)
      ec_evt%tnext = interval
      ec_evt%occurred = .false.
      ! don't call update_parameters here; it needs to be called after the
      ! flow solver and lpt are set up

      call ec%monitor_setup()

    end block initialize_ec
    if (amroot) call log("Initialized EC.")

    ! Read xflow parameters
    call param_read('Xflow rms mult', xflowmult)
    !call param_read('Xflow relax time', xflowrelaxtime)
    xflowvel = sqrt(2.0_WP / 3 * ec_params(3)) * xflowmult

    ! Create a single-phase flow solvers
    create_and_initialize_flow_solver: block
      use hypre_str_class, only: pcg_pfmg
      use incomp_class,    only: dirichlet, clipped_neumann

      ! Cube
      cube%fs = incomp(cfg=cube%cfg,name='NS solver')
      cube%fs%visc(:,:,:) = ec_params(1) * ec_params(5)
      cube%ps = fft3d(cfg=cube%cfg, name='Pressure', nst=7)
      call cube%fs%setup(pressure_solver=cube%ps)
      cube%fs%rho = ec_params(1)
      cube%rho(:,:,:) = ec_params(1)

      ! XFlow
      xflw%fs = incomp(cfg=xflw%cfg, name='NS solver')
      xflw%fs%visc(:,:,:) = ec_params(1) * ec_params(5)
      xflw%fs%rho = ec_params(1)
      xflw%rho(:,:,:) = ec_params(1)
      call xflw%fs%add_bcond(name='inflow', type=dirichlet, face='x',         &
        dir=-1, canCorrect=.false., locator=xflow_inflow_locator)
      call xflw%fs%add_bcond(name='outflow', type=clipped_neumann, face='x',  &
        dir=+1, canCorrect=.true., locator=xflow_outflow_locator)
      xflw%ps = fft2d(cfg=xflw%cfg, name='Pressure', nst=7)
      call xflw%fs%setup(pressure_solver=xflw%ps)

    end block create_and_initialize_flow_solver
    if (amroot) call log("Initialized flow solver.")

    ! Initialize LPT in cube
    initialize_lpt_cube: block
      use random, only: random_uniform, random_normal
      integer :: i, Np
      integer, dimension(3) :: idx0, idx
      real(WP) :: xl, xh, yl, yh, zl, zh

      ! Create solver
      cube%lp = lpt(cfg=cube%cfg, name='Cube LPT')

      ! Get drag model from the input
      call param_read('Drag model', cube%lp%drag_model, default='Schiller-Naumann')

      ! Get number of particles
      call param_read('Number of particles', Np)

      ! set filter width
      cube%lp%filter_width = FILTER_MESH_RATIO * cube%cfg%min_meshsize

      ! set density
      cube%lp%rho = ec_params(2)

      ! Root process initializes Np uniformly
      xl = cube%lp%cfg%x(cube%lp%cfg%imin)
      xh = cube%lp%cfg%x(cube%lp%cfg%imax+1)
      yl = cube%lp%cfg%y(cube%lp%cfg%jmin)
      yh = cube%lp%cfg%y(cube%lp%cfg%jmax+1)
      zl = cube%lp%cfg%z(cube%lp%cfg%kmin)
      zh = cube%lp%cfg%z(cube%lp%cfg%kmax+1)
      idx0(:) = (/ cube%lp%cfg%imin, cube%lp%cfg%jmin, cube%lp%cfg%kmin /)

      if (cube%lp%cfg%amRoot) then
        call cube%lp%resize(Np)
        do i = 1, Np
          ! Give id
          cube%lp%p(i)%id = int(i, 8)
          ! Set the diameter
          cube%lp%p(i)%d = ec_params(6)
          ! Assign random position in the domain
          cube%lp%p(i)%pos = (/ random_uniform(xl, xh),                       &
            random_uniform(yl, yh), random_uniform(zl, zh) /)
          ! Locate the particle on the mesh
          cube%lp%p(i)%ind = cube%lp%cfg%get_ijk_global(cube%lp%p(i)%pos, idx0)
          ! Give equilibrium velocity
          idx(:) = cube%lp%p(i)%ind(:)
          cube%lp%p(i)%vel(:) = cube%cfg%get_velocity(pos=cube%lp%p(i)%pos,   &
            i0=idx(1), j0=idx(2), k0=idx(3), U=cube%fs%U, V=cube%fs%V,        &
            W=cube%fs%W)
          ! Give zero velocity
          !cube%lp%p(i)%vel(:) = 0.0_WP
          ! Activate the particle
          cube%lp%p(i)%flag = 0
        end do
      else
        call cube%lp%resize(0)
      end if

      ! Distribute particles
      call cube%lp%sync()

      ! Get initial particle volume fraction
      call cube%lp%update_VF()

    end block initialize_lpt_cube
    if (amroot) call log("Initialized cube LPT.")

    ! Initialize LPT in xflow
    ! we put one particle toward the end so that the first ensight save isn't
    ! empty
    initialize_lpt_xflow: block

      xflw%lp = lpt(cfg=xflw%cfg, name='Xflow LPT')

      xflw%lp%drag_model = cube%lp%drag_model
      xflw%lp%filter_width = cube%lp%filter_width
      xflw%lp%rho = cube%lp%rho

      if (xflw%lp%cfg%amRoot) then
        call xflw%lp%resize(1)
        xflw%lp%p(1)%id = int(666, 8)
        xflw%lp%p(1)%d = ec_params(6)
        xflw%lp%p(1)%pos = (/ 0.8_WP * xflw%cfg%xL, 0.3_WP * xflw%cfg%yL,       &
          0.3_WP * xflw%cfg%zL /)
        xflw%lp%p(1)%ind = (/ xflw%cfg%imax-2, xflw%cfg%jmin+1, xflw%cfg%kmin+1 /)
        xflw%lp%p(1)%ind = xflw%lp%cfg%get_ijk_global(xflw%lp%p(1)%pos,         &
          xflw%lp%p(1)%ind)
        xflw%lp%p(1)%vel(:) = xflw%cfg%get_velocity(pos=xflw%lp%p(1)%pos,       &
          i0=xflw%lp%p(1)%ind(1), j0=xflw%lp%p(1)%ind(2),                       &
          k0=xflw%lp%p(1)%ind(3), U=xflw%fs%U, V=xflw%fs%V, W=xflw%fs%W)
        xflw%lp%p(1)%flag = 3
      else
        call xflw%lp%resize(0)
      end if
      call xflw%lp%sync()
      call xflw%lp%update_VF()

    end block initialize_lpt_xflow
    if (amroot) call log("Initialized xflow LPT.")

    ! Set up cubestats and xflwstats
    !TODO re-enable
    !initialize_statsobjs: block
    !  character(len=str_medium) :: filterfile
    !  integer, dimension(3) :: hitFFTN, xflowFFTN
    !  call param_read('Filter list', filterfile)
    !  call param_read('Cube FFT mesh', hitFFTN)
    !  call param_read('XFlow FFT mesh', xflowFFTN)
    !  call cube%stats%init(cube%cfg, filterfile, FFTN, cube%rho,              &
    !    cube%fs%visc, cube%fs%U, cube%fs%V, cube%fs%W, cube%ps)
    !  call xflow%stats%init(cube%cfg, filterfile, FFTN, cube%rho,             &
    !    cube%fs%visc, cube%fs%U, cube%fs%V, cube%fs%W, cube%ps)
    !  call cube%stats%init_filters()
    !  call xflow%stats%init_filters()
    !  call cube%stats%setup_sliceio()
    !  call xflow%stats%setup_sliceio()
    !end block initialize_statsobjs
    if (amroot) call log("Initialized statistics.")

    ! Initialize forcing
    initialize_forcing: block
      call param_read('Forcing constant', G)
    end block initialize_forcing

    ! Set up domain couplers
    dom_cpl = coupler(src_grp=cube%cfg%group, dst_grp=xflw%cfg%group,         &
      name='Domain Coupler')
    dom_lptcpl = lptcoupler(src_grp=cube%cfg%group, dst_grp=xflw%cfg%group,   &
      name='Domain LPT Coupler')
    dom_lptcpl%srcflag = 2
    dom_lptcpl%dstflag = 0
    allocate(dom_lptcpl%pushflagfilter(1))
    dom_lptcpl%pushflagfilter(1) = 2
    call dom_cpl%set_src(cube%cfg)
    call dom_lptcpl%set_src(cube%lp)
    if (.true.) then    ! all tasks are in the destination group
      call dom_cpl%set_dst(xflw%cfg)
      call dom_lptcpl%set_dst(xflw%lp)
    end if
    call dom_cpl%initialize()
    call dom_lptcpl%initialize()
    if (amroot) call log("Initialized domain couplers.")

    ! Prepare initial velocity fields
    call init_vel(cube%fs, ec_params(3), cube%resU, cube%resV, cube%resW)
    !call init_vel(xflw%fs, ec_params(3), xflw%resU, xflw%resV, xflw%resW,     &
    !  mask=dom_cpl%overlap)
    xflw%fs%U(:,:,:) = 0.0_WP
    xflw%fs%V(:,:,:) = 0.0_WP
    xflw%fs%W(:,:,:) = 0.0_WP
    if (amroot) call log("Initial velocity fields prepared.")

    ! Calculate divergence
    call cube%fs%get_div()

    ! update parameters to print logging information
    call update_parameters()

    ! Compute initial turbulence stats
    call compute_cube_stats()

    ! Set up ensight output
    call setup_ens(cube)
    if (amroot) call log("Initialized ensight for cube.")
    call setup_ens(xflw)
    if (amroot) call log("Initialized ensight for xflow.")

    ! Create event for Ensight output
    ens_evt = periodic_event(time=time, name='Ensight output')
    call param_read('Ensight at intervals', ens_at_ints)
    if (.not. ens_at_ints) then
      call param_read('Ensight output period',ens_evt%tper)
    else
      ens_evt%tper = 6e66_WP
    end if

    ! Initialize timers
    wt_total%time = 0.0_WP    ; wt_total%percent = 0.0_WP;
    wt_cube_vel%time = 0.0_WP ; wt_cube_vel%percent = 0.0_WP;
    wt_xflw_vel%time = 0.0_WP ; wt_xflw_vel%percent = 0.0_WP;
    wt_cube_pres%time = 0.0_WP; wt_cube_pres%percent = 0.0_WP;
    wt_xflw_pres%time = 0.0_WP; wt_xflw_pres%percent = 0.0_WP;
    wt_cube_lpt%time = 0.0_WP ; wt_cube_lpt%percent = 0.0_WP;
    wt_xflw_lpt%time = 0.0_WP ; wt_xflw_lpt%percent = 0.0_WP;
    wt_stat%time = 0.0_WP     ; wt_stat%percent = 0.0_WP;
    wt_force%time = 0.0_WP    ; wt_force%percent = 0.0_WP;
    wt_cpl%time = 0.0_WP      ; wt_cpl%percent = 0.0_WP;
    wt_rest%time = 0.0_WP     ; wt_rest%percent = 0.0_WP;


    ! Create a monitor file
    create_monitor: block

      ! Prepare some info about fields
      call cube%fs%get_cfl(time%dt,time%cfl)
      call cube%fs%get_max()
      call xflw%fs%get_cfl(time%dt,time%cfl)
      call xflw%fs%get_max()

      ! Initialize monitor files
      call setup_gendom_monitors(cube)
      call setup_hitcube_monitors(cube)
      call setup_gendom_monitors(xflw)
      tfile = monitor(amroot=amroot, name='timing')
      call tfile%add_column(time%n, 'Timestep number')
      call tfile%add_column(time%t, 'Time')
      call tfile%add_column(wt_total%time, 'Total [s]')
      call tfile%add_column(wt_cube_vel%time, 'Cube Vel [s]')
      call tfile%add_column(wt_cube_vel%percent, 'Cube Vel [%]')
      call tfile%add_column(wt_cube_pres%time, 'Cube Pres [s]')
      call tfile%add_column(wt_cube_pres%percent, 'Cube Pres [%]')
      call tfile%add_column(wt_cube_lpt%time, 'Cube LPT [s]')
      call tfile%add_column(wt_cube_lpt%percent, 'Cube LPT [%]')
      call tfile%add_column(wt_xflw_vel%time, 'Xflow Vel [s]')
      call tfile%add_column(wt_xflw_vel%percent, 'Xflow Vel [%]')
      call tfile%add_column(wt_xflw_pres%time, 'Xflow Pres [s]')
      call tfile%add_column(wt_xflw_pres%percent, 'Xflow Pres [%]')
      call tfile%add_column(wt_xflw_lpt%time, 'Xflow LPT [s]')
      call tfile%add_column(wt_xflw_lpt%percent, 'Xflow LPT [%]')
      call tfile%add_column(wt_stat%time, 'Stats [s]')
      call tfile%add_column(wt_stat%percent, 'Stats [%]')
      call tfile%add_column(wt_force%time, 'Forcing [s]')
      call tfile%add_column(wt_force%percent, 'Forcing [%]')
      call tfile%add_column(wt_cpl%time, 'Coupling [s]')
      call tfile%add_column(wt_cpl%percent, 'Coupling [%]')
      call tfile%add_column(wt_rest%time, 'Rest [s]')
      call tfile%add_column(wt_rest%percent, 'Rest [%]')
      call tfile%write()

    end block create_monitor
    if (amroot) call log("Primary monitor files created.")

    ! allocate filter memory
    !TODO

    ! write zero time statistics
    !TODO
    !call ec%compute_statistics(Re_lambda, Stk, phiinf, Wovk, urms, eta, time%t, time%n)
    !call ec%io_write(time%t)

    ! Output initial ensight
    if (.not. ens_at_ints .and. ens_evt%occurs()) then
      call write_ens(cube)
      call write_ens(xflw)
    end if

    ! Announce success
    if (amroot) call log("Finished simulation initialization.")

  end subroutine simulation_init

  !> Time integrate our problem
  subroutine simulation_run
    use messager,            only: die
    use parallel,            only: parallel_time
    use eccontroller_class,  only: FORCE_TIMESCALE
    implicit none
    integer :: it

    ! Perform time integration
    do while (.not. (time%done() .or. ec_done))

      do while (.not. ec_evt%occurs())

        !TODO subiterate the cube and the crossflow separately
        !TODO the right way to do this is probably just to write a step routine that take domain objects (and has its own
        !subiteration counter), then call it for both objects

        !! Wallclock
        wt_total%time_in=parallel_time()

        !! Increment time
        call cube%fs%get_cfl(time%dt, cube%cfl)
        call xflw%fs%get_cfl(time%dt, xflw%cfl)
        time%cfl = max(cube%cfl, xflw%cfl)
        call time%adjust_dt()
        time%dt = min(time%dt, FORCE_TIMESCALE / G)
        !time%dt = min(time%dt, sqrt(nu / max(cube%EPS, cube%EPSp)))
        time%dt = min(time%dt, sqrt(nu / cube%EPS))
        !TODO re-enable when controller is re-integrated
        !time%dt = min(time%dt, ec_evt%tnext - time%t)
        call time%increment()

        !! Cube step
        ! Remember old velocity
        cube%fs%Uold = cube%fs%U
        cube%fs%Vold = cube%fs%V
        cube%fs%Wold = cube%fs%W
        ! advance particles
        wt_cube_lpt%time_in = parallel_time()
        !call cube%lp%collide(dt=time%dtmid)
        call cube%lp%advance(dt=time%dtmid, U=cube%fs%U, V=cube%fs%V,         &
          W=cube%fs%W, rho=cube%rho, visc=cube%fs%visc)
        wt_cube_lpt%time = wt_cube_lpt%time + parallel_time() - wt_cube_lpt%time_in
        ! re-flag particles within the acceptable range as copyable
        reset_particle_flags: block
          integer :: i
          do i = 1, cube%lp%np_
            if (part_reset_bdrys(1) .lt. cube%lp%p(i)%pos(1) .and.            &
                part_reset_bdrys(2) .gt. cube%lp%p(i)%pos(1) .and.            &
                cube%lp%p(i)%flag .eq. 2) then
              cube%lp%p(i)%flag = 0
            end if
          end do
        end block reset_particle_flags
        ! Perform sub-iterations
        do it = 1, time%itmax
          ! Velocity step
          wt_cube_vel%time_in = parallel_time()
          cube%fs%U = 0.5_WP * (cube%fs%U + cube%fs%Uold)
          cube%fs%V = 0.5_WP * (cube%fs%V + cube%fs%Vold)
          cube%fs%W = 0.5_WP * (cube%fs%W + cube%fs%Wold)
          call cube%fs%get_dmomdt(cube%resU, cube%resV, cube%resW)
          cube%resU = -2.0_WP*(cube%fs%rho*cube%fs%U-cube%fs%rho*cube%fs%Uold)+time%dt*cube%resU
          cube%resV = -2.0_WP*(cube%fs%rho*cube%fs%V-cube%fs%rho*cube%fs%Vold)+time%dt*cube%resV
          cube%resW = -2.0_WP*(cube%fs%rho*cube%fs%W-cube%fs%rho*cube%fs%Wold)+time%dt*cube%resW
          wt_cube_vel%time = wt_cube_vel%time + parallel_time() - wt_cube_vel%time_in
          ! Linear forcing term (Bassenne et al. 2016)
          wt_force%time_in = parallel_time()
          ! note we have stats from the previous step; they haven't changed
          !TODO DEBUG --- recompute statistics to see if this is the issue
          call compute_cube_stats()
          linear_forcing: block
            real(WP) :: A
            ! - Eq. (7) (forcing constant TKE)
            A = (cube%EPSp - (G / FORCE_TIMESCALE) * (cube%TKE -              &
              TKE_target)) / (2.0_WP * cube%TKE) * cube%fs%rho
            !A = (cube%EPS - (G / FORCE_TIMESCALE) * (cube%TKE -               &
            !  TKE_target)) / (2.0_WP * cube%TKE) * cube%fs%rho
            ! update residuals
            cube%resU = cube%resU + time%dt * A * (cube%fs%U - cube%meanU)
            cube%resV = cube%resV + time%dt * A * (cube%fs%V - cube%meanV)
            cube%resW = cube%resW + time%dt * A * (cube%fs%W - cube%meanW)
          end block linear_forcing
          wt_force%time = wt_force%time + parallel_time() - wt_force%time_in
          ! Apply residuals
          wt_cube_vel%time_in = parallel_time()
          cube%fs%U = 2.0_WP * cube%fs%U - cube%fs%Uold + cube%resU / cube%fs%rho
          cube%fs%V = 2.0_WP * cube%fs%V - cube%fs%Vold + cube%resV / cube%fs%rho
          cube%fs%W = 2.0_WP * cube%fs%W - cube%fs%Wold + cube%resW / cube%fs%rho
          wt_cube_vel%time = wt_cube_vel%time + parallel_time() - wt_cube_vel%time_in
          ! Solve Poisson equation, correct velocity
          wt_cube_pres%time_in = parallel_time()
          call cube%fs%correct_mfr()
          call cube%fs%get_div()
          cube%fs%psolv%rhs = -cube%fs%cfg%vol*cube%fs%div*cube%fs%rho/time%dt
          call cube%fs%psolv%solve()
          call cube%fs%shift_p(cube%fs%psolv%sol)
          call cube%fs%get_pgrad(cube%fs%psolv%sol, cube%resU, cube%resV, cube%resW)
          cube%fs%P = cube%fs%P + cube%fs%psolv%sol
          cube%fs%U = cube%fs%U - time%dt * cube%resU / cube%fs%rho
          cube%fs%V = cube%fs%V - time%dt * cube%resV / cube%fs%rho
          cube%fs%W = cube%fs%W - time%dt * cube%resW / cube%fs%rho
          wt_cube_pres%time = wt_cube_pres%time + parallel_time() - wt_cube_pres%time_in
        end do

        !! Apply xflow bcs
        ! update inflow bc using right side of cube 
        ! here we use resU, resV, resW as pieces of memory that share a
        ! compatible layout and pgrid with both geometries and contain
        ! information at the inflow boundary (because it's a staggered mesh)
        ! the overlap region needs to be as small as possible (for
        ! performance) while being large enough to not cause any funky
        ! edge cases or interpolation issues
        call dom_cpl%push(cube%fs%U)
        call dom_cpl%transfer()
        call dom_cpl%pull(xflw%resU)
        call dom_cpl%push(cube%fs%V)
        call dom_cpl%transfer()
        call dom_cpl%pull(xflw%resV)
        call dom_cpl%push(cube%fs%W)
        call dom_cpl%transfer()
        call dom_cpl%pull(xflw%resW)
        call dom_lptcpl%push()
        call dom_lptcpl%transfer()
        call dom_lptcpl%pull()
        apply_inflow_vals: block
          if (xflw%cfg%imin .eq. xflw%cfg%imin_) then
            xflw%fs%U(xflw%cfg%imino_:xflw%cfg%imin_,:,:) = xflw%resU(xflw%cfg%imino_:xflw%cfg%imin_,:,:)
            xflw%fs%V(xflw%cfg%imino_:xflw%cfg%imin_-1,:,:) = xflw%resV(xflw%cfg%imino_:xflw%cfg%imin_-1,:,:)
            xflw%fs%W(xflw%cfg%imino_:xflw%cfg%imin_-1,:,:) = xflw%resW(xflw%cfg%imino_:xflw%cfg%imin_-1,:,:)
          end if
        end block apply_inflow_vals
        call xflw%fs%apply_bcond(time%t, time%dt)

        !! Xflw step
        ! Remember old velocity
        xflw%fs%Uold = xflw%fs%U
        xflw%fs%Vold = xflw%fs%V
        xflw%fs%Wold = xflw%fs%W
        ! advance particles
        wt_xflw_lpt%time_in = parallel_time()
        !call xflw%lp%collide(dt=time%dtmid)
        call xflw%lp%advance(dt=time%dtmid, U=xflw%fs%U, V=xflw%fs%V,    &
          W=xflw%fs%W, rho=xflw%rho, visc=xflw%fs%visc)
        wt_xflw_lpt%time = wt_xflw_lpt%time + parallel_time() - wt_xflw_lpt%time_in
        ! Perform sub-iterations
        do it = 1, time%itmax
          ! Velocity step
          wt_xflw_vel%time_in = parallel_time()
          xflw%fs%U = 0.5_WP * (xflw%fs%U + xflw%fs%Uold)
          xflw%fs%V = 0.5_WP * (xflw%fs%V + xflw%fs%Vold)
          xflw%fs%W = 0.5_WP * (xflw%fs%W + xflw%fs%Wold)
          call xflw%fs%get_dmomdt(xflw%resU, xflw%resV, xflw%resW)
          xflw%resU = -2.0_WP*(xflw%fs%rho*xflw%fs%U-xflw%fs%rho*xflw%fs%Uold)+time%dt*xflw%resU
          xflw%resV = -2.0_WP*(xflw%fs%rho*xflw%fs%V-xflw%fs%rho*xflw%fs%Vold)+time%dt*xflw%resV
          xflw%resW = -2.0_WP*(xflw%fs%rho*xflw%fs%W-xflw%fs%rho*xflw%fs%Wold)+time%dt*xflw%resW
          wt_xflw_vel%time = wt_xflw_vel%time + parallel_time() - wt_xflw_vel%time_in
          ! Apply residuals
          wt_xflw_vel%time_in = parallel_time()
          xflw%fs%U = 2.0_WP * xflw%fs%U - xflw%fs%Uold + xflw%resU / xflw%fs%rho
          xflw%fs%V = 2.0_WP * xflw%fs%V - xflw%fs%Vold + xflw%resV / xflw%fs%rho
          xflw%fs%W = 2.0_WP * xflw%fs%W - xflw%fs%Wold + xflw%resW / xflw%fs%rho
          call xflw%fs%apply_bcond(time%t, time%dt)
          wt_xflw_vel%time = wt_xflw_vel%time + parallel_time() - wt_xflw_vel%time_in
          ! Solve Poisson equation, correct velocity
          wt_xflw_pres%time_in = parallel_time()
          call xflw%fs%correct_mfr()
          call xflw%fs%get_div()
          xflw%fs%psolv%rhs = -xflw%fs%cfg%vol*xflw%fs%div*xflw%fs%rho/time%dt
          call xflw%fs%psolv%solve()
          call xflw%fs%shift_p(xflw%fs%psolv%sol)
          call xflw%fs%get_pgrad(xflw%fs%psolv%sol, xflw%resU, xflw%resV, xflw%resW)
          xflw%fs%P = xflw%fs%P + xflw%fs%psolv%sol
          xflw%fs%U = xflw%fs%U - time%dt * xflw%resU / xflw%fs%rho
          xflw%fs%V = xflw%fs%V - time%dt * xflw%resV / xflw%fs%rho
          xflw%fs%W = xflw%fs%W - time%dt * xflw%resW / xflw%fs%rho
          wt_xflw_pres%time = wt_xflw_pres%time + parallel_time() - wt_xflw_pres%time_in
        end do

        ! Recompute divergence
        wt_cube_vel%time_in = parallel_time()
        call cube%fs%get_div()
        wt_cube_vel%time = wt_cube_vel%time + parallel_time() - wt_cube_vel%time_in
        wt_xflw_vel%time_in = parallel_time()
        call xflw%fs%get_div()
        wt_xflw_vel%time = wt_xflw_vel%time + parallel_time() - wt_xflw_vel%time_in

        ! Output to ensight
        if (.not. ens_at_ints .and. ens_evt%occurs()) then
          call write_ens(cube)
          call write_ens(xflw)
        end if

        ! Statistics and output monitoring
        wt_stat%time_in = parallel_time()
        call compute_cube_stats()
        call write_gendom_monitors(cube)
        call write_hitcube_monitors(cube)
        wt_stat%time = wt_stat%time + parallel_time() - wt_stat%time_in
        wt_stat%time_in = parallel_time()
        call write_gendom_monitors(xflw)
        wt_stat%time = wt_stat%time + parallel_time() - wt_stat%time_in

        ! Monitor timing
        wt_total%time = parallel_time() - wt_total%time_in
        wt_rest%time = wt_total%time - wt_cube_vel%time - wt_cube_pres%time &
          - wt_cube_lpt%time - wt_xflw_vel%time - wt_xflw_pres%time         &
          - wt_xflw_lpt%time - wt_stat%time - wt_force%time - wt_cpl%time
        wt_total%percent = 100.0_WP
        wt_cube_vel%percent = wt_cube_vel%time / wt_total%time * 100.0_WP
        wt_cube_pres%percent = wt_cube_pres%time / wt_total%time * 100.0_WP
        wt_cube_lpt%percent = wt_cube_lpt%time / wt_total%time * 100.0_WP
        wt_xflw_vel%percent = wt_xflw_vel%time / wt_total%time * 100.0_WP
        wt_xflw_pres%percent = wt_xflw_pres%time / wt_total%time * 100.0_WP
        wt_xflw_lpt%percent = wt_xflw_lpt%time / wt_total%time * 100.0_WP
        wt_stat%percent = wt_stat%time / wt_total%time * 100.0_WP
        wt_force%percent = wt_force%time / wt_total%time * 100.0_WP
        wt_cpl%percent = wt_cpl%time / wt_total%time * 100.0_WP
        wt_rest%percent = wt_rest%time / wt_total%time * 100.0_WP
        call tfile%write()
        wt_total%time = 0.0_WP
        wt_cube_vel%time = 0.0_WP
        wt_cube_pres%time = 0.0_WP
        wt_cube_lpt%time = 0.0_WP
        wt_xflw_vel%time = 0.0_WP
        wt_xflw_pres%time = 0.0_WP
        wt_xflw_lpt%time = 0.0_WP
        wt_stat%time = 0.0_WP
        wt_force%time = 0.0_WP
        wt_cpl%time = 0.0_WP

      end do

      ! move to next parameters
      ec_next: block
        real(WP) :: interval

        !TODO re-enable
        !call cube%stats%compute_stats(time%n)
        !call xflow%stats%compute_stats(time%n)
        call ec%update_write(cube%Re_lambda, cube%Stk, cube%phiinf,           &
          cube%Wovk, cube%urms, cube%eta, time%t, time%n)
        !call cube%stats%write_sliceio(time%t)
        !call xflow%stats%write_sliceio(time%t)

        if (ens_at_ints) then
          call write_ens(cube)
          call write_ens(xflw)
        end if

        call ec%get_next_params(ec_params, ec_done)
        call ec%get_interval(interval)
        ec_evt%tnext = ec_evt%tnext + interval
        ec_evt%occurred = .false.
        call update_parameters()

      end block ec_next

    end do

    if (time%done() .and. .not. ec_done) call die("timer stopped before &
      &parameter sweep was complete")

  end subroutine simulation_run

  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    !TODO

  end subroutine simulation_final


end module simulation

