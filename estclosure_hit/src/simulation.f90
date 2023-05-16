!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use fftsolver3d_class,    only: fftsolver3d
  use incomp_class,         only: incomp
  use timetracker_class,    only: timetracker
  use ensight_class,        only: ensight
  use partmesh_class,       only: partmesh
  use event_class,          only: periodic_event, threshold_event
  use monitor_class,        only: monitor
  use datafile_class,       only: datafile
  use estimclosures_class,  only: estimclosures_mesh
  use lpt_class,            only: lpt
  use string,               only: str_medium
  implicit none
  private

  real(WP), parameter :: FILTER_MESH_RATIO = 3.5_WP

  !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
  type(fftsolver3d), public :: ps
  type(incomp),      public :: fs
  type(timetracker), public :: time

  ! particles
  type(lpt), target, public :: lp
  real(WP), dimension(:,:,:), pointer :: rho, sx, sy, sz

  !> Ensight postprocessing
  type(partmesh)       :: pmesh
  type(ensight)        :: ens_out
  type(periodic_event) :: ens_evt
  logical              :: ens_at_ints

  !> Closure Estimation
  type(estimclosures_mesh) :: ec
  type(threshold_event)    :: ec_evt
  logical                  :: ec_done

  !> Turbulent Parameters
  real(WP), dimension(7)   :: ec_params   ! rhof, rhop, ktarget, epstarget, nu, dp, g
  real(WP) :: TKE_target, EPS_target, nu, dp

  !> Simulation monitor file
  type(monitor) :: mfile, cflfile, hitfile, lptfile, tfile, ssfile

  public :: simulation_init, simulation_run, simulation_final

  !> Private work arrays
  real(WP), dimension(:,:,:), allocatable     :: resU, resV, resW
  real(WP), dimension(:,:,:), allocatable     :: Ui, Vi, Wi
  real(WP), dimension(:,:,:,:), allocatable   :: SR       ! strain rate
  real(WP), dimension(:,:,:,:,:), allocatable :: gradu

  !> Turbulent parameters
  real(WP) :: meanU, meanV, meanW, EPS, EPSp, TKE, urms

  !> Forcing constant (and burn in constant)
  real(WP) :: G

  !> For monitoring
  real(WP) :: Re_lambda, EPS_ratio, TKE_ratio, eta, ell_ratio, dx_eta,        &
    nondimtime, permissible_eps_err, Stk, phiinf, Wovk

  !> Wallclock time for monitoring
  type :: timer; real(WP) :: time_in, time, percent; end type timer
  type(timer) :: wt_total,wt_vel,wt_pres,wt_lpt,wt_rest,wt_stat,wt_force

contains

  !> Compute turbulence stats
  subroutine compute_stats()
    use mathtools, only: PI
    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
    use parallel,  only: MPI_REAL_WP
    real(WP) :: myTKE, myEPS, myEPSp
    integer :: i, j, k, ierr

    ! Compute mean velocities
    call fs%cfg%integrate(A=fs%U,integral=meanU); meanU = meanU / fs%cfg%vol_total;
    call fs%cfg%integrate(A=fs%V,integral=meanV); meanV = meanV / fs%cfg%vol_total;
    call fs%cfg%integrate(A=fs%W,integral=meanW); meanW = meanW / fs%cfg%vol_total;

    ! Compute strainrate and grad(u)
    call fs%get_strainrate(SR=SR)
    call fs%get_gradu(gradu=gradu)

    ! compute TKE, EPS, EPSp, forcing_offset
    myTKE = 0.0_WP; myEPS = 0.0_WP; myEPSp = 0.0_WP;
    do k=fs%cfg%kmin_,fs%cfg%kmax_
      do j=fs%cfg%jmin_,fs%cfg%jmax_
        do i=fs%cfg%imin_,fs%cfg%imax_
          myTKE = myTKE + 0.5_WP*((fs%U(i,j,k)-meanU)**2 +                    &
            (fs%V(i,j,k)-meanV)**2 + (fs%W(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)
          myEPS = myEPS + 2.0_WP*fs%visc(i,j,k)*fs%cfg%vol(i,j,k)*(           &
            SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(             &
            SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))/fs%rho
          myEPSp = myEPSp + fs%cfg%vol(i,j,k)*fs%visc(i,j,k)*(                &
            gradu(1,1,i,j,k)**2 + gradu(1,2,i,j,k)**2 + gradu(1,3,i,j,k)**2 + &
            gradu(2,1,i,j,k)**2 + gradu(2,2,i,j,k)**2 + gradu(2,3,i,j,k)**2 + &
            gradu(3,1,i,j,k)**2 + gradu(3,2,i,j,k)**2 + gradu(3,3,i,j,k)**2)
        end do
      end do
    end do
    call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
    TKE = TKE / fs%cfg%vol_total
    call MPI_ALLREDUCE(myEPS,EPS,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
    EPS = EPS / fs%cfg%vol_total
    call MPI_ALLREDUCE(myEPSp,EPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
    EPSp = EPSp / (fs%cfg%vol_total * fs%rho)

    ! urms
    urms = sqrt(2.0_WP / 3 * TKE)

    ! additional monitoring quantities
    eta = (nu**3 / EPS)**0.25_WP
    Re_lambda = sqrt(20.0_WP / 3) * TKE * (eta / nu)**2
    EPS_ratio = EPS / EPS_target
    TKE_ratio = TKE / TKE_target
    ell_ratio = (2.0_WP / 3 * TKE)**1.5_WP / EPS / (cfg%vol_total)**(1.0_WP/3)
    dx_eta = cfg%max_meshsize / eta
    nondimtime = time%t / (2 * TKE_target) * (3 * EPS_target)
    permissible_eps_err = (TKE - TKE_target) / TKE_target * G * EPS_target

    ! nondimensional particle quantities
    Stk = ec_params(2) / ec_params(1) * ec_params(6)**2 / (18 * eta**2)
    phiinf = lp%np * PI * ec_params(6)**3 / (6 * cfg%vol_total)
    Wovk = ec_params(7) * ec_params(2) / ec_params(1) * ec_params(6)**2 * eta &
      / (18 * nu**2)

  end subroutine compute_stats

  !> params_primit - 7 items - rhof, rhop, ktarget, epstarget, nu, dp, g
  subroutine update_parameters()
    implicit none
    integer :: i

    ! fluid parameters
    fs%rho = ec_params(1)
    rho(:,:,:) = fs%rho
    fs%visc(:,:,:) = ec_params(1) * ec_params(5)
    sx(:,:,:) = 0.0_WP; sy(:,:,:) = 0.0_WP; sz(:,:,:) = 0.0_WP;

    ! particle parameters
    lp%filter_width = FILTER_MESH_RATIO * cfg%min_meshsize
    lp%rho = ec_params(2)
    do i = 1, lp%np_
      lp%p(i)%d = ec_params(6)
    end do
    call lp%sync()

    ! gravity
    lp%gravity(:) = (/ -ec_params(7), 0.0_WP, 0.0_WP /)

    ! store parameters that are difficult to access for reference elsewhere
    TKE_target = ec_params(3); EPS_target = ec_params(4);
    nu = ec_params(5); dp = ec_params(6);

    ! print current parameters
    print_statistics: block
      use string, only: str_long
      use messager, only: log
      character(len=str_long) :: message

      if (cfg%amroot) then
        write(message,'("At t = ",e16.8," updated turbulent parameters to:")') time%t; call log(message);
        write(message,'("    Fluid density: ",e12.6)') fs%rho; call log(message);
        write(message,'("    Target TKE: ",e12.6)') TKE_target; call log(message);
        write(message,'("    Target EPS: ",e12.6)') EPS_target; call log(message);
        write(message,'("    Kinematic Viscosity: ",e12.6)') nu; call log(message);
      end if

    end block print_statistics

  end subroutine update_parameters

  subroutine update_pmesh()
    implicit none
    integer :: i

    call pmesh%reset()
    call pmesh%set_size(lp%np_)
    do i = 1, lp%np_
      pmesh%pos(:,i) = lp%p(i)%pos
      pmesh%var(1,i) = lp%p(i)%id
      pmesh%var(2,i) = lp%p(i)%d
      pmesh%vec(:,1,i) = lp%p(i)%vel
      lp%p(i)%ind = cfg%get_ijk_global(lp%p(i)%pos, lp%p(i)%ind)
      pmesh%vec(:,2,i) = lp%cfg%get_velocity(pos=lp%p(i)%pos,                 &
        i0=lp%p(i)%ind(1), j0=lp%p(i)%ind(2), k0=lp%p(i)%ind(3), U=fs%U,      &
        V=fs%V, W=fs%W)
    end do

  end subroutine update_pmesh

  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none

    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(SR  (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(gradu(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays

    ! Initialize time tracker with 2 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max cfl number',time%cflmax)
      call param_read('Max iter',time%nmax)
      time%dtmin=1e3*epsilon(0.0_WP)
      time%dt=time%dtmax
      time%itmax=2
    end block initialize_timetracker

    ! Initialize timers
    initialize_timers: block
      wt_total%time=0.0_WP; wt_total%percent=0.0_WP;
      wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP;
      wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP;
      wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP;
      wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP;
      wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP;
      wt_force%time=0.0_WP; wt_force%percent=0.0_WP;
    end block initialize_timers

    ! Initialize EC and get parameters
    initialize_ec: block
      use messager, only: die
      real(WP) :: interval

      ec = estimclosures_mesh(cfg)
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

    ! Create a single-phase flow solver without bconds
    create_and_initialize_flow_solver: block
      use hypre_str_class, only: pcg_pfmg

      ! Create flow solver
      fs=incomp(cfg=cfg,name='NS solver')
      ! Assign constant viscosity
      fs%visc(:,:,:) = ec_params(1) * ec_params(5)
      ! Prepare and configure pressure solver
      ps=fftsolver3d(cfg=cfg,name='Pressure',nst=7)
      ! Setup the solver
      call fs%setup(pressure_solver=ps)
      fs%rho = ec_params(1)
      ! Set up constant arrays for lpt
      allocate(rho(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(sx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(sy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(sz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      rho(:,:,:) = ec_params(1)

    end block create_and_initialize_flow_solver

    ! Initialize LPT
    initialize_lpt: block
      use random, only: random_uniform, random_normal
      integer :: i, Np

      ! Create solver
      lp = lpt(cfg=cfg, name='LPT')

      ! Get drag model from the input
      call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')

      ! Get number of particles
      call param_read('Number of particles', Np)

      ! set filter width
      lp%filter_width = FILTER_MESH_RATIO * cfg%min_meshsize

      ! set density
      lp%rho = ec_params(2)

      ! Root process initializes np particles randomly
      if (lp%cfg%amRoot) then
        call lp%resize(np)
        do i=1,np
          ! Give id
          lp%p(i)%id=int(i,8)
          ! Set the diameter
          lp%p(i)%d=ec_params(6)
          ! Assign random position in the domain
          lp%p(i)%pos=(/random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),  &
            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),              &
            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))/)
          ! Locate the particle on the mesh
          lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,(/lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin/))
          ! Give equilibrium velocity
          !lp%p(i)%vel=cfg%get_velocity(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
          ! Give zero velocity
          lp%p(i)%vel(:) = 0.0_WP
          ! Activate the particle
          lp%p(i)%flag=0
        end do
      else
        call lp%resize(0)
      end if

      ! Distribute particles
      call lp%sync()

      ! Get initial particle volume fraction
      call lp%update_VF()

    end block initialize_lpt

    ! Initialize forcing
    initialize_forcing: block
      call param_read('Forcing constant', G)
    end block initialize_forcing

    ! Prepare initial velocity field
    initialize_velocity: block
      use random,    only: random_normal
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use mathtools, only: pi
      real(WP) :: Urms0
      integer :: i,j,k

      Urms0 = sqrt(2.0_WP/3*ec_params(3))
      ! Gaussian initial field
      do k=fs%cfg%kmin_,fs%cfg%kmax_
        do j=fs%cfg%jmin_,fs%cfg%jmax_
          do i=fs%cfg%imin_,fs%cfg%imax_
            fs%U(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
            fs%V(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
            fs%W(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
          end do
        end do
      end do
      call fs%cfg%sync(fs%U)
      call fs%cfg%sync(fs%V)
      call fs%cfg%sync(fs%W)

      ! Compute mean and remove it from the velocity field to obtain <U>=0
      call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total;
      call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total;
      call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total;
      fs%U = fs%U - meanU; fs%V = fs%V - meanV; fs%W = fs%W - meanW;

      ! Project to ensure divergence-free
      call fs%get_div()
      fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
      fs%psolv%sol=0.0_WP
      call fs%psolv%solve()
      call fs%shift_p(fs%psolv%sol)
      call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
      fs%P=fs%P+fs%psolv%sol
      fs%U=fs%U-time%dt*resU/fs%rho
      fs%V=fs%V-time%dt*resV/fs%rho
      fs%W=fs%W-time%dt*resW/fs%rho

      ! Calculate divergence
      call fs%get_div()

      ! update parameters to print logging information
      call update_parameters()

    end block initialize_velocity

    ! Compute initial turbulence stats
    call compute_stats()

    ! initialize hitstats objects now that we have pointers to fs and lp
    call ec%init(lp, rho, fs%visc, fs%U, fs%V, fs%W, sx, sy, sz)

    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      pmesh=partmesh(nvar=2,nvec=2,name='lpt')
      pmesh%varname(1)="id"
      pmesh%varname(2)="dp"
      pmesh%vecname(1)="vel"
      pmesh%vecname(2)="fld"
      call update_pmesh()
    end block create_pmesh

    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='HIT')
      ! Create event for Ensight output
      ens_evt=periodic_event(time=time,name='Ensight output')
      call param_read('Ensight at intervals', ens_at_ints)
      if (.not. ens_at_ints) then
        call param_read('Ensight output period',ens_evt%tper)
      else
        ens_evt%tper = 6e66_WP
      end if
      ! Add variables to output
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('pressure',fs%P)
      call ens_out%add_particle('particles',pmesh)
      ! Set up ensight output for filtered quantities
      call ec%ensight_setup()
      call ec%compute_statistics(Re_lambda, Stk, phiinf, Wovk, urms, ETA, nu, time%t, time%n)
      call ec%ensight_write(time%t)
      ! Output to ensight
      if (.not. ens_at_ints .and. ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight

    ! Create a monitor file
    create_monitor: block
      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_max()
      ! Create simulation monitor
      mfile=monitor(fs%cfg%amRoot,'simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(time%cfl,'Maximum CFL')
      call mfile%add_column(fs%Umax,'Umax')
      call mfile%add_column(meanU,'Umean')
      call mfile%add_column(fs%Vmax,'Vmax')
      call mfile%add_column(meanV,'Vmean')
      call mfile%add_column(fs%Wmax,'Wmax')
      call mfile%add_column(meanW,'Wmean')
      call mfile%add_column(fs%Pmax,'Pmax')
      call mfile%add_column(fs%divmax,'Maximum divergence')
      call mfile%add_column(fs%psolv%it,'Pressure iteration')
      call mfile%add_column(fs%psolv%rerr,'Pressure error')
      call mfile%write()
      ! Create CFL monitor
      cflfile=monitor(fs%cfg%amRoot,'cfl')
      call cflfile%add_column(time%n,'Timestep number')
      call cflfile%add_column(time%t,'Time')
      call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
      call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
      call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
      call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
      call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
      call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
      call cflfile%write()
      ! Create hit monitor
      hitfile=monitor(fs%cfg%amRoot,'hit')
      call hitfile%add_column(time%n,'Timestep number')
      call hitfile%add_column(time%t,'Time')
      call hitfile%add_column(Re_lambda,'Re_lambda')
      call hitfile%add_column(TKE_ratio,'TKEratio')
      call hitfile%add_column(EPS_ratio,'EPSratio')
      call hitfile%add_column(ell_ratio,'ell_ratio')
      call hitfile%add_column(dx_eta,'dx_eta')
      call hitfile%add_column(permissible_eps_err,'permepserr')
      call hitfile%add_column(eta,'eta')
      call hitfile%add_column(TKE,'TKE')
      call hitfile%add_column(EPS,'EPS')
      !call hitfile%add_column(urms,'Urms')
      !call hitfile%add_column(ell,'L')
      call hitfile%write()
      ! Create hit convergence monitor
      ssfile=monitor(fs%cfg%amRoot,'convergence')
      call ssfile%add_column(time%n,'Timestep number')
      call ssfile%add_column(time%t,'Time')
      call ssfile%add_column(nondimtime,'Time/t_int')
      call ssfile%add_column(eps_ratio,'EPS_ratio')
      call ssfile%add_column(tke_ratio,'TKE_ratio')
      call ssfile%add_column(dx_eta,'dx/eta')
      call ssfile%add_column(ell_ratio,'ell_ratio')
      call ssfile%write()
      ! Create LPT monitor
      call lp%get_max()
      lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
      call lptfile%add_column(time%n,'Timestep number')
      call lptfile%add_column(time%t,'Time')
      call lptfile%add_column(lp%np,'Particle number')
      call lptfile%add_column(lp%Umin,'Particle Umin')
      call lptfile%add_column(lp%Umax,'Particle Umax')
      call lptfile%add_column(lp%Vmin,'Particle Vmin')
      call lptfile%add_column(lp%Vmax,'Particle Vmax')
      call lptfile%add_column(lp%Wmin,'Particle Wmin')
      call lptfile%add_column(lp%Wmax,'Particle Wmax')
      call lptfile%add_column(lp%dmin,'Particle dmin')
      call lptfile%add_column(lp%dmax,'Particle dmax')
      call lptfile%write()
      ! Create timing monitor
      tfile=monitor(amroot=fs%cfg%amRoot,name='timing')
      call tfile%add_column(time%n,'Timestep number')
      call tfile%add_column(time%t,'Time')
      call tfile%add_column(wt_total%time,'Total [s]')
      call tfile%add_column(wt_vel%time,'Velocity [s]')
      call tfile%add_column(wt_vel%percent,'Velocity [%]')
      call tfile%add_column(wt_pres%time,'Pressure [s]')
      call tfile%add_column(wt_pres%percent,'Pressure [%]')
      call tfile%add_column(wt_stat%time,'Stats [s]')
      call tfile%add_column(wt_stat%percent,'Stats [%]')
      call tfile%add_column(wt_force%time,'Forcing [s]')
      call tfile%add_column(wt_force%percent,'Forcing [%]')
      call tfile%add_column(wt_rest%time,'Rest [s]')
      call tfile%add_column(wt_rest%percent,'Rest [%]')
      call tfile%write()
    end block create_monitor

  end subroutine simulation_init

  !> Time integrate our problem
  subroutine simulation_run
    use messager,            only: die
    use parallel,            only: parallel_time
    use estimclosures_class, only: FORCE_TIMESCALE
    implicit none

    ! Perform time integration
    do while (.not. (time%done() .or. ec_done))

      do while (.not. ec_evt%occurs())

        ! init wallclock
        wt_total%time_in=parallel_time()

        ! Increment time
        call fs%get_cfl(time%dt,time%cfl)
        call time%adjust_dt()
        time%dt = min(time%dt, FORCE_TIMESCALE / G)
        time%dt = min(time%dt, sqrt(nu / max(EPS, EPSp)))
        time%dt = min(time%dt, ec_evt%tnext - time%t)
        call time%increment()

        ! Remember old velocity
        fs%Uold=fs%U; fs%Vold=fs%V; fs%Wold=fs%W;

        ! advance particles
        wt_lpt%time_in=parallel_time()
        !call lp%collide(dt=time%dtmid)
        call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=rho,visc=fs%visc)
        wt_lpt%time=wt_lpt%time+parallel_time()-wt_lpt%time_in

        ! Perform sub-iterations
        do while (time%it.le.time%itmax)

          wt_vel%time_in=parallel_time()

          ! Build mid-time velocity
          fs%U=0.5_WP*(fs%U+fs%Uold)
          fs%V=0.5_WP*(fs%V+fs%Vold)
          fs%W=0.5_WP*(fs%W+fs%Wold)

          ! Explicit calculation of drho*u/dt from NS
          call fs%get_dmomdt(resU,resV,resW)

          ! Assemble explicit residual
          resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
          resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
          resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
          wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

          ! Add linear forcing term: Bassenne et al. (2016)
          wt_force%time_in=parallel_time()
          ! we have stats from the previous step; they haven't changed
          !call compute_stats()
          linear_forcing: block
            real(WP) :: A
            ! - Eq. (7) (forcing constant TKE)
            A = (EPSp - (G / FORCE_TIMESCALE) * (TKE - TKE_target)) / (2.0_WP * TKE) * fs%rho
            ! update residuals
            resU = resU + time%dt * (fs%U - meanU) * A
            resV = resV + time%dt * (fs%V - meanV) * A
            resW = resW + time%dt * (fs%W - meanW) * A
          end block linear_forcing
          wt_force%time=wt_force%time+parallel_time()-wt_force%time_in

          ! Apply these residuals, apply bcs
          wt_vel%time_in=parallel_time()
          fs%U = 2.0_WP * fs%U - fs%Uold + resU / fs%rho
          fs%V = 2.0_WP * fs%V - fs%Vold + resV / fs%rho
          fs%W = 2.0_WP * fs%W - fs%Wold + resW / fs%rho
          call fs%apply_bcond(time%t,time%dt)
          wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

          ! Solve Poisson equation, correct velocity
          wt_pres%time_in=parallel_time()
          call fs%correct_mfr()
          call fs%get_div()
          fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
          call fs%psolv%solve()
          call fs%shift_p(fs%psolv%sol)
          call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
          fs%P = fs%P + fs%psolv%sol
          fs%U = fs%U - time%dt * resU / fs%rho
          fs%V = fs%V - time%dt * resV / fs%rho
          fs%W = fs%W - time%dt * resW / fs%rho
          wt_pres%time=wt_pres%time+parallel_time()-wt_pres%time_in

          ! Increment sub-iteration counter
          time%it=time%it+1

        end do

        ! Recompute divergence
        wt_vel%time_in=parallel_time()
        call fs%get_div()
        wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

        ! Compute turbulence stats
        wt_stat%time_in=parallel_time()
        call compute_stats()
        wt_stat%time=wt_stat%time+parallel_time()-wt_stat%time_in

        ! Output to ensight
        if (.not. ens_at_ints .and. ens_evt%occurs()) then
          call update_pmesh()
          call fs%interp_vel(Ui, Vi, Wi)
          call ens_out%write_data(time%t)
        end if

        ! Output monitoring
        call fs%get_max()
        call lp%get_max()
        call mfile%write()
        call cflfile%write()
        call hitfile%write()
        call lptfile%write()
        call ssfile%write()

        ! Monitor timing
        wt_total%time=parallel_time()-wt_total%time_in
        wt_vel%percent=wt_vel%time/wt_total%time*100.0_WP
        wt_pres%percent=wt_pres%time/wt_total%time*100.0_WP
        wt_lpt%percent=wt_lpt%time/wt_total%time*100.0_WP
        wt_stat%percent=wt_stat%time/wt_total%time*100.0_WP
        wt_force%percent=wt_force%time/wt_total%time*100.0_WP
        wt_rest%time=wt_total%time-wt_vel%time-wt_pres%time-wt_stat%time      &
          -wt_force%time
        wt_rest%percent=wt_rest%time/wt_total%time*100.0_WP
        call tfile%write()
        wt_total%time=0.0_WP; wt_total%percent=0.0_WP
        wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
        wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
        wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
        wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP
        wt_force%time=0.0_WP; wt_force%percent=0.0_WP
        wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP

      end do

      ! move to next parameters
      ec_next: block
        real(WP) :: interval

        call ec%compute_statistics(Re_lambda, Stk, phiinf, Wovk, urms, ETA, nu, time%t, time%n)
        call ec%ensight_write(time%t)

        if (ens_at_ints) then
          call update_pmesh()
          call ens_out%write_data(time%t)
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

    ! Get rid of all objects - need destructors
    ! monitor
    ! ensight
    ! bcond
    ! timetracker

    ! Deallocate work arrays
    deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradu,rho)

  end subroutine simulation_final

end module simulation
