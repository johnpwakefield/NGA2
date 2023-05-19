!> Definition for a crossflow class
module crossflow_class
   use precision,         only: WP
   use config_class,      only: config
   use iterator_class,    only: iterator
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: crossflow
   
   !> Crossflow object
   type :: crossflow
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(incomp)      :: fs    !< NS flow solver
      type(lpt)         :: lp    !< LPT solver
      type(fft2d)       :: ps    !< 2D FFT solver for pressure
      type(ddadi)       :: vs    !< DDADI solver for velocity
      type(timetracker) :: time  !< Time info
      
      !> Ensight postprocessing
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      type(partmesh) :: pmesh    !< Particle mesh for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      type(monitor) :: lptfile  !< LPT monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities

      !> Flow properties
      real(WP) :: visc,Ubulk
      
   contains
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type crossflow
   

contains
   
   
   !> Initialization of crossflow simulation
   subroutine init(this)
      implicit none
      class(crossflow), intent(inout) :: this
      
      
      ! Create the crossflow mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         ! Read in grid definition
         call param_read('2 Lx',Lx); call param_read('2 nx',nx); allocate(x(nx+1))
         call param_read('2 Ly',Ly); call param_read('2 ny',ny); allocate(y(ny+1))
         call param_read('2 Lz',Lz); call param_read('2 nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='Crossflow')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
         this%cfg%VF=1.0_WP
      end block create_config
      

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         call param_read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use param,           only: param_read
         use incomp_class,    only: dirichlet,clipped_neumann,bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k      
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='NS solver')
         call param_read('Density',this%fs%rho)
         call param_read('Dynamic viscosity',this%visc)
         this%fs%visc=this%visc
         call param_read('Bulk velocity',this%Ubulk)
         ! Define inflow boundary condition on the left
         call this%fs%add_bcond(name='inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Define outflow boundary condition on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         ! Configure pressure solver
         this%ps=fft2d(cfg=this%cfg,name='Pressure 2',nst=7)
         ! Configure implicit velocity solver
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply convective velocity
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=this%Ubulk
         end do
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
       end block create_flow_solver


       ! Initialize LPT solver
       initialize_lpt: block
         use param, only: param_read
         ! Create solver
         this%lp=lpt(cfg=this%cfg,name='LPT')
         ! Set various particle parameters
         call param_read('Gravity',this%lp%gravity)
         call param_read('Particle density',this%lp%rho)
         this%lp%drag_model='Schiller-Naumann'
         ! Initialize with zero particles
         call this%lp%resize(0)
         ! Injection parameters
         call param_read('Particle mass flow rate',this%lp%mfr)
         call param_read('Particle mean diameter',this%lp%inj_dmean)
         this%lp%inj_dsd=0.0_WP; this%lp%inj_dshift=0.0_WP; this%lp%inj_dmin=this%lp%inj_dmean; this%lp%inj_dmax=this%lp%inj_dmean
         this%lp%inj_pos(1)=this%lp%cfg%x(this%lp%cfg%imin)+0.5_WP*this%lp%inj_dmax
         this%lp%inj_vel=(/ this%Ubulk,0.0_WP, 0.0_WP /)
         this%lp%inj_T=300.0_WP
         ! Get initial particle volume fraction
         call this%lp%update_VF()
       end block initialize_lpt


       ! Create partmesh object for Lagrangian particle output
       create_pmesh: block
         integer :: i
         this%pmesh=partmesh(nvar=0,nvec=1,name='lpt')
         this%pmesh%vecname(1)='velocity'
         call this%lp%update_partmesh(this%pmesh)
         do i=1,this%lp%np_
            this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
         end do
       end block create_pmesh


      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='crossflow')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_particle('particles',this%pmesh)
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_xflow')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
         ! Create LPT monitor
         call this%lp%get_max()
         this%lptfile=monitor(this%lp%cfg%amRoot,'lpt')
         call this%lptfile%add_column(this%time%n,'Timestep number')
         call this%lptfile%add_column(this%time%t,'Time')
         call this%lptfile%add_column(this%lp%np,'Particle number')
         call this%lptfile%add_column(this%lp%Umin,'Particle Umin')
         call this%lptfile%add_column(this%lp%Umax,'Particle Umax')
         call this%lptfile%add_column(this%lp%Vmin,'Particle Vmin')
         call this%lptfile%add_column(this%lp%Vmax,'Particle Vmax')
         call this%lptfile%add_column(this%lp%Wmin,'Particle Wmin')
         call this%lptfile%add_column(this%lp%Wmax,'Particle Wmax')
         call this%lptfile%write()
      end block create_monitor
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(crossflow), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()

      ! Inject and advance particles
      call this%lp%inject(dt=this%time%dt)
      this%resU=this%fs%rho; this%resV=this%fs%visc
      call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*(this%fs%U-this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%V-this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%W-this%fs%Wold)+this%time%dt*this%resW 
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Solve Poisson equation
         call this%fs%correct_mfr()
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU
         this%fs%V=this%fs%V-this%time%dt*this%resV
         this%fs%W=this%fs%W-this%time%dt*this%resW
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         update_pmesh: block
           integer :: i
           call this%lp%update_partmesh(this%pmesh)
           do i=1,this%lp%np_
              this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
           end do
         end block update_pmesh
         ! Perform ensight output
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%lp%get_max()
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%lptfile%write()
      
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(crossflow), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
   end subroutine final
   
   
   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   
   
end module crossflow_class
