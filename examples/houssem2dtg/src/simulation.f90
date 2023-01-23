!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use muscl_class,          only: muscl
  use hyperbolic_houssem2d, only: make_houssem2d_muscl, houssem2d_rhs
  use timetracker_class,    only: timetracker
  use ensight_class,        only: ensight
  use partmesh_class,       only: partmesh
  use event_class,          only: event
  use monitor_class,        only: monitor
  use string,               only: str_medium
  implicit none
  private

  !> Flow solver and a time tracker
  type(muscl), public :: fs
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> physical value arrays for ensight output
  !> up, vp
  real(WP), dimension(:,:,:,:), pointer :: p_vel

  !> tg parameters
  real(WP) :: tg_nu
  integer :: tg_a, tg_b

  !> src params and storage arrays
  real(WP) :: taup
  real(WP), dimension(2) :: gvec
  real(WP), dimension(:,:,:,:), allocatable :: src_rhs, src_mid

  !> simulation monitor file
  type(monitor) :: mfile, cflfile, consfile

  !> simulation control functions
  public :: simulation_init, simulation_run, simulation_final

contains

  !> update TG vortex
  subroutine update_tg(t)
    implicit none
    real(WP), intent(in) :: t
    integer :: i, j
    real(WP) :: sy, cy

    do j = cfg%jmin_, cfg%jmax_
      sy = cos(tg_b * cfg%ym(j)); cy = sin(tg_b * cfg%ym(j));
      do i = cfg%imin_, cfg%imax_
        fs%params(1,i,j,:) = sin(tg_a * cfg%xm(i)) * sy / (+tg_a)
        fs%params(2,i,j,:) = cos(tg_a * cfg%xm(i)) * cy / (-tg_b)
      end do
    end do

    if (tg_nu .ne. 0.0_WP) then
      fs%params(:,:,:,:) = exp(-2 * tg_nu * t) * fs%params(:,:,:,:)
    end if

  end subroutine update_tg

  subroutine rk2src()
    implicit none
    integer i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call  houssem2d_rhs(gvec, taup, fs%params(:,i,j,k), fs%Uc(:,i,j,k), &
            & src_rhs(:,i,j,k))
        end do
      end do
    end do

    src_rhs = src_rhs * (0.5_WP * time%dt)
    src_mid = fs%Uc + src_rhs

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call  houssem2d_rhs(gvec, taup, fs%params(:,i,j,k),                 &
            & src_mid(:,i,j,k), src_rhs(:,i,j,k))
        end do
      end do
    end do

    fs%dU = time%dt * src_rhs

  end subroutine

  !> update phys
  subroutine update_p_vel()
    implicit none

    p_vel(1,:,:,:) = fs%Uc(2,:,:,:) / fs%Uc(1,:,:,:)
    p_vel(2,:,:,:) = fs%Uc(3,:,:,:) / fs%Uc(1,:,:,:)

    call cfg%sync(p_vel)

  end subroutine

  !> initialization of problem solver
  subroutine simulation_init()
    use param, only: param_read
    use messager, only: die
    implicit none

    !TODO need to figure out what the right interaction with this is
    initialize_timetracker: block

      time = timetracker(amRoot=cfg%amRoot)
      call param_read('Max cfl number', time%cflmax)
      call param_read('Max timestep', time%dtmax, 'Max dt', huge(1.0_WP))
      time%dt = 1e-3_WP
      time%itmax = 1

    end block initialize_timetracker

    ! create a single-phase flow solver
    create_and_initialize_flow_solver: block

      ! call constructor
      fs = make_houssem2d_muscl(cfg)

      ! get TG params
      call param_read('Num vortices x', tg_a)
      call param_read('Num vortices x', tg_b)
      call param_read('Fluid kinematic viscosity', tg_nu, 'Fl kin vel', 0.0_WP)

      ! allocate source arrays
      allocate(src_rhs(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,         &
        & cfg%kmino_:cfg%kmaxo_), src_mid(6,cfg%imino_:cfg%imaxo_,            &
        & cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      ! get source params
      call param_read('Gravity', gvec(1))
      gvec(2) = 0.0_WP
      call param_read('Particle relax time', taup)

    end block create_and_initialize_flow_solver

    ! prepare initial fields
    initialize_fields: block
      real(WP) :: rhopn_init, v1_init, v2_init, e11_init, e12_init, e22_init

      call update_tg(time%t)

      call param_read('Initial rhopn', rhopn_init)
      call param_read('Initial x velocity', v1_init)
      call param_read('Initial y velocity', v2_init)
      call param_read('Initial e11', e11_init)
      call param_read('Initial e12', e12_init)
      call param_read('Initial e22', e22_init)

      fs%Uc(1,:,:,:) = rhopn_init
      fs%Uc(2,:,:,:) = rhopn_init * v1_init
      fs%Uc(3,:,:,:) = rhopn_init * v2_init
      fs%Uc(4,:,:,:) = e11_init
      fs%Uc(5,:,:,:) = e12_init
      fs%Uc(6,:,:,:) = e22_init

      call fs%recalc_cfl()

    end block initialize_fields

    ! Add Ensight output
    create_ensight: block
      use ensight_class, only: add_rscalar
      use string, only: str_short
      real(WP), dimension(:,:,:), pointer :: scl_ptr

      ! create array to hold particle velocities
      allocate(p_vel(2,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,        &
        & cfg%kmino_:cfg%kmaxo_))

      ! Create Ensight output from cfg
      ens_out = ensight(cfg=cfg, name='houssem2dtg')

      ! Create event for Ensight output
      ens_evt = event(time=time, name='Ensight output')
      call param_read('Ensight output period', ens_evt%tper)

      ! Add variables to output
      scl_ptr => fs%Uc(1,:,:,:)
      call add_rscalar(ens_out, 'pt_density', scl_ptr)
      scl_ptr => p_vel(1,:,:,:)
      call add_rscalar(ens_out, 'pt_x_velocity', scl_ptr)
      scl_ptr => p_vel(2,:,:,:)
      call add_rscalar(ens_out, 'pt_y_velocity', scl_ptr)
      scl_ptr => fs%params(1,:,:,:)
      call add_rscalar(ens_out, 'fluid_x_velocity', scl_ptr)
      scl_ptr => fs%params(2,:,:,:)
      call add_rscalar(ens_out, 'fluid_y_velocity', scl_ptr)
      scl_ptr => fs%Uc(4,:,:,:)
      call add_rscalar(ens_out, 'e11', scl_ptr)
      scl_ptr => fs%Uc(5,:,:,:)
      call add_rscalar(ens_out, 'e12', scl_ptr)
      scl_ptr => fs%Uc(6,:,:,:)
      call add_rscalar(ens_out, 'e22', scl_ptr)

      ! Output to ensight
      if (ens_evt%occurs()) then
        call ens_out%write_data(time%t)
      end if

    end block create_ensight

    ! Create a monitor file
    create_monitor: block
      use monitor_class, only: add_column_real
      use string, only: str_short
      real(WP), pointer :: real_ptr

      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_range()

      ! Create simulation monitor
      mfile = monitor(fs%cfg%amRoot, 'simulation')
      call mfile%add_column(time%n, 'Timestep number')
      call mfile%add_column(time%t, 'Time')
      call mfile%add_column(time%dt, 'Timestep size')
      call mfile%add_column(time%cfl, 'Maximum CFL')
      real_ptr => fs%Umin(1)
      call add_column_real(mfile, real_ptr, 'pt_density_min')
      real_ptr => fs%Umax(1)
      call add_column_real(mfile, real_ptr, 'pt_density_max')
      real_ptr => fs%Umin(2)
      call add_column_real(mfile, real_ptr, 'pt_momx_min')
      real_ptr => fs%Umax(2)
      call add_column_real(mfile, real_ptr, 'pt_momx_max')
      real_ptr => fs%Umin(3)
      call add_column_real(mfile, real_ptr, 'pt_momx_min')
      real_ptr => fs%Umax(3)
      call add_column_real(mfile, real_ptr, 'pt_momx_max')
      real_ptr => fs%Umin(4)
      call add_column_real(mfile, real_ptr, 'e11_min')
      real_ptr => fs%Umax(4)
      call add_column_real(mfile, real_ptr, 'e11_max')
      real_ptr => fs%Umin(5)
      call add_column_real(mfile, real_ptr, 'e12_min')
      real_ptr => fs%Umax(5)
      call add_column_real(mfile, real_ptr, 'e12_max')
      real_ptr => fs%Umin(6)
      call add_column_real(mfile, real_ptr, 'e22_min')
      real_ptr => fs%Umax(6)
      call add_column_real(mfile, real_ptr, 'e22_max')
      call mfile%write()

      ! Create CFL monitor
      cflfile = monitor(fs%cfg%amRoot, 'cfl')
      call cflfile%add_column(time%n, 'Timestep number')
      call cflfile%add_column(time%t, 'Time')
      call cflfile%add_column(fs%CFL_x, 'CFLx')
      call cflfile%add_column(fs%CFL_y, 'CFLy')
      call cflfile%add_column(fs%CFL_z, 'CFLz')
      call cflfile%write()

      ! Create conservation monitor
      consfile = monitor(fs%cfg%amRoot, 'conservation')
      call consfile%add_column(time%n, 'Timestep number')
      call consfile%add_column(time%t, 'Time')
      real_ptr => fs%Uint(1)
      call add_column_real(consfile, real_ptr, 'pt_density_int')
      real_ptr => fs%Uint(2)
      call add_column_real(consfile, real_ptr, 'pt_momx_int')
      real_ptr => fs%Uint(3)
      call add_column_real(consfile, real_ptr, 'pt_momy_int')
      real_ptr => fs%Uint(4)
      call add_column_real(consfile, real_ptr, 'e11_int')
      real_ptr => fs%Uint(5)
      call add_column_real(consfile, real_ptr, 'e12_int')
      real_ptr => fs%Uint(6)
      call add_column_real(consfile, real_ptr, 'e22_int')
      call consfile%write()

    end block create_monitor

  end subroutine simulation_init

  !> do the thing
  subroutine simulation_run
    implicit none

    ! Perform time integration
    do while (.not.time%done())

      ! Increment time
      call fs%get_cfl(time%dt, time%cfl)
      call time%adjust_dt()
      call time%increment()

      ! take step (Strang)
      call update_tg(time%t)
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_x(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_y(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call rk2src()
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_y(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_x(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU

      ! Output to ensight
      if (ens_evt%occurs()) then
        call update_p_vel()
        call ens_out%write_data(time%t)
      end if

      ! Perform and output monitoring
      call fs%get_range()
      call mfile%write()
      call cflfile%write()
      call consfile%write()

    end do

  end subroutine simulation_run

  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    ! Get rid of all objects - need destructors
    !TODO
    ! monitor
    ! ensight
    ! bcond
    ! timetracker


  end subroutine simulation_final

end module simulation

