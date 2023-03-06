!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use muscl_class,          only: muscl
  use hyperbolic_houssem2d, only: make_houssem2d_muscl, houssem2d_rhs, houssem2d_backeuler
  use timetracker_class,    only: timetracker
  use ensight_class,        only: ensight
  use partmesh_class,       only: partmesh
  use event_class,          only: event
  use monitor_class,        only: monitor
  use string,               only: str_medium
  implicit none
  private

  real(WP), parameter :: pi = 4.0_WP * atan(1.0_WP)

  !> Flow solver and a time tracker
  type(muscl), public :: fs
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> physical value arrays for ensight output
  !> up, vp
  real(WP), dimension(:,:,:,:), pointer :: p_vel
  real(WP), dimension(:,:,:), pointer :: scl_zeros

  !> debug values for ensight output
  real(WP), dimension(:,:,:), pointer :: csquared

  !> tg parameters
  real(WP) :: tg_nu
  integer :: tg_a, tg_b

  !> src params and storage arrays
  real(WP), dimension(3) :: eqPp
  real(WP) :: taup
  real(WP), dimension(2) :: gvec
  real(WP), dimension(:,:,:,:), allocatable :: src_rhs, src_mid

  !> simulation monitor file
  type(monitor) :: mfile, cflfile, consfile, rangefile

  !> simulation control functions
  public :: simulation_init, simulation_run, simulation_final

contains

  !> update TG vortex
  subroutine update_tg(t)
    implicit none
    real(WP), intent(in) :: t
    integer :: i, j
    real(WP) :: sy, cy

    do j = cfg%jmino_, cfg%jmaxo_
      sy = sin(2 * pi * tg_b * cfg%ym(j)); cy = cos(2 * pi * tg_b * cfg%ym(j));
      do i = cfg%imino_, cfg%imaxo_
        fs%params(1,i,j,:) = sin(2 * pi * tg_a * cfg%xm(i)) * cy / (+tg_a)
        fs%params(2,i,j,:) = cos(2 * pi * tg_a * cfg%xm(i)) * sy / (-tg_b)
      end do
    end do

    if (tg_nu .ne. 0.0_WP) then
      fs%params(:,:,:,:) = exp(-2 * tg_nu * t) * fs%params(:,:,:,:)
    end if

  end subroutine update_tg

  subroutine rk2src()
    implicit none
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call  houssem2d_rhs(eqPp, gvec, taup, fs%params(:,i,j,k),           &
            fs%Uc(:,i,j,k), src_rhs(:,i,j,k))
        end do
      end do
    end do

    src_rhs = src_rhs * (0.5_WP * time%dt)
    src_mid = fs%Uc + src_rhs

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call  houssem2d_rhs(eqPp, gvec, taup, fs%params(:,i,j,k),           &
            src_mid(:,i,j,k), src_rhs(:,i,j,k))
        end do
      end do
    end do

    fs%dU = time%dt * src_rhs

  end subroutine rk2src

  subroutine eulersrc()
    implicit none
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call  houssem2d_rhs(eqPp, gvec, taup, fs%params(:,i,j,k),           &
            fs%Uc(:,i,j,k), src_rhs(:,i,j,k))
        end do
      end do
    end do

    fs%dU = time%dt * src_rhs

  end subroutine eulersrc

  subroutine backeulersrc(dt)
    implicit none
    real(WP), intent(in) :: dt
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call  houssem2d_backeuler(eqPp, gvec, taup, fs%params(:,i,j,k), dt, &
            fs%Uc(:,i,j,k), fs%dU(:,i,j,k))
        end do
      end do
    end do

  end subroutine backeulersrc

  subroutine halfhalfsrc(dt)
    implicit none
    real(WP), intent(in) :: dt
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call houssem2d_rhs(eqPp, gvec, taup, fs%params(:,i,j,k),            &
            fs%Uc(:,i,j,k), src_rhs(:,i,j,k))
        end do
      end do
    end do

    src_rhs = src_rhs * (0.5 * time%dt)
    src_mid = fs%Uc + src_rhs

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call  houssem2d_backeuler(eqPp, gvec, taup, fs%params(:,i,j,k),     &
            0.5_WP * dt, src_mid(:,i,j,k), fs%dU(:,i,j,k))
        end do
      end do
    end do

    fs%dU = fs%dU + src_rhs

  end subroutine halfhalfsrc(dt)

  !> update phys
  subroutine update_p_vel()
    implicit none

    p_vel(1,:,:,:) = fs%Uc(2,:,:,:) / fs%Uc(1,:,:,:)
    p_vel(2,:,:,:) = fs%Uc(3,:,:,:) / fs%Uc(1,:,:,:)

    call cfg%sync(p_vel)

  end subroutine update_p_vel

  !> update debug
  subroutine update_debug()
    implicit none

    csquared(:,:,:) = fs%Uc(4,:,:,:) - fs%Uc(2,:,:,:)**2 / fs%Uc(1,:,:,:)
    csquared(:,:,:) = csquared(:,:,:) / fs%Uc(1,:,:,:)

    call cfg%sync(csquared)

  end subroutine update_debug

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
      real(WP) :: eqPpii

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
      call param_read('Equilibrium particle pressure', eqPpii)
      eqPp(:) = (/ eqPpii, 0.0_WP, eqPpii /)

    end block create_and_initialize_flow_solver

    ! prepare initial fields
    initialize_fields: block
      real(WP) :: rhopn_init, v_init_scale, e11_init, e12_init, e22_init

      call update_tg(time%t)

      call param_read('Initial rhopn', rhopn_init)
      call param_read('Initial velocity scale', v_init_scale)
      call param_read('Initial e11', e11_init)
      call param_read('Initial e12', e12_init)
      call param_read('Initial e22', e22_init)

      fs%Uc(1,:,:,:) = rhopn_init
      fs%Uc(2,:,:,:) = rhopn_init * v_init_scale * fs%params(1,:,:,:)
      fs%Uc(3,:,:,:) = rhopn_init * v_init_scale * fs%params(2,:,:,:)
      fs%Uc(4,:,:,:) = rhopn_init * e11_init
      fs%Uc(5,:,:,:) = rhopn_init * e12_init
      fs%Uc(6,:,:,:) = rhopn_init * e22_init

      call fs%recalc_cfl()

    end block initialize_fields

    ! Add Ensight output
    create_ensight: block
      use string, only: str_short
      real(WP), dimension(:,:,:), pointer :: scl_ptr_1, scl_ptr_2

      ! create array to hold particle velocities
      allocate(p_vel(2,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,           &
        & cfg%kmino_:cfg%kmaxo_))
      allocate(scl_zeros(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,         &
        cfg%kmino_:cfg%kmaxo_))

      ! create debug value arrays
      allocate(csquared(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,          &
        & cfg%kmino_:cfg%kmaxo_))

      ! Create Ensight output from cfg
      ens_out = ensight(cfg=cfg, name='houssem2dtg')

      ! Create event for Ensight output
      ens_evt = event(time=time, name='Ensight output')
      call param_read('Ensight output period', ens_evt%tper)

      ! Add variables to output
      scl_ptr_1 => fs%Uc(1,:,:,:)
      call ens_out%add_scalar('pt_density', scl_ptr_1)
      scl_ptr_1 => p_vel(1,:,:,:)
      scl_ptr_2 => p_vel(2,:,:,:)
      call ens_out%add_vector('pt_vel', scl_ptr_1, scl_ptr_2, scl_zeros)
      scl_ptr_1 => fs%params(1,:,:,:)
      scl_ptr_2 => fs%params(2,:,:,:)
      call ens_out%add_vector('fl_vel', scl_ptr_1, scl_ptr_2, scl_zeros)
      scl_ptr_1 => fs%Uc(4,:,:,:)
      call ens_out%add_scalar('E11', scl_ptr_1)
      scl_ptr_1 => fs%Uc(5,:,:,:)
      call ens_out%add_scalar('E12', scl_ptr_1)
      scl_ptr_1 => fs%Uc(6,:,:,:)
      call ens_out%add_scalar('E22', scl_ptr_1)
      scl_ptr_1 => csquared(:,:,:)
      call ens_out%add_scalar('c2', scl_ptr_1)

      ! Output to ensight
      if (ens_evt%occurs()) then
        call ens_out%write_data(time%t)
      end if

    end block create_ensight

    ! Create a monitor file
    create_monitor: block
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
      call mfile%write()

      ! Create range monitor
      rangefile = monitor(fs%cfg%amRoot, 'range')
      call rangefile%add_column(time%n, 'Timestep number')
      call rangefile%add_column(time%t, 'Time')
      real_ptr => fs%Umin(1)
      call rangefile%add_column(real_ptr, 'pt_density_min')
      real_ptr => fs%Umax(1)
      call rangefile%add_column(real_ptr, 'pt_density_max')
      real_ptr => fs%Umin(2)
      call rangefile%add_column(real_ptr, 'pt_momx_min')
      real_ptr => fs%Umax(2)
      call rangefile%add_column(real_ptr, 'pt_momx_max')
      real_ptr => fs%Umin(3)
      call rangefile%add_column(real_ptr, 'pt_momx_min')
      real_ptr => fs%Umax(3)
      call rangefile%add_column(real_ptr, 'pt_momx_max')
      real_ptr => fs%Umin(4)
      call rangefile%add_column(real_ptr, 'e11_min')
      real_ptr => fs%Umax(4)
      call rangefile%add_column(real_ptr, 'e11_max')
      real_ptr => fs%Umin(5)
      call rangefile%add_column(real_ptr, 'e12_min')
      real_ptr => fs%Umax(5)
      call rangefile%add_column(real_ptr, 'e12_max')
      real_ptr => fs%Umin(6)
      call rangefile%add_column(real_ptr, 'e22_min')
      real_ptr => fs%Umax(6)
      call rangefile%add_column(real_ptr, 'e22_max')
      call rangefile%write()

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
      call consfile%add_column(real_ptr, 'pt_density_int')
      real_ptr => fs%Uint(2)
      call consfile%add_column(real_ptr, 'pt_momx_int')
      real_ptr => fs%Uint(3)
      call consfile%add_column(real_ptr, 'pt_momy_int')
      real_ptr => fs%Uint(4)
      call consfile%add_column(real_ptr, 'e11_int')
      real_ptr => fs%Uint(5)
      call consfile%add_column(real_ptr, 'e12_int')
      real_ptr => fs%Uint(6)
      call consfile%add_column(real_ptr, 'e22_int')
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
      time%cfl = max(time%cfl, 0.5_WP / (time%cflmax * taup) * time%dt)
      call time%adjust_dt()
      call time%increment()

      ! update vortex
      call update_tg(time%t)

      ! take step (Strang)
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_x(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_y(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call backeulersrc(time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_y(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_x(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU

      ! take step (Lie backward source, unsplit space)
      !fs%dU(:,:,:,:) = 0.0_WP
      !call backeulersrc(time%dt)
      !fs%Uc = fs%Uc + fs%dU
      !fs%dU(:,:,:,:) = 0.0_WP
      !call fs%compute_dU_x(time%dt)
      !call fs%compute_dU_y(time%dt)
      !fs%Uc = fs%Uc + fs%dU

      ! take step (unsplit)
      !fs%dU(:,:,:,:) = 0.0_WP
      !call backeulersrc(time%dt)
      !call fs%compute_dU_x(time%dt)
      !call fs%compute_dU_y(time%dt)
      !fs%Uc = fs%Uc + fs%dU

      ! Output to ensight
      if (ens_evt%occurs()) then
        call update_p_vel()
        call update_debug()
        call ens_out%write_data(time%t)
      end if

      ! Perform and output monitoring
      call fs%get_range()
      call mfile%write()
      call rangefile%write()
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

