!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use muscl_class,          only: muscl
  use hyperbolic_euler,     only: make_euler_muscl, euler_tocons,             &
    & euler_tophys, SOD_PHYS_L, SOD_PHYS_R, DIATOMIC_GAMMA
  use timetracker_class,    only: timetracker
  use ensight_class,        only: ensight
  use partmesh_class,       only: partmesh
  use event_class,          only: event
  use monitor_class,        only: monitor
  implicit none
  private

  !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
  type(muscl), public :: fs
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> physical value arrays for ensight output
  real(WP), dimension(:,:,:,:), pointer :: phys_out

  !> simulation monitor file
  type(monitor) :: mfile, cflfile

  !> simulation control functions
  public :: simulation_init, simulation_run, simulation_final

contains

  !> update phys
  subroutine update_phys_out()
    implicit none
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call euler_tophys(DIATOMIC_GAMMA, fs%Uc(:,i,j,k), phys_out(:,i,j,k))
        end do
      end do
    end do

    call cfg%sync(phys_out)

  end subroutine

  !> initialization of problem solver
  subroutine simulation_init()
    use param, only: param_read
    use string, only: str_medium
    use messager, only: die
    implicit none

    !TODO need to figure out what the right interaction with this is
    initialize_timetracker: block

      time = timetracker(amRoot=cfg%amRoot)
      call param_read('Max cfl number', time%cflmax)
      time%dt = 1e-3_WP
      time%itmax = 1

    end block initialize_timetracker

    ! create a single-phase flow solver without bconds
    create_and_initialize_flow_solver: block

      ! call constructor
      fs = make_euler_muscl(cfg)

    end block create_and_initialize_flow_solver

    ! prepare initial fields
    initialize_fields: block
      real(WP) :: r2, oR
      real(WP), dimension(3) :: x, c
      real(WP), dimension(5) :: insval, outval
      integer :: i, j, k

      call param_read('Initial diameter', oR)
      oR = 0.5_WP * oR

      call euler_tocons(DIATOMIC_GAMMA, SOD_PHYS_L, insval)
      call euler_tocons(DIATOMIC_GAMMA, SOD_PHYS_R, outval)

      c = 0.5 * (/ cfg%xL, cfg%yL, cfg%zL /)

      do k = cfg%kmin_, cfg%kmax_
        x(3) = cfg%zm(k)
        do j = cfg%jmin_, cfg%jmax_
          x(2) = cfg%ym(j)
          do i = cfg%imin_, cfg%imax_
            x(1) = cfg%xm(i)
            r2 = sum((c - x)**2)
            if (r2 .lt. oR**2) then
              fs%Uc(:,i,j,k) = insval
            else
              fs%Uc(:,i,j,k) = outval
            end if
          end do
        end do
      end do

      call fs%recalc_cfl()

    end block initialize_fields

    ! Add Ensight output
    create_ensight: block
      use ensight_class, only: add_rscalar
      use string, only: str_short
      real(WP), dimension(:,:,:), pointer :: scl_ptr

      ! create array to hold physical coordinates
      allocate(phys_out(5,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,        &
        & cfg%kmino_:cfg%kmaxo_))

      ! Create Ensight output from cfg
      ens_out = ensight(cfg=cfg, name='SodSphereTest')

      ! Create event for Ensight output
      ens_evt = event(time=time, name='Ensight output')
      call param_read('Ensight output period', ens_evt%tper)

      ! Add variables to output
      scl_ptr => phys_out(1,:,:,:)
      call add_rscalar(ens_out, 'density', scl_ptr)
      scl_ptr => phys_out(2,:,:,:)
      call add_rscalar(ens_out, 'x_velocity', scl_ptr)
      scl_ptr => phys_out(3,:,:,:)
      call add_rscalar(ens_out, 'y_velocity', scl_ptr)
      scl_ptr => phys_out(4,:,:,:)
      call add_rscalar(ens_out, 'z_velocity', scl_ptr)
      scl_ptr => phys_out(5,:,:,:)
      call add_rscalar(ens_out, 'pressure', scl_ptr)
      scl_ptr => fs%params(1,:,:,:)
      call add_rscalar(ens_out, 'gamma', scl_ptr)

      ! Output to ensight
      if (ens_evt%occurs()) then
        call update_phys_out()
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
      ! not sure why this doesn't work?
      !call mfile%add_column(real_ptr, fields(i:i)//'min')
      !call mfile%add_column(real_ptr, fields(i:i)//'max')
      real_ptr => fs%Umin(1)
      call add_column_real(mfile, real_ptr, 'density_min')
      real_ptr => fs%Umax(1)
      call add_column_real(mfile, real_ptr, 'density_max')
      real_ptr => fs%Umin(2)
      call add_column_real(mfile, real_ptr, 'momentum_x_min')
      real_ptr => fs%Umax(2)
      call add_column_real(mfile, real_ptr, 'momentum_x_max')
      real_ptr => fs%Umin(3)
      call add_column_real(mfile, real_ptr, 'momentum_y_min')
      real_ptr => fs%Umax(3)
      call add_column_real(mfile, real_ptr, 'momentum_y_max')
      real_ptr => fs%Umin(4)
      call add_column_real(mfile, real_ptr, 'momentum_z_min')
      real_ptr => fs%Umax(4)
      call add_column_real(mfile, real_ptr, 'momentum_z_max')
      real_ptr => fs%Umin(5)
      call add_column_real(mfile, real_ptr, 'totalenergy_min')
      real_ptr => fs%Umax(5)
      call add_column_real(mfile, real_ptr, 'totalenergy_max')
      call mfile%write()

      ! Create CFL monitor
      cflfile = monitor(fs%cfg%amRoot, 'cfl')
      call cflfile%add_column(time%n, 'Timestep number')
      call cflfile%add_column(time%t, 'Time')
      call cflfile%add_column(fs%CFL_x, 'CFLx')
      call cflfile%add_column(fs%CFL_y, 'CFLy')
      call cflfile%add_column(fs%CFL_z, 'CFLz')
      call cflfile%write()

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
      fs%dU(:, :, :, :) = 0.0_WP
      call fs%compute_dU_x(0.5 * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:, :, :, :) = 0.0_WP
      call fs%compute_dU_y(time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:, :, :, :) = 0.0_WP
      call fs%compute_dU_x(0.5 * time%dt)
      fs%Uc = fs%Uc + fs%dU

      ! apply boundary conditions
      !call fs%apply_bcond(time%t, time%dt)

      ! Output to ensight
      if (ens_evt%occurs()) then
        call update_phys_out()
        call ens_out%write_data(time%t)
      end if

      ! Perform and output monitoring
      call fs%get_range()
      call mfile%write()
      call cflfile%write()

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

