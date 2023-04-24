!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use viehlli2d_class,      only: viehlli2d, viehlli2d_cons2prim, viehlli2d_prim2cons, viehlli2d_isvalid
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
  type(viehlli2d), target, public :: fs
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> physical (primitive) value arrays for ensight output
  real(WP), dimension(:,:,:,:), pointer :: physout
  real(WP), dimension(:,:,:), pointer :: scl_zeros

  !> fluid velocity array
  real(WP), dimension(:,:,:,:), pointer :: flvel

  !> tg parameters
  real(WP) :: tg_nu
  integer :: tg_a, tg_b
  real(WP), dimension(2) :: tg_offset

  !> src params and storage arrays
  real(WP) :: taup
  real(WP), dimension(2) :: gvec
  real(WP), dimension(:,:,:,:), allocatable :: rhs, src_mid

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

    if (tg_a .eq. 0 .or. tg_b .eq. 0) then
      flvel(:,:,:,:) = 0.0_WP
    else
      do j = cfg%jmino_, cfg%jmaxo_
        sy = sin(2 * pi * tg_b * cfg%ym(j)); cy = cos(2 * pi * tg_b * cfg%ym(j));
        do i = cfg%imino_, cfg%imaxo_
          flvel(1,i,j,:) = sin(2 * pi * tg_a * cfg%xm(i)) * cy / (+tg_a)
          flvel(2,i,j,:) = cos(2 * pi * tg_a * cfg%xm(i)) * sy / (-tg_b)
        end do
      end do
    end if

    if (tg_nu .ne. 0.0_WP) then
      flvel(:,:,:,:) = exp(-2 * tg_nu * t) * flvel(:,:,:,:)
    end if

    flvel(1,:,:,:) = tg_offset(1) + flvel(1,:,:,:)
    flvel(2,:,:,:) = tg_offset(2) + flvel(2,:,:,:)

  end subroutine update_tg

  !subroutine rk2src()
  !  implicit none
  !  integer :: i, j, k
!
  !  do k = cfg%kmin_, cfg%kmax_
  !    do j = cfg%jmin_, cfg%jmax_
  !      do i = cfg%imin_, cfg%imax_
  !        call  viehlli2d_rhs(gvec, taup, flvel(:,i,j,k), fs%Uc(:,i,j,k), rhs(:,i,j,k))
  !      end do
  !    end do
  !  end do
!
  !  rhs = rhs * (0.5_WP * time%dt)
  !  src_mid = fs%Uc + rhs
!
  !  do k = cfg%kmin_, cfg%kmax_
  !    do j = cfg%jmin_, cfg%jmax_
  !      do i = cfg%imin_, cfg%imax_
  !        call fs%rhs(gvec, taup, flvel(:,i,j,k), src_mid(:,i,j,k), rhs(:,i,j,k))
  !      end do
  !    end do
  !  end do
!
  !  fs%dU = time%dt * rhs
!
  !end subroutine rk2src

  subroutine eulersrc(dt)
    implicit none
    real(WP), intent(in) :: dt
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call fs%rhs(gvec, taup, flvel(:,i,j,k), fs%Uc(:,i,j,k), rhs(:,i,j,k))
        end do
      end do
    end do

    fs%dU(:,:,:,:) = fs%dU(:,:,:,:) + dt * rhs(:,:,:,:)

  end subroutine eulersrc

  subroutine backeulersrc(dt)
    implicit none
    real(WP), intent(in) :: dt
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call fs%backeuler(gvec, taup, flvel(:,i,j,k), dt, fs%Uc(:,i,j,k), rhs(:,i,j,k))
        end do
      end do
    end do

    fs%dU(:,:,:,:) = fs%dU(:,:,:,:) + dt * rhs(:,:,:,:)

  end subroutine backeulersrc

  !> update phys
  subroutine update_physout()
    implicit none
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call viehlli2d_cons2prim(fs%Uc(:,i,j,k), physout(:,i,j,k))
        end do
      end do
    end do

    call cfg%sync(physout)

  end subroutine update_physout

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
      fs = viehlli2d(cfg, .false.)

      ! get TG params
      call param_read('Num vortices x', tg_a)
      call param_read('Num vortices y', tg_b)
      call param_read('Fluid kinematic viscosity', tg_nu, 'Fl kin vel', 0.0_WP)
      call param_read('Fluid offset velocity x', tg_offset(1), 'Fl off x', 0.0_WP)
      call param_read('Fluid offset velocity y', tg_offset(2), 'Fl off y', 0.0_WP)

      ! allocate source arrays
      allocate(rhs(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,             &
        & cfg%kmino_:cfg%kmaxo_), src_mid(6,cfg%imino_:cfg%imaxo_,            &
        & cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      ! allocate fluid velocity array
      allocate(flvel(2,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,           &
        & cfg%kmino_:cfg%kmaxo_))

      ! get source params
      call param_read('Gravity', gvec(1))
      gvec(2) = 0.0_WP
      call param_read('Particle relax time', taup)

    end block create_and_initialize_flow_solver

    ! prepare initial fields
    initialize_fields: block
      real(WP) :: rhopn_in, rhopn_out, v_init_scale_in, v_init_scale_out,     &
        Pp11_in, Pp12_in, Pp22_in, Pp11_out, Pp12_out, Pp22_out, r, x, y
      real(WP), dimension(6) :: U
      integer :: i, j, k

      call update_tg(time%t)

      call param_read('Initial rhopn inside', rhopn_in)
      call param_read('Initial rhopn outside', rhopn_out)
      call param_read('Initial velocity scale inside', v_init_scale_in)
      call param_read('Initial velocity scale outside', v_init_scale_out)
      call param_read('Initial Pp11 inside', Pp11_in)
      call param_read('Initial Pp12 inside', Pp12_in)
      call param_read('Initial Pp22 inside', Pp22_in)
      call param_read('Initial Pp11 outside', Pp11_out)
      call param_read('Initial Pp12 outside', Pp12_out)
      call param_read('Initial Pp22 outside', Pp22_out)
      call param_read('Initial circle radius', r)

      do k = cfg%kmin_, cfg%kmax_
        do j = cfg%jmin_, cfg%jmax_
          y = cfg%ym(j)
          do i = cfg%imin_, cfg%imax_
            x = cfg%xm(i)
            if ((x - 0.5_WP * cfg%xL)**2 + (y - 0.5_WP * cfg%yL)**2 < r**2) then
              U(1) = rhopn_in
              U(2) = v_init_scale_in * flvel(1,i,j,k)
              U(3) = v_init_scale_in * flvel(2,i,j,k)
              U(4) = Pp11_in; U(5) = Pp12_in; U(6) = Pp22_in;
            else
              U(1) = rhopn_out
              U(2) = v_init_scale_out * flvel(1,i,j,k)
              U(3) = v_init_scale_out * flvel(2,i,j,k)
              U(4) = Pp11_out; U(5) = Pp12_out; U(6) = Pp22_out;
            end if
            call viehlli2d_prim2cons(U, fs%Uc(:,i,j,k))
          end do
        end do
      end do

      call cfg%sync(fs%Uc)

    end block initialize_fields

    ! Add Ensight output
    create_ensight: block
      use string, only: str_short
      real(WP), dimension(:,:,:), pointer :: scl_ptr_1, scl_ptr_2

      ! create array to hold primitive variables
      allocate(physout(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,         &
        & cfg%kmino_:cfg%kmaxo_))

      ! Create Ensight output from cfg
      ens_out = ensight(cfg=cfg, name='viehlli2dtg')

      ! Create event for Ensight output
      ens_evt = event(time=time, name='Ensight output')
      call param_read('Ensight output period', ens_evt%tper)

      ! allocate and set scl_zeros to zero
      allocate(scl_zeros(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,         &
        cfg%kmino_:cfg%kmaxo_))
      scl_zeros(:,:,:) = 0.0_WP

      ! Add variables to output
      scl_ptr_1 => physout(1,:,:,:)
      call ens_out%add_scalar('pt_density', scl_ptr_1)
      scl_ptr_1 => physout(2,:,:,:)
      scl_ptr_2 => physout(3,:,:,:)
      call ens_out%add_vector('pt_vel', scl_ptr_1, scl_ptr_2, scl_zeros)
      scl_ptr_1 => flvel(1,:,:,:)
      scl_ptr_2 => flvel(2,:,:,:)
      call ens_out%add_vector('fl_vel', scl_ptr_1, scl_ptr_2, scl_zeros)
      scl_ptr_1 => physout(4,:,:,:)
      call ens_out%add_scalar('Pp11', scl_ptr_1)
      scl_ptr_1 => physout(5,:,:,:)
      call ens_out%add_scalar('Pp12', scl_ptr_1)
      scl_ptr_1 => physout(6,:,:,:)
      call ens_out%add_scalar('Pp22', scl_ptr_1)

      ! Output to ensight
      if (ens_evt%occurs()) then
        call update_physout()
        call ens_out%write_data(time%t)
      end if

    end block create_ensight

    ! Create a monitor file
    create_monitor: block
      use string, only: str_short
      real(WP), pointer :: real_ptr

      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_max()

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
      call rangefile%add_column(real_ptr, 'E11_min')
      real_ptr => fs%Umax(4)
      call rangefile%add_column(real_ptr, 'E11_max')
      real_ptr => fs%Umin(5)
      call rangefile%add_column(real_ptr, 'E12_min')
      real_ptr => fs%Umax(5)
      call rangefile%add_column(real_ptr, 'E12_max')
      real_ptr => fs%Umin(6)
      call rangefile%add_column(real_ptr, 'E22_min')
      real_ptr => fs%Umax(6)
      call rangefile%add_column(real_ptr, 'E22_max')
      call rangefile%write()

      ! Create CFL monitor
      cflfile = monitor(fs%cfg%amRoot, 'cfl')
      call cflfile%add_column(time%n, 'Timestep number')
      call cflfile%add_column(time%t, 'Time')
      call cflfile%add_column(fs%CFL_x, 'CFLx')
      call cflfile%add_column(fs%CFL_y, 'CFLy')
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
      call consfile%add_column(real_ptr, 'E11_int')
      real_ptr => fs%Uint(5)
      call consfile%add_column(real_ptr, 'E12_int')
      real_ptr => fs%Uint(6)
      call consfile%add_column(real_ptr, 'E22_int')
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
      time%cfl = max(time%cfl, time%dt / (taup * time%cflmax))
      call time%adjust_dt()
      call time%increment()

      ! update vortex
      call update_tg(time%t)

      ! take step (Strang)
      ! there were some maybe not ideal choices made about how ghost cells
      ! are handled
      ! after a sync then x step, cells in the boundary region needed for the
      ! y step are accurate, as compute_dU_x operates in ghost cells as well
      ! after the y step, the cells needed in the x direction are not accurate
      ! because the ghost corners were not updated between these stages
      call cfg%sync(fs%Uc)
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_x(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      call cfg%sync(fs%Uc)
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_y(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      call cfg%sync(fs%Uc)
      fs%dU(:,:,:,:) = 0.0_WP
      call eulersrc()
      fs%Uc = fs%Uc + fs%dU
      call cfg%sync(fs%Uc)
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_y(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU
      call cfg%sync(fs%Uc)
      fs%dU(:,:,:,:) = 0.0_WP
      call fs%compute_dU_x(0.5_WP * time%dt)
      fs%Uc = fs%Uc + fs%dU

      ! Output to ensight
      if (ens_evt%occurs()) then
        call cfg%sync(fs%Uc)
        call update_physout()
        call ens_out%write_data(time%t)
      end if

      ! Perform and output monitoring
      call fs%get_max()
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

