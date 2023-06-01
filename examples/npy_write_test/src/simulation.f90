!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use ensight_class,        only: ensight
  use npy_class,            only: npy
  use partmesh_class,       only: partmesh
  use lpt_class,            only: lpt
  use string,               only: str_medium
  implicit none
  private

  !> particles
  type(lpt), target, public :: lp

  !> Partmesh
  type(partmesh)       :: pmesh

  !> Ensight postprocessing
  type(ensight)        :: ens_out

  !> Npy postprocessing
  type(npy)            :: npy_out

  public :: simulation_init, simulation_run, simulation_final

  !> Arrays
  real(WP), dimension(:,:,:), target, allocatable :: rho, sine, sqpsq, sqmsq
  real(WP), dimension(:,:,:,:), target, allocatable :: vfield

contains

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
      pmesh%vec(:,2,i) =  0.0_WP
    end do

  end subroutine update_pmesh

  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    use mathtools, only: pi
    implicit none
    integer :: i, j, k, Np
    real(WP) :: r, theta
    real(WP), dimension(3) :: c

    ! Allocate arrays
    allocate(sine(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    allocate(sqpsq(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    allocate(sqmsq(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    allocate(vfield(3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    allocate(rho(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    rho(:,:,:) = 1.0_WP

    ! Initialize LPT
    lp = lpt(cfg=cfg, name='LPT')
    call param_read('Number of particles', Np)
    lp%filter_width = cfg%min_meshsize
    lp%rho = 2.0_WP
    if (lp%cfg%amRoot) then
      call lp%resize(Np)
      r = 0.3_WP * max(lp%cfg%xL, lp%cfg%yL)
      c = (/ 0.5_WP * lp%cfg%xL, 0.5_WP * lp%cfg%yL, 0.5_WP * lp%cfg%zL /)
      do i = 1, Np
        ! Give id
        lp%p(i)%id = int(i,8)
        ! Set the diameter
        lp%p(i)%d = 0.04_WP
        ! Put them in a circle
        theta = 2 * pi * (i - 1) / Np
        lp%p(i)%pos = c + (/ r * cos(theta), r * sin(theta), 0.0_WP /)
        ! Locate the particle on the mesh
        lp%p(i)%ind = lp%cfg%get_ijk_global(lp%p(i)%pos, (/lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin/))
        ! Give zero velocity
        lp%p(i)%vel(:) = 0.0_WP
        ! Activate the particle
        lp%p(i)%flag=0
      end do
    else
      call lp%resize(0)
    end if
    call lp%sync()
    call lp%update_VF()

    ! Prepare arrays
    do k = cfg%kmino_, cfg%kmaxo_
      do j = cfg%jmino_, cfg%jmaxo_
        do i = cfg%imino_, cfg%imaxo_
          sine(i,j,k) = sin(2.0_WP * pi * i / cfg%nx) * sin(2.0_WP * pi * j / cfg%ny)
          sqpsq(i,j,k) = cfg%xm(i)**2 + cfg%ym(j)**2 + cfg%zm(k)
          sqmsq(i,j,k) = cfg%xm(i)**2 - cfg%ym(j)**2 + cfg%zm(k)
          vfield(:,i,j,k) = (/ sin(2.0_WP * pi * i / cfg%nx), cos(2.0_WP * pi * j / cfg%ny), 0.0_WP /)
        end do
      end do
    end do
    call cfg%sync(sine)
    call cfg%sync(sqpsq)
    call cfg%sync(sqmsq)
    call cfg%sync(vfield)

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
      ens_out=ensight(cfg=cfg,name='TEST')
      ! Add variables to output
      call ens_out%add_scalar('sine', sine)
      call ens_out%add_scalar('sqpsq', sqpsq)
      call ens_out%add_scalar('sqmsq', sqmsq)
      call ens_out%add_vector('vfield', vfield(1,:,:,:), vfield(2,:,:,:), vfield(3,:,:,:))
      call ens_out%add_particle('particles',pmesh)
    end block create_ensight

    ! Add Npy output
    create_npy: block
      ! Create Ensight output from cfg
      npy_out=npy(pg=cfg, folder='TEST')
      ! Add variables to output
      call npy_out%add_scalar('sine', sine)
      call npy_out%add_scalar('sqpsq', sqpsq)
      call npy_out%add_scalar('sqmsq', sqmsq)
      call npy_out%add_vector('vfield', vfield(1,:,:,:), vfield(2,:,:,:), vfield(3,:,:,:))
      call npy_out%add_particle('particles',pmesh)
    end block create_npy

  end subroutine simulation_init

  !> Write data
  subroutine simulation_run
    implicit none

    if (cfg%amRoot) write(*,*) "Writing files..."

    call update_pmesh()

    ! Output to ensight
    call ens_out%write_data(0.0_WP)
    if (cfg%rank .eq. 0) write(*,*) "Ensight output written"

    ! Output Npy
    call npy_out%write_data(0.0_WP)
    if (cfg%rank .eq. 0) write(*,*) "Npy output written"

  end subroutine simulation_run

  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    ! Deallocate work arrays
    !TODO

  end subroutine simulation_final

end module simulation
