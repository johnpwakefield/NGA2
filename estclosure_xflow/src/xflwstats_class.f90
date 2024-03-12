!> Statistics for Crossflow Domain
!>
!> John P Wakefield, September 2023
!>
!> This class needs to collect two different types of statistics; those of
!> the type collected by hitstats but without the benefits of periodicity and
!> simple statistics down the crossflow domain as a function of x.
!>
!> Note slice is in the yz plane
!>
!> For crossflow statistics, the convolution in the x direction is done
!> explicitly onto the 2d plane, then the convolution in the cross-wise
!> directions are done in real space.  For real values this is straightforward,
!> but for gradients it should be added that this means we need a separate
!> set of grids in which gradients of filtered quantities in the x-direction
!> may be explicitly projected, then spectral derivatives are taken in the y
!> and z directions.
!>
module xflwstats_class
  use precision,        only: WP
  use mpi_f08,          only: MPI_GROUP
  use cfourier_class,   only: cfourier
  use monitor_class,    only: monitor
  use string,           only: str_medium, str_long
  use sgrid_class,      only: sgrid
  use pgrid_class,      only: pgrid
  use config_class,     only: config
  use lpt_class,        only: lpt
  use npy_class,        only: npy
  use coupler_class,    only: coupler
  use lptcoupler_class, only: lptcoupler
  use partmesh_class,   only: partmesh
  use filterstats
  implicit none
  private


  public :: xflwstats


  !> microscopic statistics along domain (saved as monitor files)
  ! data only allocated / updated on root
  type :: xdepstats
    integer :: step
    integer, dimension(24) :: units
    ! these arrays are only imin to imax and don't include overlap
    integer, dimension(:), allocatable :: num
    real(WP), dimension(:), allocatable :: nden, ndensq, p, pke, fke
    real(WP), dimension(:,:), allocatable :: vel, velsq, slp, slpsq, drg, fvel  ! 3 by nx
    integer, dimension(:,:,:), allocatable :: tmpmesh
  contains
    procedure :: init => xdepstats_init
    procedure :: write => xdepstats_write
    final :: xdepstats_destroy
  end type xdepstats


  !> filtered fluctuation statistics
  type :: xplanestats
    character(len=str_medium) :: out_fname
    type(monitor) :: mon
    integer :: step
    real(WP) :: PB, PC2
    real(WP), dimension(3) :: UB, SB, PCUB, UBPCG  ! x, y, z
    real(WP), dimension(3,3) :: UC2   ! only upper triangle used
  contains
    procedure :: init => xplanestats_init
    final :: xplanestats_destroy
  end type xplanestats


  !> filter information
  type :: filter
    character(len=str_medium) :: name
    integer :: type
    real(WP), dimension(FLT_NUM_PARAMS) :: params
    ! we make this 3d but 1 cell wide in x so that we can use existing fft
    complex(WP), dimension(:,:,:), allocatable :: flt_f
  contains
    procedure :: init => filter_init
    final :: filter_destroy
  end type filter


  !> observation plane in xflow region
  ! each observation plane requires its own fft grids and data storage arrays
  type :: obsplane

    real(WP) :: offset

    class(config), pointer :: sim_cfg

    type(MPI_GROUP) :: grp
    logical :: in_grp
    type(pgrid) :: pg
    type(cfourier), pointer :: fft
    type(lptcoupler) :: lptcpl

    ! these can't be shared across obsplanes because the z indices have to
    ! match, but they can be shared across filters
    complex(WP), dimension(:,:,:), allocatable :: work_c
    real(WP), dimension(:,:,:), allocatable :: work_r
    ! even though the memory for these is the same for each obsplane, we have
    ! to re-project for each filter
    real(WP) :: phimean
    real(WP), dimension(:,:,:), pointer :: phi, phi_x
    real(WP), dimension(:,:,:,:), pointer :: up, uf

    ! io data, set up be xplanestats if needed
    character(len=str_medium) :: fn
    type(npy) :: io_f
    type(partmesh) :: io_pm
    real(WP), dimension(:,:,:,:), allocatable :: phi_io, phi_x_io
    real(WP), dimension(:,:,:,:,:), allocatable :: up_io, uf_io

  contains

    procedure :: init => obsplane_init
    procedure :: apply_filter => obsplane_apply_filter
    procedure :: destroy => obsplane_destroy

  end type obsplane


  !> top level xflow statistics object (interfaced with in simulation.f90)
  type :: xflwstats

    class(config), pointer :: sim_cfg

    real(WP), dimension(:,:,:), pointer :: vf, rhof, p, visc, U, V, W, zeros
    type(lpt), pointer :: ps

    integer :: num_filters, num_planes
    type(xdepstats) :: mstats
    type(obsplane), dimension(:), allocatable :: planes
    type(filter), dimension(:,:), allocatable :: filters      ! num_planes by num_filters
    type(xplanestats), dimension(:,:), allocatable :: fstats  ! num_planes by num_filters
    logical, dimension(:), allocatable :: sliceio_setup       ! num_planes

  contains

    procedure :: init => xflwstats_init
    procedure :: sliceio_init => xflwstats_sliceio_init
    procedure :: sliceio_destroy => xflwstats_sliceio_destroy
    procedure :: compute_xdep_stats => xflwstats_xdep_stats ! x-dependent stats
    procedure :: compute_stats => xflwstats_opln_stats      ! observation plane stats
    procedure :: destroy => xflwstats_destroy

  end type xflwstats


contains


  !> xdepstats (simple stats across xflow)

  subroutine xdepstats_init(this, g, mkdir)
    implicit none
    class(xdepstats), intent(inout) :: this
    class(pgrid), target, intent(in) :: g
    logical, intent(in) :: mkdir
    integer :: i

    allocate(this%num(g%imin:g%imax),                                         &
      this%nden(g%imin:g%imax), this%ndensq(g%imin:g%imax),                   &
      this%vel(3,g%imin:g%imax), this%velsq(3,g%imin:g%imax),                 &
      this%slp(3,g%imin:g%imax), this%slpsq(3,g%imin:g%imax),                 &
      this%drg(3,g%imin:g%imax), this%p(g%imin:g%imax),                       &
      this%fvel(3,g%imin:g%imax),                                             &
      this%pke(g%imin:g%imax), this%fke(g%imin:g%imax))
    allocate(this%tmpmesh(g%imin:g%imax,g%jmin:g%jmax,g%kmin:g%kmax))

    ! just use raw files for all the x dependent statistics
    if (g%rank .eq. 0) then
      if (mkdir) call execute_command_line('mkdir -p xflw_stat')
      open(newunit=this%units( 1), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/num')
      open(newunit=this%units( 2), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/nden')
      open(newunit=this%units( 3), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/ndensq')
      open(newunit=this%units( 4), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/velx')
      open(newunit=this%units( 5), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/vely')
      open(newunit=this%units( 6), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/velz')
      open(newunit=this%units( 7), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/velsqx')
      open(newunit=this%units( 8), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/velsqy')
      open(newunit=this%units( 9), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/velsqz')
      open(newunit=this%units(10), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/slpx')
      open(newunit=this%units(11), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/slpy')
      open(newunit=this%units(12), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/slpz')
      open(newunit=this%units(13), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/slpsqx')
      open(newunit=this%units(14), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/slpsqy')
      open(newunit=this%units(15), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/slpsqz')
      open(newunit=this%units(16), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/drgx')
      open(newunit=this%units(17), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/drgy')
      open(newunit=this%units(18), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/drgz')
      open(newunit=this%units(19), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/press')
      open(newunit=this%units(20), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/fvelx')
      open(newunit=this%units(21), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/fvely')
      open(newunit=this%units(22), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/fvelz')
      open(newunit=this%units(23), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/pke')
      open(newunit=this%units(24), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/fke')
      !TODO split this so that x and xm are used appropriately
      do i = 1, 24
        write(this%units(i),'(A8,*(E20.6))') ' x ', g%x(g%imin:g%imax)
      end do
    end if

  end subroutine xdepstats_init

  subroutine xdepstats_write(this, g)
    implicit none
    class(xdepstats), intent(inout) :: this
    class(pgrid), intent(in) :: g
    integer :: n

    if (g%rank .eq. 0) then

      n = this%step

      write(this%units( 1),'(I8,*(I20))')   n, this%num(       g%imin:g%imax)
      write(this%units( 2),'(I8,*(E20.8E3))') n, this%nden(    g%imin:g%imax)
      write(this%units( 3),'(I8,*(E20.8E3))') n, this%ndensq(  g%imin:g%imax)
      write(this%units( 4),'(I8,*(E20.8E3))') n, this%vel(1,   g%imin:g%imax)
      write(this%units( 5),'(I8,*(E20.8E3))') n, this%vel(2,   g%imin:g%imax)
      write(this%units( 6),'(I8,*(E20.8E3))') n, this%vel(3,   g%imin:g%imax)
      write(this%units( 7),'(I8,*(E20.8E3))') n, this%velsq(1, g%imin:g%imax)
      write(this%units( 8),'(I8,*(E20.8E3))') n, this%velsq(2, g%imin:g%imax)
      write(this%units( 9),'(I8,*(E20.8E3))') n, this%velsq(3, g%imin:g%imax)
      write(this%units(10),'(I8,*(E20.8E3))') n, this%slp(1,   g%imin:g%imax)
      write(this%units(11),'(I8,*(E20.8E3))') n, this%slp(2,   g%imin:g%imax)
      write(this%units(12),'(I8,*(E20.8E3))') n, this%slp(3,   g%imin:g%imax)
      write(this%units(13),'(I8,*(E20.8E3))') n, this%slpsq(1, g%imin:g%imax)
      write(this%units(14),'(I8,*(E20.8E3))') n, this%slpsq(2, g%imin:g%imax)
      write(this%units(15),'(I8,*(E20.8E3))') n, this%slpsq(3, g%imin:g%imax)
      write(this%units(16),'(I8,*(E20.8E3))') n, this%drg(1,   g%imin:g%imax)
      write(this%units(17),'(I8,*(E20.8E3))') n, this%drg(2,   g%imin:g%imax)
      write(this%units(18),'(I8,*(E20.8E3))') n, this%drg(3,   g%imin:g%imax)
      write(this%units(19),'(I8,*(E20.8E3))') n, this%p(       g%imin:g%imax)
      write(this%units(20),'(I8,*(E20.8E3))') n, this%fvel(1,  g%imin:g%imax)
      write(this%units(21),'(I8,*(E20.8E3))') n, this%fvel(2,  g%imin:g%imax)
      write(this%units(22),'(I8,*(E20.8E3))') n, this%fvel(3,  g%imin:g%imax)
      write(this%units(23),'(I8,*(E20.8E3))') n, this%pke(     g%imin:g%imax)
      write(this%units(24),'(I8,*(E20.8E3))') n, this%fke(     g%imin:g%imax)

    end if

  end subroutine xdepstats_write

  subroutine xdepstats_destroy(this)
    implicit none
    type(xdepstats), intent(inout) :: this
    integer :: i

    ! monitor falls out of scope and has its own destructor

    do i = 1, 24
      close(unit=this%units(i))
    end do

    deallocate(this%num, this%nden, this%ndensq, this%vel, this%velsq,        &
      this%slp, this%slpsq, this%drg, this%p, this%fvel, this%pke, this%fke)

  end subroutine xdepstats_destroy


  !> observation plane statistics object

  subroutine xplanestats_init(this, amroot)
    implicit none
    class(xplanestats), intent(inout) :: this
    logical, intent(in) :: amroot

    this%mon = monitor(amroot, this%out_fname)
    call this%mon%add_column(this%step,     'step')
    call this%mon%add_column(this%PB,       'PB')
    call this%mon%add_column(this%PC2,      'PC2')
    call this%mon%add_column(this%UB(1),    'UBx')
    call this%mon%add_column(this%UB(2),    'UBy')
    call this%mon%add_column(this%UB(3),    'UBz')
    call this%mon%add_column(this%SB(1),    'SBx')
    call this%mon%add_column(this%SB(2),    'SBy')
    call this%mon%add_column(this%SB(3),    'SBz')
    call this%mon%add_column(this%PCUB(1),  'PCUBx')
    call this%mon%add_column(this%PCUB(2),  'PCUBy')
    call this%mon%add_column(this%PCUB(3),  'PCUBz')
    call this%mon%add_column(this%UBPCG(1), 'UBPCGx')
    call this%mon%add_column(this%UBPCG(2), 'UBPCGy')
    call this%mon%add_column(this%UBPCG(3), 'UBPCGz')
    call this%mon%add_column(this%UC2(1,1), 'UC2xx')
    call this%mon%add_column(this%UC2(1,2), 'UC2xy')
    call this%mon%add_column(this%UC2(1,3), 'UC2xz')
    call this%mon%add_column(this%UC2(2,2), 'UC2yy')
    call this%mon%add_column(this%UC2(2,3), 'UC2yz')
    call this%mon%add_column(this%UC2(3,3), 'UC2zz')

  end subroutine xplanestats_init

  subroutine xplanestats_destroy(this)
    implicit none
    type(xplanestats), intent(inout) :: this

    ! nothing to do (monitor falls out of scope and has its own destructor)

  end subroutine xplanestats_destroy


  !> filter object

  ! assumes name, params, and type are already set
  subroutine filter_init(this, fft)
    use mpi_f08, only:   mpi_allreduce, MPI_SUM
    use parallel, only:  MPI_REAL_WP
    use mathtools, only: pi
    use messager, only: die
    implicit none
    class(filter), intent(inout) :: this
    class(cfourier), intent(inout) :: fft
    real(WP) :: flt_int_l, flt_int_g

    ! allocate Fourier space filter
    allocate(this%flt_f(fft%pg%imin_:fft%pg%imin_,fft%pg%jmin_:fft%pg%jmax_,  &
      fft%pg%kmin_:fft%pg%kmax_))   !       ^^^^^ not a typo

    ! fill real filter; we build fft%pg such that the corner is at origin
    select case (this%type)
    case (FLT_GAUSSIAN)
      init_gaussian: block
        real(WP) :: b, r2_z, r2_y, yoff, zoff
        integer :: j, k, jshift, kshift
        this%flt_f(:,:,:) = (0.0_WP, 0.0_WP)
        b = -6.0_WP / this%params(1)**2
        do jshift = -2, 2
          yoff = jshift * fft%pg%yL
          do kshift = -2, 2
            zoff = kshift * fft%pg%zL
            do k = fft%pg%kmin_, fft%pg%kmax_
              r2_z = (fft%pg%zm(k) + zoff)**2
              do j = fft%pg%jmin_, fft%pg%jmax_
                r2_y = (fft%pg%ym(j) + yoff)**2
                this%flt_f(fft%pg%imin_,j,k) = this%flt_f(fft%pg%imin_,j,k) + &
                  exp(b * (r2_z + r2_y))
              end do
            end do
          end do
        end do
      end block init_gaussian
    case default
      call die("[xflwstats] invalid filter type")
    end select

    ! normalize filter
    flt_int_l = sum(realpart(this%flt_f))
    call mpi_allreduce(flt_int_l, flt_int_g, 1, MPI_REAL_WP, MPI_SUM,         &
      fft%pg%comm)
    this%flt_f = this%flt_f / flt_int_g

    ! transform filter
    call fft%ytransform_forward(this%flt_f)
    call fft%ztransform_forward(this%flt_f)
    ! this shouldn't be necessary since we normalized in real space, but
    ! whatever
    if (fft%oddball) this%flt_f(fft%pg%imin_,fft%pg%jmin_,fft%pg%kmin_) =     &
      (1.0_WP, 0.0_WP)

  end subroutine filter_init

  subroutine filter_destroy(this)
    implicit none
    type(filter), intent(inout) :: this

    deallocate(this%flt_f)

  end subroutine filter_destroy


  !> observation plane

  subroutine build_slice_pg(sim_pg, ny, nz, offset, width, obs_pg, obs_grp, in_grp)
    use mpi_f08, only: mpi_allreduce, MPI_LOGICAL, MPI_LOR, mpi_group_range_incl
    use sgrid_class, only: cartesian
    use messager, only: die
    implicit none
    class(pgrid), intent(in) :: sim_pg
    integer :: ny, nz
    real(WP), intent(in) :: offset, width
    type(pgrid), intent(out) :: obs_pg
    type(MPI_GROUP), intent(out) :: obs_grp
    logical, intent(out) :: in_grp
    real(WP) :: dy, dz
    integer, dimension(3,sim_pg%nproc) :: ranges
    logical, dimension(sim_pg%nproc) :: in_group_loc, in_group
    type(sgrid) :: obs_sg
    integer :: i, j, k, ierr
    integer, dimension(3) :: ind
    logical :: prev
    real(WP), dimension(:), allocatable :: xs, ys, zs

    ! allocate xs, ys, zs arrays
    allocate(xs(2), ys(ny+1), zs(nz+1))

    ! build sgrid
    xs(1) = offset - 0.5_WP * width
    xs(2) = offset + 0.5_WP * width
    dy = sim_pg%yL / ny; dz = sim_pg%zL / nz;
    ys(:) = (/ (real(i,WP) * dy + sim_pg%y(sim_pg%jmin), i = 0, ny) /)
    zs(:) = (/ (real(i,WP) * dz + sim_pg%z(sim_pg%kmin), i = 0, nz) /)
    obs_sg = sgrid(coord=cartesian, no=0, x=xs, y=ys, z=zs,                   &
      xper=.false., yper=.true., zper=.true., name='XF_OBS_SG')

    ! find middle process
    ind(:) = sim_pg%get_ijk_global((/ offset, ys(2), zs(2) /))
    in_group_loc(:) = .false.
    in_group_loc(sim_pg%rank+1) = sim_pg%imin_ .le. ind(1)                    &
      .and. ind(1) .le. sim_pg%imax_
    call mpi_allreduce(in_group_loc, in_group, sim_pg%nproc,                  &
      MPI_LOGICAL, MPI_LOR, sim_pg%comm, ierr)

    ! build the processor group
    if (in_group(1)) then; k = 1; else; k = 0; end if;
    do i = 2, sim_pg%nproc
      if (in_group(i) .and. .not. in_group(i-1)) k = k + 1
    end do
    j = 0
    prev = .false.
    do i = 1, sim_pg%nproc
      if (in_group(i) .and. .not. prev) then
        j = j + 1
        ranges(1,j) = i - 1
      end if
      if (.not. in_group(i) .and. prev) then
        ranges(2,j) = i - 2
      end if
      prev = in_group(i)
    end do
    if (prev) ranges(2,j) = sim_pg%nproc - 1
    ranges(3,:) = 1
    if (j .ne. k) call die("[XS] error determining obs plane group")
    call mpi_group_range_incl(sim_pg%group, k, ranges(1:3,1:k), obs_grp, ierr)
    in_grp = in_group(sim_pg%rank+1)

    ! free memory
    deallocate(xs, ys, zs)

    ! set up pgrid
    if (in_grp) obs_pg = pgrid(obs_sg, obs_grp, (/ 1, sim_pg%npy, sim_pg%npz /))

  end subroutine build_slice_pg

  subroutine obsplane_init(this, sim_cfg, FFTN, ps, offset, width)
    use messager, only: die
    use mpi_f08, only: mpi_bcast, mpi_allreduce, MPI_LOGICAL, MPI_LOR, mpi_group_range_incl
    use sgrid_class, only: cartesian
    implicit none
    class(obsplane), intent(inout) :: this
    class(config), target, intent(in) :: sim_cfg
    integer, dimension(3), intent(in) :: FFTN
    class(lpt), target, intent(in) :: ps            ! lpt with xflow particles
    real(WP), intent(in) :: offset, width
    character(len=10) :: offsetstr
    integer :: li, hi, lj, hj, lk, hk

    ! set basic properties
    this%sim_cfg => sim_cfg
    this%offset = offset
    write(offsetstr, '(F10.6)') offset
    this%fn = "obsplane_"//trim(adjustl(offsetstr))

    ! check one cell in x
    if (FFTN(1) .ne. 1) call die("[XFS] FFTN(1) must be 1 for xflwstats.")

    ! build obs slice
    call build_slice_pg(sim_cfg, FFTN(2), FFTN(3), offset, width, this%pg,    &
      this%grp, this%in_grp)

    ! set up fft
    if (this%in_grp) allocate(this%fft, source=cfourier(this%pg))

    ! setup particle holder for projection onto plane and particle coupler
    this%lptcpl = lptcoupler(src_grp=sim_cfg%group, dst_grp=this%grp,         &
      name='OBS_LPT_CPL')
    call this%lptcpl%set_src(ps)
    if (this%in_grp) call this%lptcpl%set_dst(this%pg)
    call this%lptcpl%initialize()

    ! allocate arrays
    if (this%in_grp) then
      li = this%pg%imino_; lj = this%pg%jmino_; lk = this%pg%kmino_;
      hi = this%pg%imaxo_; hj = this%pg%jmaxo_; hk = this%pg%kmaxo_;
      allocate(this%work_c(li:hi,lj:hj,lk:hk), this%work_r(li:hi,lj:hj,lk:hk),&
        this%phi(li:hi,lj:hj,lk:hk), this%phi_x(li:hi,lj:hj,lk:hk),           &
        this%up(li:hi,lj:hj,lk:hk,1:3), this%uf(li:hi,lj:hj,lk:hk,1:3))
    end if

  end subroutine obsplane_init

  ! assumes obs_lpt is up to date
  ! U, V, W, Ux, Vx, Wx are from flow solver and on flow solver pgrid
  subroutine obsplane_apply_filter(this, flt, ps, U, V, W)
    use mpi_f08, only: mpi_allreduce, MPI_SUM
    use parallel, only: MPI_REAL_WP
    use lpt_class, only: part
    use mathtools, only: pi
    implicit none
    class(obsplane), intent(inout) :: this
    type(filter), intent(in) :: flt
    type(part), dimension(:), intent(inout) :: ps
    real(WP), dimension(this%sim_cfg%imin_:,this%sim_cfg%jmin_:,              &
      this%sim_cfg%kmin_:), intent(in) :: U, V, W
    integer :: n, i, j, k
    integer, dimension(3) :: ind
    real(WP), dimension(3) :: fvel
    real(WP) :: b, cellvoli, dx, vp, vpexpbdx2

    if (.not. this%in_grp) return

    ! project particles
    this%phi(:,:,:) = 0.0_WP; this%phi_x(:,:,:) = 0.0_WP;
    this%up(:,:,:,:) = 0.0_WP; this%uf(:,:,:,:) = 0.0_WP;
    do n = 1, size(ps)
      vp = pi * ps(n)%d**3 / 6.0_WP
      ps(n)%ind = this%sim_cfg%get_ijk_global(ps(n)%pos)
      fvel = this%sim_cfg%get_velocity(pos=ps(n)%pos, i0=ps(n)%ind(1),        &
        j0=ps(n)%ind(2), k0=ps(n)%ind(3), U=U, V=V, W=W)
      ind(:) = this%pg%get_ijk_global(ps(n)%pos)
      i = this%pg%imin_; j = ind(2); k = ind(3);
      dx = ps(n)%pos(1) - this%offset
      b = -6.0_WP / flt%params(1)**2
      vpexpbdx2 = vp * exp(b * dx**2)
      this%phi(i,j,k) = this%phi(i,j,k) + vpexpbdx2
      this%phi_x(i,j,k) = this%phi_x(i,j,k) + 2 * b * dx * vpexpbdx2
      this%up(i,j,k,:) = this%up(i,j,k,:) + ps(n)%vel * vpexpbdx2
      this%uf(i,j,k,:) = this%uf(i,j,k,:) + fvel * vpexpbdx2
    end do

    ! normalize phi
    cellvoli = this%pg%nx * this%pg%ny * this%pg%nz / this%pg%vol_total
    this%phi(:,:,:) = cellvoli * this%phi(:,:,:)
    b = sum(this%phi(this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_, &
      this%pg%kmin_:this%pg%kmax_))
    call mpi_allreduce(b, this%phimean, 1, MPI_REAL_WP, MPI_SUM, this%pg%comm)
    this%phimean = this%phimean / (this%pg%nx * this%pg%ny * this%pg%nz)

    ! we normalize everything by phimean to avoid roundoff error
    this%phi(:,:,:) = this%phi(:,:,:) / this%phimean
    this%phi_x(:,:,:) = this%phi_x(:,:,:) / this%phimean
    this%up(:,:,:,:) = this%up(:,:,:,:) / this%phimean
    this%uf(:,:,:,:) = this%uf(:,:,:,:) / this%phimean
    this%phi(:,:,:) = this%phi(:,:,:) - 1.0

    ! filter in normal directions
    this%work_c(:,:,:) = this%phi(:,:,:)
    call this%fft%ytransform_forward(this%work_c(:,:,:))
    call this%fft%ztransform_forward(this%work_c(:,:,:))
    this%work_c(:,:,:) = this%work_c(:,:,:) * flt%flt_f(:,:,:)
    call this%fft%ztransform_backward(this%work_c(:,:,:))
    call this%fft%ytransform_backward(this%work_c(:,:,:))
    this%phi(:,:,:) = real(this%work_c(:,:,:), WP)
    this%work_c(:,:,:) = cellvoli * this%phi_x(:,:,:)
    call this%fft%ytransform_forward(this%work_c(:,:,:))
    call this%fft%ztransform_forward(this%work_c(:,:,:))
    this%work_c(:,:,:) = this%work_c(:,:,:) * flt%flt_f(:,:,:)
    call this%fft%ztransform_backward(this%work_c(:,:,:))
    call this%fft%ytransform_backward(this%work_c(:,:,:))
    this%phi_x(:,:,:) = real(this%work_c(:,:,:), WP)
    do n = 1, 3
      this%work_c(:,:,:) = cellvoli * this%up(:,:,:,n)
      call this%fft%ytransform_forward(this%work_c(:,:,:))
      call this%fft%ztransform_forward(this%work_c(:,:,:))
      this%work_c(:,:,:) = this%work_c(:,:,:) * flt%flt_f(:,:,:)
      call this%fft%ztransform_backward(this%work_c(:,:,:))
      call this%fft%ytransform_backward(this%work_c(:,:,:))
      this%up(:,:,:,n) = real(this%work_c(:,:,:), WP)
      this%work_c(:,:,:) = cellvoli * this%uf(:,:,:,n)
      call this%fft%ytransform_forward(this%work_c(:,:,:))
      call this%fft%ztransform_forward(this%work_c(:,:,:))
      this%work_c(:,:,:) = this%work_c(:,:,:) * flt%flt_f(:,:,:)
      call this%fft%ztransform_backward(this%work_c(:,:,:))
      call this%fft%ytransform_backward(this%work_c(:,:,:))
      this%uf(:,:,:,n) = real(this%work_c(:,:,:), WP)
    end do

  end subroutine obsplane_apply_filter

  subroutine obsplane_destroy(this)
    implicit none
    class(obsplane), intent(inout) :: this

    deallocate(this%fft, this%work_c, this%work_r, this%phi, this%phi_x,      &
      this%up, this%uf)

  end subroutine obsplane_destroy


  !> xflwstats

  subroutine xflwstats_init(this, sim_cfg, filterfile, FFTN, vf, rhof, p, visc, U, V, W, ps)
    use param,       only: param_read
    use mpi_f08,     only: mpi_bcast, MPI_REAL, MPI_INTEGER, MPI_CHARACTER,   &
      MPI_LOGICAL, MPI_COMM_WORLD
    use parallel,    only: MPI_REAL_WP
    use sgrid_class, only: cartesian
    use messager,    only: die, log
    implicit none
    class(xflwstats), intent(inout) :: this
    class(config), target, intent(in) :: sim_cfg
    character(len=str_medium), intent(in) :: filterfile
    integer, dimension(2), intent(in) :: FFTN
    type(lpt), target, intent(in) :: ps
    real(WP), dimension(sim_cfg%imino_:,sim_cfg%jmino_:,sim_cfg%kmino_:),     &
      intent(in), target :: vf, rhof, p, visc, U, V, W
    real(WP), dimension(:), allocatable :: planelocs, lptwidths
    type(filter_info_row) :: f_info_raw
    character(len=8) :: istr, nstr, ostr
    integer :: i, j, fh, ierr
    logical :: use_slice_io

    ! store pointer to simulation config and data
    this%sim_cfg => sim_cfg; this%ps => ps;
    this%vf => vf; this%rhof => rhof; this%p => p; this%visc => visc;
    this%U => U; this%V => V; this%W => W;

    ! allocate array of zeros
    allocate(this%zeros(sim_cfg%imino_:sim_cfg%imaxo_,                        &
      sim_cfg%jmino_:sim_cfg%jmaxo_,sim_cfg%kmino_:sim_cfg%kmaxo_))
    this%zeros(:,:,:) = 0.0_WP

    ! set up microscopic statistics
    !TODO set flag somewhere else for mkdir arg
    call this%mstats%init(ps%cfg, .true.)

    ! read observation plane info
    if (sim_cfg%rank .eq. 0) then
      call param_read('HS obs plane count', this%num_planes)
      allocate(planelocs(this%num_planes),lptwidths(this%num_planes))
      call param_read('HS obs plane locations', planelocs)
      call param_read('HS obs plane lpt widths', lptwidths)
    end if
    call mpi_bcast(this%num_planes, 1, MPI_INTEGER, 0, sim_cfg%comm, ierr)
    if (sim_cfg%rank .ne. 0)                                                  &
      allocate(planelocs(this%num_planes),lptwidths(this%num_planes))
    call mpi_bcast(planelocs, this%num_planes, MPI_REAL_WP, 0, sim_cfg%comm, ierr)
    call mpi_bcast(lptwidths, this%num_planes, MPI_REAL_WP, 0, sim_cfg%comm, ierr)

    ! read filter info (as root)
    if (sim_cfg%rank .eq. 0) then
      call param_read('HS write slice', use_slice_io)
      open(newunit=fh, file=filterfile, action='READ', access='SEQUENTIAL')
      this%num_filters = 0
      do
        read(fh,*,iostat=ierr) f_info_raw
        if (ierr .ne. 0) exit
        this%num_filters = this%num_filters + 1
      end do
      allocate(this%filters(this%num_planes,this%num_filters))
      rewind(fh)
      do i = 1, this%num_filters
        read(fh,*,iostat=ierr) f_info_raw
        if (ierr .ne. 0) call die("[EC] error reading filterfile row")
        do j = 1, this%num_planes
          select case (f_info_raw%typename)
          case ('gaussian', 'GAUSSIAN', 'g', 'G')
            this%filters(j,i)%type = FLT_GAUSSIAN
          case ('box', 'BOX', 'b', 'B')
            this%filters(j,i)%type = FLT_BOX
          case default
            call die("[EC] unknown filter type '"//f_info_raw%typename//"'")
          end select
          this%filters(j,i)%name = 'xflw_'//trim(adjustl(f_info_raw%out_fname))
          this%filters(j,i)%params(:) = f_info_raw%params(:)
        end do
      end do
      close(fh)
    end if

    ! broadcast filter info
    call mpi_bcast(this%num_filters, 1, MPI_INTEGER, 0, sim_cfg%comm, ierr)
    if (sim_cfg%rank .ne. 0) allocate(this%filters(this%num_planes,this%num_filters))
    do j = 1, this%num_planes
      do i = 1, this%num_filters
        call mpi_bcast(this%filters(j,i)%name, str_medium, MPI_CHARACTER, 0,  &
          sim_cfg%comm, ierr)
        this%filters(j,i)%name = trim(this%filters(j,i)%name)
        call mpi_bcast(this%filters(j,i)%type, 1, MPI_INTEGER, 0,             &
          sim_cfg%comm, ierr)
        call mpi_bcast(this%filters(j,i)%params, FLT_NUM_PARAMS, MPI_REAL_WP, &
          0, sim_cfg%comm, ierr)
        !TODO check if this is synced at init, else sync it here
        !call mpi_bcast(this%filters(j,i)%use_slice_io, 1, MPI_LOGICAL, 0,         &
        !  sim_cfg%comm, ierr)
      end do
    end do

    ! init observation planes
    allocate(this%planes(this%num_planes))
    write(nstr, '(I8)') this%num_planes
    do i = 1, this%num_planes
      call this%planes(i)%init(sim_cfg, (/ 1, FFTN(1), FFTN(2) /), ps,        &
        planelocs(i), lptwidths(i))
      write(istr, '(I8)') i
      if (sim_cfg%rank .eq. 0) call log(                                      &
        "Initialized XF observation plane "//istr//" of "//nstr//".")
    end do

    ! init filters
    if (this%num_planes .eq. 0 .and. this%num_filters .gt. 0) call die(       &
      "Must have at least one observation plane if any filters are used.")
    write(nstr, '(I8)') this%num_filters
    do i = 1, this%num_filters
      do j = 1, this%num_planes
        if (this%planes(j)%in_grp)                                            &
          call this%filters(j,i)%init(this%planes(j)%fft)
      end do
      write(istr, '(I8)') i
      if (sim_cfg%rank .eq. 0) call log(                                      &
        "Initialized XF filter "//istr//" of "//nstr//".")
    end do

    ! initialize fstats objects
    allocate(this%fstats(this%num_planes,this%num_filters))
    do i = 1, this%num_filters
      do j = 1, this%num_planes
        write(ostr, '(F8.4)') planelocs(j)
        write(nstr, '(F8.4)') this%filters(j,i)%params(1)
        this%fstats(j,i)%out_fname = 'xflw_obs_'//trim(adjustl(ostr))//'_'//trim(adjustl(nstr))
        if (this%planes(j)%in_grp)                                            &
          call this%fstats(j,i)%init(this%planes(j)%fft%pg%rank .eq. 0)
      end do
    end do

    ! allocate sliceio setup array
    allocate(this%sliceio_setup(this%num_planes))
    this%sliceio_setup(:) = .false.

    ! cleanup
    deallocate(planelocs)

  end subroutine xflwstats_init

  subroutine xflwstats_sliceio_init(this)
    implicit none
    class(xflwstats), intent(inout) :: this
    real(WP), dimension(:,:,:), pointer :: vx, vy, vz
    character(len=8) :: offstr
    integer :: m, n, li, hi, lj, hj, lk, hk

    do m = 1, this%num_planes

      if (.not. this%planes(m)%in_grp) cycle

      write(offstr, '(F8.6)') this%planes(m)%offset

      ! set up partmesh
      this%planes(m)%io_pm = partmesh(nvar=2, nvec=2, name='lpt')
      this%planes(m)%io_pm%varname(1) = "id"
      this%planes(m)%io_pm%varname(2) = "dp"
      this%planes(m)%io_pm%vecname(1) = "vel"
      this%planes(m)%io_pm%vecname(2) = "fld"

      ! set up npy writer
      this%planes(m)%io_f = npy(pg=this%planes(m)%pg, folder=this%planes(m)%fn)
      do n = 1, this%num_filters
        call this%planes(m)%io_f%add_scalar(                                  &
          trim(this%filters(m,n)%name)//"_"//offstr//"_phi",                  &
          this%planes(m)%phi)
        call this%planes(m)%io_f%add_scalar(                                  &
          trim(this%filters(m,n)%name)//"_"//offstr//"_phi_x",                &
          this%planes(m)%phi_x)
        vx => this%planes(m)%up(:,:,:,1)
        vy => this%planes(m)%up(:,:,:,2)
        vz => this%planes(m)%up(:,:,:,3)
        call this%planes(m)%io_f%add_vector(                                  &
          trim(this%filters(m,n)%name)//"_"//offstr//"_up", vx, vy, vz)
        vx => this%planes(m)%uf(:,:,:,1)
        vy => this%planes(m)%uf(:,:,:,2)
        vz => this%planes(m)%uf(:,:,:,3)
        call this%planes(m)%io_f%add_vector(                                  &
          trim(this%filters(m,n)%name)//"_"//offstr//"_uf", vx, vy, vz)
      end do
      call this%planes(m)%io_f%add_particle('partslice', this%planes(m)%io_pm)

      li = this%planes(m)%pg%imino_; hi = this%planes(m)%pg%imaxo_;
      lj = this%planes(m)%pg%jmino_; hj = this%planes(m)%pg%jmaxo_;
      lk = this%planes(m)%pg%kmino_; hk = this%planes(m)%pg%kmaxo_;
      allocate(this%planes(m)%phi_io(li:hi,lj:hj,lk:hk,1:this%num_filters),   &
             this%planes(m)%phi_x_io(li:hi,lj:hj,lk:hk,1:this%num_filters),   &
             this%planes(m)%up_io(li:hi,lj:hj,lk:hk,1:3,1:this%num_filters),  &
             this%planes(m)%uf_io(li:hi,lj:hj,lk:hk,1:3,1:this%num_filters))

      this%sliceio_setup(m) = .true.

    end do

  end subroutine xflwstats_sliceio_init

  subroutine xflwstats_sliceio_destroy(this)
    implicit none
    class(xflwstats), intent(inout) :: this
    integer :: m

    do m = 1, this%num_planes

      if (.not. this%planes(m)%in_grp) cycle

      !TODO call npy destructor
      !TODO call partmesh destructor

      ! deallocate io arrays
      deallocate(this%planes(m)%phi_io, this%planes(m)%phi_x_io,              &
        this%planes(m)%up_io, this%planes(m)%uf_io)

    end do

    this%sliceio_setup(:) = .false.

  end subroutine xflwstats_sliceio_destroy

  ! computes microscopic statistics along xflow region
  subroutine xflwstats_xdep_stats(this, step)
    use mathtools, only:  pi
    use mpi_f08, only:    mpi_reduce, MPI_SUM, MPI_INTEGER
    use parallel, only:   MPI_REAL_WP
    use messager, only:   die
    implicit none
    class(xflwstats), intent(inout) :: this
    integer, intent(in) :: step
    integer :: n, i, j, k
    real(WP) :: opt_dt, slicevol, mp
    integer, dimension(3) :: indmin, indmax
    real(WP), dimension(3) :: acc, fvel, velmean
    integer, dimension(this%ps%cfg%imin:this%ps%cfg%imax) :: num
    real(WP), dimension(this%ps%cfg%imin:this%ps%cfg%imax) :: press, pke, fke
    real(WP), dimension(3,this%ps%cfg%imin:this%ps%cfg%imax) :: vel, velsq,   &
      slp, slpsq, drg

    indmin(:) = (/ this%ps%cfg%imin, this%ps%cfg%jmin, this%ps%cfg%kmin /)
    indmax(:) = (/ this%ps%cfg%imax, this%ps%cfg%jmax, this%ps%cfg%kmax /)

    ! store current step
    this%mstats%step = step

    ! compute mean quantities
    call this%ps%cfg%integrate(this%U, velmean(1))
    call this%ps%cfg%integrate(this%V, velmean(2))
    call this%ps%cfg%integrate(this%W, velmean(3))
    velmean(:) = velmean(:) / this%ps%cfg%vol_total

    ! zero out arrays
    num(:) = 0; this%mstats%tmpmesh(:,:,:) = 0.0_WP;
    vel(:,:) = 0.0_WP; velsq(:,:) = 0.0_WP;
    slp(:,:) = 0.0_WP; slpsq(:,:) = 0.0_WP;
    drg(:,:) = 0.0_WP; pke(:) = 0.0_WP;

    ! compute microscopic statistics
    do n = 1, this%ps%np_
      this%ps%p(n)%ind = this%ps%cfg%get_ijk_global(this%ps%p(n)%pos,         &
        this%ps%p(n)%ind)
      i = this%ps%p(n)%ind(1)
      j = this%ps%p(n)%ind(2)
      k = this%ps%p(n)%ind(3)
      num(i) = num(i) + 1
      this%mstats%tmpmesh(i,j,k) = this%mstats%tmpmesh(i,j,k) + 1
      !vf(i) = vf(i) + pi * this%ps%p(n)%d**3 / 6
      vel(:,i) = vel(:,i) + (this%ps%p(n)%vel(:) - velmean(:))
      velsq(:,i) = velsq(:,i) + (this%ps%p(n)%vel(:) - velmean(:))**2
      fvel = this%ps%cfg%get_velocity(pos=this%ps%p(n)%pos, i0=i, j0=j, k0=k, &
        U=this%U, V=this%V, W=this%W)
      slp(:,i) = slp(:,i) + fvel - this%ps%p(n)%vel
      slpsq(:,i) = slpsq(:,i) + (fvel - this%ps%p(n)%vel)**2
      !TODO check this is the right viscosity
      call this%ps%get_rhs(U=this%U, V=this%V, W=this%W, rho=this%rhof,       &
        visc=this%visc, stress_x=this%zeros, stress_y=this%zeros,             &
        stress_z=this%zeros, p=this%ps%p(n), acc=acc, opt_dt=opt_dt)
      drg(:,i) = drg(:,i) + acc(:)
      mp = pi * this%ps%p(n)%d**3 / 6.0_WP
      pke(i) = pke(i) + 0.5_WP * mp * sum((this%ps%p(n)%vel - velmean)**2)
    end do

    ! collect statistics to root
    call mpi_reduce(num,   this%mstats%num,   this%ps%cfg%nx,   MPI_INTEGER,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(vel,   this%mstats%vel,   3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(velsq, this%mstats%velsq, 3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(slp,   this%mstats%slp,   3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(slpsq, this%mstats%slpsq, 3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(drg,   this%mstats%drg,   3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(pke,   this%mstats%pke,   this%ps%cfg%nx,   MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)

    ! divide as appropriate
    if (this%ps%cfg%rank .eq. 0) then
      do i = this%ps%cfg%imin, this%ps%cfg%imax
        slicevol = this%ps%cfg%dx(i) * this%ps%cfg%yL * this%ps%cfg%zL
        this%mstats%nden(i) = this%mstats%num(i) / slicevol
        this%mstats%ndensq(i) = 0.0_WP
        do k = this%ps%cfg%kmin, this%ps%cfg%kmax
          do j = this%ps%cfg%jmin, this%ps%cfg%jmax
            this%mstats%ndensq(i) = this%mstats%ndensq(i)                     &
              + real(this%mstats%tmpmesh(i,j,k)**2,WP)
          end do
        end do
        this%mstats%ndensq(i) = this%mstats%ndensq(i) / slicevol
        if (this%mstats%num(i) .ne. 0) then
          this%mstats%vel(:,i)   = this%mstats%vel(:,i)   / this%mstats%num(i)
          this%mstats%velsq(:,i) = this%mstats%velsq(:,i) / this%mstats%num(i)
          this%mstats%slp(:,i)   = this%mstats%slp(:,i)   / this%mstats%num(i)
          this%mstats%slpsq(:,i) = this%mstats%slpsq(:,i) / this%mstats%num(i)
          this%mstats%drg(:,i)   = this%mstats%drg(:,i)   / this%mstats%num(i)
        end if
      end do
    end if

    ! compute continuum stats along x direction
    vel(:,:) = 0.0_WP; press(:) = 0.0_WP; fke(:) = 0.0_WP;
    do k = this%ps%cfg%kmin_, this%ps%cfg%kmax_
      do j = this%ps%cfg%jmin_, this%ps%cfg%jmax_
        do i = this%ps%cfg%imin_, this%ps%cfg%imax_
          vel(1,i) = vel(1,i) + this%U(i,j,k)
          vel(2,i) = vel(2,i) + this%V(i,j,k)
          vel(3,i) = vel(3,i) + this%W(i,j,k)
          press(i) = press(i) + this%p(i,j,k)
          mp = this%rhof(i,j,k) * this%ps%cfg%dx(i) * this%ps%cfg%dy(j)       &
            * this%ps%cfg%dz(k)
          fke(i) = fke(i) + 0.5_WP * mp * sum(((/ this%U(i,j,k),              &
            this%V(i,j,k), this%W(i,j,k) /) - velmean)**2)
        end do
      end do
    end do
    call mpi_reduce(vel,   this%mstats%fvel, 3*this%ps%cfg%nx, MPI_REAL_WP,   &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(press, this%mstats%p,    this%ps%cfg%nx,   MPI_REAL_WP,   &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(fke, this%mstats%fke,    this%ps%cfg%nx,   MPI_REAL_WP,   &
      MPI_SUM, 0, this%ps%cfg%comm)
    if (this%ps%cfg%rank .eq. 0) then
      this%mstats%fvel = this%mstats%fvel / (this%ps%cfg%ny * this%ps%cfg%nz)
      this%mstats%p    = this%mstats%p    / (this%ps%cfg%ny * this%ps%cfg%nz)
    end if

    ! write to monitor file
    call this%mstats%write(this%ps%cfg)

  end subroutine xflwstats_xdep_stats

  ! observation plane
  ! will mostly just call single below, but needs to be its own function to
  ! iterate through planes and filters, along with calling slice io and doing
  ! any needed setup work
  subroutine xflwstats_opln_stats(this, step, t)
    use mpi_f08, only: MPI_UNDEFINED
    implicit none
    class(xflwstats), intent(inout) :: this
    integer, intent(in) :: step
    real(WP), intent(in) :: t
    integer :: m, n, i

    do m = 1, this%num_planes

      ! update obsplane coupler
      !this%planes(m)%obs_lpt%p(:)%flag = 1   using array mode
      !call this%planes(m)%obs_lpt%recycle()  using array mode
      call this%planes(m)%lptcpl%push()
      call this%planes(m)%lptcpl%transfer()

      ! move to next if not in destination group
      if (.not. this%planes(m)%in_grp) cycle
      if (this%planes(m)%lptcpl%drank .eq. MPI_UNDEFINED) cycle

      ! pull from obsplane coupler
      call this%planes(m)%lptcpl%pull()

      ! update io pmesh
      if (this%sliceio_setup(m)) then
        call this%planes(m)%io_pm%reset()
        call this%planes(m)%io_pm%set_size(this%planes(m)%lptcpl%pulledparticlecount)
        do i = 1, this%planes(m)%lptcpl%pulledparticlecount
          this%planes(m)%io_pm%pos(:,i) = this%planes(m)%lptcpl%pulledparticles(i)%pos
          this%planes(m)%io_pm%var(1,i) = this%planes(m)%lptcpl%pulledparticles(i)%id
          this%planes(m)%io_pm%var(2,i) = this%planes(m)%lptcpl%pulledparticles(i)%d
          this%planes(m)%io_pm%vec(:,1,i) = this%planes(m)%lptcpl%pulledparticles(i)%vel
          !TODO fluid vel?
          !ind = this%ps%cfg%get_ijk_global(this%io_lpt%p(i)%pos)
          !      this%io_pm%vec(:,2,i) = this%ps%cfg%get_velocity(             &
          !        pos=this%io_lpt%p(i)%pos, i0=ind(1), j0=ind(2), k0=ind(3),     &
          !        U=this%U, V=this%V, W=this%W)
        end do
      end if

      do n = 1, this%num_filters

        ! apply filter
        call this%planes(m)%apply_filter(this%filters(m,n),                   &
          this%planes(m)%lptcpl%pulledparticles, this%U, this%V, this%W)

        ! compute relevant stats and write monitor
        this%fstats(m,n)%step = step
        call compute_opln_stats_single(this%fstats(m,n),                      &
          this%planes(m)%fft, this%planes(m)%phimean, this%planes(m)%phi,     &
          this%planes(m)%phi_x, this%planes(m)%up, this%planes(m)%uf,         &
          this%planes(m)%work_r, this%planes(m)%work_c)
        if (this%planes(m)%fft%pg%rank .eq. 0)                                &
          call this%fstats(m,n)%mon%write()

        ! update Eulerian fields
        if (this%sliceio_setup(m)) then
          this%planes(m)%phi_io(:,:,:,n)   = this%planes(m)%phi(:,:,:)
          this%planes(m)%phi_x_io(:,:,:,n) = this%planes(m)%phi_x(:,:,:)
          this%planes(m)%up_io(:,:,:,:,n)  = this%planes(m)%up(:,:,:,:)
          this%planes(m)%uf_io(:,:,:,:,n)  = this%planes(m)%uf(:,:,:,:)
        end if

      end do

      ! write npy output for slice
      if (this%planes(m)%in_grp .and. this%sliceio_setup(m))                  &
        call this%planes(m)%io_f%write_data(t)

    end do

  end subroutine xflwstats_opln_stats

  ! computes statistics from filtered fields and puts them in the stats object
  subroutine compute_opln_stats_single(stats, fft, phimean, phi, phi_x,       &
      up, uf, work_r, work_c)
    use mpi_f08, only:  mpi_allreduce, MPI_SUM
    use parallel, only: MPI_REAL_WP
    implicit none
    type(xplanestats), intent(inout) :: stats
    type(cfourier), intent(inout) :: fft
    real(WP), intent(in) :: phimean
    real(WP), dimension(fft%pg%imino_:,fft%pg%jmino_:,fft%pg%kmino_:),        &
      intent(in) :: phi, phi_x
    real(WP), dimension(fft%pg%imino_:,fft%pg%jmino_:,fft%pg%kmino_:,1:),     &
      intent(in) :: up, uf
    real(WP), dimension(fft%pg%imino_:,fft%pg%jmino_:,fft%pg%kmino_:),        &
      intent(out) ::  work_r
    complex(WP), dimension(fft%pg%imino_:,fft%pg%jmino_:,fft%pg%kmino_:),     &
      intent(out) :: work_c
    real(WP) :: PC2_loc
    real(WP), dimension(3) :: UB_loc, SB_loc, PCUB_loc, UBPCG_loc
    real(WP), dimension(3,3) :: UC2_loc
    integer :: i, j, k, n, m, N2, ierr, lj, hj, lk, hk

    N2 = fft%pg%ny * fft%pg%nz
    i = fft%pg%imin_
    lj = fft%pg%jmin_; hj = fft%pg%jmax_;
    lk = fft%pg%kmin_; hk = fft%pg%kmax_;

    ! store PB
    stats%PB = phimean

    ! compute UB and UC2
    UC2_loc(:,:) = 0.0_WP
    do n = 1, 3; UB_loc(n) = sum(up(i,lj:hj,lk:hk,n)); end do
    call mpi_allreduce(UB_loc, stats%UB, 3, MPI_REAL_WP, MPI_SUM,             &
      fft%pg%comm, ierr)
    do n = 1, 3
      do m = n, 3
        work_r(i,:,:) = (up(i,:,:,n) - stats%UB(n)) * (up(i,:,:,m) - stats%UB(m))
        UC2_loc(n,m) = sum(work_r(i,lj:hj,lk:hk))
      end do
    end do
    call mpi_allreduce(UC2_loc, stats%UC2, 9, MPI_REAL_WP, MPI_SUM,           &
      fft%pg%comm, ierr)
    stats%UB = stats%UB * (phimean / N2)
    stats%UC2 = stats%UC2 * (phimean**2 / N2)

    ! compute SB
    do n = 1, 3
      work_r(i,:,:) = uf(i,:,:,n) - up(i,:,:,n)
      SB_loc(n) = sum(work_r(i,lj:hj,lk:hk))
    end do
    call mpi_allreduce(SB_loc, stats%SB, 3, MPI_REAL_WP, MPI_SUM,             &
      fft%pg%comm, ierr)
    stats%SB = stats%SB * (phimean / N2)

    ! compute PC2
    work_r(i,:,:) = phi(i,:,:)**2
    PC2_loc = sum(work_r(i,lj:hj,lk:hk))
    call mpi_allreduce(PC2_loc, stats%PC2, 1, MPI_REAL_WP, MPI_SUM,           &
      fft%pg%comm, ierr)
    stats%PC2 = stats%PC2 * (phimean**2 / N2)

    ! compute PCUB
    do n = 1, 3
      work_r(i,:,:) = phi(i,:,:) * up(i,:,:,n)
      PCUB_loc(n) = sum(work_r(i,lj:hj,lk:hk))
    end do
    call mpi_allreduce(PCUB_loc, stats%PCUB, 3, MPI_REAL_WP, MPI_SUM,         &
      fft%pg%comm, ierr)
    stats%PCUB = stats%PCUB * (phimean**2 / N2)

    ! store/compute products of upbar with gradients of phi check
    work_r(i,:,:) = up(i,:,:,1) * phi_x(i,:,:)
    UBPCG_loc(1) = sum(work_r(i,lj:hj,lk:hk))
    work_c(:,:,:) = phi(:,:,:)
    call fft%ytransform_forward(work_c)
    do k = fft%pg%kmin_, fft%pg%kmax_
      work_c(i,:,k) = im * fft%ky(:) * work_c(i,:,k)
    end do
    call fft%ytransform_backward(work_c)
    work_r(i,:,:) = up(i,:,:,2) * realpart(work_c(i,:,:))
    UBPCG_loc(2) = sum(work_r(i,lj:hj,lk:hk))
    work_c(:,:,:) = phi(:,:,:)
    call fft%ztransform_forward(work_c)
    do j = fft%pg%jmin_, fft%pg%jmax_
      work_c(i,j,:) = im * fft%kz(:) * work_c(i,j,:)
    end do
    call fft%ztransform_backward(work_c)
    work_r(i,:,:) = up(i,:,:,3) * realpart(work_c(i,:,:))
    UBPCG_loc(3) = sum(work_r(i,lj:hj,lk:hk))
    call mpi_allreduce(UBPCG_loc, stats%UBPCG, 3, MPI_REAL_WP, MPI_SUM,     &
      fft%pg%comm)
    stats%UBPCG = stats%UBPCG * (phimean**2 / N2)

  end subroutine compute_opln_stats_single

  subroutine xflwstats_destroy(this)
    use messager, only: die
    implicit none
    class(xflwstats), intent(inout) :: this

    deallocate(this%planes, this%filters, this%fstats, this%sliceio_setup,    &
      this%zeros)

  end subroutine xflwstats_destroy


end module xflwstats_class

