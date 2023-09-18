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
  type :: xdepstats
    type(monitor) :: globalmon
    integer :: step
    integer, dimension(18) :: units
    ! data only allocated / updated on root
    real(WP) :: rhof, muf, dp, ntot
    ! these arrays are only imin to imax and don't include overlap
    real(WP), dimension(:), allocatable :: num, vf, vfsq
    real(WP), dimension(:,:), allocatable :: vel, velsq, slp, slpsq, drg   ! 3 by nx
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
    procedure :: write => xplanestats_write
    final :: xplanestats_destroy
  end type xplanestats


  !> filter information
  type :: filter
    character(len=str_medium) :: filtername
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

    character(len=str_medium) :: label
    real(WP) :: offset

    type(MPI_GROUP) :: grp
    logical :: in_grp

    type(MPI_GROUP) :: fft_grp
    class(config), pointer :: sim_cfg
    type(pgrid) :: fft_pg
    type(cfourier), pointer :: fft
    logical :: in_fft_grp

    ! we use a coupler to subsample and relocate particles to the appropriate
    ! processes within some cutoff; as a result we end up needing a second grid
    ! for each plane; this second grid uses a processor group that is
    ! consistent with that of fft_pg, but is larger
    type(MPI_GROUP) :: obs_lpt_grp
    type(pgrid) :: obs_lpt_pg
    type(lptcoupler) :: obs_lptcpl
    logical :: in_lpt_grp

    ! these can't be shared across obsplanes because the z indices have to
    ! match, but they can be shared across filters
    complex(WP), dimension(:,:,:), allocatable :: work_c
    real(WP), dimension(:,:,:), allocatable :: work_r
    ! even though the memory for these is the same for each obsplane, we have
    ! to re-project for each filter
    real(WP) :: phimean
    real(WP), dimension(:,:,:), allocatable :: phi, phi_x
    real(WP), dimension(:,:,:,:), allocatable :: up, uf

    ! io
    logical :: io_setup = .false.
    logical :: write_slices = .false.
    type(npy) :: io_out
    type(partmesh) :: io_pmesh

  contains

    procedure :: init => obsplane_init
    procedure :: apply_filter => obsplane_apply_filter
    procedure :: destroy => obsplane_destroy
    !procedure :: io_setup => obsplane_io_setup
    !procedure :: io_write => obsplane_io_write
    !procedure :: io_destroy => obsplane_io_destroy

  end type obsplane


  !> top level xflow statistics object (interfaced with in simulation.f90)
  type :: xflwstats

    class(config), pointer :: sim_cfg

    real(WP), dimension(:,:,:), pointer :: vf, rhof, visc, U, V, W
    type(lpt), pointer :: ps

    integer :: num_filters, num_planes
    type(xdepstats) :: mstats
    type(obsplane), dimension(:), pointer :: planes
    type(filter), dimension(:), pointer :: filters
    type(xplanestats), dimension(:,:), pointer :: fstats  ! num_planes by num_filters

  contains

    procedure :: init => xflwstats_init
    !procedure :: setup_sliceio => xflwstats_setup_sliceio
    !procedure :: write_sliceio => xflwstats_write_sliceio
    procedure :: compute_xdep_stats => xflwstats_xdep_stats ! x-dependent
    procedure :: compute_opln_stats => xflwstats_opln_stats ! observation plane
    procedure :: destroy => xflwstats_destroy

  end type xflwstats


contains


  !> xdepstats (simple stats across xflow)

  subroutine xdepstats_init(this, mkdir, g)
    implicit none
    class(xdepstats), intent(inout) :: this
    logical, intent(in) :: mkdir
    class(pgrid), intent(in) :: g
    integer :: i

    allocate(this%num(g%imin:g%imax), this%vf(g%imin:g%imax),                 &
      this%vfsq(g%imin:g%imax), this%vel(3,g%imin:g%imax),                    &
      this%velsq(3,g%imin:g%imax), this%slp(3,g%imin:g%imax),                 &
      this%slpsq(3,g%imin:g%imax), this%drg(3,g%imin:g%imax))

    this%globalmon = monitor(g%rank .eq. 0, 'xflw_stat_global')
    call this%globalmon%add_column(this%step, 'step')
    call this%globalmon%add_column(this%rhof, 'rhof')
    call this%globalmon%add_column(this%muf,  'muf')
    call this%globalmon%add_column(this%dp,   'dp')
    call this%globalmon%add_column(this%ntot, 'ntot')

    ! just use raw files for all the x dependent statistics
    if (g%rank .eq. 0) then
      if (mkdir) call execute_command_line('mkdir xflw_stat')
      open(newunit=this%units( 1), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/num')
      open(newunit=this%units( 2), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/vf')
      open(newunit=this%units( 3), action='WRITE', form='FORMATTED',          &
        file='xflw_stat/vfsq')
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
      do i = 1, 18
        write(i,'(A8,*(e16.6))') ' x ', g%x(g%imin:g%imax)
      end do
    end if

  end subroutine xdepstats_init

  subroutine xdepstats_write(this, g)
    implicit none
    class(xdepstats), intent(inout) :: this
    class(pgrid), intent(in) :: g

    if (g%rank .eq. 0) then

      call this%globalmon%write()

      write( 1,'(I8,*(e16.6))') this%step, this%num(     g%imin:g%imax)
      write( 2,'(I8,*(e16.6))') this%step, this%vf(      g%imin:g%imax)
      write( 3,'(I8,*(e16.6))') this%step, this%vfsq(    g%imin:g%imax)
      write( 4,'(I8,*(e16.6))') this%step, this%vel(1,   g%imin:g%imax)
      write( 5,'(I8,*(e16.6))') this%step, this%vel(2,   g%imin:g%imax)
      write( 6,'(I8,*(e16.6))') this%step, this%vel(3,   g%imin:g%imax)
      write( 7,'(I8,*(e16.6))') this%step, this%velsq(1, g%imin:g%imax)
      write( 8,'(I8,*(e16.6))') this%step, this%velsq(2, g%imin:g%imax)
      write( 9,'(I8,*(e16.6))') this%step, this%velsq(3, g%imin:g%imax)
      write(10,'(I8,*(e16.6))') this%step, this%slp(1,   g%imin:g%imax)
      write(11,'(I8,*(e16.6))') this%step, this%slp(2,   g%imin:g%imax)
      write(12,'(I8,*(e16.6))') this%step, this%slp(3,   g%imin:g%imax)
      write(13,'(I8,*(e16.6))') this%step, this%slpsq(1, g%imin:g%imax)
      write(14,'(I8,*(e16.6))') this%step, this%slpsq(2, g%imin:g%imax)
      write(15,'(I8,*(e16.6))') this%step, this%slpsq(3, g%imin:g%imax)
      write(16,'(I8,*(e16.6))') this%step, this%drg(1,   g%imin:g%imax)
      write(17,'(I8,*(e16.6))') this%step, this%drg(2,   g%imin:g%imax)
      write(18,'(I8,*(e16.6))') this%step, this%drg(3,   g%imin:g%imax)

    end if

  end subroutine xdepstats_write

  subroutine xdepstats_destroy(this)
    implicit none
    type(xdepstats), intent(inout) :: this
    integer :: i

    ! monitor falls out of scope and has its own destructor

    do i = 1, 18
      close(unit=this%units(i))
    end do

    deallocate(this%num, this%vf, this%vfsq, this%vel, this%velsq, this%slp,  &
      this%slpsq, this%drg)

  end subroutine xdepstats_destroy


  !> macroscopic statistics object

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

  subroutine xplanestats_write(this)
    implicit none
    class(xplanestats), intent(inout) :: this

    call this%mon%write()

  end subroutine xplanestats_write

  subroutine xplanestats_destroy(this)
    implicit none
    type(xplanestats), intent(inout) :: this

    ! nothing to do (monitor falls out of scope and has its own destructor)

  end subroutine xplanestats_destroy


  !> filter object

  ! assumes filtername, params, and type are already set
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

  subroutine build_slice_pg(sim_pg, N, offset, width, obs_pg, obs_grp, in_grp)
    use mpi_f08, only: mpi_allreduce, MPI_LOGICAL, MPI_LOR, mpi_group_range_incl
    use sgrid_class, only: cartesian
    use messager, only: die
    implicit none
    class(pgrid), intent(in) :: sim_pg
    integer, dimension(3) :: N
    real(WP), intent(in) :: offset, width
    type(pgrid), intent(out) :: obs_pg
    type(MPI_GROUP), intent(out) :: obs_grp
    logical, intent(out) :: in_grp
    real(WP), dimension(3) :: dx
    integer, dimension(3,sim_pg%nproc) :: ranges
    logical, dimension(sim_pg%nproc) :: in_group_loc, in_group
    type(sgrid) :: obs_sg
    integer :: i, j, k, ierr
    logical :: prev
    real(WP), dimension(:), allocatable :: xs, ys, zs

    ! check plane is one cell wide and allocate xs, ys, zs arrays
    if (N(1) .ne. 1) call die("Observation plane must be one cell wide.")
    allocate(xs(N(1)+1), ys(N(2)+1), zs(N(3)+1))

    ! build sgrid and find overlapping processes
    dx(:) = (/ width, sim_pg%yL, sim_pg%zL /) / N(:)
    do i = 1, N(1) + 1; xs(i) = (real(i,WP) - 1.0_WP) * dx(1); end do;
    do i = 1, N(2) + 1; ys(i) = (real(i,WP) - 1.0_WP) * dx(2); end do;
    do i = 1, N(3) + 1; zs(i) = (real(i,WP) - 1.0_WP) * dx(3); end do;
    xs(:) = xs(:) + offset - 0.5_WP * width
    ys(:) = ys(:) + sim_pg%y(sim_pg%jmin)
    zs(:) = zs(:) + sim_pg%z(sim_pg%kmin)
    obs_sg = sgrid(coord=cartesian, no=0, x=xs, y=ys, z=zs,                   &
      xper=.true., yper=.true., zper=.true., name='XF_OBS_SG')
    in_group_loc(:) = .false.
    in_group_loc(sim_pg%rank+1) =                                             &
      offset + 0.5_WP * width .gt. sim_pg%x(sim_pg%imin_) .and.               &
      offset - 0.5_WP * width .lt. sim_pg%x(sim_pg%imax_+1)
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

    ! set up pgrid
    obs_pg = pgrid(obs_sg, obs_grp, N)

    ! free memory
    deallocate(xs, ys, zs)

  end subroutine build_slice_pg

  subroutine obsplane_init(this, sim_cfg, FFTN, ps, offset, fft_width, lpt_width)
    use messager, only: die
    use mpi_f08, only: mpi_allreduce, MPI_LOGICAL, MPI_LOR, mpi_group_range_incl
    use sgrid_class, only: cartesian
    implicit none
    class(obsplane), intent(inout) :: this
    class(config), target, intent(in) :: sim_cfg
    integer, dimension(3), intent(in) :: FFTN
    class(lpt), target, intent(in) :: ps            ! lpt with xflow particles
    real(WP), intent(in) :: offset, fft_width, lpt_width
    integer :: li, hi, lj, hj, lk, hk

    ! set basic properties
    this%sim_cfg => sim_cfg
    this%offset = offset
    write(this%label, '(A,F10.6)') 'obsplane_', offset

    ! set up lpt slice
    call build_slice_pg(sim_cfg, FFTN, offset, lpt_width, this%obs_lpt_pg,    &
      this%obs_lpt_grp, this%in_lpt_grp)

    ! set up fft slice
    call build_slice_pg(sim_cfg, FFTN, offset, fft_width, this%fft_pg,        &
      this%fft_grp, this%in_fft_grp)

    ! set up fft
    allocate(this%fft, source=cfourier(this%fft_pg))

    ! setup particle holder for projection onto plane and particle coupler
    this%obs_lptcpl = lptcoupler(src_grp=sim_cfg%group,                       &
      dst_grp=this%obs_lpt_grp, name='OBS_LPT_CPL')
    call this%obs_lptcpl%set_src(ps)
    if (this%in_lpt_grp) call this%obs_lptcpl%set_dst(this%obs_lpt_pg)
    call this%obs_lptcpl%initialize()

    ! allocate arrays
    li = this%fft_pg%imin_; lj = this%fft_pg%jmin_; lk = this%fft_pg%kmin_;
    hi = this%fft_pg%imax_; hj = this%fft_pg%jmax_; hk = this%fft_pg%kmax_;
    allocate(this%work_c(li:hi,lj:hj,lk:hk), this%work_r(li:hi,lj:hj,lk:hk),  &
      this%phi(li:hi,lj:hj,lk:hk), this%phi_x(li:hi,lj:hj,lk:hk),             &
      this%up(li:hi,lj:hj,lk:hk,3), this%uf(li:hi,lj:hj,lk:hk,3))

  end subroutine obsplane_init

  ! assumes obs_lpt is up to date
  ! U, V, W, Ux, Vx, Wx are from flow solver and on flow solver pgrid
  subroutine obsplane_apply_filter(this, flt, ps, U, V, W, Ux, Vx, Wx)
    use mpi_f08, only: mpi_allreduce, MPI_SUM
    use parallel, only: MPI_REAL_WP
    use lpt_class, only: part
    use mathtools, only: pi
    implicit none
    class(obsplane), intent(inout) :: this
    type(filter), intent(in) :: flt
    type(part), dimension(:), intent(inout) :: ps
    real(WP), dimension(this%sim_cfg%imin_:,this%sim_cfg%jmin_:,              &
      this%sim_cfg%kmin_:), intent(in) :: U, V, W, Ux, Vx, Wx
    integer :: n, i, j, k
    integer, dimension(3) :: ind, dom_n, dom_imin
    real(WP), dimension(3) :: fvel, dom_c, dom_d
    real(WP) :: b, cellvoli, dx, pvol

    if (.not. this%in_grp) return

    ! project particles 
    this%phi(:,:,:) = 0.0_WP; this%phi_x(:,:,:) = 0.0_WP;
    this%up(:,:,:,:) = 0.0_WP; this%uf(:,:,:,:) = 0.0_WP;
    dom_c = (/ this%fft_pg%x(this%fft_pg%imin),                               &
      this%fft_pg%y(this%fft_pg%jmin), this%fft_pg%z(this%fft_pg%kmin) /)
    dom_n = (/ this%fft_pg%nx, this%fft_pg%ny, this%fft_pg%nz /)
    dom_d = (/ this%fft_pg%xL, this%fft_pg%yL, this%fft_pg%zL /) / dom_n
    dom_imin = (/ this%fft_pg%imin, this%fft_pg%jmin, this%fft_pg%kmin /)
    do n = 1, size(ps)
      pvol = pi * ps(n)%d**3 / 6.0_WP
      ps(n)%ind = this%sim_cfg%get_ijk_global(ps(n)%pos, ps(n)%ind)
      fvel = this%sim_cfg%get_velocity(pos=ps(n)%pos, i0=ps(n)%ind(1),         &
        j0=ps(n)%ind(2), k0=ps(n)%ind(3), U=U, V=V, W=W)
      ind = int((ps(n)%pos - dom_c) / dom_d) + dom_imin
      ind(:) = this%fft_pg%get_ijk_global(ps(n)%pos, ind)
      ind(1) = this%fft_pg%imin_
      dx = ps(n)%pos(1) - this%offset
      b = -6.0_WP / flt%params(1)**2
      this%phi(i,j,k) = this%phi(i,j,k) + pvol * exp(b * dx**2)
      this%phi_x(i,j,k) = this%phi_x(i,j,k) +                                 &
        2 * b * dx * pvol * exp(b * dx**2)
      this%up(i,j,k,:) = this%up(i,j,k,:) +                                   &
        ps(n)%vel * pvol * exp(b * dx**2)
      this%uf(i,j,k,:) = this%uf(i,j,k,:) + fvel * pvol * exp(b * dx**2)
    end do

    ! normalize phi
    cellvoli = this%fft_pg%nx * this%fft_pg%ny * this%fft_pg%nz /             &
      this%fft_pg%vol_total
    this%phi(:,:,:) = cellvoli * this%phi(:,:,:)
    b = sum(this%phi(this%fft_pg%imin_:this%fft_pg%imax_,                     &
                     this%fft_pg%jmin_:this%fft_pg%jmax_,                     &
                     this%fft_pg%kmin_:this%fft_pg%kmax_))
    call mpi_allreduce(b, this%phimean, 1, MPI_REAL_WP, MPI_SUM,              &
      this%fft_pg%comm)
    this%phimean = this%phimean / (this%fft_pg%nx * this%fft_pg%ny *          &
      this%fft_pg%nz)
    this%phi(:,:,:) = this%phi(:,:,:) - this%phimean

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

    !TODO enable this once written
    !if (this%io_setup) call this%io_destroy()

  end subroutine obsplane_destroy

  !TODO this method hasn't been touched since being pulled from the cube
  !TODO this is actually not necessary/much easier, since we already have a single obsplane, we just need to write them out if the
  !flag is switched
  !TODO however, note that this is filter dependent, and must be set up for each filter
  !subroutine obsplane_io_setup(this)
  !  use mpi_f08    !,  only: mpi_group_range_incl, MPI_LOGICAL, MPI_LOR,           &
  !    !mpi_allreduce, mpi_bcast
  !  use parallel, only: MPI_REAL_WP
  !  use param,    only: param_read
  !  use messager, only: die
  !  use sgrid_class, only: cartesian
  !  implicit none
  !  class(obsplane), intent(inout) :: this
  !  type(filter), pointer :: flt
  !  integer :: n, li, hi, lj, hj, lk, hk, i, j, ierr
  !  real(WP), dimension(:,:,:), pointer :: vx, vy, vz
  !  real(WP) :: thickness, height
  !  real(WP), dimension(2) :: z
  !  integer, dimension(3) :: partition
  !  logical, dimension(this%sim_pg%nproc) :: in_group_loc, in_group
  !  integer, dimension(3,this%sim_pg%nproc) :: ranges
  !  logical :: prev

  !  if (.not. this%in_grp) return

  !  this%io_pmesh = partmesh(nvar=2, nvec=2, name='lpt')
  !  this%io_pmesh%varname(1) = "id"
  !  this%io_pmesh%varname(2) = "dp"
  !  this%io_pmesh%vecname(1) = "vel"
  !  this%io_pmesh%vecname(2) = "fld"
  !  this%io_out = npy(pg=this%obs_lpt_pg, folder='obsvis'//this%label)
  !  call this%io_out%add_particle('partslice', this%io_pmesh)

  !  ! set up output scalars for each interesting filter quantity
  !  !TODO
  !  do n = 1, this%num_filters
  !    flt => this%filters(n)
  !    if (.not. flt%use_slice_io .or. .not. this%in_io_grp) cycle
  !    li = this%io_cfg%imino_; hi = this%io_cfg%imaxo_;
  !    lj = this%io_cfg%jmino_; hj = this%io_cfg%jmaxo_;
  !    lk = this%io_cfg%kmino_; hk = this%io_cfg%kmaxo_;
  !    allocate(flt%io_phi(li:hi,lj:hj,lk:hk),                                &
  !      flt%io_up(li:hi,lj:hj,lk:hk,3), flt%io_uf(li:hi,lj:hj,lk:hk,3))
  !    call this%io_out%add_scalar(trim(flt%filtername) // 'phi', flt%io_phi)
  !    vx => flt%io_up(:,:,:,1); vy => flt%io_up(:,:,:,2); vz => flt%io_up(:,:,:,3);
  !    call this%io_out%add_vector(trim(flt%filtername) // 'up', vx, vy, vz)
  !    vx => flt%io_uf(:,:,:,1); vy => flt%io_uf(:,:,:,2); vz => flt%io_uf(:,:,:,3);
  !    call this%io_out%add_vector(trim(flt%filtername) // 'uf', vx, vy, vz)
  !  end do

  !end subroutine obsplane_io_setup

  !!TODO this method hasn't been touched since being pulled from the cube
  !subroutine obsplane_io_write(this, t)
  !  implicit none
  !  class(xflwstats), intent(inout) :: this
  !  real(WP), intent(in) :: t

  !  if (this%in_io_grp) call this%io_out%write_data(t)

  !end subroutine obsplane_io_write

  !TODO obsplane_io_destroy


  !> xflwstats

  subroutine xflwstats_init(this, sim_cfg, filterfile, FFTN, vf, rhof, visc, U, V, W, ps)
    use param,       only: param_read
    use mpi_f08,     only: mpi_bcast, MPI_REAL, MPI_INTEGER, MPI_CHARACTER,   &
      MPI_LOGICAL, MPI_COMM_WORLD
    use parallel,    only: MPI_REAL_WP
    use sgrid_class, only: cartesian
    use messager,    only: die
    implicit none
    class(xflwstats), intent(inout) :: this
    class(config), target, intent(in) :: sim_cfg
    character(len=str_medium), intent(in) :: filterfile
    integer, dimension(2), intent(in) :: FFTN
    type(lpt), target, intent(in) :: ps
    real(WP), dimension(sim_cfg%imino_:,sim_cfg%jmino_:,sim_cfg%kmino_:),     &
      intent(in), target :: vf, rhof, visc, U, V, W
    real(WP), dimension(:), allocatable :: planelocs, lptwidths
    type(filter_info_row) :: f_info_raw
    integer :: i, j, fh, ierr
    real(WP) :: fftwidth
    logical :: use_slice_io

    ! store pointer to simulation config and data
    this%sim_cfg => sim_cfg; this%ps => ps;
    this%vf => vf; this%rhof => rhof; this%visc => visc;
    this%U => U; this%V => V; this%W => W;

    ! set up microscopic statistics
    call this%mstats%init(sim_cfg%amroot, sim_cfg)

    ! read observation plane info
    if (sim_cfg%rank .eq. 0) then
      call param_read('HS obs plane count', this%num_planes)
      allocate(planelocs(this%num_planes),lptwidths(this%num_planes))
      call param_read('HS obs plane locations', planelocs)
      call param_read('HS obs plane lpt widths', lptwidths)
      call param_read('HS obs plane fft width', fftwidth)
    end if

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
      allocate(this%filters(this%num_filters))
      rewind(fh)
      do i = 1, this%num_filters
        read(fh,*,iostat=ierr) f_info_raw
        if (ierr .ne. 0) call die("[EC] error reading filterfile row")
        select case (f_info_raw%typename)
        case ('gaussian', 'GAUSSIAN', 'g', 'G')
          this%filters(i)%type = FLT_GAUSSIAN
        case ('box', 'BOX', 'b', 'B')
          this%filters(i)%type = FLT_BOX
        case default
          call die("[EC] unknown filter type '"//f_info_raw%typename//"'")
        end select
        this%filters(i)%filtername = f_info_raw%out_fname
        this%filters(i)%params(:) = f_info_raw%params(:)
        !TODO adapt sliceio to new object layout
        !this%filters(i)%use_slice_io = use_slice_io
      end do
      close(fh)
    end if

    ! broadcast filter info
    call mpi_bcast(this%num_filters, 1, MPI_INTEGER, 0, sim_cfg%comm, ierr)
    if (sim_cfg%rank .ne. 0) allocate(this%filters(this%num_filters))
    do i = 1, this%num_filters
      call mpi_bcast(this%filters(i)%filtername, str_medium, MPI_CHARACTER, 0,&
        sim_cfg%comm, ierr)
      this%filters(i)%filtername = trim(this%filters(i)%filtername)
      call mpi_bcast(this%filters(i)%type, 1, MPI_INTEGER, 0, sim_cfg%comm,   &
        ierr)
      call mpi_bcast(this%filters(i)%params, FLT_NUM_PARAMS, MPI_REAL_WP, 0,  &
        sim_cfg%comm, ierr)
      !call mpi_bcast(this%filters(i)%use_slice_io, 1, MPI_LOGICAL, 0,         &
      !  sim_cfg%comm, ierr)
    end do

    ! init observation planes
    allocate(this%planes(this%num_planes))
    do i = 1, this%num_planes
      call this%planes(i)%init(sim_cfg, (/ 1, FFTN(1), FFTN(2) /), ps,        &
        planelocs(i), fftwidth, lptwidths(i))
    end do

    ! init filters
    if (this%num_planes .eq. 0 .and. this%num_filters .gt. 0) call die(       &
      "Must have at least one observation plane if any filters are used.")
    do i = 1, this%num_filters
      call this%filters(i)%init(this%planes(1)%fft)
    end do

    ! initialize fstats objects
    allocate(this%fstats(this%num_planes,this%num_filters))
    do j = 1, this%num_filters
      do i = 1, this%num_planes
        !TODO check this is compatible with the root chosen at call / write time
        call this%fstats(i,j)%init(this%planes(1)%fft%pg%rank .eq. 0)
      end do
    end do

    ! set all processes to not be io processes until it's set up
    !TODO adapt this to new layout
    !this%in_io_grp = .false.

    ! cleanup
    deallocate(planelocs)

  end subroutine xflwstats_init

  !TODO xflwstats_setup_sliceio

  !TODO xflwstats_write_sliceio

  ! computes microscopic statistics along xflow region
  subroutine xflwstats_xdep_stats(this, step)
    use mpi_f08, only:   mpi_reduce, MPI_SUM, MPI_INTEGER
    use parallel, only:  MPI_REAL_WP
    implicit none
    class(xflwstats), intent(inout) :: this
    integer, intent(in) :: step
    integer :: n, i, j, k
    !real(WP) :: opt_dt TODO add this back in
    integer, dimension(3) :: indmin, indmax
    real(WP), dimension(3) :: acc, fvel
    integer, dimension(this%ps%cfg%imin:this%ps%cfg%imax) :: num
    real(WP), dimension(this%ps%cfg%imin:this%ps%cfg%imax) :: vf, vfsq
    real(WP), dimension(3,this%ps%cfg%imin:this%ps%cfg%imax) :: vel, velsq,   &
      slp, slpsq, drg

    ! to be computed:
    !real(WP) :: rhof, muf, dp, ntot
    !! these arrays are only imin to imax and don't include overlap
    !real(WP), dimension(:), allocatable :: num, vf, vfsq

    indmin(:) = (/ this%ps%cfg%imin, this%ps%cfg%jmin, this%ps%cfg%kmin /)
    indmax(:) = (/ this%ps%cfg%imax, this%ps%cfg%jmax, this%ps%cfg%kmax /)

    ! zero out arrays
    num(:) = 0.0_WP; vf(:) = 0.0_WP; vfsq(:) = 0.0_WP;
    vel(:,:) = 0.0_WP; velsq(:,:) = 0.0_WP;
    slp(:,:) = 0.0_WP; slpsq(:,:) = 0.0_WP;
    drg(:,:) = 0.0_WP;

    ! compute microscopic statistics
    do n = 1, this%ps%np_
      this%ps%p(n)%ind = this%ps%cfg%get_ijk_global(this%ps%p(n)%pos,         &
        this%ps%p(n)%ind)
      i = this%ps%p(n)%ind(1); j = this%ps%p(n)%ind(2); k = this%ps%p(n)%ind(3);
      if (i .lt. this%ps%cfg%imin_) cycle
      if (i .gt. this%ps%cfg%imax_) cycle
      num(i) = num(i) + 1
      vel(:,i) = vel(:,i) + this%ps%p(n)%vel
      velsq(:,i) = velsq(:,i) + this%ps%p(n)%vel**2
      fvel = this%ps%cfg%get_velocity(pos=this%ps%p(n)%pos, i0=i, j0=j, k0=k, &
        U=this%U, V=this%V, W=this%W)
      slp(:,i) = slp(:,i) + fvel - this%ps%p(n)%vel
      slpsq(:,i) = slpsq(:,i) + (fvel - this%ps%p(n)%vel)**2
      !TODO check this is the right viscosity
      !call this%ps%get_rhs(U=this%U, V=this%V, W=this%W, rho=this%rhof,       &
      !  visc=this%visc, p=this%ps%p(n), acc=acc, opt_dt=opt_dt)
      !TODO repopulate fluid stress back through these routines
      acc(:) = 0.0_WP
      drg(:,i) = drg(:,i) + acc
    end do

    ! compute coarse statistics across simulation grid
    do k = this%ps%cfg%kmin_, this%ps%cfg%kmax_
      do j = this%ps%cfg%jmin_, this%ps%cfg%jmax_
        do i = this%ps%cfg%imin_, this%ps%cfg%imax_
          vf(i) = vf(i) + this%vf(i,j,k)
          vfsq(i) = vfsq(i) + this%vf(i,j,k)**2
        end do
      end do
    end do

    ! collect statistics to root
    call mpi_reduce(this%mstats%num,   num,   this%ps%cfg%nx,   MPI_INTEGER,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(this%mstats%vel,   vel,   3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(this%mstats%velsq, velsq, 3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(this%mstats%slp,   slp,   3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(this%mstats%slpsq, slpsq, 3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)
    call mpi_reduce(this%mstats%drg,   drg,   3*this%ps%cfg%nx, MPI_REAL_WP,  &
      MPI_SUM, 0, this%ps%cfg%comm)

    ! sum particles in xflow region
    if (this%ps%cfg%rank .eq. 0) this%mstats%ntot = sum(this%mstats%num)

    ! note that currently writes are done in a separate method

  end subroutine xflwstats_xdep_stats

  !TODO xflwstats_opln_stats ! observation plane
  ! will mostly just call single below, but needs to be its own function to iterate through planes and filters, along with calling
  ! slice io and doing any needed setup work
  !TODO this has not been edited since being copied from the old version; none of the micro stats/setup stuff is
  !correct, the iteration through filters has changed, couplers need to be added for particles, etc
  subroutine xflwstats_opln_stats(this, step, U, V, W, Ux, Vx, Wx)
    implicit none
    class(xflwstats), intent(inout) :: this
    integer, intent(in) :: step
    real(WP), dimension(this%sim_cfg%imin_:,this%sim_cfg%jmin_:,              &
      this%sim_cfg%kmin_:), intent(in) :: U, V, W, Ux, Vx, Wx
    integer :: m, n

    do m = 1, this%num_planes

      ! update obsplane coupler
      !this%planes(m)%obs_lpt%p(:)%flag = 1   using array mode
      !call this%planes(m)%obs_lpt%recycle()  using array mode
      call this%planes(m)%obs_lptcpl%push()
      call this%planes(m)%obs_lptcpl%transfer()
      if (this%planes(m)%in_grp) call this%planes(m)%obs_lptcpl%pull()

      do n = 1, this%num_filters

        ! apply filter
        call this%planes(m)%apply_filter(this%filters(n),                     &
          this%planes(m)%obs_lptcpl%pulledparticles, U, V, W, Ux, Vx, Wx)

        ! compute relevant stats and write monitor
        this%fstats(m,n)%step = step
        call compute_opln_stats_single(this%fstats(m,n),                      &
          this%planes(m)%fft, this%planes(m)%phimean, this%planes(m)%phi,     &
          this%planes(m)%phi_x, this%planes(m)%up, this%planes(m)%uf,         &
          this%planes(m)%work_r, this%planes(m)%work_c)
        call this%fstats(m,n)%mon%write()

        !TODO update IO (call obsplane routines or move them here, but only
        ! keep one of the two
        ! xy slice first
        !if (this%filters(n)%use_slice_io) then
        !  ! transfer continuum properties
        !  call this%io_cpl%push(this%phi)
        !  call this%io_cpl%transfer()
        !  if (this%in_io_grp) call this%io_cpl%pull(this%filters(n)%io_phi)
        !  do m = 1, 3
        !    call this%io_cpl%push(this%up(:,:,:,m))
        !    call this%io_cpl%transfer()
        !    if (this%in_io_grp) call this%io_cpl%pull(this%filters(n)%io_up(:,:,:,m))
        !    call this%io_cpl%push(this%uf(:,:,:,m))
        !    call this%io_cpl%transfer()
        !    if (this%in_io_grp) call this%io_cpl%pull(this%filters(n)%io_uf(:,:,:,m))
        !  end do
        !  ! transfer particle properties to lpt
        !  call this%io_lptcpl%push()
        !  call this%io_lptcpl%transfer()
        !  if (this%in_io_grp) call this%io_lptcpl%pull()
        !  ! transfer particles to pmesh
        !  if (this%in_io_grp) then
        !    call this%io_pmesh%reset()
        !    call this%io_pmesh%set_size(this%io_lpt%np_)
        !    do i = 1, this%io_lpt%np_
        !      this%io_pmesh%pos(:,i) = this%io_lpt%p(i)%pos
        !      this%io_pmesh%var(1,i) = this%io_lpt%p(i)%id
        !      this%io_pmesh%var(2,i) = this%io_lpt%p(i)%d
        !      this%io_pmesh%vec(:,1,i) = this%io_lpt%p(i)%vel
        !      ind = this%ps%cfg%get_ijk_global(this%io_lpt%p(i)%pos)
        !      this%io_pmesh%vec(:,2,i) = this%ps%cfg%get_velocity(             &
        !        pos=this%io_lpt%p(i)%pos, i0=ind(1), j0=ind(2), k0=ind(3),     &
        !        U=this%U, V=this%V, W=this%W)
        !    end do
        !  end if
        !end if

      end do

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
    real(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:),           &
      intent(in) :: phi, phi_x
    real(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:,1:),        &
      intent(in) :: up, uf
    real(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:),           &
      intent(out) ::  work_r
    complex(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:),        &
      intent(out) :: work_c
    real(WP) :: PC2_loc
    real(WP), dimension(3) :: UB_loc, SB_loc, PCUB_loc, UBPCG_loc
    real(WP), dimension(3,3) :: UC2_loc
    integer :: i, j, k, n, m, N2, ierr

    !TODO check fft_pg has no overlap, or fix indices

    N2 = fft%pg%ny * fft%pg%nz
    i = fft%pg%imin_

    ! store PB
    stats%PB = phimean

    ! compute UB
    do n = 1, 3; UB_loc(n) = sum(up(i,:,:,n)); end do
    call mpi_allreduce(UB_loc, stats%UB, 3, MPI_REAL_WP, MPI_SUM,             &
      fft%pg%comm, ierr)
    stats%UB = stats%UB / N2

    ! compute Reynolds-stress-type terms
    !TODO check this has the correct mean
    UC2_loc(n,m) = 0.0_WP
    do n = 1, 3
      do m = n, 3
        work_r(i,:,:) = up(i,:,:,n) * up(i,:,:,m)
        UC2_loc(n,m) = sum(work_r(i,:,:))
      end do
    end do
    call mpi_allreduce(UC2_loc, stats%UC2, 9, MPI_REAL_WP, MPI_SUM,           &
      fft%pg%comm, ierr)
    stats%UC2 = stats%UC2 / N2

    ! compute SB
    do n = 1, 3
      work_r(i,:,:) = uf(i,:,:,n) - up(i,:,:,n)
      SB_loc(n) = sum(work_r(i,:,:))
    end do
    call mpi_allreduce(SB_loc, stats%SB, 3, MPI_REAL_WP, MPI_SUM,             &
      fft%pg%comm, ierr)
    stats%SB = stats%SB / N2

    ! compute PC2
    work_r(i,:,:) = phi(i,:,:)**2
    PC2_loc = sum(work_r(i,:,:))
    call mpi_allreduce(PC2_loc, stats%PC2, 1, MPI_REAL_WP, MPI_SUM,           &
      fft%pg%comm, ierr)
    stats%PC2 = stats%PC2 / N2

    ! compute UB
    do n = 1, 3
      UB_loc(n) = sum(up(i,:,:,n))
    end do
    call mpi_allreduce(UB_loc, stats%UB, 3, MPI_REAL_WP, MPI_SUM,             &
      fft%pg%comm, ierr)
    stats%UB = stats%UB / N2

    ! compute PCUB
    do n = 1, 3
      work_r(i,:,:) = phi(i,:,:) * up(i,:,:,n)
      PCUB_loc(n) = sum(work_r(i,:,:))
    end do
    call mpi_allreduce(PCUB_loc, stats%PCUB, 3, MPI_REAL_WP, MPI_SUM,         &
      fft%pg%comm, ierr)
    stats%PCUB = stats%PCUB / N2

    ! store/compute products of upbar with gradients of phi check
    UBPCG_loc(1) = sum(up(i,:,:,1) * phi_x(i,:,:))
    work_c(:,:,:) = phi(:,:,:)
    call fft%ytransform_forward(work_c)
    do k = fft%pg%kmin_, fft%pg%kmax_
      do i = fft%pg%imin_, fft%pg%imax_
        work_c(i,:,k) = im * fft%ky(:) * work_c(i,:,k)
      end do
    end do
    call fft%ytransform_backward(work_c)
    UBPCG_loc(2) = sum(up(i,:,:,2) * realpart(work_c(i,:,:)))
    work_c(:,:,:) = phi(:,:,:)
    call fft%ztransform_forward(work_c)
    do j = fft%pg%jmin_, fft%pg%jmax_
      do i = fft%pg%imin_, fft%pg%imax_
        work_c(i,j,:) = im * fft%kz(:) * work_c(i,j,:)
      end do
    end do
    call fft%ztransform_backward(work_c)
    UBPCG_loc(3) = sum(up(i,:,:,3) * realpart(work_c(i,:,:)))
    call mpi_allreduce(UBPCG_loc, stats%UBPCG, 3, MPI_REAL_WP, MPI_SUM,     &
      fft%pg%comm)
    stats%UBPCG = stats%UBPCG / N2

  end subroutine compute_opln_stats_single

  subroutine xflwstats_destroy(this)
    use messager, only: die
    implicit none
    class(xflwstats), intent(inout) :: this
    type(filter), pointer :: flt
    integer :: n

    !TODO this

    do n = 1, this%num_filters
      flt => this%filters(n)
      !TODO this
    end do

  end subroutine xflwstats_destroy


end module xflwstats_class

