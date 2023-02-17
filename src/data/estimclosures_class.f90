!>
!> Module for estimation of closures in Shankar's two-fluid models.
!>
!> This module serves two purposes:
!>    1) estimation of statistics and io
!>    2) management of parameters and integral length scales
!>
!> This depends on p3dfft; even though it's a mess, I don't want to write my
!> own 3d fft parallelization.
!>
!> Originally written by John P Wakefield in Feb 2023
!>
module estimclosures_class
  use precision,      only: WP
  use mathtools,      only: pi
  use p3dfft_plus_plus
  use monitor_class,  only: monitor
  use string,         only: str_medium, str_long
  use pgrid_class,    only: pgrid
  use config_class,   only: config
  use lpt_class,      only: lpt
  implicit none
  private

  public :: estimclosures

  !> give all filters six parameters, not all of which may be used
  integer, parameter :: FLT_NUM_PARAMS = 6

  !> filter types
  integer, parameter :: FLT_GAUSSIAN = 1

  !> sampling modes
  integer, parameter :: SEQ_MESH = 1

  !> structure containing statistics to store
  type :: statpoint
    ! time
    real(WP) :: time
    integer :: step
    ! particle parameters
    real(WP) :: dp, volfrac
    ! Two-Fluid State
    real(WP) :: meanphibar, meanphicheck2
    real(WP), dimension(3) :: meanmesopvel
    ! Closure Values
    real(WP), dimension(4) :: meanVV        ! xx, xy, xz, yz
    real(WP), dimension(3) :: meanphicheckupmcheck, meanphicheck2upmcheck
  end type statpoint

  !> filter information and corresponding monitor file
  type :: filterstat
    character(len=str_medium) :: out_fname
    type(monitor) :: mon
    type(statpoint) :: d
    integer :: type
    real(WP), dimension(FLT_NUM_PARAMS) :: params
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: flt_f
    class(filterstat), pointer :: next
  contains
    procedure :: setup => filter_setup
    procedure :: destruct => filter_destruct
  end type filterstat
  type :: filter_info_row
    character(len=str_medium) :: out_fname
    character(len=str_medium) :: typename
    real(WP), dimension(FLT_NUM_PARAMS) :: params
  end type

  !> object managing computation of closures for each filter and navigation
  !> of paramter space
  type :: estimclosures
    ! configs for simulation and fft
    class(config), pointer :: sim_cfg
    type(pgrid) :: fft_pg
    ! fft plans and arrays
    integer :: trans_f, trans_b
    logical :: has_zero_wavenumber
    integer(C_INT), dimension(3) :: gdims_f, ldims_f
    real(C_DOUBLE), dimension(:,:,:), allocatable :: phi_r, work_r_x,         &
      work_r_y, work_r_z
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: phi_f,        &
      vel_f_x, vel_f_y, vel_f_z, work_f_x, work_f_y, work_f_z
    ! list of filters
    integer :: num_filters
    class(filterstat), pointer :: first_filter
    contains
      procedure :: get_done
      procedure :: get_next_parameters
      procedure :: get_next_time
      procedure :: compute_statistics
      procedure :: destruct => estimclosures_destruct
  end type estimclosures

  interface estimclosures; procedure estimclosures_from_args; end interface

contains

  subroutine fft_setup(pg, trans_f, trans_b, gdims_f, ldims_f, has_zero_wavenumber)
    !TODO pull changes to this from Olivier's version
    !use mpi_f08,  only: MPI_COMM_WORLD
    use messager, only: die
    implicit none
    include "mpif.h"
    class(pgrid), intent(in) :: pg
    integer, intent(out) :: trans_f, trans_b
    logical, intent(out) :: has_zero_wavenumber
    integer(C_INT), dimension(3), intent(out) :: gdims_f, ldims_f
    integer :: type_ccr, type_rcc, ndcp, ndims
    integer(C_INT) :: p3dfft_grid, grid_r, grid_f
    integer, dimension(3) :: type_ids_f, type_ids_b, glob_start_r, glob_start_f
    integer(C_INT), dimension(3) :: pdims, ldims_r, mem_order_r, mem_order_f, &
      dmap_r, dmap_f, gdims_r

    ! Set up global dimensions of the grid
    gdims_r(:) = (/ pg%nx, pg%ny, pg%nz /)

    ! Define the processor grid
    pdims(:) = (/ pg%npx, pg%npy, pg%npz /)

    ! Check periodicity and uniformity of mesh
    if (.not.(pg%xper.and.pg%uniform_x)) call die('[estimclosure] &
      &Need x-direction needs to be periodic and uniform')
    if (.not.(pg%yper.and.pg%uniform_y)) call die('[estimclosure] &
      &Need y-direction needs to be periodic and uniform')
    if (.not.(pg%zper.and.pg%uniform_z)) call die('[estimclosure] &
      &Need z-direction needs to be periodic and uniform')

    ! Ensure that we have at least one non-decomposed direction
    ndims = count(gdims_r .gt. 1); ndcp = count(pdims .gt. 1);
    if (ndcp .ge. ndims) call die('[pdfft3] Need at least one NON-decomposed direction')

    ! Set up work structures for P3DFFT
    call p3dfft_setup

    ! Set up 2 transform types for 3D transforms
    type_ids_f(:) = (/ P3DFFT_R2CFFT_D, P3DFFT_CFFT_FORWARD_D, P3DFFT_CFFT_FORWARD_D /)
    type_ids_b(:) = (/ P3DFFT_C2RFFT_D, P3DFFT_CFFT_BACKWARD_D, P3DFFT_CFFT_BACKWARD_D /)

    ! Now initialize 3D transforms (forward and backward) with these types
    call p3dfft_init_3Dtype(type_rcc, type_ids_f)
    call p3dfft_init_3Dtype(type_ccr, type_ids_b)

    ! set fourier dimensions, first index is halved due to real symmetry
    gdims_f(:) = gdims_r(:)
    gdims_f(1) = gdims_f(1) / 2 + 1

    ! set memory access order; cannot have a singleton dimension be decomposed,
    ! and mem_order / dmap permute processor coordinates; for 2d we just
    ! use 0, 1, 2 and ensure nx is not one
    mem_order_r(:) = (/ 0, 1, 2 /); dmap_r(:) = (/ 0, 1, 2 /);
    if (any(gdims_r .eq. 1)) then                      ! 2d case
      mem_order_f(:) = (/ 0, 1, 2 /); dmap_f(:) = (/ 0, 1, 2 /);
      if (gdims_r(1) .eq. 1) call die("[pf3dft] cannot run 2d with nx = 1")
    else                                                  ! 3d case
      mem_order_f(:) = (/ 1, 2, 0 /); dmap_f(:) = (/ 1, 2, 0 /);
    end if

    ! Specify the default communicator for P3DFFT++. This can be different
    ! from your program default communicator.
    p3dfft_grid = p3dfft_init_proc_grid(pdims, MPI_COMM_WORLD)

    ! Initialize initial grid, no conjugate symmetry (-1)
    call p3dfft_init_data_grid(grid_r, ldims_r, glob_start_r,                 &
      gdims_r, -1, p3dfft_grid, dmap_r, mem_order_r)

    ! Check ldims_real against pgrid
    if (.not. all((/ pg%nx_, pg%ny_, pg%nz_ /) .eq. ldims_r))                 &
      call die("[pfftd3] parallel fft decomposition does not match pg")

    ! Final grid has conjugate symmetry in X dimension (0)
    call p3dfft_init_data_grid(grid_f, ldims_f, glob_start_f,&
      gdims_f, 0, p3dfft_grid, dmap_f, mem_order_f)

    ! Plan transforms
    call p3dfft_plan_3Dtrans(trans_f, grid_r, grid_f, type_rcc)
    call p3dfft_plan_3Dtrans(trans_b, grid_f, grid_r, type_ccr)

    ! mark zero wavenumber
    has_zero_wavenumber = all((/ pg%iproc, pg%jproc, pg%kproc /)  &
      .eq. 1)

  end subroutine fft_setup

  ! assumes out_fname, type, num_params, params are set
  ! assumes fft has been setup
  subroutine filter_setup(flt, fft_pg, trans_f, ldims_f, work_r)
    use messager, only: die
    use mpi_f08, only: MPI_SUM, mpi_allreduce
    use parallel, only: MPI_REAL_WP
    implicit none
    class(filterstat), intent(inout) :: flt
    class(pgrid), intent(in) :: fft_pg
    integer :: trans_f
    integer, dimension(3), intent(in) :: ldims_f
    real(WP), dimension(fft_pg%imin_:,fft_pg%jmin_:,fft_pg%kmin_:), intent(inout) :: work_r
    real(WP) :: r2_y, r2_z, a, b, flt_int_l, flt_int_g
    integer :: j, k

    ! setup monitor file
    flt%mon = monitor(fft_pg%amRoot, flt%out_fname)
    call flt%mon%add_column(flt%d%time, 'time')
    call flt%mon%add_column(flt%d%step, 'step')
    call flt%mon%add_column(flt%d%dp, 'dp')
    call flt%mon%add_column(flt%d%volfrac, 'volfrac')
    call flt%mon%add_column(flt%d%meanphibar, 'phibar')
    call flt%mon%add_column(flt%d%meanphicheck2, 'phichk2')
    call flt%mon%add_column(flt%d%meanmesopvel(1), 'mesopvelx')
    call flt%mon%add_column(flt%d%meanmesopvel(2), 'mesopvely')
    call flt%mon%add_column(flt%d%meanmesopvel(3), 'mesopvelz')
    call flt%mon%add_column(flt%d%meanVV(1), 'VVxx')
    call flt%mon%add_column(flt%d%meanVV(2), 'VVxy')
    call flt%mon%add_column(flt%d%meanVV(3), 'VVxz')
    call flt%mon%add_column(flt%d%meanVV(4), 'VVyz')
    call flt%mon%add_column(flt%d%meanphicheckupmcheck(1), 'phichkupmchkx')
    call flt%mon%add_column(flt%d%meanphicheckupmcheck(2), 'phichkupmchky')
    call flt%mon%add_column(flt%d%meanphicheckupmcheck(3), 'phichkupmchkz')
    call flt%mon%add_column(flt%d%meanphicheck2upmcheck(1), 'phichkupmchkx')
    call flt%mon%add_column(flt%d%meanphicheck2upmcheck(2), 'phichkupmchky')
    call flt%mon%add_column(flt%d%meanphicheck2upmcheck(3), 'phichkupmchkz')

    ! allocate Fourier space filter
    allocate(flt%flt_f(ldims_f(1),ldims_f(2),ldims_f(3)))

    ! fill real filter; we build fft_pg such that the corner is at origin
    select case (flt%type)
    case (FLT_GAUSSIAN)
      a = (6.0_WP / (pi * flt%params(1)**2))**1.5_WP
      b = -6.0_WP / flt%params(1)**2
      do k = fft_pg%kmin_, fft_pg%kmax_
        r2_z = fft_pg%z(k)**2
        do j = fft_pg%jmin_, fft_pg%jmax_
          r2_y = fft_pg%y(j)**2
          work_r(:,j,k) = a * exp(b * (r2_z + r2_y +                          &
            fft_pg%x(fft_pg%imin_:fft_pg%imax_)**2))
        end do
      end do
    case default
      call die("[estimclosures] invalid filter type")
    end select

    ! check filter is normalized
    flt_int_l = sum(work_r)
    call mpi_allreduce(flt_int_l, flt_int_g, 1, MPI_REAL_WP, MPI_SUM, fft_pg%comm)
    work_r = work_r / flt_int_g

    ! transform filter
    call p3dfft_3Dtrans_double(trans_f, work_r, flt%flt_f, 0)
    if (all((/ fft_pg%iproc, fft_pg%jproc, fft_pg%kproc /) .eq. 1)) then
      flt%flt_f(1,1,1) = 1.0_C_DOUBLE
    end if

  end subroutine filter_setup

  function estimclosures_from_args(sim_cfg) result(ec)
    use param,       only: param_read
    use mpi_f08,     only: mpi_bcast, MPI_REAL, MPI_INTEGER, MPI_CHARACTER
    use parallel,    only: MPI_REAL_WP
    use sgrid_class, only: sgrid, cartesian
    use parallel,    only: group
    use messager,    only: die
    implicit none
    type(estimclosures) :: ec
    class(config), target, intent(in) :: sim_cfg
    real(WP), dimension(3) :: L, dx
    integer, dimension(3) :: FFTN
    real(WP), dimension(:), allocatable :: x, y, z
    type(sgrid) :: fft_sg
    character(len=str_medium) :: filterfile
    class(filterstat), pointer :: currfilter, prevfilter
    type(filter_info_row) :: f_info_raw
    integer :: i, fh, ierr

    ! store pointer to simulation config
    ec%sim_cfg => sim_cfg

    ! setup ffts
    L(:) = (/ sim_cfg%xL, sim_cfg%yL, sim_cfg%zL /)
    call param_read('EC FFT Mesh', FFTN)
    dx(:) = L / FFTN
    allocate(x(FFTN(1)+1), y(FFTN(2)+1), z(FFTN(3)+1))
    do i = 1, FFTN(1) + 1; x(i) = real(i-1,WP) * dx(1); end do;
    do i = 1, FFTN(2) + 1; y(i) = real(i-1,WP) * dx(2); end do;
    do i = 1, FFTN(3) + 1; z(i) = real(i-1,WP) * dx(3); end do;
    fft_sg = sgrid(coord=cartesian, no=0, x=x, y=y, z=z, xper=.true.,         &
      yper=.true., zper=.true., name='EC_FFT_CFG')
    ec%fft_pg = pgrid(grp=group, decomp=(/ sim_cfg%npx, sim_cfg%npy,          &
      sim_cfg%npz /), grid=fft_sg)
    call fft_setup(ec%fft_pg, ec%trans_f, ec%trans_b, ec%gdims_f, ec%ldims_f, &
      ec%has_zero_wavenumber)

    ! allocate arrays
    allocate(ec%phi_r(ec%fft_pg%imin_:ec%fft_pg%imax_,                        &
      ec%fft_pg%jmin_:ec%fft_pg%jmax_, ec%fft_pg%kmin_:ec%fft_pg%kmax_))
    allocate(ec%work_r_x(ec%fft_pg%imin_:ec%fft_pg%imax_,                     &
      ec%fft_pg%jmin_:ec%fft_pg%jmax_, ec%fft_pg%kmin_:ec%fft_pg%kmax_))
    allocate(ec%work_r_y(ec%fft_pg%imin_:ec%fft_pg%imax_,                     &
      ec%fft_pg%jmin_:ec%fft_pg%jmax_, ec%fft_pg%kmin_:ec%fft_pg%kmax_))
    allocate(ec%work_r_z(ec%fft_pg%imin_:ec%fft_pg%imax_,                     &
      ec%fft_pg%jmin_:ec%fft_pg%jmax_, ec%fft_pg%kmin_:ec%fft_pg%kmax_))
    allocate(ec%phi_f(ec%ldims_f(1),ec%ldims_f(2),ec%ldims_f(3)))
    allocate(ec%vel_f_x(ec%ldims_f(1),ec%ldims_f(2),ec%ldims_f(3)))
    allocate(ec%vel_f_y(ec%ldims_f(1),ec%ldims_f(2),ec%ldims_f(3)))
    allocate(ec%vel_f_z(ec%ldims_f(1),ec%ldims_f(2),ec%ldims_f(3)))
    allocate(ec%work_f_x(ec%ldims_f(1),ec%ldims_f(2),ec%ldims_f(3)))
    allocate(ec%work_f_y(ec%ldims_f(1),ec%ldims_f(2),ec%ldims_f(3)))
    allocate(ec%work_f_z(ec%ldims_f(1),ec%ldims_f(2),ec%ldims_f(3)))

    ! read filter info (as root)
    if (sim_cfg%rank .eq. 0) then
      call param_read('EC Filter List', filterfile)
      open(newunit=fh, file=filterfile, action='READ', access='SEQUENTIAL')
      ec%num_filters = 0
      nullify(currfilter)
      do
        read(fh,*,iostat=ierr) f_info_raw
        if (ierr .ne. 0) exit
        prevfilter => currfilter
        nullify(currfilter)
        allocate(currfilter)
        currfilter%next => prevfilter
        ec%num_filters = ec%num_filters + 1
        select case (f_info_raw%typename)
        case ('gaussian', 'GAUSSIAN', 'g', 'G')
          currfilter%type = FLT_GAUSSIAN
        case default
          call die("[EC] unknown filter type '"//f_info_raw%typename//"'")
        end select
        currfilter%out_fname = f_info_raw%out_fname
        currfilter%params(:) = f_info_raw%params
      end do
      ec%first_filter => currfilter
      close(fh)
      write(*,'("[EC] read in ",I2," filters")') ec%num_filters
    end if

    ! broadcast filter info
    call mpi_bcast(ec%num_filters, 1, MPI_INTEGER, 0, sim_cfg%comm, ierr)
    if (sim_cfg%rank .ne. 0) then
      nullify(currfilter)
      do i = 1, ec%num_filters
        prevfilter => currfilter
        nullify(currfilter)
        allocate(currfilter)
        currfilter%next => prevfilter
      end do
      ec%first_filter => currfilter
    end if
    currfilter => ec%first_filter
    do i = 1, ec%num_filters
      call mpi_bcast(currfilter%out_fname, str_medium, MPI_CHARACTER, 0,      &
        sim_cfg%comm, ierr)
      call mpi_bcast(currfilter%type, 1, MPI_INTEGER, 0, sim_cfg%comm, ierr)
      call mpi_bcast(currfilter%params, FLT_NUM_PARAMS, MPI_REAL_WP, 0,       &
        sim_cfg%comm, ierr)
      currfilter => currfilter%next
    end do

    ! setup filters
    currfilter => ec%first_filter
    do while (associated(currfilter))
      call currfilter%setup(ec%fft_pg, ec%trans_f, ec%ldims_f, ec%work_r_x)
      currfilter => currfilter%next
    end do

  end function estimclosures_from_args

  subroutine get_done(ec, done)
    implicit none
    class(estimclosures), intent(inout) :: ec
    logical, intent(out) :: done




  end subroutine get_done

  subroutine get_next_parameters(ec)
    implicit none
    class(estimclosures), intent(inout) :: ec




  end subroutine get_next_parameters

  subroutine get_next_time(ec)
    implicit none
    class(estimclosures), intent(inout) :: ec







  end subroutine get_next_time

  subroutine filter_destruct(flt)
    implicit none
    class(filterstat), intent(inout) :: flt

    deallocate(flt%flt_f)

    if (associated(flt%next)) call flt%next%destruct()

  end subroutine filter_destruct

  subroutine estimclosures_destruct(ec, dealloc_sim_cfg)
    implicit none
    class(estimclosures), intent(inout) :: ec
    logical, intent(in), optional :: dealloc_sim_cfg

    if (present(dealloc_sim_cfg) .and. dealloc_sim_cfg) deallocate(ec%sim_cfg)

    deallocate(ec%phi_r, ec%work_r_x, ec%work_r_y, ec%work_r_z,               &
      ec%phi_f, ec%vel_f_x, ec%vel_f_y, ec%vel_f_z, ec%work_f_x, ec%work_f_y, &
      ec%work_f_z)

    ! recursively deletes list
    call ec%first_filter%destruct()

  end subroutine estimclosures_destruct

  subroutine compute_statistics(ec, time, step, ps)
    use mpi_f08, only:  MPI_CHARACTER, MPI_INTEGER, mpi_reduce, MPI_SUM
    use parallel, only: MPI_REAL_WP
    use messager, only: die
    implicit none
    class(estimclosures), intent(inout) :: ec
    real(WP), intent(in) :: time
    integer, intent(in) :: step
    class(lpt), intent(in) :: ps
    class(filterstat), pointer :: flt
    real(WP) :: fft_rescale, dp, volfrac
    real(WP) :: pvol, dp_local, volfrac_local, phibar_local
    real(WP), dimension(4) :: meanVV_local, meanVV
    real(WP), dimension(3) :: tmpvec
    integer :: n, i, j, k

    ! rescaling values
    fft_rescale = 1.0_WP / (ec%fft_pg%nx * ec%fft_pg%ny * ec%fft_pg%nz)

    ! project particle volumes and velocities to mesh, compute unfiltered
    ! quantities
    dp_local = 0.0_WP; volfrac_local = 0.0_WP; meanVV_local(:) = 0.0_WP;
    ec%phi_r(:,:,:) = 0.0_WP
    ec%work_r_x(:,:,:) = 0.0_WP
    ec%work_r_y(:,:,:) = 0.0_WP
    ec%work_r_z(:,:,:) = 0.0_WP
    do n = 1, ps%np_
      dp_local = dp_local + ps%p(n)%d
      pvol = pi * ps%p(n)%d**3 / 6.0_WP
      volfrac_local = volfrac_local + pvol
      meanVV_local(1) = meanVV_local(1) + ps%p(n)%vel(1) * ps%p(n)%vel(1)
      meanVV_local(2) = meanVV_local(2) + ps%p(n)%vel(1) * ps%p(n)%vel(2)
      meanVV_local(3) = meanVV_local(3) + ps%p(n)%vel(1) * ps%p(n)%vel(3)
      meanVV_local(4) = meanVV_local(4) + ps%p(n)%vel(2) * ps%p(n)%vel(3)
      i = ps%p(n)%ind(1); j = ps%p(n)%ind(2); k = ps%p(n)%ind(3);
      ! why is the necessary?
      i = i + ec%fft_pg%imin_ - 1;
      j = j + ec%fft_pg%jmin_ - 1;
      k = k + ec%fft_pg%kmin_ - 1;
      ec%phi_r(i,j,k) = ec%phi_r(i,j,k) + pvol
      ec%work_r_x(i,j,k) = ec%work_r_x(i,j,k) + pvol * ps%p(n)%vel(1)
      ec%work_r_y(i,j,k) = ec%work_r_y(i,j,k) + pvol * ps%p(n)%vel(2)
      ec%work_r_z(i,j,k) = ec%work_r_z(i,j,k) + pvol * ps%p(n)%vel(3)
    end do
    ec%phi_r = ec%phi_r / ec%sim_cfg%vol_total
    ec%work_r_x = ec%work_r_x / ec%sim_cfg%vol_total
    ec%work_r_y = ec%work_r_y / ec%sim_cfg%vol_total
    ec%work_r_z = ec%work_r_z / ec%sim_cfg%vol_total
    call mpi_reduce(ps%np_, n, 1, MPI_INTEGER, MPI_SUM, 0, ec%fft_pg%comm)
    call mpi_reduce(dp_local, dp, 1, MPI_REAL_WP, MPI_SUM, 0, ec%fft_pg%comm)
    call mpi_reduce(volfrac_local, volfrac, 1, MPI_REAL_WP, MPI_SUM, 0,       &
      ec%fft_pg%comm)
    call mpi_reduce(meanVV_local, meanVV, 4, MPI_REAL_WP, MPI_SUM, 0,         &
      ec%fft_pg%comm)
    if (ec%fft_pg%rank .eq. 0) then
      dp = dp / n; volfrac = volfrac / ec%sim_cfg%vol_total;
      meanVV = meanVV / n;
    end if

    ! do forward transforms to get variables in fourier space
    call p3dfft_3Dtrans_double(ec%trans_f, ec%phi_r, ec%phi_f, 0)
    call p3dfft_3Dtrans_double(ec%trans_f, ec%work_r_x, ec%vel_f_x, 0)
    call p3dfft_3Dtrans_double(ec%trans_f, ec%work_r_y, ec%vel_f_y, 0)
    call p3dfft_3Dtrans_double(ec%trans_f, ec%work_r_z, ec%vel_f_z, 0)

    ! for each filter, convolve and compute means, write result
    flt => ec%first_filter
    do while (associated(flt))

      ! compute filtered phi
      ec%work_f_x = ec%phi_f * flt%flt_f
      call p3dfft_3Dtrans_double(ec%trans_b, ec%work_f_x, ec%phi_r, 0)
      ec%phi_r = fft_rescale * ec%phi_r
      tmpvec(1) = sum(ec%phi_r)
      call mpi_allreduce(tmpvec(1), flt%d%meanphibar, 1, MPI_REAL_WP,         &
        MPI_SUM, 0, ec%fft_pg%comm)
      ec%phi_r = ec%phi_r - flt%d%meanphibar

      ! compute filtered vel
      ec%work_f_x = ec%vel_f_x * flt%flt_f
      ec%work_f_y = ec%vel_f_y * flt%flt_f
      ec%work_f_z = ec%vel_f_z * flt%flt_f
      call p3dfft_3Dtrans_double(ec%trans_b, ec%work_f_x, ec%work_r_x, 0)
      call p3dfft_3Dtrans_double(ec%trans_b, ec%work_f_y, ec%work_r_y, 0)
      call p3dfft_3Dtrans_double(ec%trans_b, ec%work_f_z, ec%work_r_z, 0)
      ec%work_r_x = fft_rescale * ec%work_r_x
      ec%work_r_y = fft_rescale * ec%work_r_y
      ec%work_r_z = fft_rescale * ec%work_r_z
      tmpvec(:) = (/ sum(ec%work_r_x), sum(ec%work_r_y), sum(ec%work_r_z) /)
      call mpi_allreduce(tmpvec, flt%d%meanmesopvel, 3, MPI_REAL_WP, MPI_SUM, &
        0, ec%fft_pg%comm)
      ec%work_r_x = ec%work_r_x - flt%d%meanmesopvel(1)
      ec%work_r_y = ec%work_r_y - flt%d%meanmesopvel(2)
      ec%work_r_z = ec%work_r_z - flt%d%meanmesopvel(3)

      ! compute stats (phi_r and work_r_? now contain phicheck and upcheck)
      ! not all of these values are correct off of the root processor, but
      ! that's fine
      flt%d%time = time; flt%d%step = step;
      flt%d%dp = dp; flt%d%volfrac = volfrac;
      ec%work_r_x = ec%work_r_x * ec%phi_r
      ec%work_r_y = ec%work_r_y * ec%phi_r
      ec%work_r_z = ec%work_r_z * ec%phi_r
      tmpvec(:) = (/ sum(ec%work_r_x), sum(ec%work_r_y), sum(ec%work_r_z) /)
      call mpi_reduce(tmpvec, flt%d%meanphicheckupmcheck, 3, MPI_REAL_WP,     &
        MPI_SUM, 0, ec%fft_pg%comm)
      ec%work_r_x = ec%work_r_x * ec%phi_r
      ec%work_r_y = ec%work_r_y * ec%phi_r
      ec%work_r_z = ec%work_r_z * ec%phi_r
      tmpvec(:) = (/ sum(ec%work_r_x), sum(ec%work_r_y), sum(ec%work_r_z) /)
      call mpi_reduce(tmpvec, flt%d%meanphicheck2upmcheck, 3, MPI_REAL_WP,    &
        MPI_SUM, 0, ec%fft_pg%comm)
      ec%phi_r = ec%phi_r * ec%phi_r  ! in place should be faster than work_r
      tmpvec(1) = sum(ec%phi_r)
      call mpi_reduce(tmpvec(1), flt%d%meanphicheck2, 1, MPI_REAL_WP,         &
        MPI_SUM, 0, ec%fft_pg%comm)
      ! save stats
      if (ec%fft_pg%rank .eq. 0) write(*,*) "writing..."
      if (ec%fft_pg%rank .eq. 0) call flt%mon%write()
      ! move to next filter
      flt => flt%next
    end do

  end subroutine compute_statistics

end module estimclosures_class

