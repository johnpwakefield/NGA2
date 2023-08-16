!>
!> Originally written by John P Wakefield in May 2023
!> Reworked by John P Wakefield in August 2023
!>
!> Note slice is in the yz plane
!>
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
  implicit none
  private

  public :: xflwstats

  !> filtered fluctuation statistics
  type :: xflwstats
    character(len=str_medium) :: out_fname
    type(monitor) :: mon
    integer :: step
    real(WP) :: PB, PC2
    real(WP), dimension(3) :: UB, SB, PCUB, UBPCG, GPCUB  ! x, y, z
    real(WP), dimension(3,3) :: UC2  ! only upper triangle used
  contains
    procedure :: init => xflwstats_init_mon
    final :: xflwstats_destroy_mon
  end type xflwstats

  !> filter information
  type :: filter
    character(len=str_medium) :: filtername
    type(xflwstats) :: fstats
    integer :: type
    real(WP), dimension(FLT_NUM_PARAMS) :: params
    complex(WP), dimension(:,:,:), allocatable :: flt_f
    ! io optional for each filter
    logical :: use_slice_io
    real(WP), dimension(:,:,:), pointer :: io_phi
    real(WP), dimension(:,:,:,:), pointer :: io_up
  contains
    procedure :: init => filter_init
    final :: filter_destruct
  end type filter

  ! contains both simulation aligned (*_aln) and standard grids
  ! aligned grids use the sim_pg decomposition, but are unallocated where the
  ! simulation region does not overlap with the fft region
  ! the fft region is centered around the observation plane, but is much larger
  ! in particular, the 3D FFT is still used for filtering; this window must be
  !   large enough to contain the filter support on each side so periodicity
  !   does not effect the results
  !TODO left off here
  type :: filterset
    class(pgrid), pointer :: sim_pg
    real(WP), dimension(:,:,:), pointer :: rhof, visc, U, V, W, sx, sy, sz
    type(lpt), pointer :: ps
    type(sgrid) :: fft_sg
    type(pgrid) :: fft_pg
    type(cfourier), pointer :: fft
    integer :: num_filters
    type(microstats) :: mstats
    type(filter), dimension(:), pointer :: filters
    complex(WP), dimension(:,:,:), allocatable :: work_c, phi_in_f
    complex(WP), dimension(:,:,:,:), allocatable :: up_in_f, uf_in_f
    real(WP), dimension(:,:,:), allocatable :: work_r
    real(WP), dimension(:,:,:), pointer :: phi_aln, phi, phicheck
    real(WP), dimension(:,:,:,:), pointer :: up_aln, uf_aln, up, uf
    ! coupler moving from sim_pg to fft_pg arrays; note no particles are moved
    ! as the projection occurs before this step
    type(coupler), pointer :: cpl
    ! only used for slice io
    type(sgrid) :: io_sg
    type(config) :: io_cfg
    type(MPI_GROUP) :: io_grp
    logical :: in_io_grp
    type(coupler), pointer :: io_cpl
    type(npy) :: io_out
    type(lpt) :: io_lpt
    type(lptcoupler), pointer :: io_lptcpl
    type(partmesh) :: io_pmesh
  contains
    procedure :: init => filterset_init
    procedure :: init_filters => filterset_init_filters
    procedure :: setup_sliceio => filterset_setup_sliceio
    procedure :: write_sliceio => filterset_write_sliceio
    procedure :: compute_stats => filterset_compute_stats
    final :: filterset_destruct
  end type filterset


contains


  ! projects particles to the aligned arrays if within the support cutoff region
  ! computes microscopic statistics if in the micro region
  subroutine compute_micro_stats(sim_pg, fft_pg, stats, ps, step, rhof, visc,  &
    U, V, W, sx, sy, sz, indfield, pvelfield, fvelfield)
    use mpi_f08, only:   mpi_allreduce, mpi_reduce, MPI_SUM, MPI_INTEGER
    use parallel, only:  MPI_REAL_WP
    use mathtools, only: pi
    implicit none
    class(pgrid), intent(in) :: sim_pg, fft_pg
    type(microstats), intent(inout) :: stats
    type(lpt), intent(inout) :: ps
    integer, intent(in) :: step
    real(WP), dimension(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:), intent(in) :: rhof, visc, U, V, W, sx, sy, sz
    real(WP), dimension(fft_pg%imin_:,fft_pg%jmin_:,fft_pg%kmin_:), optional, intent(out) :: indfield
    real(WP), dimension(fft_pg%imin_:,fft_pg%jmin_:,fft_pg%kmin_:,1:), optional, intent(out) :: pvelfield, fvelfield
    real(WP) :: dp_loc, VF_loc, taup_loc, pvol, cellvoli
    real(WP), dimension(3) :: VB_loc, slp_loc, drg_loc, slp, drg, fvel, dxdydz
    real(WP), dimension(6) :: VV_loc, junk
    integer, dimension(3) :: ind, indmin, indmax
    integer :: n

    ! copy over step and rho
    stats%step = step
    stats%rhop = ps%rho

    dp_loc = 0.0_WP; VF_loc = 0.0_WP; taup_loc = 0.0_WP;
    VB_loc(:) = 0.0_WP; slp_loc(:) = 0.0_WP; drg_loc(:) = 0.0_WP;
    VV_loc(:) = 0.0_WP;
    if (present(indfield)) indfield(:,:,:) = 0.0_WP
    if (present(pvelfield)) pvelfield(:,:,:,:) = 0.0_WP
    if (present(fvelfield)) fvelfield(:,:,:,:) = 0.0_WP

    indmin(:) = (/ fft_pg%imin, fft_pg%jmin, fft_pg%kmin /)
    indmax(:) = (/ fft_pg%imax, fft_pg%jmax, fft_pg%kmax /)
    dxdydz(:) = (/ fft_pg%dx(indmin(1)), fft_pg%dy(indmin(2)),                &
      fft_pg%dz(indmin(3)) /)

    do n = 1, ps%np_
      ps%p(n)%ind = ps%cfg%get_ijk_global(ps%p(n)%pos, ps%p(n)%ind)
      ind = min(max(floor(ps%p(n)%pos / dxdydz), indmin), indmax)
      ind = fft_pg%get_ijk_global(ps%p(n)%pos, ind)
      dp_loc = dp_loc + ps%p(n)%d
      pvol = pi * ps%p(n)%d**3 / 6.0_WP
      VF_loc = VF_loc + pvol
      VB_loc = VB_loc + ps%p(n)%vel
      VV_loc(1) = VV_loc(1) + ps%p(n)%vel(1) * ps%p(n)%vel(1)
      VV_loc(2) = VV_loc(2) + ps%p(n)%vel(1) * ps%p(n)%vel(2)
      VV_loc(3) = VV_loc(3) + ps%p(n)%vel(1) * ps%p(n)%vel(3)
      VV_loc(4) = VV_loc(4) + ps%p(n)%vel(2) * ps%p(n)%vel(2)
      VV_loc(5) = VV_loc(5) + ps%p(n)%vel(2) * ps%p(n)%vel(3)
      VV_loc(6) = VV_loc(6) + ps%p(n)%vel(3) * ps%p(n)%vel(3)
      fvel = ps%cfg%get_velocity(pos=ps%p(n)%pos, i0=ps%p(n)%ind(1),          &
        j0=ps%p(n)%ind(2), k0=ps%p(n)%ind(3), U=U, V=V, W=W)
      slp = fvel - ps%p(n)%vel
      slp_loc = slp_loc + slp
      call ps%get_rhs(U=U, V=V, W=W, rho=rhof, visc=visc, stress_x=sx,        &
        stress_y=sy, stress_z=sz, p=ps%p(n), acc=drg,                         &
        opt_dt=junk(4))
      drg_loc = drg_loc + drg
      taup_loc = taup_loc + sum(slp * drg) / sum(drg**2)
      if (present(indfield)) indfield(ind(1),ind(2),ind(3)) =                 &
        indfield(ind(1),ind(2),ind(3)) + pvol
      if (present(pvelfield)) pvelfield(ind(1),ind(2),ind(3),:) =             &
        pvelfield(ind(1),ind(2),ind(3),:) + pvol * ps%p(n)%vel
      if (present(fvelfield)) fvelfield(ind(1),ind(2),ind(3),:) =             &
        fvelfield(ind(1),ind(2),ind(3),:) + pvol * fvel
    end do
    cellvoli = fft_pg%nx * fft_pg%ny * fft_pg%nz / fft_pg%vol_total
    if (present(indfield)) indfield = indfield * cellvoli
    if (present(pvelfield)) pvelfield = pvelfield * cellvoli
    if (present(fvelfield)) fvelfield = fvelfield * cellvoli
    call mpi_allreduce(ps%np_,   stats%Np,   1, MPI_INTEGER, MPI_SUM, sim_pg%comm)
    call mpi_allreduce(dp_loc,   stats%dp,   1, MPI_REAL_WP, MPI_SUM, sim_pg%comm)
    call mpi_allreduce(VF_loc,   stats%VF,   1, MPI_REAL_WP, MPI_SUM, sim_pg%comm)
    call mpi_allreduce(taup_loc, stats%taup, 1, MPI_REAL_WP, MPI_SUM, sim_pg%comm)
    call mpi_allreduce(VB_loc,   stats%VB,   3, MPI_REAL_WP, MPI_SUM, sim_pg%comm)
    call mpi_allreduce(slp_loc,  stats%slp,  3, MPI_REAL_WP, MPI_SUM, sim_pg%comm)
    call mpi_allreduce(drg_loc,  stats%drg,  3, MPI_REAL_WP, MPI_SUM, sim_pg%comm)
    call mpi_allreduce(VV_loc,   stats%VV,   6, MPI_REAL_WP, MPI_SUM, sim_pg%comm)
    stats%dp = stats%dp / stats%Np; stats%VF = stats%VF / sim_pg%vol_total;
    stats%taup = stats%taup / stats%Np; stats%VB = stats%VB / stats%Np;
    stats%slp = stats%slp / stats%Np; stats%drg = stats%drg / stats%Np;
    stats%VV = stats%VV / stats%Np;

    call stats%mon%write()

  end subroutine compute_micro_stats


  !> filter objects

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
    allocate(this%flt_f(fft%pg%imin_:fft%pg%imax_,fft%pg%jmin_:fft%pg%jmax_,  &
      fft%pg%kmin_:fft%pg%kmax_))

    ! fill real filter; we build fft%pg such that the corner is at origin
    select case (this%type)
    case (FLT_GAUSSIAN)
      init_gaussian: block
        real(WP) :: b, r2_z, r2_y, r2_x, xoff, yoff, zoff
        integer :: i, j, k, ishift, jshift, kshift
        this%flt_f(:,:,:) = (0.0_WP, 0.0_WP)
        b = -6.0_WP / this%params(1)**2
        do ishift = -2, 2
          xoff = ishift * fft%pg%xL
          do jshift = -2, 2
            yoff = jshift * fft%pg%yL
            do kshift = -2, 2
              zoff = kshift * fft%pg%zL
              do k = fft%pg%kmin_, fft%pg%kmax_
                r2_z = (fft%pg%zm(k) + zoff)**2
                do j = fft%pg%jmin_, fft%pg%jmax_
                  r2_y = (fft%pg%ym(j) + yoff)**2
                  do i = fft%pg%imin_, fft%pg%imax_
                    r2_x = (fft%pg%xm(i) + xoff)**2
                    this%flt_f(i,j,k) = this%flt_f(i,j,k) + exp(b * (r2_z + r2_y + r2_x))
                  end do
                end do
              end do
            end do
          end do
        end do
      end block init_gaussian
    case default
      call die("[filterset] invalid filter type")
    end select

    ! normalize filter
    flt_int_l = sum(realpart(this%flt_f))
    call mpi_allreduce(flt_int_l, flt_int_g, 1, MPI_REAL_WP, MPI_SUM,         &
      fft%pg%comm)
    this%flt_f = this%flt_f / flt_int_g

    ! transform filter
    call fft%xtransform_forward(this%flt_f)
    call fft%ytransform_forward(this%flt_f)
    call fft%ztransform_forward(this%flt_f)
    if (fft%oddball) this%flt_f(fft%pg%imin_,fft%pg%jmin_,fft%pg%kmin_) =     &
      (1.0_WP, 0.0_WP)

    this%fstats%out_fname = this%filtername
    call this%fstats%init(fft%pg%amroot)

  end subroutine filter_init

  subroutine filter_destruct(this)
    implicit none
    type(filter), intent(inout) :: this

    ! nothing to do

  end subroutine filter_destruct


  !> macroscopic statistics

  subroutine xflwstats_init_mon(this, amroot)
    implicit none
    class(xflwstats_cube), intent(inout) :: this
    logical, intent(in) :: amroot

    this%mon = monitor(amroot, this%out_fname)
    call this%mon%add_column(this%step,      'step')
    call this%mon%add_column(this%PB,        'PB')
    call this%mon%add_column(this%PC2,       'PC2')
    call this%mon%add_column(this%UB(1),     'UBx')
    call this%mon%add_column(this%UB(2),     'UBy')
    call this%mon%add_column(this%UB(3),     'UBz')
    call this%mon%add_column(this%SB(1),     'SBx')
    call this%mon%add_column(this%SB(2),     'SBy')
    call this%mon%add_column(this%SB(3),     'SBz')
    call this%mon%add_column(this%PCUB(1),   'PCUBx')
    call this%mon%add_column(this%PCUB(2),   'PCUBy')
    call this%mon%add_column(this%PCUB(3),   'PCUBz')
    call this%mon%add_column(this%UBPCG(1),  'UBPCGx')
    call this%mon%add_column(this%UBPCG(2),  'UBPCGy')
    call this%mon%add_column(this%UBPCG(3),  'UBPCGz')
    call this%mon%add_column(this%UBPCG(3),  'UBPCGz')
    call this%mon%add_column(this%UC2(1,1), 'UC2xx')
    call this%mon%add_column(this%UC2(1,2), 'UC2xy')
    call this%mon%add_column(this%UC2(1,3), 'UC2xz')
    call this%mon%add_column(this%UC2(2,2), 'UC2yy')
    call this%mon%add_column(this%UC2(2,3), 'UC2yz')
    call this%mon%add_column(this%UC2(3,3), 'UC2zz')

  end subroutine xflwstats_init_mon

  subroutine xflwstats_destroy_mon(this)
    implicit none
    type(xflwstats_cube), intent(inout) :: this

    ! nothing to do

  end subroutine xflwstats_destroy_mon

  subroutine filter_ffts_forward(fft, phi, up, uf, phi_in_f, up_in_f, uf_in_f)
    implicit none
    type(cfourier), intent(inout) :: fft
    real(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:), intent(in) :: phi
    real(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:,1:), intent(in) :: up, uf
    complex(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:), intent(out) :: phi_in_f
    complex(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:,1:), intent(out) :: up_in_f, uf_in_f
    integer :: n

    phi_in_f(:,:,:) = phi(:,:,:)
    call fft%xtransform_forward(phi_in_f)
    call fft%ytransform_forward(phi_in_f)
    call fft%ztransform_forward(phi_in_f)

    do n = 1, 3
      up_in_f(:,:,:,n) = up(:,:,:,n)
      call fft%xtransform_forward(up_in_f(:,:,:,n))
      call fft%ytransform_forward(up_in_f(:,:,:,n))
      call fft%ztransform_forward(up_in_f(:,:,:,n))
      uf_in_f(:,:,:,n) = uf(:,:,:,n)
      call fft%xtransform_forward(uf_in_f(:,:,:,n))
      call fft%ytransform_forward(uf_in_f(:,:,:,n))
      call fft%ztransform_forward(uf_in_f(:,:,:,n))
      if (fft%oddball) then
        ! we make the decision here to fix the particle velocity at zero and
        ! allow the fluid velocity to have nonzero velocity so that computing
        ! the Reynolds stress terms is easy
        uf_in_f(fft%pg%imin_,fft%pg%jmin_,fft%pg%kmin_,n) =                   &
          uf_in_f(fft%pg%imin_,fft%pg%jmin_,fft%pg%kmin_,n) -                 &
          up_in_f(fft%pg%imin_,fft%pg%jmin_,fft%pg%kmin_,n)
        up_in_f(fft%pg%imin_,fft%pg%jmin_,fft%pg%kmin_,n) =                   &
          (0.0_WP, 0.0_WP)
      end if
    end do

  end subroutine filter_ffts_forward

  ! leaving this function, phi, phicheck, up, uf contain the appropriate
  ! (filtered) quantities
  subroutine compute_filter_stats(stats, fft, flt, step, phi_in_f,       &
    up_in_f, uf_in_f, phi, phicheck, up, uf, work_r, work_c)
    use mpi_f08, only:  mpi_allreduce, MPI_SUM
    use parallel, only: MPI_REAL_WP
    implicit none
    type(xflwstats_cube), intent(inout) :: stats
    type(cfourier), intent(inout) :: fft
    type(filter), intent(in) :: flt
    integer, intent(in) :: step
    complex(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:), intent(in) :: phi_in_f
    complex(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:,1:), intent(in) :: up_in_f, uf_in_f
    real(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:), intent(out) :: phi, phicheck, work_r
    real(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:,1:), intent(out) :: up, uf
    complex(WP), dimension(fft%pg%imin_:,fft%pg%jmin_:,fft%pg%kmin_:), intent(out) :: work_c
    real(WP) :: PB_loc, PC2_loc
    real(WP), dimension(3) :: UB_loc, SB_loc, PCUB_loc, UBPCG_loc
    real(WP), dimension(3,3) :: UC2_loc
    integer :: i, j, k, n, m, N3, ierr

    N3 = fft%pg%nx * fft%pg%ny * fft%pg%nz

    ! copy over step
    stats%step = step

    ! convolve and compute backward transform
    work_c(:,:,:) = phi_in_f(:,:,:) * flt%flt_f(:,:,:)
    call fft%xtransform_backward(work_c)
    call fft%ytransform_backward(work_c)
    call fft%ztransform_backward(work_c)
    phi = realpart(work_c)
    do n = 1, 3
      work_c(:,:,:) = up_in_f(:,:,:,n) * flt%flt_f(:,:,:)
      call fft%ztransform_backward(work_c)
      call fft%ytransform_backward(work_c)
      call fft%xtransform_backward(work_c)
      up(:,:,:,n) = realpart(work_c)
      work_c(:,:,:) = uf_in_f(:,:,:,n) * flt%flt_f(:,:,:)
      call fft%ztransform_backward(work_c)
      call fft%ytransform_backward(work_c)
      call fft%xtransform_backward(work_c)
      uf(:,:,:,n) = realpart(work_c)
    end do

    ! compute PB
    PB_loc = sum(phi)
    call mpi_allreduce(PB_loc, stats%PB, 1, MPI_REAL_WP, MPI_SUM, fft%pg%comm, ierr)
    stats%PB = stats%PB / N3

    ! compute UB
    !TODO this will be exactly zero from above; remove this and make sure the correct value is output to monitor files elsewhere
    do n = 1, 3; UB_loc(n) = sum(up(:,:,:,n)); end do
    call mpi_allreduce(UB_loc, stats%UB, 3, MPI_REAL_WP, MPI_SUM, fft%pg%comm, ierr)
    stats%UB = stats%UB / N3

    ! compute Reynolds-stress-type terms
    UC2_loc(n,m) = 0.0_WP
    do n = 1, 3
      do m = n, 3
        work_r = up(:,:,:,n) * up(:,:,:,m)
        UC2_loc(n,m) = sum(work_r)
      end do
    end do
    call mpi_allreduce(UC2_loc, stats%UC2, 9, MPI_REAL_WP, MPI_SUM, fft%pg%comm, ierr)
    stats%UC2 = stats%UC2 / N3

    ! compute SB
    do n = 1, 3
      work_r = uf(:,:,:,n) - up(:,:,:,n)
      SB_loc(n) = sum(work_r)
    end do
    call mpi_allreduce(SB_loc, stats%SB, 3, MPI_REAL_WP, MPI_SUM, fft%pg%comm, ierr)
    stats%SB = stats%SB / N3

    ! compute PC, PC2
    phicheck = phi - stats%PB
    work_r = phicheck**2
    PC2_loc = sum(work_r)
    call mpi_allreduce(PC2_loc, stats%PC2, 1, MPI_REAL_WP, MPI_SUM, fft%pg%comm, ierr)
    stats%PC2 = stats%PC2 / N3

    ! compute PCUB
    do n = 1, 3
      work_r = phicheck * up(:,:,:,n)
      PCUB_loc(n) = sum(work_r)
    end do
    call mpi_allreduce(PCUB_loc, stats%PCUB, 3, MPI_REAL_WP, MPI_SUM, fft%pg%comm, ierr)
    stats%PCUB = stats%PCUB / N3

    ! take gradients and compute UBPCG
    ! (the gradient of phi and phicheck is the same for HIT, but since we have
    ! phi check we will use it)
    ! x direction
    do k = fft%pg%kmin_, fft%pg%kmax_
      do j = fft%pg%jmin_, fft%pg%jmax_
        work_c(:,j,k) = im * fft%kx(:) * phicheck(:,j,k)
      end do
    end do
    call fft%ztransform_backward(work_c)
    call fft%ytransform_backward(work_c)
    call fft%xtransform_backward(work_c)
    UBPCG_loc(1) = sum(up(:,:,:,1) * realpart(work_c))
    ! y direction
    do k = fft%pg%kmin_, fft%pg%kmax_
      do i = fft%pg%imin_, fft%pg%imax_
        work_c(i,:,k) = im * fft%ky(:) * phicheck(i,:,k)
      end do
    end do
    call fft%ztransform_backward(work_c)
    call fft%ytransform_backward(work_c)
    call fft%xtransform_backward(work_c)
    UBPCG_loc(2) = sum(up(:,:,:,2) * realpart(work_c))
    ! z direction
    do j = fft%pg%jmin_, fft%pg%jmax_
      do i = fft%pg%imin_, fft%pg%imax_
        work_c(i,j,:) = im * fft%kz(:) * phicheck(i,j,:)
      end do
    end do
    call fft%ztransform_backward(work_c)
    call fft%ytransform_backward(work_c)
    call fft%xtransform_backward(work_c)
    UBPCG_loc(3) = sum(up(:,:,:,3) * realpart(work_c))
    call mpi_allreduce(UBPCG_loc, stats%UBPCG, 3, MPI_REAL_WP, MPI_SUM,     &
      fft%pg%comm)
    stats%UBPCG = stats%UBPCG / N3

    ! write
    call stats%mon%write()

  end subroutine compute_filter_stats


  !> filterset

  subroutine filterset_init(this, sim_pg, filterfile, FFTN, rhof, visc, U, V, W, sx, sy, sz, ps)
    use param,       only: param_read
    use mpi_f08,     only: mpi_bcast, MPI_REAL, MPI_INTEGER, MPI_CHARACTER,   &
      MPI_LOGICAL, MPI_COMM_WORLD
    use parallel,    only: MPI_REAL_WP, GROUP
    use sgrid_class, only: cartesian
    use messager,    only: die
    implicit none
    class(filterset), intent(inout) :: this
    class(pgrid), target, intent(in) :: sim_pg
    character(len=str_medium), intent(in) :: filterfile
    integer, dimension(3), intent(in) :: FFTN
    type(lpt), target, intent(in) :: ps
    real(WP), intent(in), target ::                                           &
      rhof(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:),                     &
      visc(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:),                     &
      U(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:),                        &
      V(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:),                        &
      W(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:),                        &
      sx(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:),                       &
      sy(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:),                       &
      sz(sim_pg%imino_:,sim_pg%jmino_:,sim_pg%kmino_:)
    real(WP), dimension(3) :: dx
    real(WP), dimension(:), allocatable :: fftxs, fftys, fftzs
    type(filter_info_row) :: f_info_raw
    integer :: li, hi, lj, hj, lk, hk, i, fh, ierr
    logical :: use_slice_io, use_slice_xz_io

    ! store pointer to simulation config and data
    this%sim_pg => sim_pg; this%ps => ps;
    this%rhof => rhof; this%visc => visc;
    this%U => U; this%V => V; this%W => W;
    this%sx => sx; this%sy => sy; this%sz => sz;

    ! set up microscopic statistics
    this%mstats%out_fname = 'microstats'
    call this%mstats%init(sim_pg%amroot)

    ! setup ffts (they need their own, much finer, pgrid)
    dx(:) = (/ sim_pg%xL, sim_pg%yL, sim_pg%zL /) / FFTN
    allocate(fftxs(FFTN(1)+1), fftys(FFTN(2)+1), fftzs(FFTN(3)+1))
    do i = 1, FFTN(1) + 1; fftxs(i) = real(i-1,WP) * dx(1); end do;
    do i = 1, FFTN(2) + 1; fftys(i) = real(i-1,WP) * dx(2); end do;
    do i = 1, FFTN(3) + 1; fftzs(i) = real(i-1,WP) * dx(3); end do;
    this%fft_sg = sgrid(coord=cartesian, no=0, x=fftxs, y=fftys, z=fftzs,     &
      xper=.true., yper=.true., zper=.true., name='EC_FFT_G')
    this%fft_pg = pgrid(this%fft_sg, GROUP, (/ sim_pg%npx, sim_pg%npy,        &
      sim_pg%npz /))
    allocate(this%fft, source=cfourier(this%fft_pg))

    ! allocate arrays
    li = this%fft_pg%imin_; lj = this%fft_pg%jmin_; lk = this%fft_pg%kmin_;
    hi = this%fft_pg%imax_; hj = this%fft_pg%jmax_; hk = this%fft_pg%kmax_;
    allocate(this%work_r(li:hi,lj:hj,lk:hk), this%work_c(li:hi,lj:hj,lk:hk),  &
      this%phi_in_f(li:hi,lj:hj,lk:hk), this%up_in_f(li:hi,lj:hj,lk:hk,3),    &
      this%uf_in_f(li:hi,lj:hj,lk:hk,3), this%phi(li:hi,lj:hj,lk:hk),         &
      this%phicheck(li:hi, lj:hj, lk:hk), this%up(li:hi,lj:hj,lk:hk,3),       &
      this%uf(li:hi,lj:hj,lk:hk,3))

    ! read filter info (as root)
    if (sim_pg%rank .eq. 0) then
      call param_read('HS write slice', use_slice_io)
      call param_read('HS write slice', use_slice_xz_io)
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
        this%filters(i)%use_slice_io = use_slice_io
        this%filters(i)%use_slice_xz_io = use_slice_xz_io
      end do
      close(fh)
    end if

    ! broadcast filter info and setup filters
    call mpi_bcast(this%num_filters, 1, MPI_INTEGER, 0, sim_pg%comm, ierr)
    if (sim_pg%rank .ne. 0) allocate(this%filters(this%num_filters))
    do i = 1, this%num_filters
      call mpi_bcast(this%filters(i)%filtername, str_medium, MPI_CHARACTER, 0,&
        sim_pg%comm, ierr)
      this%filters(i)%filtername = trim(this%filters(i)%filtername)
      call mpi_bcast(this%filters(i)%type, 1, MPI_INTEGER, 0, sim_pg%comm, ierr)
      call mpi_bcast(this%filters(i)%params, FLT_NUM_PARAMS, MPI_REAL_WP, 0,  &
        sim_pg%comm, ierr)
      call mpi_bcast(this%filters(i)%use_slice_io, 1, MPI_LOGICAL, 0,      &
        sim_pg%comm, ierr)
      call mpi_bcast(this%filters(i)%use_slice_xz_io, 1, MPI_LOGICAL, 0,      &
        sim_pg%comm, ierr)
    end do

    ! set all processes to not be io processes until it's set up
    this%in_io_grp = .false.; this%in_io_xz_grp = .false.;

  end subroutine filterset_init

  subroutine filterset_init_filters(this)
    implicit none
    class(filterset), intent(inout) :: this
    integer :: i

    do i = 1, this%num_filters
      call this%filters(i)%init(this%fft)
    end do

  end subroutine filterset_init_filters

  subroutine filterset_setup_sliceio(this)
    use mpi_f08    !,  only: mpi_group_range_incl, MPI_LOGICAL, MPI_LOR,           &
      !mpi_allreduce, mpi_bcast
    use parallel, only: MPI_REAL_WP
    use param,    only: param_read
    use messager, only: die
    use sgrid_class, only: cartesian
    implicit none
    class(filterset), intent(inout) :: this
    type(filter), pointer :: flt
    integer :: n, li, hi, lj, hj, lk, hk, i, j, ierr
    real(WP), dimension(:,:,:), pointer :: vx, vy, vz
    real(WP) :: thickness, height
    real(WP), dimension(2) :: z
    integer, dimension(3) :: partition
    logical, dimension(this%sim_pg%nproc) :: in_group_loc, in_group
    integer, dimension(3,this%sim_pg%nproc) :: ranges
    logical :: prev

    ! the the io group we use the parallel decomposition of the first
    ! two indices used for the simulation and index 1 in the third dimension;
    ! this means things will run slightly faster if io_idx is in the first
    ! 1/npz-th of the domain and io_xz_idx in in the first 1/npy-th of the
    ! domain

    !TODO add option to save slices at lower resolution than the fft grid

    ! xy geometry
    if (this%sim_pg%amroot)                                                  &
      call param_read('HS slice thickness proportion', thickness)
    call mpi_bcast(thickness, 1, MPI_REAL_WP, 0, this%sim_pg%comm, ierr)
    thickness = this%sim_pg%zL * thickness
    height = 0.5_WP * this%sim_pg%zL / this%sim_pg%npz
    height = height + this%sim_pg%x(this%sim_pg%imin)
    z(1) = height; z(2) = height + thickness;
    this%io_sg = sgrid(coord=cartesian, no=1,                              &
      x=this%fft_pg%x(this%fft_pg%imin:this%fft_pg%imax+1),                   &
      y=this%fft_pg%y(this%fft_pg%jmin:this%fft_pg%jmax+1),                   &
      z=z, xper=.true., yper=.true., zper=.false., name='filtslice_xy')

    ! xz geometry
    if (this%sim_pg%amroot)                                                  &
      call param_read('HS slice height proportion', height)
    call mpi_bcast(height, 1, MPI_REAL_WP, 0, this%sim_pg%comm, ierr)
    height = this%sim_pg%yL * height
    height = height + this%sim_pg%y(this%sim_pg%jmin)
    z(1) = height; z(2) = height + thickness;
    this%io_xz_sg = sgrid(coord=cartesian, no=1,                              &
      x=this%fft_pg%x(this%fft_pg%imin:this%fft_pg%imax+1),                   &
      y=z, z=this%fft_pg%z(this%fft_pg%kmin:this%fft_pg%kmax+1),              &
      xper=.true., yper=.false., zper=.true., name='filtslice_xz')

    ! xy group
    in_group_loc(:) = .false.
    in_group_loc(this%sim_pg%rank + 1) = this%sim_pg%kproc .eq. 1
    call mpi_allreduce(in_group_loc, in_group, this%sim_pg%nproc,             &
      MPI_LOGICAL, MPI_LOR, this%sim_pg%comm, ierr)
    if (in_group(1)) then; n = 1; else; n = 0; end if;
    do i = 2, this%sim_pg%nproc
      if (in_group(i) .and. .not. in_group(i-1)) n = n + 1
    end do
    j = 0
    prev = .false.
    do i = 1, this%sim_pg%nproc
      if (in_group(i) .and. .not. prev) then
        j = j + 1
        ranges(1,j) = i - 1
      end if
      if (.not. in_group(i) .and. prev) then
        ranges(2,j) = i - 2
      end if
      prev = in_group(i)
    end do
    if (prev) ranges(2,j) = this%sim_pg%nproc - 1
    ranges(3,:) = 1
    if (j .ne. n) call die("[EC] error determining io xy group")
    call mpi_group_range_incl(this%sim_pg%group, n, ranges(1:3,1:n),          &
      this%io_grp, ierr)
    this%in_io_grp = in_group(this%sim_pg%rank+1)

    ! xz group
    in_group_loc(:) = .false.
    in_group_loc(this%sim_pg%rank + 1) = this%sim_pg%jproc .eq. 1
    call mpi_allreduce(in_group_loc, in_group, this%sim_pg%nproc,             &
      MPI_LOGICAL, MPI_LOR, this%sim_pg%comm, ierr)
    if (in_group(1)) then; n = 1; else; n = 0; end if;
    do i = 2, this%sim_pg%nproc
      if (in_group(i) .and. .not. in_group(i-1)) n = n + 1
    end do
    j = 0
    prev = .false.
    do i = 1, this%sim_pg%nproc
      if (in_group(i) .and. .not. prev) then
        j = j + 1
        ranges(1,j) = i - 1
      end if
      if (.not. in_group(i) .and. prev) then
        ranges(2,j) = i - 2
      end if
      prev = in_group(i)
    end do
    if (prev) ranges(2,j) = this%sim_pg%nproc - 1
    ranges(3,:) = 1
    if (j .ne. n) call die("[EC] error determining io xz group")
    call mpi_group_range_incl(this%sim_pg%group, n, ranges(1:3,1:n),          &
      this%io_xz_grp, ierr)
    this%in_io_xz_grp = in_group(this%sim_pg%rank+1)

    ! xy partition
    partition(:) = (/ this%sim_pg%npx, this%sim_pg%npy, 1 /)
    if (this%in_io_grp) then
      this%io_cfg = config(grp=this%io_grp, decomp=partition, grid=this%io_sg)
      this%io_lpt = lpt(cfg=this%io_cfg, name='hitslicexy')
      this%io_cfg%VF = 1.0_WP
      this%io_pmesh = partmesh(nvar=2, nvec=2, name='lpt')
      this%io_pmesh%varname(1) = "id"
      this%io_pmesh%varname(2) = "dp"
      this%io_pmesh%vecname(1) = "vel"
      this%io_pmesh%vecname(2) = "fld"
    end if
    allocate(this%io_cpl, source=coupler(src_grp=this%sim_pg%group,           &
      dst_grp=this%io_grp, name='hitfiltslicexy'))
    allocate(this%io_lptcpl, source=lptcoupler(src_grp=this%sim_pg%group,     &
      dst_grp=this%io_grp, name='hitpartslicexy'))
    call this%io_cpl%set_src(this%fft_pg)
    call this%io_lptcpl%set_src(this%ps)
    if (this%in_io_grp) then
      call this%io_cpl%set_dst(this%io_cfg)
      call this%io_lptcpl%set_dst(this%io_lpt)
    end if
    call this%io_cpl%initialize()
    call this%io_lptcpl%initialize()
    if (this%in_io_grp) then
      this%io_out = npy(pg=this%io_cfg, folder='hitslicexy')
      call this%io_out%add_particle('partslice', this%io_pmesh)
    end if

    ! xz partition
    partition(:) = (/ this%sim_pg%npx, 1, this%sim_pg%npz /)
    if (this%in_io_xz_grp) then
      this%io_xz_cfg = config(grp=this%io_xz_grp, decomp=partition, grid=this%io_xz_sg)
      this%io_xz_lpt = lpt(cfg=this%io_xz_cfg, name='hitslicexz')
      this%io_xz_cfg%VF = 1.0_WP
      this%io_xz_pmesh = partmesh(nvar=2, nvec=2, name='lpt')
      this%io_xz_pmesh%varname(1) = "id"
      this%io_xz_pmesh%varname(2) = "dp"
      this%io_xz_pmesh%vecname(1) = "vel"
      this%io_xz_pmesh%vecname(2) = "fld"
    end if
    allocate(this%io_xz_cpl, source=coupler(src_grp=this%sim_pg%group,           &
      dst_grp=this%io_xz_grp, name='hitfiltslicexz'))
    allocate(this%io_xz_lptcpl, source=lptcoupler(src_grp=this%sim_pg%group,     &
      dst_grp=this%io_xz_grp, name='hitpartslicexz'))
    call this%io_xz_cpl%set_src(this%fft_pg)
    call this%io_xz_lptcpl%set_src(this%ps)
    if (this%in_io_xz_grp) then
      call this%io_xz_cpl%set_dst(this%io_xz_cfg)
      call this%io_xz_lptcpl%set_dst(this%io_xz_lpt)
    end if
    call this%io_xz_cpl%initialize()
    call this%io_xz_lptcpl%initialize()
    if (this%in_io_xz_grp) then
      this%io_xz_out = npy(pg=this%io_xz_cfg, folder='hitslicexz')
      call this%io_xz_out%add_particle('partslice', this%io_xz_pmesh)
    end if

    ! xy init filter memory
    do n = 1, this%num_filters
      flt => this%filters(n)
      if (.not. flt%use_slice_io .or. .not. this%in_io_grp) cycle
      li = this%io_cfg%imino_; hi = this%io_cfg%imaxo_;
      lj = this%io_cfg%jmino_; hj = this%io_cfg%jmaxo_;
      lk = this%io_cfg%kmino_; hk = this%io_cfg%kmaxo_;
      allocate(flt%io_phi(li:hi,lj:hj,lk:hk),                                &
        flt%io_up(li:hi,lj:hj,lk:hk,3), flt%io_uf(li:hi,lj:hj,lk:hk,3))
      call this%io_out%add_scalar(trim(flt%filtername) // 'phi', flt%io_phi)
      vx => flt%io_up(:,:,:,1); vy => flt%io_up(:,:,:,2); vz => flt%io_up(:,:,:,3);
      call this%io_out%add_vector(trim(flt%filtername) // 'up', vx, vy, vz)
      vx => flt%io_uf(:,:,:,1); vy => flt%io_uf(:,:,:,2); vz => flt%io_uf(:,:,:,3);
      call this%io_out%add_vector(trim(flt%filtername) // 'uf', vx, vy, vz)
    end do

    ! xz init filter memory
    do n = 1, this%num_filters
      flt => this%filters(n)
      if (.not. flt%use_slice_xz_io .or. .not. this%in_io_xz_grp) cycle
      li = this%io_xz_cfg%imino_; hi = this%io_xz_cfg%imaxo_;
      lj = this%io_xz_cfg%jmino_; hj = this%io_xz_cfg%jmaxo_;
      lk = this%io_xz_cfg%kmino_; hk = this%io_xz_cfg%kmaxo_;
      allocate(flt%io_xz_phi(li:hi,lj:hj,lk:hk),                                &
        flt%io_xz_up(li:hi,lj:hj,lk:hk,3), flt%io_xz_uf(li:hi,lj:hj,lk:hk,3))
      call this%io_xz_out%add_scalar(trim(flt%filtername) // 'phi', flt%io_xz_phi)
      vx => flt%io_xz_up(:,:,:,1); vy => flt%io_xz_up(:,:,:,2); vz => flt%io_xz_up(:,:,:,3);
      call this%io_xz_out%add_vector(trim(flt%filtername) // 'up', vx, vy, vz)
      vx => flt%io_xz_uf(:,:,:,1); vy => flt%io_xz_uf(:,:,:,2); vz => flt%io_xz_uf(:,:,:,3);
      call this%io_xz_out%add_vector(trim(flt%filtername) // 'uf', vx, vy, vz)
    end do

  end subroutine filterset_setup_sliceio

  subroutine filterset_write_sliceio(this, t)
    implicit none
    class(filterset), intent(inout) :: this
    real(WP), intent(in) :: t

    if (this%in_io_grp) call this%io_out%write_data(t)
    if (this%in_io_xz_grp) call this%io_xz_out%write_data(t)

  end subroutine filterset_write_sliceio

  subroutine filterset_compute_stats(this, step)
    implicit none
    class(filterset), intent(inout) :: this
    integer, intent(in) :: step
    integer, dimension(3) :: ind
    integer :: m, n, i

    call compute_micro_stats(this%sim_pg, this%fft%pg, this%mstats, this%ps,  &
      step, this%rhof, this%visc, this%U, this%V, this%W, this%sx, this%sy,   &
      this%sz, this%phi, this%up, this%uf)

    call filter_ffts_forward(this%fft, this%phi, this%up, this%uf,            &
      this%phi_in_f, this%up_in_f, this%uf_in_f)

    do n = 1, this%num_filters

      call compute_filter_stats(this%filters(n)%fstats, this%fft,             &
        this%filters(n), step, this%phi_in_f, this%up_in_f, this%uf_in_f,     &
        this%phi, this%phicheck, this%up, this%uf, this%work_r, this%work_c)

      ! xy slice first
      if (this%filters(n)%use_slice_io) then
        ! transfer continuum properties
        call this%io_cpl%push(this%phi)
        call this%io_cpl%transfer()
        if (this%in_io_grp) call this%io_cpl%pull(this%filters(n)%io_phi)
        do m = 1, 3
          call this%io_cpl%push(this%up(:,:,:,m))
          call this%io_cpl%transfer()
          if (this%in_io_grp) call this%io_cpl%pull(this%filters(n)%io_up(:,:,:,m))
          call this%io_cpl%push(this%uf(:,:,:,m))
          call this%io_cpl%transfer()
          if (this%in_io_grp) call this%io_cpl%pull(this%filters(n)%io_uf(:,:,:,m))
        end do
        ! transfer particle properties to lpt
        call this%io_lptcpl%push()
        call this%io_lptcpl%transfer()
        if (this%in_io_grp) call this%io_lptcpl%pull()
        ! transfer particles to pmesh
        if (this%in_io_grp) then
          call this%io_pmesh%reset()
          call this%io_pmesh%set_size(this%io_lpt%np_)
          do i = 1, this%io_lpt%np_
            this%io_pmesh%pos(:,i) = this%io_lpt%p(i)%pos
            this%io_pmesh%var(1,i) = this%io_lpt%p(i)%id
            this%io_pmesh%var(2,i) = this%io_lpt%p(i)%d
            this%io_pmesh%vec(:,1,i) = this%io_lpt%p(i)%vel
            ind = this%ps%cfg%get_ijk_global(this%io_lpt%p(i)%pos)
            this%io_pmesh%vec(:,2,i) = this%ps%cfg%get_velocity(             &
              pos=this%io_lpt%p(i)%pos, i0=ind(1), j0=ind(2), k0=ind(3),     &
              U=this%U, V=this%V, W=this%W)
          end do
        end if
      end if

      ! xz slice second
      if (this%filters(n)%use_slice_xz_io) then
        ! transfer coninuum properties
        call this%io_xz_cpl%push(this%phi)
        call this%io_xz_cpl%transfer()
        if (this%in_io_xz_grp) call this%io_xz_cpl%pull(this%filters(n)%io_xz_phi)
        do m = 1, 3
          call this%io_xz_cpl%push(this%up(:,:,:,m))
          call this%io_xz_cpl%transfer()
          if (this%in_io_xz_grp) call this%io_xz_cpl%pull(this%filters(n)%io_xz_up(:,:,:,m))
          call this%io_xz_cpl%push(this%uf(:,:,:,m))
          call this%io_xz_cpl%transfer()
          if (this%in_io_xz_grp) call this%io_xz_cpl%pull(this%filters(n)%io_xz_uf(:,:,:,m))
        end do
        ! transfer particle properties to lpt
        call this%io_xz_lptcpl%push()
        call this%io_xz_lptcpl%transfer()
        if (this%in_io_xz_grp) call this%io_xz_lptcpl%pull()
        ! transfer particles to pmesh
        if (this%in_io_xz_grp) then
          call this%io_xz_pmesh%reset()
          call this%io_xz_pmesh%set_size(this%io_xz_lpt%np_)
          do i = 1, this%io_xz_lpt%np_
            this%io_xz_pmesh%pos(:,i) = this%io_xz_lpt%p(i)%pos
            this%io_xz_pmesh%var(1,i) = this%io_xz_lpt%p(i)%id
            this%io_xz_pmesh%var(2,i) = this%io_xz_lpt%p(i)%d
            this%io_xz_pmesh%vec(:,1,i) = this%io_xz_lpt%p(i)%vel
            ind = this%ps%cfg%get_ijk_global(this%io_xz_lpt%p(i)%pos)
            this%io_xz_pmesh%vec(:,2,i) = this%ps%cfg%get_velocity(               &
              pos=this%io_xz_lpt%p(i)%pos, i0=ind(1), j0=ind(2), k0=ind(3),       &
              U=this%U, V=this%V, W=this%W)
          end do
        end if
      end if

    end do

  end subroutine filterset_compute_stats

  subroutine filterset_destruct(this)
    use messager, only: die
    implicit none
    type(filterset), intent(inout) :: this
    type(filter), pointer :: flt
    integer :: n

    deallocate(this%phi, this%phicheck, this%up, this%uf, this%phi_in_f,       &
      this%up_in_f, this%uf_in_f, this%work_r, this%work_c)

    do n = 1, this%num_filters
      flt => this%filters(n)
      if (flt%use_slice_io) deallocate(flt%io_phi, flt%io_up, flt%io_uf)
      if (flt%use_slice_xz_io) deallocate(flt%io_xz_phi, flt%io_xz_up, flt%io_xz_uf)
    end do

  end subroutine filterset_destruct

end module xflwstats_class

