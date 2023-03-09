!> FFT for 3D periodic uniform computational domains decomposed in at most 2
!> directions. Makes use of FFTW and in-house parallel transpose operations.
!>
!> The FFT-based pressure solver that uses this class internally may be found
!> in fftsolver_clas.f90
!>
!> Unlike FFTW (and several other libaries), the transforms provided by this
!> class have the correct scaling and sign.  `forward_transform` and
!> `backward_transform` are truly inverses; multiplication of the inverse
!> transform result by -Nx*Ny*Nz is not necessary.  In addition, the real-space
!> integral/Fourier-space zero value is preserved through forward and inverse
!> transforms.
!>
module fft3d_class
  use precision,    only: WP
  use pgrid_class,  only: pgrid
  use string,       only: str_short
  use, intrinsic :: iso_c_binding
  implicit none
  private


  ! Expose type/constructor/methods
  public :: fft3d


  !> fft3d object definition
  type :: fft3d

    ! pgrid
    class(pgrid), pointer :: pg

    ! FFT's oddball
    logical :: oddball

    ! Data storage for FFTW plans
    complex(WP), dimension(:), allocatable :: in_x,out_x
    complex(WP), dimension(:), allocatable :: in_y,out_y
    complex(WP), dimension(:), allocatable :: in_z,out_z

    ! FFTW plans
    type(C_PTR) :: fplan_x,bplan_x
    type(C_PTR) :: fplan_y,bplan_y
    type(C_PTR) :: fplan_z,bplan_z

    !> Unstrided arrays
    complex(WP), dimension(:,:,:), allocatable :: factored_operator
    complex(WP), dimension(:,:,:), allocatable :: transformed_rhs

    ! Storage for transposed data
    complex(WP), dimension(:,:,:), allocatable :: xtrans, ytrans, ztrans

    ! Transpose partition - X
    integer, dimension(:), allocatable :: imin_x,imax_x
    integer, dimension(:), allocatable :: jmin_x,jmax_x
    integer, dimension(:), allocatable :: kmin_x,kmax_x
    integer, dimension(:), allocatable :: nx_x,ny_x,nz_x
    complex(WP), dimension(:,:,:,:), allocatable :: sendbuf_x,recvbuf_x
    integer :: sendcount_x,recvcount_x
    character :: xdir

    ! Transpose partition - Y
    integer, dimension(:), allocatable :: imin_y,imax_y
    integer, dimension(:), allocatable :: jmin_y,jmax_y
    integer, dimension(:), allocatable :: kmin_y,kmax_y
    integer, dimension(:), allocatable :: nx_y,ny_y,nz_y
    complex(WP), dimension(:,:,:,:), allocatable :: sendbuf_y,recvbuf_y
    integer :: sendcount_y,recvcount_y
    character :: ydir

    ! Transpose partition - Z
    integer, dimension(:), allocatable :: imin_z,imax_z
    integer, dimension(:), allocatable :: jmin_z,jmax_z
    integer, dimension(:), allocatable :: kmin_z,kmax_z
    integer, dimension(:), allocatable :: nx_z,ny_z,nz_z
    complex(WP), dimension(:,:,:,:), allocatable :: sendbuf_z,recvbuf_z
    integer :: sendcount_z,recvcount_z
    character :: zdir

  contains

    procedure :: print_short=>fft3d_print_short   !< One-line printing of solver status
    procedure :: print=>fft3d_print               !< Long-form printing of solver status
    procedure :: log=>fft3d_log                   !< Long-form logging of solver status
    final :: fft3d_destroy                        !< Deconstructor

    procedure :: forward_transform => fft3d_fourier_transform
    procedure :: backward_transform => fft3d_inverse_transform

    procedure, private :: fft3d_xtranspose_init
    procedure, private :: fft3d_ytranspose_init
    procedure, private :: fft3d_ztranspose_init

    procedure, private :: fft3d_xtranspose_forward
    procedure, private :: fft3d_ytranspose_forward
    procedure, private :: fft3d_ztranspose_forward

    procedure, private :: fft3d_xtranspose_backward
    procedure, private :: fft3d_ytranspose_backward
    procedure, private :: fft3d_ztranspose_backward

  end type fft3d


  !> Declare fft3d constructor
  interface fft3d; procedure fft3d_from_args; end interface fft3d;


contains


  !> Constructor for a fft3d object
  !>  (pgrid is a subclass of pgrid; a full cfg object can be passed here
  !>  without issue)
  function fft3d_from_args(pg) result(self)
    use messager, only: die
    implicit none
    type(fft3d) :: self
    class(pgrid), target, intent(in) :: pg
    include 'fftw3.f03'

    ! Link the config and store the name
    self%pg=>pg

    ! Various checks to ensure we can use this solver
    check_solver_is_useable: block
      integer :: ndim,ndcp
      ! Periodicity and uniformity of mesh
      if (self%pg%nx.gt.1.and..not.(self%pg%xper.and.self%pg%uniform_x)) &
        call die('[fft3d constructor] Need x-direction needs to be &
        &periodic and uniform')
      if (self%pg%ny.gt.1.and..not.(self%pg%yper.and.self%pg%uniform_y)) &
        call die('[fft3d constructor] Need y-direction needs to be &
        &periodic and uniform')
      if (self%pg%nz.gt.1.and..not.(self%pg%zper.and.self%pg%uniform_z)) &
        call die('[fft3d constructor] Need z-direction needs to be &
        &periodic and uniform')
      ! Ensure that we have at least one non-decomposed direction
      ndim = count((/ self%pg%nx, self%pg%ny, self%pg%nz /) .gt. 1)
      ndcp = count((/ self%pg%npx, self%pg%npy, self%pg%npz /) .gt. 1)
      if (ndcp.ge.ndim) call die('[fft3d constructor] Need at least one &
        &NON-decomposed direction')
    end block check_solver_is_useable

    ! Initialize transpose and FFTW plans in x
    if (self%pg%nx.gt.1) then
      call self%fft3d_xtranspose_init()
      allocate(self%in_x(self%pg%nx),self%out_x(self%pg%nx))
      self%fplan_x=fftw_plan_dft_1d(self%pg%nx,self%in_x,self%out_x,FFTW_FORWARD,FFTW_MEASURE)
      self%bplan_x=fftw_plan_dft_1d(self%pg%nx,self%in_x,self%out_x,FFTW_BACKWARD,FFTW_MEASURE)
    end if

    ! Initialize transpose and FFTW plans in y
    if (self%pg%ny.gt.1) then
      call self%fft3d_ytranspose_init()
      allocate(self%in_y(self%pg%ny),self%out_y(self%pg%ny))
      self%fplan_y=fftw_plan_dft_1d(self%pg%ny,self%in_y,self%out_y,FFTW_FORWARD,FFTW_MEASURE)
      self%bplan_y=fftw_plan_dft_1d(self%pg%ny,self%in_y,self%out_y,FFTW_BACKWARD,FFTW_MEASURE)
    end if

    ! Initialize transpose and FFTW plans in z
    if (self%pg%nz.gt.1) then
      call self%fft3d_ztranspose_init()
      allocate(self%in_z(self%pg%nz),self%out_z(self%pg%nz))
      self%fplan_z=fftw_plan_dft_1d(self%pg%nz,self%in_z,self%out_z,FFTW_FORWARD,FFTW_MEASURE)
      self%bplan_z=fftw_plan_dft_1d(self%pg%nz,self%in_z,self%out_z,FFTW_BACKWARD,FFTW_MEASURE)
    end if

    ! Find who owns the oddball
    self%oddball=all((/ self%pg%iproc, self%pg%jproc, self%pg%kproc /) .eq. 1)

  end function fft3d_from_args


  !> Initialize transpose tool in x
  subroutine fft3d_xtranspose_init(this)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    integer :: ierr,ip,q,r

    ! Determine non-decomposed direction to use for transpose
    if      (this%pg%npx.eq.1.and.this%pg%nx.gt.1) then
      this%xdir='x'
    else if (this%pg%npy.eq.1.and.this%pg%ny.gt.1) then
      this%xdir='y'
    else if (this%pg%npz.eq.1.and.this%pg%nz.gt.1) then
      this%xdir='z'
    end if

    ! Allocate global partitions
    allocate(  this%nx_x(this%pg%npx))
    allocate(  this%ny_x(this%pg%npx))
    allocate(  this%nz_x(this%pg%npx))
    allocate(this%imin_x(this%pg%npx))
    allocate(this%imax_x(this%pg%npx))
    allocate(this%jmin_x(this%pg%npx))
    allocate(this%jmax_x(this%pg%npx))
    allocate(this%kmin_x(this%pg%npx))
    allocate(this%kmax_x(this%pg%npx))

    ! Partition
    select case (trim(this%xdir))
    case ('x')

      ! No transpose required, use local partition
      this%nx_x=this%pg%nx_
      this%ny_x=this%pg%ny_
      this%nz_x=this%pg%nz_
      this%imin_x=this%pg%imin_
      this%imax_x=this%pg%imax_
      this%jmin_x=this%pg%jmin_
      this%jmax_x=this%pg%jmax_
      this%kmin_x=this%pg%kmin_
      this%kmax_x=this%pg%kmax_

    case ('y')

      ! Store old local indices from each processor
      call MPI_AllGather(this%pg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
      call MPI_AllGather(this%pg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
      this%nx_x=this%imax_x-this%imin_x+1

      ! Partition new local indices
      do ip=1,this%pg%npx
        q=this%pg%ny/this%pg%npx
        r=mod(this%pg%ny,this%pg%npx)
        if (ip.le.r) then
          this%ny_x(ip)  =q+1
          this%jmin_x(ip)=this%pg%jmin+(ip-1)*(q+1)
        else
          this%ny_x(ip)  =q
          this%jmin_x(ip)=this%pg%jmin+r*(q+1)+(ip-r-1)*q
        end if
        this%jmax_x(ip)=this%jmin_x(ip)+this%ny_x(ip)-1
      end do
      this%nz_x=this%pg%nz_
      this%kmin_x=this%pg%kmin_
      this%kmax_x=this%pg%kmax_

      ! Variables for AllToAll communication
      this%sendcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%pg%nz_
      this%recvcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%pg%nz_
      allocate(this%sendbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%pg%kmin_:this%pg%kmax_,this%pg%npx))
      allocate(this%recvbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%pg%kmin_:this%pg%kmax_,this%pg%npx))

      ! Zero out buffers
      this%sendbuf_x=0.0_WP
      this%recvbuf_x=0.0_WP

    case ('z')

      ! Store old local indices from each processor
      call MPI_AllGather(this%pg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
      call MPI_AllGather(this%pg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
      this%nx_x=this%imax_x-this%imin_x+1

      ! Partition new local indices
      do ip=1,this%pg%npx
        q=this%pg%nz/this%pg%npx
        r=mod(this%pg%nz,this%pg%npx)
        if (ip.le.r) then
          this%nz_x(ip)  =q+1
          this%kmin_x(ip)=this%pg%kmin+(ip-1)*(q+1)
        else
          this%nz_x(ip)  =q
          this%kmin_x(ip)=this%pg%kmin+r*(q+1)+(ip-r-1)*q
        end if
        this%kmax_x(ip)=this%kmin_x(ip)+this%nz_x(ip)-1
      end do
      this%ny_x=this%pg%ny_
      this%jmin_x=this%pg%jmin_
      this%jmax_x=this%pg%jmax_

      ! Variables for AllToAll communication
      this%sendcount_x=maxval(this%nx_x)*this%pg%ny_*maxval(this%nz_x)
      this%recvcount_x=maxval(this%nx_x)*this%pg%ny_*maxval(this%nz_x)
      allocate(this%sendbuf_x(maxval(this%nx_x),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_x),this%pg%npx))
      allocate(this%recvbuf_x(maxval(this%nx_x),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_x),this%pg%npx))

      ! Zero out buffers
      this%sendbuf_x=0.0_WP
      this%recvbuf_x=0.0_WP

    end select

    ! Allocate storage
    allocate(this%xtrans(this%pg%imin:this%pg%imax,this%jmin_x(this%pg%iproc):this%jmax_x(this%pg%iproc),this%kmin_x(this%pg%iproc):this%kmax_x(this%pg%iproc)))

  end subroutine fft3d_xtranspose_init


  !> Initialize transpose tool in y
  subroutine fft3d_ytranspose_init(this)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    integer :: ierr,jp,q,r

    ! Determine non-decomposed direction to use for transpose
    if      (this%pg%npy.eq.1.and.this%pg%ny.gt.1) then
      this%ydir='y'
    else if (this%pg%npz.eq.1.and.this%pg%nz.gt.1) then
      this%ydir='z'
    else if (this%pg%npx.eq.1.and.this%pg%nx.gt.1) then
      this%ydir='x'
    end if

    ! Allocate global partitions
    allocate(  this%nx_y(this%pg%npy))
    allocate(  this%ny_y(this%pg%npy))
    allocate(  this%nz_y(this%pg%npy))
    allocate(this%imin_y(this%pg%npy))
    allocate(this%imax_y(this%pg%npy))
    allocate(this%jmin_y(this%pg%npy))
    allocate(this%jmax_y(this%pg%npy))
    allocate(this%kmin_y(this%pg%npy))
    allocate(this%kmax_y(this%pg%npy))

    ! Partition
    select case (trim(this%ydir))
    case ('x')

      ! Store old local indices from each processor
      call MPI_AllGather(this%pg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
      call MPI_AllGather(this%pg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
      this%ny_y=this%jmax_y-this%jmin_y+1

      ! Partition new local indices
      do jp=1,this%pg%npy
        q=this%pg%nx/this%pg%npy
        r=mod(this%pg%nx,this%pg%npy)
        if (jp.le.r) then
          this%nx_y(jp)  =q+1
          this%imin_y(jp)=this%pg%imin+(jp-1)*(q+1)
        else
          this%nx_y(jp)  =q
          this%imin_y(jp)=this%pg%imin+r*(q+1)+(jp-r-1)*q
        end if
        this%imax_y(jp)=this%imin_y(jp)+this%nx_y(jp)-1
      end do
      this%nz_y=this%pg%nz_
      this%kmin_y=this%pg%kmin_
      this%kmax_y=this%pg%kmax_

      ! Variables for AllToAll communication
      this%sendcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%pg%nz_
      this%recvcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%pg%nz_
      allocate(this%sendbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%pg%kmin_:this%pg%kmax_,this%pg%npy))
      allocate(this%recvbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%pg%kmin_:this%pg%kmax_,this%pg%npy))

      ! Zero out buffers
      this%sendbuf_y=0.0_WP
      this%recvbuf_y=0.0_WP

    case ('y')

      ! No transpose required, use local partition
      this%nx_y=this%pg%nx_
      this%ny_y=this%pg%ny_
      this%nz_y=this%pg%nz_
      this%imin_y=this%pg%imin_
      this%imax_y=this%pg%imax_
      this%jmin_y=this%pg%jmin_
      this%jmax_y=this%pg%jmax_
      this%kmin_y=this%pg%kmin_
      this%kmax_y=this%pg%kmax_

    case ('z')

      ! Store old local indices from each processor
      call MPI_AllGather(this%pg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
      call MPI_AllGather(this%pg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
      this%ny_y=this%jmax_y-this%jmin_y+1

      ! Partition new local indices
      do jp=1,this%pg%npy
        q=this%pg%nz/this%pg%npy
        r=mod(this%pg%nz,this%pg%npy)
        if (jp.le.r) then
          this%nz_y(jp)  =q+1
          this%kmin_y(jp)=this%pg%kmin+(jp-1)*(q+1)
        else
          this%nz_y(jp)  =q
          this%kmin_y(jp)=this%pg%kmin+r*(q+1)+(jp-r-1)*q
        end if
        this%kmax_y(jp)=this%kmin_y(jp)+this%nz_y(jp)-1
      end do
      this%nx_y=this%pg%nx_
      this%imin_y=this%pg%imin_
      this%imax_y=this%pg%imax_

      ! Variables for AllToAll communication
      this%sendcount_y=this%pg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
      this%recvcount_y=this%pg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
      allocate(this%sendbuf_y(this%pg%imin_:this%pg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%pg%npy))
      allocate(this%recvbuf_y(this%pg%imin_:this%pg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%pg%npy))

      ! Zero out buffers
      this%sendbuf_y=0.0_WP
      this%recvbuf_y=0.0_WP

    end select

    ! Allocate storage
    allocate(this%ytrans(this%imin_y(this%pg%jproc):this%imax_y(this%pg%jproc),this%pg%jmin:this%pg%jmax,this%kmin_y(this%pg%jproc):this%kmax_y(this%pg%jproc)))

  end subroutine fft3d_ytranspose_init


  !> Initialize transpose tool in z
  subroutine fft3d_ztranspose_init(this)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    integer :: ierr,kp,q,r

    ! Determine non-decomposed direction to use for transpose
    if      (this%pg%npz.eq.1.and.this%pg%nz.gt.1) then
      this%zdir='z'
    else if (this%pg%npx.eq.1.and.this%pg%nx.gt.1) then
      this%zdir='x'
    else if (this%pg%npy.eq.1.and.this%pg%ny.gt.1) then
      this%zdir='y'
    end if

    ! Allocate global partitions
    allocate(  this%nx_z(this%pg%npz))
    allocate(  this%ny_z(this%pg%npz))
    allocate(  this%nz_z(this%pg%npz))
    allocate(this%imin_z(this%pg%npz))
    allocate(this%imax_z(this%pg%npz))
    allocate(this%jmin_z(this%pg%npz))
    allocate(this%jmax_z(this%pg%npz))
    allocate(this%kmin_z(this%pg%npz))
    allocate(this%kmax_z(this%pg%npz))

    ! Partition
    select case (trim(this%zdir))
    case ('x')

      ! Store old local indices from each processor
      call MPI_AllGather(this%pg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
      call MPI_AllGather(this%pg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
      this%nz_z=this%kmax_z-this%kmin_z+1

      ! Partition new local indices
      do kp=1,this%pg%npz
        q=this%pg%nx/this%pg%npz
        r=mod(this%pg%nx,this%pg%npz)
        if (kp.le.r) then
          this%nx_z(kp)  =q+1
          this%imin_z(kp)=this%pg%imin+(kp-1)*(q+1)
        else
          this%nx_z(kp)  =q
          this%imin_z(kp)=this%pg%imin+r*(q+1)+(kp-r-1)*q
        end if
        this%imax_z(kp)=this%imin_z(kp)+this%nx_z(kp)-1
      end do
      this%ny_z=this%pg%ny_
      this%jmin_z=this%pg%jmin_
      this%jmax_z=this%pg%jmax_

      ! Variables for AllToAll communication
      this%sendcount_z=maxval(this%nx_z)*this%pg%ny_*maxval(this%nz_z)
      this%recvcount_z=maxval(this%nx_z)*this%pg%ny_*maxval(this%nz_z)
      allocate(this%sendbuf_z(maxval(this%nx_z),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_z),this%pg%npz))
      allocate(this%recvbuf_z(maxval(this%nx_z),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_z),this%pg%npz))

      ! Zero out buffers
      this%sendbuf_z=0.0_WP
      this%recvbuf_z=0.0_WP

    case ('y')

      ! Store old local indices from each processor
      call MPI_AllGather(this%pg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
      call MPI_AllGather(this%pg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
      this%nz_z=this%kmax_z-this%kmin_z+1

      ! Partition new local indices
      do kp=1,this%pg%npz
        q=this%pg%ny/this%pg%npz
        r=mod(this%pg%ny,this%pg%npz)
        if (kp.le.r) then
          this%ny_z(kp)  =q+1
          this%jmin_z(kp)=this%pg%jmin+(kp-1)*(q+1)
        else
          this%ny_z(kp)  =q
          this%jmin_z(kp)=this%pg%jmin+r*(q+1)+(kp-r-1)*q
        end if
        this%jmax_z(kp)=this%jmin_z(kp)+this%ny_z(kp)-1
      end do
      this%nx_z=this%pg%nx_
      this%imin_z=this%pg%imin_
      this%imax_z=this%pg%imax_

      ! Variables for AllToAll communication
      this%sendcount_z=this%pg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
      this%recvcount_z=this%pg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
      allocate(this%sendbuf_z(this%pg%imin_:this%pg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%pg%npz))
      allocate(this%recvbuf_z(this%pg%imin_:this%pg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%pg%npz))

      ! Zero out buffers
      this%sendbuf_z=0.0_WP
      this%recvbuf_z=0.0_WP

    case ('z')

      ! No transpose required, use local partition
      this%nx_z=this%pg%nx_
      this%ny_z=this%pg%ny_
      this%nz_z=this%pg%nz_
      this%imin_z=this%pg%imin_
      this%imax_z=this%pg%imax_
      this%jmin_z=this%pg%jmin_
      this%jmax_z=this%pg%jmax_
      this%kmin_z=this%pg%kmin_
      this%kmax_z=this%pg%kmax_

    end select

    ! Allocate storage
    allocate(this%ztrans(this%imin_z(this%pg%kproc):this%imax_z(this%pg%kproc),this%jmin_z(this%pg%kproc):this%jmax_z(this%pg%kproc),this%pg%kmin:this%pg%kmax))

  end subroutine fft3d_ztranspose_init


  !> Perform forward transpose in x
  subroutine fft3d_xtranspose_forward(this,A,At)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(in) :: A
    complex(WP), dimension(this%pg%imin :,this%jmin_x(this%pg%iproc):,this%kmin_x(this%pg%iproc):), intent(out) :: At
    integer :: i,j,k,ip,ii,jj,kk,ierr

    select case (trim(this%xdir))
    case ('x')
      ! No transpose required
      At=A
    case ('y')
      ! Transpose x=>y
      do ip=1,this%pg%npx
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_x(ip),this%jmax_x(ip)
            do i=this%pg%imin_,this%pg%imax_
              jj=j-this%jmin_x(ip)+1
              ii=i-this%pg%imin_+1
              this%sendbuf_x(ii,jj,k,ip)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%pg%xcomm,ierr)
      do ip=1,this%pg%npx
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
            do i=this%imin_x(ip),this%imax_x(ip)
              jj=j-this%jmin_x(this%pg%iproc)+1
              ii=i-this%imin_x(ip)+1
              At(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
            end do
          end do
        end do
      end do
    case ('z')
      ! Transpose x=>z
      do ip=1,this%pg%npx
        do k=this%kmin_x(ip),this%kmax_x(ip)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%pg%imin_,this%pg%imax_
              kk=k-this%kmin_x(ip)+1
              ii=i-this%pg%imin_+1
              this%sendbuf_x(ii,j,kk,ip)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%pg%xcomm,ierr)
      do ip=1,this%pg%npx
        do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_x(ip),this%imax_x(ip)
              kk=k-this%kmin_x(this%pg%iproc)+1
              ii=i-this%imin_x(ip)+1
              At(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
            end do
          end do
        end do
      end do
    end select

  end subroutine fft3d_xtranspose_forward


  !> Perform forward transpose in y
  subroutine fft3d_ytranspose_forward(this,A,At)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(in) :: A
    complex(WP), dimension(this%imin_y(this%pg%jproc):,this%pg%jmin:,this%kmin_y(this%pg%jproc):), intent(out) :: At
    integer :: i,j,k,jp,ii,jj,kk,ierr

    select case (trim(this%ydir))
    case ('x')
      ! Transpose y=>x
      do jp=1,this%pg%npy
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_y(jp),this%imax_y(jp)
              ii=i-this%imin_y(jp)+1
              jj=j-this%pg%jmin_+1
              this%sendbuf_y(ii,jj,k,jp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%pg%ycomm,ierr)
      do jp=1,this%pg%npy
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
              ii=i-this%imin_y(this%pg%jproc)+1
              jj=j-this%jmin_y(jp)+1
              At(i,j,k)=this%recvbuf_y(ii,jj,k,jp)
            end do
          end do
        end do
      end do
    case ('y')
      ! No transpose required
      At=A
    case ('z')
      ! Transpose y=>z
      do jp=1,this%pg%npy
        do k=this%kmin_y(jp),this%kmax_y(jp)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%pg%imin_,this%pg%imax_
              kk=k-this%kmin_y(jp)+1
              jj=j-this%pg%jmin_+1
              this%sendbuf_y(i,jj,kk,jp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%pg%ycomm,ierr)
      do jp=1,this%pg%npy
        do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%pg%imin_,this%pg%imax_
              kk=k-this%kmin_y(this%pg%jproc)+1
              jj=j-this%jmin_y(jp)+1
              At(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
            end do
          end do
        end do
      end do
    end select

  end subroutine fft3d_ytranspose_forward


  !> Perform forward transpose in z
  subroutine fft3d_ztranspose_forward(this,A,At)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(in) :: A
    complex(WP), dimension(this%imin_z(this%pg%kproc):,this%jmin_z(this%pg%kproc):,this%pg%kmin:), intent(out) :: At
    integer :: i,j,k,kp,ii,jj,kk,ierr

    select case (trim(this%zdir))
    case ('x')
      ! Transpose z=>x
      do kp=1,this%pg%npz
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_z(kp),this%imax_z(kp)
              ii=i-this%imin_z(kp)+1
              kk=k-this%pg%kmin_+1
              this%sendbuf_z(ii,j,kk,kp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%pg%zcomm,ierr)
      do kp=1,this%pg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
              ii=i-this%imin_z(this%pg%kproc)+1
              kk=k-this%kmin_z(kp)+1
              At(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
            end do
          end do
        end do
      end do
    case ('y')
      ! Transpose z=>y
      do kp=1,this%pg%npz
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_z(kp),this%jmax_z(kp)
            do i=this%pg%imin_,this%pg%imax_
              jj=j-this%jmin_z(kp)+1
              kk=k-this%pg%kmin_+1
              this%sendbuf_z(i,jj,kk,kp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%pg%zcomm,ierr)
      do kp=1,this%pg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
            do i=this%pg%imin_,this%pg%imax_
              jj=j-this%jmin_z(this%pg%kproc)+1
              kk=k-this%kmin_z(kp)+1
              At(i,j,k)=this%recvbuf_z(i,jj,kk,kp)
            end do
          end do
        end do
      end do
    case ('z')
      ! No transpose required
      At=A
    end select

  end subroutine fft3d_ztranspose_forward


  !> Perform backward transpose in x
  subroutine fft3d_xtranspose_backward(this,At,A)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%pg%imin :,this%jmin_x(this%pg%iproc):,this%kmin_x(this%pg%iproc):), intent(in) :: At
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(out) :: A
    integer :: i,j,k,ip,ii,jj,kk,ierr

    select case (trim(this%xdir))
    case ('x')
      ! No transpose required
      A=At
    case ('y')
      ! Transpose y=>x
      do ip=1,this%pg%npx
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
            do i=this%imin_x(ip),this%imax_x(ip)
              jj=j-this%jmin_x(this%pg%iproc)+1
              ii=i-this%imin_x(ip)+1
              this%sendbuf_x(ii,jj,k,ip)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%pg%xcomm,ierr)
      do ip=1,this%pg%npx
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_x(ip),this%jmax_x(ip)
            do i=this%imin_x(this%pg%iproc),this%imax_x(this%pg%iproc)
              jj=j-this%jmin_x(ip)+1
              ii=i-this%imin_x(this%pg%iproc)+1
              A(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
            end do
          end do
        end do
      end do
    case ('z')
      ! Transpose z=>x
      do ip=1,this%pg%npx
        do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_x(ip),this%imax_x(ip)
              kk=k-this%kmin_x(this%pg%iproc)+1
              ii=i-this%imin_x(ip)+1
              this%sendbuf_x(ii,j,kk,ip)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%pg%xcomm,ierr)
      do ip=1,this%pg%npx
        do k=this%kmin_x(ip),this%kmax_x(ip)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_x(this%pg%iproc),this%imax_x(this%pg%iproc)
              kk=k-this%kmin_x(ip)+1
              ii=i-this%imin_x(this%pg%iproc)+1
              A(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
            end do
          end do
        end do
      end do
    end select

  end subroutine fft3d_xtranspose_backward


  !> Perform backward transpose in y
  subroutine fft3d_ytranspose_backward(this,At,A)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%imin_y(this%pg%jproc):,this%pg%jmin:,this%kmin_y(this%pg%jproc):), intent(in) :: At
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(out) :: A
    integer :: i,j,k,jp,ii,jj,kk,ierr

    select case (trim(this%ydir))
    case ('x')
      ! Transpose x=>y
      do jp=1,this%pg%npy
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
              ii=i-this%imin_y(this%pg%jproc)+1
              jj=j-this%jmin_y(jp)+1
              this%sendbuf_y(ii,jj,k,jp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%pg%ycomm,ierr)
      do jp=1,this%pg%npy
        do k=this%pg%kmin_,this%pg%kmax_
          do j=this%jmin_y(this%pg%jproc),this%jmax_y(this%pg%jproc)
            do i=this%imin_y(jp),this%imax_y(jp)
              ii=i-this%imin_y(jp)+1
              jj=j-this%jmin_y(this%pg%jproc)+1
              A(i,j,k)=this%recvbuf_y(ii,jj,k,jp)
            end do
          end do
        end do
      end do
    case ('y')
      ! No transpose required
      A=At
    case ('z')
      ! Transpose z=>y
      do jp=1,this%pg%npy
        do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%pg%imin_,this%pg%imax_
              kk=k-this%kmin_y(this%pg%jproc)+1
              jj=j-this%jmin_y(jp)+1
              this%sendbuf_y(i,jj,kk,jp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%pg%ycomm,ierr)
      do jp=1,this%pg%npy
        do k=this%kmin_y(jp),this%kmax_y(jp)
          do j=this%jmin_y(this%pg%jproc),this%jmax_y(this%pg%jproc)
            do i=this%pg%imin_,this%pg%imax_
              kk=k-this%kmin_y(jp)+1
              jj=j-this%jmin_y(this%pg%jproc)+1
              A(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
            end do
          end do
        end do
      end do
    end select

  end subroutine fft3d_ytranspose_backward


  !> Perform backward transpose in z
  subroutine fft3d_ztranspose_backward(this,At,A)
    use mpi_f08
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%imin_z(this%pg%kproc):,this%jmin_z(this%pg%kproc):,this%pg%kmin:), intent(in) :: At
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(out) :: A
    integer :: i,j,k,kp,ii,jj,kk,ierr

    select case (trim(this%zdir))
    case ('x')
      ! Transpose x=>z
      do kp=1,this%pg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
              ii=i-this%imin_z(this%pg%kproc)+1
              kk=k-this%kmin_z(kp)+1
              this%sendbuf_z(ii,j,kk,kp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%pg%zcomm,ierr)
      do kp=1,this%pg%npz
        do k=this%kmin_z(this%pg%kproc),this%kmax_z(this%pg%kproc)
          do j=this%pg%jmin_,this%pg%jmax_
            do i=this%imin_z(kp),this%imax_z(kp)
              ii=i-this%imin_z(kp)+1
              kk=k-this%kmin_z(this%pg%kproc)+1
              A(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
            end do
          end do
        end do
      end do
    case ('y')
      ! Transpose y=>z
      do kp=1,this%pg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
            do i=this%pg%imin_,this%pg%imax_
              jj=j-this%jmin_z(this%pg%kproc)+1
              kk=k-this%kmin_z(kp)+1
              this%sendbuf_z(i,jj,kk,kp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%pg%zcomm,ierr)
      do kp=1,this%pg%npz
        do k=this%kmin_z(this%pg%kproc),this%kmax_z(this%pg%kproc)
          do j=this%jmin_z(kp),this%jmax_z(kp)
            do i=this%pg%imin_,this%pg%imax_
              jj=j-this%jmin_z(kp)+1
              kk=k-this%kmin_z(this%pg%kproc)+1
              A(i,j,k)=this%recvbuf_z(i,jj,kk,kp)
            end do
          end do
        end do
      end do
    case ('z')
      ! No transpose required
      A=At
    end select

  end subroutine fft3d_ztranspose_backward


  !> Transpose A and perform Fourier transform
  subroutine fft3d_fourier_transform(this,A)
    use messager, only: die
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
    integer :: i,j,k
    include 'fftw3.f03'

    if (this%pg%nx.gt.1) then
      ! Transpose in X
      call this%fft3d_xtranspose_forward(A,this%xtrans)
      ! Forward transform - X
      do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
        do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
          this%in_x=this%xtrans(:,j,k)
          call fftw_execute_dft(this%fplan_x,this%in_x,this%out_x)
          this%xtrans(:,j,k)=this%out_x
        end do
      end do
      ! Transpose back
      call this%fft3d_xtranspose_backward(this%xtrans,A)
    end if

    if (this%pg%ny.gt.1) then
      ! Transpose in Y
      call this%fft3d_ytranspose_forward(A,this%ytrans)
      ! Forward transform - Y
      do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
        do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
          this%in_y=this%ytrans(i,:,k)
          call fftw_execute_dft(this%fplan_y,this%in_y,this%out_y)
          this%ytrans(i,:,k)=this%out_y
        end do
      end do
      ! Transpose back
      call this%fft3d_ytranspose_backward(this%ytrans,A)
    end if

    if (this%pg%nz.gt.1) then
      ! Transpose in Z
      call this%fft3d_ztranspose_forward(A,this%ztrans)
      ! Forward transform - Z
      do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
        do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
          this%in_z=this%ztrans(i,j,:)
          call fftw_execute_dft(this%fplan_z,this%in_z,this%out_z)
          this%ztrans(i,j,:)=this%out_z
        end do
      end do
      ! Transpose back
      call this%fft3d_ztranspose_backward(this%ztrans,A)
    end if

  end subroutine fft3d_fourier_transform


  !> FFT -> real and transpose back
  subroutine fft3d_inverse_transform(this,A)
    use messager, only: die
    implicit none
    class(fft3d), intent(inout) :: this
    complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
    integer :: i,j,k
    include 'fftw3.f03'

    if (this%pg%nx.gt.1) then
      ! Transpose in X
      call this%fft3d_xtranspose_forward(A,this%xtrans)
      ! Inverse transform
      do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
        do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
          this%in_x=this%xtrans(:,j,k)
          call fftw_execute_dft(this%bplan_x,this%in_x,this%out_x)
          this%xtrans(:,j,k)=this%out_x/this%pg%nx
        end do
      end do
      ! Transpose back
      call this%fft3d_xtranspose_backward(this%xtrans,A)
    end if

    if (this%pg%ny.gt.1) then
      ! Transpose in Y
      call this%fft3d_ytranspose_forward(A,this%ytrans)
      ! Inverse transform
      do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
        do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
          this%in_y=this%ytrans(i,:,k)
          call fftw_execute_dft(this%bplan_y,this%in_y,this%out_y)
          this%ytrans(i,:,k)=this%out_y/this%pg%ny
        end do
      end do
      ! Transpose back
      call this%fft3d_ytranspose_backward(this%ytrans,A)
    end if

    if (this%pg%nz.gt.1) then
      ! Transpose in Z
      call this%fft3d_ztranspose_forward(A,this%ztrans)
      ! Inverse transform
      do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
        do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
          this%in_z=this%ztrans(i,j,:)
          call fftw_execute_dft(this%bplan_z,this%in_z,this%out_z)
          this%ztrans(i,j,:)=this%out_z/this%pg%nz
        end do
      end do
      ! Transpose back
      call this%fft3d_ztranspose_backward(this%ztrans,A)
    end if

  end subroutine fft3d_inverse_transform


  !> Destroy FFT object
  subroutine fft3d_destroy(this)
    use messager, only: die
    implicit none
    type(fft3d), intent(inout) :: this
    include 'fftw3.f03'

    ! Destroy our plans
    call fftw_destroy_plan(this%fplan_x); call fftw_destroy_plan(this%bplan_x);
    call fftw_destroy_plan(this%fplan_y); call fftw_destroy_plan(this%bplan_y);
    call fftw_destroy_plan(this%fplan_z); call fftw_destroy_plan(this%bplan_z);

  end subroutine fft3d_destroy


  !> Log fft3d info
  subroutine fft3d_log(this)
    use string,   only: str_long
    use messager, only: log
    implicit none
    class(fft3d), intent(in) :: this
    character(len=str_long) :: message

    if (this%pg%amRoot) then
      write(message,'("fft3d for config [",a,"]")') trim(this%pg%name)
      call log(message)
    end if

  end subroutine fft3d_log


  !> Print fft3d info to the screen
  subroutine fft3d_print(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(fft3d), intent(in) :: this

    if (this%pg%amRoot) then
      write(output_unit,'("fft3d for config [",a,"]")') trim(this%pg%name)
    end if

  end subroutine fft3d_print


  !> Short print of fft3d info to the screen
  subroutine fft3d_print_short(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(fft3d), intent(in) :: this

    if (this%pg%amRoot) then
      write(output_unit,'("fft3d for config [",a16,"]")') trim(this%pg%name)
    end if

  end subroutine fft3d_print_short


end module fft3d_class

