!> FFT for 3D periodic uniform computational domains decomposed in at most 2
!> directions. Makes use of FFTW and in-house parallel transpose operations.
module fourier3d_class
  use precision,    only: WP
  use config_class, only: config
  use string,       only: str_short
  use, intrinsic :: iso_c_binding
  implicit none
  private


  ! Expose type/constructor/methods
  public :: fourier3d


  !> fourier3d object definition
  type :: fourier3d

    ! cfg
    class(config), pointer :: cfg

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

    procedure :: print_short=>fourier3d_print_short   !< One-line printing of solver status
    procedure :: print=>fourier3d_print               !< Long-form printing of solver status
    procedure :: log=>fourier3d_log                   !< Long-form logging of solver status
    final :: fourier3d_destroy                        !< Deconstructor

    procedure :: fourier3d_fourier_transform
    procedure :: fourier3d_inverse_transform

    procedure, private :: fourier3d_xtranspose_init
    procedure, private :: fourier3d_ytranspose_init
    procedure, private :: fourier3d_ztranspose_init

    procedure, private :: fourier3d_xtranspose_forward
    procedure, private :: fourier3d_ytranspose_forward
    procedure, private :: fourier3d_ztranspose_forward

    procedure, private :: fourier3d_xtranspose_backward
    procedure, private :: fourier3d_ytranspose_backward
    procedure, private :: fourier3d_ztranspose_backward

  end type fourier3d


  !> Declare fourier3d constructor
  interface fourier3d; procedure fourier3d_from_args; end interface fourier3d;


contains


  !> Constructor for a fourier3d object
  function fourier3d_from_args(cfg) result(self)
    use messager, only: die
    implicit none
    type(fourier3d) :: self
    class(config), target, intent(in) :: cfg
    integer :: ierr
    integer, dimension(3) :: periodicity,offset
    include 'fftw3.f03'

    ! Link the config and store the name
    self%cfg=>cfg

    ! Various checks to ensure we can use this solver
    check_solver_is_useable: block
      integer :: ndim,ndcp
      ! Periodicity and uniformity of mesh
      if (self%cfg%nx.gt.1.and..not.(self%cfg%xper.and.self%cfg%uniform_x)) &
        call die('[fourier3d constructor] Need x-direction needs to be &
        &periodic and uniform')
      if (self%cfg%ny.gt.1.and..not.(self%cfg%yper.and.self%cfg%uniform_y)) &
        call die('[fourier3d constructor] Need y-direction needs to be &
        &periodic and uniform')
      if (self%cfg%nz.gt.1.and..not.(self%cfg%zper.and.self%cfg%uniform_z)) &
        call die('[fourier3d constructor] Need z-direction needs to be &
        &periodic and uniform')
      ! Ensure that we have at least one non-decomposed direction
      ndim = count((/ self%cfg%nx, self%cfg%ny, self%cfg%nz /) .gt. 1)
      ndcp = count((/ self%cfg%npx, self%cfg%npy, self%cfg%npz /) .gt. 1)
      if (ndcp.ge.ndim) call die('[fourier3d constructor] Need at least one &
        &NON-decomposed direction')
    end block check_solver_is_useable

    ! Initialize transpose and FFTW plans in x
    if (self%cfg%nx.gt.1) then
      call self%fourier3d_xtranspose_init()
      allocate(self%in_x(self%cfg%nx),self%out_x(self%cfg%nx))
      self%fplan_x=fftw_plan_dft_1d(self%cfg%nx,self%in_x,self%out_x,FFTW_FORWARD,FFTW_MEASURE)
      self%bplan_x=fftw_plan_dft_1d(self%cfg%nx,self%in_x,self%out_x,FFTW_BACKWARD,FFTW_MEASURE)
    end if

    ! Initialize transpose and FFTW plans in y
    if (self%cfg%ny.gt.1) then
      call self%fourier3d_ytranspose_init()
      allocate(self%in_y(self%cfg%ny),self%out_y(self%cfg%ny))
      self%fplan_y=fftw_plan_dft_1d(self%cfg%ny,self%in_y,self%out_y,FFTW_FORWARD,FFTW_MEASURE)
      self%bplan_y=fftw_plan_dft_1d(self%cfg%ny,self%in_y,self%out_y,FFTW_BACKWARD,FFTW_MEASURE)
    end if

    ! Initialize transpose and FFTW plans in z
    if (self%cfg%nz.gt.1) then
      call self%fourier3d_ztranspose_init()
      allocate(self%in_z(self%cfg%nz),self%out_z(self%cfg%nz))
      self%fplan_z=fftw_plan_dft_1d(self%cfg%nz,self%in_z,self%out_z,FFTW_FORWARD,FFTW_MEASURE)
      self%bplan_z=fftw_plan_dft_1d(self%cfg%nz,self%in_z,self%out_z,FFTW_BACKWARD,FFTW_MEASURE)
    end if

    ! Find who owns the oddball
    self%oddball=all((/ self%cfg%iproc, self%cfg%jproc, self%cfg%kproc /) .eq. 1)

  end function fourier3d_from_args


  !> Initialize transpose tool in x
  subroutine fourier3d_xtranspose_init(this)
    use mpi_f08
    implicit none
    class(fourier3d), intent(inout) :: this
    integer :: ierr,ip,q,r

    ! Determine non-decomposed direction to use for transpose
    if      (this%cfg%npx.eq.1.and.this%cfg%nx.gt.1) then
      this%xdir='x'
    else if (this%cfg%npy.eq.1.and.this%cfg%ny.gt.1) then
      this%xdir='y'
    else if (this%cfg%npz.eq.1.and.this%cfg%nz.gt.1) then
      this%xdir='z'
    end if

    ! Allocate global partitions
    allocate(  this%nx_x(this%cfg%npx))
    allocate(  this%ny_x(this%cfg%npx))
    allocate(  this%nz_x(this%cfg%npx))
    allocate(this%imin_x(this%cfg%npx))
    allocate(this%imax_x(this%cfg%npx))
    allocate(this%jmin_x(this%cfg%npx))
    allocate(this%jmax_x(this%cfg%npx))
    allocate(this%kmin_x(this%cfg%npx))
    allocate(this%kmax_x(this%cfg%npx))

    ! Partition
    select case (trim(this%xdir))
    case ('x')

      ! No transpose required, use local partition
      this%nx_x=this%cfg%nx_
      this%ny_x=this%cfg%ny_
      this%nz_x=this%cfg%nz_
      this%imin_x=this%cfg%imin_
      this%imax_x=this%cfg%imax_
      this%jmin_x=this%cfg%jmin_
      this%jmax_x=this%cfg%jmax_
      this%kmin_x=this%cfg%kmin_
      this%kmax_x=this%cfg%kmax_

    case ('y')

      ! Store old local indices from each processor
      call MPI_AllGather(this%cfg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
      call MPI_AllGather(this%cfg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
      this%nx_x=this%imax_x-this%imin_x+1

      ! Partition new local indices
      do ip=1,this%cfg%npx
        q=this%cfg%ny/this%cfg%npx
        r=mod(this%cfg%ny,this%cfg%npx)
        if (ip.le.r) then
          this%ny_x(ip)  =q+1
          this%jmin_x(ip)=this%cfg%jmin+(ip-1)*(q+1)
        else
          this%ny_x(ip)  =q
          this%jmin_x(ip)=this%cfg%jmin+r*(q+1)+(ip-r-1)*q
        end if
        this%jmax_x(ip)=this%jmin_x(ip)+this%ny_x(ip)-1
      end do
      this%nz_x=this%cfg%nz_
      this%kmin_x=this%cfg%kmin_
      this%kmax_x=this%cfg%kmax_

      ! Variables for AllToAll communication
      this%sendcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%cfg%nz_
      this%recvcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%cfg%nz_
      allocate(this%sendbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npx))
      allocate(this%recvbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npx))

      ! Zero out buffers
      this%sendbuf_x=0.0_WP
      this%recvbuf_x=0.0_WP

    case ('z')

      ! Store old local indices from each processor
      call MPI_AllGather(this%cfg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
      call MPI_AllGather(this%cfg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
      this%nx_x=this%imax_x-this%imin_x+1

      ! Partition new local indices
      do ip=1,this%cfg%npx
        q=this%cfg%nz/this%cfg%npx
        r=mod(this%cfg%nz,this%cfg%npx)
        if (ip.le.r) then
          this%nz_x(ip)  =q+1
          this%kmin_x(ip)=this%cfg%kmin+(ip-1)*(q+1)
        else
          this%nz_x(ip)  =q
          this%kmin_x(ip)=this%cfg%kmin+r*(q+1)+(ip-r-1)*q
        end if
        this%kmax_x(ip)=this%kmin_x(ip)+this%nz_x(ip)-1
      end do
      this%ny_x=this%cfg%ny_
      this%jmin_x=this%cfg%jmin_
      this%jmax_x=this%cfg%jmax_

      ! Variables for AllToAll communication
      this%sendcount_x=maxval(this%nx_x)*this%cfg%ny_*maxval(this%nz_x)
      this%recvcount_x=maxval(this%nx_x)*this%cfg%ny_*maxval(this%nz_x)
      allocate(this%sendbuf_x(maxval(this%nx_x),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_x),this%cfg%npx))
      allocate(this%recvbuf_x(maxval(this%nx_x),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_x),this%cfg%npx))

      ! Zero out buffers
      this%sendbuf_x=0.0_WP
      this%recvbuf_x=0.0_WP

    end select

    ! Allocate storage
    allocate(this%xtrans(this%cfg%imin:this%cfg%imax,this%jmin_x(this%cfg%iproc):this%jmax_x(this%cfg%iproc),this%kmin_x(this%cfg%iproc):this%kmax_x(this%cfg%iproc)))

  end subroutine fourier3d_xtranspose_init


  !> Initialize transpose tool in y
  subroutine fourier3d_ytranspose_init(this)
    use mpi_f08
    implicit none
    class(fourier3d), intent(inout) :: this
    integer :: ierr,jp,q,r

    ! Determine non-decomposed direction to use for transpose
    if      (this%cfg%npy.eq.1.and.this%cfg%ny.gt.1) then
      this%ydir='y'
    else if (this%cfg%npz.eq.1.and.this%cfg%nz.gt.1) then
      this%ydir='z'
    else if (this%cfg%npx.eq.1.and.this%cfg%nx.gt.1) then
      this%ydir='x'
    end if

    ! Allocate global partitions
    allocate(  this%nx_y(this%cfg%npy))
    allocate(  this%ny_y(this%cfg%npy))
    allocate(  this%nz_y(this%cfg%npy))
    allocate(this%imin_y(this%cfg%npy))
    allocate(this%imax_y(this%cfg%npy))
    allocate(this%jmin_y(this%cfg%npy))
    allocate(this%jmax_y(this%cfg%npy))
    allocate(this%kmin_y(this%cfg%npy))
    allocate(this%kmax_y(this%cfg%npy))

    ! Partition
    select case (trim(this%ydir))
    case ('x')

      ! Store old local indices from each processor
      call MPI_AllGather(this%cfg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
      call MPI_AllGather(this%cfg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
      this%ny_y=this%jmax_y-this%jmin_y+1

      ! Partition new local indices
      do jp=1,this%cfg%npy
        q=this%cfg%nx/this%cfg%npy
        r=mod(this%cfg%nx,this%cfg%npy)
        if (jp.le.r) then
          this%nx_y(jp)  =q+1
          this%imin_y(jp)=this%cfg%imin+(jp-1)*(q+1)
        else
          this%nx_y(jp)  =q
          this%imin_y(jp)=this%cfg%imin+r*(q+1)+(jp-r-1)*q
        end if
        this%imax_y(jp)=this%imin_y(jp)+this%nx_y(jp)-1
      end do
      this%nz_y=this%cfg%nz_
      this%kmin_y=this%cfg%kmin_
      this%kmax_y=this%cfg%kmax_

      ! Variables for AllToAll communication
      this%sendcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%cfg%nz_
      this%recvcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%cfg%nz_
      allocate(this%sendbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npy))
      allocate(this%recvbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npy))

      ! Zero out buffers
      this%sendbuf_y=0.0_WP
      this%recvbuf_y=0.0_WP

    case ('y')

      ! No transpose required, use local partition
      this%nx_y=this%cfg%nx_
      this%ny_y=this%cfg%ny_
      this%nz_y=this%cfg%nz_
      this%imin_y=this%cfg%imin_
      this%imax_y=this%cfg%imax_
      this%jmin_y=this%cfg%jmin_
      this%jmax_y=this%cfg%jmax_
      this%kmin_y=this%cfg%kmin_
      this%kmax_y=this%cfg%kmax_

    case ('z')

      ! Store old local indices from each processor
      call MPI_AllGather(this%cfg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
      call MPI_AllGather(this%cfg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
      this%ny_y=this%jmax_y-this%jmin_y+1

      ! Partition new local indices
      do jp=1,this%cfg%npy
        q=this%cfg%nz/this%cfg%npy
        r=mod(this%cfg%nz,this%cfg%npy)
        if (jp.le.r) then
          this%nz_y(jp)  =q+1
          this%kmin_y(jp)=this%cfg%kmin+(jp-1)*(q+1)
        else
          this%nz_y(jp)  =q
          this%kmin_y(jp)=this%cfg%kmin+r*(q+1)+(jp-r-1)*q
        end if
        this%kmax_y(jp)=this%kmin_y(jp)+this%nz_y(jp)-1
      end do
      this%nx_y=this%cfg%nx_
      this%imin_y=this%cfg%imin_
      this%imax_y=this%cfg%imax_

      ! Variables for AllToAll communication
      this%sendcount_y=this%cfg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
      this%recvcount_y=this%cfg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
      allocate(this%sendbuf_y(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%cfg%npy))
      allocate(this%recvbuf_y(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%cfg%npy))

      ! Zero out buffers
      this%sendbuf_y=0.0_WP
      this%recvbuf_y=0.0_WP

    end select

    ! Allocate storage
    allocate(this%ytrans(this%imin_y(this%cfg%jproc):this%imax_y(this%cfg%jproc),this%cfg%jmin:this%cfg%jmax,this%kmin_y(this%cfg%jproc):this%kmax_y(this%cfg%jproc)))

  end subroutine fourier3d_ytranspose_init


  !> Initialize transpose tool in z
  subroutine fourier3d_ztranspose_init(this)
    use mpi_f08
    implicit none
    class(fourier3d), intent(inout) :: this
    integer :: ierr,kp,q,r

    ! Determine non-decomposed direction to use for transpose
    if      (this%cfg%npz.eq.1.and.this%cfg%nz.gt.1) then
      this%zdir='z'
    else if (this%cfg%npx.eq.1.and.this%cfg%nx.gt.1) then
      this%zdir='x'
    else if (this%cfg%npy.eq.1.and.this%cfg%ny.gt.1) then
      this%zdir='y'
    end if

    ! Allocate global partitions
    allocate(  this%nx_z(this%cfg%npz))
    allocate(  this%ny_z(this%cfg%npz))
    allocate(  this%nz_z(this%cfg%npz))
    allocate(this%imin_z(this%cfg%npz))
    allocate(this%imax_z(this%cfg%npz))
    allocate(this%jmin_z(this%cfg%npz))
    allocate(this%jmax_z(this%cfg%npz))
    allocate(this%kmin_z(this%cfg%npz))
    allocate(this%kmax_z(this%cfg%npz))

    ! Partition
    select case (trim(this%zdir))
    case ('x')

      ! Store old local indices from each processor
      call MPI_AllGather(this%cfg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
      call MPI_AllGather(this%cfg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
      this%nz_z=this%kmax_z-this%kmin_z+1

      ! Partition new local indices
      do kp=1,this%cfg%npz
        q=this%cfg%nx/this%cfg%npz
        r=mod(this%cfg%nx,this%cfg%npz)
        if (kp.le.r) then
          this%nx_z(kp)  =q+1
          this%imin_z(kp)=this%cfg%imin+(kp-1)*(q+1)
        else
          this%nx_z(kp)  =q
          this%imin_z(kp)=this%cfg%imin+r*(q+1)+(kp-r-1)*q
        end if
        this%imax_z(kp)=this%imin_z(kp)+this%nx_z(kp)-1
      end do
      this%ny_z=this%cfg%ny_
      this%jmin_z=this%cfg%jmin_
      this%jmax_z=this%cfg%jmax_

      ! Variables for AllToAll communication
      this%sendcount_z=maxval(this%nx_z)*this%cfg%ny_*maxval(this%nz_z)
      this%recvcount_z=maxval(this%nx_z)*this%cfg%ny_*maxval(this%nz_z)
      allocate(this%sendbuf_z(maxval(this%nx_z),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_z),this%cfg%npz))
      allocate(this%recvbuf_z(maxval(this%nx_z),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_z),this%cfg%npz))

      ! Zero out buffers
      this%sendbuf_z=0.0_WP
      this%recvbuf_z=0.0_WP

    case ('y')

      ! Store old local indices from each processor
      call MPI_AllGather(this%cfg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
      call MPI_AllGather(this%cfg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
      this%nz_z=this%kmax_z-this%kmin_z+1

      ! Partition new local indices
      do kp=1,this%cfg%npz
        q=this%cfg%ny/this%cfg%npz
        r=mod(this%cfg%ny,this%cfg%npz)
        if (kp.le.r) then
          this%ny_z(kp)  =q+1
          this%jmin_z(kp)=this%cfg%jmin+(kp-1)*(q+1)
        else
          this%ny_z(kp)  =q
          this%jmin_z(kp)=this%cfg%jmin+r*(q+1)+(kp-r-1)*q
        end if
        this%jmax_z(kp)=this%jmin_z(kp)+this%ny_z(kp)-1
      end do
      this%nx_z=this%cfg%nx_
      this%imin_z=this%cfg%imin_
      this%imax_z=this%cfg%imax_

      ! Variables for AllToAll communication
      this%sendcount_z=this%cfg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
      this%recvcount_z=this%cfg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
      allocate(this%sendbuf_z(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%cfg%npz))
      allocate(this%recvbuf_z(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%cfg%npz))

      ! Zero out buffers
      this%sendbuf_z=0.0_WP
      this%recvbuf_z=0.0_WP

    case ('z')

      ! No transpose required, use local partition
      this%nx_z=this%cfg%nx_
      this%ny_z=this%cfg%ny_
      this%nz_z=this%cfg%nz_
      this%imin_z=this%cfg%imin_
      this%imax_z=this%cfg%imax_
      this%jmin_z=this%cfg%jmin_
      this%jmax_z=this%cfg%jmax_
      this%kmin_z=this%cfg%kmin_
      this%kmax_z=this%cfg%kmax_

    end select

    ! Allocate storage
    allocate(this%ztrans(this%imin_z(this%cfg%kproc):this%imax_z(this%cfg%kproc),this%jmin_z(this%cfg%kproc):this%jmax_z(this%cfg%kproc),this%cfg%kmin:this%cfg%kmax))

  end subroutine fourier3d_ztranspose_init


  !> Perform forward transpose in x
  subroutine fourier3d_xtranspose_forward(this,A,At)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(in) :: A
    complex(WP), dimension(this%cfg%imin :,this%jmin_x(this%cfg%iproc):,this%kmin_x(this%cfg%iproc):), intent(out) :: At
    integer :: i,j,k,ip,ii,jj,kk,ierr

    select case (trim(this%xdir))
    case ('x')
      ! No transpose required
      At=A
    case ('y')
      ! Transpose x=>y
      do ip=1,this%cfg%npx
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_x(ip),this%jmax_x(ip)
            do i=this%cfg%imin_,this%cfg%imax_
              jj=j-this%jmin_x(ip)+1
              ii=i-this%cfg%imin_+1
              this%sendbuf_x(ii,jj,k,ip)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%cfg%xcomm,ierr)
      do ip=1,this%cfg%npx
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_x(this%cfg%iproc),this%jmax_x(this%cfg%iproc)
            do i=this%imin_x(ip),this%imax_x(ip)
              jj=j-this%jmin_x(this%cfg%iproc)+1
              ii=i-this%imin_x(ip)+1
              At(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
            end do
          end do
        end do
      end do
    case ('z')
      ! Transpose x=>z
      do ip=1,this%cfg%npx
        do k=this%kmin_x(ip),this%kmax_x(ip)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
              kk=k-this%kmin_x(ip)+1
              ii=i-this%cfg%imin_+1
              this%sendbuf_x(ii,j,kk,ip)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%cfg%xcomm,ierr)
      do ip=1,this%cfg%npx
        do k=this%kmin_x(this%cfg%iproc),this%kmax_x(this%cfg%iproc)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_x(ip),this%imax_x(ip)
              kk=k-this%kmin_x(this%cfg%iproc)+1
              ii=i-this%imin_x(ip)+1
              At(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
            end do
          end do
        end do
      end do
    end select

  end subroutine fourier3d_xtranspose_forward


  !> Perform forward transpose in y
  subroutine fourier3d_ytranspose_forward(this,A,At)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(in) :: A
    complex(WP), dimension(this%imin_y(this%cfg%jproc):,this%cfg%jmin:,this%kmin_y(this%cfg%jproc):), intent(out) :: At
    integer :: i,j,k,jp,ii,jj,kk,ierr

    select case (trim(this%ydir))
    case ('x')
      ! Transpose y=>x
      do jp=1,this%cfg%npy
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_y(jp),this%imax_y(jp)
              ii=i-this%imin_y(jp)+1
              jj=j-this%cfg%jmin_+1
              this%sendbuf_y(ii,jj,k,jp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%cfg%ycomm,ierr)
      do jp=1,this%cfg%npy
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%imin_y(this%cfg%jproc),this%imax_y(this%cfg%jproc)
              ii=i-this%imin_y(this%cfg%jproc)+1
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
      do jp=1,this%cfg%npy
        do k=this%kmin_y(jp),this%kmax_y(jp)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
              kk=k-this%kmin_y(jp)+1
              jj=j-this%cfg%jmin_+1
              this%sendbuf_y(i,jj,kk,jp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%cfg%ycomm,ierr)
      do jp=1,this%cfg%npy
        do k=this%kmin_y(this%cfg%jproc),this%kmax_y(this%cfg%jproc)
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%cfg%imin_,this%cfg%imax_
              kk=k-this%kmin_y(this%cfg%jproc)+1
              jj=j-this%jmin_y(jp)+1
              At(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
            end do
          end do
        end do
      end do
    end select

  end subroutine fourier3d_ytranspose_forward


  !> Perform forward transpose in z
  subroutine fourier3d_ztranspose_forward(this,A,At)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(in) :: A
    complex(WP), dimension(this%imin_z(this%cfg%kproc):,this%jmin_z(this%cfg%kproc):,this%cfg%kmin:), intent(out) :: At
    integer :: i,j,k,kp,ii,jj,kk,ierr

    select case (trim(this%zdir))
    case ('x')
      ! Transpose z=>x
      do kp=1,this%cfg%npz
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_z(kp),this%imax_z(kp)
              ii=i-this%imin_z(kp)+1
              kk=k-this%cfg%kmin_+1
              this%sendbuf_z(ii,j,kk,kp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%cfg%zcomm,ierr)
      do kp=1,this%cfg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_z(this%cfg%kproc),this%imax_z(this%cfg%kproc)
              ii=i-this%imin_z(this%cfg%kproc)+1
              kk=k-this%kmin_z(kp)+1
              At(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
            end do
          end do
        end do
      end do
    case ('y')
      ! Transpose z=>y
      do kp=1,this%cfg%npz
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_z(kp),this%jmax_z(kp)
            do i=this%cfg%imin_,this%cfg%imax_
              jj=j-this%jmin_z(kp)+1
              kk=k-this%cfg%kmin_+1
              this%sendbuf_z(i,jj,kk,kp)=A(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%cfg%zcomm,ierr)
      do kp=1,this%cfg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%jmin_z(this%cfg%kproc),this%jmax_z(this%cfg%kproc)
            do i=this%cfg%imin_,this%cfg%imax_
              jj=j-this%jmin_z(this%cfg%kproc)+1
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

  end subroutine fourier3d_ztranspose_forward


  !> Perform backward transpose in x
  subroutine fourier3d_xtranspose_backward(this,At,A)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%cfg%imin :,this%jmin_x(this%cfg%iproc):,this%kmin_x(this%cfg%iproc):), intent(in) :: At
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(out) :: A
    integer :: i,j,k,ip,ii,jj,kk,ierr

    select case (trim(this%xdir))
    case ('x')
      ! No transpose required
      A=At
    case ('y')
      ! Transpose y=>x
      do ip=1,this%cfg%npx
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_x(this%cfg%iproc),this%jmax_x(this%cfg%iproc)
            do i=this%imin_x(ip),this%imax_x(ip)
              jj=j-this%jmin_x(this%cfg%iproc)+1
              ii=i-this%imin_x(ip)+1
              this%sendbuf_x(ii,jj,k,ip)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%cfg%xcomm,ierr)
      do ip=1,this%cfg%npx
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_x(ip),this%jmax_x(ip)
            do i=this%imin_x(this%cfg%iproc),this%imax_x(this%cfg%iproc)
              jj=j-this%jmin_x(ip)+1
              ii=i-this%imin_x(this%cfg%iproc)+1
              A(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
            end do
          end do
        end do
      end do
    case ('z')
      ! Transpose z=>x
      do ip=1,this%cfg%npx
        do k=this%kmin_x(this%cfg%iproc),this%kmax_x(this%cfg%iproc)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_x(ip),this%imax_x(ip)
              kk=k-this%kmin_x(this%cfg%iproc)+1
              ii=i-this%imin_x(ip)+1
              this%sendbuf_x(ii,j,kk,ip)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_DOUBLE_COMPLEX,this%recvbuf_x,this%recvcount_x,MPI_DOUBLE_COMPLEX,this%cfg%xcomm,ierr)
      do ip=1,this%cfg%npx
        do k=this%kmin_x(ip),this%kmax_x(ip)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_x(this%cfg%iproc),this%imax_x(this%cfg%iproc)
              kk=k-this%kmin_x(ip)+1
              ii=i-this%imin_x(this%cfg%iproc)+1
              A(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
            end do
          end do
        end do
      end do
    end select

  end subroutine fourier3d_xtranspose_backward


  !> Perform backward transpose in y
  subroutine fourier3d_ytranspose_backward(this,At,A)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%imin_y(this%cfg%jproc):,this%cfg%jmin:,this%kmin_y(this%cfg%jproc):), intent(in) :: At
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(out) :: A
    integer :: i,j,k,jp,ii,jj,kk,ierr

    select case (trim(this%ydir))
    case ('x')
      ! Transpose x=>y
      do jp=1,this%cfg%npy
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%imin_y(this%cfg%jproc),this%imax_y(this%cfg%jproc)
              ii=i-this%imin_y(this%cfg%jproc)+1
              jj=j-this%jmin_y(jp)+1
              this%sendbuf_y(ii,jj,k,jp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%cfg%ycomm,ierr)
      do jp=1,this%cfg%npy
        do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%jmin_y(this%cfg%jproc),this%jmax_y(this%cfg%jproc)
            do i=this%imin_y(jp),this%imax_y(jp)
              ii=i-this%imin_y(jp)+1
              jj=j-this%jmin_y(this%cfg%jproc)+1
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
      do jp=1,this%cfg%npy
        do k=this%kmin_y(this%cfg%jproc),this%kmax_y(this%cfg%jproc)
          do j=this%jmin_y(jp),this%jmax_y(jp)
            do i=this%cfg%imin_,this%cfg%imax_
              kk=k-this%kmin_y(this%cfg%jproc)+1
              jj=j-this%jmin_y(jp)+1
              this%sendbuf_y(i,jj,kk,jp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_DOUBLE_COMPLEX,this%recvbuf_y,this%recvcount_y,MPI_DOUBLE_COMPLEX,this%cfg%ycomm,ierr)
      do jp=1,this%cfg%npy
        do k=this%kmin_y(jp),this%kmax_y(jp)
          do j=this%jmin_y(this%cfg%jproc),this%jmax_y(this%cfg%jproc)
            do i=this%cfg%imin_,this%cfg%imax_
              kk=k-this%kmin_y(jp)+1
              jj=j-this%jmin_y(this%cfg%jproc)+1
              A(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
            end do
          end do
        end do
      end do
    end select

  end subroutine fourier3d_ytranspose_backward


  !> Perform backward transpose in z
  subroutine fourier3d_ztranspose_backward(this,At,A)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%imin_z(this%cfg%kproc):,this%jmin_z(this%cfg%kproc):,this%cfg%kmin:), intent(in) :: At
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(out) :: A
    integer :: i,j,k,kp,ii,jj,kk,ierr

    select case (trim(this%zdir))
    case ('x')
      ! Transpose x=>z
      do kp=1,this%cfg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_z(this%cfg%kproc),this%imax_z(this%cfg%kproc)
              ii=i-this%imin_z(this%cfg%kproc)+1
              kk=k-this%kmin_z(kp)+1
              this%sendbuf_z(ii,j,kk,kp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%cfg%zcomm,ierr)
      do kp=1,this%cfg%npz
        do k=this%kmin_z(this%cfg%kproc),this%kmax_z(this%cfg%kproc)
          do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%imin_z(kp),this%imax_z(kp)
              ii=i-this%imin_z(kp)+1
              kk=k-this%kmin_z(this%cfg%kproc)+1
              A(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
            end do
          end do
        end do
      end do
    case ('y')
      ! Transpose y=>z
      do kp=1,this%cfg%npz
        do k=this%kmin_z(kp),this%kmax_z(kp)
          do j=this%jmin_z(this%cfg%kproc),this%jmax_z(this%cfg%kproc)
            do i=this%cfg%imin_,this%cfg%imax_
              jj=j-this%jmin_z(this%cfg%kproc)+1
              kk=k-this%kmin_z(kp)+1
              this%sendbuf_z(i,jj,kk,kp)=At(i,j,k)
            end do
          end do
        end do
      end do
      call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_DOUBLE_COMPLEX,this%recvbuf_z,this%recvcount_z,MPI_DOUBLE_COMPLEX,this%cfg%zcomm,ierr)
      do kp=1,this%cfg%npz
        do k=this%kmin_z(this%cfg%kproc),this%kmax_z(this%cfg%kproc)
          do j=this%jmin_z(kp),this%jmax_z(kp)
            do i=this%cfg%imin_,this%cfg%imax_
              jj=j-this%jmin_z(kp)+1
              kk=k-this%kmin_z(this%cfg%kproc)+1
              A(i,j,k)=this%recvbuf_z(i,jj,kk,kp)
            end do
          end do
        end do
      end do
    case ('z')
      ! No transpose required
      A=At
    end select

  end subroutine fourier3d_ztranspose_backward


  !> Transpose A and perform Fourier transform
  subroutine fourier3d_fourier_transform(this,A)
    use messager, only: die
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
    integer :: i,j,k
    include 'fftw3.f03'

    if (this%cfg%nx.gt.1) then
      ! Transpose in X
      call this%fourier3d_xtranspose_forward(A,this%xtrans)
      ! Forward transform - X
      do k=this%kmin_x(this%cfg%iproc),this%kmax_x(this%cfg%iproc)
        do j=this%jmin_x(this%cfg%iproc),this%jmax_x(this%cfg%iproc)
          this%in_x=this%xtrans(:,j,k)
          call fftw_execute_dft(this%fplan_x,this%in_x,this%out_x)
          this%xtrans(:,j,k)=this%out_x
        end do
      end do
      ! Transpose back
      call this%fourier3d_xtranspose_backward(this%xtrans,A)
    end if

    if (this%cfg%ny.gt.1) then
      ! Transpose in Y
      call this%fourier3d_ytranspose_forward(A,this%ytrans)
      ! Forward transform - Y
      do k=this%kmin_y(this%cfg%jproc),this%kmax_y(this%cfg%jproc)
        do i=this%imin_y(this%cfg%jproc),this%imax_y(this%cfg%jproc)
          this%in_y=this%ytrans(i,:,k)
          call fftw_execute_dft(this%fplan_y,this%in_y,this%out_y)
          this%ytrans(i,:,k)=this%out_y
        end do
      end do
      ! Transpose back
      call this%fourier3d_ytranspose_backward(this%ytrans,A)
    end if

    if (this%cfg%nz.gt.1) then
      ! Transpose in Z
      call this%fourier3d_ztranspose_forward(A,this%ztrans)
      ! Forward transform - Z
      do j=this%jmin_z(this%cfg%kproc),this%jmax_z(this%cfg%kproc)
        do i=this%imin_z(this%cfg%kproc),this%imax_z(this%cfg%kproc)
          this%in_z=this%ztrans(i,j,:)
          call fftw_execute_dft(this%fplan_z,this%in_z,this%out_z)
          this%ztrans(i,j,:)=this%out_z
        end do
      end do
      ! Transpose back
      call this%fourier3d_ztranspose_backward(this%ztrans,A)
    end if

    ! Oddball
    if (this%oddball) A(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_)=0.0_WP

  end subroutine fourier3d_fourier_transform


  !> FFT -> real and transpose back
  subroutine fourier3d_inverse_transform(this,A)
    use messager, only: die
    implicit none
    class(fourier3d), intent(inout) :: this
    complex(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
    integer :: i,j,k
    include 'fftw3.f03'

    if (this%cfg%nx.gt.1) then
      ! Transpose in X
      call this%fourier3d_xtranspose_forward(A,this%xtrans)
      ! Inverse transform
      do k=this%kmin_x(this%cfg%iproc),this%kmax_x(this%cfg%iproc)
        do j=this%jmin_x(this%cfg%iproc),this%jmax_x(this%cfg%iproc)
          this%in_x=this%xtrans(:,j,k)
          call fftw_execute_dft(this%bplan_x,this%in_x,this%out_x)
          this%xtrans(:,j,k)=this%out_x/this%cfg%nx
        end do
      end do
      ! Transpose back
      call this%fourier3d_xtranspose_backward(this%xtrans,A)
    end if

    if (this%cfg%ny.gt.1) then
      ! Transpose in Y
      call this%fourier3d_ytranspose_forward(A,this%ytrans)
      ! Inverse transform
      do k=this%kmin_y(this%cfg%jproc),this%kmax_y(this%cfg%jproc)
        do i=this%imin_y(this%cfg%jproc),this%imax_y(this%cfg%jproc)
          this%in_y=this%ytrans(i,:,k)
          call fftw_execute_dft(this%bplan_y,this%in_y,this%out_y)
          this%ytrans(i,:,k)=this%out_y/this%cfg%ny
        end do
      end do
      ! Transpose back
      call this%fourier3d_ytranspose_backward(this%ytrans,A)
    end if

    if (this%cfg%nz.gt.1) then
      ! Transpose in Z
      call this%fourier3d_ztranspose_forward(A,this%ztrans)
      ! Inverse transform
      do j=this%jmin_z(this%cfg%kproc),this%jmax_z(this%cfg%kproc)
        do i=this%imin_z(this%cfg%kproc),this%imax_z(this%cfg%kproc)
          this%in_z=this%ztrans(i,j,:)
          call fftw_execute_dft(this%bplan_z,this%in_z,this%out_z)
          this%ztrans(i,j,:)=this%out_z/this%cfg%nz
        end do
      end do
      ! Transpose back
      call this%fourier3d_ztranspose_backward(this%ztrans,A)
    end if

  end subroutine fourier3d_inverse_transform


  !> Destroy FFT object
  subroutine fourier3d_destroy(this)
    use messager, only: die
    implicit none
    type(fourier3d), intent(inout) :: this
    include 'fftw3.f03'

    ! Destroy our plans
    call fftw_destroy_plan(this%fplan_x); call fftw_destroy_plan(this%bplan_x);
    call fftw_destroy_plan(this%fplan_y); call fftw_destroy_plan(this%bplan_y);
    call fftw_destroy_plan(this%fplan_z); call fftw_destroy_plan(this%bplan_z);

  end subroutine fourier3d_destroy


  !> Log fourier3d info
  subroutine fourier3d_log(this)
    use string,   only: str_long
    use messager, only: log
    implicit none
    class(fourier3d), intent(in) :: this
    character(len=str_long) :: message

    if (this%cfg%amRoot) then
      write(message,'("fourier3d for config [",a,"]")') trim(this%cfg%name)
      call log(message)
    end if

  end subroutine fourier3d_log


  !> Print fourier3d info to the screen
  subroutine fourier3d_print(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(fourier3d), intent(in) :: this

    if (this%cfg%amRoot) then
      write(output_unit,'("fourier3d for config [",a,"]")') trim(this%cfg%name)
    end if

  end subroutine fourier3d_print


  !> Short print of fourier3d info to the screen
  subroutine fourier3d_print_short(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(fourier3d), intent(in) :: this

    if (this%cfg%amRoot) then
      write(output_unit,'("fourier3d for config [",a16,"]")') trim(this%cfg%name)
    end if

  end subroutine fourier3d_print_short


end module fourier3d_class

