!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use string,        only: str_medium
   use geometry,      only: cfg1,grp1,isInGrp1
   use geometry,      only: cfg2,grp2,isInGrp2
   use ensight_class, only: ensight
   use coupler_class, only: coupler
   use datafile_class,only: datafile
   use incomp_class,  only: incomp
   implicit none
   private

   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final

   !> Ensight postprocessing
   type(ensight) :: ens1,ens2

   !> Give ourselves work arrays
   real(WP), dimension(:,:,:), allocatable :: Uid,Vid,Wid
   real(WP), dimension(:,:,:), allocatable :: Uil,Vil,Wil

   !> Datafile declaration
   type(datafile) :: dfDNS
   type(datafile) :: dfLES

   character(len=str_medium) :: dir_write
   real(WP) :: dt

   !> Solver declaration
   type(incomp),      public :: fsd
   type(incomp),      public :: fsl

   !> Coupler
   type(coupler) :: cpl

contains

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      InitDummySolv: block
         fsd=incomp(cfg=cfg1,name='DNS')
         fsl=incomp(cfg=cfg2,name='LES')
      end block InitDummySolv

      allocate(Uid  (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
      allocate(Vid  (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
      allocate(Wid  (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))

      allocate(Uil  (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
      allocate(Vil  (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
      allocate(Wil  (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))

      ! Read in DNS datafile
      readDNS: block
         character(len=str_medium) :: dir_read
         call param_read('DNS file (i)', dir_read)
         ! Read the datafile
         dfDNS=datafile(pg=cfg1,fdata=trim(adjustl(dir_read)))
         call dfDNS%pullvar(name='U',var=fsd%U)
         call dfDNS%pullvar(name='V',var=fsd%V)
         call dfDNS%pullvar(name='W',var=fsd%W)
         call dfDNS%pullvar(name='P',var=fsd%P)
      end block readDNS

      ! Init LES datafile
      initLES: block
         call param_read('LES file (o)', dir_write)
         ! Init the datafile to write
         dfLES=datafile(pg=cfg2,filename=trim(dir_write),nval=2,nvar=4)
         dfLES%valname(1)='t'
         dfLES%valname(2)='dt'
         dfLES%varname(1)='U'
         dfLES%varname(2)='V'
         dfLES%varname(3)='W'
         dfLES%varname(4)='P'
      end block initlES

      ! Both groups prepare the coupler
      if (isInGrp1.or.isInGrp2) then
         ! Create the coupler
         cpl=coupler(src_grp=grp1,dst_grp=grp2,name='LESfromDNS')
         ! Set the grids
         if (isInGrp1) call cpl%set_src(cfg1)
         if (isInGrp2) call cpl%set_dst(cfg2)
         ! Initialize the metrics
         call cpl%initialize()
      end if

      ! Group1 does its initialization work here
      if (isInGrp1) then

         ! Group1 also outputs to Ensight
         ens1=ensight(cfg=cfg1,name='DNS')
         call ens1%add_scalar('U',fsd%U)
         call ens1%add_scalar('V',fsd%V)
         call ens1%add_scalar('W',fsd%W)
         call ens1%add_scalar('P',fsd%P)
         call ens1%write_data(0.0_WP)

      end if

      ! Group2 allocates the field to receive
      if (isInGrp2) then
         fsl%U=0.0_WP
         fsl%V=0.0_WP
         fsl%W=0.0_WP
         fsl%P=0.0_WP
      end if

      if (isInGrp2) then
         ens2=ensight(cfg=cfg2,name='LES')
      end if

      ! Both groups work on the coupling
      ! NOTE: Data to be transfered is attached to the coupler object one at a time.
      !       Each interpolated variable must then be treated separately unless
      !       the coupler class source code is modified.

      coupling: block
         if (isInGrp1) call cpl%push(fsd%U)
         call cpl%transfer()
         if (isInGrp2) call cpl%pull(fsl%U)
         call ens2%add_scalar('U',fsl%U)
         call dfLES%pushvar(name='U' ,var=fsl%U   )

         if (isInGrp1) call cpl%push(fsd%V)
         call cpl%transfer()
         if (isInGrp2) call cpl%pull(fsl%V)
         call ens2%add_scalar('V',fsl%V)
         call dfLES%pushvar(name='V' ,var=fsl%V   )

         if (isInGrp1) call cpl%push(fsd%W)
         call cpl%transfer()
         if (isInGrp2) call cpl%pull(fsl%W)
         call ens2%add_scalar('W',fsl%W)
         call dfLES%pushvar(name='W' ,var=fsl%W   )

         if (isInGrp1) call cpl%push(fsd%P)
         call cpl%transfer()
         if (isInGrp2) call cpl%pull(fsl%P)
         call ens2%add_scalar('P',fsl%P)
         call dfLES%pushvar(name='P' ,var=fsl%P   )

         call dfDNS%pullval(name='dt' ,val=dt   )
         call dfLES%pushval(name='dt' ,val=dt   )
      end block coupling

      ! Group2 outputs its received data
      if (isInGrp2) then
         call ens2%add_scalar('overlap',cpl%overlap)
         call ens2%write_data(0.0_WP)

         call dfLES%write(fdata=trim(adjustl(dir_write)))
      end if

   end subroutine simulation_init


   !> Check TKE of LES and DNS
   subroutine simulation_run
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM
      use parallel,only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,ierr
      real(WP) :: myKE,KELES,KEDNS,ratio

      call dfLES%pullvar(name='U',var=fsl%U)
      call dfLES%pullvar(name='V',var=fsl%V)
      call dfLES%pullvar(name='W',var=fsl%W)

      call fsd%interp_vel(Uid,Vid,Wid)
      call fsl%interp_vel(Uil,Vil,Wil)

      myKE=0.0_WP
      do k=fsd%cfg%kmin_,fsd%cfg%kmax_
         do j=fsd%cfg%jmin_,fsd%cfg%jmax_
            do i=fsd%cfg%imin_,fsd%cfg%imax_
               myKE=myKE+0.5_WP*(Uid(i,j,k)**2+Vid(i,j,k)**2+Wid(i,j,k)**2)*fsd%cfg%vol(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(myKE,KEDNS,1,MPI_REAL_WP,MPI_SUM,fsd%cfg%comm,ierr); KEDNS = KEDNS/fsd%cfg%vol_total

      myKE=0.0_WP
      do k=fsl%cfg%kmin_,fsl%cfg%kmax_
         do j=fsl%cfg%jmin_,fsl%cfg%jmax_
            do i=fsl%cfg%imin_,fsl%cfg%imax_
               myKE=myKE+0.5_WP*(Uil(i,j,k)**2+Vil(i,j,k)**2+Wil(i,j,k)**2)*fsl%cfg%vol(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(myKE,KELES,1,MPI_REAL_WP,MPI_SUM,fsl%cfg%comm,ierr); KELES = KELES/fsl%cfg%vol_total

      ratio = KELES / KEDNS
      if (fsd%cfg%amRoot) print *, "TKE_interp / TKE_input = ", ratio

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      deallocate(Uid,Vid,Wid,Uil,Vil,Wil)
   end subroutine simulation_final

end module simulation
