!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,       only: WP
   use hit_class,       only: hit
   use crossflow_class, only: crossflow
   implicit none
   private
   
   !> HIT simulation
   type(hit) :: turb
   logical :: isInHITGrp
   
   !> Crossflow simulation
   type(crossflow) :: xflow
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      type(MPI_Group) :: hit_group
      
      ! Initialize crossflow simulation
      call xflow%init()
      
      ! Create an MPI group using leftmost processors only
      create_hit_group: block
         use parallel, only: group
         use mpi_f08,  only: MPI_Group_incl
         integer, dimension(:), allocatable :: ranks
         integer, dimension(3) :: coord
         integer :: ngrp,ierr,ny,nz
         ngrp=xflow%cfg%npy*xflow%cfg%npz
         allocate(ranks(ngrp))
         ngrp=0
         do nz=1,xflow%cfg%npz
            do ny=1,xflow%cfg%npy
               ngrp=ngrp+1
               coord=[0,ny-1,nz-1]
               call MPI_CART_RANK(xflow%cfg%comm,coord,ranks(ngrp),ierr)
            end do
         end do
         call MPI_Group_incl(group,ngrp,ranks,hit_group,ierr)
         if (xflow%cfg%iproc.eq.1) then
            isInHITGrp=.true.
         else
            isInHITGrp=.false.
         end if
       end block create_hit_group
      
      ! Prepare HIT simulation
      if (isInHITGrp) then
         prepare_hit: block
            real(WP) :: dt
            ! Initialize HIT
            call turb%init(group=hit_group)
            ! Run HIT until t/tau_eddy=20
            dt=0.15_WP*turb%cfg%min_meshsize/turb%Urms_tgt !< Estimate maximum stable dt
            do while (turb%time%t.lt.2.0_WP*turb%tau_tgt)
               call turb%step(dt)
            end do
            ! Reset time
            turb%time%t=0.0_WP
         end block prepare_hit
      end if      
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
     implicit none
      
      ! Crossflow drives overall time integration
      do while (.not.xflow%time%done())
         
         ! Advance crossflow simulation
         call xflow%step()
         
         ! Advance HIT simulation and transfer velocity info
         if (isInHITGrp) then
            ! Advance HIT
            call turb%step(xflow%time%dt)
            ! Transfer turbulent velocity from hit to xflow
            apply_boundary_condition: block
               use incomp_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k,ihit
               real(WP) :: rescaling
               rescaling=1.0_WP!turb%ti/turb%Urms_tgt
               call xflow%fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  ihit=i-xflow%fs%cfg%imin+turb%fs%cfg%imax+1
                  xflow%fs%U(i  ,j,k)=xflow%Ubulk+turb%fs%U(ihit  ,j,k)*rescaling
                  xflow%fs%V(i-1,j,k)=       turb%fs%V(ihit-1,j,k)*rescaling
                  xflow%fs%W(i-1,j,k)=       turb%fs%W(ihit-1,j,k)*rescaling
               end do
            end block apply_boundary_condition
         end if
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize crossflow simulation
      call xflow%final()
      
      ! Finalize HIT simulation
      if (isInHITGrp) call turb%final()
      
   end subroutine simulation_final
   
   
end module simulation
