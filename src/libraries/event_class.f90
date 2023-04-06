module event_class
   use precision,         only: WP
   use string,            only: str_medium
   use timetracker_class, only: timetracker
   implicit none
   private

   ! Expose type/constructor/methods
   public :: event, periodic_event, threshold_event

   !> Safety coefficient to avoid barely missing an occurence due to round off
   real(WP), parameter :: CSAFE=100.0_WP

   !> Event objects
   type, abstract :: event
      character(len=str_medium)   :: name='UNNAMED_EVENT' !< Name for event
      class(timetracker), pointer :: time                 !< Timetracker for event
   contains
      procedure(occur_ftype), deferred :: occurs          !< Check if event is occuring
   end type event
   interface
     logical function occur_ftype(this)
       import :: event
       implicit none
       class(event), intent(inout) :: this
     end function occur_ftype
   end interface

   type, extends(event) :: periodic_event
      integer  :: nper                        !< Period in elapsed number of time steps
      real(WP) :: tper                        !< Period in elapsed time
   contains
      procedure :: occurs => periodic_occurs  !< Check if event is occuring
   end type periodic_event
   interface periodic_event
      procedure periodic_event_constructor
   end interface periodic_event

   type, extends(event) :: threshold_event
     integer  :: nnext
     real(WP) :: tnext
     logical  :: occurred
   contains
     procedure :: occurs => threshold_occurs  !< Check if event is occuring
   end type threshold_event
   interface threshold_event
      procedure threshold_event_constructor
   end interface threshold_event


contains


   !> Constructor for periodic event
   function periodic_event_constructor(time,name) result(self)
      implicit none
      type(periodic_event) :: self
      class(timetracker), target, intent(in) :: time
      character(len=*), optional :: name

      ! Point to timetracker object
      self%time=>time
      ! Set the event name
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to 0 periods for event
      self%nper=-1
      self%tper=-1.0_WP

   end function periodic_event_constructor


   !> Occurence check for periodic event
   logical function periodic_occurs(this) result(occurs)
      implicit none
      class(periodic_event), intent(inout) :: this

      ! Assume not occuring
      occurs = .false.
      ! Go through standard occurence tests
      occurs = this%nper.gt.0 .and. mod(this%time%n,this%nper).eq.0
      occurs = occurs .or. (this%tper.gt.0.0_WP .and.                         &
        mod(this%time%t+CSAFE*epsilon(this%time%t),this%tper).lt.this%time%dt)

   end function periodic_occurs


   !> Constructor for threshold event
   function threshold_event_constructor(time,name) result(self)
      implicit none
      type(threshold_event) :: self
      class(timetracker), target, intent(in) :: time
      character(len=*), optional :: name

      ! Point to timetracker object
      self%time=>time
      ! Set the event name
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to never occuring
      self%nnext = -1
      self%tnext = -1.0_WP
      self%occurred = .false.

   end function threshold_event_constructor


   !> Occurence check for threshold event
   logical function threshold_occurs(this) result(occurs)
      implicit none
      class(threshold_event), intent(inout) :: this

      if (this%occurred) then
        occurs = .false.
      else
        occurs = this%nnext .gt. 0 .and. this%nnext .le. this%time%n
        occurs = (this%tnext .gt. 0.0_WP .and. this%tnext .le. this%time%t) .or. occurs
        this%occurred = this%occurred .or. occurs
      end if

   end function threshold_occurs


end module event_class
