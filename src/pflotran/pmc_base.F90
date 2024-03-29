! Process Model Coupler Base class
module PMC_Base_class

  use PM_Base_class
  use Timestepper_Base_class
  use Option_module
  use Waypoint_module
  use PM_Base_Pointer_module
  use Output_module, only : Output
  use Simulation_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
  
  private
  
  ! process model coupler type
  type, public :: pmc_base_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: stage
    PetscBool :: is_master
    type(option_type), pointer :: option
    class(timestepper_base_type), pointer :: timestepper
    class(pm_base_type), pointer :: pms
    type(waypoint_list_type), pointer :: waypoint_list
    class(pmc_base_type), pointer :: below
    class(pmc_base_type), pointer :: next
    type(pm_base_pointer_type), pointer :: pm_ptr
    type(simulation_aux_type),pointer :: sim_aux
    procedure(Output), nopass, pointer :: Output
  contains
    procedure, public :: Init => PMCBaseInit
    procedure, public :: InitializeRun
    procedure, public :: CastToBase => PMCCastToBase
    procedure, public :: SetTimestepper => PMCBaseSetTimestepper
    procedure, public :: RunToTime => PMCBaseRunToTime
    procedure, public :: Checkpoint => PMCBaseCheckpoint
    procedure, public :: Restart => PMCBaseRestart
    procedure, public :: FinalizeRun
    procedure, public :: OutputLocal
    procedure, public :: UpdateSolution => PMCBaseUpdateSolution
    procedure, public :: Destroy => PMCBaseDestroy
    procedure, public :: AccumulateAuxData
    procedure, public :: GetAuxData
    procedure, public :: SetAuxData
  end type pmc_base_type
  
  abstract interface
    subroutine Synchronize(pmc)
      import pmc_base_type
      implicit none
        class(pmc_base_type) :: pmc
    end subroutine Synchronize
  end interface

  ! For checkpointing  
  type, public :: pmc_base_header_type
    integer*8 :: plot_number      ! in the checkpoint file format
    integer*8 :: times_per_h5_file! in the checkpoint file format
  end type pmc_base_header_type
  PetscSizeT, parameter, private :: bagsize = 16

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: pmc_base_header_type
      implicit none
#include "finclude/petscbag.h"      
      PetscBag :: bag
      class(pmc_base_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData   
    
  public :: PMCBaseCreate, &
            PMCBaseInit, &
            PMCBaseStrip, &
            SetOutputFlags
  
contains

! ************************************************************************** !

function PMCBaseCreate()
  ! 
  ! Allocates and initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_base_type), pointer :: PMCBaseCreate
  
  class(pmc_base_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCBase%Create()'
#endif
  
  allocate(pmc)
  call pmc%Init()

  PMCBaseCreate => pmc  
  
end function PMCBaseCreate

! ************************************************************************** !

subroutine PMCBaseInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
#ifdef DEBUG
  print *, 'PMCBase%Init()'
#endif
  
  this%name = 'PMCBase'
  this%stage = 0
  this%is_master = PETSC_FALSE
  nullify(this%option)
  nullify(this%timestepper)
  nullify(this%pms)
  nullify(this%waypoint_list)
  nullify(this%below)
  nullify(this%next)
  nullify(this%sim_aux)
  this%Output => Null()
  
  allocate(this%pm_ptr)
  nullify(this%pm_ptr%ptr)
  
end subroutine PMCBaseInit

! ************************************************************************** !

function PMCCastToBase(this)
  ! 
  ! PMCBaseCastToBase: Casts an extended PMC to a pointer to the base class
  !                    in order to avoid a 'select type' statement when
  !                    pointing a pmc_base_type pointer to an extended class.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_base_type), target :: this
  
  class(pmc_base_type), pointer :: PMCCastToBase

  PMCCastToBase => this
  
end function PMCCastToBase

! ************************************************************************** !

subroutine PMCBaseSetTimestepper(this,timestepper)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_Base_class
  
  implicit none
  
  class(pmc_base_type) :: this
  class(timestepper_base_type), pointer :: timestepper

#ifdef DEBUG
  call printMsg(this%option,'PMCBase%SetTimestepper()')
#endif
  
  this%timestepper => timestepper
  
end subroutine PMCBaseSetTimestepper

! ************************************************************************** !

recursive subroutine InitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
  class(pm_base_type), pointer :: cur_pm
  
#ifdef DEBUG
  call printMsg(this%option,'PMCBase%InitializeRun()')
#endif
  
  if (associated(this%timestepper)) then
    this%option%time = this%timestepper%target_time
  endif
  cur_pm => this%pms
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%InitializeRun()
    cur_pm => cur_pm%next
  enddo
  
  if (associated(this%below)) then
    call this%below%InitializeRun()
  endif
  
  if (associated(this%next)) then
    call this%next%InitializeRun()
  endif

end subroutine InitializeRun

! ************************************************************************** !

recursive subroutine PMCBaseRunToTime(this,sync_time,stop_flag)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_Base_class
  
  implicit none
  
#include "finclude/petscviewer.h"  

  class(pmc_base_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  PetscInt :: local_stop_flag
  PetscBool :: failure
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag
  PetscBool :: checkpoint_flag
  class(pm_base_type), pointer :: cur_pm
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  if (this%stage /= 0) then
    call PetscLogStagePush(this%stage,ierr);CHKERRQ(ierr)
  endif
  this%option%io_buffer = trim(this%name) // ':' // trim(this%pms%name)  
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  call this%GetAuxData()
  
  local_stop_flag = TS_CONTINUE
  do
    if (local_stop_flag /= TS_CONTINUE) exit ! end simulation
    if (this%timestepper%target_time >= sync_time) exit
    
    call SetOutputFlags(this)
    plot_flag = PETSC_FALSE
    transient_plot_flag = PETSC_FALSE
    checkpoint_flag = PETSC_FALSE
    
    call this%timestepper%SetTargetTime(sync_time,this%option, &
                                        local_stop_flag,plot_flag, &
                                        transient_plot_flag,checkpoint_flag)
    call this%timestepper%StepDT(this%pms,local_stop_flag)

    if (local_stop_flag == TS_STOP_FAILURE) exit ! failure
    ! Have to loop over all process models coupled in this object and update
    ! the time step size.  Still need code to force all process models to
    ! use the same time step size if tightly or iteratively coupled.
    cur_pm => this%pms
    do
      if (.not.associated(cur_pm)) exit
      ! have to update option%time for conditions
      this%option%time = this%timestepper%target_time
      call cur_pm%UpdateSolution()
      call this%timestepper%UpdateDT(cur_pm)
      cur_pm => cur_pm%next
    enddo

    ! Accumulate data needed by process-model
    call this%AccumulateAuxData()

    ! Run underlying process model couplers
    if (associated(this%below)) then
      ! Set data needed by process-models
      call this%SetAuxData()
      call this%below%RunToTime(this%timestepper%target_time,local_stop_flag)
      ! Get data from other process-models
      call this%GetAuxData()
    endif
    
    ! only print output for process models of depth 0
    if (associated(this%Output)) then
      if (this%timestepper%time_step_cut_flag) then
        plot_flag = PETSC_FALSE
      endif
      ! however, if we are using the modulus of the output_option%imod, we may
      ! still print
      if (mod(this%timestepper%steps, &
              this%pms% &
                output_option%periodic_output_ts_imod) == 0) then
        plot_flag = PETSC_TRUE
      endif
      if (plot_flag .or. mod(this%timestepper%steps, &
                             this%pms%output_option% &
                               periodic_tr_output_ts_imod) == 0) then
        transient_plot_flag = PETSC_TRUE
      endif
      call this%Output(this%pms%realization_base,plot_flag, &
                       transient_plot_flag)
    endif
    
    if (this%is_master) then
      if (this%option%checkpoint_flag .and. this%option%checkpoint_frequency > 0) then
        if (mod(this%timestepper%steps,this%option%checkpoint_frequency) == 0) then
           checkpoint_flag = PETSC_TRUE
        endif
      endif
    else
      checkpoint_flag = PETSC_FALSE
    endif

    if (checkpoint_flag) then
      ! if checkpointing, need to sync all other PMCs.  Those "below" are
      ! already in sync, but not those "next".
      ! Set data needed by process-model
      call this%SetAuxData()
      ! Run neighboring process model couplers
      if (associated(this%next)) then
        call this%next%RunToTime(this%timestepper%target_time,local_stop_flag)
      endif
      call this%Checkpoint(viewer,this%timestepper%steps)
    endif
    
    if (this%is_master) then
      if (this%timestepper%WallClockStop(this%option)) then
         local_stop_flag = TS_STOP_WALLCLOCK_EXCEEDED
      endif
    endif

  enddo
  
  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%next)) then
    call this%next%RunToTime(sync_time,local_stop_flag)
  endif
  
  stop_flag = max(stop_flag,local_stop_flag)
  
  if (this%stage /= 0) then
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMCBaseRunToTime

! ************************************************************************** !

recursive subroutine PMCBaseUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

  class(pm_base_type), pointer :: cur_pm
  
#ifdef DEBUG
  call printMsg(this%option,'PMCBase%UpdateSolution()')
#endif
  
  cur_pm => this%pms
  do
    if (.not.associated(cur_pm)) exit
    ! have to update option%time for conditions
    this%option%time = this%timestepper%target_time
    call cur_pm%UpdateSolution()
    cur_pm => cur_pm%next
  enddo  
  
end subroutine PMCBaseUpdateSolution

! ************************************************************************** !

recursive subroutine FinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
  character(len=MAXSTRINGLENGTH) :: string
  
#ifdef DEBUG
  call printMsg(this%option,'PMCBase%FinalizeRun()')
#endif
  
  if (associated(this%timestepper)) then
    call this%timestepper%FinalizeRun(this%option)
  endif

  if (associated(this%below)) then
    call this%below%FinalizeRun()
  endif
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif
  
end subroutine FinalizeRun

! ************************************************************************** !

subroutine SetOutputFlags(this)
  ! 
  ! Toggles flags that determine whether output is printed
  ! to the screen and output file during a time step.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/29/13
  ! 

  use Option_module
  use Output_Aux_module
  
  implicit none
  
  class(pmc_base_type) :: this
  
  type(output_option_type), pointer :: output_option
  
  output_option => this%pms%output_option

  if (OptionPrintToScreen(this%option) .and. &
      mod(this%timestepper%steps,output_option%screen_imod) == 0) then
    this%option%print_screen_flag = PETSC_TRUE
  else
    this%option%print_screen_flag = PETSC_FALSE
  endif

  if (OptionPrintToFile(this%option) .and. &
      mod(this%timestepper%steps,output_option%output_file_imod) == 0) then
    this%option%print_file_flag = PETSC_TRUE
  else
    this%option%print_file_flag = PETSC_FALSE
      
  endif
  
end subroutine SetOutputFlags

! ************************************************************************** !

recursive subroutine OutputLocal(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
  class(pm_base_type), pointer :: cur_pm
  
#ifdef DEBUG
  call printMsg(this%option,'PMC%Output()')
#endif
  
  cur_pm => this%pms
  do
    if (.not.associated(cur_pm)) exit
!    call Output(cur_pm%realization,plot_flag,transient_plot_flag)
    cur_pm => cur_pm%next
  enddo
    
end subroutine OutputLocal

! ************************************************************************** !

recursive subroutine PMCBaseCheckpoint(this,viewer,id,id_stamp)
  ! 
  ! Checkpoints PMC timestepper and state variables.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Logging_module
  use Checkpoint_module, only : OpenCheckpointFile, CloseCheckpointFile

  implicit none
  
#include "finclude/petscviewer.h"

  class(pmc_base_type) :: this
  PetscViewer :: viewer
  PetscInt :: id
  character(len=MAXWORDLENGTH), optional, intent(in) :: id_stamp
  
  class(pm_base_type), pointer :: cur_pm
  class(pmc_base_header_type), pointer :: header
  PetscBag :: bag
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  ! if the top PMC, 
  if (this%is_master) then
    call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)
    call PetscLogEventBegin(logging%event_checkpoint,ierr);CHKERRQ(ierr)
    call PetscTime(tstart,ierr);CHKERRQ(ierr)
    if (present(id_stamp)) then
       call OpenCheckpointFile(viewer,id,this%option,id_stamp)
    else
       call OpenCheckpointFile(viewer,id,this%option)
    endif
    ! create header for storing local information specific to PMc
    call PetscBagCreate(this%option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
    call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
    call PMCBaseRegisterHeader(this,bag,header)
    call PMCBaseSetHeader(this,bag,header)
    call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
    call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  endif
  
  call this%timestepper%Checkpoint(viewer,this%option)
  cur_pm => this%pms
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%Checkpoint(viewer)
    cur_pm => cur_pm%next
  enddo
  
  if (associated(this%below)) then
    call this%below%Checkpoint(viewer,UNINITIALIZED_INTEGER)
  endif
  
  if (associated(this%next)) then
    call this%next%Checkpoint(viewer,UNINITIALIZED_INTEGER)
  endif
  
  if (this%is_master) then
    call CloseCheckpointFile(viewer)
    call PetscTime(tend,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer, &
          '("      Seconds to write to checkpoint file: ", f10.2)') &
      tend-tstart
    call printMsg(this%option)
    call PetscLogEventEnd(logging%event_checkpoint,ierr);CHKERRQ(ierr)
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif
    
end subroutine PMCBaseCheckpoint

! ************************************************************************** !

subroutine PMCBaseRegisterHeader(this,bag,header)
  ! 
  ! Register header entries.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/13
  ! 

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

  class(pmc_base_type) :: this
  class(pmc_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  ! bagsize = 2 * 8 bytes = 16 bytes
  call PetscBagRegisterInt(bag,header%plot_number,0, &
                           "plot number","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%times_per_h5_file,0, &
                           "times_per_h5_file","",ierr);CHKERRQ(ierr)

end subroutine PMCBaseRegisterHeader

! ************************************************************************** !

subroutine PMCBaseSetHeader(this,bag,header)
  ! 
  ! Sets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/13
  ! 

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

  class(pmc_base_type) :: this
  class(pmc_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  header%plot_number = &
    this%pms%realization_base%output_option%plot_number

  header%times_per_h5_file = &
    this%pms%realization_base%output_option%times_per_h5_file

end subroutine PMCBaseSetHeader

! ************************************************************************** !

recursive subroutine PMCBaseRestart(this,viewer)
  ! 
  ! Restarts PMC timestepper and state variables.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Logging_module
  use Checkpoint_module, only : OpenCheckpointFile, CloseCheckpointFile

  implicit none
  
#include "finclude/petscviewer.h"

  class(pmc_base_type) :: this
  PetscViewer :: viewer
  character(len=MAXWORDLENGTH) :: filename
  
  class(pm_base_type), pointer :: cur_pm
  class(pmc_base_header_type), pointer :: header
  PetscBag :: bag
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  ! if the top PMC, 
  if (this%is_master) then
    call PetscLogEventBegin(logging%event_restart,ierr);CHKERRQ(ierr)
    call PetscTime(tstart,ierr);CHKERRQ(ierr)
    call PetscViewerBinaryOpen(this%option%mycomm, &
                               this%option%restart_filename, &
                               FILE_MODE_READ,viewer,ierr);CHKERRQ(ierr)
    ! skip reading info file when loading, but not working
    call PetscViewerBinarySetSkipOptions(viewer,PETSC_TRUE,ierr);CHKERRQ(ierr)

    ! read pmc header
    call PetscBagCreate(this%option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
    call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
    call PMCBaseRegisterHeader(this,bag,header)
    call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
    call PMCBaseGetHeader(this,header)
    if (Initialized(this%option%restart_time)) then
      this%pms%realization_base%output_option%plot_number = 0
    endif
    call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  endif
  
  call this%timestepper%Restart(viewer,this%option)
  if (Initialized(this%option%restart_time)) then
    ! simply a flag to set time back to zero, no matter what the restart
    ! time is set to.
    call this%timestepper%Reset()
    ! note that this sets the target time back to zero.
  endif
  
  ! Point cur_waypoint to the correct waypoint.
  !geh: there is a problem here in that the timesteppers "prev_waypoint"
  !     may not be set correctly if the time step does not converge. See
  !     top of TimestepperBaseSetTargetTime().
  call WaypointSkipToTime(this%timestepper%cur_waypoint, &
                          this%timestepper%target_time)
  !geh: this is a bit of a kludge.  Need to use the timestepper target time
  !     directly.  Time is needed to update boundary conditions within 
  !     this%UpdateSolution
  this%option%time = this%timestepper%target_time

  cur_pm => this%pms
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%Restart(viewer)
    cur_pm => cur_pm%next
  enddo
  
  if (associated(this%below)) then
    call this%below%Restart(viewer)
  endif
  
  if (associated(this%next)) then
    call this%next%Restart(viewer)
  endif
  
  if (this%is_master) then
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    call PetscTime(tend,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer, &
          '("      Seconds to read from restart file: ", f10.2)') &
      tend-tstart
    call printMsg(this%option)
    call PetscLogEventEnd(logging%event_restart,ierr);CHKERRQ(ierr)
  endif
    
end subroutine PMCBaseRestart

! ************************************************************************** !

subroutine PMCBaseGetHeader(this,header)
  ! 
  ! Gets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/13
  ! 

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

  class(pmc_base_type) :: this
  class(pmc_base_header_type) :: header
  
  character(len=MAXSTRINGLENGTH) :: string

  this%pms%realization_base%output_option%plot_number = &
    header%plot_number

  ! Check the value of 'times_per_h5_file'
  if (header%times_per_h5_file /= &
      this%pms%realization_base%output_option%times_per_h5_file) then
    write(string,*) header%times_per_h5_file
    this%option%io_buffer = 'From checkpoint file: times_per_h5_file ' // trim(string)
    call printMsg(this%option)
    write(string,*) this%pms%realization_base%output_option%times_per_h5_file
    this%option%io_buffer = 'From inputdeck      : times_per_h5_file ' // trim(string)
    call printMsg(this%option)
    this%option%io_buffer = 'times_per_h5_file specified in inputdeck does not ' // &
      'match that stored in checkpoint file. Correct the inputdeck.'
    call printErrMsg(this%option)
  endif

  this%pms%realization_base%output_option%times_per_h5_file = &
    header%times_per_h5_file

end subroutine PMCBaseGetHeader

! ************************************************************************** !

subroutine AccumulateAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/21/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

end subroutine AccumulateAuxData

! ************************************************************************** !

subroutine GetAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/21/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

end subroutine GetAuxData

! ************************************************************************** !

subroutine SetAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/21/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

end subroutine SetAuxData

! ************************************************************************** !

subroutine PMCBaseUpdateMaterialProperties(this)
  !
  ! At a prescribed time, updates material properties based on instructions
  ! provided by a material update waypoint.
  !
  ! Author: Glenn Hammond
  ! Date: 09/18/14
  
  implicit none
  
  class(pmc_base_type) :: this

end subroutine PMCBaseUpdateMaterialProperties
! ************************************************************************** !

subroutine PMCBaseStrip(this)
  !
  ! Deallocates members of PMC Base.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_base_type) :: this

  nullify(this%option)

  if (associated(this%timestepper)) then
    call this%timestepper%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%timestepper)
    nullify(this%timestepper)
  endif
  if (associated(this%pms)) then
    ! destroy does not currently destroy; it strips
    call this%pms%Destroy()
    deallocate(this%pms)
    nullify(this%pms)
  endif
  nullify(this%waypoint_list) ! deleted in realization
!  call WaypointListDestroy(this%waypoint_list)
  if (associated(this%pm_ptr)) then
    nullify(this%pm_ptr%ptr) ! solely a pointer
    deallocate(this%pm_ptr)
    nullify(this%pm_ptr)
  endif
  nullify(this%sim_aux)

end subroutine PMCBaseStrip

! ************************************************************************** !

recursive subroutine PMCBaseDestroy(this)
  ! 
  ! Deallocates a pmc object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_base_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMC%Destroy()')
#endif
  
  if (associated(this%below)) then
    call this%below%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%below)
    nullify(this%below)
  endif 
  
  if (associated(this%next)) then
    call this%next%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%next)
    nullify(this%next)
  endif 
  
!  deallocate(pmc)
!  nullify(pmc)
  
end subroutine PMCBaseDestroy
  
end module PMC_Base_class
