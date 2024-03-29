module Hydrogeophysics_Simulation_class
  
  use Option_module
  use Subsurface_Simulation_class
  use PMC_Hydrogeophysics_class

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type, public, extends(subsurface_simulation_type) :: &
    hydrogeophysics_simulation_type
    ! pointer to hydrogeophysics coupler
    class(pmc_hydrogeophysics_type), pointer :: hydrogeophysics_coupler
    PetscMPIInt :: pf_e4d_scatter_comm
    PetscMPIInt :: pf_e4d_scatter_grp
    PetscMPIInt :: pf_e4d_scatter_size
    PetscMPIInt :: pf_e4d_scatter_rank
    PetscMPIInt :: pf_e4d_master_comm
    PetscMPIInt :: pf_e4d_master_grp
    PetscMPIInt :: pf_e4d_master_size
    PetscMPIInt :: pf_e4d_master_rank
    PetscBool :: pflotran_process
    Vec :: solution_mpi
    ! these PetscMPIInts are used to save the process decomposition that 
    ! enters hydrogeophysics_simulation.F90:HydrogeophysicsInitialize()
    ! - they are set in HydrogeophysicsInitialize() 
    ! - their counterparts in option are set back in HydrogeophysicsDestroy()
    PetscMPIInt :: mycomm_save
    PetscMPIInt :: myrank_save
    PetscMPIInt :: mycommsize_save
    PetscMPIInt :: mygroup_save
    PetscMPIInt :: mygroup_id_save
  contains
    procedure, public :: Init => HydrogeophysicsInit
    procedure, public :: ExecuteRun => HydrogeophysicsExecuteRun
!    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => HydrogeophysicsFinalizeRun
    procedure, public :: Strip => HydrogeophysicsStrip
  end type hydrogeophysics_simulation_type
  
  public :: HydrogeophysicsCreate, &
            HydrogeophysicsDestroy
  
contains

! ************************************************************************** !

function HydrogeophysicsCreate(option)
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/17/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(hydrogeophysics_simulation_type), pointer :: HydrogeophysicsCreate
  
  call printMsg(option,'HydrogeophysicsCreate()')
  
  allocate(HydrogeophysicsCreate)
  call HydrogeophysicsCreate%Init(option)
  
end function HydrogeophysicsCreate

! ************************************************************************** !

subroutine HydrogeophysicsInit(this,option)
  ! 
  ! Initializes simulation values
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/17/13
  ! 

  use Option_module
  
  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  type(option_type), pointer :: option
  
  call SubsurfaceSimulationInit(this,option)
  nullify(this%hydrogeophysics_coupler)
  this%solution_mpi = 0
  ! UNINITIALIZED_INTEGER denotes uninitialized
  this%pf_e4d_scatter_comm = MPI_COMM_NULL
  this%pf_e4d_scatter_grp = UNINITIALIZED_INTEGER
  this%pf_e4d_scatter_size = UNINITIALIZED_INTEGER
  this%pf_e4d_scatter_rank = UNINITIALIZED_INTEGER
  this%pf_e4d_master_comm = MPI_COMM_NULL
  this%pf_e4d_master_grp = UNINITIALIZED_INTEGER
  this%pf_e4d_master_size = UNINITIALIZED_INTEGER
  this%pf_e4d_master_rank = UNINITIALIZED_INTEGER
  this%pflotran_process = PETSC_FALSE
  this%mycomm_save = 0
  this%myrank_save = 0
  this%mycommsize_save = 0
  this%mygroup_save = 0
  this%mygroup_id_save = 0
   
end subroutine HydrogeophysicsInit

! ************************************************************************** !

subroutine HydrogeophysicsExecuteRun(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Simulation_Base_class

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  PetscReal :: final_time
  PetscReal :: dt
  PetscReal :: current_time
  
  call printMsg(this%option,'Hydrogeophysics%ExecuteRun()')

  if (this%pflotran_process) then
    final_time = SimulationGetFinalWaypointTime(this)
    ! take hourly steps until final time
    current_time = 0.d0
    dt = 365.d0*24.d0*3600.d0
    do
      current_time = min(current_time + dt,final_time)
      call this%RunToTime(current_time)
      if (this%stop_flag > 0) exit
    enddo
  else
    ! do nothing for E4D as it is waiting to receive instructions
  endif
  
end subroutine HydrogeophysicsExecuteRun

! ************************************************************************** !

subroutine HydrogeophysicsFinalizeRun(this)
  ! 
  ! Finalizes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/17/13
  ! 

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'Hydrogeophysics%FinalizeRun()')
  
  if (this%pflotran_process) then
    call SubsurfaceFinalizeRun(this)
  endif
  
end subroutine HydrogeophysicsFinalizeRun

! ************************************************************************** !

subroutine HydrogeophysicsStrip(this)
  ! 
  ! Deallocates members of hydrogeophysics simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Hydrogeophysics_Wrapper_module, only : HydrogeophysicsWrapperDestroy

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'Hydrogeophysics%Strip()')
  
  ! we place a barrier here to ensure that the pflotran processes do not
  ! rate ahead.
  call MPI_Barrier(this%mycomm_save,ierr)
  if (this%pflotran_process) then
    call SubsurfaceSimulationStrip(this)
  else
    call HydrogeophysicsWrapperDestroy(this%option)
  endif
  ! created in HydrogeophysicsInitialize()
  if (this%solution_mpi /= 0) then
    call VecDestroy(this%solution_mpi ,ierr);CHKERRQ(ierr)
  endif
  this%solution_mpi = 0

  if (this%pf_e4d_scatter_comm /= MPI_COMM_NULL)  then
    call MPI_Comm_free(this%pf_e4d_scatter_comm,ierr)
  endif
  this%pf_e4d_scatter_comm = 0
  if (this%pf_e4d_master_comm /= MPI_COMM_NULL)  then
    call MPI_Comm_free(this%pf_e4d_master_comm,ierr)
  endif
  this%pf_e4d_master_comm = 0

  ! reset original communicator back to initial state  
  this%option%mycomm = this%mycomm_save
  this%option%myrank = this%myrank_save
  this%option%mycommsize = this%mycommsize_save
  this%option%mygroup = this%mygroup_save
  this%option%mygroup_id = this%mygroup_id_save
  
end subroutine HydrogeophysicsStrip

! ************************************************************************** !

subroutine HydrogeophysicsDestroy(simulation)
  ! 
  ! Deallocates a simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/17/13
  ! 

  implicit none
  
  class(hydrogeophysics_simulation_type), pointer :: simulation

  call printMsg(simulation%option,'Hydrogeophysics%Destroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine HydrogeophysicsDestroy
  
end module Hydrogeophysics_Simulation_class
