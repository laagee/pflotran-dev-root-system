module Root_module

  use Dataset_Common_HDF5_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public :: root_property_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: ssd
    PetscReal :: Kcomp
    PetscReal :: Hthreshold
    character(len=MAXWORDLENGTH) :: ssd_dataset_name
    class(dataset_common_hdf5_type), pointer :: ssd_dataset
    type(root_property_type), pointer :: next
  end type root_property_type

  type, public :: root_property_ptr_type
    type(root_property_type), pointer :: ptr
  end type root_property_ptr_type

  public :: RootPropertyCreate, &
            RootPropertyDestroy, &
            RootPropertyAddToList, &
            RootPropertyConvertListToArray, &
            RootPropertyGetPtrFromList, &
            RootPropertyGetPtrFromArray, &
            RootPropertyRead

contains

! ************************************************************************** !
!> This subroutine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/19/13
! ************************************************************************** !
function RootPropertyCreate()

  implicit none

  type(root_property_type), pointer :: RootPropertyCreate

  type(root_property_type), pointer :: root_property

  allocate(root_property)
  root_property%id = 0
  root_property%name = ''
  root_property%ssd = 0.d0
  root_property%Kcomp = 0.d0
  root_property%Hthreshold = 0.d0
  root_property%ssd_dataset_name = ''
  nullify(root_property%ssd_dataset)
  nullify(root_property%next)
  RootPropertyCreate => root_property

end function RootPropertyCreate

! ************************************************************************** !
!> This subroutine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/19/13
! ************************************************************************** !
subroutine RootPropertyRead(root_property,input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(root_property_type) :: root_property
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: length

  input%ierr = 0
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','ROOT_PROPERTY')
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('NAME')
        call InputReadWord(input,option,root_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','root_property')
      case('ID')
        call InputReadInt(input,option,root_property%id)
        call InputErrorMsg(input,option,'id','root_property')
      case('COMPENSATORY_CONDUCATNACE')
        call InputReadDouble(input,option,root_property%Kcomp)
        call InputErrorMsg(input,option,'Kcomp','root_property')
      case('STANDARD_WATER_UPDATE_DENISTY')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'dataset or value','SSD')
        call StringToLower(word)
        length = len_trim(word)
        if (length == SEVEN_INTEGER .and. StringCompare(word,'dataset')) then
          call InputReadWord(input,option,root_property%ssd_dataset_name,PETSC_TRUE)
          call InputErrorMsg(input,option,'dataset','root_property')
        else
          input%buf = trim(word)
          call InputReadDouble(input,option,root_property%ssd)
          call InputErrorMsg(input,option,'value','SSD')
        endif
      case('COLLAR_WATER_POT_THRESHOLD')
        call InputReadDouble(input,option,root_property%Hthreshold)
        call InputErrorMsg(input,option,'Hthreshold','root_property')
      case default
        option%io_buffer = 'Keyword (' // trim(keyword) // &
                           ') not recognized in root_property'
        call printErrMsg(option)
    end select
  enddo

end subroutine

! ************************************************************************** !
!> This subroutine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/19/13
! ************************************************************************** !
subroutine RootPropertyAddToList(root_property,list)

  implicit none

  type(root_property_type), pointer :: root_property
  type(root_property_type), pointer :: list

  type(root_property_type), pointer :: cur_root_property

  if (associated(list)) then
    cur_root_property => list
    do
      if (.not.associated(cur_root_property%next)) exit
      cur_root_property => cur_root_property%next
    enddo
    cur_root_property%next => root_property
  else
    list => root_property
  endif

end subroutine RootPropertyAddToList

! ************************************************************************** !
!> This subroutine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/19/13
! ************************************************************************** !
subroutine RootPropertyConvertListToArray(list,array,option)

  use Option_module
  use String_module

  implicit none

  type(root_property_type), pointer :: list
  type(root_property_ptr_type), pointer :: array(:)
  type(option_type) :: option

  type(root_property_type), pointer :: cur_root_property
  type(root_property_type), pointer :: prev_root_property
  type(root_property_type), pointer :: next_root_property
  PetscInt :: i, j, length1,length2, max_id
  PetscInt, allocatable :: id_count(:)
  PetscBool :: error_flag
  character(len=MAXSTRINGLENGTH) :: string

  max_id = 0
  cur_root_property => list
  do 
    if (.not.associated(cur_root_property)) exit
    max_id = max(max_id,cur_root_property%id)
    cur_root_property => cur_root_property%next
  enddo
  
  allocate(array(max_id))
  do i = 1, max_id
    nullify(array(i)%ptr)
  enddo
  
  ! use id_count to ensure that an id is not duplicated
  allocate(id_count(max_id))
  id_count = 0
  
  cur_root_property => list
  do 
    if (.not.associated(cur_root_property)) exit
    id_count(cur_root_property%id) = &
      id_count(cur_root_property%id) + 1
    array(cur_root_property%id)%ptr => cur_root_property
    cur_root_property => cur_root_property%next
  enddo

  ! check to ensure that an id is not duplicated
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (id_count(i) > 1) then
      write(string,*) i
      option%io_buffer = 'Root ID ' // trim(adjustl(string)) // &
        ' is duplicated in input file.'
      call printMsg(option)
      error_flag = PETSC_TRUE
    endif
  enddo

  deallocate(id_count)

  if (error_flag) then
    option%io_buffer = 'Duplicate Root IDs.'
    call printErrMsg(option)
  endif
  
  ! ensure unique material names
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (associated(array(i)%ptr)) then
      length1 = len_trim(array(i)%ptr%name)
      do j = 1, i-1
        if (associated(array(j)%ptr)) then
          length2 = len_trim(array(j)%ptr%name)
          if (length1 /= length2) cycle
          if (StringCompare(array(i)%ptr%name,array(j)%ptr%name,length1)) then
            option%io_buffer = 'Root name "' // &
              trim(adjustl(array(i)%ptr%name)) // &
              '" is duplicated in input file.'
            call printMsg(option)
            error_flag = PETSC_TRUE
          endif
        endif
      enddo
    endif
  enddo

  if (error_flag) then
    option%io_buffer = 'Duplicate Root names.'
    call printErrMsg(option)
  endif
  
end subroutine RootPropertyConvertListToArray

! ************************************************************************** !
!> This subroutine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/19/13
! ************************************************************************** !
function RootPropertyGetPtrFromList(root_property_name,root_property_list)

  use String_module

  implicit none

  type(root_property_type), pointer :: RootPropertyGetPtrFromList
  character(len=MAXWORDLENGTH) :: root_property_name
  type(root_property_type), pointer :: root_property_list

  PetscInt :: length
  type(root_property_type), pointer :: root_property

  nullify(RootPropertyGetPtrFromList)
  root_property => root_property_list

  length = len_trim(root_property_name)
  do
    if (.not.associated(root_property)) exit
    if (length == len_trim(root_property%name) .and. &
        StringCompare(root_property%name,root_property_name,length)) then
        RootPropertyGetPtrFromList => root_property
      return
    endif
  enddo

end function RootPropertyGetPtrFromList

! ************************************************************************** !
!> This subroutine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/19/13
! ************************************************************************** !
function RootPropertyGetPtrFromArray(root_property_name,root_property_array)

  use String_module

  implicit none

  type(root_property_type), pointer :: RootPropertyGetPtrFromArray
  character(len=MAXWORDLENGTH) :: root_property_name
  type(root_property_ptr_type), pointer :: root_property_array(:)

  PetscInt :: length
  PetscInt :: iroot_property

  nullify(RootPropertyGetPtrFromArray)

  length = len_trim(root_property_name)
  do iroot_property = 1, size(root_property_array)
    if (.not.associated(root_property_array(iroot_property)%ptr)) cycle
    if (length == len_trim(root_property_array(iroot_property)%ptr%name) .and. &
        StringCompare(root_property_array(iroot_property)%ptr%name, &
                      root_property_name,length)) then
        RootPropertyGetPtrFromArray => root_property_array(iroot_property)%ptr
      return
    endif
  enddo

end function RootPropertyGetPtrFromArray

! ************************************************************************** !
!> This subroutine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/19/13
! ************************************************************************** !
recursive subroutine RootPropertyDestroy(root_property)

  implicit none

  type(root_property_type), pointer :: root_property

  if (.not.associated(root_property)) return

  call RootPropertyDestroy(root_property%next)

  deallocate(root_property)
  nullify(root_property)

end subroutine

end module Root_module
