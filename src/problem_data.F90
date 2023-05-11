!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_data

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private
  
  public :: problem_init

  real(r8), public :: visc

contains

  subroutine problem_init(params, stat, errmsg)
    use parameter_list_type
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call params%get('visc', visc, stat=stat, errmsg=errmsg)
  end subroutine

end module problem_data
