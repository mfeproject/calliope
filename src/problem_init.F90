!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_init

  use problem_data
  use common_io, only: read_tagged_data
  implicit none
  private
  
  public :: read_problem_data
  
contains

  subroutine read_problem_data(params)
    use parameter_list_type
    type(parameter_list), intent(inout) :: params
    call read_tagged_data(visc, 'PDE viscosity coefficient')
    call params%set('visc', visc)
  end subroutine

end module problem_init
