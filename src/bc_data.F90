!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bc_data

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants, only: NEQNS
  implicit none
  private

  type, public :: NodeBC
    integer :: x_type
    integer :: u_type(NEQNS)
    real(r8) :: x_value
    real(r8) :: u_value(NEQNS)
  end type NodeBC

  integer, parameter, public :: FREE  = 0
  integer, parameter, public :: FIXED = 1

  type(NodeBC), save, public :: bc_left, bc_right

end module bc_data
