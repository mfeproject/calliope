!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_constants

  use kind_parameters, only: r8
  implicit none
  private

  integer, parameter, public :: wp = r8                 ! Use 8-byte reals
  
  integer, parameter, public :: NEQNS = 3               ! Number equations in the PDE system

  !!! DO NOT MODIFY THESE PARAMETERS !!!
  integer, parameter, public :: DIMEN = 1               ! Spatial dimension
  integer, parameter, public :: NVARS = NEQNS + DIMEN   ! Total number of variables per node
  integer, parameter, public :: NVERT = DIMEN + 1       ! Number of vertices per cell
  integer, parameter, public :: NFACE = NVERT           ! Number of faces per cell

end module mfe_constants
