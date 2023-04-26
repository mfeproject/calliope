!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_constants

  implicit none
  public

  integer, parameter :: NEQNS = 3               ! Number equations in the PDE system

  !!! DO NOT MODIFY THESE PARAMETERS !!!
  integer, parameter :: DIMEN = 1               ! Spatial dimension
  integer, parameter :: NVARS = NEQNS + DIMEN   ! Total number of variables per node
  integer, parameter :: NVERT = DIMEN + 1       ! Number of vertices per cell
  integer, parameter :: NFACE = NVERT           ! Number of faces per cell

end module mfe_constants
