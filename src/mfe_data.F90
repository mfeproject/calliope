!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_data

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants
  implicit none
  private

  integer,  public :: kreg
  real(r8), public :: eqw(NEQNS), eltvsc(NEQNS), segspr(NEQNS)
  real(r8), public :: fdinc, dxmin

end module mfe_data
