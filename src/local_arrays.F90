!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_arrays

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  integer, public :: ncell, ix

  ! Local solution arrays.
  real(r8), allocatable, public :: u(:,:,:), udot(:,:,:)

  ! Local result arrays.
  real(r8), allocatable, public :: r(:,:,:)
  real(r8), allocatable, public :: mtx(:,:,:,:,:)

  ! Intermediate data arrays.
  real(r8), allocatable, public :: l(:,:)
  real(r8), allocatable, public :: n(:,:,:)
  real(r8), allocatable, public :: dudx(:,:)
  real(r8), allocatable, public :: dx(:)
  real(r8), allocatable, public :: du(:,:)

end module local_arrays
