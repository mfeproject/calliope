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
  use mfe_types, only: NodeVar, NodeMtx
  implicit none
  private

  type, public :: TwoVec
    real(r8) :: x, u
  end type TwoVec

  integer, public :: ncell

  ! Local solution arrays.
  type(NodeVar), allocatable, public :: u(:,:), udot(:,:)

  ! Local result arrays.
  type(NodeVar), allocatable, public :: r(:,:)
  type(NodeMtx), allocatable, public :: mtx(:,:,:)

  ! Intermediate data arrays.
  real(r8), allocatable, public :: l(:,:)
  type(TwoVec), allocatable, public :: n(:,:)
  real(r8), allocatable, public :: dudx(:,:)
  real(r8), allocatable, public :: dx(:)
  real(r8), allocatable, public :: du(:,:)

end module local_arrays
