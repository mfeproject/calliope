!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_arrays

  use mfe_constants, only: wp
  use mfe_types, only: NodeVar, NodeMtx
  implicit none
  private

  type, public :: TwoVec
    real(kind=wp) :: x, u
  end type TwoVec

  integer, save, public :: ncell

  ! Local solution arrays.
  type(NodeVar), dimension(:,:),   allocatable, save, public :: u, udot

  ! Local result arrays.
  type(NodeVar), dimension(:,:),   allocatable, save, public :: r
  type(NodeMtx), dimension(:,:,:), allocatable, save, public :: mtx

  ! Intermediate data arrays.
  real(kind=wp), dimension(:,:),   allocatable, save, public :: l
  type(TwoVec),  dimension(:,:),   allocatable, save, public :: n
  real(kind=wp), dimension(:,:),   allocatable, save, public :: dudx
  real(kind=wp), dimension(:),     allocatable, save, public :: dx
  real(kind=wp), dimension(:,:),   allocatable, save, public :: du

end module local_arrays
