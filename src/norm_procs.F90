!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module norm_procs

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants
  use mfe_types, only: NodeVar
  use mfe_data, only: dxmin
  use norm_data
  use common_io, only: element_info
  implicit none
  private

  public :: eval_norm, check_soln

  real(r8), allocatable, save, public :: dx(:)

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  EVAL_NORM
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function eval_norm(du, key) result(norm)

    type(NodeVar(*)), intent(in) :: du(:)
    integer, intent(in) :: key
    real(r8) :: norm

    integer :: k
    real(r8) :: del_dx(size(du)-1)

    ! Weighted max-norm.
    norm = maxval(abs(du%x)) / ptol%x
    do k = 1, NEQNS
      norm = max(norm, maxval(abs(du%u(k))) / ptol%u(k))
    end do

    del_dx = du(2:size(du))%x - du(1:size(du)-1)%x

    ! Relative norm on element lengths.
    norm = max(norm, maxval(abs(del_dx) / dx) / rtol)

  end function eval_norm

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  CHECK_SOLN
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_soln(u, key, errc)

    type(NodeVar(*)), intent(in) :: u(:)
    integer, intent(in)  :: key
    integer, intent(out) :: errc

    integer :: loc

    if (.not. allocated(dx)) then
      allocate (dx(size(u)-1))
    end if

    dx = u(2:size(u))%x - u(1:size(u)-1)%x
    loc = minloc(dx, dim=1)

    if (dx(loc) < dxmin) then
      call element_info('CHECK_SOLN: BAD ELEMENT', loc)
      errc = 1
    else
      errc = 0
    end if

  end subroutine check_soln

end module norm_procs
