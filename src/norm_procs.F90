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
  use mfe1_vector_type
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

    type(mfe1_vector), intent(in) :: du
    integer, intent(in) :: key
    real(r8) :: norm

    integer :: k
    real(r8) :: del_dx(du%nnode-1)

    ! Weighted max-norm.
    norm = maxval(abs(du%array(NEQNS+1,:))) / ptol(NEQNS+1)
    do k = 1, NEQNS
      norm = max(norm, maxval(abs(du%array(k,:))) / ptol(k))
    end do

    del_dx = du%array(NEQNS+1,2:du%nnode) - du%array(NEQNS+1,1:du%nnode-1)

    ! Relative norm on element lengths.
    norm = max(norm, maxval(abs(del_dx) / dx) / rtol)

  end function eval_norm

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  CHECK_SOLN
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_soln(u, key, errc)

    type(mfe1_vector), intent(in) :: u
    integer, intent(in)  :: key
    integer, intent(out) :: errc

    integer :: loc

    if (.not.allocated(dx)) allocate(dx(u%nnode-1))

    associate (x => u%array(NEQNS+1,:))
      dx = x(2:size(x)) - x(1:size(x)-1)
    end associate
    loc = minloc(dx, dim=1)

    if (dx(loc) < dxmin) then
      call element_info('CHECK_SOLN: BAD ELEMENT', loc)
      errc = 1
    else
      errc = 0
    end if

  end subroutine check_soln

end module norm_procs
