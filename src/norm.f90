!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module norm_data

  use mfe_constants
  use mfe_types, only: NodeVar
  implicit none
  private

  type(NodeVar), save, public :: ptol
  real(kind=wp), save, public :: rtol

end module norm_data

module norm_procs

  use mfe_constants
  use mfe_types, only: NodeVar
  use mfe_data, only: dxmin
  use norm_data
  use common_io, only: element_info
  implicit none
  private

  public :: eval_norm, check_soln

  real(kind=wp), dimension(:), allocatable, save, public :: dx

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_NORM
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function eval_norm (du, key) result (norm)

      type(NodeVar), dimension(:), intent(in) :: du
      integer, intent(in) :: key
      real(kind=wp) :: norm

      integer :: k
      real(kind=wp), dimension(size(du)-1) :: del_dx

      ! Weighted max-norm.
      norm = maxval(abs(du % x)) / ptol % x
      do k = 1, NEQNS
        norm = max (norm, maxval(abs(du % u(k))) / ptol % u(k))
      end do

      del_dx = du(2:size(du)) % x - du(1:size(du)-1) % x

      ! Relative norm on element lengths.
      norm = max (norm, maxval(abs(del_dx) / dx) / rtol)

    end function eval_norm

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  CHECK_SOLN
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine check_soln (u, key, errc)

      type(NodeVar), dimension(:), intent(in) :: u
      integer, intent(in)  :: key
      integer, intent(out) :: errc

      integer, dimension(1) :: loc

      if (.not. allocated(dx)) then
        allocate (dx(size(u)-1))
      end if

      dx = u(2:size(u)) % x - u(1:size(u)-1) % x
      loc = minloc(dx)

      if (dx(loc(1)) < dxmin) then
        call element_info ("CHECK_SOLN: BAD ELEMENT", loc(1))
        errc = 1
      else
        errc = 0
      end if

    end subroutine check_soln

end module norm_procs
