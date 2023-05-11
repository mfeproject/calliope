!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module initialize

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants
  implicit none
  private

  public :: refine

contains

  subroutine refine(x0, n, x)

    use mfe1_vector_type

    integer, intent(in) :: n(:)
    real(r8), intent(in) :: x0(:,:)
    type(mfe1_vector), intent(inout) :: x

    integer :: node, j, k
    real(r8) :: dx, du(NEQNS)

    node = 1

    associate (xx => x%array(NEQNS+1,:), u => x%array(:NEQNS,:))
      do k = 1, size(n)
        u(:,node) = x0(2:NVARS,k)
        xx(node) = x0(1,k)
        node = node + 1

        dx = (x0(1,k+1) - x0(1,k)) / n(k)
        du = (x0(2:NVARS,k+1) - x0(2:NVARS,k)) / n(k)

        do j = 1, n(k) - 1
          xx(node) = xx(node-1) + dx
          u(:,node) = u(:,node-1) + du
          node = node + 1
        end do
      end do
      u(:,node) = x0(2:NVARS,size(n)+1)
      xx(node) = x0(1,size(n)+1)
    end associate

  end subroutine refine

end module initialize

