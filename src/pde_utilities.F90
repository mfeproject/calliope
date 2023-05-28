#include "f90_assert.fpp"

module pde_utilities

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use cell_data_type
  implicit none
  private

  public :: add_lapl

  real(r8), parameter :: ETA = 0.01_r8, &
    C3 = 1.0_r8 / 3.0_r8, C5 = 1.0_r8 / 5.0_r8, C7 = 1.0_r8 / 7.0_r8

contains

  pure subroutine add_lapl(cdata, i, c, gx, gu)

    type(cell_data(*)), intent(in) :: cdata
    integer, intent(in) :: i
    real(r8), intent(in) :: c(:)
    real(r8), intent(inout) :: gx(:), gu(:,:)

    real(r8) :: r1, m, e, s1, s2

    r1 = 1.0_r8 / cdata%nu(i)
    m = cdata%dudx(i)

    e = m / r1
    if (abs(e) > ETA) then
      s1 = -sign(log(abs(m)+r1), m)
    else
      s1 = -e*(1 + (e**2)*(C3 + (e**2)*(C5 + C7*(e**2))))
    end if

    s2 = m**2/(1 + r1)

    gx(1)   = gx(1)   - (c(1)*s2)
    gu(i,1) = gu(i,1) - (c(1)*s1)

    gx(2)   = gx(2)   + (c(2)*s2)
    gu(i,2) = gu(i,2) + (c(2)*s1)

  end subroutine

end module pde_utilities
