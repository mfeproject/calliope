!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_laplacian

  use mfe_constants
  use mfe_data, only: eqw
  use local_arrays
  implicit none
  private

  public  :: laplacian
  private :: lapl_const_coef, lapl_var_coef

  interface laplacian
    module procedure lapl_const_coef, lapl_var_coef
  end interface

  real(kind=wp), parameter :: ETA = 0.01_wp, &
    C3 = 1.0_wp / 3.0_wp, C5 = 1.0_wp / 5.0_wp, C7 = 1.0_wp / 7.0_wp

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! LAPL_CONST_COEF (LAPLACIAN)
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lapl_const_coef (eqno, coef)

      integer, intent(in) :: eqno
      real(kind=wp), intent(in) :: coef

      integer :: i
      real(kind=wp) :: c, r1, m, e, s1, s2

      c = eqw(eqno) * coef

      do i = 1, ncell

        r1 = 1.0_wp / n(eqno,i) % u
        m = dudx(eqno,i)

       !!!
       !!! FIRST- AND SECOND-KIND INTEGRALS

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_wp + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_wp + r1)

       !!!
       !!! LOAD THE INNER PRODUCTS

        r(1,i) % x       = r(1,i) % x       - (c * s2)
        r(1,i) % u(eqno) = r(1,i) % u(eqno) - (c * s1)

        r(2,i) % x       = r(2,i) % x       + (c * s2)
        r(2,i) % u(eqno) = r(2,i) % u(eqno) + (c * s1)

      end do

    end subroutine lapl_const_coef

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! LAPL_VAR_COEF (LAPLACIAN)
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lapl_var_coef (eqno, coef)

      integer, intent(in) :: eqno
      real(kind=wp), dimension(:,:), intent(in) :: coef

      integer :: i
      real(kind=wp) :: c, r1, m, e, s1, s2

      do i = 1, ncell

        r1 = 1.0_wp / n(eqno,i) % u
        m = dudx(eqno,i)

       !!!
       !!! FIRST- AND SECOND-KIND INTEGRALS

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_wp + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_wp + r1)

       !!!
       !!! LOAD THE INNER PRODUCTS

        c = eqw(eqno) * coef(1,i)
        r(1,i) % x       = r(1,i) % x       - c * s2
        r(1,i) % u(eqno) = r(1,i) % u(eqno) - c * s1

        c = eqw(eqno) * coef(2,i)
        r(2,i) % x       = r(2,i) % x       + c * s2
        r(2,i) % u(eqno) = r(2,i) % u(eqno) + c * s1

      end do

    end subroutine lapl_var_coef

end module local_laplacian
