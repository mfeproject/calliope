!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_laplacian

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_data, only: eqw
  use local_arrays
  implicit none
  private

  public  :: laplacian
  private :: lapl_const_coef, lapl_var_coef

  interface laplacian
    module procedure lapl_const_coef, lapl_var_coef
  end interface

  real(r8), parameter :: ETA = 0.01_r8, &
    C3 = 1.0_r8 / 3.0_r8, C5 = 1.0_r8 / 5.0_r8, C7 = 1.0_r8 / 7.0_r8

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! LAPL_CONST_COEF (LAPLACIAN)
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lapl_const_coef(eqno, coef)

    integer, intent(in) :: eqno
    real(r8), intent(in) :: coef

    integer :: i
    real(r8) :: c, r1, m, e, s1, s2

    c = eqw(eqno) * coef

    do i = 1, ncell

      r1 = 1.0_r8 / n(2,eqno,i)
      m = dudx(eqno,i)

     !!!
     !!! FIRST- AND SECOND-KIND INTEGRALS

      e = m / r1
      if (abs(e) > ETA) then
        s1 = - sign ( log(abs(m) + r1), m )
      else
        s1 = - e * (1.0_r8 + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
      end if

      s2 = m ** 2 / (1.0_r8 + r1)

     !!!
     !!! LOAD THE INNER PRODUCTS

      r(ix,1,i)   = r(ix,1,i)   - (c * s2)
      r(eqno,1,i) = r(eqno,1,i) - (c * s1)

      r(ix,2,i)   = r(ix,2,i)   + (c * s2)
      r(eqno,2,i) = r(eqno,2,i) + (c * s1)

    end do

  end subroutine lapl_const_coef

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! LAPL_VAR_COEF (LAPLACIAN)
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lapl_var_coef(eqno, coef)

    integer, intent(in) :: eqno
    real(r8), intent(in) :: coef(:,:)

    integer :: i
    real(r8) :: c, r1, m, e, s1, s2

    do i = 1, ncell

      r1 = 1.0_r8 / n(2,eqno,i)
      m = dudx(eqno,i)

     !!!
     !!! FIRST- AND SECOND-KIND INTEGRALS

      e = m / r1
      if (abs(e) > ETA) then
        s1 = - sign ( log(abs(m) + r1), m )
      else
        s1 = - e * (1.0_r8 + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
      end if

      s2 = m ** 2 / (1.0_r8 + r1)

     !!!
     !!! LOAD THE INNER PRODUCTS

      c = eqw(eqno) * coef(1,i)
      r(ix,1,i)   = r(ix,1,i)   - c * s2
      r(eqno,1,i) = r(eqno,1,i) - c * s1

      c = eqw(eqno) * coef(2,i)
      r(ix,2,i)   = r(ix,2,i)   + c * s2
      r(eqno,2,i) = r(eqno,2,i) + c * s1

    end do

  end subroutine lapl_var_coef

end module local_laplacian
