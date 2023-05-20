#include "f90_assert.fpp"

module mfe1_disc_core_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: mfe1_disc_core
    integer :: neqns, nvars, ncell
    real(r8), allocatable :: u(:,:,:), udot(:,:,:)
    real(r8), allocatable :: dx(:), du(:,:)
    real(r8), allocatable :: dudx(:,:)
    real(r8), allocatable :: l(:,:)
    real(r8), allocatable :: n(:,:,:)
    real(r8), allocatable :: r(:,:,:)
    real(r8), allocatable :: mtx(:,:,:,:,:)
    real(r8), allocatable :: eqw(:)
  contains
    generic   :: laplacian => lapl_const_coef, lapl_var_coef
    procedure, private :: lapl_const_coef, lapl_var_coef
  end type

  real(r8), parameter :: ETA = 0.01_r8, &
    C3 = 1.0_r8 / 3.0_r8, C5 = 1.0_r8 / 5.0_r8, C7 = 1.0_r8 / 7.0_r8

contains

  subroutine lapl_const_coef(this, eqno, coef)

    class(mfe1_disc_core), intent(inout) :: this
    integer, intent(in) :: eqno
    real(r8), intent(in) :: coef

    integer :: i
    real(r8) :: c, r1, m, e, s1, s2

    associate (rx => this%r(this%nvars,:,:))
      c = this%eqw(eqno) * coef
      do i = 1, this%ncell
        r1 = 1.0_r8 / this%n(2,eqno,i)
        m = this%dudx(eqno,i)

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_r8 + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_r8 + r1)

        rx(1,i) = rx(1,i) - (c * s2)
        this%r(eqno,1,i) = this%r(eqno,1,i) - (c * s1)

        rx(2,i) = rx(2,i) + (c * s2)
        this%r(eqno,2,i) = this%r(eqno,2,i) + (c * s1)
      end do
    end associate

  end subroutine


  subroutine lapl_var_coef(this, eqno, coef)

    class(mfe1_disc_core), intent(inout) :: this
    integer, intent(in) :: eqno
    real(r8), intent(in) :: coef(:,:)

    integer :: i
    real(r8) :: c, r1, m, e, s1, s2

    associate (rx => this%r(this%nvars,:,:))
      do i = 1, this%ncell
        r1 = 1.0_r8 / this%n(2,eqno,i)
        m = this%dudx(eqno,i)

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_r8 + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_r8 + r1)

        c = this%eqw(eqno) * coef(1,i)
        rx(1,i) = rx(1,i) - c * s2
        this%r(eqno,1,i) = this%r(eqno,1,i) - c * s1

        c = this%eqw(eqno) * coef(2,i)
        rx(2,i) = rx(2,i) + c * s2
        this%r(eqno,2,i) = this%r(eqno,2,i) + c * s1

      end do
    end associate

  end subroutine

end module
