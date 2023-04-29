!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_pde

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use problem_data
  use mfe_data, only: eqw
  use local_arrays
  use local_laplacian
  implicit none
  private

  public  :: pde_rhs
  private :: flux

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  LOAD_PDE_RHS -- Inner products for gasdynamics.
 !!
 !!        Equation 1 == mass density,
 !!        Equation 2 == momentum density,
 !!        Equation 3 == total energy density (internal plus kinetic).
 !!
 !!     Each equation is in conservation law form
 !!
 !!        du/dt = -df/dx + visc* d2u/dx2.
 !!
 !!  FLUX -- Gas dynamics flux f.
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pde_rhs(t)

    real(r8), intent(in) :: t

    ! Local variables
    integer :: i, k
    real(r8) :: rx1, rx2, term
    real(r8) :: umid(3), f1(3), f2(3), favg(3)
    real(r8), parameter :: c1 = 1.0_r8 / 6.0_r8, c2 = 4.0_r8 / 6.0_r8

    do i = 1, ncell

      f1 = flux(u(:3,1,i))                      ! Flux at left endpoint.
      f2 = flux(u(:3,2,i))                      ! Flux at right endpoint.
      umid = 0.5_r8 * (u(:3,1,i) + u(:3,2,i))   ! Variables at midpoint.
      favg = c1 * (f1 + f2) + c2 * flux(umid)   ! Average flux (Simpson).

      rx1 = 0.0_r8
      rx2 = 0.0_r8

      do k = 1, 3

        term = (eqw(k) * (favg(k) - f1(k)))
        rx1 = rx1 - term * n(1,k,i)
        r(k,1,i) = - term * n(2,k,i)

        term = (eqw(k) * (f2(k) - favg(k)))
        rx2 = rx2 - term * n(1,k,i)
        r(k,2,i) = - term * n(2,k,i)

      end do

      r(ix,1,i) = rx1
      r(ix,2,i) = rx2

    end do

    call laplacian(eqno=1, coef=visc)
    call laplacian(eqno=2, coef=visc)
    call laplacian(eqno=3, coef=visc)

  end subroutine pde_rhs


  function flux(u) result(f)

    real(r8), intent(in) :: u(:)
    real(r8) :: f(3)

    real(r8), parameter :: gamma = 1.4_r8,                     &
                                c1 = 0.5_r8 * (3.0_r8 - gamma),     &
                                c2 = gamma - 1.0_r8,                &
                                c3 = 0.5_r8 * (1.0_r8 - gamma),     &
                                c4 = gamma

    ! Mass flux (equation 1).
    f(1) = u(2)

    ! Momentum flux (equation 2).
    f(2) = c1 * (u(2) * (u(2) / u(1))) + c2 * u(3)

    ! Energy flux (equation 3).
    f(3) = (c3 * (u(2) * (u(2) / u(1))) + c4 * u(3)) * (u(2) / u(1))

  end function flux

end module problem_pde
