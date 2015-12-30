!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module problem_data

  use mfe_constants, only: wp
  implicit none
  private
  
  real(kind=wp), save, public :: visc
  
end module problem_data
  
module problem_init

  use problem_data
  use common_io
  implicit none
  private
  
  public :: read_problem_data
  
  contains
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_PROBLEM_DATA
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    subroutine read_problem_data ()

      call read_tagged_data (visc, "PDE viscosity coefficient")

    end subroutine read_problem_data

end module problem_init

module problem_pde

  use problem_data
  use mfe_constants, only: wp
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

     subroutine pde_rhs (t)

       real(kind=wp), intent(in) :: t

       ! Local variables
       integer :: i, k
       real(kind=wp) :: rx1, rx2, term
       real(kind=wp), dimension(3) :: umid, f1, f2, favg
       real(kind=wp), parameter :: c1 = 1.0_wp / 6.0_wp, c2 = 4.0_wp / 6.0_wp

       do i = 1, ncell

         f1 = flux (u(1,i) % u)                      ! Flux at left endpoint.
         f2 = flux (u(2,i) % u)                      ! Flux at right endpoint.
         umid = 0.5_wp * (u(1,i) % u + u(2,i) % u)   ! Variables at midpoint.
         favg = c1 * (f1 + f2) + c2 * flux (umid)    ! Average flux (Simpson).

         rx1 = 0.0_wp
         rx2 = 0.0_wp

         do k = 1, 3

           term = (eqw(k) * (favg(k) - f1(k)))         
           rx1 = rx1 - term * n(k,i) % x
           r(1,i) % u(k) = - term * n(k,i) % u

           term = (eqw(k) * (f2(k) - favg(k)))
           rx2 = rx2 - term * n(k,i) % x
           r(2,i) % u(k) = - term * n(k,i) % u

         end do

         r(1,i) % x = rx1
         r(2,i) % x = rx2

       end do

       call laplacian (eqno=1, coef=visc)
       call laplacian (eqno=2, coef=visc)
       call laplacian (eqno=3, coef=visc)

     end subroutine pde_rhs


     function flux (u) result (f)

       real(kind=wp), dimension(:), intent(in) :: u
       real(kind=wp), dimension(3) :: f

       real(kind=wp), parameter :: gamma = 1.4_wp,                     &
                                   c1 = 0.5_wp * (3.0_wp - gamma),     &
                                   c2 = gamma - 1.0_wp,                &
                                   c3 = 0.5_wp * (1.0_wp - gamma),     &
                                   c4 = gamma

       ! Mass flux (equation 1).
       f(1) = u(2)

       ! Momentum flux (equation 2).
       f(2) = c1 * (u(2) * (u(2) / u(1))) + c2 * u(3)

       ! Energy flux (equation 3).
       f(3) = (c3 * (u(2) * (u(2) / u(1))) + c4 * u(3)) * (u(2) / u(1))

     end function flux

end module problem_pde
