#include "f90_assert.fpp"

module navier_stokes_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pde_class
  use mfe1_disc_core_type
  use parameter_list_type
  implicit none
  private

  type, extends(pde), public :: navier_stokes
    type(mfe1_disc_core), pointer :: disc => null() ! reference only
    real(r8) :: visc
  contains
    procedure, nopass :: neqns
    procedure :: init
    procedure :: rhs
  end type

  public :: alloc_pde

contains

  pure integer function neqns()
    neqns = 3
  end function

  subroutine alloc_pde(cp) bind(c)
    use,intrinsic :: iso_c_binding
    type(c_ptr), intent(out) :: cp
    type(pde_box), pointer :: box
    allocate(box)
    allocate(navier_stokes :: box%p)
    cp = c_loc(box)
  end subroutine

  subroutine init(this, disc, params, stat, errmsg)
    use parameter_list_type
    class(navier_stokes), intent(out) :: this
    type(mfe1_disc_core), intent(in), target :: disc
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%disc => disc
    call params%get('visc', this%visc, stat=stat, errmsg=errmsg)
  end subroutine

  subroutine rhs(this, t)

    class(navier_stokes), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, k
    real(r8) :: rx1, rx2, term
    real(r8) :: umid(3), f1(3), f2(3), favg(3)
    real(r8), parameter :: c1 = 1.0_r8 / 6.0_r8, c2 = 4.0_r8 / 6.0_r8

    associate (disc => this%disc)

    do i = 1, disc%ncell

      f1 = flux(disc%u(:3,1,i))                      ! Flux at left endpoint.
      f2 = flux(disc%u(:3,2,i))                      ! Flux at right endpoint.
      umid = 0.5_r8 * (disc%u(:3,1,i) + disc%u(:3,2,i))   ! Variables at midpoint.
      favg = c1 * (f1 + f2) + c2 * flux(umid)   ! Average flux (Simpson).

      rx1 = 0.0_r8
      rx2 = 0.0_r8

      do k = 1, 3

        term = (disc%eqw(k) * (favg(k) - f1(k)))
        rx1 = rx1 - term * disc%n(1,k,i)
        disc%r(k,1,i) = - term * disc%n(2,k,i)

        term = (disc%eqw(k) * (f2(k) - favg(k)))
        rx2 = rx2 - term * disc%n(1,k,i)
        disc%r(k,2,i) = - term * disc%n(2,k,i)

      end do

      disc%r(4,1,i) = rx1
      disc%r(4,2,i) = rx2

    end do

    call disc%laplacian(eqno=1, coef=this%visc)
    call disc%laplacian(eqno=2, coef=this%visc)
    call disc%laplacian(eqno=3, coef=this%visc)

    end associate

  contains

    pure function flux(u) result(f)

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

    end function

  end subroutine rhs

end module
