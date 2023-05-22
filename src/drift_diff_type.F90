#include "f90_assert.fpp"

!! Single carrier semiconductor equations where the hole density is given by
!! p = exp(USCF*U) and the voltage is v = VSCF*V. The equations for U and V are
!!
!!    dU/dt = Grad(U)*(VSCF*Grad(V) + LAMBDA*USCF*Grad(U))
!!               + (1 - exp(USCF*U))/USCF + LAMBDA*Lapl(U)
!!
!!    dV/dt = (1/EPS)*(Lapl(V) - (1 - exp(USCF*U))/VSCF)
!!
!!   LAMBDA and EPS are parameters and USCF and VSCF are scaling factors.

module drift_diff_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pde_class
  use mfe1_disc_core_type
  use parameter_list_type
  implicit none
  private

  type, extends(pde), public :: drift_diff
    private
    real(r8) :: lambda, eps, uscf, vscf
  contains
    procedure, nopass :: neqns
    procedure :: init
    procedure :: rhs
  end type

  public :: alloc_pde

contains

  pure integer function neqns()
    neqns = 2
  end function

  subroutine alloc_pde(cp) bind(c)
    use,intrinsic :: iso_c_binding
    type(c_ptr), intent(out) :: cp
    type(pde_box), pointer :: box
    allocate(box)
    allocate(drift_diff :: box%p)
    cp = c_loc(box)
  end subroutine

  subroutine init(this, disc, params, stat, errmsg)
    use parameter_list_type
    class(drift_diff), intent(out) :: this
    type(mfe1_disc_core), intent(in), target :: disc
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%disc => disc
    call params%get('lambda', this%lambda, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('eps', this%eps, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('u-scale-factor', this%uscf, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('v-scale-factor', this%vscf, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
  end subroutine

  !! Copied from gwmfe1ds with minimal changes

  subroutine rhs(this, t)

    class(drift_diff), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: j
    real(r8) :: c1, c2, c3, c4, pa0, pa1, term1, term2
    
    real(r8), parameter :: g3w(3) = [5.0_r8, 8.0_r8, 5.0_r8] / 18.0_r8
    real(r8), parameter :: a = (1-sqrt(0.6_r8))/2, b = 0.5_r8, c = (1+sqrt(0.6_r8))/2
    real(r8), parameter :: g3a0(3) = [a,b,c], g3a1(3) = [c,b,a]

    c1 = 0.5_r8*this%disc%eqw(1)*this%vscf
    c2 = 0.5_r8*this%disc%eqw(1)*this%uscf*this%lambda
    c3 = this%disc%eqw(1)/this%uscf
    c4 = -this%disc%eqw(2)/(this%eps*this%vscf)

    associate (u0 => this%disc%u(:,1,:), u1 => this%disc%u(:,2,:), &
               gx0 => this%disc%r(3,1,:), gx1 => this%disc%r(3,2,:), &
               gu0 => this%disc%r(1:2,1,:), gu1 => this%disc%r(1:2,2,:), &
               du => this%disc%du, dx => this%disc%dx, &
               n1 => this%disc%n(1,:,:), n2 => this%disc%n(2,:,:))

      do j = 1, this%disc%ncell

         pa0 = g3w(1)*g3a0(1)*exp(this%uscf*(g3a0(1)*u0(1,j)+g3a1(1)*u1(1,j))) &
             + g3w(2)*g3a0(2)*exp(this%uscf*(g3a0(2)*u0(1,j)+g3a1(2)*u1(1,j))) &
             + g3w(3)*g3a0(3)*exp(this%uscf*(g3a0(3)*u0(1,j)+g3a1(3)*u1(1,j)))

         pa1 = g3w(1)*g3a1(1)*exp(this%uscf*(g3a0(1)*u0(1,j)+g3a1(1)*u1(1,j))) &
             + g3w(2)*g3a1(2)*exp(this%uscf*(g3a0(2)*u0(1,j)+g3a1(2)*u1(1,j))) &
             + g3w(3)*g3a1(3)*exp(this%uscf*(g3a0(3)*u0(1,j)+g3a1(3)*u1(1,j)))

         term1 = du(1,j)*(c1*du(2,j) + c2*du(1,j))/dx(j)
         term2 = dx(j)*(.5 - pa0)
         gx0(j) = (term1 + c3*term2)*n1(1,j) + (c4*term2)*n1(2,j)
         gu0(1,j) = (term1 + c3*term2)*n2(1,j)
         gu0(2,j) = (c4*term2)*n2(2,j)


         term2 = dx(j)*(.5 - pa1)
         gx1(j) = (term1 + c3*term2)*n1(1,j) + (c4*term2)*n1(2,j)
         gu1(1,j) = (term1 + c3*term2)*n2(1,j)
         gu1(2,j) = (c4*term2)*n2(2,j)

      end do
    end associate

    call this%disc%laplacian(eqno=1, coef=this%lambda)
    call this%disc%laplacian(eqno=2, coef=1.0_r8/this%eps)

  end subroutine rhs

end module
