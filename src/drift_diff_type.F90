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
  use parameter_list_type
  implicit none
  private

  type, extends(pde), public :: drift_diff
    private
    real(r8) :: lambda, eps, uscf, vscf
    real(r8) :: lapl_coef(2,2)
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

  subroutine init(this, eqw, params, stat, errmsg)
    use parameter_list_type
    class(drift_diff), intent(out) :: this
    real(r8), intent(in) :: eqw(:)
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%eqw = eqw(1:2)
    call params%get('u-scale-factor', this%uscf, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('v-scale-factor', this%vscf, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('lambda', this%lambda, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('eps', this%eps, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    this%lapl_coef = spread([this%eqw(1)*this%lambda, this%eqw(2)/this%eps], dim=1, ncopies=2)
  end subroutine

  !! Copied from gwmfe1ds with minimal changes

  pure subroutine rhs(this, t, cdata, gx, gu)

    use cell_data_type
    use pde_utilities, only: add_lapl

    class(drift_diff), intent(inout) :: this
    real(r8), intent(in) :: t
    type(cell_data), intent(in) :: cdata
    real(r8), intent(out) :: gx(:), gu(:,:)

    integer :: j
    real(r8) :: c1, c2, c3, c4, pa0, pa1, term1, term2

    real(r8), parameter :: g3w(3) = [5.0_r8, 8.0_r8, 5.0_r8] / 18.0_r8
    real(r8), parameter :: a = (1-sqrt(0.6_r8))/2, b = 0.5_r8, c = (1+sqrt(0.6_r8))/2
    real(r8), parameter :: g3a0(3) = [a,b,c], g3a1(3) = [c,b,a]

    c1 = 0.5_r8*this%eqw(1)*this%vscf
    c2 = 0.5_r8*this%eqw(1)*this%uscf*this%lambda
    c3 = this%eqw(1)/this%uscf
    c4 = -this%eqw(2)/(this%eps*this%vscf)

    associate (u0 => cdata%u(:,1), u1 => cdata%u(:,2), &
               du => cdata%du, dx => cdata%dx, nx => cdata%nx, nu => cdata%nu)

      pa0 = g3w(1)*g3a0(1)*exp(this%uscf*(g3a0(1)*u0(1)+g3a1(1)*u1(1))) &
          + g3w(2)*g3a0(2)*exp(this%uscf*(g3a0(2)*u0(1)+g3a1(2)*u1(1))) &
          + g3w(3)*g3a0(3)*exp(this%uscf*(g3a0(3)*u0(1)+g3a1(3)*u1(1)))

      pa1 = g3w(1)*g3a1(1)*exp(this%uscf*(g3a0(1)*u0(1)+g3a1(1)*u1(1))) &
          + g3w(2)*g3a1(2)*exp(this%uscf*(g3a0(2)*u0(1)+g3a1(2)*u1(1))) &
          + g3w(3)*g3a1(3)*exp(this%uscf*(g3a0(3)*u0(1)+g3a1(3)*u1(1)))

      term1 = du(1)*(c1*du(2) + c2*du(1))/dx
      term2 = dx*(0.5_r8 - pa0)
      gx(1) = (term1 + c3*term2)*nx(1) + (c4*term2)*nx(2)
      gu(1,1) = (term1 + c3*term2)*nu(1)
      gu(2,1) = (c4*term2)*nu(2)


      term2 = dx*(0.5_r8 - pa1)
      gx(2) = (term1 + c3*term2)*nx(1) + (c4*term2)*nx(2)
      gu(1,2) = (term1 + c3*term2)*nu(1)
      gu(2,2) = (c4*term2)*nu(2)
    end associate

    call add_lapl(cdata, 1, this%lapl_coef(:,1), gx, gu)
    call add_lapl(cdata, 2, this%lapl_coef(:,2), gx, gu)

  end subroutine rhs

end module
