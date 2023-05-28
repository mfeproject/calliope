#include "f90_assert.fpp"

module navier_stokes_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pde_class
  use parameter_list_type
  implicit none
  private

  type, extends(pde), public :: navier_stokes
    private
    real(r8) :: lapl_coef(2,3)
  contains
    procedure :: init
    procedure :: rhs
  end type

  public :: alloc_pde

contains

  subroutine alloc_pde(cp) bind(c)
    use,intrinsic :: iso_c_binding
    type(c_ptr), intent(out) :: cp
    type(pde_box), pointer :: box
    allocate(box)
    allocate(navier_stokes(npde=3) :: box%p)
    cp = c_loc(box)
  end subroutine

  subroutine init(this, eqw, params, stat, errmsg)
    use parameter_list_type
    class(navier_stokes(*)), intent(out) :: this
    real(r8), intent(in) :: eqw(:)
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8) :: visc
    this%eqw = eqw(1:3)
    call params%get('visc', visc, stat=stat, errmsg=errmsg)
    this%lapl_coef = spread(this%eqw*visc, dim=1, ncopies=2)
  end subroutine

  pure subroutine rhs(this, t, cdata, gx, gu)

    use cell_data_type
    use pde_utilities, only: add_lapl

    class(navier_stokes(*)), intent(inout) :: this
    real(r8), intent(in) :: t
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(out) :: gx(:), gu(:,:)

    integer :: k
    real(r8) :: term
    real(r8) :: umid(3), f1(3), f2(3), favg(3)
    real(r8), parameter :: c1 = 1.0_r8 / 6.0_r8, c2 = 4.0_r8 / 6.0_r8

    f1 = flux(cdata%u(:,1)) ! flux at left endpoint
    f2 = flux(cdata%u(:,2)) ! flux at right endpoint
    umid = 0.5_r8 * (cdata%u(:,1) + cdata%u(:,2)) ! unknowns at midpoint
    favg = c1*(f1 + f2) + c2*flux(umid) ! average flux (Simpson)

    gx(1) = 0.0_r8
    gx(2) = 0.0_r8

    do k = 1, 3
      term = this%eqw(k)*(favg(k) - f1(k))
      gx(1)   = gx(1) - term*cdata%nx(k)
      gu(k,1) =       - term*cdata%nu(k)

      term = this%eqw(k)*(f2(k) - favg(k))
      gx(2)   = gx(2) - term*cdata%nx(k)
      gu(k,2) =       - term*cdata%nu(k)
    end do

    do k = 1, 3
      call add_lapl(cdata, k, this%lapl_coef(:,k), gx, gu)
    end do

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
