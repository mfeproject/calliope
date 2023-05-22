module burgers_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pde_class
  use mfe1_disc_core_type
  use parameter_list_type
  implicit none
  private

  type, extends(pde), public :: burgers
    private
    real(r8) :: visc
  contains
    procedure, nopass :: neqns
    procedure :: init
    procedure :: rhs
  end type

  public :: alloc_pde

contains

  pure integer function neqns()
    neqns = 1
  end function

  subroutine alloc_pde(cp) bind(c)
    use,intrinsic :: iso_c_binding
    type(c_ptr), intent(out) :: cp
    type(pde_box), pointer :: box
    allocate(box)
    allocate(burgers :: box%p)
    cp = c_loc(box)
  end subroutine

  subroutine init(this, disc, params, stat, errmsg)
    use parameter_list_type
    class(burgers), intent(out) :: this
    type(mfe1_disc_core), intent(in), target :: disc
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%disc => disc
    call params%get('visc', this%visc, stat=stat, errmsg=errmsg)
  end subroutine

  subroutine rhs(this, t)

    class(burgers), intent(inout) :: this
    real(r8), intent(in) :: t

    integer  :: j
    real(r8) :: c

    associate (rx1 => this%disc%r(2,1,:), rx2 => this%disc%r(2,2,:), &
               ru1 => this%disc%r(1,1,:), ru2 => this%disc%r(1,2,:), &
               u1 => this%disc%u(1,1,:), u2 => this%disc%u(1,2,:), &
        du => this%disc%du(1,:), n1 => this%disc%n(1,1,:), n2 => this%disc%n(2,1,:))
      do j = 1, this%disc%ncell
        c = -(du(j)/6.0_r8)*(2*u1(j)+u2(j))
        rx1(j) = c * n1(j)
        ru1(j) = c * n2(j)
        c = -(du(j)/6.0_r8)*(2*u2(j)+u1(j))
        rx2(j) = c * n1(j)
        ru2(j) = c * n2(j)
      end do
    end associate

    call this%disc%laplacian(eqno=1, coef=this%visc)

  end subroutine

end module burgers_type
