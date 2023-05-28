module cell_data_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: cell_data(npde)
    integer,len :: npde
    !integer :: npde, nvar
    real(r8) :: x(2), dx
    real(r8) :: u(npde,2), du(npde), dudx(npde), l(npde), nx(npde), nu(npde)
  contains
    !procedure :: init
    generic   :: update => update1, update2
    procedure, private :: update1, update2
    procedure :: set_val
  end type

contains

!  subroutine init(this, npde)
!    class(cell_data), intent(out) :: this
!    integer, intent(in) :: npde
!    this%npde = npde
!    this%nvar = npde + 1
!    allocate(this%u(npde,2))
!    allocate(this%du(npde), this%dudx(npde), this%l(npde), this%nx(npde), this%nu(npde))
!  end subroutine

  pure subroutine update1(this)
    class(cell_data(*)), intent(inout) :: this
    this%dx = this%x(2) - this%x(1)
    this%du = this%u(:,2) - this%u(:,1)
    this%dudx = this%du / this%dx
    this%l = sqrt(this%dx**2 + this%du**2)
    this%nx = -this%du/this%l
    this%nu =  this%dx/this%l
  end subroutine

  pure subroutine update2(this, y)
    class(cell_data(*)), intent(inout) :: this
    real(r8), intent(in) :: y(:,:)
    this%u = y(:this%npde,:)
    this%x = y(this%npde+1,:)
    call update1(this)
  end subroutine

  pure subroutine set_val(this, i, k, val, update)
    class(cell_data(*)), intent(inout) :: this
    integer, intent(in) :: i, k ! pde, vertex indices
    real(r8), intent(in) :: val
    logical, intent(in) :: update
    if (i == this%npde + 1) then
      this%x(k) = val
    else
      this%u(i,k) = val
    end if
    if (update) call update1(this)
  end subroutine

end module cell_data_type
