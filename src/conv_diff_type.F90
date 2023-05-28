module conv_diff_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pde_class
  use parameter_list_type
  implicit none
  private

  type, extends(pde), public :: conv_diff
    private
    real(r8) :: lapl_coef(2)
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
    allocate(conv_diff(1) :: box%p)
    cp = c_loc(box)
  end subroutine

  subroutine init(this, eqw, params, stat, errmsg)
    use parameter_list_type
    class(conv_diff(*)), intent(out) :: this
    real(r8), intent(in) :: eqw(:)
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8) :: visc
    this%eqw = eqw(1:1) ! unused
    call params%get('visc', visc, stat=stat, errmsg=errmsg)
    this%lapl_coef = [visc, visc]
  end subroutine

  pure subroutine rhs(this, t, cdata, gx, gu)

    use cell_data_type
    use pde_utilities, only: add_lapl

    class(conv_diff(*)), intent(inout) :: this
    real(r8), intent(in) :: t
    type(cell_data), intent(in) :: cdata
    real(r8), intent(out) :: gx(:), gu(:,:)

    integer  :: j
    real(r8) :: c

    c = -0.5_r8*cdata%du(1)
    gx(1)   = c * cdata%nx(1)
    gu(1,1) = c * cdata%nu(1)
    gx(2)   = c * cdata%nx(1)
    gu(1,2) = c * cdata%nu(1)

    call add_lapl(cdata, 1, this%lapl_coef, gx, gu)

  end subroutine

end module conv_diff_type
