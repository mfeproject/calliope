module pde_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use cell_data_type
  use parameter_list_type
  implicit none
  private

  type, abstract, public :: pde(npde)
    integer,len :: npde
    real(r8), allocatable :: eqw(:)
  contains
    !procedure(neqns), nopass, deferred :: neqns
    procedure(init), deferred :: init
    procedure(rhs), deferred :: rhs
  end type

  abstract interface
    !pure integer function neqns()
    !end function
    subroutine init(this, eqw, params, stat, errmsg)
      import pde, parameter_list, r8
      class(pde(*)), intent(out) :: this
      real(r8), intent(in) :: eqw(:)
      type(parameter_list), intent(inout) :: params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
    pure subroutine rhs(this, t, cdata, gx, gu)
      import pde, cell_data, r8
      class(pde(*)), intent(inout) :: this
      real(r8), intent(in) :: t
      type(cell_data(*)), intent(in) :: cdata
      real(r8), intent(out) :: gx(:), gu(:,:)
    end subroutine
  end interface

  type, public :: pde_box
    class(pde(:)), allocatable :: p
  end type

end module
