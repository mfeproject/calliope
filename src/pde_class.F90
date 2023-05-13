module pde_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe1_disc_core_type
  use parameter_list_type
  implicit none
  private
  
  type, abstract, public :: pde
  contains
    procedure(neqns), nopass, deferred :: neqns
    procedure(init), deferred :: init
    procedure(rhs), deferred :: rhs
  end type
  
  abstract interface
    pure integer function neqns()
    end function
    subroutine init(this, disc, params, stat, errmsg)
      import pde, mfe1_disc_core, parameter_list
      class(pde), intent(out) :: this
      type(mfe1_disc_core), intent(in), target :: disc
      type(parameter_list), intent(inout) :: params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
    subroutine rhs(this, t)
      import pde, r8
      class(pde), intent(inout) :: this
      real(r8), intent(in) :: t
    end subroutine
  end interface

  type, public :: pde_box
    class(pde), allocatable :: p
  end type

end module
