module index_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: index_func
    integer,  allocatable :: index(:)
    real(r8), allocatable :: value(:)
  contains
    procedure :: compute
  end type

contains

  subroutine compute(this, t)
    class(index_func), intent(inout) :: this
    real(r8), intent(in) :: t
    ! a no-op for the time being
  end subroutine

end module index_func_type
