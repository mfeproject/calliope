#include "f90_assert.fpp"

module btd_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use block_linear_solver
  implicit none
  private

  type, public :: btd_matrix
    real(r8), allocatable :: l(:,:,:), d(:,:,:), u(:,:,:)
  contains
    procedure :: init
    procedure :: factor
    procedure :: solve
    procedure :: set_dir_var
  end type

contains

  subroutine init(this, m, n)
    class(btd_matrix), intent(out) :: this
    integer, intent(in) :: m, n
    allocate(this%l(m,m,n))
    allocate(this%d, this%u, mold=this%l)
  end subroutine

  subroutine factor(this)

    class(btd_matrix), intent(inout) :: this

    integer :: i

    call fct(this%d(:,:,1))
    do i = 2, size(this%d,dim=3)
      call mslv(this%d(:,:,i-1), this%u(:,:,i-1))
      call cmab(this%d(:,:,i), this%l(:,:,i), this%u(:,:,i-1))
      call fct(this%d(:,:,i))
    end do

  end subroutine


  subroutine solve(this, b)

    class(btd_matrix), intent(in) :: this
    real(r8), intent(inout) :: b(:,:)

    integer :: i

    ASSERT(size(b,1) == size(this%d,1))
    ASSERT(size(b,2) == size(this%d,3))

    !! Forward substitution
    call slv(this%d(:,:,1), b(:,1))
    do i = 2, size(b,2)
      call ymax(b(:,i), this%l(:,:,i), b(:,i-1))
      call slv(this%d(:,:,i), b(:,i))
    end do

    !! Backward substitution
    do i = size(b,2)-1, 1, -1
      call ymax(b(:,i), this%u(:,:,i), b(:,i+1))
    end do

  end subroutine

  subroutine set_dir_var(this, i, j)
    class(btd_matrix), intent(inout) :: this
    integer, intent(in) :: i, j
    this%l(i,:,j) = 0.0_r8
    this%d(i,:,j) = 0.0_r8
    this%u(i,:,j) = 0.0_r8
    this%d(i,i,j) = 1.0_r8
  end subroutine

end module btd_matrix_type
