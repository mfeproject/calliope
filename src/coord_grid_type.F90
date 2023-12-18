!!
!! COORD_GRID_TYPE
!!

module coord_grid_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: coord_grid
    private
    real(r8), allocatable :: grid(:,:)
    integer,  allocatable :: nint(:)
    real(r8), allocatable :: ratio(:)
  contains
    procedure :: init
    procedure :: get_grid
  end type

contains

  subroutine init(this, params, grid_key, nint_key, ratio_key, stat, errmsg)

    use parameter_list_type
    use string_utilities, only: i_to_c

    class(coord_grid), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: grid_key, nint_key, ratio_key
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: rarray(:)
    integer,  allocatable :: iarray(:)
    integer :: n

    !! The coarse grid points (required)
    if (.not.params%is_parameter(grid_key)) then
      stat = 1
      errmsg = 'no such parameter: "' // grid_key // '"'
      return
    end if
    if (params%is_vector(grid_key)) then
      call params%get(grid_key, rarray, stat, errmsg)
      if (stat /= 0) return
      allocate(this%grid(1,size(rarray)))
      this%grid(1,:) = rarray
    else if (params%is_matrix(grid_key)) then
      call params%get(grid_key, this%grid, stat, errmsg)
      if (stat /= 0) return
    else
      stat = 1
      errmsg = 'not a vector or matrix parameter: "' // grid_key // '"'
      return
    end if
    n = size(this%grid,dim=2)
    if (n < 2) then
      stat = 1
      errmsg = '"' // grid_key // '" requires at least two values'
      return
    end if

    n = n - 1 ! expected size of the following arrays

    !! The number of subintervals for each coarse grid segment (required)
    call params%get(nint_key, iarray, stat, errmsg)
    if (stat /= 0) return
    if (size(iarray) /= n) then
      stat = 1
      errmsg = i_to_c(n) // ' values required for "' // nint_key // '"'
      return
    else if (any(iarray < 1)) then
      stat = 1
      errmsg = '"' // nint_key // '" values must be > 0'
      return
    end if
    call move_alloc(iarray, this%nint)

    !! The ratio for biased subdivision for each segment (optional)
    if (params%is_parameter(ratio_key)) then
      call params%get(ratio_key, rarray, stat, errmsg)
      if (stat /= 0) return
      if (size(rarray) /= n) then
        stat = 1
        errmsg = i_to_c(n) // ' values required for "' // ratio_key // '"'
        return
      else if (any(rarray <= 0)) then
        stat = 1
        errmsg = '"' // ratio_key // '" values must be > 0'
        return
      end if
      call move_alloc(rarray, this%ratio)
    else
      allocate(this%ratio(n))
      this%ratio = 1.0_r8
    end if

  end subroutine init

  !! Using the components GRID(1:N+1), NINT(1:N), and RATIO(1:N) this method
  !! partitions the interval [GRID(1), GRID(N+1)], returning the results in
  !! the array X(0:M), where M = SUM(NINT). The interval [GRID(j), GRID(j+1)]
  !! is partitioned into NINT(j) subintervals, with the ratio of the lengths
  !! of successive subintervals equal to RATIO(j).

  subroutine get_grid(this, x)

    class(coord_grid), intent(in) :: this
    real(r8), allocatable, intent(out) :: x(..)

    integer :: i, j, n1, n2

    select rank (x)
    rank (1)

      allocate(x(1+sum(this%nint)))
      n2 = 1
      x(1) = this%grid(1,1)
      do j = 1, size(this%nint)
        n1 = n2
        n2 = n1 + this%nint(j)
        x(n2) = this%grid(1,j+1)
        call partition_interval(this%ratio(j), x(n1:n2))
      end do

    rank (2)

      allocate(x(size(this%grid,dim=1),1+sum(this%nint)))
      n2 = 1
      x(:,1) = this%grid(:,1)
      do j = 1, size(this%nint)
        n1 = n2
        n2 = n1 + this%nint(j)
        x(:,n2) = this%grid(:,j+1)
        do i = 1, size(x,dim=1)
          call partition_interval(this%ratio(j), x(i,n1:n2))
        end do
      end do

    rank default

      error stop "coord_grid%get_grid: invalid rank for x"

    end select

  end subroutine get_grid

  !! Given an array X(0:N) this auxiliary subroutine partitions the interval
  !! [X(0), X(N)], defining the intermediate values X(j) for j = 1 to N-1.
  !! The ratio of the lengths of successive subintervals is RATIO.

  subroutine partition_interval(ratio, x)
    real(r8), intent(in) :: ratio
    real(r8), intent(inout) :: x(0:)
    integer :: j, n
    real(r8) :: a, dx
    n = ubound(x,1)
    a = 1
    do j = 1, n-1
      a = ratio*a + 1
    end do
    dx = (x(n)-x(0))/a
    a = 1
    do j = 1, n-1
      x(j) = x(0) + a*dx
      a = ratio*a + 1
    end do
  end subroutine

end module coord_grid_type
