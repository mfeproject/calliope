!!
!! MFE1_VECTOR_TYPE
!!
!! An implementation of the VECTOR base class that stores the 1D MFE unknowns
!! as a rank-2 array in node major order.
!!

module mfe1_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  !use timer_tree_type
  implicit none
  private

  type, extends(vector), public :: mfe1_vector
    integer :: neqns, nnode
    real(r8), allocatable :: array(:,:)
  contains
    !! Deferred base class procedures
    procedure :: clone1
    procedure :: clone2
    procedure :: copy_
    procedure :: setval
    procedure :: scale
    procedure :: update1_
    procedure :: update2_
    procedure :: update3_
    procedure :: update4_
    procedure :: dot_
    procedure :: norm2 => norm2_
    procedure :: checksum
    !! Additional procedures specific to this type
    generic :: init => init_dim, init_mold
    procedure, private :: init_dim, init_mold
  end type

contains

  subroutine init_dim(this, neqns, nnode)
    class(mfe1_vector), intent(out) :: this
    integer, intent(in) :: neqns, nnode
    this%neqns = neqns
    this%nnode = nnode
    allocate(this%array(1+neqns,nnode))
  end subroutine

  subroutine init_mold(this, mold)
    class(mfe1_vector), intent(out) :: this
    class(mfe1_vector), intent(in)  :: mold
    this%neqns = mold%neqns
    this%nnode = mold%nnode
    allocate(this%array, mold=mold%array)
  end subroutine

  subroutine clone1(this, clone)
    class(mfe1_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    allocate(clone, source=this) ! easy, but with an undesired copy of array data too
  end subroutine

  subroutine clone2(this, clone, n)
    class(mfe1_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    allocate(clone(n), source=this) ! easy, but with an undesired copy of array data too
  end subroutine

  subroutine copy_(dest, src)
    class(mfe1_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    !call start_timer('vector%copy_')
    select type (src)
    type is (mfe1_vector)
      !dest%array(:,:) = src%array
      call cpy(dest%array, src%array, size(dest%array))
    end select
    !call stop_timer('vector%copy_')
  contains
    pure subroutine cpy(lhs, rhs, n)
      integer, intent(in) :: n
      real(r8), intent(out) :: lhs(n)
      real(r8), intent(in)  :: rhs(n)
      integer :: j
      do j = 1, n
        lhs(j) = rhs(j)
      end do
    end subroutine
  end subroutine

  subroutine setval(this, val)
    class(mfe1_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    !call start_timer('vector%setval')
    !this%array = val
    call cpy(this%array, val, size(this%array))
    !call stop_timer('vector%setval')
  contains
    pure subroutine cpy(lhs, val, n)
      integer, intent(in) :: n
      real(r8), intent(out) :: lhs(n)
      real(r8), intent(in) :: val
      integer :: j
      do j = 1, n
        lhs(j) = val
      end do
    end subroutine
  end subroutine

  subroutine scale(this, a)
    class(mfe1_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    !call start_timer('vector%scale')
    !this%array = a*this%array
    call scl(this%array, a, size(this%array))
    !call stop_timer('vector%scale')
  contains
    pure subroutine scl(lhs, a, n)
      integer, intent(in) :: n
      real(r8), intent(inout) :: lhs(n)
      real(r8), intent(in) :: a
      integer :: j
      do j = 1, n
        lhs(j) = a*lhs(j)
      end do
    end subroutine
  end subroutine

  !! Conventional SAXPY procedure: y <-- a*x + y
  subroutine update1_(this, a, x)
    class(mfe1_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    !call start_timer('vector%update1_')
    select type (x)
    type is (mfe1_vector)
      !this%array = a * x%array + this%array
      call aux(this%array, a, x%array, size(this%array))
    end select
    !call stop_timer('vector%update1_')
  contains
    pure subroutine aux(y, a, x, n)
      integer, intent(in) :: n
      real(r8), intent(inout) :: y(n)
      real(r8), intent(in) :: a, x(n)
      integer :: j
      do j = 1, n
        y(j) = y(j) + a*x(j)
      end do
    end subroutine
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  subroutine update2_(this, a, x, b)
    class(mfe1_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    !call start_timer('vector%update2_')
    select type (x)
    class is (mfe1_vector)
      !this%array = a * x%array + b * this%array
      call aux(this%array, a, x%array, b, size(this%array))
    end select
    !call stop_timer('vector%update2_')
  contains
    pure subroutine aux(y, a, x, b, n)
      integer, intent(in) :: n
      real(r8), intent(inout) :: y(n)
      real(r8), intent(in) :: a, x(n), b
      integer :: j
      do j = 1, n
        y(j) = b*y(j) + a*x(j)
      end do
    end subroutine
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  subroutine update3_(this, a, x, b, y)
    class(mfe1_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    !call start_timer('vector%update3_')
    select type (x)
    class is (mfe1_vector)
      select type (y)
      class is (mfe1_vector)
        !this%array = a * x%array + b * y%array + this%array
        call aux(this%array, a, x%array, b, y%array, size(this%array))
      end select
    end select
    !call stop_timer('vector%update3_')
  contains
    pure subroutine aux(z, a, x, b, y, n)
      integer, intent(in) :: n
      real(r8), intent(inout) :: z(n)
      real(r8), intent(in) :: a, x(n), b, y(n)
      integer :: j
      do j = 1, n
        z(j) = z(j) + a*x(j) + b*y(j)
      end do
    end subroutine
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  subroutine update4_(this, a, x, b, y, c)
    class(mfe1_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    !call start_timer('vector%update4_')
    select type (x)
    class is (mfe1_vector)
      select type (y)
      class is (mfe1_vector)
        !this%array = a * x%array + b * y%array + c * this%array
        call aux(this%array, a, x%array, b, y%array, c, size(this%array))
      end select
    end select
    !call stop_timer('vector%update4_')
  contains
    pure subroutine aux(z, a, x, b, y, c, n)
      integer, intent(in) :: n
      real(r8), intent(inout) :: z(n)
      real(r8), intent(in) :: a, x(n), b, y(n), c
      integer :: j
      do j = 1, n
        z(j) = c*z(j) + a*x(j) + b*y(j)
      end do
    end subroutine
  end subroutine

  function dot_(x, y) result(dp)
    class(mfe1_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    !call start_timer('vector%dot_')
    select type (y)
    class is (mfe1_vector)
      !dp = 0.0_r8
      !do j = 1, size(x%array,dim=2)
      !  do i = 1, size(x%array,dim=1)
      !    dp = dp + x%array(i,j) * y%array(i,j)
      !  end do
      !end do
      dp = aux(x%array, y%array, size(x%array))
    end select
    !call stop_timer('vector%dot_')
  contains
    pure function aux(x, y, n) result(xdoty)
      integer, intent(in) :: n
      real(r8), intent(in) :: x(n), y(n)
      real(r8) :: xdoty
      integer :: j
      xdoty = 0.0_r8
      do j = 1, n
        xdoty = xdoty + x(j)*y(j)
      end do
      !xdoty = dot_product(x, y)
    end function
  end function

  function norm2_(this)
    class(mfe1_vector), intent(in) :: this
    real(r8) :: norm2_
    !integer :: i, j
    !call start_timer('vector%norm2')
    norm2_ = norm2(this%array)
    !call stop_timer('vector%norm2')
    !norm2_ = 0.0_r8
    !do j = 1, this%ny
    !  do i = 1, this%nx
    !    norm2_ = norm2_ + this%array(i,j)**2
    !  end do
    !end do
    !norm2_ = sqrt(norm2_)
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(mfe1_vector), intent(in) :: this
    logical, intent(in), optional :: full ! default is FALSE
    character(:), allocatable :: string
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full ! no distinction for now
    call hash%update(this%array)
    string = hash%hexdigest()
  end function

end module mfe1_vector_type
