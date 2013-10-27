#include "f90_assert.fpp"

module nodal_FEM_1D

  use bdf2_dae, only: bdf2_model_evaluator 
  use kinds
  implicit none
  private
  
  type, extends(bdf2_model_evaluator), public :: ht_model
    private
    real(r8), public :: d
    real(r8), public :: atol, rtol

    integer :: n, prob = 0
    real(r8) :: u_left, u_right
    real(r8), allocatable, public :: x(:)
    real(r8), allocatable :: dx(:), jac(:,:)
  contains
    procedure :: size => model_size
    procedure :: init => init_simple
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
    procedure :: schk
    procedure :: udot
    procedure :: user
  end type
  
contains

  integer function model_size (this)
    class(ht_model), intent(in) :: this
    model_size = this%n
  end function model_size

  subroutine schk (this, u, stage, errc)
    class(ht_model) :: this
    real(r8), intent(in)  :: u(:)
    integer,  intent(in)  :: stage
    integer,  intent(out) :: errc
    errc = 0
  end subroutine schk

  subroutine init_simple (this, nnode, x0, x1, u0, u1, pnum)
  
    class(ht_model), intent(out) :: this
    integer,  intent(in) :: nnode   ! Number of mesh nodes.
    real(r8), intent(in) :: x0, x1  ! End points of the domain interval.
    real(r8), intent(in) :: u0, u1  ! Dirichlet boundary values.
    integer,  intent(in) :: pnum    ! Problem number (for PDE parameters).
    
    integer :: j
    real(r8) :: h
    
    ASSERT( nnode > 2 )
    ASSERT( x0 < x1 )
    
    this%n = nnode
    this%prob = pnum  ! set the problem number.
    
    this%u_left = u0
    this%u_right = u1
    
    allocate(this%x(nnode), this%dx(nnode-1), this%jac(3,nnode))
    
    !! Create the equi-spaced mesh.
    h = (x1 - x0) / real(nnode-1,r8)
    this%dx = h
    this%x(1) = x0
    do j = 2, nnode-1
      this%x(j) = this%x(j-1) + h
    end do
    this%x(nnode) = x1
    
  end subroutine init_simple

  subroutine eval_diff_coef (this, u, a)
  
    class(ht_model), intent(in) :: this
    real(r8), intent(in)  :: u(:)
    real(r8), intent(out) :: a(:)
    
    integer :: j
    
    ASSERT(size(u) == this%n)
    ASSERT(size(a) == this%n-1)
  
    select case (this%prob)
    case (1)  ! constant coefficient
      a = this%d
    case (2)
      do j = 1, this%n-1
        a(j) = this%d + max(0.0_r8, 0.5_r8*(u(j) + u(j+1)))
      end do
    case default
      INSIST( .false. )
    end select

  end subroutine eval_diff_coef

  subroutine eval_mass_matrix (dx, m)
  
    real(r8), intent(in)  :: dx(:)
    real(r8), intent(out) :: m(:,:)
    
    integer :: j
    
    ASSERT( size(m,1) == 3 )
    ASSERT( size(m,2) == size(dx)+1 )
    
    m(2,1) = 0.0_r8
    do j = 1, size(dx)
      m(2,j) = m(2,j) + dx(j)/3.0_r8
      m(3,j) = dx(j)/6.0_r8
      m(1,j+1) = dx(j)/6.0_r8
      m(2,j+1) = dx(j)/3.0_r8
    end do
    
  end subroutine eval_mass_matrix
  
  subroutine tdfactor (a)
  
    real(r8), intent(inout) :: a(:,:)
    
    integer :: j
    
    ASSERT( size(a,1) == 3 )
    
    do j = 2, size(a,2)
      a(3,j-1) = a(3,j-1)/a(2,j-1)
      a(2,j) = a(2,j) - a(1,j)*a(3,j-1)
    end do
    
  end subroutine tdfactor

  subroutine tdsolve (a, x)
  
    real(r8), intent(in) :: a(:,:)
    real(r8), intent(inout) :: x(:)
    
    integer :: j
    
    ASSERT( size(a,1) == 3 )
    ASSERT( size(a,2) == size(x) )
    
    !! Forward substitution.
    x(1) = x(1)/a(2,1)
    do j = 2, size(x)
      x(j) = (x(j) - a(1,j)*x(j-1))/a(2,j)
    end do
    
    !! Backward substitution.
    do j = size(x)-1, 1, -1
      x(j) = x(j) - a(3,j)*x(j+1)
    end do
    
  end subroutine tdsolve

  subroutine compute_f (this, t, u, udot, f)
  
    class(ht_model) :: this
    real(r8), intent(in)  :: t, u(:), udot(:)
    real(r8), intent(out) :: f(:)
    
    integer :: j, n
    real(r8) :: a(this%n-1)
    
    ASSERT( size(u) == this%n )
    ASSERT( size(udot) == this%n )
    ASSERT( size(f) == this%n )
    
    n = this%n
    call eval_diff_coef (this, u, a)
    f(2) = (this%dx(1)/3.0_r8)*udot(2) + (a(1)/this%dx(1))*(u(2) - this%u_left)
    do j = 2, n-2
      f(j)   = f(j) + (this%dx(j)/6.0_r8)*(2.0_r8*udot(j) + udot(j+1)) - (a(j)/this%dx(j))*(u(j+1) - u(j))
      f(j+1) =        (this%dx(j)/6.0_r8)*(udot(j) + 2.0_r8*udot(j+1)) + (a(j)/this%dx(j))*(u(j+1) - u(j))
    end do
    f(n-1) = f(n-1) + (this%dx(n-1)/3.0_r8)*udot(n-1) - (a(n-1)/this%dx(n-1))*(this%u_right - u(n-1))
    
    f(1) = u(1) - this%u_left
    f(n) = u(n) - this%u_right
    
  end subroutine compute_f
  
  subroutine apply_precon (this, t, u, f)
  
    class(ht_model) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), intent(inout) :: f(:)
    
    call tdsolve (this%jac, f)
    
  end subroutine apply_precon
  
  function udot (this, t, u)
  
    class(ht_model), intent(in) :: this
    real(r8), intent(in)  :: t, u(:)
    real(r8) :: udot(this%n)
    
    integer :: n, j
    real(r8) :: m(3,this%n), a(this%n-1)
    
    ASSERT( size(u) == this%n )
    
    n = this%n
    
    !! Evaluate the rhs of the linear system into UDOT.
    call eval_diff_coef (this, u, a)
    udot(2) = - (a(1)/this%dx(1))*(u(2) - this%u_left)
    do j = 2, n-2
      udot(j)   = udot(j) + (a(j)/this%dx(j))*(u(j+1) - u(j))
      udot(j+1) =         - (a(j)/this%dx(j))*(u(j+1) - u(j))
    end do
    udot(n-1) = udot(n-1) + (a(n-1)/this%dx(n-1))*(this%u_right - u(n-1))
    
    !! Evaluate the mass matrix.
    call eval_mass_matrix (this%dx, m)
    
    !! Solve for UDOT on the interior nodes only.
    call tdfactor (m(:,2:n-1))
    call tdsolve (m(:,2:n-1), udot(2:n-1))
    
    !! Time independent Dirichlet BV.
    udot(1) = 0.0_r8
    udot(n) = 0.0_r8
      
  end function udot
  
  subroutine compute_precon (this, t, u, dt)
  
    class(ht_model) :: this
    real(r8), intent(in) :: t, u(:), dt
    
    integer :: n, j
    real(r8) :: a(this%n-1), tmp
    
    n = this%n
    
    !! Jacobian of the linear term in udot.
    call eval_mass_matrix (this%dx, this%jac)
    this%jac = this%jac / dt
    
    !! Jacobian of the linear diffusion term in u.
    call eval_diff_coef (this, u, a)
    do j = 1, n-1
      tmp = a(j)/this%dx(j)
      this%jac(2,j) = this%jac(2,j) + tmp
      this%jac(3,j) = this%jac(3,j) - tmp
      this%jac(1,j+1) = this%jac(1,j+1) - tmp
      this%jac(2,j+1) = this%jac(2,j+1) + tmp
    end do
    
    !! Dirichlet BC at the left-most node.
    this%jac(2,1) = 1.0_r8
    this%jac(3,1) = 0.0_r8
    this%jac(1,2) = 0.0_r8
    
    !! Dirichlet BC at the right-most node.
    this%jac(2,n) = 1.0_r8
    this%jac(1,n) = 0.0_r8
    this%jac(3,n-1) = 0.0_r8
    
    call tdfactor (this%jac) 

  end subroutine compute_precon
  
  subroutine user (this, t, u)
    class(ht_model), intent(in) :: this
    real(r8), intent(in) :: t, u(:)
    integer :: j
    print *, 'T=', t
    print '(2es13.5)', (this%x(j), u(j), j=1,this%n)
  end subroutine user
  
  subroutine du_norm (this, u, du, error)
    class(ht_model) :: this
    real(r8), intent(in) :: u(:), du(:)
    real(r8), intent(out) :: error
    error = maxval(abs(du)/(this%atol + this%rtol*abs(u)))
  end subroutine du_norm
  
end module nodal_FEM_1D
