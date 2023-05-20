module mfe_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_env_type
  use mfe_model_type
  use btd_matrix_type
  use mfe1_vector_type
  implicit none
  private

  type, public :: mfe_precon
    private
    type(mfe_env), pointer :: env => null() ! reference only
    type(mfe_model), pointer :: model => null() ! reference only
    type(btd_matrix) :: jac
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type

contains

  subroutine init(this, env, model)
    class(mfe_precon), intent(out) :: this
    type(mfe_env), intent(in), target :: env
    type(mfe_model), intent(in), target :: model
    this%env => env
    this%model => model
    call this%jac%init(this%model%nvars, this%model%nnode)
  end subroutine

  subroutine compute(this, u, udot, t, dt, stat)
    class(mfe_precon), intent(inout) :: this
    type(mfe1_vector), intent(in) :: u, udot
    real(r8), intent(in) :: t, dt
    integer, intent(out) :: stat
    call eval_jacobian(this, u, udot, t, dt, stat)
  end subroutine

  subroutine apply(this, u)
    class(mfe_precon), intent(in) :: this
    type(mfe1_vector), intent(inout) :: u
    call this%jac%solve(u%array)
  end subroutine

  subroutine eval_jacobian(this, u, udot, t, h, stat)

    use block_linear_solver

    type(mfe_precon), intent(inout) :: this
    type(mfe1_vector), intent(in) :: u, udot
    real(r8), intent(in) :: t, h
    integer, intent(out) :: stat

    integer :: i, j, n
    real(r8) :: diag(u%neqns+1,u%neqns+1,u%nnode)
    character(:), allocatable :: errmsg

    call this%env%timer%start('preprocessing')
    call this%model%disc%update(u, udot)
    call this%env%timer%stop('preprocessing')

    call this%env%timer%start('mass-matrix')
    call this%model%disc%eval_mass_matrix(factor = -1.0_r8 / h)
    call this%env%timer%stop('mass-matrix')

    ! Capture the unscaled diagonal for preconditioning.
    call this%env%timer%start('diag-pc')
    call this%model%disc%assemble_diagonal(diag)
    diag = (-h) * diag
    do i = 1, size(this%model%bc_dir)
      do j = 1, size(this%model%bc_dir(i)%index)
        n = this%model%bc_dir(i)%index(j)
        diag(:,i,n) = 0.0_r8
        diag(i,:,n) = 0.0_r8
        diag(i,i,n) = 1.0_r8
      end do
    end do
    call this%env%timer%stop('diag-pc')

    call this%env%timer%start('eval_dfdy')
    call this%model%disc%eval_dfdy(t, stat, errmsg)
    call this%env%timer%stop('eval_dfdy')

    if (stat /= 0) then
      call this%env%log%info('eval_dfdy: ' // errmsg)
      stat = 1
      return
    end if

    call this%env%timer%start('assembly')
    call this%model%disc%assemble_matrix(this%jac)
    do i = 1, size(this%model%bc_dir)
      do j = 1, size(this%model%bc_dir(i)%index)
        n = this%model%bc_dir(i)%index(j)
        call this%jac%set_dir_var(i, n)
      end do
    end do
    call this%env%timer%stop('assembly')

    ! Diagonal preconditioning.
    call this%env%timer%start('diag-pc')
    call vfct(diag)
    call vmslv(diag, this%jac%l)
    call vmslv(diag, this%jac%d)
    call vmslv(diag, this%jac%u)
    call this%env%timer%stop('diag-pc')

    ! Factorize the Jacobian.
    call this%env%timer%start('factorization')
    call this%jac%factor
    call this%env%timer%stop('factorization')

  end subroutine eval_jacobian

end module mfe_precon_type
