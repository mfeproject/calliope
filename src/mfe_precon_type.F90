module mfe_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_env_type
  use mfe_model_type
  use btd_matrix_type
  use mfe1_vector_type
  implicit none
  private

  type, public :: mfe_precon
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
    call this%jac%init(model%nvars, model%nnode)
  end subroutine

  subroutine compute(this, u, udot, t, dt, stat)
    class(mfe_precon), intent(inout) :: this
    type(mfe1_vector), intent(in) :: u, udot
    real(r8), intent(in) :: t, dt
    integer, intent(out) :: stat
    call eval_jacobian(this%env, this%model, this%jac, u, udot, t, dt, stat)
  end subroutine

  subroutine apply(this, u)
    class(mfe_precon), intent(in) :: this
    type(mfe1_vector), intent(inout) :: u
    call this%jac%solve(u%array)
  end subroutine

  subroutine eval_jacobian(env, model, jac, u, udot, t, h, stat)

    use mfe_model_type
    use block_linear_solver
    use timer_tree_type

    type(mfe_env), intent(in) :: env
    type(mfe_model), intent(inout) :: model
    type(btd_matrix), intent(inout) :: jac
    type(mfe1_vector), intent(in) :: u, udot
    real(r8), intent(in) :: t, h
    integer, intent(out) :: stat

    integer :: i, j, n
    real(r8) :: diag(u%neqns+1,u%neqns+1,u%nnode)
    character(:), allocatable :: errmsg

    call start_timer('preprocessing')
    call model%disc%update(u, udot)
    call stop_timer('preprocessing')

    call start_timer('mass-matrix')
    call model%disc%eval_mass_matrix(factor = -1.0_r8 / h)
    call stop_timer('mass-matrix')

    ! Capture the unscaled diagonal for preconditioning.
    call start_timer('diag-pc')
    call model%disc%assemble_diagonal(diag)
    diag = (-h) * diag
    do i = 1, size(model%bc_dir)
      do j = 1, size(model%bc_dir(i)%index)
        n = model%bc_dir(i)%index(j)
        diag(:,i,n) = 0.0_r8
        diag(i,:,n) = 0.0_r8
        diag(i,i,n) = 1.0_r8
      end do
    end do
    call stop_timer('diag-pc')

    call start_timer('eval_dfdy')
    call model%disc%eval_dfdy(t, stat, errmsg)
    call stop_timer('eval_dfdy')

    if (stat /= 0) then
      call env%log%info('eval_dfdy: ' // errmsg)
      stat = 1
      return
    end if

    call start_timer('assembly')
    call model%disc%assemble_matrix(jac)
    do i = 1, size(model%bc_dir)
      do j = 1, size(model%bc_dir(i)%index)
        n = model%bc_dir(i)%index(j)
        call jac%set_dir_var(i, n)
      end do
    end do
    call stop_timer('assembly')

    ! Diagonal preconditioning.
    call start_timer('diag-pc')
    call vfct(diag)
    call vmslv(diag, jac%l)
    call vmslv(diag, jac%d)
    call vmslv(diag, jac%u)
    call stop_timer('diag-pc')

    ! Factorize the Jacobian.
    call start_timer('factorization')
    call jac%factor
    call stop_timer('factorization')

  end subroutine eval_jacobian

end module mfe_precon_type
