#include "f90_assert.fpp"

module mfe_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use idaesol_type, only: idaesol_model
  use mfe_env_type
  use mfe_constants, only: NEQNS
  use mfe_model_type
  use mfe_precon_type
  use vector_class
  use mfe1_vector_type
  use parameter_list_type
  use timer_tree_type
  implicit none
  private

  type, extends(idaesol_model), public :: mfe_idaesol_model
    integer :: nnode
    type(mfe_env), pointer :: env => null() ! reference only
    type(mfe_model), pointer :: model => null() ! reference only
    type(mfe_precon) :: precon
    real(r8) :: dxmin ! for check_state
    real(r8) :: rtol
    real(r8), allocatable :: ptol(:)
    real(r8), allocatable :: dx(:)  ! saved state; FIXME!
  contains
    procedure :: init
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
    procedure :: check_state
  end type mfe_idaesol_model

contains

  subroutine init(this, env, model, params, stat, errmsg)

    class(mfe_idaesol_model), intent(out) :: this
    type(mfe_env), target, intent(in) :: env
    type(mfe_model), target, intent(in) :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    this%env => env

    this%nnode = model%nnode
    this%model => model
    allocate(this%dx(this%nnode-1))

    call this%precon%init(this%env, this%model)

    call params%get('dxmin', this%dxmin, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%dxmin < 0.0_r8) then !TODO: <= 0.0?
      stat = 1
      errmsg = 'dxmin is <= 0.0'
      return
    end if

    call params%get('rtol', this%rtol, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%rtol <= 0.0_r8) then
      stat = 1
      errmsg = 'rtol is <= 0.0'
      return
    end if

    call params%get('ptol', this%ptol, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (size(this%ptol) /= this%model%neqns+1) then
      stat = 1
      errmsg = 'wrong number of ptol values'
      return
    else if (any(this%ptol <= 0.0)) then
      stat = 1
      errmsg = 'ptol is <= 0'
      return
    end if

  end subroutine

  subroutine alloc_vector(this, vec)
    use mfe1_vector_type
    class(mfe_idaesol_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec
    type(mfe1_vector), allocatable :: tmp
    allocate(tmp)
    call tmp%init(this%model%neqns, this%nnode)
    call move_alloc(tmp, vec)
  end subroutine

  subroutine compute_f(this, t, u, udot, f)
    class(mfe_idaesol_model) :: this
    real(r8), intent(in)  :: t
    class(vector), intent(inout) :: u, udot
    class(vector), intent(inout) :: f
    call stop_timer('integration')
    call start_timer('compute_f')
    select type (u)
    class is (mfe1_vector)
      select type (udot)
      class is (mfe1_vector)
        select type (f)
        class is (mfe1_vector)
          call this%model%eval_residual(u, udot, t, f)
        end select
      end select
    end select
    call stop_timer('compute_f')
    call start_timer('integration')
  end subroutine

  subroutine apply_precon(this, t, u, f)
    class(mfe_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u
    class(vector), intent(inout) :: f
    call stop_timer('integration')
    call start_timer('apply_precon')
    select type (f)
    class is (mfe1_vector)
      call this%precon%apply(f)
    end select
    call stop_timer('apply_precon')
    call start_timer('integration')
  end subroutine

  subroutine compute_precon(this, t, u, udot, dt)
    class(mfe_idaesol_model) :: this
    real(r8), intent(in)  :: t, dt
    class(vector), intent(inout) :: u, udot
    integer :: stat
    call stop_timer('integration')
    call start_timer('compute_precon')
    select type (u)
    class is (mfe1_vector)
      select type (udot)
      class is (mfe1_vector)
        call this%precon%compute(u, udot, t, dt, stat)
        INSIST(stat == 0)
      end select
    end select
    call stop_timer('compute_precon')
    call start_timer('integration')
  end subroutine


!FIXME: du_norm relies on check_state having cached the relevant dx

  subroutine du_norm(this, u, du, error)

    class(mfe_idaesol_model) :: this
    class(vector), intent(in) :: u, du
    real(r8), intent(out) :: error

    call stop_timer('integration')
    call start_timer('du_norm')

    select type (u)
    class is (mfe1_vector)
      select type (du)
      class is (mfe1_vector)

        ! Weighted max-norm.
        error = maxval(maxval(abs(du%array),dim=2) / this%ptol)

        ! Relative norm on element lengths.
        associate (dx => du%array(du%neqns+1,:))
          error = max(error, maxval(abs(dx(2:size(dx))-dx(1:size(dx)-1))/this%dx) / this%rtol)
        end associate

      end select
    end select

    call stop_timer('du_norm')
    call start_timer('integration')

  end subroutine

  subroutine check_state(this, u, stage, stat)

    use string_utilities, only: i_to_c

    class(mfe_idaesol_model) :: this
    class(vector), intent(in) :: u
    integer, intent(in)  :: stage
    integer, intent(out) :: stat

    integer :: loc

    call stop_timer('integration')
    call start_timer('check_state')

    select type (u)
    class is (mfe1_vector)

      associate (x => u%array(NEQNS+1,:))
        this%dx = x(2:size(x)) - x(1:size(x)-1)
      end associate
      loc = minloc(this%dx, dim=1)

      !TODO: return errmsg and let caller handle output?
      if (this%dx(loc) < this%dxmin) then
        call this%env%log%info('check_state: bad element ' // i_to_c(loc))
        stat = 1
      else
        stat = 0
      end if

    end select

    call stop_timer('check_state')
    call start_timer('integration')

  end subroutine

end module mfe_idaesol_model_type
