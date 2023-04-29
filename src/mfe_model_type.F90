#include "f90_assert.fpp"

module mfe_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use idaesol_type, only: idaesol_model
  use mfe_constants, only: NEQNS
  use vector_class
  use mfe1_vector_type
  use timer_tree_type
  implicit none
  private
  
  type, extends(idaesol_model), public :: mfe_model
    integer :: nnode
  contains
    procedure :: init
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
    procedure :: check_state
  end type mfe_model

contains

  subroutine init (this, nnode)
    class(mfe_model), intent(out) :: this
    integer, intent(in) :: nnode
    this%nnode = nnode
  end subroutine

  subroutine alloc_vector(this, vec)
    use mfe1_vector_type
    class(mfe_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec
    type(mfe1_vector), allocatable :: tmp
    allocate(tmp)
    call tmp%init(NEQNS, this%nnode)
    call move_alloc(tmp, vec)
  end subroutine
  
  subroutine compute_f(this, t, u, udot, f)
    use mfe_procs, only: eval_residual
    class(mfe_model) :: this
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
          call eval_residual(u, udot, t, f)
        end select
      end select
    end select
    call stop_timer('compute_f')
    call start_timer('integration')
  end subroutine
  
  subroutine apply_precon(this, t, u, f)
    use mfe_procs, only: apply_inverse_jacobian
    class(mfe_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u
    class(vector), intent(inout) :: f
    call stop_timer('integration')
    call start_timer('apply_precon')
    select type (f)
    class is (mfe1_vector)
      call apply_inverse_jacobian(f)
    end select
    call stop_timer('apply_precon')
    call start_timer('integration')
  end subroutine
  
  subroutine compute_precon(this, t, u, udot, dt)
    use mfe_procs, only: eval_jacobian
    class(mfe_model) :: this
    real(r8), intent(in)  :: t, dt
    class(vector), intent(inout) :: u, udot
    integer :: stat
    call stop_timer('integration')
    call start_timer('compute_precon')
    select type (u)
    class is (mfe1_vector)
      select type (udot)
      class is (mfe1_vector)
        call eval_jacobian(u, udot, t, dt, stat)
        INSIST(stat == 0)
      end select
    end select
    call stop_timer('compute_precon')
    call start_timer('integration')
  end subroutine
  
  subroutine du_norm(this, u, du, error)
    use norm_procs, only: eval_norm
    class(mfe_model) :: this
    class(vector), intent(in) :: u, du
    real(r8), intent(out) :: error
    call stop_timer('integration')
    call start_timer('du_norm')
    select type (u)
    class is (mfe1_vector)
      select type (du)
      class is (mfe1_vector)
        error = eval_norm(du, 0)
      end select
    end select
    call stop_timer('du_norm')
    call start_timer('integration')
  end subroutine
  
  subroutine check_state (this, u, stage, stat)
    use norm_procs, only: check_soln
    class(mfe_model) :: this
    class(vector), intent(in) :: u
    integer, intent(in)  :: stage
    integer, intent(out) :: stat
    call stop_timer('integration')
    call start_timer('check_state')
    select type (u)
    class is (mfe1_vector)
      call check_soln(u, stage, stat)
    end select
    call stop_timer('check_state')
    call start_timer('integration')
  end subroutine
  
end module mfe_model_type
