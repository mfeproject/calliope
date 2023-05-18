module mfe1_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_env_type
  use idaesol_type
  use mfe_idaesol_model_type
  use mfe_model_type
  use mfe1_vector_type
  use parameter_list_type
  implicit none
  private

  type, public :: mfe1_solver
    type(mfe_env), pointer :: env => null() ! reference only
    type(mfe_idaesol_model) :: integ_model
    type(idaesol) :: integ
    type(mfe_model), pointer :: model => null()
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: advance
    procedure :: last_time
    procedure :: get_last_state_copy
    procedure :: get_interpolated_state
    procedure :: write_metrics
  end type

contains

  subroutine init(this, env, model, params, stat, errmsg)

    class(mfe1_solver), intent(out), target :: this
    type(mfe_env), intent(in), target :: env
    type(mfe_model), intent(inout), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    logical :: flag
    integer :: lun

    this%env => env
    this%model => model

    call this%integ_model%init(this%env, this%model, params, stat, errmsg)
    if (stat /= 0) return
    call this%integ%init(this%integ_model, params, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'integrator input error: ' // errmsg
      return
    end if

    call params%get('verbose-stepping', flag, default=.false.)
    if (flag) then
      open(newunit=lun,file='bdfout',action='write',status='replace')
      call this%integ%set_verbose_stepping(lun)
    end if

  end subroutine

  subroutine set_initial_state(this, t, u, udot)
    class(mfe1_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(in) :: u, udot
    call this%integ%set_initial_state(t, u, udot)
  end subroutine

  subroutine advance(this, t, u, hnext, stat, errmsg)
    class(mfe1_solver), intent(inout) :: this
    real(r8), intent(inout) :: t
    type(mfe1_vector), intent(inout) :: u
    real(r8), intent(out)   :: hnext
    integer,  intent(out)   :: stat
    character(:), allocatable :: errmsg
    call this%integ%advance(t, u, hnext, stat, errmsg)
  end subroutine

  function last_time(this) result(t)
    class(mfe1_solver), intent(in) :: this
    real(r8) :: t
    t = this%integ%last_time()
  end function

  subroutine get_last_state_copy(this, copy)
    class(mfe1_solver), intent(in) :: this
    type(mfe1_vector), intent(inout) :: copy
    call this%integ%get_last_state_copy(copy)
  end subroutine

  subroutine get_interpolated_state(this, t, u)
    class(mfe1_solver), intent(in) :: this
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(inout) :: u
    call this%integ%get_interpolated_state(t, u)
  end subroutine

  subroutine write_metrics(this, unit)
    class(mfe1_solver), intent(in) :: this
    integer, intent(in) :: unit
    call this%integ%write_metrics(unit)
  end subroutine

end module
