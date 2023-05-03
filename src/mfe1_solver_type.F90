module mfe1_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use idaesol_type
  use mfe_idaesol_model_type
  use mfe_model_type
  use mfe1_vector_type
  use parameter_list_type
  implicit none
  private

  type, public :: mfe1_solver
    type(mfe_idaesol_model) :: integ_model
    type(idaesol) :: integ
    type(mfe_model), pointer :: model => null()
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: integrate
    procedure :: advance
    procedure :: last_time
    procedure :: get_last_state_copy
    procedure :: get_interpolated_state
    procedure :: get_metrics
    procedure :: write_metrics
  end type

  public :: SOLVED_TO_TOUT
  public :: SOLVED_TO_NSTEP
  public :: STEP_FAILED
  public :: STEP_SIZE_TOO_SMALL
  public :: BAD_INPUT

contains

  subroutine init(this, model, params, stat, errmsg)

    class(mfe1_solver), intent(out), target :: this
    type(mfe_model), intent(inout), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    logical :: flag
    integer :: lun

    this%model => model

    call this%integ_model%init(this%model, params, stat, errmsg)
    if (stat /= 0) return
    call this%integ%init(this%integ_model, params)

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

  subroutine integrate(this, hnext, status, nstep, tout, hmin, hmax, mtry)

    class(mfe1_solver), intent(inout) :: this
    real(r8), intent(inout) :: hnext
    integer,  intent(out)   :: status
    integer,  intent(in) :: nstep, mtry
    real(r8), intent(in) :: tout, hmin, hmax
    optional :: nstep, tout, hmin, hmax, mtry

    call this%integ%integrate(hnext, status, nstep, tout, hmin, hmax, mtry)

  end subroutine

  subroutine advance(this, hmin, mtry, t, u, hnext, stat)

    class(mfe1_solver), intent(inout) :: this
    real(r8), intent(in)    :: hmin
    integer,  intent(in)    :: mtry
    real(r8), intent(inout) :: t
    type(mfe1_vector), intent(inout) :: u
    real(r8), intent(out)   :: hnext
    integer,  intent(out)   :: stat

    call this%integ%advance(hmin, mtry, t, u, hnext, stat)

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

  subroutine get_metrics(this, nstep) ! TEMPORARILY
    class(mfe1_solver), intent(in) :: this
    integer, intent(out) :: nstep
    call this%integ%get_metrics(nstep=nstep)
  end subroutine

  subroutine write_metrics(this, unit)
    class(mfe1_solver), intent(in) :: this
    integer, intent(in) :: unit
    call this%integ%write_metrics(unit)
  end subroutine

end module
