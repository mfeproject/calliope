!!
!! IDAESOL_TYPE
!!
!! A solver for index-1 DAE in implicit form using the BDF2 method.
!! This F2008 version is adapted from much earlier F77 and F95 implementations.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2013 Neil N. Carlson
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module idaesol_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use state_history_type
  use nka_type
  use vector_class
  use string_utilities, only: i_to_c
  implicit none
  private

  type, public :: idaesol
    private
    class(idaesol_model), pointer :: model => null()  ! reference only -- do not own
    integer  :: n                   ! number of unknowns
    integer  :: seq = -1            ! number of steps taken
    real(r8) :: hlast               ! last step size
    real(r8) :: hpc                 ! step size built into the current preconditioner
    logical  :: usable_pc = .false. ! whether the current preconditioner is usable
    integer  :: freeze_count = 0    ! don't increase step size for this number of steps
    integer  :: mitr = 5            ! maximum number of nonlinear iterations
    real(r8) :: ntol = 0.1_r8       ! nonlinear solver error tolerance (relative to 1)
    type(nka), allocatable :: nka   ! nonlinear Krylov accelerator
    type(state_history) :: uhist    ! solution history structure
    real(r8) :: h_min, h_max  !TODO: change to hmin, hmax after counters refactored
    integer  :: max_try

    !! Persistent temporary workspace
    class(vector), allocatable :: u, u0, up
    class(vector), allocatable :: du, udot  ! local to bce_step

    !! Perfomance counters
    integer :: pcfun_calls = 0      ! number of calls to PCFUN
    integer :: updpc_calls = 0      ! number of calls to UPDPC
    integer :: updpc_failed = 0     ! number of UPDPC calls returning an error
    integer :: retried_bce = 0      ! number of retried BCE steps
    integer :: failed_bce = 0       ! number of completely failed BCE steps
    integer :: rejected_steps = 0   ! number of steps rejected on error tolerance
    real(r8) :: hmin = huge(1.0_r8) ! minimum step size used on a successful step
    real(r8) :: hmax = tiny(1.0_r8) ! maximum step size used on a successful step

    !! Diagnostics
    integer :: unit = 0
    logical :: verbose = .false.
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: advance
    procedure :: step
    procedure :: commit_state
    procedure :: get_interpolated_state
    procedure :: get_last_state_copy
    procedure :: get_last_state_view
    procedure :: last_time
    procedure :: last_step_size
    procedure :: set_verbose_stepping
    procedure :: set_quiet_stepping
    procedure :: write_metrics
    procedure :: get_metrics
  end type idaesol

  type, abstract, public :: idaesol_model
  contains
    procedure(alloc_vector), deferred :: alloc_vector
    procedure(compute_f), deferred :: compute_f
    procedure(apply_precon), deferred :: apply_precon
    procedure(compute_precon), deferred :: compute_precon
    procedure(du_norm), deferred :: du_norm
    procedure(check_state), deferred :: check_state
  end type

  abstract interface
    subroutine alloc_vector(this, vec)
      import idaesol_model, vector
      class(idaesol_model), intent(in) :: this
      class(vector), allocatable, intent(out) :: vec
    end subroutine
    subroutine compute_f(this, t, u, udot, f)
      import idaesol_model, vector, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t
      class(vector), intent(inout) :: u, udot ! TODO: why inout and not in?
      class(vector), intent(inout) :: f
    end subroutine
    subroutine apply_precon(this, t, u, f)
      import idaesol_model, vector, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t
      class(vector), intent(inout) :: u ! TODO: why inout and not in?
      class(vector), intent(inout) :: f
    end subroutine
    subroutine compute_precon(this, t, u, udot, dt)
      import idaesol_model, vector, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t, dt
      class(vector), intent(inout) :: u, udot ! TODO: why inout and not in?
    end subroutine
    subroutine du_norm(this, u, du, error)
      import :: idaesol_model, vector, r8
      class(idaesol_model) :: this
      class(vector), intent(in) :: u, du
      real(r8), intent(out) :: error
    end subroutine
    subroutine check_state(this, u, stage, stat)
      import :: idaesol_model, vector
      class(idaesol_model) :: this
      class(vector), intent(in) :: u
      integer, intent(in) :: stage
      integer, intent(out) :: stat
    end subroutine
  end interface

  real(r8), parameter, private :: RMIN = 0.25_r8
  real(r8), parameter, private :: RMAX = 4.0_r8
  real(r8), parameter, private :: MARGIN = 3.0_r8

  !! Successful STATUS return codes:
  integer, parameter, public :: SOLVED_TO_TOUT = 1
  integer, parameter, public :: SOLVED_TO_NSTEP = 2

  !! Unsuccessful STATUS return codes:
  integer, parameter, public :: BAD_INPUT = -1
  integer, parameter, public :: STEP_FAILED = -2
  integer, parameter, public :: STEP_SIZE_TOO_SMALL = -3

contains

  !! Return a pointer to the current system state.
  subroutine get_last_state_view (this, view)
    class(idaesol), intent(in) :: this
    class(vector), pointer, intent(out) :: view
    call this%uhist%get_last_state_view(view)
  end subroutine get_last_state_view

  !! Return a copy of the current system state.
  subroutine get_last_state_copy(this, copy)
    class(idaesol), intent(in) :: this
    class(vector), intent(inout) :: copy
    call this%uhist%get_last_state_copy(copy)
  end subroutine get_last_state_copy

  !! Return the current system time.
  function last_time (this) result (t)
    class(idaesol), intent(in) :: this
    real(r8) :: t
    t = this%uhist%last_time()
  end function last_time

  !! Return the last time step size used (successfully).
  function last_step_size (this) result (h)
    class(idaesol), intent(in) :: this
    real(r8) :: h
    h = this%hlast
  end function last_step_size

  !! Computes the interpolated system state vector at time T.  The specified
  !! time should generally be contained in the interval between the current
  !! time and the previous time.  Typically the array U would have the same
  !! size as the system state vector, but more generally U may return any
  !! contiguous segment of the interpolated state starting at index FIRST
  !! (default 1) and length the size of U.

  subroutine get_interpolated_state(this, t, u)
    class(idaesol), intent(in) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u
    call this%uhist%interp_state(t, u)
  end subroutine get_interpolated_state

  subroutine write_metrics (this, unit)
    class(idaesol), intent(in) :: this
    integer, intent(in) :: unit
    write(unit,fmt='(/,a,i0,a,g0.5,a,g0.5)') &
      'STEP=', this%seq, ', T=', this%uhist%last_time(), ', H=', this%hlast
    write(unit,fmt='(a,i0,":",i0,a,5(i0,:,":"))') &
      'NFUN:NPC=', this%pcfun_calls, this%updpc_calls, &
      ', NPCF:NNR:NNF:NSR=', this%updpc_failed, &
      this%retried_bce, this%failed_bce, this%rejected_steps
  end subroutine write_metrics

  subroutine get_metrics (this, nstep, hmin, hmax, counters)
    class(idaesol), intent(in) :: this
    integer, intent(out), optional :: nstep, counters(:)
    real(r8), intent(out), optional :: hmin, hmax
    if (present(nstep)) nstep = this%seq
    if (present(hmin)) hmin = this%hmin
    if (present(hmax)) hmax = this%hmax
    if (present(counters)) then
      ASSERT(size(counters) == 6)
      counters(1) = this%pcfun_calls
      counters(2) = this%updpc_calls
      counters(3) = this%updpc_failed
      counters(4) = this%retried_bce
      counters(5) = this%failed_bce
      counters(6) = this%rejected_steps
    end if
  end subroutine get_metrics

  subroutine set_verbose_stepping (this, unit)
    class(idaesol), intent(inout) :: this
    integer, intent(in) :: unit
    this%unit = unit
    this%verbose = .true.
  end subroutine set_verbose_stepping

  subroutine set_quiet_stepping (this)
    class(idaesol), intent(inout) :: this
    this%verbose = .false.
  end subroutine set_quiet_stepping

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(idaesol), intent(out) :: this
    class(idaesol_model), intent(in), target :: model
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: maxv
    real(r8) :: vtol

    this%model => model

    call params%get('nlk-max-iter', this%mitr, default=5, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%mitr < 1) then
      stat = -1
      errmsg = '"nlk-max-iter" must be >= 1'
      return
    end if

    call params%get('nlk-tol', this%ntol, default=0.1_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%ntol <= 0 .or. this%ntol > 1) then
      stat = -1
      errmsg = '"nlk-tol" must be in (0,1]'
      return
    end if

    call params%get('nlk-max-vec', maxv, default=this%mitr-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (maxv <= 0) then
      stat = -1
      errmsg = '"nlk-max-vec" must be > 0'
      return
    end if
    maxv = min(maxv, this%mitr-1)

    call params%get('nlk-vec-tol', vtol, default=0.01_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (vtol <= 0) then
      stat = -1
      errmsg = '"nlk-vec-tol" must be > 0'
      return
    end if

    call params%get('h-max', this%h_max, default=huge(this%h_max), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%h_max <= 0) then
      stat = -1
      errmsg = '"h-max" must be > 0'
      return
    end if

    call params%get('h-min', this%h_min, default=tiny(this%h_min), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%h_min < 0) then
      stat = -1
      errmsg = '"h-min" must be >= 0'
      return
    else if (this%h_min > this%h_max) then
      stat = -1
      errmsg = '"h-min" must be <= "h-max"'
      return
    end if

    call params%get('max-try', this%max_try, default=10, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%max_try < 1) then
      stat = -1
      errmsg = '"max-try" must be > 0'
      return
    end if

    call model%alloc_vector(this%u)
    call model%alloc_vector(this%u0)
    call model%alloc_vector(this%up)
    call model%alloc_vector(this%du)
    call model%alloc_vector(this%udot)

    !! Initialize the NKA structure.
    if (maxv > 0) then
      allocate(this%nka)
      call this%nka%init(this%u, maxv)
      call this%nka%set_vec_tol(vtol)
    end if

    !! We need to maintain 3 solution vectors for quadratic extrapolation.
    call this%uhist%init(3, this%u)

  end subroutine init

  !! Sets the initial state (t, u, du/dt) of the DAE system.  This must be
  !! called before beginning integration.  It is important that the passed
  !! state satisfy the DAE, f(t,u,du/dt) = 0.  This is trivial for explicit
  !! DE, du/dt = f(t,u), but is often not for DAEs.

  subroutine set_initial_state(this, t, u, udot)
    class(idaesol), intent(inout) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u, udot
    call this%uhist%flush(t, u, udot)
    this%seq = 0
  end subroutine set_initial_state

  !! This advances the state by a single "robust" time step. The target time
  !! for the step is T, but in the case of failure, the step size will be
  !! repeatedly reduced until a successful step is achieved, or until a bound
  !! is exceeded: either the minimum allowed step size or the maximum of
  !! allowed attempts. If successful, STAT returns 0 and the time and solution
  !! are committed as the current state, and also returned in T and U.  Note
  !! that the returned value of T may differ from its input value.  HNEXT
  !! returns the requested next step size. If the step is unsuccesful, STAT
  !! returns a negative value and ERRMSG is assigned an explanatory message.

  subroutine advance(this, t, u, hnext, stat, errmsg)

    class(idaesol), intent(inout) :: this
    real(r8), intent(inout) :: t
    class(vector), intent(inout) :: u
    real(r8), intent(out)   :: hnext
    integer,  intent(out)   :: stat
    character(:), allocatable :: errmsg

    integer :: try
    real(r8) :: h

    do try = 1, this%max_try
      !! Check for a too-small step size.
      h = t - this%uhist%last_time()
      if (h < this%h_min) then
        hnext = h
        stat = STEP_SIZE_TOO_SMALL
        errmsg = 'next time step size is too small'
        return
      end if

      !! Attempt a BDF2 step.
      call step (this, t, u, hnext, stat)
      if (stat == 0) then ! successful
        hnext = min(hnext, this%h_max)
        call commit_state (this, t, u)  !TODO: Should the client have the responsibility of doing this?
        return
      else  ! failed; try again with the suggested step size.
        if (this%verbose) write(this%unit,fmt=1) hnext/h
        t = this%uhist%last_time() + hnext
      end if
    end do

    stat = STEP_FAILED
    errmsg = 'time step unsuccessful after repeated attempts ("max-try"=' // i_to_c(this%max_try) // ')'

    1 format(2x,'Changing H by a factor of ',f6.3)

  end subroutine advance

  !! This subroutine commits the passed state as the new current state.  The
  !! state should be that returned by the STEP procedure.  The reason for not
  !! having STEP automatically commit the state is for use cases where other
  !! PDEs are being simultaneously integrated (time split).  This allows the
  !! driving procedure control over acceptance of the overarching time step.

  subroutine commit_state (this, t, u)

    class(idaesol), intent(inout) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u

    real(r8) :: h

    h = t - this%uhist%last_time()
    ASSERT(h > 0.0_r8)
    call this%uhist%record_state (t, u)

    this%hlast = h
    this%seq = this%seq + 1
    this%freeze_count = max(0, this%freeze_count - 1)

    this%hmin = min(h, this%hmin)
    this%hmax = max(h, this%hmax)

  end subroutine commit_state

  !! Starting from the current state, this subroutine takes a single step to
  !! compute the solution at time T.  If successful (STAT==0), the advanced
  !! solution is returned in U, and HNEXT returns the next step size to use.
  !! Note that this does not advance the current state; COMMIT_STATE must be
  !! called with the computed solution to advance the state.  Normally a BDF2
  !! step is taken, except at the start of integration when insufficient state
  !! history is available and a BDF1 step is taken instead.  This subroutine
  !! manages the update of the preconditioner, and will retry a failed step
  !! with a fresh preconditioner as needed.  If the specified step to T is
  !! ultimately unsuccessful, STAT returns a nonzero value and HNEXT returns
  !! a (reduced) step size to attempt next.  Failure can occur for several
  !! reasons: failure to solve the nonlinear system (STAT==1), predictor
  !! error exceeded tolerance (STAT==2), and inadmissable predicted solution
  !! (STAT==3).

  subroutine step (this, t, u, hnext, stat)

    class(idaesol), intent(inout) :: this
    real(r8), intent(in)  :: t
    class(vector), intent(inout) :: u ! data is intent(out)
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: stat

    real(r8) :: eta, etah, h, t0, tlast, perr, dt(3)
    logical  :: fresh_pc, predictor_error

    ASSERT(this%seq >= 0)

    tlast = this%uhist%last_time()
    h = t - tlast
    INSIST(h > 0)

    if (this%uhist%depth() == 2) then ! trapezoid step to bootstrap
      call trap_step (this, t, u, stat)
      if (stat /= 0) then
        hnext = 0.1_r8 * h
      else
        hnext = h
      end if
      return
    end if

    if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h

    !! Predicted solution and base point for BCE step.
    if (this%uhist%depth() == 2) then ! BDF1
      etah = h
      t0 = tlast
      call this%uhist%interp_state(t,  this%up, order=1)
      call this%uhist%interp_state(t0, this%u0, order=0)
    else  ! BDF2
      eta = (this%hlast + h) / (this%hlast + 2.0_r8 * h)
      etah = eta * h
      t0 = tlast + (1.0_r8 - eta)*h
      call this%uhist%interp_state(t,  this%up, order=2)
      call this%uhist%interp_state(t0, this%u0, order=1)
    end if

    !! Check the predicted solution for admissibility.
    call this%model%check_state (this%up, 0, stat)
    if (stat /= 0) then ! it's bad; cut h and retry.
      this%rejected_steps = this%rejected_steps + 1
      if (this%verbose) write(this%unit,fmt=7)
      hnext = 0.25_r8 * h
      this%freeze_count = 1
      this%usable_pc = .false.
      stat = 3
      return
    end if

    fresh_pc = .false.

    !! If the PC step size is too different than the current step size we tag
    !! it as unusable in order to preempt a possible nonlinear solve failure.
    if (this%usable_pc) then
      if (this%hpc/etah > 1.0_r8 + MARGIN) this%usable_pc = .false.
      if (etah/this%hpc > 1.0_r8 + MARGIN) this%usable_pc = .false.
    end if

    BCE: do

      !! Update the preconditioner if necessary.
      if (.not.this%usable_pc) then
        this%updpc_calls = this%updpc_calls + 1
        !call this%model%compute_precon(t, up, (up-u0)/etah, etah)
        call this%udot%copy(this%up)
        call this%udot%update(-1.0_r8, this%u0)
        call this%udot%scale(1.0_r8/etah)
        call this%model%compute_precon(t, this%up, this%udot, etah)
        if (this%verbose) write(this%unit,fmt=3) t
        this%hpc = etah
        this%usable_pc = .true.
        fresh_pc = .true.
      end if

      !! Solve the nonlinear BCE system.
      call u%copy(this%up) ! Initial solution guess is the predictor.
      call bce_step(this, t, etah, this%u0, u, stat)
      if (stat == 0) exit BCE ! the BCE step was successful.

      if (fresh_pc) then ! preconditioner was fresh; cut h and return error condition.
        if (stat == 2) then ! inadmissible iterate generated
          this%rejected_steps = this%rejected_steps + 1
          hnext = 0.25_r8 * h
          this%freeze_count = 1
          this%usable_pc = .false.
        else
          this%failed_bce = this%failed_bce + 1
          hnext = 0.5_r8 * h
          this%freeze_count = 2
        end if
        stat = 1
        return
      else ! update the preconditioner and retry the nonlinear solve.
        this%retried_bce = this%retried_bce + 1
        this%usable_pc = .false.
        cycle BCE
      end if

    end do BCE

    predictor_error = (this%seq >= 2)

    if (predictor_error) then

      !! Predictor error control.
      !du = u - up
      call this%du%copy(u)
      call this%du%update(-1.0_r8, this%up)
      call this%model%du_norm(u, this%du, perr)
      if (perr < 4.0_r8) then ! accept the step.
        if (this%verbose) write(this%unit,fmt=4) perr
        stat = 0
      else ! reject the step; cut h and return error condition.
        this%rejected_steps = this%rejected_steps + 1
        if (this%verbose) write(this%unit,fmt=5) perr
        hnext = 0.5_r8 * h
        this%freeze_count = 1
        stat = 2
        return
      end if

      !! Select the next step size based on the predictor error and past step
      !! size history, but don't change the step size by too great a factor.
      dt(1) = h
      dt(2:) = h + this%uhist%time_deltas()
      call select_step_size(dt, perr, hnext)
      hnext = max(RMIN*h, min(RMAX*h, hnext))
      if (this%freeze_count /= 0) hnext = min(h, hnext)

    else

      if (this%verbose) write(this%unit,fmt=6)
      hnext = h
      stat = 0

    end if

    1 format(/,'BDF2 step ',i6,': T=',es12.5,', H=',es12.5)
    3 format(2x,'Preconditioner updated at T=',es12.5)
    4 format(2x,'Step accepted: perr=',es12.5)
    5 format(2x,'Step REJECTED: perr=',es12.5)
    6 format(2x,'Step accepted: no local error control')
    7 format(2x,'Step REJECTED: inadmissible predicted solution')

  end subroutine step

  !! This subroutine selects the next step size based on the predictor error
  !! for the current step and the recent history of step sizes. The step size
  !! is chosen so that the estimated predictor error on the next step is 1/2.
  !! Recall that the correction norm computed by the model is a scaled norm
  !! with value 1 separating 'small' from 'large'. See [1] for a description
  !! of the scheme, which leads to a problem of finding the root of a third
  !! order polynomial that is solved here using Newton iteration.

  subroutine select_step_size (dt, perr, h)

    real(r8), intent(in)  :: dt(:), perr
    real(r8), intent(out) :: h

    real(r8), parameter :: tol = 0.001_r8
    real(r8) :: a, dh, phi, dphi

    ASSERT(size(dt) == 3)

    a = 0.5_r8*dt(1)*dt(2)*dt(3)/max(perr,0.001_r8)
    h = dt(1)
    do ! until converged
      phi  = h*(h + dt(1))*(h + dt(2)) - a
      dphi = (2.0_r8*h + dt(1))*(h + dt(2)) + h*(h + dt(1))
      dh = phi / dphi
      h = h - dh
      if (abs(dh) / h < tol) exit
    end do

  end subroutine select_step_size

  subroutine trap_step(this, t, u, stat)

    class(idaesol), intent(inout) :: this
    real(r8), intent(in)  :: t
    class(vector), intent(inout) :: u
    integer,  intent(out) :: stat

    real(r8) :: tlast, t0, h, etah

    tlast = this%uhist%last_time()
    h = t - tlast
    etah = 0.5_r8 * h
    t0 = tlast + etah

    if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h

    call this%uhist%interp_state(t,  u,  order=1)
    call this%uhist%interp_state(t0, this%u0, order=1)

    !! Check the predicted solution for admissibility.
    call this%model%check_state(u, 0, stat)
    if (stat /= 0) then ! it's bad; cut h and retry.
      this%rejected_steps = this%rejected_steps + 1
      if (this%verbose) write(this%unit,fmt=7)
      stat = 3
      return
    end if

    !! Update the preconditioner.
    this%updpc_calls = this%updpc_calls + 1
    call this%udot%copy(u)
    call this%udot%update(-1.0_r8, this%u0)
    call this%udot%scale(1.0_r8/etah)
    call this%model%compute_precon(t, u, this%udot, etah)
    if (this%verbose) write(this%unit,fmt=3) t
    this%hpc = etah

    !! Solve the nonlinear BCE system.
    call bce_step(this, t, etah, this%u0, u, stat)
    if (stat /= 0) then
      this%failed_bce = this%failed_bce + 1
      this%freeze_count = 1
      stat = 1
    else
      if (this%verbose) write(this%unit,fmt=6)
      stat = 0
    end if

    1 format(/,'TRAP step ',i6,': T=',es12.5,', H=',es12.5)
    3 format(2x,'Preconditioner updated at T=',es12.5)
    6 format(2x,'Step accepted: no local error control')
    7 format(2x,'Step REJECTED: inadmissible predicted solution')

  end subroutine trap_step

  !! The backward Cauchy-Euler (BCE) method applied to the implicit DAE
  !! f(t,u,u') = 0 yields the nonlinear system f(t,u,(u-u0)/h) = 0 for
  !! advancing the solution from a given state u0 at time t - h to the
  !! unknown state u at time t. This subroutine solves that nonlinear system
  !! using an accelerated fixed point iteration [1] for the preconditioned
  !! system g(u) = pc(f(t,u,(u-u0)/h)) = 0:
  !!
  !!    u given
  !!    do until converged:
  !!      du <-- g(u)
  !!      du <-- NKA(du)
  !!      u  <-- u - du
  !!    end do
  !!
  !! The nonlinear Krylov acceleration (NKA) procedure uses information about
  !! g' gleaned from the unaccelerated correction du=g(u) and previous g
  !! values to compute an improved correction.  The preconditioning function
  !! pc() is typically an approximate solution of the Newton correction
  !! equation J*du = f(t,u,(u-u0)/h) where J is an approximation to the
  !! Jacobian of f(t,u,(u-u0)/h) as a function of u.
  !!
  !! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
  !!     weighted moving finite element code I: in one dimension", SIAM J.
  !!     Sci. Comput;, 19 (1998), pp. 728-765..

  subroutine bce_step(this, t, h, u0, u, stat)

    type(idaesol), intent(inout) :: this
    real(r8), intent(in)    :: t, h
    class(vector), intent(in) :: u0
    class(vector), intent(inout) :: u
    integer,  intent(out)   :: stat

    integer  :: itr
    real(r8) :: error

    call this%nka%restart

    itr = 0
    do

      if (itr >= this%mitr) then  ! too many nonlinear iterations
        if (this%verbose) write(this%unit,fmt=1) itr, error
        stat = 1
        exit
      end if

      itr = itr + 1

      !! Evaluate the nonlinear function and precondition it.
      this%pcfun_calls = this%pcfun_calls + 1

      !call this%model%compute_f(t, u, (u-u0)/h, du)
      call this%udot%copy(u)
      call this%udot%update(-1.0_r8, u0)
      call this%udot%scale(1.0_r8/h)
      call this%model%compute_f(t, u, this%udot, this%du)
      call this%model%apply_precon(t, u, this%du)

      !! NKA accelerated correction.
      if (allocated(this%nka)) call this%nka%accel_update(this%du)

      !! Next solution iterate.
      !u = u - du
      call u%update(-1.0_r8, this%du)

      !! Check the solution iterate for admissibility.
      call this%model%check_state (u, 1, stat)
      if (stat /= 0) then ! iterate is bad; bail.
        if (this%verbose) write(this%unit,fmt=4) itr
        stat = 2
        exit
      end if

      !! Error estimate.
      call this%model%du_norm(u, this%du, error)
      if (this%verbose) write(this%unit,fmt=3) itr, error

      !! Check for convergence.
      if (((error < this%ntol) .and. (itr > 1)) .or. (error < 0.01_r8 * this%ntol)) then
        if (this%verbose) write(this%unit,fmt=2) itr, error
        stat = 0
        exit
      end if

    end do

    1 format(2x,'NLK BCE solve FAILED: ',i3,' iterations (max), error=',es12.5)
    2 format(2x,'NLK BCE solve succeeded: ',i3,' iterations, error=',es12.5)
    3 format(2x,i3,': error=',es12.5)
    4 format(2x,'NLK BCE solve FAILED: inadmissible solution iterate: itr=',i3)

  end subroutine bce_step

end module idaesol_type
