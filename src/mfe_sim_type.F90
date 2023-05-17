#include "f90_assert.fpp"

module mfe_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_env_type
  use mfe_model_type
  use mfe1_solver_type
  use mfe1_vector_type
  use sim_event_queue_type
  use time_step_sync_type
  use timer_tree_type
  implicit none
  private

  type, public :: mfe_sim
    type(mfe_env), pointer :: env => null() ! reference only
    type(mfe_model) :: model
    type(mfe1_solver) :: solver
    type(mfe1_vector) :: u
    type(sim_event_queue) :: hard_eventq, soft_eventq
    type(time_step_sync) :: ts_sync
    real(r8) :: tlast, hlast
    real(r8) :: dt_init
    integer :: ofreq, mstep
  contains
    procedure :: init
    procedure :: run
  end type

  type, extends(event_action) :: output_event
  end type

contains

  subroutine init(this, env, params, stat, errmsg)

    use parameter_list_type
    use coord_grid_type

    class(mfe_sim), intent(out), target :: this
    type(mfe_env), intent(in), target :: env
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: nnode, itmp
    real(r8) :: t_init
    real(r8), allocatable :: u_array(:,:)
    type(coord_grid) :: grid
    type(mfe1_vector) :: udot

    this%env => env

    call grid%init(params, 'useg', 'niseg', 'ratio', stat, errmsg)
    if (stat /= 0) return

    call grid%get_grid(u_array)
    nnode = size(u_array,dim=2)

    call this%model%init(nnode, params, stat, errmsg)
    if (stat /= 0) return

    !! Ensure the initial solution is properly sized for the model
    if (size(u_array,dim=1) /= this%model%nvars) then
      stat = -1
      errmsg = 'wrong number of columns for useg'
      return
    end if

    !TODO: this should be an optional check (e.g., motion by curvature)
    if (any(u_array(1,2:) <= u_array(1,1:nnode-1))) then
      stat = -1
      errmsg = 'x coordinates not strictly increasing'
      return
    end if

    !TODO: should the model have a method for creating a vector? and then set its values?
    call this%u%init(x = u_array(1,:), u = u_array(2:,:))

    call this%model%set_boundary_values(this%u) ! get BV values from the initial solution

    block ! Solution output times
      real(r8), allocatable :: rarray(:)
      integer :: j
      call params%get('output-times', rarray, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (size(rarray) < 1) then
        stat = -1
        errmsg = '"output-times" must be assigned at least 1 values='
        return
      end if

      !! Add the output times to the soft event queue. The integrator will not
      !! hit these times exactly but will output an interpolated solution
      do j = 1, size(rarray)
        call this%soft_eventq%add_event(rarray(j), output_event())
        !call this%hard_eventq%add_event(rarray(j), output_event()) ! WORKS AS EXPECTED
      end do
    end block

    call params%get('mstep', this%mstep, default=huge(this%mstep), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%mstep < 0) then
      stat = -1
      errmsg = 'mstep must be >= 0'
      return
    end if
    call params%get('ofreq', this%ofreq, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%ofreq <= 0) this%ofreq = this%mstep

    call params%get('t-init', t_init, stat=stat, errmsg=errmsg)
    if (stat /= 0) return

    call params%get('h-init', this%dt_init, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%dt_init <= 0) then
      stat = -1
      errmsg = '"h-init" must be > 0'
      return
    end if

    !! Look-ahead for hard events and associated time_step_sync object
    call params%get('hard-event-lookahead', itmp, default=2, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (itmp < 1) then
      stat = -1
      errmsg = '"hard-event-lookahead" must be > 0'
      return
    end if
    this%ts_sync = time_step_sync(itmp)

    call this%solver%init(this%env, this%model, params, stat, errmsg)
    if (stat /= 0) return

    call udot%init(this%u)
    call this%model%eval_udot(this%u, t_init, udot, stat)
    if (stat /= 0) then
      stat = 1
      errmsg = 'failed to solve for the initial time derivative'
      return
    end if
    call this%solver%set_initial_state(t_init, this%u, udot)

  end subroutine init


  subroutine run(this, stat, errmsg)

    use,intrinsic :: iso_fortran_env, only: output_unit

    class(mfe_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, nstep
    real(r8) :: t, hnext, t_write, t_hard, t_soft
    character(16) :: string
    logical :: write_diag
    type(action_list), allocatable :: actions
    class(event_action), allocatable :: action

    call start_timer('integration')

    t = this%solver%last_time()
    this%tlast = t
    hnext = this%dt_init
    this%hlast = hnext

    call this%env%log%info('')
    write(string,'(g0.6)') t
    call this%env%log%info('Starting integration at t = ' // trim(string))

    call write_soln(this%env%grf_unit, t, this%u)
    t_write = t
    write_diag = .false.

    call this%hard_eventq%fast_forward(t)
    call this%soft_eventq%fast_forward(t)
    t_soft = this%soft_eventq%next_time()

    stat = 0
    nstep = 0
    do

      !! PER STEP ACTIONS
      if (nstep > 0) then
        if (mod(nstep,this%ofreq) == 0 .and. t_write /= t) then
          call this%solver%get_last_state_copy(this%u)  !TODO: get view?
          call write_soln(this%env%grf_unit, t, this%u)
          t_write = t
        end if
        if (write_diag) then
          !TODO: use log%info for this
          call this%solver%write_metrics(this%env%log_unit)
          call this%solver%write_metrics(output_unit)
          write_diag = .false.
        end if
      end if

      if (nstep == this%mstep) then
        stat = 0
        errmsg = 'completed maximum number of time steps'
        exit
      end if

      if (stat < 0) then
        if (this%tlast /= t_write) then
          call this%solver%get_last_state_copy(this%u)
          call write_soln(this%env%grf_unit, this%tlast, this%u)
          call this%solver%write_metrics(this%env%log_unit)
          call this%solver%write_metrics(output_unit)
        end if
        exit
      end if

      if (this%hard_eventq%is_empty() .and. this%soft_eventq%is_empty()) then
        stat = 0
        write(string,'(g0.6)') t
        call this%env%log%info('')
        call this%env%log%info('completed integration to final time t = ' // trim(string))
        exit
      end if

      !! Time for next step; nominally TLAST+HNEXT but possibly adjusted
      t_hard = this%hard_eventq%next_time()
      if (this%hard_eventq%is_empty()) then
        t = this%tlast + hnext
      else
        t = this%ts_sync%next_time(t_hard, this%tlast, this%hlast, hnext)
      end if
      nstep = nstep + 1

      call this%solver%advance(t, this%u, hnext, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'integration failure: ' // errmsg
        cycle
      end if

      this%hlast = t - this%tlast
      this%tlast = t

      do while (t >= t_soft) ! HANDLE SOFT EVENT ACTIONS AT T_SOFT
        call this%soft_eventq%pop_actions(actions)
        do
          call actions%get_next_action(action)
          if (.not.allocated(action)) exit
          select type (action)
          type is (output_event)
            call this%solver%get_interpolated_state(t_soft, this%u)
            call write_soln(this%env%grf_unit, t_soft, this%u)
            t_write = t_soft
            write_diag = .true.
          class default
            INSIST(.false.)
          end select
        end do
        t_soft = this%soft_eventq%next_time()
      end do

      if (t == t_hard) then ! HANDLE HARD EVENT ACTIONS AT T_HARD
        call this%hard_eventq%pop_actions(actions)
        do
          call actions%get_next_action(action)
          if (.not.allocated(action)) exit
          select type (action)
          type is (output_event)
            if (t /= t_write) then
              call this%solver%get_last_state_copy(this%u)
              call write_soln(this%env%grf_unit, t, this%u)
              t_write = t
              write_diag = .true.
            end if
          class default
            INSIST(.false.)
          end select
        end do
      end if

    end do

    call stop_timer('integration')

  end subroutine run

  subroutine write_soln(lun, t, u)

    integer, intent(in) :: lun
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(in) :: u

    integer :: j
    character(16) :: fmt

    associate (x => u%array(u%neqns+1,:), uu => u%array(:u%neqns,:))
      write(lun,'(a,es13.5)') 'TIME = ', t
      write(fmt,'(a,i0,a)') '(', u%neqns+1, 'es17.8)'
      do j = 1, u%nnode
        write(lun,fmt) x(j), uu(:,j)
      end do
    end associate

  end subroutine write_soln

end module mfe_sim_type
