module mfe_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_model_type
  use mfe1_solver_type
  use mfe1_vector_type
  use timer_tree_type
  implicit none
  private

  type, public :: mfe_sim
    type(mfe_model) :: model
    type(mfe1_solver) :: solver
    type(mfe1_vector) :: u
    real(r8), allocatable :: tout(:)
    real(r8) :: tlast, hlast
    real(r8) :: dt_init, hlb, hub
    integer :: mtry, ofreq, mstep
  contains
    procedure :: init
    procedure :: run
  end type

contains

  subroutine init(this, neqns, params, stat, errmsg)

    use parameter_list_type

    class(mfe_sim), intent(out), target :: this
    integer, intent(in) :: neqns
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: nnode
    type(mfe1_vector) :: udot

    !! Generate the initial solution (includes the mesh)
    block
      use initialize, only: refine  !TODO: move into this module
      integer,  allocatable :: niseg(:)
      real(r8), allocatable :: useg(:,:)
      call params%get('niseg', niseg, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (size(niseg) < 1) then
        stat = 1
        errmsg = 'niseg must be assigned at least 1 value'
        return
      else if (any(niseg < 1)) then
        stat = 1
        errmsg = 'niseg values must be > 0'
        return
      end if
      call params%get('useg', useg, stat=stat, errmsg=errmsg)
      if (size(useg,1) /= neqns+1) then
        stat = 1
        errmsg = 'wrong number of columns for useg'
        return
      else if (size(useg,2) /= size(niseg)+1) then
        stat = 1
        errmsg = 'wrong number of rows for useg'
        return
      else if (any(useg(1,2:) <= useg(1,1:size(niseg)))) then
        stat = 1
        errmsg = 'useg x coordinates not strictly increasing'
        return
      end if
      nnode = sum(niseg) + 1
      call this%u%init(neqns, nnode)
      call refine(useg, niseg, this%u)
    end block

    call params%get('tout', this%tout, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (size(this%tout) < 2) then
      stat = 1
      errmsg = 'tout must be assigned at least 2 values'
      return
    else if (any(this%tout(2:) <= this%tout(1:size(this%tout)-1))) then
      stat = 1
      errmsg = 'tout values must be strictly increasing'
      return
    end if

    call this%model%init(neqns, nnode, params, stat, errmsg)
    if (stat /= 0) return
    call this%model%set_boundary_values(this%u) ! get BV values from the initial solution

    call this%solver%init(this%model, params, stat, errmsg)
    if (stat /= 0) return

    call udot%init(this%u)
    call this%model%eval_udot(this%u, this%tout(1), udot, stat)
    if (stat /= 0) then
      stat = 1
      errmsg = 'failed to solve for the initial time derivative'
      return
    end if
    call this%solver%set_initial_state(this%tout(1), this%u, udot)

    call params%get('mstep', this%mstep, default=huge(this%mstep), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%mstep < 0) then
      stat = 1
      errmsg = 'mstep must be >= 0'
      return
    end if
    call params%get('ofreq', this%ofreq, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%ofreq <= 0) this%ofreq = this%mstep

    call params%get('hlb', this%hlb, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%hlb <= 0) then
      stat = 1
      errmsg = 'hlb must be > 0.0'
      return
    end if
    call params%get('hub', this%hub, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%hub <= this%hlb) then
      stat = 1
      errmsg = 'hub must be > hlb'
      return
    end if
    call params%get('h', this%dt_init, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%dt_init < this%hlb) then
      stat = 1
      errmsg = 'h must be >= hlb'
      return
    else if (this%dt_init > this%hub) then
      stat = 1
      errmsg = 'h must be <= hub'
      return
    end if
    call params%get('mtry', this%mtry, default=9, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%mtry < 1) then
      stat = 1
      errmsg = 'mtry must be >= 1'
      return
    end if

  end subroutine init


  subroutine run(this) !, stat, errmsg)

    use output, only: write_soln ! TEMPORARILY
    use common_io ! TEMPORARILY
    use,intrinsic :: iso_fortran_env, only: error_unit, output_unit

    class(mfe_sim), intent(inout) :: this
    !integer, intent(out) :: stat
    !character(:), allocatable, intent(out) :: errmsg

    integer :: j, nstep, stat
    real(r8) :: t, hnext
    character(16) :: string
    character(:), allocatable :: errmsg

    t = this%solver%last_time()
    this%tlast = t
    hnext = this%dt_init
    this%hlast = hnext

    write(string,'(es11.3)') t
    call info([log_unit, output_unit], 'BEGINNING SOLUTION AT T = ' // string)

    call write_soln(this%u, t)

    j = 2
    do

      call start_timer('integration')
      call this%solver%get_metrics(nstep=nstep)
      !call this%solver%integrate(hnext, rtype, this%ofreq-modulo(nstep,this%ofreq), &
      !    this%tout(j), this%hlb, this%hub, this%mtry)
      call integrate(this, this%ofreq-modulo(nstep,this%ofreq), this%tout(j), hnext, stat, errmsg)
      call stop_timer('integration')

      call start_timer('output')
      t = this%solver%last_time()
      call this%solver%write_metrics(log_unit)
      call this%solver%write_metrics(output_unit)

      select case (stat)
      case (SOLVED_TO_TOUT)       ! Integrated to TOUT.
        call this%solver%get_interpolated_state(this%tout(j), this%u)
        call write_soln(this%u, this%tout(j))
        j = j + 1

      case (SOLVED_TO_NSTEP)       ! Integrated OFREQ more steps.
        call this%solver%get_last_state_copy(this%u)
        call write_soln(this%u, t)

      case (STEP_FAILED, STEP_SIZE_TOO_SMALL)
        call this%solver%get_last_state_copy(this%u)
        call write_soln(this%u, t)
        call abort([log_unit, error_unit], errmsg)

      case default
        call this%solver%get_last_state_copy(this%u)
        call write_soln(this%u, t)
        call abort([log_unit, error_unit], 'Unknown return type!')

      end select
      call stop_timer('output')

      if (j > size(this%tout)) then
        call info([log_unit, output_unit], 'Integrated to final TOUT.  Done.')
        exit
      end if

      call this%solver%get_metrics(nstep=nstep)
      if (nstep >= this%mstep) then
        call info([log_unit, output_unit], 'Maximum number of steps taken.  Done.')
        exit
      end if

    end do

  end subroutine run

!TODO: handle ofreq output in integrate

  subroutine integrate(this, nstep, tout, hnext, stat, errmsg)

    class(mfe_sim), intent(inout) :: this
    integer, intent(in) :: nstep
    real(r8), intent(in) :: tout
    real(r8), intent(inout) :: hnext
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: step
    real(r8) :: t

    step = 0
    do

      !! Check for a normal return !TODO: move to end of loop?
      if (tout <= this%tlast) then
        stat = SOLVED_TO_TOUT
        return
      else if (step >= nstep) then
        stat = SOLVED_TO_NSTEP
        return
      end if

      step = step + 1

      t = this%tlast + hnext
      call this%solver%advance(this%hlb, this%mtry, t, this%u, hnext, stat)
      select case (stat)
      case (STEP_SIZE_TOO_SMALL)
        errmsg = 'next time step size is too small'
        return
      case (STEP_FAILED)
        errmsg = 'time step failed'
        return
      end select
      !TODO: idaessol should return the message string itself

      this%hlast = t - this%tlast
      this%tlast = t
      hnext = min(hnext, this%hub)

    end do

  end subroutine

end module mfe_sim_type
