!!
!! SIMLOG_TYPE
!!
!! This module provides a common set of procedures for writing log/trace
!! messages -- the types of messages typically written to the terminal or
!! simulation log file.  All application code should use the facilities
!! provided here, and eschew use of raw Fortran writes for this purpose.
!!
!! Note: See truchas-pbf for a mpi parallel-aware version
!!

module simlog_type

  implicit none
  private

  type, public :: simlog
    private
    integer, allocatable :: units(:)
    integer :: verbosity
  contains
    procedure :: init
    procedure :: info
    procedure :: warn
    procedure :: error
    procedure :: fatal
    procedure :: close
    procedure :: exit
  end type simlog

  !! Verbosity level named constants.
  integer, parameter, public :: VERB_SILENT = 0
  integer, parameter, public :: VERB_NORMAL = 1
  integer, parameter, public :: VERB_NOISY  = 2

contains

  subroutine init(this, units, verbosity)
    class(simlog), intent(out) :: this
    integer, intent(in) :: units(:)
    integer, intent(in), optional :: verbosity
    this%units = units
    if (present(verbosity)) then
      this%verbosity = verbosity
    else
      this%verbosity = VERB_NORMAL
    end if
  end subroutine

  subroutine close(this)
    use,intrinsic :: iso_fortran_env, only: output_unit, error_unit
    class(simlog), intent(in) :: this
    character(:), allocatable :: date, time, zone
    integer :: n
    call timestamp(date, time, zone)
    call this%info('Normal terminatation on ' // date // ' at '// time // ' ' // zone)
    do n = 1, size(this%units)
      if (all(this%units(n) /= [error_unit, output_unit])) close(this%units(n))
    end do
  end subroutine

  subroutine info(this, message, verbosity)
    class(simlog), intent(in) :: this
    character(*), intent(in) :: message
    integer, intent(in), optional :: verbosity
    integer :: n, level
    if (present(verbosity)) then
      level = max(1, verbosity)
    else
      level = VERB_NORMAL
    end if
    if (level <= this%verbosity) then
      do n = 1, size(this%units)
        write(this%units(n),'(a)') message(:len_trim(message))
      end do
    end if
  end subroutine

  subroutine warn(this, message)
    class(simlog), intent(in) :: this
    character(*),  intent(in) :: message
    call labeled_message(this, 'Warning: ', message)
  end subroutine

  subroutine error(this, message)
    class(simlog), intent(in) :: this
    character(*),  intent(in) :: message
    call labeled_message(this, 'ERROR: ', message)
  end subroutine

  subroutine labeled_message(this, label, message)
    class(simlog), intent(in) :: this
    character(*),  intent(in) :: label, message
    integer :: n
    do n = 1, size(this%units)
      write(this%units(n),'(2a)') label, message(:len_trim(message))
    end do
  end subroutine

  subroutine exit(this)
    class(simlog), intent(in) :: this
    character(:), allocatable :: date, time, zone
    call timestamp(date, time, zone)
    call info(this, 'Normal terminatation on ' // date // ' at '// time // ' ' // zone)
    stop
  end subroutine

  subroutine fatal(this, message)
    class(simlog), intent(in) :: this
    character(*), intent(in) :: message
    character(:), allocatable :: date, time, zone
    call labeled_message(this, 'FATAL: ', message)
    call timestamp(date, time, zone)
    call info(this, 'Abnormal termination on ' // date // ' at ' // time // ' ' // zone)
    stop 1
  end subroutine

  subroutine timestamp(date, time, zone)
    character(:), allocatable, intent(out) :: date, time, zone
    character(8)  :: d
    character(10) :: t
    character(5)  :: z
    call date_and_time(date=d, time=t, zone=z)
    date = d(1:4) // '-' // d(5:6) // '-' // d(7:8)
    time = t(1:2) // ':' // t(3:4) // ':' // t(5:6)
    zone = z(1:5)
  end subroutine

end module simlog_type
