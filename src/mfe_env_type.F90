module mfe_env_type

  use simlog_type
  implicit none
  private

  type, public :: mfe_env
    type(simlog) :: log
    integer :: log_unit, grf_unit
  end type

  public :: VERB_SILENT, VERB_NORMAL, VERB_NOISY  ! export parameters from simlog_type

end module mfe_env_type