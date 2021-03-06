!!
!! KINDS
!!
!! This module provides kind parameter values for integer and real types
!! that are identified by byte sizes.
!!

module kinds

  public

  integer, parameter :: r4 = selected_real_kind(6)   ! 4-byte IEEE float
  integer, parameter :: r8 = selected_real_kind(15)  ! 8-byte IEEE float

  integer, parameter :: i4 = selected_int_kind(9)    ! 4-byte IEEE integer
  integer, parameter :: i8 = selected_int_kind(18)   ! 8-byte IEEE integer
  
end module kinds
