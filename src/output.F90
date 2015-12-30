!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module output

  use mfe_constants
  use mfe_types
  use common_io
  implicit none
  private

  !public :: write_status, write_soln, write_vels
  public :: write_soln, write_vels

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  WRITE_STATUS
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    subroutine write_status (units, cpusec, time)
!    
!      use mfe_ode_solver, only: bdf2_inquire
!    
!      integer, dimension(:), intent(in) :: units
!      real, intent(in) :: cpusec
!      real(kind=wp), intent(in), optional :: time
!      
!      real(kind=wp) :: t, h
!      integer :: j
!      integer, dimension(5) :: n
!      
!      character(len=*), parameter :: fmt1 = "(/, i5, a, es11.4, a, es10.3, a, es10.3, a)"
!      character(len=*), parameter :: fmt2 = "(t8, a, i5.5, ':', i4.4)"
!      character(len=*), parameter :: fmt3 = "(a, 5(i3.3, :, ':'))"
!      
!      if(present(time)) then
!        call bdf2_inquire( nstep=n(1), h_last=h )
!        do j = 1, size(units)
!          write( unit=units(j), fmt=fmt1 ) n(1), ": T =", time, ", H =", h, ", CPU =", cpusec, " SEC"
!        end do
!        
!      else
!        call bdf2_inquire( nstep=n(1), t=t, h_last=h )
!        do j = 1, size(units)
!          write( unit=units(j), fmt=fmt1 ) n(1), ": T =", t, ", H =", h, ", CPU =", cpusec, " SEC"
!        end do
!      end if
!      
!      call bdf2_inquire( nres=n(1), njac=n(2) )
!      do j = 1, size(units)
!        write( unit=units(j), fmt=fmt2, advance="no" ) "NRES:NJAC = ", n(1:2)
!      end do
!      
!      call bdf2_inquire( nbp=n(1), njf=n(2), nnr=n(3), nnf=n(4), nsr=n(5) )
!      do j = 1, size(units)
!        write( unit=units(j), fmt=fmt3 ) ",  NBP:NJF:NNR:NNF:NSR = ", n(1:5)
!      end do
!      
!    end subroutine write_status
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  WRITE_SOLN
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine write_soln (u, t)

      type(NodeVar), dimension(:), intent(in) :: u
      real(kind=wp), intent(in) :: t

      integer :: j
      character(len=16) :: fmt

      write (unit=out_unit, fmt="(a,es13.5)") "TIME = ", t
      write (unit=fmt,fmt="(a,i2,a)") "(", NVARS, "es17.8)"
      do j = 1, size(u)
        write (unit=out_unit,fmt=fmt) u(j) % x, u(j) % u
      end do

    end subroutine write_soln

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  WRITE_VELS
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine write_vels (unit, u, udot, t)

      integer, intent(in) :: unit
      type(NodeVar), dimension(:), intent(in) :: u, udot
      real(kind=wp), intent(in) :: t

      integer :: j
      character(len=16) :: fmt

      write (unit=unit, fmt="(a,es13.5)") "TIME = ", t
      write (unit=fmt,fmt="(a,i2,a)") "(", 1+NVARS, "es17.8)"
      do j = 1, size(u)
        write (unit=unit,fmt=fmt) u(j) % x, udot(j) % x, udot(j) % u
      end do

    end subroutine write_vels

end module output
