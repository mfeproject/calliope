!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module output

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  !use mfe_constants
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
!      integer, intent(in) :: units(:)
!      real, intent(in) :: cpusec
!      real(r8), intent(in), optional :: time
!      
!      real(r8) :: t, h
!      integer :: j
!      integer :: n(5)
!      
!      character(*), parameter :: fmt1 = "(/, i5, a, es11.4, a, es10.3, a, es10.3, a)"
!      character(*), parameter :: fmt2 = "(t8, a, i5.5, ':', i4.4)"
!      character(*), parameter :: fmt3 = "(a, 5(i3.3, :, ':'))"
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

    type(NodeVar(*)), intent(in) :: u(:)
    real(r8), intent(in) :: t

    integer :: j
    character(16) :: fmt

    write(out_unit,'(a,es13.5)') 'TIME = ', t
    write(fmt,'(a,i0,a)') '(', 1+u%npde, 'es17.8)'
    do j = 1, size(u)
      write(out_unit,fmt) u(j)%x, u(j)%u
    end do

  end subroutine write_soln

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  WRITE_VELS
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_vels(unit, u, udot, t)

    integer, intent(in) :: unit
    type(NodeVar(*)), intent(in) :: u(:), udot(:)
    real(r8), intent(in) :: t

    integer :: j
    character(16) :: fmt

    write(unit, '(a,es13.5)') 'TIME = ', t
    write(fmt,'(a,i0,a)') '(', 2+u%npde, 'es17.8)'
    do j = 1, size(u)
      write(unit,fmt) u(j)%x, udot(j)%x, udot(j)%u
    end do

  end subroutine write_vels

end module output
