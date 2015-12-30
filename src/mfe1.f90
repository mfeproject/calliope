!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  MFE1 --
!!
!!  This is a Fortran 90 implementation of the piecewise-linear
!!  Gradient-Weighted Moving Finite Element Method (GWMFE) for
!!  time-dependent systems of partial differential equations in
!!  one space dimension.
!!
!!  Version 0.3, 12 August 1997
!!
!!  Neil N. Carlson, Dept of Math, Purdue University
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997, Neil N. Carlson
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

program mfe1

  use mfe_constants
  use mfe_types
  use common_io
  use output
  use initialize
  !use mfe_ode_solver
  use idaesol_type
  use mfe_model_type
  use parameter_list_type
  use mfe_procs, only: eval_udot
  use secure_hash_factory
  use,intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none

  real(kind=wp) :: t
  real(kind=wp), dimension(6) :: rvar
  integer :: mode, rtype, j, errc, debug_unit, nstep
  integer, dimension(3) :: ivar
  type(NodeVar), dimension(:), pointer :: u, udot
  real(wp), allocatable :: uflat(:), udotflat(:)
  type(idaesol) :: solver
  type(mfe_model), target :: model
  type(parameter_list) :: params
  class(secure_hash), allocatable :: hash

  character(len=16) :: string
  
  real :: cpusec, cpusec0

  call cpu_time (cpusec0)
  
  ! Open the input and output files.

  open (newunit=input_unit, file="mfein", position="rewind", action="read", status="old")

  open (newunit=log_unit, file="mfelog", position="rewind", action="write", status="replace")

  open (newunit=out_unit, file="mfegrf", position="rewind", action="write", status="replace")

  call read_soln (u, udot)
  call read_data ()

do j = 1, size(u)
  write(*,*) u(j)%u, u(j)%x
end do

  t = tout(1)
  !mode = START_SOLN
  call eval_udot (u, t, udot, errc)
  call write_soln (u, t)
  j = 2

  if (errc /= 0) then
    call abort ( (/log_unit,error_unit/), "Bad initial solution.")
  end if
  
  if (mstep <= 0) then
    stop
  end if
  call checksum_nodevar (u, 'U')
  call checksum_nodevar (udot, 'UDOT')
  allocate(uflat(NVARS*size(u)), udotflat(NVARS*size(udot)))
  call copy_from_nodevar (u, uflat)
  call copy_from_nodevar (udot, udotflat)
  
  call new_secure_hash (hash, 'md5')
  call hash%update (uflat)
  write(output_unit,'(a)') 'U: ' // hash%hexdigest()
  call hash%update (udotflat)
  write(output_unit,'(a)') 'UDOT: ' // hash%hexdigest()
  
  call params%set ('nlk-max-iter', mitr)
  call params%set ('nlk-ntol', ntol)
  call params%set ('nlk-max-vec', mvec)
  call params%set ('nlk-vec-tol', vtol)
  
  call model%init (size(u))
  call solver%init (model, params)
  call solver%set_initial_state (t, uflat, udotflat)

  if (debug /= 0) then
    open (newunit=debug_unit, file="bdfout", action="write", status="replace")
    !call set_solver_messages (debug, debug_unit)
    call solver%set_verbose_stepping (debug_unit)
  end if

  rvar(1) = h
  rvar(2) = hlb
  rvar(3) = hub
  !rvar(4) = ntol
  !rvar(5) = margin ! not used
  !rvar(6) = vtol

  ivar(1) = mtry
  !ivar(2) = mitr
  !ivar(3) = mvec

  write(unit=string, fmt="(es11.3)") t
  call info ([log_unit, output_unit], "BEGINNING SOLUTION AT T = " // string)

  if (ofreq <= 0) then
    ofreq = mstep
  end if

  do

    !call bdf2_solver (mode, rvar, ivar, tout(j), ofreq, rtype, u, udot, t)
    
    call solver%get_metrics (nstep=nstep)
    call solver%integrate (h, rtype, ofreq-modulo(nstep,ofreq), tout(j), hlb, hub, mtry)

    !mode = RESUME_SOLN

    t = solver%last_time()
    call cpu_time (cpusec)
    cpusec = cpusec - cpusec0
    !call write_soln (u, t)
        call solver%write_metrics (log_unit)
        call solver%write_metrics (output_unit)

    select case (rtype)

      !case (SOLN_AT_TOUT)       ! Integrated to TOUT.
      case (SOLVED_TO_TOUT)       ! Integrated to TOUT.
        call solver%get_interpolated_state (tout(j), uflat)
        call copy_to_nodevar (uflat, u)
        call write_soln (u, tout(j))
        !call write_status ([log_unit, output_unit], cpusec, t)
        j = j + 1

      !case (SOLN_AT_STEP)       ! Integrated OFREQ more steps.
      case (SOLVED_TO_NSTEP)       ! Integrated OFREQ more steps.
        call solver%get_last_state_copy (uflat)
        call copy_to_nodevar (uflat, u)
        call write_soln (u, t)
        !call write_status ([log_unit, output_unit], cpusec)

      !case (FAIL_ON_STEP)
      case (STEP_FAILED)
        call solver%get_last_state_copy (uflat)
        call copy_to_nodevar (uflat, u)
        call write_soln (u, t)
        !call write_status ([log_unit, output_unit], cpusec)
        call abort ([log_unit, error_unit], "Repeated failure at a step.")

      !case (SMALL_H_FAIL)
      case (STEP_SIZE_TOO_SMALL)
        call solver%get_last_state_copy (uflat)
        call copy_to_nodevar (uflat, u)
        call write_soln (u, t)
        !call write_status ([log_unit, output_unit], cpusec)
        call abort ([log_unit, error_unit], "Next time step is too small.")

      !case (FAIL_ON_START)
      case (BAD_INPUT)
        call solver%get_last_state_copy (uflat)
        call copy_to_nodevar (uflat, u)
        call write_soln (u, t)
        !call write_status ([log_unit, output_unit], cpusec)
        call abort ([log_unit, error_unit], "Bad input.")

      case default
        call solver%get_last_state_copy (uflat)
        call copy_to_nodevar (uflat, u)
        call write_soln (u, t)
        !call write_status ([log_unit, output_unit], cpusec)
        call abort ([log_unit, error_unit], "Unknown return type!")

    end select

    if (j > size(tout)) then
      call info ([log_unit, output_unit], "Integrated to final TOUT.  Done.")
      exit
    end if

    call solver%get_metrics (nstep=nstep)
    if (nstep >= mstep) then
      call info ([log_unit, output_unit], "Maximum number of steps taken.  Done.")
      exit
    end if

  end do

contains

  subroutine copy_to_nodevar (u, ustruct)
    use mfe_types, only: NodeVar
    real(wp), intent(in), target :: u(:)
    type(NodeVar), intent(out) :: ustruct(:)
    integer :: j, k
    real(wp), pointer :: u2(:,:)
    u2(1:NVARS,1:size(ustruct)) => u
    do j = 1, size(ustruct)
      do k = 1, NEQNS
        ustruct(j)%u(k) = u2(k,j)
      end do
      ustruct(j)%x = u2(k,j)
    end do
  end subroutine
  
  
  subroutine copy_from_nodevar (ustruct, u)
    use mfe_types, only: NodeVar
    type(NodeVar), intent(in) :: ustruct(:)
    real(wp), intent(out), target :: u(:)
    integer :: j, k
    real(wp), pointer :: u2(:,:)
    u2(1:NVARS,1:size(ustruct)) => u
    do j = 1, size(ustruct)
      do k = 1, NEQNS
        u2(k,j) = ustruct(j)%u(k)
      end do
      u2(k,j) = ustruct(j)%x
    end do
  end subroutine

  subroutine checksum_nodevar (u, name)
    use secure_hash_factory
    type(NodeVar), intent(in) :: u(:)
    character(*), intent(in) :: name
    integer :: j
    class(secure_hash), allocatable :: hash
    call new_secure_hash (hash, 'md5')
    do j = 1, size(u)
      call hash%update(u(j)%u)
      call hash%update(u(j)%x)
    end do
    write(output_unit,'(a)') name // ': ' // hash%hexdigest()
  end subroutine

end program mfe1
