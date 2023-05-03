!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module initialize

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants
  use problem_init
  use common_io
  implicit none
  private

  public :: read_input
  public :: refine

  integer :: ofreq, mstep, debug
  real(r8), allocatable :: tout(:)

  ! MFE ODE solver parameters.
  integer  :: mtry, mitr, mvec
  real(r8) :: h, hlb, hub, ntol, margin, vtol

contains

  subroutine read_input(params)

    use parameter_list_type

    type(parameter_list), intent(inout) :: params

    integer :: nseg, i, n, nout
    integer, allocatable  :: niseg(:)
    real(r8), allocatable :: useg(:,:)
    integer :: bc_left_u_type(NEQNS), bc_right_u_type(NEQNS), bc_left_x_type, bc_right_x_type

    integer :: kreg
    real(r8) :: eqw(NEQNS), eltvsc(NEQNS), segspr(NEQNS)
    real(r8) :: fdinc, dxmin

    integer, parameter :: FREE  = 0
    integer, parameter :: FIXED = 1

    write(log_unit,'(//a/)') 'I N I T I A L   C O N D I T I O N S'

    call read_tagged_data(nseg, 'Mesh segments (NSEG)')

    call params%set('nseg', nseg)

    allocate(niseg(nseg), useg(NVARS,nseg+1))

    call read_tagged_data(niseg, 'Elements per segment (NISEG)')
    call read_tagged_data(useg, 'Initial solution (USEG)')
    call params%set('niseg', niseg)
    call params%set('useg', useg)

    write(log_unit,*)
    call read_tagged_data(bc_left_u_type, 'Left-end boundary conditions')
    call read_tagged_data(bc_right_u_type, 'Right-end boundary conditions')
    call read_tagged_data(bc_left_x_type, 'Left-end node type')
    call read_tagged_data(bc_right_x_type, 'Right-end node type')

    call params%set('bc-left-u-type', bc_left_u_type == FIXED)
    call params%set('bc-right-u-type', bc_right_u_type == FIXED)
    call params%set('bc-left-x-type', bc_left_x_type == FIXED)
    call params%set('bc-right-x-type', bc_right_x_type == FIXED)

    write(log_unit,'(//a/)') 'M F E   P A R A M E T E R S'

    if (NEQNS > 1) then
      call read_tagged_data(eqw, 'Equation weights (EQW)')
    else
      eqw(1) = 1.0_r8
    end if
    call params%set('eqw', eqw)

    call read_tagged_data(kreg, 'Segment visc type (KREG)')
    call read_tagged_data(eltvsc, 'Segment visc coef (ELTVSC)')
    call read_tagged_data(segspr, 'Segment spring coef (SEGSPR)')
    call read_tagged_data(dxmin, 'Minimum element size (DXMIN)')
    call read_tagged_data(fdinc, 'Jacob FD increment (FDINC)')

    call params%set('kreg', kreg)
    call params%set('eltvsc', eltvsc)
    call params%set('segspr', segspr)
    call params%set('dxmin', dxmin)
    call params%set('fdinc', fdinc)

    write(log_unit,*)
    call read_tagged_data(debug, 'Solver debug level (DEBUG)')
    call read_tagged_data(ofreq, 'Output frequency (OFREQ)')
    call read_tagged_data(mstep, 'Maximum steps (MSTEP)')

    call params%set('verbose-stepping', (debug/=0))
    call params%set('ofreq', ofreq)
    call params%set('mstep', mstep)

    call read_tagged_data(nout)
    allocate(tout(nout))
    call read_tagged_data(tout, 'Output times (TOUT)')

    call params%set('tout', tout)

    write(log_unit,*)
    block
      real(r8) :: rtol, ptol(NVARS)
      call read_tagged_data(rtol, 'Relative dx tolerance (RTOL)')
      call read_tagged_data(ptol(:neqns), 'U predictor tolerance (PTOL)')
      call read_tagged_data(ptol(neqns+1), 'X predictor tolerance (PTOL)')

      call params%set('rtol', rtol)
      call params%set('ptol', ptol)
    end block

    write(log_unit,*)
    call read_tagged_data(h, 'Initial time step (H)')
    call read_tagged_data(hlb, 'Time step lower bound (HLB)')
    call read_tagged_data(hub, 'Time step upper bound (HUB)')
    call read_tagged_data(margin, 'Jacob update margin (MARGIN)')
    call read_tagged_data(ntol, 'FPI tolerance (NTOL)')
    call read_tagged_data(mitr, 'Max FP iterations (MITR)')

    call params%set('h', h)
    call params%set('hlb', hlb)
    call params%set('hub', hub)
    call params%set('margin', margin)
    call params%set('nlk-ntol', ntol)
    call params%set('nlk-max-iter', mitr)

    write(log_unit,*)
    call read_tagged_data(vtol, 'FPA vector tolerance (VTOL)')
    call read_tagged_data(mvec, 'Max FPA vectors (MVEC)')

    call params%set('nlk-vec-tol', vtol)
    call params%set('nlk-max-vec', mvec)

    write(log_unit,'(//a/)') 'P R O B L E M   S P E C I F I C   P A R A M E T E R S'

    call read_problem_data

  end subroutine read_input

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  REFINE
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine refine(x0, n, x)

    use mfe1_vector_type

    integer, intent(in) :: n(:)
    real(r8), intent(in) :: x0(:,:)
    type(mfe1_vector), intent(inout) :: x

    integer :: node, j, k
    real(r8) :: dx, du(NEQNS)

    node = 1

    associate (xx => x%array(NEQNS+1,:), u => x%array(:NEQNS,:))
      do k = 1, size(n)
        u(:,node) = x0(2:NVARS,k)
        xx(node) = x0(1,k)
        node = node + 1

        dx = (x0(1,k+1) - x0(1,k)) / n(k)
        du = (x0(2:NVARS,k+1) - x0(2:NVARS,k)) / n(k)

        do j = 1, n(k) - 1
          xx(node) = xx(node-1) + dx
          u(:,node) = u(:,node-1) + du
          node = node + 1
        end do
      end do
      u(:,node) = x0(2:NVARS,size(n)+1)
      xx(node) = x0(1,size(n)+1)
    end associate

  end subroutine refine

end module initialize

