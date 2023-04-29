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
  use mfe1_vector_type
  use mfe_data
  use bc_data
  use norm_data
  use problem_init
  use common_io
  implicit none
  private

  public  :: read_soln, read_data
  private :: refine

  integer, save :: nnod, nelt

  integer, save, public :: ofreq, mstep, debug
  real(r8), allocatable, save, public :: tout(:)

  ! MFE ODE solver parameters.
  integer, save, public  :: mtry, mitr, mvec
  real(r8), save, public :: h, hlb, hub, ntol, margin, vtol

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  READ_SOLN
  !!
  !!    This procedure should:
  !!    1) complete the initialization of the data in mfe_extents;
  !!    2) generate the initial solution u and allocate storage for udot; and
  !!    3) initialize all data in the module bc_data.
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_soln(u, udot)

    type(mfe1_vector), intent(out) :: u, udot

    integer :: nseg
    integer, allocatable  :: niseg(:)
    real(r8), allocatable :: useg(:,:)

    write(log_unit,'(//a/)') 'I N I T I A L   C O N D I T I O N S'

    call read_tagged_data(nseg, 'Mesh segments (NSEG)')

    allocate(niseg(nseg), useg(NVARS,nseg+1))

    call read_tagged_data(niseg, 'Elements per segment (NISEG)')
    call read_tagged_data(useg, 'Initial solution (USEG)')

    nelt = sum(niseg)
    nnod = nelt + 1

    call u%init(NEQNS, nnod)
    call udot%init(u)

    call refine(useg, niseg, u)

    deallocate(niseg, useg)

    write(log_unit,*)
    call read_tagged_data(bc_left%u_type, 'Left-end boundary conditions')
    call read_tagged_data(bc_right%u_type, 'Right-end boundary conditions')
    call read_tagged_data(bc_left%x_type, 'Left-end node type')
    call read_tagged_data(bc_right%x_type, 'Right-end node type')

    ! Save the initial boundary values.
    bc_left%x_value = u%array(neqns+1,1)
    bc_left%u_value = u%array(:neqns,1)
    bc_right%x_value = u%array(neqns+1,nnod)
    bc_right%u_value = u%array(:neqns,nnod)

  end subroutine read_soln

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  READ_DATA
 !!
 !!    This procedure should:
 !!    1) initialize the parameters in module mfe_data;
 !!    2) initialize the parameters in module norm_data;
 !!    3) ...
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_data

    integer :: nout

    write(log_unit,'(//a/)') 'M F E   P A R A M E T E R S'

    if (NEQNS > 1) then
      call read_tagged_data(eqw, 'Equation weights (EQW)')
    else
      eqw(1) = 1.0_r8
    end if

    call read_tagged_data(kreg, 'Segment visc type (KREG)')
    call read_tagged_data(eltvsc, 'Segment visc coef (ELTVSC)')
    call read_tagged_data(segspr, 'Segment spring coef (SEGSPR)')
    call read_tagged_data(dxmin, 'Minimum element size (DXMIN)')
    call read_tagged_data(fdinc, 'Jacob FD increment (FDINC)')

    write(log_unit,*)
    call read_tagged_data(debug, 'Solver debug level (DEBUG)')
    call read_tagged_data(ofreq, 'Output frequency (OFREQ)')
    call read_tagged_data(mstep, 'Maximum steps (MSTEP)')

    call read_tagged_data(nout)
    allocate(tout(nout))
    call read_tagged_data(tout, 'Output times (TOUT)')

    write(log_unit,*)
    call read_tagged_data(rtol, 'Relative dx tolerance (RTOL)')
    call read_tagged_data(ptol(:neqns), 'U predictor tolerance (PTOL)')
    call read_tagged_data(ptol(neqns+1), 'X predictor tolerance (PTOL)')

    write(log_unit,*)
    call read_tagged_data(h, 'Initial time step (H)')
    call read_tagged_data(hlb, 'Time step lower bound (HLB)')
    call read_tagged_data(hub, 'Time step upper bound (HUB)')
    call read_tagged_data(margin, 'Jacob update margin (MARGIN)')
    call read_tagged_data(ntol, 'FPI tolerance (NTOL)')
    call read_tagged_data(mitr, 'Max FP iterations (MITR)')

    write(log_unit,*)
    call read_tagged_data(vtol, 'FPA vector tolerance (VTOL)')
    call read_tagged_data(mvec, 'Max FPA vectors (MVEC)')

    mtry = 9

    write(log_unit,'(//a/)') 'P R O B L E M   S P E C I F I C   P A R A M E T E R S'

    call read_problem_data

  end subroutine read_data

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  REFINE
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine refine(x0, n, x)

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

        !dx%u = (x0(2:NVARS,k+1) - x0(2:NVARS,k)) / n(k)
        !dx%x = (x0(1,k+1) - x0(1,k)) / n(k)

        dx = (x0(1,k+1) - x0(1,k)) / n(k)
        du = (x0(2:NVARS,k+1) - x0(2:NVARS,k)) / n(k)

        do j = 1, n(k) - 1
          !x(node)%x = x(node-1)%x + dx%x
          !x(node)%u = x(node-1)%u + dx%u
          xx(node) = xx(node-1) + dx
          u(:,node) = u(:,node-1) + du
          node = node + 1
        end do
      end do
      !x(node)%u = x0(2:NVARS,size(n)+1)
      !x(node)%x = x0(1,size(n)+1)
      u(:,node) = x0(2:NVARS,size(n)+1)
      xx(node) = x0(1,size(n)+1)
    end associate

  end subroutine refine

end module initialize

