!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_procs

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants, only: NVARS
  use mfe1_vector_type
  use local_mfe_procs
  use bc_procs
  use block_linear_solver
  use common_io, only: element_info
  use timer_tree_type
  implicit none
  private

  ! Publically accessible procedures.
  public :: eval_residual, eval_jacobian, apply_inverse_jacobian, eval_udot

  ! Storage for the Jacobian matrix
  real(r8), allocatable :: jac_lowr(:,:,:), jac_diag(:,:,:), jac_uppr(:,:,:)

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  EVAL_RESIDUAL
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eval_residual(u, udot, t, r)

    type(mfe1_vector), intent(in) :: u, udot
    type(mfe1_vector), intent(inout) :: r
    real(r8), intent(in) :: t

    real(r8) :: diag(NVARS,NVARS,u%nnode)

   !!!
   !!! EVALUATE THE RESIDUAL

    call gather_local_solution(u, udot)

    call preprocessor

    call pde_rhs(t)
    call reg_rhs
    call res_mass_matrix
    call assemble_vector(r)
    call force_bc(u, r)

   !!!
   !!! DIAGONAL PRECONDITIONING

    call eval_mass_matrix(diag_only=.true.)
    call assemble_diagonal(diag)
    call force_bc(diag)

    call vfct(diag)
    call vslv(diag, r%array)

    call free_local_arrays
    !call apply_inverse_jacobian(r)

  end subroutine eval_residual


  subroutine apply_inverse_jacobian(r)
    type(mfe1_vector), intent(inout) :: r
    call btslv(jac_lowr, jac_diag, jac_uppr, r%array)
  end subroutine apply_inverse_jacobian

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  EVAL_JACOBIAN
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eval_jacobian(u, udot, t, h, errc)

    type(mfe1_vector), intent(in) :: u, udot
    real(r8), intent(in) :: t, h
    integer, intent(out) :: errc

    ! local variables.
    real(r8) :: diag(NVARS,NVARS,u%nnode)

    call start_timer('preprocessing')
    call gather_local_solution(u, udot)

    call preprocessor
    call stop_timer('preprocessing')

    call start_timer('mass-matrix')
    call eval_mass_matrix(factor = -1.0_r8 / h)
    call stop_timer('mass-matrix')

    ! Capture the unscaled diagonal for preconditioning.
    call start_timer('diag-pc')
    call assemble_diagonal(diag)
    diag = (-h) * diag
    call force_bc(diag)
    call stop_timer('diag-pc')

    call start_timer('eval_dfdy')
    call eval_dfdy(t, errc)
    call stop_timer('eval_dfdy')

    if (errc /= 0) then
      call element_info('EVAL_DFDY: BAD ELEMENT', errc)
      call free_local_arrays
      errc = 1
      return
    end if

    if (.not.allocated(jac_lowr)) allocate(jac_lowr, jac_diag, jac_uppr, mold=diag)

    call start_timer('assembly')
    call assemble_matrix(jac_lowr, jac_diag, jac_uppr)
    call force_bc(jac_lowr, jac_diag, jac_uppr)
    call stop_timer('assembly')

    ! Diagonal preconditioning.
    call start_timer('diag-pc')
    call vfct(diag)
    call vmslv(diag, jac_lowr)
    call vmslv(diag, jac_diag)
    call vmslv(diag, jac_uppr)
    call stop_timer('diag-pc')

    ! Factorize the Jacobian.
    call start_timer('factorization')
    call btfct(jac_lowr, jac_diag, jac_uppr)
    call stop_timer('factorization')

    call free_local_arrays

    errc = 0

  end subroutine eval_jacobian

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  EVAL_UDOT
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eval_udot(u, t, udot, errc)

    type(mfe1_vector), intent(in) :: u
    type(mfe1_vector), intent(inout) :: udot
    real(r8), intent(in) :: t
    integer, intent(out) :: errc

    call gather_local_solution(u)

    call preprocessor(errc)

    if (errc /= 0) then
      call element_info('EVAL_UDOT: BAD ELEMENT', errc)
      call free_local_arrays
      errc = 1
      return
    end if

    call eval_mass_matrix

    if (.not.allocated(jac_lowr)) then
      allocate(jac_lowr(NVARS,NVARS,u%nnode), jac_diag(NVARS,NVARS,u%nnode), jac_uppr(NVARS,NVARS,u%nnode))
    end if

    call assemble_matrix(jac_lowr, jac_diag, jac_uppr)

    call pde_rhs(t)
    call reg_rhs

    ! Assemble the RHS of the MFE equations...
    call assemble_vector(udot)

    call force_bc(jac_lowr, jac_diag, jac_uppr, udot)

    ! Solve for du/dt...
    call btfct(jac_lowr, jac_diag, jac_uppr)
    call btslv(jac_lowr, jac_diag, jac_uppr, udot%array)

    deallocate(jac_lowr, jac_diag, jac_uppr)
    call free_local_arrays

    errc = 0

  end subroutine eval_udot

end module mfe_procs
