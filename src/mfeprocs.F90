!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_procs

  use mfe_constants, only: wp
  use mfe_types
  use local_mfe_procs
  use bc_procs
  use bls_solver
  use common_io, only: element_info
  implicit none
  private

  ! Publically accessible procedures.
  public :: eval_residual, eval_jacobian, apply_inverse_jacobian, eval_udot

  ! Storage for the Jacobian matrix
  type(NodeMtx), dimension(:), allocatable, save :: jac_lowr, jac_diag, jac_uppr

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_RESIDUAL
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_residual (u, udot, t, r)

      type(NodeVar), dimension(:), intent(in)  :: u, udot
      type(NodeVar), dimension(:), intent(out) :: r
      real(kind=wp), intent(in) :: t

      ! Local variables.
      type(NodeMtx), dimension(size(u)) :: diag

     !!!
     !!! EVALUATE THE RESIDUAL

      call gather_local_solution (u, udot)

      call preprocessor ()

      call pde_rhs (t)
      call reg_rhs ()
      call res_mass_matrix ()
      call assemble_vector (r)
      call force_bc (u, r)

     !!!
     !!! DIAGONAL PRECONDITIONING

      call eval_mass_matrix (diag_only=.true.)
      call assemble_diagonal (diag)
      call force_bc (diag)

      call factor (diag)
      call solve (diag, r)

      call free_local_arrays ()
      !call apply_inverse_jacobian (r)

    end subroutine eval_residual

    subroutine apply_inverse_jacobian (r)

      type(NodeVar), dimension(:), intent(inout) :: r

      call solve (jac_lowr, jac_diag, jac_uppr, r)

    end subroutine apply_inverse_jacobian

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_JACOBIAN
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_jacobian (u, udot, t, h, errc)

      type(NodeVar), dimension(:), intent(in) :: u, udot
      real(kind=wp), intent(in) :: t, h
      integer, intent(out) :: errc

      ! local variables.
      type(NodeMtx), dimension(size(u)) :: diag

      call gather_local_solution (u, udot)

      call preprocessor ()

      call eval_mass_matrix (factor = -1.0_wp / h)

      ! Capture the unscaled diagonal for preconditioning.
      call assemble_diagonal (diag)
      diag = (-h) * diag
      call force_bc (diag)

      call eval_dfdy (t, errc)

      if (errc /= 0) then
        call element_info ("EVAL_DFDY: BAD ELEMENT", errc)
        call free_local_arrays ()
        errc = 1
        return
      end if

      if (.not. allocated(jac_lowr)) then
        allocate (jac_lowr(size(u)), jac_diag(size(u)), jac_uppr(size(u)))
      end if

      call assemble_matrix (jac_lowr, jac_diag, jac_uppr)
      call force_bc (jac_lowr, jac_diag, jac_uppr)

      ! Diagonal preconditioning.
      call factor (diag)
      call solve (diag, jac_lowr)
      call solve (diag, jac_diag)
      call solve (diag, jac_uppr)

      ! Factorize the Jacobian.
      call factor (jac_lowr, jac_diag, jac_uppr)

      call free_local_arrays ()

      errc = 0

    end subroutine eval_jacobian

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  EVAL_UDOT
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_udot (u, t, udot, errc)

      type(NodeVar), dimension(:), intent(in)  :: u
      type(NodeVar), dimension(:), intent(out) :: udot
      real(kind=wp), intent(in) :: t
      integer, intent(out) :: errc

      call gather_local_solution (u)

      call preprocessor (errc)

      if (errc /= 0) then
        call element_info ("EVAL_UDOT: BAD ELEMENT", errc)
        call free_local_arrays ()
        errc = 1
        return
      end if

      call eval_mass_matrix ()

      if (.not. allocated(jac_lowr)) then
        allocate (jac_lowr(size(u)), jac_diag(size(u)), jac_uppr(size(u)))
      end if

      call assemble_matrix (jac_lowr, jac_diag, jac_uppr)

      call pde_rhs (t)
      call reg_rhs ()

      ! Assemble the RHS of the MFE equations...
      call assemble_vector (udot)
      
      call force_bc( jac_lowr, jac_diag, jac_uppr, udot )

      ! Solve for du/dt...
      call factor (jac_lowr, jac_diag, jac_uppr)
      call solve  (jac_lowr, jac_diag, jac_uppr, udot)

      deallocate (jac_lowr, jac_diag, jac_uppr)
      call free_local_arrays ()

      errc = 0

    end subroutine eval_udot

end module mfe_procs
