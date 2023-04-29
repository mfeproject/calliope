!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_mfe

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants, only: NVERT, NEQNS, NVARS
  use mfe1_vector_type
  use mfe_data
  use local_arrays
  use problem_pde
  implicit none
  private

  public :: preprocessor, res_mass_matrix, eval_mass_matrix, reg_rhs, eval_dfdy

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! PREPROCESSOR --
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine preprocessor(errc)

    integer, intent(out), optional :: errc

    integer :: j
    logical :: check_dx

    if (present(errc)) then
      check_dx = .true.
      errc = 0
    else
      check_dx = .false.
    end if

!NB: This establishes the order of unknowns at a node: x, u1, u2, ...
    do j = 1, ncell

      ! CELL LENGTH
      dx(j) = u(NEQNS+1,2,j) - u(NEQNS+1,1,j)

      if (check_dx .and. dx(j) < dxmin) then
        errc = j
        return
      end if

      ! U DIFFERENCE ACROSS CELL
      du(:,j) = u(:NEQNS,2,j) - u(:NEQNS,1,j)

      ! CELL LENGTH ON SOLUTION MANIFOLD
      l(:,j) = sqrt( dx(j)**2 + du(:,j)**2 )

      ! UNIT NORMAL TO THE SOLUTION MANIFOLD
      n(1,:,j) = - du(:,j) / l(:,j)
      n(2,:,j) =     dx(j) / l(:,j)

      ! SOLUTION GRADIENT
      dudx(:,j) = du(:,j) / dx(j)

    end do

  end subroutine preprocessor

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RES_MASS_MATRIX
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine res_mass_matrix

    integer :: i, j
    real(r8) :: c(neqns)
    real(r8) :: term2, dxdot, dudot
    real(r8) :: ndot(NVERT), term(NVERT)

   !!!
   !!!  Pure MFE mass matrix terms.

    c = eqw / 6.0_r8

    associate (rx => r(ix,:,:), xdot => udot(ix,:,:))
      do j = 1, ncell
        do i = 1, neqns
          ! Normal velocity at each vertex.
          ndot(:) = n(1,i,j) * xdot(:,j) + n(2,i,j) * udot(i,:,j)
          term(:) = (c(i) * l(i,j)) * (sum(ndot) + ndot(:))

          rx(:,j)  = rx(:,j)  - term * n(1,i,j)
          r(i,:,j) = r(i,:,j) - term * n(2,i,j)
        end do
      end do

     !!!
     !!!  Regularization contribution to the mass matrix.

      select case (kreg)
      case (1)

       !!!
       !!! RATE OF DEFORMATION PENALIZATION

        c = eqw * eltvsc
        do j = 1, ncell
          dxdot = xdot(2,j) - xdot(1,j)
          do i = 1, neqns
            dudot = udot(i,2,j) - udot(i,1,j)
            term2 = (c(i) / l(i,j)) * (n(2,i,j) * dxdot - n(1,i,j) * dudot)

            rx(1,j) = rx(1,j) + (term2 * n(2,i,j))
            rx(2,j) = rx(2,j) - (term2 * n(2,i,j))

            r(i,1,j) = r(i,1,j) - (term2 * n(1,i,j))
            r(i,2,j) = r(i,2,j) + (term2 * n(1,i,j))
          end do
        end do

      case (2)

       !!!
       !!! TOTAL GRADIENT PENALIZATION

        c = eqw * eltvsc
        do j = 1, ncell
          dxdot = xdot(2,j) - xdot(1,j)
          do i = 1, neqns
            dudot = udot(i,2,j) - udot(i,1,j)

            rx(1,j) = rx(1,j) + ((c(i) / l(i,j)) * dxdot)
            rx(2,j) = rx(2,j) - ((c(i) / l(i,j)) * dxdot)

            r(i,1,j) = r(i,1,j) + ((c(i) / l(i,j)) * dudot)
            r(i,2,j) = r(i,2,j) - ((c(i) / l(i,j)) * dudot)
          end do
        end do
      end select
    end associate

  end subroutine res_mass_matrix

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EVAL_MASS_MATRIX
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eval_mass_matrix(factor, diag_only)

    real(r8), intent(in), optional :: factor
    logical, intent(in), optional :: diag_only

    integer :: i, j
    logical :: diagonal
    real(r8) :: blk(NVARS,NVARS)
    real(r8) :: fac, term, term_xx, term_xu, term_uu
    real(r8) :: c(neqns)

    if (present(factor)) then
      fac = factor
    else
      fac = 1.0_r8
    end if

    if (present(diag_only)) then
      diagonal = diag_only
    else
      diagonal = .false.
    end if

   !!!
   !!! PURE MFE MASS MATRIX

    c = (fac / 3.0_r8) * eqw
    do j = 1, ncell
      blk = 0.0_r8
      do i = 1, neqns
        associate (n1 => n(1,i,j), n2 => n(2,i,j))
          blk(ix,ix)   =  (c(i) * l(i,j)) * n1 * n1 + blk(ix,ix)
          blk(ix,i)    =  (c(i) * l(i,j)) * n1 * n2
          blk(i,ix)    =  (c(i) * l(i,j)) * n1 * n2
          blk(i,i)     =  (c(i) * l(i,j)) * n2 * n2
        end associate
      end do

      ! COPY THE BASIC blk

      mtx(:,:,1,1,j) = blk
      mtx(:,:,2,2,j) = blk

      if (diagonal) cycle

      blk = 0.5_r8 * blk

      mtx(:,:,2,1,j) = blk
      mtx(:,:,1,2,j) = blk
    end do

   !!!
   !!!  Regularization contribution to the mass matrix.

    select case (kreg)
    case (1)

     !!!
     !!! RATE OF DEFORMATION PENALIZATION

      c = fac * eqw * eltvsc
      do j = 1, ncell
        do i = 1, neqns
          associate (n1 => n(1,i,j), n2 => n(2,i,j))
            term = c(i) / l(i,j)
            term_xx =   term * n2 * n2
            term_xu = - term * n1 * n2
            term_uu =   term * n1 * n1

            mtx(ix,ix,1,1,j) = mtx(ix,ix,1,1,j) + term_xx
            mtx(ix,i,1,1,j)  = mtx(ix,i,1,1,j)  + term_xu
            mtx(i,ix,1,1,j)  = mtx(i,ix,1,1,j)  + term_xu
            mtx(i,i,1,1,j)   = mtx(i,i,1,1,j)   + term_uu

            mtx(ix,ix,2,2,j) = mtx(ix,ix,2,2,j) + term_xx
            mtx(ix,i,2,2,j)  = mtx(ix,i,2,2,j)  + term_xu
            mtx(i,ix,2,2,j)  = mtx(i,ix,2,2,j)  + term_xu
            mtx(i,i,2,2,j)   = mtx(i,i,2,2,j)   + term_uu
          end associate

          if (diagonal) then
            cycle
          end if

          mtx(ix,ix,2,1,j) = mtx(ix,ix,2,1,j) - term_xx
          mtx(ix,i,2,1,j)  = mtx(ix,i,2,1,j)  - term_xu
          mtx(i,ix,2,1,j)  = mtx(i,ix,2,1,j)  - term_xu
          mtx(i,i,2,1,j)   = mtx(i,i,2,1,j)   - term_uu

          mtx(ix,ix,1,2,j) = mtx(ix,ix,1,2,j) - term_xx
          mtx(ix,i,1,2,j)  = mtx(ix,i,1,2,j)  - term_xu
          mtx(i,ix,1,2,j)  = mtx(i,ix,1,2,j)  - term_xu
          mtx(i,i,1,2,j)   = mtx(i,i,1,2,j)   - term_uu
        end do
      end do

    case (2)

     !!!
     !!! TOTAL GRADIENT PENALIZATION

      c = fac * eqw * eltvsc
      do j = 1, ncell
        do i = 1, neqns
          term = c(i) / l(i,j)

          mtx(ix,ix,1,1,j) = mtx(ix,ix,1,1,j) + term
          mtx(i,i,1,1,j)   = mtx(i,i,1,1,j)   + term

          mtx(ix,ix,2,2,j) = mtx(ix,ix,2,2,j) + term
          mtx(i,i,2,2,j)   = mtx(i,i,2,2,j)   + term

          if (diagonal) then
            cycle
          end if

          mtx(ix,ix,2,1,j) = mtx(ix,ix,2,1,j) - term
          mtx(i,i,2,1,j)   = mtx(i,i,2,1,j)   - term

          mtx(ix,ix,1,2,j) = mtx(ix,ix,1,2,j) - term
          mtx(i,i,1,2,j)   = mtx(i,i,1,2,j)   - term
        end do
      end do
    end select

  end subroutine eval_mass_matrix

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! REG_RHS --
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine reg_rhs

    integer :: i, j
    real(r8) :: term, term_x, term_u
    real(r8) :: c(neqns)

    associate (rx => r(ix,:,:))
      c = eqw * segspr
      do j = 1, ncell
        do i = 1, neqns
          term = c(i) / l(i,j)**2
          term_x =   term * n(2,i,j)
          term_u = - term * n(1,i,j)

          rx(1,j)  = rx(1,j)  - term_x
          r(i,1,j) = r(i,1,j) - term_u

          rx(2,j)  = rx(2,j)  + term_x
          r(i,2,j) = r(i,2,j) + term_u
        end do
      end do
    end associate

  end subroutine reg_rhs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EVAL_DFDY
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eval_dfdy(t, errc)

    real(r8), intent(in) :: t
    integer, intent(out) :: errc

    integer :: i, j, k, l
    real(r8) :: rh
    real(r8), allocatable :: r_save(:,:,:)

    rh = 1.0_r8 / fdinc

   !!!
   !!!  Compute and save the unperturbed residual.

    call pde_rhs(t)
    call reg_rhs
    call res_mass_matrix

    r_save = r

    do k = 1, NVERT

      do i = 1, neqns+1
        u(i,k,:) = u(i,k,:) + fdinc
        call preprocessor(errc)
        if (errc /= 0) return
        call pde_rhs(t)
        call reg_rhs
        call res_mass_matrix
        u(i,k,:) = u(i,k,:) - fdinc

        do j = 1, ncell
          do l = 1, NVERT
            mtx(:,i,l,k,j) = mtx(:,i,l,k,j) + rh * (r(:,l,j) - r_save(:,l,j))
          end do
        end do
      end do

    end do

    errc = 0

  end subroutine eval_dfdy

end module local_mfe

