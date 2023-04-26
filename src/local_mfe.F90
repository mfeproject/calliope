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
  use mfe_constants
  use mfe_types, only: NodeVar, NodeMtx
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

    do j = 1, ncell

      ! CELL LENGTH
      dx(j) = u(2,j)%x - u(1,j)%x

      if (check_dx .and. dx(j) < dxmin) then
        errc = j
        return
      end if

      ! U DIFFERENCE ACROSS CELL
      du(:,j) = u(2,j)%u - u(1,j)%u

      ! CELL LENGTH ON SOLUTION MANIFOLD
      l(:,j) = sqrt( dx(j)**2 + du(:,j)**2 )

      ! UNIT NORMAL TO THE SOLUTION MANIFOLD
      n(:,j)%x = - du(:,j) / l(:,j)
      n(:,j)%u =     dx(j) / l(:,j)

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
    real(r8) :: c(NEQNS)
    real(r8) :: term2, dxdot, dudot
    real(r8) :: ndot(NVERT), term(NVERT)

   !!!
   !!!  Pure MFE mass matrix terms.

    c = eqw / 6.0_r8

    do j = 1, ncell

      do i = 1, NEQNS

        ! Normal velocity at each vertex.
        ndot(:) = n(i,j)%x * udot(:,j)%x + n(i,j)%u * udot(:,j)%u(i)
        term(:) = (c(i) * l(i,j)) * (sum(ndot) + ndot(:))

        r(:,j)%x    = r(:,j)%x    - term * n(i,j)%x
        r(:,j)%u(i) = r(:,j)%u(i) - term * n(i,j)%u

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

        dxdot = udot(2,j)%x - udot(1,j)%x

        do i = 1, NEQNS

          dudot = udot(2,j)%u(i) - udot(1,j)%u(i)
          term2 = (c(i) / l(i,j)) * (n(i,j)%u * dxdot - n(i,j)%x * dudot)

          r(1,j)%x    = r(1,j)%x    + (term2 * n(i,j)%u)
          r(2,j)%x    = r(2,j)%x    - (term2 * n(i,j)%u)

          r(1,j)%u(i) = r(1,j)%u(i) - (term2 * n(i,j)%x)
          r(2,j)%u(i) = r(2,j)%u(i) + (term2 * n(i,j)%x)

        end do

      end do

    case (2)

     !!!
     !!! TOTAL GRADIENT PENALIZATION

      c = eqw * eltvsc

      do j = 1, ncell

        dxdot = udot(2,j)%x - udot(1,j)%x

        do i = 1, NEQNS

          dudot = udot(2,j)%u(i) - udot(1,j)%u(i)

          r(1,j)%x    = r(1,j)%x    + ((c(i) / l(i,j)) * dxdot)
          r(2,j)%x    = r(2,j)%x    - ((c(i) / l(i,j)) * dxdot)

          r(1,j)%u(i) = r(1,j)%u(i) + ((c(i) / l(i,j)) * dudot)
          r(2,j)%u(i) = r(2,j)%u(i) - ((c(i) / l(i,j)) * dudot)

        end do

      end do

    end select

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
    type(NodeMtx) :: block
    real(r8) :: fac, term, term_xx, term_xu, term_uu
    real(r8) :: c(NEQNS)

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

      block%xx = 0.0_r8
      block%uu = 0.0_r8

      do i = 1, NEQNS

        block%xx      =  (c(i) * l(i,j)) * n(i,j)%x * n(i,j)%x + block%xx
        block%xu(i)   =  (c(i) * l(i,j)) * n(i,j)%x * n(i,j)%u
        block%ux(i)   =  (c(i) * l(i,j)) * n(i,j)%x * n(i,j)%u
        block%uu(i,i) =  (c(i) * l(i,j)) * n(i,j)%u * n(i,j)%u

      end do

      ! COPY THE BASIC BLOCK

      mtx(1,1,j) = block
      mtx(2,2,j) = block

      if (diagonal) cycle

      block%xx = 0.5_r8 * block%xx
      block%xu = 0.5_r8 * block%xu
      block%ux = 0.5_r8 * block%ux
      block%uu = 0.5_r8 * block%uu

      mtx(2,1,j) = block
      mtx(1,2,j) = block

    end do

   !!!
   !!!  Regularization contribution to the mass matrix.

    select case (kreg)
    case (1)

     !!!
     !!! RATE OF DEFORMATION PENALIZATION

      c = fac * eqw * eltvsc

      do j = 1, ncell

        do i = 1, NEQNS

          term = c(i) / l(i,j)
          term_xx =   term * n(i,j)%u * n(i,j)%u
          term_xu = - term * n(i,j)%x * n(i,j)%u
          term_uu =   term * n(i,j)%x * n(i,j)%x

          mtx(1,1,j)%xx      = mtx(1,1,j)%xx      + term_xx
          mtx(1,1,j)%xu(i)   = mtx(1,1,j)%xu(i)   + term_xu
          mtx(1,1,j)%ux(i)   = mtx(1,1,j)%ux(i)   + term_xu
          mtx(1,1,j)%uu(i,i) = mtx(1,1,j)%uu(i,i) + term_uu

          mtx(2,2,j)%xx      = mtx(2,2,j)%xx      + term_xx
          mtx(2,2,j)%xu(i)   = mtx(2,2,j)%xu(i)   + term_xu
          mtx(2,2,j)%ux(i)   = mtx(2,2,j)%ux(i)   + term_xu
          mtx(2,2,j)%uu(i,i) = mtx(2,2,j)%uu(i,i) + term_uu

          if (diagonal) then
            cycle
          end if

          mtx(2,1,j)%xx      = mtx(2,1,j)%xx      - term_xx
          mtx(2,1,j)%xu(i)   = mtx(2,1,j)%xu(i)   - term_xu
          mtx(2,1,j)%ux(i)   = mtx(2,1,j)%ux(i)   - term_xu
          mtx(2,1,j)%uu(i,i) = mtx(2,1,j)%uu(i,i) - term_uu

          mtx(1,2,j)%xx      = mtx(1,2,j)%xx      - term_xx
          mtx(1,2,j)%xu(i)   = mtx(1,2,j)%xu(i)   - term_xu
          mtx(1,2,j)%ux(i)   = mtx(1,2,j)%ux(i)   - term_xu
          mtx(1,2,j)%uu(i,i) = mtx(1,2,j)%uu(i,i) - term_uu

        end do

      end do

    case (2)

     !!!
     !!! TOTAL GRADIENT PENALIZATION

      c = fac * eqw * eltvsc

      do j = 1, ncell

        do i = 1, NEQNS

          term = c(i) / l(i,j)

          mtx(1,1,j)%xx      = mtx(1,1,j)%xx      + term
          mtx(1,1,j)%uu(i,i) = mtx(1,1,j)%uu(i,i) + term

          mtx(2,2,j)%xx      = mtx(2,2,j)%xx      + term
          mtx(2,2,j)%uu(i,i) = mtx(2,2,j)%uu(i,i) + term

          if (diagonal) then
            cycle
          end if

          mtx(2,1,j)%xx      = mtx(2,1,j)%xx      - term
          mtx(2,1,j)%uu(i,i) = mtx(2,1,j)%uu(i,i) - term

          mtx(1,2,j)%xx      = mtx(1,2,j)%xx      - term
          mtx(1,2,j)%uu(i,i) = mtx(1,2,j)%uu(i,i) - term

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
    real(r8) :: c(NEQNS)

    c = eqw * segspr
    do j = 1, ncell
      do i = 1, NEQNS
        term = c(i) / l(i,j)**2
        term_x =   term * n(i,j)%u
        term_u = - term * n(i,j)%x

        r(1,j)%x    = r(1,j)%x    - term_x
        r(1,j)%u(i) = r(1,j)%u(i) - term_u

        r(2,j)%x    = r(2,j)%x    + term_x
        r(2,j)%u(i) = r(2,j)%u(i) + term_u
      end do
    end do

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
    type(NodeVar) :: r_save(NVERT,ncell)

    rh = 1.0_r8 / fdinc

   !!!
   !!!  Compute and save the unperturbed residual.

    call pde_rhs(t)
    call reg_rhs
    call res_mass_matrix

    r_save = r

    do k = 1, NVERT

     !!!
     !!! PARTIALS WITH RESPECT TO X AT VERTEX K

      u(k,:)%x = u(k,:)%x + fdinc

      call preprocessor(errc)
      if (errc /= 0) return
      call pde_rhs(t)
      call reg_rhs
      call res_mass_matrix

      u(k,:)%x = u(k,:)%x - fdinc

      do j = 1, ncell
        do l = 1, NVERT
          mtx(l,k,j)%xx = mtx(l,k,j)%xx + rh * (r(l,j)%x - r_save(l,j)%x)
          mtx(l,k,j)%ux = mtx(l,k,j)%ux + rh * (r(l,j)%u - r_save(l,j)%u)
        end do
      end do

     !!!
     !!! PARTIALS WITH RESPECT TO THE DEPENDENT VARIABLES AT VERTEX K

      do i = 1, NEQNS

        u(k,:)%u(i) = u(k,:)%u(i) + fdinc

        call preprocessor(errc)
        if (errc /= 0) return
        call pde_rhs(t)
        call reg_rhs
        call res_mass_matrix

        u(k,:)%u(i) = u(k,:)%u(i) - fdinc

        do j = 1, ncell
          do l = 1, NVERT
            mtx(l,k,j)%xu(i)   = mtx(l,k,j)%xu(i)   + rh * (r(l,j)%x - r_save(l,j)%x)
            mtx(l,k,j)%uu(:,i) = mtx(l,k,j)%uu(:,i) + rh * (r(l,j)%u - r_save(l,j)%u)
          end do
        end do

      end do

    end do

    errc = 0

  end subroutine eval_dfdy

end module local_mfe

