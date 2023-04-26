!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bc_procs

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants
  use bc_data
  use mfe_types, only: NodeVar, NodeMtx
  implicit none
  private

  public :: force_bc

  interface force_bc
    module procedure bc_res, bc_diag, bc_jac, bc_udot
  end interface

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RESIDUAL_BC
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bc_res(u, r)

    type(NodeVar), intent(in) :: u(:)
    type(NodeVar), intent(inout) :: r(:)

    integer :: n

    n = size(r)

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      r(1)%x = u(1)%x - bc_left%x_value
    end if

    where (bc_left%u_type == FIXED)
      r(1)%u = u(1)%u - bc_left%u_value
    end where

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      r(n)%x = u(n)%x - bc_right%x_value
    end if

    where (bc_right%u_type == FIXED)
      r(n)%u = u(n)%u - bc_right%u_value
    end where

  end subroutine bc_res


  subroutine bc_diag(diag)

    type(NodeMtx), intent(inout) :: diag(:)

    integer :: n, i

    n = size(diag)

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      diag(1)%xu = 0.0_r8
      diag(1)%ux = 0.0_r8
      diag(1)%xx = 1.0_r8
    end if

    do i = 1, NEQNS
      if (bc_left%u_type(i) == FIXED) then
        diag(1)%uu(i,:) = 0.0_r8
        diag(1)%uu(:,i) = 0.0_r8
        diag(1)%uu(i,i) = 1.0_r8
        diag(1)%ux(i)   = 0.0_r8
        diag(1)%xu(i)   = 0.0_r8
      end if
    end do

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      diag(n)%xu = 0.0_r8
      diag(n)%ux = 0.0_r8
      diag(n)%xx = 1.0_r8
    end if

    do i = 1, NEQNS
      if (bc_right%u_type(i) == FIXED) then
        diag(n)%uu(i,:) = 0.0_r8
        diag(n)%uu(:,i) = 0.0_r8
        diag(n)%uu(i,i) = 1.0_r8
        diag(n)%ux(i)   = 0.0_r8
        diag(n)%xu(i)   = 0.0_r8
      end if
    end do

  end subroutine bc_diag


  subroutine bc_jac(jac_l, jac_d, jac_u)

    type(NodeMtx), intent(inout) :: jac_l(:), jac_d(:), jac_u(:)

    integer :: n, i

    n = size(jac_d)

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      jac_d(1)%xu = 0.0_r8
      jac_d(1)%xx = 1.0_r8
      jac_u(1)%xu = 0.0_r8
      jac_u(1)%xx = 0.0_r8
    end if

    do i = 1, NEQNS
      if (bc_left%u_type(i) == FIXED) then
        jac_d(1)%uu(i,:) = 0.0_r8
        jac_d(1)%uu(i,i) = 1.0_r8
        jac_d(1)%ux(i)   = 0.0_r8
        jac_u(1)%uu(i,:) = 0.0_r8
        jac_u(1)%ux(i)   = 0.0_r8
      end if
    end do

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      jac_l(n)%xu = 0.0_r8
      jac_l(n)%xx = 0.0_r8
      jac_d(n)%xu = 0.0_r8
      jac_d(n)%xx = 1.0_r8
    end if

    do i = 1, NEQNS
      if (bc_right%u_type(i) == FIXED) then
        jac_l(n)%uu(i,:) = 0.0_r8
        jac_l(n)%ux(i)   = 0.0_r8
        jac_d(n)%uu(i,:) = 0.0_r8
        jac_d(n)%uu(i,i) = 1.0_r8
        jac_d(n)%ux(i)   = 0.0_r8
      end if
    end do

  end subroutine bc_jac

  subroutine bc_udot(a_l, a_d, a_u, g)

    type(NodeMtx), intent(inout) :: a_l(:), a_d(:), a_u(:)
    type(NodeVar), intent(inout) :: g(:)

    integer :: n, i

    n = size(g)

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      a_d(1)%xu = 0.0_r8
      a_d(1)%xx = 1.0_r8
      a_u(1)%xu = 0.0_r8
      a_u(1)%xx = 0.0_r8
      g(1)%x = 0.0_r8
    end if

    do i = 1, NEQNS
      if (bc_left%u_type(i) == FIXED) then
        a_d(1)%uu(i,:) = 0.0_r8
        a_d(1)%uu(i,i) = 1.0_r8
        a_d(1)%ux(i)   = 0.0_r8
        a_u(1)%uu(i,:) = 0.0_r8
        a_u(1)%ux(i)   = 0.0_r8
        g(1)%u(i) = 0.0_r8
      end if
    end do

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      a_l(n)%xu = 0.0_r8
      a_l(n)%xx = 0.0_r8
      a_d(n)%xu = 0.0_r8
      a_d(n)%xx = 1.0_r8
      g(n)%x = 0.0_r8
    end if

    do i = 1, NEQNS
      if (bc_right%u_type(i) == FIXED) then
        a_l(n)%uu(i,:) = 0.0_r8
        a_l(n)%ux(i)   = 0.0_r8
        a_d(n)%uu(i,:) = 0.0_r8
        a_d(n)%uu(i,i) = 1.0_r8
        a_d(n)%ux(i)   = 0.0_r8
        g(n)%u(i) = 0.0_r8
      end if
    end do

  end subroutine bc_udot

end module bc_procs
