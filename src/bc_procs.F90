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
  use mfe1_vector_type
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

    type(mfe1_vector), intent(in) :: u
    type(mfe1_vector), intent(inout) :: r

    integer :: n, neqns

    n = r%nnode
    neqns = r%neqns

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      r%array(neqns+1,1) = u%array(neqns+1,1) - bc_left%x_value
    end if

    where (bc_left%u_type == FIXED)
      r%array(:neqns,1) = u%array(:neqns,1) - bc_left%u_value
    end where

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      r%array(neqns+1,n) = u%array(neqns+1,n) - bc_right%x_value
    end if

    where (bc_right%u_type == FIXED)
      r%array(:neqns,n) = u%array(:neqns,n) - bc_right%u_value
    end where

  end subroutine bc_res


  subroutine bc_diag(diag)

    real(r8), intent(inout) :: diag(:,:,:)

    integer :: n, i, ix

    n = size(diag,dim=3)
    ix = size(diag,dim=1)

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      diag(ix,:,1) = 0.0_r8
      diag(:,ix,1) = 0.0_r8
      diag(ix,ix,1) = 1.0_r8
    end if

    do i = 1, NEQNS
      if (bc_left%u_type(i) == FIXED) then
        diag(i,:,1) = 0.0_r8
        diag(:,i,1) = 0.0_r8
        diag(i,i,1) = 1.0_r8
      end if
    end do

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      diag(ix,:,n) = 0.0_r8
      diag(:,ix,n) = 0.0_r8
      diag(ix,ix,n) = 1.0_r8
    end if

    do i = 1, NEQNS
      if (bc_right%u_type(i) == FIXED) then
        diag(i,:,n) = 0.0_r8
        diag(:,i,n) = 0.0_r8
        diag(i,i,n) = 1.0_r8
      end if
    end do

  end subroutine bc_diag


  subroutine bc_jac(jac_l, jac_d, jac_u)

    real(r8), intent(inout) :: jac_l(:,:,:), jac_d(:,:,:), jac_u(:,:,:)

    integer :: n, i, ix

    n = size(jac_d,dim=3)
    ix = size(jac_d,dim=1)

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      jac_d(ix,:,1) = 0.0_r8
      jac_d(ix,ix,1) = 1.0_r8
      jac_u(ix,:,1) = 0.0_r8
    end if

    do i = 1, NEQNS
      if (bc_left%u_type(i) == FIXED) then
        jac_d(i,:,1) = 0.0_r8
        jac_d(i,i,1) = 1.0_r8
        jac_u(i,:,1) = 0.0_r8
      end if
    end do

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      jac_l(ix,:,n) = 0.0_r8
      jac_d(ix,:,n) = 0.0_r8
      jac_d(ix,ix,n) = 1.0_r8
    end if

    do i = 1, NEQNS
      if (bc_right%u_type(i) == FIXED) then
        jac_l(i,:,n) = 0.0_r8
        jac_d(i,:,n) = 0.0_r8
        jac_d(i,i,n) = 1.0_r8
      end if
    end do

  end subroutine bc_jac

  subroutine bc_udot(a_l, a_d, a_u, g)

    real(r8), intent(inout) :: a_l(:,:,:), a_d(:,:,:), a_u(:,:,:)
    type(mfe1_vector), intent(inout) :: g

    integer :: n, i, ix

    n = g%nnode
    ix = size(a_d,dim=1)

    !!! LEFT ENDPOINT !!!

    if (bc_left%x_type == FIXED) then
      a_d(ix,:,1) = 0.0_r8
      a_d(ix,ix,1) = 1.0_r8
      a_u(ix,:,1) = 0.0_r8
      g%array(neqns+1,1) = 0.0_r8
    end if

    do i = 1, NEQNS
      if (bc_left%u_type(i) == FIXED) then
        a_d(i,:,1) = 0.0_r8
        a_d(i,i,1) = 1.0_r8
        a_u(i,:,1) = 0.0_r8
        g%array(i,1)   = 0.0_r8
      end if
    end do

    !!! RIGHT ENDPOINT !!!

    if (bc_right%x_type == FIXED) then
      a_l(ix,:,n) = 0.0_r8
      a_d(ix,:,n) = 0.0_r8
      a_d(ix,ix,n) = 1.0_r8
      g%array(neqns+1,n) = 0.0_r8
    end if

    do i = 1, NEQNS
      if (bc_right%u_type(i) == FIXED) then
        a_l(i,:,n) = 0.0_r8
        a_d(i,:,n) = 0.0_r8
        a_d(i,i,n) = 1.0_r8
        g%array(i,n) = 0.0_r8
      end if
    end do

  end subroutine bc_udot

end module bc_procs
