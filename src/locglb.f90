!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_global

  use mfe_constants
  use mfe_types
  use local_arrays
  implicit none
  private

  public :: gather_local_solution, free_local_arrays, assemble_vector, &
            assemble_matrix, assemble_diagonal

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  GATHER_LOCAL_SOLUTION
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine gather_local_solution (global_u, global_udot)

      type(NodeVar), dimension(:), intent(in)           :: global_u
      type(NodeVar), dimension(:), intent(in), optional :: global_udot

      integer :: j

      ncell = size(global_u) - 1

     !!!
     !!! ALLOCATE WORKING STORAGE

      allocate (u(NVERT,ncell), r(NVERT,ncell), mtx(NVERT,NVERT,ncell), &
                l(NEQNS,ncell), n(NEQNS,ncell), dudx(NEQNS,ncell), dx(ncell), du(NEQNS,ncell))

      if (present(global_udot)) then
        allocate (udot(NVERT,ncell))
      end if

     !!!
     !!! GATHER-UP LOCAL COPIES OF THE SOLUTION VECTORS

      do j = 1, ncell

        u(1,j) % x = global_u(j) % x
        u(1,j) % u = global_u(j) % u

        u(2,j) % x = global_u(j+1) % x
        u(2,j) % u = global_u(j+1) % u

      end do

      if (present(global_udot)) then

        do j = 1, ncell

          udot(1,j) % x = global_udot(j) % x
          udot(1,j) % u = global_udot(j) % u

          udot(2,j) % x = global_udot(j+1) % x
          udot(2,j) % u = global_udot(j+1) % u

        end do

      end if

    end subroutine gather_local_solution

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  FREE_LOCAL_ARRAYS
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine free_local_arrays ()

      if (allocated(u)) then
        deallocate (u, r, mtx, l, n, dudx, dx, du)
      end if

      if (allocated(udot)) then
        deallocate (udot)
      end if

    end subroutine free_local_arrays

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ASSEMBLE_VECTOR
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine assemble_vector (global_r)

      type(NodeVar), dimension(:), intent(out) :: global_r

      integer :: j

      global_r(1) % x = r(1,1) % x
      global_r(1) % u = r(1,1) % u

      do j = 2, ncell
        global_r(j) % x = r(1,j) % x + r(2,j-1) % x
        global_r(j) % u = r(1,j) % u + r(2,j-1) % u
      end do

      global_r(ncell+1) % x = r(2,ncell) % x
      global_r(ncell+1) % u = r(2,ncell) % u

    end subroutine assemble_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ASSEMBLE_MATRIX
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine assemble_matrix (lowr, diag, uppr)

      type(NodeMtx), dimension(:), intent(out) :: lowr, diag, uppr

      integer :: j

      diag(1) = mtx(1,1,1)
      do j = 2, ncell
        diag(j) % xx = mtx(1,1,j) % xx + mtx(2,2,j-1) % xx
        diag(j) % xu = mtx(1,1,j) % xu + mtx(2,2,j-1) % xu
        diag(j) % ux = mtx(1,1,j) % ux + mtx(2,2,j-1) % ux
        diag(j) % uu = mtx(1,1,j) % uu + mtx(2,2,j-1) % uu
      end do
      diag(ncell+1) = mtx(2,2,ncell)

      do j = 1, ncell
        uppr(j) = mtx(1,2,j)
      end do
      uppr(ncell+1) = 0.0_wp

      lowr(1) = 0.0_wp
      do j = 2, ncell + 1
        lowr(j) = mtx(2,1,j-1)
      end do

    end subroutine assemble_matrix

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ASSEMBLE_DIAGONAL
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine assemble_diagonal (diag)

      type(NodeMtx), dimension(:), intent(out) :: diag

      integer :: j

      diag(1) = mtx(1,1,1)
      do j = 2, ncell
        diag(j) % xx = mtx(1,1,j) % xx + mtx(2,2,j-1) % xx
        diag(j) % xu = mtx(1,1,j) % xu + mtx(2,2,j-1) % xu
        diag(j) % ux = mtx(1,1,j) % ux + mtx(2,2,j-1) % ux
        diag(j) % uu = mtx(1,1,j) % uu + mtx(2,2,j-1) % uu
      end do
      diag(ncell+1) = mtx(2,2,ncell)

    end subroutine assemble_diagonal

end module local_global
