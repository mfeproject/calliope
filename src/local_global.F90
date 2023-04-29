!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module local_global

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants, only: NVERT, NVARS
  use mfe1_vector_type
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

  subroutine gather_local_solution(global_u, global_udot)

    type(mfe1_vector), intent(in) :: global_u
    type(mfe1_vector), intent(in), optional :: global_udot

    integer :: j, neqns

    ncell = global_u%nnode - 1
    neqns = global_u%neqns
    ix = neqns + 1

   !!!
   !!! ALLOCATE WORKING STORAGE

    allocate(u(neqns+1,NVERT,ncell), r(neqns+1,NVERT,ncell), &
             l(neqns,ncell), n(2,neqns,ncell), dudx(neqns,ncell), dx(ncell), du(neqns,ncell))
    allocate(mtx(NVARS,NVARS,NVERT,NVERT,ncell))

    if (present(global_udot)) allocate(udot(neqns+1,NVERT,ncell))

   !!!
   !!! GATHER-UP LOCAL COPIES OF THE SOLUTION VECTORS

    do j = 1, ncell
      u(:,1,j) = global_u%array(:,j)
      u(:,2,j) = global_u%array(:,j+1)
    end do

    if (present(global_udot)) then
      do j = 1, ncell
        udot(:,1,j) = global_udot%array(:,j)
        udot(:,2,j) = global_udot%array(:,j+1)
      end do
    end if

  end subroutine gather_local_solution

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  FREE_LOCAL_ARRAYS
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_local_arrays
    if (allocated(u)) deallocate(u, r, mtx, l, n, dudx, dx, du)
    if (allocated(udot)) deallocate (udot)
  end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  ASSEMBLE_VECTOR
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine assemble_vector(global_r)

    type(mfe1_vector), intent(inout) :: global_r

    integer :: j

    global_r%array(:,1) = r(:,1,1)
    do j = 2, ncell
      global_r%array(:,j) = r(:,2,j-1) + r(:,1,j)
    end do
    global_r%array(:,ncell+1) = r(:,2,ncell)

  end subroutine assemble_vector

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  ASSEMBLE_MATRIX
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine assemble_matrix(lowr, diag, uppr)

    real(r8), intent(out) :: lowr(:,:,:), diag(:,:,:), uppr(:,:,:)

    integer :: j

    diag(:,:,1) = mtx(:,:,1,1,1)
    do j = 2, ncell
      diag(:,:,j) = mtx(:,:,1,1,j) + mtx(:,:,2,2,j-1)
    end do
    diag(:,:,ncell+1) = mtx(:,:,2,2,ncell)

    do j = 1, ncell
      uppr(:,:,j) = mtx(:,:,1,2,j)
    end do
    uppr(:,:,ncell+1) = 0.0_r8

    lowr(:,:,1) = 0.0_r8
    do j = 2, ncell + 1
      lowr(:,:,j) = mtx(:,:,2,1,j-1)
    end do

  end subroutine assemble_matrix

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  ASSEMBLE_DIAGONAL
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine assemble_diagonal(diag)

    real(r8), intent(out) :: diag(:,:,:)

    integer :: j

    diag(:,:,1) = mtx(:,:,1,1,1)
    do j = 2, ncell
      diag(:,:,j) = mtx(:,:,1,1,j) + mtx(:,:,2,2,j-1)
    end do
    diag(:,:,ncell+1) = mtx(:,:,2,2,ncell)

  end subroutine assemble_diagonal

end module local_global
