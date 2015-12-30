!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bls_solver

  use mfe_constants
  use mfe_types
  implicit none
  private

  public  :: factor, solve

  interface factor
    module procedure fct, vfct, btfct
  end interface

  interface solve
    module procedure slv, mslv, vslv, vmslv, btslv
  end interface

  interface update
    module procedure cmab, ymax
  end interface

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! BTFCT (FACTOR) -- Factor a (block) tridiagonal matrix
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine btfct (l, d, u)

      type(NodeMtx), dimension(:), intent(inout) :: l, d, u

      integer :: i

      call factor (d(1))

      do i = 2, size(d)

        call solve (d(i-1), u(i-1))
        call update (d(i), l(i), u(i-1))
        call factor (d(i))

      end do

    end subroutine btfct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! BTSLV (SOLVE) -- Solve a (block) tridiagonal linear system
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine btslv (l, d, u, b)

      type(NodeMtx), dimension(:), intent(in)    :: l, d, u
      type(NodeVar), dimension(:), intent(inout) :: b

      integer :: i

      call solve (d(1),  b(1))      !!! FORWARD SUBSTITUTION !!!

      do i = 2, size(d)

        call update (b(i), l(i), b(i-1))
        call solve (d(i), b(i))

      end do

      do i = size(d) - 1, 1, -1     !!! BACKWARD SUBSTITUTION !!!

        call update (b(i), u(i), b(i+1))

      end do

    end subroutine btslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! FCT (FACTOR) -- LU-factor a nodal matrix
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fct (a)

      type(NodeMtx), intent(inout) :: a

      integer :: i, j, k
      real(kind=wp) :: lkk, lkj, ujk, lxx, lxj, ujx

      do k = 1, NEQNS
        lkk = a % uu(k,k)
        do j = 1, k - 1
          lkj = a % uu(k,j)
          ujk = a % uu(j,k)
          do i = 1, j - 1
            lkj = lkj - a % uu(k,i) * a % uu(i,j)
            ujk = ujk - a % uu(j,i) * a % uu(i,k)
          end do
          ujk = a % uu(j,j) * ujk
          lkk = lkk - lkj * ujk
          a % uu(k,j) = lkj
          a % uu(j,k) = ujk
        end do
        a % uu(k,k) = 1.0_wp / lkk
      end do

      lxx = a % xx
      do j = 1, NEQNS
        lxj = a % xu(j)
        ujx = a % ux(j)
        do i = 1, j - 1
          lxj = lxj - a % xu(i) * a % uu(i,j)
          ujx = ujx - a % uu(j,i) * a % ux(i)
        end do
        ujx = a % uu(j,j) * ujx
        lxx = lxx - lxj * ujx
        a % xu(j) = lxj
        a % ux(j) = ujx
      end do
      a % xx = 1.0_wp / lxx

    end subroutine fct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! VFCT (FACTOR) -- LU-factor a rank-1 array of nodal matrices
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vfct (a)

      type(NodeMtx), dimension(:), intent(inout) :: a

      integer :: i, j, k, l
      real(kind=wp) :: lkk, lkj, ujk, lxx, lxj, ujx

      do l = 1, size(a)

        do k = 1, NEQNS
          lkk = a(l) % uu(k,k)
          do j = 1, k - 1
            lkj = a(l) % uu(k,j)
            ujk = a(l) % uu(j,k)
            do i = 1, j - 1
              lkj = lkj - a(l) % uu(k,i) * a(l) % uu(i,j)
              ujk = ujk - a(l) % uu(j,i) * a(l) % uu(i,k)
            end do
            ujk = a(l) % uu(j,j) * ujk
            lkk = lkk - lkj * ujk
            a(l) % uu(k,j) = lkj
            a(l) % uu(j,k) = ujk
          end do
          a(l) % uu(k,k) = 1.0_wp / lkk
        end do

        lxx = a(l) % xx
        do j = 1, NEQNS
          lxj = a(l) % xu(j)
          ujx = a(l) % ux(j)
          do i = 1, j - 1
            lxj = lxj - a(l) % xu(i) * a(l) % uu(i,j)
            ujx = ujx - a(l) % uu(j,i) * a(l) % ux(i)
          end do
          ujx = a(l) % uu(j,j) * ujx
          lxx = lxx - lxj * ujx
          a(l) % xu(j) = lxj
          a(l) % ux(j) = ujx
        end do
        a(l) % xx = 1.0_wp / lxx

      end do

    end subroutine vfct

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! SLV (SOLVE) -- Solve a nodal matrix-vector linear system
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine slv (a, b)

      type(NodeMtx), intent(in)    :: a
      type(NodeVar), intent(inout) :: b

      integer  :: i, j
      real(kind=wp) :: bj, bx

      do j = 1, NEQNS
        bj = b % u(j)
        do i = 1, j - 1
          bj = bj - a % uu(j,i) * b % u(i)
        end do
        b % u(j) = bj * a % uu(j,j)
      end do

      bx = b % x
      do i = 1, NEQNS
        bx = bx - a % xu(i) * b % u(i)
      end do
      b % x = bx * a % xx

      do j = NEQNS, 1, -1
        bj = b % u(j)
        do i = j + 1, NEQNS
          bj = bj - a % uu(j,i) * b % u(i)
        end do
        b % u(j) = bj - a % ux(j) * b % x
      end do

    end subroutine slv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! MSLV (SOLVE) -- Solve a nodal matrix-matrix linear system
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mslv (a, b)

      type(NodeMtx), intent(in)    :: a
      type(NodeMtx), intent(inout) :: b

      integer :: i, j, k
      real(kind=wp) :: bj, bx

      do k = 1, NEQNS

        do j = 1, NEQNS
          bj = b % uu(j,k)
          do i = 1, j - 1
            bj = bj - a % uu(j,i) * b % uu(i,k)
          end do
          b % uu(j,k) = bj * a % uu(j,j)
        end do

        bx = b % xu(k)
        do i = 1, NEQNS
          bx = bx - a % xu(i) * b % uu(i,k)
        end do
        b % xu(k) = bx * a % xx

        do j = NEQNS, 1, -1
          bj = b % uu(j,k)
          do i = j + 1, NEQNS
            bj = bj - a % uu(j,i) * b % uu(i,k)
          end do
          b % uu(j,k) = bj - a % ux(j) * b % xu(k)
        end do

      end do

      do j = 1, NEQNS
        bj = b % ux(j)
        do i = 1, j - 1
          bj = bj - a % uu(j,i) * b % ux(i)
        end do
        b % ux(j) = bj * a % uu(j,j)
      end do

      bx = b % xx
      do i = 1, NEQNS
        bx = bx - a % xu(i) * b % ux(i)
      end do
      b % xx = bx * a % xx

      do j = NEQNS, 1, -1
        bj = b % ux(j)
        do i = j + 1, NEQNS
          bj = bj - a % uu(j,i) * b % ux(i)
        end do
        b % ux(j) = bj - a % ux(j) * b % xx
      end do

    end subroutine mslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! VSLV (SOLVE) -- Solve a rank-1 array of nodal matrix-vector linear systems
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vslv (a, b)

      type(NodeMtx), dimension(:), intent(in)    :: a
      type(NodeVar), dimension(:), intent(inout) :: b

      integer :: i, j, l
      real(kind=wp) :: bj, bx

      do l = 1, size(a)

        do j = 1, NEQNS
          bj = b(l) % u(j)
          do i = 1, j - 1
            bj = bj - a(l) % uu(j,i) * b(l) % u(i)
          end do
          b(l) % u(j) = bj * a(l) % uu(j,j)
        end do

        bx = b(l) % x
        do i = 1, NEQNS
          bx = bx - a(l) % xu(i) * b(l) % u(i)
        end do
        b(l) % x = bx * a(l) % xx

        do j = NEQNS, 1, -1
          bj = b(l) % u(j)
          do i = j + 1, NEQNS
            bj = bj - a(l) % uu(j,i) * b(l) % u(i)
          end do
          b(l) % u(j) = bj - a(l) % ux(j) * b(l) % x
        end do

      end do

    end subroutine vslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! VMSLV (SOLVE) -- Solve a rank-1 array of nodal matrix-matrix linear systems
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine vmslv (a, b)

      type(NodeMtx), dimension(:), intent(in)    :: a
      type(NodeMtx), dimension(:), intent(inout) :: b

      integer :: i, j, k, l
      real(kind=wp) :: bj, bx

      do l = 1, size(a)

        do k = 1, NEQNS

          do j = 1, NEQNS
            bj = b(l) % uu(j,k)
            do i = 1, j - 1
              bj = bj - a(l) % uu(j,i) * b(l) % uu(i,k)
            end do
            b(l) % uu(j,k) = bj * a(l) % uu(j,j)
          end do

          bx = b(l) % xu(k)
          do i = 1, NEQNS
            bx = bx - a(l) % xu(i) * b(l) % uu(i,k)
          end do
          b(l) % xu(k) = bx * a(l) % xx

          do j = NEQNS, 1, -1
            bj = b(l) % uu(j,k)
            do i = j + 1, NEQNS
              bj = bj - a(l) % uu(j,i) * b(l) % uu(i,k)
            end do
            b(l) % uu(j,k) = bj - a(l) % ux(j) * b(l) % xu(k)
          end do

        end do

        do j = 1, NEQNS
          bj = b(l) % ux(j)
          do i = 1, j - 1
            bj = bj - a(l) % uu(j,i) * b(l) % ux(i)
          end do
          b(l) % ux(j) = bj * a(l) % uu(j,j)
        end do

        bx = b(l) % xx
        do i = 1, NEQNS
          bx = bx - a(l) % xu(i) * b(l) % ux(i)
        end do
        b(l) % xx = bx * a(l) % xx

        do j = NEQNS, 1, -1
          bj = b(l) % ux(j)
          do i = j + 1, NEQNS
            bj = bj - a(l) % uu(j,i) * b(l) % ux(i)
          end do
          b(l) % ux(j) = bj - a(l) % ux(j) * b(l) % xx
        end do

      end do

    end subroutine vmslv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! YMAX (UPDATE) -- Perform a nodal "vector-minus-matrix-vector" update step
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ymax (c, a, b)

      type(NodeMtx), intent(in)    :: a
      type(NodeVar), intent(in)    :: b
      type(NodeVar), intent(inout) :: c

      integer :: i, j
      real(kind=wp) :: cj, cx

      do j = 1, NEQNS
        cj = c % u(j)
        do i = 1, NEQNS
          cj = cj - a % uu(j,i) * b % u(i)
        end do
        c % u(j) = cj - a % ux(j) * b % x
      end do

      cx = c % x
      do i = 1, NEQNS
        cx = cx - a % xu(i) * b % u(i)
      end do
      c % x = cx - a % xx * b % x

    end subroutine ymax

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! CMAB (UPDATE) -- Perform a nodal "matrix-minus-matrix-matrix" update step
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cmab (c, a, b)

      type(NodeMtx), intent(in)    :: a, b
      type(NodeMtx), intent(inout) :: c

      integer :: i, j, k
      real(kind=wp) :: cjk, cjx, cxj, cxx

      do k = 1, NEQNS
        do j = 1, NEQNS
          cjk = c % uu(j,k)
          do i = 1, NEQNS
            cjk = cjk - a % uu(j,i) * b % uu(i,k)
          end do
          c % uu(j,k) = cjk - a % ux(j) * b % xu(k)
        end do
      end do

      do j = 1, NEQNS
        cjx = c % ux(j)
        do i = 1, NEQNS
          cjx = cjx - a % uu(j,i) * b % ux(i)
        end do
        c % ux(j) = cjx - a % ux(j) * b % xx
      end do

      do j = 1, NEQNS
        cxj = c % xu(j)
        do i = 1, NEQNS
          cxj = cxj - a % xu (i) * b % uu(i,j)
        end do
        c % xu(j) = cxj - a % xx * b % xu(j)
      end do

      cxx = c % xx
      do i = 1, NEQNS
        cxx = cxx - a % xu(i) * b % ux(i)
      end do
      c % xx = cxx - a % xx * b % xx

    end subroutine cmab

end module bls_solver
