!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_types

  use mfe_constants
  implicit none
  private

  public  :: assignment(=), operator(*), operator(+), operator(-), &
             operator (.Dot.), Norm, operator(/)


  private :: PlusNodeVar, PlusNodeVar1, MinusNodeVar, MinusNodeVar1, &
             ScaleNodeVar, ScaleNodeVar1, ScaleNodeMtx, ScaleNodeMtx1, &
             NodeMtxGetsReal, DotNodeVar1, NormNodeVar1

  type, public :: NodeVar
    real(kind=wp)                 :: x
    real(kind=wp), dimension(NEQNS) :: u
  end type NodeVar

  type, public :: NodeMtx
    real(kind=wp)                     :: xx
    real(kind=wp), dimension(NEQNS)     :: xu, ux
    real(kind=wp), dimension(NEQNS,NEQNS) :: uu
  end type NodeMtx

  interface operator(*)
    module procedure ScaleNodeVar, ScaleNodeVar1, ScaleNodeMtx, ScaleNodeMtx1
  end interface

  interface operator(/)
    module procedure DivNodeVar1
  end interface

  interface operator(+)
    module procedure PlusNodeVar, PlusNodeVar1
  end interface

  interface operator(-)
    module procedure MinusNodeVar, MinusNodeVar1
  end interface

  interface assignment(=)
    module procedure NodeMtxGetsReal
  end interface

  interface operator(.Dot.)
    module procedure DotNodeVar1
  end interface

  interface Norm
    module procedure NormNodeVar1
  end interface

  contains

    function ScaleNodeVar (s, a) result (b)

      real(kind=wp), intent(in) :: s
      type(NodeVar), intent(in) :: a
      type(NodeVar) :: b

      b % x = s * a % x
      b % u = s * a % u

    end function ScaleNodeVar

    function ScaleNodeVar1 (s, a) result (b)

      real(kind=wp), intent(in) :: s
      type(NodeVar), dimension(:), intent(in) :: a
      type(NodeVar), dimension(size(a)) :: b

      integer :: i

      do i = 1, size(a)
        b(i) % x = s * a(i) % x
        b(i) % u = s * a(i) % u
      end do

    end function ScaleNodeVar1

    function DivNodeVar1 (a, s) result (b)

      real(kind=wp), intent(in) :: s
      type(NodeVar), dimension(:), intent(in) :: a
      type(NodeVar), dimension(size(a)) :: b

      integer :: i

      do i = 1, size(a)
        b(i) % x = a(i) % x / s
        b(i) % u = a(i) % u / s
      end do

    end function DivNodeVar1

    function PlusNodeVar (a, b) result (c)

      type(NodeVar), intent(in) :: a, b
      type(NodeVar) :: c

      c % x = a % x + b % x
      c % u = a % u + b % u

    end function PlusNodeVar

    function PlusNodeVar1 (a, b) result (c)

      type(NodeVar), dimension(:), intent(in) :: a, b
      type(NodeVar), dimension(size(a)) :: c

      integer :: i

      do i = 1, size(a)
        c(i) % x = a(i) % x + b(i) % x
        c(i) % u = a(i) % u + b(i) % u
      end do

    end function PlusNodeVar1

    function MinusNodeVar (a, b) result (c)

      type(NodeVar), intent(in) :: a, b
      type(NodeVar) :: c

      c % x = a % x - b % x
      c % u = a % u - b % u

    end function MinusNodeVar

    function MinusNodeVar1 (a, b) result (c)

      type(NodeVar), dimension(:), intent(in) :: a, b
      type(NodeVar), dimension(size(a)) :: c

      integer :: i

      do i = 1, size(a)
        c(i) % x = a(i) % x - b(i) % x
        c(i) % u = a(i) % u - b(i) % u
      end do

    end function MinusNodeVar1

    function DotNodeVar1 (a, b) result (a_dot_b)

      type(NodeVar), dimension(:), intent(in) :: a, b
      real(kind=wp) :: a_dot_b

      integer :: i, k

      a_dot_b = 0.0_wp
      do i = 1, size(a)
        a_dot_b = a_dot_b + a(i) % x * b(i) % x
        do k = 1, NEQNS
          a_dot_b = a_dot_b + a(i) % u(k) * b(i) % u(k)
        end do
      end do

    end function DotNodeVar1

    function NormNodeVar1 (a) result (norm_a)

      type(NodeVar), dimension(:), intent(in) :: a
      real(kind=wp) :: norm_a

      integer :: i, k

      norm_a = 0.0_wp
      do i = 1, size(a)
        norm_a = norm_a + a(i) % x ** 2
        do k = 1, NEQNS
          norm_a = norm_a + a(i) % u(k) ** 2
        end do
      end do
      norm_a = sqrt(norm_a)

    end function NormNodeVar1

    function ScaleNodeMtx (s, a) result (b)

      real(kind=wp), intent(in) :: s
      type(NodeMtx), intent(in) :: a
      type(NodeMtx) :: b

      b % xx = s * a % xx
      b % xu = s * a % xu
      b % ux = s * a % ux
      b % uu = s * a % uu

    end function ScaleNodeMtx

    function ScaleNodeMtx1 (s, a) result (b)

      real(kind=wp), intent(in) :: s
      type(NodeMtx), dimension(:), intent(in) :: a
      type(NodeMtx), dimension(size(a)) :: b

      integer :: i

      do i = 1, size(a)
        b(i) % xx = s * a(i) % xx
        b(i) % xu = s * a(i) % xu
        b(i) % ux = s * a(i) % ux
        b(i) % uu = s * a(i) % uu
      end do

    end function ScaleNodeMtx1

    subroutine NodeMtxGetsReal (a, b)

      type(NodeMtx), intent(out) :: a
      real(kind=wp), intent(in)  :: b

      a % xx = b
      a % xu = b
      a % ux = b
      a % uu = b

    end subroutine NodeMtxGetsReal

end module mfe_types
