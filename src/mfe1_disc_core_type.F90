#include "f90_assert.fpp"

module mfe1_disc_core_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe1_vector_type
  use parameter_list_type
  implicit none
  private

  type, public :: mfe1_disc_core
    integer :: neqns, ncell
    real(r8), allocatable :: u(:,:,:), udot(:,:,:)
    real(r8), allocatable :: dx(:), du(:,:)
    real(r8), allocatable :: dudx(:,:)
    real(r8), allocatable :: l(:,:)
    real(r8), allocatable :: n(:,:,:)
    real(r8), allocatable :: r(:,:,:)
    real(r8), allocatable :: mtx(:,:,:,:,:)
    real(r8), allocatable :: eqw(:), eltvsc(:), segspr(:)
    integer :: kreg
    logical :: udot_valid = .false.
  contains
    procedure :: init
    procedure :: update
    procedure :: recompute
    procedure :: assemble_vector
    procedure :: assemble_matrix
    procedure :: assemble_diagonal
    procedure :: reg_rhs
    generic   :: laplacian => lapl_const_coef, lapl_var_coef
    procedure :: res_mass_matrix
    procedure :: eval_mass_matrix
    procedure, private :: lapl_const_coef, lapl_var_coef
  end type

  real(r8), parameter :: ETA = 0.01_r8, &
    C3 = 1.0_r8 / 3.0_r8, C5 = 1.0_r8 / 5.0_r8, C7 = 1.0_r8 / 7.0_r8

contains

  subroutine init(this, neqns, ncell, params, stat, errmsg)

    class(mfe1_disc_core), intent(out), target :: this
    integer, intent(in) :: neqns, ncell
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    this%neqns = neqns
    this%ncell = ncell
    allocate(this%u(neqns+1,2,this%ncell))
    allocate(this%udot, this%r, mold=this%u)
    allocate(this%l(this%neqns,this%ncell))
    allocate(this%du, this%dudx, mold=this%l)
    allocate(this%dx(this%ncell), this%n(2,this%neqns,this%ncell))
    allocate(this%mtx(neqns+1,neqns+1,2,2,this%ncell))
    allocate(this%eqw(this%neqns), this%eltvsc(this%neqns))

    !TODO: rename variables
    call params%get('eltvsc', this%eltvsc, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (size(this%eltvsc) /= this%neqns) then
      stat = 1
      errmsg = 'wrong number of values for eltvsc'
      return
    else if (any(this%eltvsc < 0.0)) then
      stat = 1
      errmsg = 'eltvsc is < 0.0'
      return
    end if

    call params%get('eqw', this%eqw, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (size(this%eqw) /= this%neqns) then
      stat = 1
      errmsg = 'wrong number of values for eqw'
      return
    else if (any(this%eqw <= 0.0)) then
      stat = 1
      errmsg = 'eqw is <= 0.0'
      return
    end if

    call params%get('kreg', this%kreg, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    select case (this%kreg)
    case (1,2) ! okay
    case default
      stat = 1
      errmsg = 'invalid value for kreg'
      return
    end select

    call params%get('segspr', this%segspr, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (size(this%segspr) /= this%neqns) then
      stat = 1
      errmsg = 'wrong number of values for segspr'
      return
    else if (any(this%segspr < 0.0)) then
      stat = 1
      errmsg = 'segspr is < 0.0'
      return
    end if

  end subroutine init

  subroutine update(this, u, udot)
    class(mfe1_disc_core), intent(inout) :: this
    type(mfe1_vector), intent(in) :: u
    type(mfe1_vector), intent(in), optional :: udot
    integer :: j
    do j = 1, this%ncell
      this%u(:,1,j) = u%array(:,j)
      this%u(:,2,j) = u%array(:,j+1)
    end do
    if (present(udot)) then
      do j = 1, this%ncell
        this%udot(:,1,j) = udot%array(:,j)
        this%udot(:,2,j) = udot%array(:,j+1)
      end do
    end if
    this%udot_valid = present(udot)
    call this%recompute
  end subroutine

  subroutine recompute(this)
    class(mfe1_disc_core), intent(inout) :: this
    integer :: i, j
    do j = 1, this%ncell
      this%dx(j) = this%u(this%neqns+1,2,j) - this%u(this%neqns+1,1,j)
      do i = 1, this%neqns
        this%du(i,j) = this%u(i,2,j) - this%u(i,1,j)
        this%l(i,j) = sqrt(this%dx(j)**2 + this%du(i,j)**2)
        this%n(1,i,j) = -this%du(i,j) / this%l(i,j)
        this%n(2,i,j) =  this%dx(j) / this%l(i,j)
        this%dudx(i,j) = this%du(i,j) / this%dx(j)
      end do
    end do
  end subroutine

  subroutine assemble_vector(this, r)
    class(mfe1_disc_core), intent(in) :: this
    type(mfe1_vector), intent(inout) :: r
    integer :: j
    r%array(:,1) = this%r(:,1,1)
    do j = 2, this%ncell
      r%array(:,j) = this%r(:,2,j-1) + this%r(:,1,j)
    end do
    r%array(:,this%ncell+1) = this%r(:,2,this%ncell)
  end subroutine

  subroutine assemble_matrix(this, a)
    use btd_matrix_type
    class(mfe1_disc_core), intent(in) :: this
    type(btd_matrix), intent(inout) :: a
    integer :: j
    a%d(:,:,1) = this%mtx(:,:,1,1,1)
    do j = 2, this%ncell
      a%d(:,:,j) = this%mtx(:,:,1,1,j) + this%mtx(:,:,2,2,j-1)
    end do
    a%d(:,:,this%ncell+1) = this%mtx(:,:,2,2,this%ncell)
    do j = 1, this%ncell
      a%u(:,:,j) = this%mtx(:,:,1,2,j)
    end do
    a%u(:,:,this%ncell+1) = 0.0_r8
    a%l(:,:,1) = 0.0_r8
    do j = 2, this%ncell + 1
      a%l(:,:,j) = this%mtx(:,:,2,1,j-1)
    end do
  end subroutine

  subroutine assemble_diagonal(this, diag)
    class(mfe1_disc_core), intent(in) :: this
    real(r8), intent(out) :: diag(:,:,:)
    integer :: j
    diag(:,:,1) = this%mtx(:,:,1,1,1)
    do j = 2, this%ncell
      diag(:,:,j) = this%mtx(:,:,1,1,j) + this%mtx(:,:,2,2,j-1)
    end do
    diag(:,:,this%ncell+1) = this%mtx(:,:,2,2,this%ncell)
  end subroutine

  subroutine res_mass_matrix(this)  !TODO: rename

    class(mfe1_disc_core), intent(inout) :: this

    integer :: i, j
    real(r8) :: c(this%neqns)
    real(r8) :: term2, dxdot, dudot
    real(r8) :: ndot(2), term(2)

    ASSERT(this%udot_valid)

    !! Pure MFE mass matrix
    c = this%eqw / 6.0_r8
    associate (rx => this%r(this%neqns+1,:,:), xdot => this%udot(this%neqns+1,:,:))
      do j = 1, this%ncell
        do i = 1, this%neqns
          ! Normal velocity at each vertex.
          ndot(:) = this%n(1,i,j) * xdot(:,j) + this%n(2,i,j) * this%udot(i,:,j)
          term(:) = (c(i) * this%l(i,j)) * (sum(ndot) + ndot(:))

          rx(:,j)  = rx(:,j)  - term * this%n(1,i,j)
          this%r(i,:,j) = this%r(i,:,j) - term * this%n(2,i,j)
        end do
      end do
    end associate

    !! Regularization contribution to the mass matrix
    select case (this%kreg)
    case (1)  ! Rate of deformation penalization

      associate (rx => this%r(this%neqns+1,:,:), xdot => this%udot(this%neqns+1,:,:))
        c = this%eqw * this%eltvsc
        do j = 1, this%ncell
          dxdot = xdot(2,j) - xdot(1,j)
          do i = 1, this%neqns
            dudot = this%udot(i,2,j) - this%udot(i,1,j)
            term2 = (c(i) / this%l(i,j)) * (this%n(2,i,j) * dxdot - this%n(1,i,j) * dudot)

            rx(1,j) = rx(1,j) + (term2 * this%n(2,i,j))
            rx(2,j) = rx(2,j) - (term2 * this%n(2,i,j))

            this%r(i,1,j) = this%r(i,1,j) - (term2 * this%n(1,i,j))
            this%r(i,2,j) = this%r(i,2,j) + (term2 * this%n(1,i,j))
          end do
        end do
      end associate

    case (2)  ! Total gradient penalization

      associate (rx => this%r(this%neqns+1,:,:), xdot => this%udot(this%neqns+1,:,:))
        c = this%eqw * this%eltvsc
        do j = 1, this%ncell
          dxdot = xdot(2,j) - xdot(1,j)
          do i = 1, this%neqns
            dudot = this%udot(i,2,j) - this%udot(i,1,j)

            rx(1,j) = rx(1,j) + ((c(i) / this%l(i,j)) * dxdot)
            rx(2,j) = rx(2,j) - ((c(i) / this%l(i,j)) * dxdot)

            this%r(i,1,j) = this%r(i,1,j) + ((c(i) / this%l(i,j)) * dudot)
            this%r(i,2,j) = this%r(i,2,j) - ((c(i) / this%l(i,j)) * dudot)
          end do
        end do
      end associate
    end select

  end subroutine res_mass_matrix

  subroutine eval_mass_matrix(this, factor, diag_only)

    class(mfe1_disc_core), intent(inout) :: this
    real(r8), intent(in), optional :: factor
    logical, intent(in), optional :: diag_only

    integer :: i, j, ix
    logical :: diagonal
    real(r8) :: blk(this%neqns+1,this%neqns+1)
    real(r8) :: fac, term, term_xx, term_xu, term_uu
    real(r8) :: c(this%neqns)

    ix = this%neqns + 1

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

    c = (fac / 3.0_r8) * this%eqw
    do j = 1, this%ncell
      blk = 0.0_r8
      do i = 1, this%neqns
        associate (n1 => this%n(1,i,j), n2 => this%n(2,i,j))
          blk(ix,ix)   =  (c(i) * this%l(i,j)) * n1 * n1 + blk(ix,ix)
          blk(ix,i)    =  (c(i) * this%l(i,j)) * n1 * n2
          blk(i,ix)    =  (c(i) * this%l(i,j)) * n1 * n2
          blk(i,i)     =  (c(i) * this%l(i,j)) * n2 * n2
        end associate
      end do

      ! COPY THE BASIC blk

      this%mtx(:,:,1,1,j) = blk
      this%mtx(:,:,2,2,j) = blk

      if (diagonal) cycle

      blk = 0.5_r8 * blk

      this%mtx(:,:,2,1,j) = blk
      this%mtx(:,:,1,2,j) = blk
    end do

   !!!
   !!!  Regularization contribution to the mass matrix.

    select case (this%kreg)
    case (1)

     !!!
     !!! RATE OF DEFORMATION PENALIZATION

      c = fac * this%eqw * this%eltvsc
      do j = 1, this%ncell
        do i = 1, this%neqns
          associate (n1 => this%n(1,i,j), n2 => this%n(2,i,j))
            term = c(i) / this%l(i,j)
            term_xx =   term * n2 * n2
            term_xu = - term * n1 * n2
            term_uu =   term * n1 * n1
          end associate

          this%mtx(ix,ix,1,1,j) = this%mtx(ix,ix,1,1,j) + term_xx
          this%mtx(ix,i,1,1,j)  = this%mtx(ix,i,1,1,j)  + term_xu
          this%mtx(i,ix,1,1,j)  = this%mtx(i,ix,1,1,j)  + term_xu
          this%mtx(i,i,1,1,j)   = this%mtx(i,i,1,1,j)   + term_uu

          this%mtx(ix,ix,2,2,j) = this%mtx(ix,ix,2,2,j) + term_xx
          this%mtx(ix,i,2,2,j)  = this%mtx(ix,i,2,2,j)  + term_xu
          this%mtx(i,ix,2,2,j)  = this%mtx(i,ix,2,2,j)  + term_xu
          this%mtx(i,i,2,2,j)   = this%mtx(i,i,2,2,j)   + term_uu

          if (diagonal) cycle

          this%mtx(ix,ix,2,1,j) = this%mtx(ix,ix,2,1,j) - term_xx
          this%mtx(ix,i,2,1,j)  = this%mtx(ix,i,2,1,j)  - term_xu
          this%mtx(i,ix,2,1,j)  = this%mtx(i,ix,2,1,j)  - term_xu
          this%mtx(i,i,2,1,j)   = this%mtx(i,i,2,1,j)   - term_uu

          this%mtx(ix,ix,1,2,j) = this%mtx(ix,ix,1,2,j) - term_xx
          this%mtx(ix,i,1,2,j)  = this%mtx(ix,i,1,2,j)  - term_xu
          this%mtx(i,ix,1,2,j)  = this%mtx(i,ix,1,2,j)  - term_xu
          this%mtx(i,i,1,2,j)   = this%mtx(i,i,1,2,j)   - term_uu
        end do
      end do

    case (2)

     !!!
     !!! TOTAL GRADIENT PENALIZATION

      c = fac * this%eqw * this%eltvsc
      do j = 1, this%ncell
        do i = 1, this%neqns
          term = c(i) / this%l(i,j)

          this%mtx(ix,ix,1,1,j) = this%mtx(ix,ix,1,1,j) + term
          this%mtx(i,i,1,1,j)   = this%mtx(i,i,1,1,j)   + term

          this%mtx(ix,ix,2,2,j) = this%mtx(ix,ix,2,2,j) + term
          this%mtx(i,i,2,2,j)   = this%mtx(i,i,2,2,j)   + term

          if (diagonal) cycle

          this%mtx(ix,ix,2,1,j) = this%mtx(ix,ix,2,1,j) - term
          this%mtx(i,i,2,1,j)   = this%mtx(i,i,2,1,j)   - term

          this%mtx(ix,ix,1,2,j) = this%mtx(ix,ix,1,2,j) - term
          this%mtx(i,i,1,2,j)   = this%mtx(i,i,1,2,j)   - term
        end do
      end do
    end select

  end subroutine eval_mass_matrix


  subroutine reg_rhs(this)

    class(mfe1_disc_core), intent(inout) :: this

    integer :: i, j
    real(r8) :: term, term_x, term_u
    real(r8) :: c(this%neqns)

    associate (rx => this%r(this%neqns+1,:,:))
      c = this%eqw * this%segspr
      do j = 1, this%ncell
        do i = 1, this%neqns
          term = c(i) / this%l(i,j)**2
          term_x =   term * this%n(2,i,j)
          term_u = - term * this%n(1,i,j)

          rx(1,j)  = rx(1,j)  - term_x
          this%r(i,1,j) = this%r(i,1,j) - term_u

          rx(2,j)  = rx(2,j)  + term_x
          this%r(i,2,j) = this%r(i,2,j) + term_u
        end do
      end do
    end associate

  end subroutine


  subroutine lapl_const_coef(this, eqno, coef)

    class(mfe1_disc_core), intent(inout) :: this
    integer, intent(in) :: eqno
    real(r8), intent(in) :: coef

    integer :: i
    real(r8) :: c, r1, m, e, s1, s2

    associate (rx => this%r(this%neqns+1,:,:))
      c = this%eqw(eqno) * coef
      do i = 1, this%ncell
        r1 = 1.0_r8 / this%n(2,eqno,i)
        m = this%dudx(eqno,i)

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_r8 + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_r8 + r1)

        rx(1,i) = rx(1,i) - (c * s2)
        this%r(eqno,1,i) = this%r(eqno,1,i) - (c * s1)

        rx(2,i) = rx(2,i) + (c * s2)
        this%r(eqno,2,i) = this%r(eqno,2,i) + (c * s1)
      end do
    end associate

  end subroutine


  subroutine lapl_var_coef(this, eqno, coef)

    class(mfe1_disc_core), intent(inout) :: this
    integer, intent(in) :: eqno
    real(r8), intent(in) :: coef(:,:)

    integer :: i
    real(r8) :: c, r1, m, e, s1, s2

    associate (rx => this%r(this%neqns+1,:,:))
      do i = 1, this%ncell
        r1 = 1.0_r8 / this%n(2,eqno,i)
        m = this%dudx(eqno,i)

        e = m / r1
        if (abs(e) > ETA) then
          s1 = - sign ( log(abs(m) + r1), m )
        else
          s1 = - e * (1.0_r8 + (e**2) * (C3 + (e**2) * (C5 + C7 * (e**2))))
        end if

        s2 = m ** 2 / (1.0_r8 + r1)

        c = this%eqw(eqno) * coef(1,i)
        rx(1,i) = rx(1,i) - c * s2
        this%r(eqno,1,i) = this%r(eqno,1,i) - c * s1

        c = this%eqw(eqno) * coef(2,i)
        rx(2,i) = rx(2,i) + c * s2
        this%r(eqno,2,i) = this%r(eqno,2,i) + c * s1

      end do
    end associate

  end subroutine

end module
