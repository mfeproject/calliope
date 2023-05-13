!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mfe_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants, only: NVARS
  use mfe1_vector_type
  use mfe1_disc_type
  use btd_matrix_type
  use block_linear_solver
  use index_func_type
  use timer_tree_type
  implicit none
  private

  type, public :: mfe_model
    integer :: neqns, nvars, ncell, nnode
    type(mfe1_disc) :: disc
    type(index_func) :: bc_dir(NVARS)
  contains
    procedure :: init
    procedure :: set_boundary_values
    procedure :: eval_residual
    procedure :: eval_udot
  end type

contains

  subroutine init(this, neqns, nnode, params, stat, errmsg)
    use parameter_list_type
    class(mfe_model), intent(out), target :: this
    integer, intent(in) :: neqns, nnode
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%neqns = neqns
    this%nvars = neqns + 1
    this%nnode = nnode
    this%ncell = nnode - 1
    call this%disc%init(this%neqns, this%ncell, params, stat, errmsg)

    block ! Initialize boundary conditions
      integer :: i, n
      logical, allocatable :: bc_left_u_type(:), bc_right_u_type(:)
      logical :: bc_left_x_type, bc_right_x_type

      do i = 1, this%neqns
        n = 0
        call params%get('bc-left-u-type', bc_left_u_type)
        call params%get('bc-right-u-type', bc_right_u_type)
        if (bc_left_u_type(i)) n = n + 1
        if (bc_right_u_type(i)) n = n + 1
        allocate(this%bc_dir(i)%index(n))
        !allocate(this%bc_dir(i)%value(n))  ! allocated by set_boundary_values
        n = 0
        if (bc_left_u_type(i)) then
          n = n + 1
          this%bc_dir(i)%index(n) = 1
        end if
        if (bc_right_u_type(i)) then
          n = n + 1
          this%bc_dir(i)%index(n) = this%nnode
        end if
      end do
      n = 0
      call params%get('bc-left-x-type', bc_left_x_type)
      call params%get('bc-right-x-type', bc_right_x_type)
      if (bc_left_x_type) n = n + 1
      if (bc_right_x_type) n = n + 1
      allocate(this%bc_dir(this%nvars)%index(n))
      !allocate(this%bc_dir(this%nvars)%value(n))  ! allocated by set_boundary_values
      n = 0
      if (bc_left_x_type) then
        n = n + 1
        this%bc_dir(this%nvars)%index(n) = 1
      end if
      if (bc_right_x_type) then
        n = n + 1
        this%bc_dir(this%nvars)%index(n) = this%nnode
      end if
    end block

  end subroutine init

  !! The Dirichlet boundary values are not specified explicitly by the BC input
  !! variables, but are instead copied from the initial conditions.

  subroutine set_boundary_values(this, u)
    class(mfe_model), intent(inout) :: this
    type(mfe1_vector), intent(in) :: u
    integer :: i, j, n
    do i = 1, size(this%bc_dir)
      allocate(this%bc_dir(i)%value(size(this%bc_dir(i)%index)))
      do j = 1, size(this%bc_dir(i)%index)
        n = this%bc_dir(i)%index(j)
        this%bc_dir(i)%value(j) = u%array(i,n)
      end do
    end do
  end subroutine

  subroutine eval_residual(this, u, udot, t, r)

    class(mfe_model), intent(inout) :: this
    type(mfe1_vector), intent(in) :: u, udot
    type(mfe1_vector), intent(inout) :: r
    real(r8), intent(in) :: t

    integer :: j, k, n
    real(r8) :: diag(this%disc%neqns+1,this%disc%neqns+1,u%nnode)

    call this%disc%update(u, udot)
    call this%disc%compute_local_residual(t)
    call this%disc%assemble_vector(r)

    !! BC modifications
    do k = 1, size(this%bc_dir)
      do j = 1, size(this%bc_dir(k)%index)
        n = this%bc_dir(k)%index(j)
        r%array(k,n) = u%array(k,n) - this%bc_dir(k)%value(j)
      end do
    end do

    !! Mass matrix block diagonal
    call this%disc%eval_mass_matrix(diag_only=.true.)
    call this%disc%assemble_diagonal(diag)

    !! BC modifications to diagonal
    do k = 1, size(this%bc_dir)
      do j = 1, size(this%bc_dir(k)%index)
        n = this%bc_dir(k)%index(j)
        diag(:,k,n) = 0.0_r8
        diag(k,:,n) = 0.0_r8
        diag(k,k,n) = 1.0_r8
      end do
    end do

    !! Multiply by the inverse of the block diagonal
    call vfct(diag)
    call vslv(diag, r%array)

  end subroutine eval_residual


  subroutine eval_udot(this, u, t, udot, errc)

    class(mfe_model), intent(inout) :: this
    type(mfe1_vector), intent(in) :: u
    type(mfe1_vector), intent(inout) :: udot
    real(r8), intent(in) :: t
    integer, intent(out) :: errc

    integer :: i, j, n
    type(btd_matrix) :: a ! mass matrix

    call this%disc%update(u)
    errc = 0

    !! Mass matrix
    call this%disc%eval_mass_matrix
    call a%init(NVARS, u%nnode)
    call this%disc%assemble_matrix(a)

    !! Right hand side
    call this%disc%pde_rhs(t)
    call this%disc%reg_rhs
    call this%disc%assemble_vector(udot)

    !! BC modifications
    do i = 1, size(this%bc_dir)
      do j = 1, size(this%bc_dir(i)%index)
        n = this%bc_dir(i)%index(j)
        call a%set_dir_var(i, n)
        udot%array(i,n) = 0.0_r8
      end do
    end do

    !! Solve for du/dt
    call a%factor
    call a%solve(udot%array)

  end subroutine eval_udot

end module mfe_model_type
