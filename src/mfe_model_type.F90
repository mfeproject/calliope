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
    type(index_func), allocatable :: bc_dir(:)
    !! Persistent temporary work space
    real(r8), allocatable, private :: diag(:,:,:)
  contains
    procedure :: init
    procedure :: set_boundary_values
    procedure :: eval_residual
    procedure :: eval_udot
  end type

contains

  subroutine init(this, nnode, params, stat, errmsg)

    use parameter_list_type
    use string_utilities, only: i_to_c

    class(mfe_model), intent(out), target :: this
    integer, intent(in) :: nnode
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist

    this%nnode = nnode
    this%ncell = nnode - 1
    call this%disc%init(this%ncell, params, stat, errmsg)
    if (stat /= 0) return
    this%neqns = this%disc%neqns
    this%nvars = this%neqns + 1

    allocate(this%diag(this%nvars,this%nvars,this%nnode))
    allocate(this%bc_dir(this%nvars))

    if (params%is_sublist('boundary-conditions')) then
      plist => params%sublist('boundary-conditions')
    else
      stat = -1
      errmsg = 'missing "boundary-conditions" sublist parameter'
      return
    end if

    !! Configure the Dirichlet boundary conditions
    select case (this%neqns)
    case (1)  ! Scalar PDE

      block
        logical :: left, right
        call plist%get('dir-left', left, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        call plist%get('dir-right', right, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        this%bc_dir(1)%index = pack([1, this%nnode], [left, right])
      end block

    case (2:) ! PDE system

      block
        integer :: i
        logical, allocatable :: left(:), right(:)
        call plist%get('dir-left', left, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        if (size(left) /= this%neqns) then
          stat = -1
          errmsg = '"dir-left" requires a vector of ' // i_to_c(this%neqns) // ' values'
          return
        end if
        call plist%get('dir-right', right, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        if (size(right) /= this%neqns) then
          stat = -1
          errmsg = '"dir-right" requires a vector of ' // i_to_c(this%neqns) // ' values'
          return
        end if
        do i = 1, this%neqns
          this%bc_dir(i)%index = pack([1, this%nnode], [left(i), right(i)])
        end do
      end block

    end select

    !! Configure fixed/free boundary nodes
    block
      logical :: left, right
      call params%get('fixed-node-left', left, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('fixed-node-right', right, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      this%bc_dir(this%nvars)%index = pack([1, this%nnode], [left, right])
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
    call this%disc%assemble_diagonal(this%diag)

    !! BC modifications to diagonal
    do k = 1, size(this%bc_dir)
      do j = 1, size(this%bc_dir(k)%index)
        n = this%bc_dir(k)%index(j)
        this%diag(:,k,n) = 0.0_r8
        this%diag(k,:,n) = 0.0_r8
        this%diag(k,k,n) = 1.0_r8
      end do
    end do

    !! Multiply by the inverse of the block diagonal
    call vfct(this%diag)
    call vslv(this%diag, r%array)

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
    call a%init(this%nvars, this%nnode)
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
