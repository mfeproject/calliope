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
  use cell_data_type
  use btd_matrix_type
  use block_linear_solver
  use index_func_type
  implicit none
  private

  type, public :: mfe_model
    integer :: neqns, nvars, ncell, nnode
    type(mfe1_disc) :: disc
    type(index_func), allocatable :: bc_dir(:)
    !! Persistent temporary work space
    real(r8), allocatable, private :: diag(:,:,:)
    type(cell_data(:)), allocatable :: cdata
  contains
    procedure :: init
    procedure :: set_boundary_values
    procedure :: eval_residual
    procedure :: eval_udot
    procedure :: compute_mass_matrix
    procedure :: add_dfdy
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
    integer :: n

    this%nnode = nnode
    this%ncell = nnode - 1
    call this%disc%init(this%ncell, params, stat, errmsg)
    if (stat /= 0) return
    this%neqns = this%disc%neqns
    this%nvars = this%neqns + 1

    n = this%neqns
    allocate(cell_data(n) :: this%cdata)
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

    call compute_res(this, t, u, udot, r)

    !! BC modifications
    do k = 1, size(this%bc_dir)
      do j = 1, size(this%bc_dir(k)%index)
        n = this%bc_dir(k)%index(j)
        r%array(k,n) = u%array(k,n) - this%bc_dir(k)%value(j)
      end do
    end do

    !! Mass matrix block diagonal
    call compute_mass_matrix_diag(this, t, u, this%diag)

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

    errc = 0

    !! Mass matrix
    call a%init(this%nvars, this%nnode)
    call compute_mass_matrix(this, t, u, a)

    !! Right hand side
    call compute_rhs(this, t, u, udot)

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


  subroutine compute_res(this, t, y, ydot, r)

    class(mfe_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(in) :: y, ydot
    type(mfe1_vector), intent(inout) :: r

    integer :: j
    !type(cell_data(this%neqns)) :: cdata
    real(r8) :: f(this%nvars,2)

    !call cdata%init(this%neqns)
    r%array(:,1) = 0.0_r8
    do j = 1, this%ncell
      call this%cdata%update(y%array(:,j:j+1))
      call this%disc%compute_cell_f(t, this%cdata, ydot%array(:,j:j+1), f)
      r%array(:,j) = r%array(:,j) + f(:,1)
      r%array(:,j+1) = f(:,2)
    end do

  end subroutine


  subroutine compute_rhs(this, t, y, g)

    class(mfe_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(in) :: y
    type(mfe1_vector), intent(inout) :: g ! data intent(out)

    integer :: j
    !type(cell_data(this%neqns)) :: cdata
    real(r8) :: gx(2), gu(this%neqns,2)

    associate (x => y%array(this%nvars,:), u => y%array(1:this%neqns,:))
      !call cdata%init(this%neqns)
      g%array(:,1) = 0.0_r8
      do j = 1, this%ncell
        call this%cdata%update(y%array(:,j:j+1))
        call this%disc%compute_cell_rhs(t, this%cdata, gx, gu)
        g%array(this%nvars,j) = g%array(this%nvars,j) + gx(1)
        g%array(1:this%neqns,j) = g%array(1:this%neqns,j) + gu(:,1)
        g%array(this%nvars,j+1) = gx(2)
        g%array(1:this%neqns,j+1) = gu(:,2)
      end do
    end associate

  end subroutine

  subroutine compute_mass_matrix(this, t, y, c)
    class(mfe_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(in) :: y
    type(btd_matrix), intent(inout) :: c  !data is intent(out)

    integer :: j
    !type(cell_data(this%neqns)) :: cdata
    real(r8) :: cell_matrix(this%nvars,this%nvars,2,2)

    associate (x => y%array(this%nvars,:), u => y%array(1:this%neqns,:))
      !call cdata%init(this%neqns)
      c%l(:,:,1) = 0.0_r8 ! unused
      c%d(:,:,1) = 0.0_r8
      do j = 1, this%ncell
        call this%cdata%update(y%array(:,j:j+1))
        call this%disc%compute_cell_mass_matrix(this%cdata, 1.0_r8, cell_matrix)
        c%d(:,:,j)   = cell_matrix(:,:,1,1) + c%d(:,:,j)
        c%l(:,:,j+1) = cell_matrix(:,:,2,1)
        c%u(:,:,j)   = cell_matrix(:,:,1,2)
        c%d(:,:,j+1) = cell_matrix(:,:,2,2)
      end do
      c%u(:,:,this%nnode) = 0.0_r8  ! unused
    end associate

  end subroutine

  subroutine compute_mass_matrix_diag(this, t, y, diag)
    class(mfe_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(in) :: y
    real(r8), intent(out) :: diag(:,:,:)

    integer :: j
    !type(cell_data(this%neqns)) :: cdata
    real(r8) :: cell_diag(this%nvars,this%nvars,2)

    associate (x => y%array(this%nvars,:), u => y%array(1:this%neqns,:))
      !call cdata%init(this%neqns)
      diag(:,:,1) = 0.0_r8
      do j = 1, this%ncell
        call this%cdata%update(y%array(:,j:j+1))
        call this%disc%compute_cell_mass_matrix_diag(this%cdata, 1.0_r8, cell_diag)
        diag(:,:,j)   = cell_diag(:,:,1) + diag(:,:,j)
        diag(:,:,j+1) = cell_diag(:,:,2)
      end do
    end associate

  end subroutine

  subroutine add_dfdy(this, t, y, ydot, dfdy)

    class(mfe_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(mfe1_vector), intent(in) :: y, ydot
    type(btd_matrix), intent(inout) :: dfdy

    integer :: j
    real(r8) :: cell_matrix(this%nvars,this%nvars,2,2)

    do j = 1, this%ncell
      call this%disc%compute_cell_dfdy(t, y%array(:,j:j+1), ydot%array(:,j:j+1), cell_matrix)
      dfdy%d(:,:,j)   = dfdy%d(:,:,j)   + cell_matrix(:,:,1,1)
      dfdy%l(:,:,j+1) = dfdy%l(:,:,j+1) + cell_matrix(:,:,2,1)
      dfdy%u(:,:,j)   = dfdy%u(:,:,j)   + cell_matrix(:,:,1,2)
      dfdy%d(:,:,j+1) = dfdy%d(:,:,j+1) + cell_matrix(:,:,2,2)
    end do

  end subroutine

end module mfe_model_type
