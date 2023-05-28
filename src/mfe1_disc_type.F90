#include "f90_assert.fpp"

module mfe1_disc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe1_vector_type
  use pde_class
  use cell_data_type
  implicit none
  private

  type, public :: mfe1_disc
    private
    integer, public :: neqns, nvars, ncell
    real(r8), allocatable :: eqw(:)
    real(r8) :: fdinc
    class(pde(:)), allocatable :: p
    integer :: kreg
    logical :: udot_valid = .false.
    real(r8), allocatable :: eltvsc(:), segspr(:)
  contains
    procedure :: init
    procedure :: compute_cell_f
    procedure :: compute_cell_rhs
    procedure :: compute_cell_mass_matrix
    procedure :: compute_cell_mass_matrix_diag
    procedure :: compute_cell_dfdy
  end type

contains

  subroutine init(this, ncell, params, stat, errmsg)

    use parameter_list_type
    use string_utilities, only: i_to_c

    class(mfe1_disc), intent(out), target :: this
    integer, intent(in) :: ncell
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: rtmp
    character(:), allocatable :: string, libdir
    type(parameter_list), pointer :: plist

    call params%get('pde-library', string, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (params%is_parameter('pde-libdir')) then
      call params%get('pde-libdir', libdir, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call load_pde(libdir // string, this%p, stat, errmsg)
    else
      call load_pde(string, this%p, stat, errmsg)
    end if
    if (stat /= 0) return

    this%neqns = this%p%npde
    this%nvars = this%neqns+1
    this%ncell = ncell

    select case (this%neqns)
    case (1)

      this%eqw = [1.0_r8] ! Not a parameter for scalar PDE

      !TODO: rename variable
      call params%get('eltvsc', rtmp, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (rtmp < 0.0) then
        stat = -1
        errmsg = '"eltvsc" must be >= 0.0'
        return
      end if
      this%eltvsc = [rtmp]

      !TODO: rename variable
      call params%get('segspr', rtmp, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (rtmp < 0.0) then
        stat = -1
        errmsg = '"segspr" must be >= 0.0'
        return
      end if
      this%segspr = [rtmp]

    case (2:)

      !TODO: rename variable
      call params%get('eqw', this%eqw, default=spread(1.0_r8,dim=1,ncopies=this%neqns), stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (size(this%eqw) /= this%neqns) then
        stat = -1
        errmsg = '"eqw" requires a vector of ' // i_to_c(this%neqns) // ' values'
        return
      else if (any(this%eqw <= 0.0)) then
        stat = -1
        errmsg = '"eqw" values must be > 0.0'
        return
      end if

      !TODO: rename variable
      call params%get('eltvsc', this%eltvsc, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (size(this%eltvsc) /= this%neqns) then
        stat = -1
        errmsg = '"eltvsc" requires a vector of ' // i_to_c(this%neqns) // ' values'
        return
      else if (any(this%eltvsc < 0.0)) then
        stat = -1
        errmsg = '"eltvsc" values must be >= 0.0'
        return
      end if

      !TODO: rename variable
      call params%get('segspr', this%segspr, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (size(this%segspr) /= this%neqns) then
        stat = -1
        errmsg = '"segspr" requires a vector of ' // i_to_c(this%neqns) // ' values'
        return
      else if (any(this%segspr < 0.0)) then
        stat = -1
        errmsg = '"segspr" values must be >= 0.0'
        return
      end if
    end select

    !TODO: rename variable
    call params%get('kreg', this%kreg, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    select case (this%kreg)
    case (1,2) ! okay
    case default
      stat = -1
      errmsg = 'invalid value for "kreg"'
      return
    end select

    !TODO: rename variable
    call params%get('fdinc', this%fdinc, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%fdinc <= 0.0) then
      stat = -1
      errmsg = '"fdinc" must be > 0.0'
      return
    end if

    if (params%is_sublist('pde-params')) then
      plist => params%sublist('pde-params')
      call this%p%init(this%eqw, plist, stat, errmsg)
      if (stat /= 0) return
    else
      stat = -1
      errmsg = 'missing "pde-params" sublist parameter'
      return
    end if

  end subroutine init

  subroutine load_pde(libname, p, stat, errmsg)

    use fortran_dynamic_loader
    use,intrinsic :: iso_c_binding

    character(*), intent(in) :: libname
    class(pde(:)), allocatable, intent(out) :: p
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(shlib) :: pde_lib
    type(pde_box), pointer :: box
    type(c_ptr) :: cp
    type(c_funptr) :: cfp

    abstract interface
      subroutine f(cp) bind(c)
        import c_ptr
        type(c_ptr), intent(out) :: cp
      end subroutine
    end interface
    procedure(f), pointer :: alloc_pde

    call pde_lib%open(libname, RTLD_LAZY, stat, errmsg)
    if (stat /= 0) return
    call pde_lib%func('alloc_pde', cfp, stat, errmsg)
    if (stat /= 0) return
    call c_f_procpointer(cfp, alloc_pde)
    call alloc_pde(cp)
    call c_f_pointer(cp, box)
    call move_alloc(box%p, p)

  end subroutine load_pde

  pure subroutine compute_cell_f(this, t, cdata, ydot, f)
    class(mfe1_disc), intent(inout) :: this
    real(r8), intent(in) :: t
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: ydot(:,:)
    real(r8), intent(out) :: f(:,:)
    associate (fx => f(this%nvars,:), fu => f(1:this%neqns,:), &
        xdot => ydot(this%nvars,:), udot => ydot(1:this%neqns,:))
      call this%compute_cell_rhs(t, cdata, fx, fu)
      call subtract_cell_lhs(this, cdata, xdot, udot, fx, fu)
    end associate
  end subroutine

  !! Computes the local RHS of the MFE equations
  pure subroutine compute_cell_rhs(this, t, cdata, gx, gu)
    class(mfe1_disc), intent(inout) :: this
    real(r8), intent(in) :: t
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(out) :: gx(:), gu(:,:)
    call this%p%rhs(t, cdata, gx, gu)
    call add_cell_press(this, cdata, gx, gu)
  end subroutine

  !! Local cell pressure regularization contribution to the MFE RHS
  pure subroutine add_cell_press(this, cdata, gx, gu)
    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(inout) :: gx(:), gu(:,:)
    integer  :: i
    real(r8) :: c, term_x, term_u
    do i = 1, this%neqns
      c = this%eqw(i) * this%segspr(i) / cdata%l(i)**2
      term_x =   c * cdata%nu(i)
      term_u = - c * cdata%nx(i)
      gx(1)   = gx(1)   - term_x
      gu(i,1) = gu(i,1) - term_u
      gx(2)   = gx(2)   + term_x
      gu(i,2) = gu(i,2) + term_u
    end do
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Mass matrix contribution to the local F
  pure subroutine subtract_cell_lhs(this, cdata, xdot, udot, rx, ru)
    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: xdot(:), udot(:,:)
    real(r8), intent(inout) :: rx(:), ru(:,:)
    call subtract_cell_mm(this, cdata, xdot, udot, rx, ru)
    select case (this%kreg)
    case (1)
      call subtract_cell_rod_reg(this, cdata, xdot, udot, rx, ru)
    case (2)
      call subtract_cell_tg_reg(this, cdata, xdot, udot, rx, ru)
    end select
  end subroutine

  !! Pure MFE mass matrix contribution to the local F
  pure subroutine subtract_cell_mm(this, cdata, xdot, udot, rx, ru)
    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: xdot(:), udot(:,:)
    real(r8), intent(inout) :: rx(:), ru(:,:)
    integer  :: i
    real(r8) :: c, ndot(2), term(2)
    do i = 1, this%neqns
      c = this%eqw(i)/6.0_r8
      ndot = cdata%nx(i) * xdot(:) + cdata%nu(i) * udot(i,:)
      term = (c * cdata%l(i)) * (sum(ndot) + ndot(:))

      rx(:)   = rx(:)   - term * cdata%nx(i)
      ru(i,:) = ru(i,:) - term * cdata%nu(i)
    end do
  end subroutine

  !! Rate-of-deformation dynamic regularization contribution to the local F
  pure subroutine subtract_cell_rod_reg(this, cdata, xdot, udot, rx, ru)
    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: xdot(:), udot(:,:)
    real(r8), intent(inout) :: rx(:), ru(:,:)
    integer  :: i
    real(r8) :: c, term
    do i = 1, this%neqns
      c = this%eqw(i)*this%eltvsc(i)/cdata%l(i)
      term = c*(cdata%nu(i)*(xdot(2)-xdot(1)) - cdata%nx(i)*(udot(i,2)-udot(i,1)))
      rx(1) = rx(1) + (term * cdata%nu(i))
      rx(2) = rx(2) - (term * cdata%nu(i))
      ru(i,1) = ru(i,1) - (term * cdata%nx(i))
      ru(i,2) = ru(i,2) + (term * cdata%nx(i))
    end do
  end subroutine

  !! Total gradient dynamic regularization contribution to the local F
  pure subroutine subtract_cell_tg_reg(this, cdata, xdot, udot, rx, ru)
    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: xdot(:), udot(:,:)
    real(r8), intent(inout) :: rx(:), ru(:,:)
    integer  :: i
    real(r8) :: c, term
    do i = 1, this%neqns
      c = this%eqw(i)*this%eltvsc(i)/cdata%l(i)
      term = c*(xdot(2) - xdot(1))
      rx(1) = rx(1) + term
      rx(2) = rx(2) - term
      c = this%eqw(i)*this%eltvsc(i)/cdata%l(i)
      term = c*(udot(i,2) - udot(i,1))
      ru(i,1) = ru(i,1) + term
      ru(i,2) = ru(i,2) - term
    end do
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Compute the local mass matrix
  subroutine compute_cell_mass_matrix(this, cdata, fct, mtx)
    class(mfe1_disc), intent(inout) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(out) :: mtx(:,:,:,:)
    call compute_cell_mfe_mass_matrix(this, cdata, fct, mtx)
    select case (this%kreg)
    case (1)
      call add_cell_rod_matrix(this, cdata, fct, mtx)
    case (2)
      call add_cell_tg_matrix(this, cdata, fct, mtx)
    end select
  end subroutine

  !! Compute the local pure MFE mass matrix
  subroutine compute_cell_mfe_mass_matrix(this, cdata, fct, mtx)
    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(out) :: mtx(:,:,:,:)
    integer :: i
    real(r8) :: c, blk(this%nvars,this%nvars)
    blk = 0.0_r8
    do i = 1, this%neqns
      c = fct*this%eqw(i)*cdata%l(i)/3.0_r8
      associate (n1 => cdata%nx(i), n2 => cdata%nu(i), ix => this%nvars)
        blk(ix,ix) = c*n1*n1 + blk(ix,ix)
        blk(ix,i)  = c*n1*n2
        blk(i,ix)  = c*n1*n2
        blk(i,i)   = c*n2*n2
      end associate
    end do
    mtx(:,:,1,1) = blk
    mtx(:,:,2,2) = blk
    blk = 0.5_r8 * blk
    mtx(:,:,2,1) = blk
    mtx(:,:,1,2) = blk
  end subroutine

  !! Rate-of-deformation dynamic regularization contribution to the local mass matrix
  subroutine add_cell_rod_matrix(this, cdata, fct, mtx)

    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(inout) :: mtx(:,:,:,:)

    integer :: i
    real(r8) :: c, term_xx, term_xu, term_uu

    do i = 1, this%neqns
      c = fct * this%eqw(i) * this%eltvsc(i) / cdata%l(i)
      associate (n1 => cdata%nu(i), n2 => cdata%nu(i), ix => this%nvars)
        term_xx =   c * n2 * n2
        term_xu = - c * n1 * n2
        term_uu =   c * n1 * n1

        mtx(ix,ix,1,1) = mtx(ix,ix,1,1) + term_xx
        mtx(ix,i,1,1)  = mtx(ix,i,1,1)  + term_xu
        mtx(i,ix,1,1)  = mtx(i,ix,1,1)  + term_xu
        mtx(i,i,1,1)   = mtx(i,i,1,1)   + term_uu

        mtx(ix,ix,2,2) = mtx(ix,ix,2,2) + term_xx
        mtx(ix,i,2,2)  = mtx(ix,i,2,2)  + term_xu
        mtx(i,ix,2,2)  = mtx(i,ix,2,2)  + term_xu
        mtx(i,i,2,2)   = mtx(i,i,2,2)   + term_uu

        mtx(ix,ix,2,1) = mtx(ix,ix,2,1) - term_xx
        mtx(ix,i,2,1)  = mtx(ix,i,2,1)  - term_xu
        mtx(i,ix,2,1)  = mtx(i,ix,2,1)  - term_xu
        mtx(i,i,2,1)   = mtx(i,i,2,1)   - term_uu

        mtx(ix,ix,1,2) = mtx(ix,ix,1,2) - term_xx
        mtx(ix,i,1,2)  = mtx(ix,i,1,2)  - term_xu
        mtx(i,ix,1,2)  = mtx(i,ix,1,2)  - term_xu
        mtx(i,i,1,2)   = mtx(i,i,1,2)   - term_uu
      end associate
    end do

  end subroutine


  !! Total gradient dynamic regularization contribution to the local mass matrix
  subroutine add_cell_tg_matrix(this, cdata, fct, mtx)

    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(inout) :: mtx(:,:,:,:)

    integer :: i
    real(r8) :: c

    do i = 1, this%neqns
      c = fct * this%eqw(i) * this%eltvsc(i) / cdata%l(i)
      associate (ix => this%nvars)
        mtx(ix,ix,1,1) = mtx(ix,ix,1,1) + c
        mtx(i,i,1,1)   = mtx(i,i,1,1)   + c

        mtx(ix,ix,2,2) = mtx(ix,ix,2,2) + c
        mtx(i,i,2,2)   = mtx(i,i,2,2)   + c

        mtx(ix,ix,2,1) = mtx(ix,ix,2,1) - c
        mtx(i,i,2,1)   = mtx(i,i,2,1)   - c

        mtx(ix,ix,1,2) = mtx(ix,ix,1,2) - c
        mtx(i,i,1,2)   = mtx(i,i,1,2)   - c
      end associate
    end do

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Compute the block diagonal of the local mass matrix
  pure subroutine compute_cell_mass_matrix_diag(this, cdata, fct, diag)
    class(mfe1_disc), intent(inout) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(out) :: diag(:,:,:)
    call compute_cell_mfe_mass_matrix_diag(this, cdata, fct, diag)
    select case (this%kreg)
    case (1)
      call add_cell_rod_matrix_diag(this, cdata, fct, diag)
    case (2)
      call add_cell_tg_matrix_diag(this, cdata, fct, diag)
    end select
  end subroutine

  !! Compute the block diagonal of the local pure MFE mass matrix
  pure subroutine compute_cell_mfe_mass_matrix_diag(this, cdata, fct, diag)

    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(out) :: diag(:,:,:)

    integer :: i
    real(r8) :: c, blk(this%nvars,this%nvars)

    blk = 0.0_r8
    do i = 1, this%neqns
      c = fct*this%eqw(i)*cdata%l(i)/3.0_r8
      associate (n1 => cdata%nx(i), n2 => cdata%nu(i), ix => this%nvars)
        blk(ix,ix) = c*n1*n1 + blk(ix,ix)
        blk(ix,i)  = c*n1*n2
        blk(i,ix)  = c*n1*n2
        blk(i,i)   = c*n2*n2
      end associate
    end do

    diag(:,:,1) = blk
    diag(:,:,2) = blk

  end subroutine

  !! Rate-of-deformation dynamic regularization contribution to the local mass matrix diagonal
  pure subroutine add_cell_rod_matrix_diag(this, cdata, fct, diag)

    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(inout) :: diag(:,:,:)

    integer :: i
    real(r8) :: c, term_xx, term_xu, term_uu

    do i = 1, this%neqns
      c = fct * this%eqw(i) * this%eltvsc(i) / cdata%l(i)
      associate (n1 => cdata%nu(i), n2 => cdata%nu(i), ix => this%nvars)
        term_xx =   c * n2 * n2
        term_xu = - c * n1 * n2
        term_uu =   c * n1 * n1

        diag(ix,ix,1) = diag(ix,ix,1) + term_xx
        diag(ix,i,1)  = diag(ix,i,1)  + term_xu
        diag(i,ix,1)  = diag(i,ix,1)  + term_xu
        diag(i,i,1)   = diag(i,i,1)   + term_uu

        diag(ix,ix,2) = diag(ix,ix,2) + term_xx
        diag(ix,i,2)  = diag(ix,i,2)  + term_xu
        diag(i,ix,2)  = diag(i,ix,2)  + term_xu
        diag(i,i,2)   = diag(i,i,2)   + term_uu
      end associate
    end do

  end subroutine

  !! Rate-of-deformation dynamic regularization contribution to the local mass matrix diagonal
  pure subroutine add_cell_tg_matrix_diag(this, cdata, fct, diag)

    class(mfe1_disc), intent(in) :: this
    type(cell_data(*)), intent(in) :: cdata
    real(r8), intent(in) :: fct
    real(r8), intent(inout) :: diag(:,:,:)

    integer :: i
    real(r8) :: c

    do i = 1, this%neqns
      c = fct * this%eqw(i) * this%eltvsc(i) / cdata%l(i)
      associate (ix => this%nvars)
        diag(ix,ix,1) = diag(ix,ix,1) + c
        diag(i,i,1)   = diag(i,i,1)   + c

        diag(ix,ix,2) = diag(ix,ix,2) + c
        diag(i,i,2)   = diag(i,i,2)   + c
      end associate
    end do

  end subroutine


  subroutine compute_cell_dfdy(this, t, y, ydot, dfdy)

    class(mfe1_disc), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: y(:,:), ydot(:,:)
    real(r8), intent(out) :: dfdy(:,:,:,:)

    real(r8) :: f0(this%nvars,2), f(this%nvars,2)

    integer :: i, k, l
    real(r8) :: rh
    type(cell_data(this%neqns)) :: cdata  !TODO? make persistent workspace? pass as argument?

    rh = 1.0_r8 / this%fdinc

    !call cdata%init(this%neqns)
    call cdata%update(y)
    call compute_cell_f(this, t, cdata, ydot, f0)

    do k = 1, 2
      do i = 1, this%nvars
        call cdata%set_val(i, k, y(i,k) + this%fdinc, update=.true.) !TODO: inverted cell check?
        call compute_cell_f(this, t, cdata, ydot, f)
        call cdata%set_val(i, k, y(i,k), update=.false.)
          do l = 1, 2
            dfdy(:,i,l,k) = rh * (f(:,l) - f0(:,l))
          end do
      end do
    end do

  end subroutine

end module
