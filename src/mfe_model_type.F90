#include "f90_assert.fpp"

module mfe_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use idaesol_type, only: idaesol_model
  use mfe_constants, only: NEQNS, NVARS
  implicit none
  private
  
  type, extends(idaesol_model), public :: mfe_model
    integer :: nnode
  contains
    procedure :: init
    procedure :: size => model_size
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: corr_norm
    procedure :: check_state
  end type mfe_model

contains

  subroutine init (this, nnode)
    class(mfe_model), intent(out) :: this
    integer, intent(in) :: nnode
    this%nnode = nnode
  end subroutine init

  integer function model_size (this)
    class(mfe_model), intent(in) :: this
    model_size = NVARS*this%nnode
  end function model_size
  
  subroutine compute_f (this, t, u, udot, f)
    use mfe_types, only: NodeVar
    use mfe_procs, only: eval_residual
    class(mfe_model) :: this
    real(r8), intent(in)  :: t, u(:), udot(:)
    real(r8), intent(out) :: f(:)
    type(NodeVar(NEQNS)) :: ustruct(this%nnode), udotstruct(this%nnode), fstruct(this%nnode)
    call copy_to_nodevar (u, ustruct)
    call copy_to_nodevar (udot, udotstruct)
    call eval_residual (ustruct, udotstruct, t, fstruct)
    call copy_from_nodevar (fstruct, f)
  end subroutine compute_f
  
  subroutine apply_precon (this, t, u, f)
    use mfe_types, only: NodeVar
    use mfe_procs, only: apply_inverse_jacobian
    class(mfe_model) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), intent(inout) :: f(:)
    type(NodeVar(NEQNS)) :: fstruct(this%nnode)
    call copy_to_nodevar (f, fstruct)
    call apply_inverse_jacobian (fstruct)
    call copy_from_nodevar (fstruct, f)
  end subroutine apply_precon
  
  subroutine compute_precon (this, t, u, udot, dt)
  
    use mfe_types, only: NodeVar
    use mfe_procs, only: eval_jacobian
  
    class(mfe_model) :: this
    real(r8), intent(in)  :: t, u(:), udot(:), dt
    
    integer :: j, stat
    type(NodeVar(NEQNS)) :: ustruct(this%nnode), udotstruct(this%nnode)
    
    call copy_to_nodevar (u, ustruct)
    call copy_to_nodevar (udot, udotstruct)
    call eval_jacobian (ustruct, udotstruct, t, dt, stat)
    INSIST(stat == 0)

  end subroutine compute_precon
  
  
  subroutine corr_norm (this, u, du, error)
  
    use norm_procs, only: eval_norm
    use mfe_types, only: NodeVar
    
    class(mfe_model) :: this
    real(r8), intent(in) :: u(:), du(:)
    real(r8), intent(out) :: error
    
    type(NodeVar(NEQNS)) :: dustruct(this%nnode)
    
    call copy_to_nodevar (du, dustruct)
    error = eval_norm(dustruct, 0)

  end subroutine corr_norm
  

  subroutine check_state (this, u, stage, stat)
  
    use norm_procs, only: check_soln
    use mfe_types, only: NodeVar
    
    class(mfe_model) :: this
    real(r8), intent(in)  :: u(:)
    integer,  intent(in)  :: stage
    integer,  intent(out) :: stat

    type(NodeVar(NEQNS)) :: ustruct(this%nnode)
    
    call copy_to_nodevar (u, ustruct)
    call check_soln (ustruct, stage, stat)
    
  end subroutine check_state
  
  
  subroutine copy_to_nodevar (u, ustruct)
    use mfe_types, only: NodeVar
    real(r8), intent(in), target :: u(:)
    type(NodeVar(*)), intent(out) :: ustruct(:)
    integer :: j, k
    real(r8), pointer :: u2(:,:)
    u2(1:NVARS,1:size(ustruct)) => u
    do j = 1, size(ustruct)
      do k = 1, NEQNS
        ustruct(j)%u(k) = u2(k,j)
      end do
      ustruct(j)%x = u2(k,j)
    end do
  end subroutine
  
  
  subroutine copy_from_nodevar (ustruct, u)
    use mfe_types, only: NodeVar
    type(NodeVar(*)), intent(in) :: ustruct(:)
    real(r8), intent(out), target :: u(:)
    integer :: j, k
    real(r8), pointer :: u2(:,:)
    u2(1:NVARS,1:size(ustruct)) => u
    do j = 1, size(ustruct)
      do k = 1, NEQNS
        u2(k,j) = ustruct(j)%u(k)
      end do
      u2(k,j) = ustruct(j)%x
    end do
  end subroutine
  
end module mfe_model_type
