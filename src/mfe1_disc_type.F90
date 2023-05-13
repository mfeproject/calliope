#include "f90_assert.fpp"

module mfe1_disc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe1_disc_core_type
  use mfe1_vector_type
  use pde_class
  use parameter_list_type
  implicit none
  private

  type, extends(mfe1_disc_core), public :: mfe1_disc
    real(r8) :: fdinc
    class(pde(:)), allocatable :: p
  contains
    procedure :: init
    procedure :: compute_local_residual
    procedure :: pde_rhs
    procedure :: eval_dfdy
  end type

contains

  subroutine init(this, neqns, ncell, params, stat, errmsg)

    use navier_stokes_type

    class(mfe1_disc), intent(out), target :: this
    integer, intent(in) :: neqns, ncell
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call this%mfe1_disc_core%init(neqns, ncell, params, stat, errmsg)
    if (stat /= 0) return

    call params%get('fdinc', this%fdinc, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%fdinc <= 0.0) then
      stat = 1
      errmsg = ' fdinc <= 0.0'
      return
    end if

    !allocate(navier_stokes :: this%p)
    call alloc_pde(this%p)
    call this%p%init(this%mfe1_disc_core, params%sublist('problem'), stat, errmsg)
    if (stat /= 0) return

  end subroutine init

  subroutine compute_local_residual(this, t)
    class(mfe1_disc), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%p%rhs(t)
    call this%reg_rhs
    call this%res_mass_matrix
  end subroutine


  subroutine pde_rhs(this, t)
    class(mfe1_disc), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%p%rhs(t)
  end subroutine


  subroutine eval_dfdy(this, t, stat, errmsg)

    class(mfe1_disc), intent(inout) :: this
    real(r8), intent(in) :: t
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j, k, l
    real(r8) :: rh
    real(r8), allocatable :: r_save(:,:,:)

    rh = 1.0_r8 / this%fdinc

    !! Compute and save the unperturbed residual
    call this%compute_local_residual(t)
    r_save = this%r

    do k = 1, 2
      do i = 1, this%neqns+1
        this%u(i,k,:) = this%u(i,k,:) + this%fdinc
        call this%recompute !TODO: should check for inverted cells (stat below)
        call this%compute_local_residual(t)
        this%u(i,k,:) = this%u(i,k,:) - this%fdinc
        do j = 1, this%ncell
          do l = 1, 2
            this%mtx(:,i,l,k,j) = this%mtx(:,i,l,k,j) + rh * (this%r(:,l,j) - r_save(:,l,j))
          end do
        end do
      end do
    end do

    stat = 0

  end subroutine eval_dfdy

end module
