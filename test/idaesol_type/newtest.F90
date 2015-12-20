program main

  use kinds
  use nodal_FEM_1D
  use idaesol_type
  use parameter_list_type
  implicit none
  
  integer :: nstep, status, nnode
  real(r8) :: tout, t, h
  real(r8), allocatable :: u(:)
  type(idaesol) :: s
  class(idaesol_model), pointer :: bdf2_prob
  type(ht_model), pointer :: prob
  type(parameter_list) :: bdf_params
  
  nnode = 201
  allocate(u(nnode))
  
  allocate(prob)
  call prob%init (nnode, x0=0.0_r8, x1=1.0_r8, u0=0.0_r8, u1=0.0_r8, pnum=2)
  bdf2_prob => prob
  
  prob%d = 0.0002_r8
  prob%rtol = 0.0
  prob%atol = 1.0d-5
  
  !call s%create (bdf2_prob, size(u), mvec=2, ntol=0.01d0)
  call bdf_params%set ('nlk-tol', 0.01d0)
  call bdf_params%set ('nlk-max-vec', 5)
  call s%init (bdf2_prob, bdf_params)
  
  u = sin(4.0_r8*atan(1.0_r8)*prob%x)
  t = 0.0d0
  
  call s%set_initial_state (t, u, prob%udot(t,u))
  call s%set_verbose_stepping (unit=10)
  
  nstep = 10
  tout = 0.2d0
  h = 1.0d-5
  
  call s%integrate (h, status, tout=tout)
  select case (status)
  case (SOLVED_TO_TOUT)
    call s%get_interpolated_state(tout, u)
    call prob%user (tout, u)
  case default
    print *, 'BDF2_STEP_DRIVER returned an unknown status: ', status
    print *, t
  end select
  
  call s%write_metrics (0)
  
end program main
