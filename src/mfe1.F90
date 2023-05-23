!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  MFE1
!!
!!  This is an implementation of the piecewise-linear Gradient-Weighted
!!  Moving Finite Element Method (GWMFE) for time-dependent systems of
!!  partial differential equations in one space dimension.
!!
!!  Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2023, Neil N. Carlson
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program mfe1

  use,intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use parameter_list_type
  use parameter_list_json
  use mfe_env_type
  use mfe_sim_type
  implicit none

  integer :: i, stat, lun, n, narg
  character(255) :: arg
  character(:), allocatable :: prog, infile, errmsg
  type(parameter_list), pointer :: params
  type(mfe_sim), target :: sim
  type(mfe_env), target :: env

  !! Initialize the simulation run environment
  open(newunit=env%log_unit, file='mfelog', position='rewind', action='write', status='replace')
  open(newunit=env%grf_unit, file='mfegrf', position='rewind', action='write', status='replace')
  call env%log%init([output_unit, env%log_unit], verbosity=VERB_NORMAL)

  call env%timer%start('simulation')
  call env%timer%start('initialization')

  !FIXME: command line parsing is fragile
  !TODO: move command line parsing to a separate procedure

  call get_command_argument(0, arg)
  i = scan(arg, '/', back=.true.)
  prog = trim(arg(i+1:))  ! remove the leading path component, if any

  narg = command_argument_count()
  if (narg < 1) then
    write(error_unit,'(3a)') 'Usage: ', prog, ' [-L <pde-lib-dir>] INFILE'
    stop 1
  end if

  n = 1
  call get_command_argument(n, arg)
  if (arg == '-L') then
    n = n + 1
    call get_command_argument(n, arg)
    i = scan(arg, '/', back=.true.)
    if (i == len_trim(arg)) i = i - 1 ! remove trailing slash, if any
    env%pde_lib_dir = trim(arg) // '/'

    n = n + 1
    call get_command_argument(n, arg)
  else
    env%pde_lib_dir = ''
  end if

  infile = trim(arg)

  if (narg /= n) then
    write(error_unit,'(3a)') 'Usage: ', prog, ' [-L <pde-lib-dir>] INFILE'
    stop 1
  end if

  open(newunit=lun, file=infile, action='read', access='stream')
  call parameter_list_from_json_stream(lun, params, errmsg)
  close(lun)
  if (.not.associated(params)) call env%log%fatal(errmsg)

  call parameter_list_to_json(params, env%log_unit, real_fmt='g0.5')

  call sim%init(env, params, stat, errmsg)
  if (stat /= 0) call env%log%fatal(errmsg)
  call env%timer%stop('initialization')

  call sim%run(stat, errmsg)

  call env%timer%stop('simulation')
  call env%timer%write(output_unit, 2)
  if (stat /= 0) call env%log%fatal(errmsg)

  call env%log%exit

end program mfe1
