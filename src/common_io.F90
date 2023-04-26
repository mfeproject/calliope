!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module common_io

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfe_constants
  use,intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private

  public :: read_tagged_data, abort, info, element_info

  interface read_tagged_data
     module procedure read_tagged_rval, read_tagged_rvec, read_tagged_rmtx, &
                      read_tagged_ival, read_tagged_ivec, read_tagged_char
  end interface

  integer, public :: input_unit, log_unit, out_unit

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  READ_TAGGED_DATA (generic)
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_tagged_rval(data, desc)

    real(r8), intent(out) :: data
    character(*), intent(in), optional :: desc

    integer :: ios
    character(16) :: tag

    read(input_unit,*,iostat=ios) tag, data
    if (ios /= 0) call read_tagged_error(ios, trim(tag))

    if (present(desc)) write(log_unit,'(t3,a,t32,"=",es12.4)') desc, data

  end subroutine read_tagged_rval


  subroutine read_tagged_rvec(data, desc)

    real(r8), intent(out) :: data(:)
    character(*), intent(in), optional :: desc

    integer :: ios
    character(16) :: tag

    read(input_unit,*,iostat=ios) tag, data
    if (ios /= 0) call read_tagged_error(ios, trim(tag))

    if (present(desc)) write(log_unit,'(t3,a,t32,"=",(t33,4es12.4))') desc, data

  end subroutine read_tagged_rvec


  subroutine read_tagged_rmtx(data, desc)

    real(r8), intent(out) :: data(:,:)
    character(*), intent(in), optional :: desc

    integer :: ios, j
    character(16) :: tag

    read(input_unit,*,iostat=ios) tag, data
    if (ios /= 0) call read_tagged_error(ios, trim(tag))

    if (present(desc)) then
      write(log_unit,'(t3,a,t32,"=",(t33,4es12.4))') desc, data(:,1)
      do j = 2, size(data,2)
        write(log_unit,'(t32,"=",(t33,4es12.4))') data(:,j)
      end do
    end if

  end subroutine read_tagged_rmtx


  subroutine read_tagged_ival(data, desc)

    integer, intent(out) :: data
    character(*), intent(in), optional :: desc

    integer :: ios
    character(16) :: tag

    read(input_unit,*,iostat=ios) tag, data
    if (ios /= 0) call read_tagged_error(ios, trim(tag))

    if (present(desc)) write(log_unit,'(t3,a,t32,"=",i6)') desc, data

  end subroutine read_tagged_ival


  subroutine read_tagged_ivec(data, desc)

    integer, intent(out) :: data(:)
    character(*), intent(in), optional :: desc

    integer :: ios
    character(16) :: tag

    read(input_unit,*,iostat=ios) tag, data
    if (ios /= 0) call read_tagged_error(ios, trim(tag))

    if (present(desc)) write(log_unit,'(t3,a,t32,"=",(t34,i5,7(tr1,i5)))') desc, data

  end subroutine read_tagged_ivec


  subroutine read_tagged_char(data, desc)

    character(*), intent(out) :: data
    character(*), intent(in), optional :: desc

    integer :: ios
    character(16) :: tag

    read(input_unit,*,iostat=ios) tag, data
    if (ios /= 0) call read_tagged_error(ios, trim(tag))

    if (present(desc)) write(log_unit,'(t3,a,t32,"=",a)') desc, trim(data)

  end subroutine read_tagged_char


  subroutine read_tagged_error(iostat, tag)

    integer, intent(in) :: iostat
    character(*), intent(in) :: tag

    character(32) :: filename

    inquire(unit=input_unit,name=filename)

    if (iostat < 0) then

      call abort([log_unit, error_unit], &
      'End-of-file encountered while reading data from ' // trim(filename) // '.')

    else

      call abort([log_unit, error_unit], &
          'Error reading data from ' // trim(filename) // &
          '.  Part of last line: ' // tag // '.')

    end if

  end subroutine read_tagged_error

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  ABORT
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine abort(units, mesg)

    integer, intent(in) :: units(:)
    character(*), intent(in) :: mesg

    integer :: j, unit

    do j = 1, size(units)
      unit = units(j)
      write(unit,'(3a)') '** ', mesg, '  Aborting.'
    end do

    stop 1

  end subroutine abort

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  INFO
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine info(units, mesg)

    integer, intent(in) :: units(:)
    character(*), intent(in) :: mesg

    integer :: j

    do j = 1, size(units)
      write(units(j),'(2a)') '** ', mesg
    end do

  end subroutine info


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  ELEMENT_INFO
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine element_info(mesg, n)

    character(*), intent(in) :: mesg
    integer, intent(in) :: n

    write(log_unit,'(/3a,i0,a)') '** ', mesg, ' (', n, ')'

  end subroutine element_info

end module common_io

