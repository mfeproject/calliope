!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1997 Neil N. Carlson
!!
!! This file is part of MFE1 which is released under the MIT license.  See the
!! file LICENSE or visit http://opensource.org/licenses/MIT for details.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module common_io

  use mfe_constants
  use,intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private

  public  :: read_tagged_data, abort, info, element_info

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

    subroutine read_tagged_rval (data, desc)

      real(kind=wp), intent(out) :: data
      character(len=*), intent(in), optional :: desc

      integer :: ios
      character(len=16) :: tag

      read (unit=input_unit, fmt=*, iostat=ios) tag, data
      if (ios /= 0) then
        call read_tagged_error (ios, trim(tag))
      end if

      if (present(desc)) then
        write (unit=log_unit, fmt="(t3,a,t32,'=',es12.4)") desc, data
      end if

    end subroutine read_tagged_rval


    subroutine read_tagged_rvec (data, desc)

      real(kind=wp), dimension(:), intent(out) :: data
      character(len=*), intent(in), optional :: desc

      integer :: ios
      character(len=16) :: tag

      read (unit=input_unit, fmt=*, iostat=ios) tag, data
      if (ios /= 0) then
        call read_tagged_error (ios, trim(tag))
      end if

      if (present(desc)) then
        write (unit=log_unit, fmt="(t3,a,t32,'=',(t33,4es12.4))") desc, data
      end if

    end subroutine read_tagged_rvec


    subroutine read_tagged_rmtx (data, desc)

      real(kind=wp), dimension(:,:), intent(out) :: data
      character(len=*), intent(in), optional :: desc

      integer :: ios, j
      character(len=16) :: tag

      read (unit=input_unit, fmt=*, iostat=ios) tag, data
      if (ios /= 0) then
        call read_tagged_error (ios, trim(tag))
      end if

      if (present(desc)) then
        write (unit=log_unit, fmt="(t3,a,t32,'=',(t33,4es12.4))") desc, data(:,1)
        do j = 2, size(data,2)
          write (unit=log_unit, fmt="(t32,'=',(t33,4es12.4))") data(:,j)
        end do
      end if

    end subroutine read_tagged_rmtx


    subroutine read_tagged_ival (data, desc)

      integer, intent(out) :: data
      character(len=*), intent(in), optional :: desc

      integer :: ios
      character(len=16) :: tag

      read (unit=input_unit, fmt=*, iostat=ios) tag, data
      if (ios /= 0) then
        call read_tagged_error (ios, trim(tag))
      end if

      if (present(desc)) then
        write (unit=log_unit, fmt="(t3,a,t32,'=',i6)") desc, data
      end if

    end subroutine read_tagged_ival


    subroutine read_tagged_ivec (data, desc)

      integer, dimension(:), intent(out) :: data
      character(len=*), intent(in), optional :: desc

      integer :: ios
      character(len=16) :: tag

      read (unit=input_unit, fmt=*, iostat=ios) tag, data
      if (ios /= 0) then
        call read_tagged_error (ios, trim(tag))
      end if

      if (present(desc)) then
        write (unit=log_unit, fmt="(t3,a,t32,'=',(t34,i5,7(tr1,i5)))") desc, data
      end if

    end subroutine read_tagged_ivec


    subroutine read_tagged_char (data, desc)

      character(len=*), intent(out) :: data
      character(len=*), intent(in), optional :: desc

      integer :: ios
      character(len=16) :: tag

      read (unit=input_unit, fmt=*, iostat=ios) tag, data
      if (ios /= 0) then
        call read_tagged_error (ios, trim(tag))
      end if

      if (present(desc)) then
        write (unit=log_unit, fmt="(t3,a,t32,'=',a)") desc, trim(data)
      end if

    end subroutine read_tagged_char


    subroutine read_tagged_error (iostat, tag)

      integer, intent(in) :: iostat
      character(len=*), intent(in) :: tag

      character(len=32) :: filename

      inquire (unit=input_unit, name=filename)

      if (iostat < 0) then

        call abort ( (/ log_unit, error_unit /), &
        "End-of-file encountered while reading data from " // trim(filename) // ".")

      else

        call abort ( (/ log_unit, error_unit /), &
        "Error reading data from " // trim(filename) // &
        ".  Part of last line: " // tag // ".")

      end if

    end subroutine read_tagged_error

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ABORT
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine abort (units, mesg)

      integer, dimension(:), intent(in) :: units
      character(len=*), intent(in) :: mesg

      integer :: j, unit

      do j = 1, size(units)
        unit = units(j)
        write (unit=unit, fmt="(3a)") "** ", mesg, "  Aborting."
      end do

      stop

    end subroutine abort

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  INFO
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine info (units, mesg)

      integer, dimension(:), intent(in) :: units
      character(len=*), intent(in) :: mesg

      integer :: j

      do j = 1, size(units)
        write (unit=units(j), fmt="(2a)") "** ", mesg
      end do

    end subroutine info


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ELEMENT_INFO
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine element_info (mesg, n)

      character(len=*), intent(in) :: mesg
      integer, intent(in) :: n

      write (unit=log_unit, fmt="(/3a,i3.3,a)") "** ", mesg, " (", n, ")"

    end subroutine element_info

end module common_io

