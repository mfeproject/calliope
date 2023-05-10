module string_utilities

  implicit none
  private

  public :: upper_case, lower_case, i_to_c

contains

  elemental function upper_case(s) result(ucs)
    character(*), intent(in) :: s
    character(len(s)) :: ucs
    integer :: i
    do i = 1, len(s)
      if (iachar(s(i:i)) >= iachar('a') .and. iachar(s(i:i)) <= iachar('z')) then
        ucs(i:i) = achar(iachar(s(i:i)) - iachar('a') + iachar('A'))
      else
        ucs(i:i) = s(i:i)
      end if
    end do
  end function

  elemental function lower_case(s) result(lcs)
    character(*), intent(in) :: s
    character(len(s)) :: lcs
    integer :: i
    do i = 1, len(s)
      if (iachar(s(i:i)) >= iachar('A') .and. iachar(s(i:i)) <= iachar('Z')) then
        lcs(i:i) = achar(iachar(s(i:i)) - iachar('A') + iachar('a'))
      else
        lcs(i:i) = s(i:i)
      end if
    end do
  end function

  pure function i_to_c(n) result(s)
    integer, intent(in) :: n
    character(:), allocatable :: s
    integer :: m, pos, digit
    character(range(n)+1) :: buffer
    m = sign(n, -1) ! fold to the larger negative range
    pos = len(buffer) + 1
    do while (m < 0)
      pos = pos - 1
      digit = 1 - mod(m, 10)
      buffer(pos:pos) = '0123456789'(digit:digit)
      m = m/10
    end do
    if (n < 0) then
      pos = pos - 1
      buffer(pos:pos) = '-'
    end if
    s = buffer(pos:)
  end function

end module string_utilities

