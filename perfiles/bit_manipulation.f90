module bit_manipulation
  ! This module provides functions that permute bits of the binary
  ! representation of a number.
  !
  ! Benjamin D. Wandelt October 1997
  ! edited by E. Hivon, October 2001 to be 'F' compatible
  use healpix_types
  implicit none
  integer(kind=i4b), private, parameter :: magic1=89478485,magic2=178956970   
  private
  public :: swapLSBMSB, invswapLSBMSB, invLSB, invMSB

contains
  !====================================================================
  function swapLSBMSB(i) result(fn_result)
  !====================================================================
  !     Returns i with even and odd bit positions interchanged
  ! Benjamin D. Wandelt October 1997
  !====================================================================
    integer(kind=i4b) :: fn_result
    integer(kind=i4b), intent(in) :: i
    integer(kind=i4b) :: msb,lsb
  !--------------------------------------------------------------------
    lsb=IAND(i,magic1)
    msb=IAND(i,magic2)

    fn_result=msb/2+lsb*2
    return
  end function swapLSBMSB

  !====================================================================
  function invswapLSBMSB(i) result(fn_result)
  !====================================================================
  !     Returns NOT(i) with even and odd bit positions interchanged
  ! Benjamin D. Wandelt October 1997
  !====================================================================
    integer(kind=i4b) :: fn_result
    integer(kind=i4b), intent(in) :: i
    integer(kind=i4b) :: msb,lsb
  !--------------------------------------------------------------------
    lsb=IAND(i,magic1)
    msb=IAND(i,magic2)

    fn_result=NOT(msb/2+lsb*2)
    return
  end function invswapLSBMSB

  !====================================================================
  function invLSB(i) result(fn_result)
  !====================================================================
  !     Returns i with even (0,2,4,...) bits inverted
  ! Benjamin D. Wandelt October 1997
  !====================================================================
    integer(kind=i4b) :: fn_result

    integer(kind=i4b), intent(in) :: i
  !--------------------------------------------------------------------
    fn_result=IEOR(i,magic1)
    return
  end function invLSB

  !====================================================================
  function invMSB(i) result(fn_result)
  !====================================================================
  !     Returns i with odd (1,3,5,...) bits inverted
  ! Benjamin D. Wandelt October 1997
  !====================================================================
    integer(kind=i4b) :: fn_result

    integer(kind=i4b), intent(in) :: i
  !--------------------------------------------------------------------
    fn_result=IEOR(i,magic2)
    return
  end function invMSB
end module bit_manipulation
