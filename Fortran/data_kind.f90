!*****************************************************************************!
!**************************     Module data_kind      ************************!
!*****************************************************************************!
module data_kind

  implicit none

  private

  integer, parameter, public :: real_4   = selected_real_kind(6)
  integer, parameter, public :: real_8   = selected_real_kind(13)
  integer, parameter, public :: real_16  = selected_real_kind(24)
  integer, parameter, public :: integer_2 = selected_int_kind(4)
  integer, parameter, public :: integer_4 = selected_int_kind(9)
  integer, parameter, public :: integer_8 = selected_int_kind(18)

  integer, parameter, public :: integer_typ = integer_8
  integer, parameter, public :: real_typ = real_8

end module data_kind

