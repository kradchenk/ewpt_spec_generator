module gwlisa__util_kinds

  implicit none

  private

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: qp = selected_real_kind(30,4931)

#ifdef QUAD
  integer, parameter, public :: wp = qp
#else
  integer, parameter, public :: wp = dp
#endif

end module gwlisa__util_kinds
