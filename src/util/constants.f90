module gwlisa__util_constants

  use gwlisa__util_kinds, only : wp

  implicit none

  private

  real(wp), public, parameter :: pi = 4.0e0_wp * atan(1.0e0_wp)
  real(wp), public, parameter :: cs_default = sqrt(1.0e0_wp / 3.0e0_wp)
  real(wp), public, parameter :: adiab_ratio = 4.0e0_wp / 3.0e0_wp
  real(wp), public, parameter :: ZP = 10.0e0_wp
  real(wp), public, parameter :: h = 0.678e0_wp
  real(wp), public, parameter :: c = 2.998e8_wp
  real(wp), public, parameter :: H0 = 1.0e2_wp * h * 1.0e3_wp / 3.086e22_wp

end module gwlisa__util_constants
