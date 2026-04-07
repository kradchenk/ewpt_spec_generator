module gwlisa__snr

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data

  implicit none

  private

  real(wp), parameter :: years_default = 3.0e0_wp

  public :: calc_snr

contains

  function calc_snr(nd, years) result(snr)

    type(noisy_data), intent(inout) :: nd
    real(wp), intent(in), optional :: years
    real(wp) :: snr

    real(wp) :: T
    real(wp) :: f1
    real(wp) :: f2
    real(wp) :: df
    real(wp) :: dx
    real(wp) :: dy
    real(wp) :: r
    integer :: i

    if (present(years)) then
      T = years
    else
      T = years_default
    end if
    T = T * 365.0e0_wp * 24.0e0_wp * 60.0e0_wp ** 2

    snr = 0.0e0_wp
    do i = 1, nd%num - 1
      f1 = nd%f(i)
      f2 = nd%f(i + 1)
      df = f2 - f1
      dx = nd%OmegaGW_hsq(f2) + nd%OmegaGW_hsq(f1)
      dy = nd%OmegaS_hsq(f2) + nd%OmegaS_hsq(f1)
      r = dx ** 2 / dy ** 2
      snr = snr + df * r
    end do
    snr = T * snr
    snr = sqrt(snr)

  end function calc_snr

end module gwlisa__snr
