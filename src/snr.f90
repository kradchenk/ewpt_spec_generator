module gwlisa__snr

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__signals_soundwave, only : spectrum_sw
  use gwlisa__signals_turbulance, only : spectrum_turb
  use gwlisa__signals_collision, only : spectrum_coll
  use gwlisa__sensitivity_curves, only : sens_curves
  use gwlisa__simulation_data, only : construct_frequencies
  use gwlisa__simulation_data, only : construct_frequencies

  implicit none

  private

  real(wp), parameter :: years_default = 3.0e0_wp
  real(wp), parameter :: A_sens = 3.0e0_wp
  real(wp), parameter :: P_sens = 15.0e0_wp

  public :: calc_snr

  interface calc_snr
    module procedure calc_snr_nd
    module procedure calc_snr_sens
  end interface calc_snr

contains

  function calc_snr_nd(nd, years) result(snr)

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

  end function calc_snr_nd

  function calc_snr_sens(sw, turb, coll, sens, years) result(snr)

    type(spectrum_sw), intent(in) :: sw
    type(spectrum_turb), intent(in), optional :: turb
    type(spectrum_coll), intent(in), optional :: coll
    type(sens_curves), intent(in), optional :: sens
    real(wp), intent(in), optional :: years
    real(wp) :: snr

    real(wp) :: T
    real(wp), allocatable :: f(:)
    integer :: i
    real(wp) :: f1
    real(wp) :: f2
    real(wp) :: df
    real(wp) :: dx
    real(wp) :: dy
    real(wp) :: r
    logical :: add_turb
    logical :: add_coll
    type(sens_curves) :: sn

    if (present(years)) then
      T = years
    else
      T = years_default
    end if
    T = T * 365.0e0_wp * 24.0e0_wp * 60.0e0_wp ** 2

    if (present(turb)) then
      add_turb = .true.
    else
      add_turb = .false.
    end if

    if (present(coll)) then
      add_coll = .true.
    else
      add_coll = .false.
    end if

    if (present(sens)) then
      sn = sens
    else
      sn = sens_curves(P_sens, A_sens)
    end if

    f = construct_frequencies()

    snr = 0.0e0_wp
    do i = 1, size(f) - 1
      f1 = f(i)
      f2 = f(i + 1)
      df = f2 - f1
      dx = signal(f2) + signal(f1)
      dy = sn%OmegaS_hsq(f2) + sn%OmegaS_hsq(f1)
      r = dx ** 2 / dy ** 2
      snr = snr + df * r
    end do
    snr = T * snr
    snr = sqrt(snr)

    contains

      function signal(f) result(y)

        real(wp), intent(in) :: f
        real(wp) :: y

        y = sw%OmegaGW_hsq(f)
        if (add_turb) y = y + turb%OmegaGW_hsq(f)
        if (add_coll) y = y + coll%OmegaGW_hsq(f)

      end function signal

  end function calc_snr_sens

end module gwlisa__snr
