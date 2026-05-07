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
  real(wp), parameter :: years_to_seconds = 365.0e0_wp * 24.0e0_wp * 60.0e0_wp ** 2
  real(wp), parameter :: A_sens = 3.0e0_wp
  real(wp), parameter :: P_sens = 15.0e0_wp

  public :: calc_snr

  interface calc_snr
    module procedure calc_snr_nd
    module procedure calc_snr_sens
    module procedure calc_snr_paras
  end interface calc_snr

contains

  function calc_snr_nd(nd, years) result(snr)

    type(noisy_data), intent(in) :: nd
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
    T = T * years_to_seconds

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
    T = T * years_to_seconds

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
      dx = sig(f2) + sig(f1)
      dy = sn%OmegaS_hsq(f2) + sn%OmegaS_hsq(f1)
      r = dx ** 2 / dy ** 2
      snr = snr + df * r
    end do
    snr = T * snr
    snr = sqrt(snr)

    contains

      function sig(f) result(y)

        real(wp), intent(in) :: f
        real(wp) :: y

        y = sw%OmegaGW_hsq(f)
        if (add_turb) y = y + turb%OmegaGW_hsq(f)
        if (add_coll) y = y + coll%OmegaGW_hsq(f)

      end function sig

  end function calc_snr_sens

  function calc_snr_paras(  &
    T, alpha, betaH, geff, vw, P, A, years,  &
    include_turbulence, include_collision) result(snr)

    real(wp), intent(in) :: T
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: betaH
    real(wp), intent(in) :: geff
    real(wp), intent(in) :: vw
    real(wp), intent(in), optional :: P
    real(wp), intent(in), optional :: A
    real(wp), intent(in), optional :: years
    logical, intent(in), optional :: include_turbulence
    logical, intent(in), optional :: include_collision
    real(wp) :: snr

    type(spectrum_sw) :: sw
    type(spectrum_turb) :: turb
    type(spectrum_coll) :: coll
    type(sens_curves) :: sens
    logical :: add_turb
    logical :: add_coll
    real(wp) :: pp
    real(wp) :: aa
    real(wp) :: yy

    if (present(P)) then
      pp = P
    else
      pp = P_sens
    end if

    if (present(A)) then
      aa = A
    else
      aa = A_sens
    end if

    if (present(years)) then
      yy = years
    else
      yy = years_default
    end if

    if (present(include_turbulence)) then
      add_turb = include_turbulence
    else
      add_turb = .true.
    end if

    if (present(include_collision)) then
      add_coll = include_collision
    else
      add_coll = .false.
    end if

    sens = sens_curves(pp, aa)

    sw = spectrum_sw(T, alpha, betaH, geff, vw)
    if (add_turb) turb = spectrum_turb(T, alpha, betaH, geff, vw)
    if (add_coll) coll = spectrum_coll(T, alpha, betaH, geff)

    if ((.not. add_turb) .and. (.not. add_coll)) then
      snr = calc_snr_sens(sw, sens=sens, years=yy)
    else if ((add_turb) .and. (.not. add_coll)) then
      snr = calc_snr_sens(sw, turb=turb, sens=sens, years=yy)
    else if ((.not. add_turb) .and. (add_coll)) then
      snr = calc_snr_sens(sw, coll=coll, sens=sens, years=yy)
    else if ((add_turb) .and. (add_coll)) then
      snr = calc_snr_sens(sw, turb=turb, coll=coll, sens=sens, years=yy)
    end if

  end function calc_snr_paras

end module gwlisa__snr
