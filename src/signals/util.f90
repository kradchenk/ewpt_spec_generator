module gwlisa__signals_util

  use gwlisa__util_kinds, only : wp
  use gwlisa__util_constants, only : pi
  use gwlisa__util_constants, only : adiab_ratio

  implicit none

  private

  public calc_kappa
  public calc_Ubarf
  public calc_HxRx
  public calc_HxtauSH
  public calc_Hx
  public calc_HxtauSW
  public calc_gxFac

contains

  ! Source: [1004.4187]
  function calc_kappa(alpha, vw, cs) result(kappa)

    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: vw
    real(wp), intent(in) :: cs
    real(wp) :: kappa

    real(wp) :: kapA
    real(wp) :: kapB
    real(wp) :: kapC
    real(wp) :: kapD
    real(wp) :: xiJ
    real(wp) :: deltaK

    if (alpha > 1.0e-10_wp) then

      kapA = vw ** (6.0e0_wp / 5.0e0_wp) * 6.9e0_wp * alpha /  &
        (1.36e0_wp - 0.037e0_wp * sqrt(alpha) + alpha)

      kapB = alpha ** (2.0e0_wp / 5.0e0_wp) /  &
        (0.017e0_wp + (0.997e0_wp + alpha) ** (2.0e0_wp / 5.0e0_wp))

      kapC = sqrt(alpha) / (0.135e0_wp + sqrt(0.98e0_wp + alpha))

      kapD = alpha / (0.73e0_wp + 0.083e0_wp *  &
        sqrt(alpha) + alpha)

      xiJ = (sqrt((2.0e0_wp / 3.0e0_wp) *alpha +  &
        alpha ** 2) + sqrt(1.0e0_wp / 3.0e0_wp)) /  &
        (1.0e0_wp + alpha)

      deltaK = -0.9e0_wp * log(sqrt(alpha) /  &
        (1.0e0_wp + sqrt(alpha)))

      if (vw < cs) then

        kappa  = cs ** (11.0e0_wp / 5.0e0_wp) * kapA * kapB /  &
          ((cs ** (11.0e0_wp / 5.0e0_wp) -  &
          vw ** (11.0e0_wp / 5.0e0_wp)) * kapB +  &
          vw * cs ** (6.0e0_wp / 5.0e0_wp) * kapA)

      else if (vw >= xiJ) then

        kappa = (xiJ - 1.0e0_wp) ** 3 * xiJ ** (5.0e0_wp / 2.0e0_wp) *  &
          vw ** (-5.0e0_wp / 2.0e0_wp) * kapC * kapD /  &
          (((xiJ - 1.0e0_wp) ** 3 - (vw - 1.0e0_wp) ** 3) *  &
          xiJ ** (5.0e0_wp / 2.0e0_wp) * kapC +  &
          (vw - 1.0e0_wp) ** 3 * kapD)

      else

          kappa = kapB + (vw - cs) * deltaK  &
            +((vw - cs) ** 3 / ((xiJ - cs) ** 3)) *  &
            (kapC - kapB - (xiJ - cs) * deltaK)

      end if

    else

        kappa = alpha ** 2 / alpha
        kappa = kappa / (0.73e0_wp + 0.083e0_wp *  &
          sqrt(alpha) + alpha)

    end if

  end function calc_kappa

  function calc_Ubarf(kappa, alpha) result(Ubarf)

    real(wp), intent(in) :: kappa
    real(wp), intent(in) :: alpha
    real(wp) :: Ubarf

    Ubarf = sqrt((0.6_wp * kappa * alpha /  &
      (1.0e0_wp + alpha)) / adiab_ratio)

  end function calc_Ubarf

  function calc_HxRx(betaH, vw, cs) result(HxRx)

    real(wp), intent(in) :: betaH
    real(wp), intent(in) :: vw
    real(wp), intent(in) :: cs
    real(wp) :: HxRx

    HxRx = (8.0e0_wp * pi) ** (1.0e0_wp / 3.0e0_wp) *  &
      max(vw, cs) / betaH

  end function calc_HxRx

  function calc_HxtauSH(HxRx, Ubarf) result(HxtauSH)

    real(wp), intent(in) :: HxRx
    real(wp), intent(in) :: Ubarf
    real(wp) :: HxtauSH

    HxtauSH = HxRx / Ubarf

  end function calc_HxtauSH

  function calc_Hx(Tx, gx) result(Hx)

    real(wp), intent(in) :: Tx
    real(wp), intent(in) :: gx
    real(wp) :: Hx

    Hx = 1.65e-5_wp * (Tx / 100.e0_wp) *  &
      (gx / 100.e0_wp) ** (1.e0_wp / 6.e0_wp)

  end function calc_Hx

  function calc_HxtauSW(HxtauSH) result(HxtauSW)

    real(wp), intent(in) :: HxtauSH
    real(wp) :: HxtauSW

    HxtauSW = min(1.0, HxtauSH)

  end function calc_HxtauSW

  function calc_gxFac(gx) result(gxFac)

    real(wp), intent(in) :: gx
    real(wp) :: gxFac

    gxFac = (100.e0_wp/ gx) ** (1.0e0_wp / 3.0e0_wp)

  end function calc_gxFac

end module gwlisa__signals_util
