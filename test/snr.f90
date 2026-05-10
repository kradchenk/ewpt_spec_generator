program gwlisa__test_snr

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__signals_soundwave, only : spectrum_sw
  use gwlisa__signals_turbulance, only : spectrum_turb
  use gwlisa__sensitivity_curves, only : sens_curves
  use gwlisa__snr, only : calc_snr

  implicit none

  ! Sensitivity curve
  real(wp), parameter :: A = 3.0e0_wp
  real(wp), parameter :: P = 15.0e0_wp

  real(wp), parameter :: T = 40.0e0_wp
  real(wp), parameter :: alpha = 0.3e0_wp
  real(wp), parameter :: betaH = 335.0e0_wp
  real(wp), parameter :: geff = 110.0e0_wp
  real(wp), parameter :: vw = 1.0e0_wp

  type(noisy_data) :: nd
  type(spectrum_sw) :: sw
  type(spectrum_turb) :: turb
  type(sens_curves) :: sens
  real(wp) :: snr

  nd = noisy_data(P, A, T, alpha, betaH, geff, vw)
  snr = calc_snr(nd)
  write(*,*) "SNR =", snr

  sw = spectrum_sw(T, alpha, betaH, geff, vw)
  turb = spectrum_turb(T, alpha, betaH, geff, vw)
  snr = calc_snr(sw, turb=turb)
  write(*,*) "SNR =", snr

  write(*,*) turb%Omega2_hsq, turb%OmegaGW_hsq(1.0e-3_wp)

  sens = sens_curves(P, A)
  snr = calc_snr(sw, turb=turb, sens=sens)
  write(*,*) "SNR =", snr

  write(*,*) turb%Omega2_hsq, turb%OmegaGW_hsq(1.0e-3_wp)

  snr = calc_snr(  &
    T, alpha, betaH, geff,  &
    vw=vw, P=P, A=A,  &
    include_turbulence=.true.,  &
    include_collision=.true.)
  write(*,*) "SNR with collision =", snr

end program gwlisa__test_snr
