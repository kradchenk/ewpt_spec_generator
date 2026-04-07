program gwlisa__test_snr

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__snr, only : calc_snr

  implicit none

  type(noisy_data) :: nd
  real(wp) :: snr

  nd = noisy_data(  &
    15.0e0_wp, 3.0e0_wp,  &
    40.0e0_wp, 0.3e0_wp, 335.0e0_wp,  &
    110.0e0_wp, 1.0e0_wp)

  snr = calc_snr(nd)

  write(*,*) "SNR =", snr

end program gwlisa__test_snr
