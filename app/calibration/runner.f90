program gwlisa__app_calibration_runner
  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_binning, only : likelihood_calibration
  implicit none

  integer, parameter :: R = 4    ! Number of realisations of noise
  real(wp), allocatable :: bk(:)
  real(wp), allocatable :: Nk(:)

  call likelihood_calibration(bk, Nk, R=R)

end program
