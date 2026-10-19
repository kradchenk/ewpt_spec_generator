program gwlisa__test_binning
  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__simulation_binning, only : create_binned_spectrum
  use csv_module,               only : csv_file
  implicit none

  ! Sensitivity curve
  real(wp), parameter :: A = 3.0e0_wp
  real(wp), parameter :: P = 15.0e0_wp

  ! Random point
  real(wp), parameter :: T = 40.0e0_wp
  real(wp), parameter :: alpha = 0.3e0_wp
  real(wp), parameter :: betaH = 335.0e0_wp
  real(wp), parameter :: geff = 110.0e0_wp
  real(wp), parameter :: vw = 1.0e0_wp
  type(noisy_data) :: nd

  ! Output file
  type(csv_file)   :: f
  logical          :: status_ok
  integer :: i

  write(*, *) "Simulating data..."
  nd = noisy_data(P, A, T, alpha, betaH, geff, vw)

  write(*, *) "Binning data..."
  call create_binned_spectrum(nd)

  write(*, *) "Writing orignal spectrum to original_spectrum.csv..."
  call f%initialize(verbose=.true.)
  call f%open( &
    "data/test/binning/orignal_spectrum.csv",  &
    n_cols=3,  &
    status_ok=status_ok)
  call f%add(["f", "D", "v"])
  call f%next_row()
  do i = 1, nd%num
    call f%add(nd%f(i), real_fmt="(es15.5)")
    call f%add(nd%Dbar(i), real_fmt="(es15.5)")
    call f%add(nd%variance(i), real_fmt="(es15.5)")
    call f%next_row()
  end do
  call f%close(status_ok)

  write(*, *) "Writing binned spectrum to binned_spectrum.csv..."
  call f%initialize(verbose=.true.)
  call f%open( &
    "data/test/binning/binned_spectrum.csv",  &
    n_cols=3,  &
    status_ok=status_ok)
  call f%add(["f", "D", "v"])
  call f%next_row()
  do i = 1, nd%num_binned
    call f%add(nd%f_binned(i), real_fmt="(es15.5)")
    call f%add(nd%Dbar_binned(i), real_fmt="(es15.5)")
    call f%add(nd%variance_binned(i), real_fmt="(es15.5)")
    call f%next_row()
  end do
  call f%close(status_ok)

end program gwlisa__test_binning
