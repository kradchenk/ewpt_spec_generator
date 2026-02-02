program gwlisa__ewpt_signal

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__signals_soundwave, only : spectrum_sw
  use gwlisa__signals_turbulance, only : spectrum_turb
  use csv_module, only : csv_file

  implicit none

  type(noisy_data) :: nd
  type(spectrum_sw) :: sw
  type(spectrum_turb) :: turb
  real(wp) :: fitness
  type(csv_file) :: f
  type(csv_file) :: g
  logical :: status_ok
  integer :: i
  integer :: j
  character(len=3) :: js

  call f%initialize(verbose=.true.)
  call f%open(  &
    "data/simulation.csv",  &
    n_cols=4, status_ok=status_ok)
  call f%add(["f0", "Db", "si", "se"])
  call f%next_row()
  nd = noisy_data(  &
    15.0e0_wp, 3.0e0_wp,  &
    60.0e0_wp, 4.0e-1_wp, 100.0e0_wp,  &
    110.0e0_wp, 0.9e0_wp)
  do i = 1, nd%num - 1
    call f%add(nd%f(i), real_fmt="(6es15.5)")
    call f%add(nd%Dbar(i), real_fmt="(6es15.5)")
    call f%add(nd%OmegaGW_hsq(nd%f(i)), real_fmt="(6es15.5)")
    call f%add(nd%OmegaS_hsq(nd%f(i)), real_fmt="(6es15.5)")
    call f%next_row()
  end do
  call f%close(status_ok)


end program gwlisa__ewpt_signal