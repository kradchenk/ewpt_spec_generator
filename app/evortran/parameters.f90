module gwlisa__app_evortran_parameters

  use gwlisa__util_kinds, only : wp

  implicit none

  ! Injected signal
  real(wp), parameter :: K_injected = 0.5188e0_wp
  real(wp), parameter :: Tx_injected = 80.94e0_wp
  real(wp), parameter :: gx_injected = 110.0e0_wp
  real(wp), parameter :: betaH_injected = 104.41e0_wp
  real(wp), parameter :: vw_injected = 0.9e0_wp

  ! Sensitivity curve
  real(wp), parameter :: A = 3.0e0_wp
  real(wp), parameter :: P = 15.0e0_wp

  ! Parameter ranges for fit
  real(wp), parameter :: Tx_lims(2) = [1.0e0_wp, 1.0e3_wp]
  real(wp), parameter :: K_lims(2) = [1.0e-10_wp, 1.0e0_wp - 1.0e-10_wp]
  real(wp), parameter :: betaH_lims(2) = [10.0e0_wp, 1.0e3_wp]

  ! Fit configurations
  integer, parameter :: npoints = 40
  integer, parameter :: popsize = 200
  integer, parameter :: max_generations = 40
  integer, parameter :: selection_size = 50
  character(len=5), parameter :: mating = 'blend'
  integer, parameter :: elite_size = 4
  character(len=7), parameter :: mutate = 'uniform'
  character(len=10), parameter :: selection = 'tournament'
  logical, parameter :: verbose = .true.
  real(wp), parameter :: hugepos = 1.0e35_wp

end module gwlisa__app_evortran_parameters
