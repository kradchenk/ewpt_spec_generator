module gwlisa__app_copa_parameters

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
  real(wp), parameter :: theta_lower(3) = [Tx_lims(1), K_lims(1), betaH_lims(1)]
  real(wp), parameter :: theta_upper(3) = [Tx_lims(2), K_lims(2), betaH_lims(2)]

  ! Fit configurations
  integer, parameter :: nthreads = 8
  integer, parameter :: nsteps = 400
  integer, parameter :: nwalkers = 20
  real(wp), parameter :: hugepos = 1.0e35_wp

end module gwlisa__app_copa_parameters
