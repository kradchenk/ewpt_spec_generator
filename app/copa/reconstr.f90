program gwlisa__app_copa_reconstr

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__signals_soundwave, only : spectrum_sw
  use gwlisa__signals_turbulance, only : spectrum_turb
  use gwlisa__sensitivity_curves, only : sens_curves
  use copa__parallel_sampler, only : run_parallel_sampler
  use copa__prior_functions, only : uniform_prior
  use copa__store, only : store_chains
  use, intrinsic :: ieee_arithmetic

implicit none

  ! Injected signal
  real(wp), parameter :: K_injected = 0.5188e0_wp
  real(wp), parameter :: Tx_injected = 80.94e0_wp
  real(wp), parameter :: gx_injected = 110.0e0_wp
  real(wp), parameter :: betaH_injected = 104.41e0_wp
  real(wp), parameter :: vw_injected = 0.9e0_wp
  type(noisy_data) :: nd_injected
  type(sens_curves) :: sens

  ! Sensitivity curve
  real(wp), parameter :: A = 3.0e0_wp
  real(wp), parameter :: P = 15.0e0_wp

  ! Parameter ranges for fit
  real(wp), parameter :: Tx_lims(2) = [1.0e0_wp, 1.0e3_wp]
  real(wp), parameter :: K_lims(2) = [1.0e-10_wp, 1.0e0_wp - 1.0e-10_wp]
  real(wp), parameter :: betaH_lims(2) = [10.0e0_wp, 1.0e3_wp]
  real(wp), parameter :: theta_lower(3) = [Tx_lims(1), K_lims(1), betaH_lims(1)]
  real(wp), parameter :: theta_upper(3) = [Tx_lims(2), K_lims(2), betaH_lims(2)]

  real(wp), parameter :: hugepos = 1.0e35_wp

! real(wp) :: chisqtest
! real(wp) :: thetatest(3)

  call init_chisq()

! thetatest = [80.93e0_wp, 0.52e0_wp, 104.41e0_wp]
! call calc_chisq(thetatest, chisqtest)
! write(*,*) chisqtest
! thetatest = [81.0e0_wp, 0.5e0_wp, 115.0e0_wp]
! call calc_chisq(thetatest, chisqtest)
! write(*,*) chisqtest

  call infer_parameters_copa()

contains

  pure function calc_alpha(K) result(a)

    real(wp), intent(in) :: K
    real(wp) :: a

    real(wp) :: d

    d = 1.0e0_wp - K

    if (d > 1.0e-10_wp) then
      a = K / d
    else
      a = 1.0e8_wp
    end if

  end function calc_alpha

  subroutine init_chisq()

    real(wp) :: alpha_injected

    alpha_injected = calc_alpha(K_injected)

    nd_injected = noisy_data(  &
      P, A,  &
      Tx_injected, alpha_injected, betaH_injected,  &
      gx_injected, vw_injected)

    sens = sens_curves(P, A)

  end subroutine init_chisq


  subroutine calc_chisq(theta, f)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(out) :: f

    real(wp) :: Tx
    real(wp) :: K
    real(wp) :: betaH
    real(wp) :: gx
    real(wp) :: vw
    real(wp) :: alpha
    type(spectrum_sw) :: sw
    type(spectrum_turb) :: turb
    integer :: i
    logical :: out_of_bounce

    call check_theta_ranges(theta, out_of_bounce)
    if (out_of_bounce) then
      f = hugepos
      return
    end if

    Tx = theta(1)
    K = theta(2)
    betaH = theta(3)
    gx = gx_injected
    vw = vw_injected
    alpha = calc_alpha(K)
    sw = spectrum_sw(Tx, alpha, betaH, gx, vw)
    turb = spectrum_turb(Tx, alpha, betaH, gx, vw)

    f = 0.0e0_wp
    do i = 1, size(nd_injected%Dbar)
      if (nd_injected%f(i) > 1.0e-2_wp) then
        exit
      end if
      f = f + 0.5e0_wp * (  &
        nd_injected%Dbar(i) -  &
        sw%OmegaGW_hsq(nd_injected%f(i)) -  &
        turb%OmegaGW_hsq(nd_injected%f(i)) -  &
        sens%OmegaS_hsq(nd_injected%f(i))) ** 2 /  &
        nd_injected%variance(i)
    end do
    if (ieee_is_nan(f)) then
      f = hugepos
    end if

  end subroutine calc_chisq

  subroutine check_theta_ranges(theta, out_of_bounce)

    real(wp), intent(in) :: theta(:)
    logical, intent(out) :: out_of_bounce

    integer :: i
    integer :: n

    n = size(theta)
    out_of_bounce = .false.
    do i = 1, n
      if ((theta(i) < theta_lower(i)) .or. (theta(i) > theta_upper(i))) then
        out_of_bounce = .true.
        return
      end if
    end do

  end subroutine check_theta_ranges

  subroutine log_like(theta, logl)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(out) :: logl

    real(wp) :: chisq

    call calc_chisq(theta, chisq)
    logl = -0.5e0_wp * chisq

  end subroutine log_like

  subroutine log_prior(theta, logp)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(out) :: logp

    call uniform_prior(theta, theta_lower, theta_upper, logp)

  end subroutine log_prior

  subroutine infer_parameters_copa()

    real(wp) :: ranges(2, 3)
    real(wp), allocatable :: walkers(:, :, :)
    real(wp), allocatable :: chains(:, :, :, :)
    real(wp), allocatable :: log_probs(:, :, :)

    ranges(1, :) = theta_lower
    ranges(2, :) = theta_upper

    call run_parallel_sampler(  &
      3, log_prior, log_like,  &
      nsteps=1000 / 10,  &
      nthreads=16,  &
      ranges=ranges,  &
      walkers=walkers,  &
      chains=chains,  &
      log_probs=log_probs)

    call store_chains(  &
      chains,  &
      "data/copa/reconstr/chains.npy",  &
      mode='machine')

  end subroutine infer_parameters_copa

end program gwlisa__app_copa_reconstr
