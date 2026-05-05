program gwlisa__app_copa_reconstr

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__signals_soundwave, only : spectrum_sw
  use gwlisa__signals_turbulance, only : spectrum_turb
  use gwlisa__sensitivity_curves, only : sens_curves
  use gwlisa__app_copa_parameters
  use copa__parallel_sampler, only : run_parallel_sampler
  use copa__prior_functions, only : uniform_prior
  use copa__prior_functions, only : log_uniform_prior
  use copa__store, only : store_chains
  use copa__store, only : store_log_probs
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(noisy_data) :: nd_injected
  type(sens_curves) :: sens
  integer :: prior_switch ! 0 = uniform, 1 = log_uniform

  call init_chisq()
  call infer_parameters_copa(0)
  call infer_parameters_copa(1)

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

    if (prior_switch == 0) then
      call uniform_prior(theta, theta_lower, theta_upper, logp)
    else
      call log_uniform_prior(theta, theta_lower, theta_upper, logp)
    end if

  end subroutine log_prior

  subroutine infer_parameters_copa(prior)

    integer, intent(in) :: prior

    real(wp) :: ranges(2, 3)
    real(wp), allocatable :: walkers(:, :, :)
    real(wp), allocatable :: chains(:, :, :, :)
    real(wp), allocatable :: log_probs(:, :, :)

    character(len=:), allocatable :: filepath_chains
    character(len=:), allocatable :: filepath_loglikes

    ranges(1, :) = theta_lower
    ranges(2, :) = theta_upper

    prior_switch = prior

    call run_parallel_sampler(  &
      3, log_prior, log_like,  &
      nsteps=nsteps,  &
      nthreads=nthreads,  &
      nwalkers=nwalkers,  &
      ranges=ranges,  &
      walkers=walkers,  &
      chains=chains,  &
      log_probs=log_probs)

    if (prior_switch == 0) then
      filepath_chains = "data/copa/reconstr/chains_uniformprior"
      filepath_loglikes = "data/copa/reconstr/loglikes_uniformprior"
    else if (prior_switch == 1) then
      filepath_chains = "data/copa/reconstr/chains_loguniformprior"
      filepath_loglikes = "data/copa/reconstr/loglikes_loguniformprior"
    end if

    call store_chains(  &
      chains,  &
      filepath_chains,  &
      mode='machine',  &
      separate=.true.)

    call store_log_probs(  &
      log_probs,  &
      filepath_loglikes,  &
      mode='machine',  &
      separate=.true.)

  end subroutine infer_parameters_copa

end program gwlisa__app_copa_reconstr
