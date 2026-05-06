program gwlisa__app_evortran_reconstr

  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__signals_soundwave, only : spectrum_sw
  use gwlisa__signals_turbulance, only : spectrum_turb
  use gwlisa__sensitivity_curves, only : sens_curves
  use gwlisa__app_evortran_parameters
  use evortran__individuals_float, only : individual
  use evortran__prng_rand, only : initialize_rands
  use evortran__evolutions_float, only : evolve_population
  use csv_module, only : csv_file
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(noisy_data) :: nd_injected
  type(sens_curves) :: sens
  integer :: prior_switch ! 0 = uniform, 1 = log_uniform

  call init_chisq()
  call infer_parameters_evortran(npoints)

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

  subroutine calc_chisq(ind, f)

    class(individual), intent(in) :: ind
    real(wp), intent(out) :: f

    real(wp) :: Tx
    real(wp) :: K
    real(wp) :: betaH
    real(wp) :: alpha
    type(spectrum_sw) :: sw
    type(spectrum_turb) :: turb
    integer :: i

    call extract_paras(ind%genes, Tx, K, betaH)
    alpha = calc_alpha(K)
    sw = spectrum_sw(Tx, alpha, betaH, gx_injected, vw_injected)
    turb = spectrum_turb(Tx, alpha, betaH, gx_injected, vw_injected)

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

  subroutine extract_paras(genes, Tx, K, betaH)

    real(wp), intent(in) :: genes(3)
    real(wp), intent(out) :: Tx
    real(wp), intent(out) :: K
    real(wp), intent(out) :: betaH

    Tx = Tx_lims(1) + genes(1) * (Tx_lims(2) - Tx_lims(1))
    K = K_lims(1) + genes(2) * (K_lims(2) - K_lims(1))
    betaH = betaH_lims(1) + genes(3) * (betaH_lims(2) - betaH_lims(1))

  end subroutine extract_paras

  subroutine infer_parameters_evortran(n)

    integer, intent(in) :: n

    type(individual) :: best_ind
    integer :: i
    real(wp) :: Tx(n)
    real(wp) :: K(n)
    real(wp) :: betaH(n)
    real(wp) :: chisq(n)

    call initialize_rands(mode='twister', seed=14)

    do i = 1, n

      write(*,*) "ITERATION: ", i

      best_ind = evolve_population(  &
        popsize, 3, calc_chisq,  &
        max_generations=max_generations,  &
        selection_size=selection_size,  &
        mating=mating,  &
        elite_size=elite_size,  &
        mutate=mutate,  &
        selection=selection,  &
        verbose=verbose)

      call extract_paras(best_ind%genes, Tx(i), K(i), betaH(i))

      chisq(i) = best_ind%get_fitness()

    end do

    call store_points(n, Tx, K, betaH, chisq)

  end subroutine infer_parameters_evortran

  subroutine store_points(n, Tx, K, betaH, chisq)

    integer, intent(in) :: n
    real(wp), intent(in) :: Tx(n)
    real(wp), intent(in) :: K(n)
    real(wp), intent(in) :: betaH(n)
    real(wp), intent(in) :: chisq(n)

    type(csv_file) :: f
    logical :: ok
    integer :: i
    character(len=:), allocatable :: fmt

    call f%initialize(verbose=.true.)
    fmt = "(6es15.5)"

    call f%open(  &
      "data/evortran/reconstr/bestinds.csv",  &
      n_cols=4, status_ok=ok)

    call f%add(["Tx", "K ", "bH", "ch"])

    call f%next_row()

    do i = 1, n
      call f%add(Tx(i), real_fmt=fmt)
      call f%add(K(i), real_fmt=fmt)
      call f%add(betaH(i), real_fmt=fmt)
      call f%add(chisq(i), real_fmt=fmt)
      call f%next_row()
    end do

    call f%close(ok)

  end subroutine store_points

end program gwlisa__app_evortran_reconstr
