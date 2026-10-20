module gwlisa__simulation_binning
  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  use gwlisa__util_progress, only : progress_bar
  implicit none

  private

  real(wp), parameter :: fcut_mod = 1.0e-4_wp
  integer, parameter :: nlog_high_mod  = 500
  integer, parameter :: R_mod  = 100
  real(wp), parameter :: A = 3.0e0_wp
  real(wp), parameter :: P = 15.0e0_wp

  public :: create_binned_spectrum
  public :: likelihood_calibration

  contains

    subroutine likelihood_calibration(bk, Nk, Vark, nbin, ink, R, fcut, nlog_high)
      real(wp), intent(out), allocatable :: bk(:)
      real(wp), intent(out), allocatable :: Nk(:)
      real(wp), intent(out), allocatable, optional :: Vark(:)
      integer, intent(out), optional :: nbin
      integer, intent(out), allocatable, optional :: ink(:)
      integer, intent(in), optional :: R
      real(wp), intent(in), optional :: fcut
      integer, intent(in), optional :: nlog_high

      integer :: R_
      real(wp) :: fcut_
      integer :: nlog_high_
      type(noisy_data) :: nd
      integer :: i
      real(wp), allocatable :: Dbar(:, :)
      integer :: nbin_
      real(wp), allocatable :: Vark_(:)
      integer :: t0, rate

      R_ = R_mod
      if (present(R)) R_ = R

      fcut_ = fcut_mod
      if (present(fcut)) fcut_ = fcut

      nlog_high_ = nlog_high_mod
      if (present(nlog_high)) nlog_high_ = nlog_higH

      ! Just to get number of bins and number of f_i per bin k (ink)
      nd = get_realisation()
      nbin_ = nd%num_binned
      allocate(Dbar(nbin_, R_))
      allocate(bk(nbin_))
      allocate(Nk(nbin_))
      allocate(Vark_(nbin_))
      if (present(ink)) ink = nd%nk_binned

      call system_clock(t0, rate)
      do i = 1, R_
        ! Simulate noise only data -> tiny alpha
        nd = get_realisation()
        Dbar(:, i) = nd%Dbar_binned
        call progress_bar(i, R_, t0, rate)
      end do

      do i = 1, nbin_
        bk(i) = sum(Dbar(i, :)) / R_
        Vark_(i) = sum((Dbar(i, :) - bk(i)) ** 2) / (R_ - 1.0e0_wp)
      end do
      Nk = bk ** 2 / Vark_

      if (present(Vark)) Vark = Vark_
      if (present(nbin)) nbin = nbin_

      contains

        function get_realisation() result(n)
          type(noisy_data) :: n

          integer :: j

          n = noisy_data(P, A, 100.0e0_wp, 1.0e-5_wp, 100.0e0_wp, 100.0e0_wp, 1.0e0_wp)

          ! Normalize to power = 1
          do j = 1, n%num
            n%Dbar(j) = n%Dbar(j) / n%OmegaS_hsq(n%f(j))
            n%variance(j) = n%variance(j) / n%OmegaS_hsq(n%f(j)) ** 2
          end do
          call create_binned_spectrum(n, fcut=fcut_, nlog_high=nlog_high_)

        end function get_realisation

    end subroutine likelihood_calibration

    subroutine create_binned_spectrum(nd, fcut, nlog_high)
      type(noisy_data), intent(inout) :: nd
      real(wp), intent(in), optional :: fcut
      integer, intent(in), optional :: nlog_high

      real(wp) :: fcut_
      integer :: nlog_high_
      integer :: i
      integer :: nlow
      integer :: nbin
      real(wp) :: logf_min_high
      real(wp) :: logf_max
      real(wp) :: dlogf
      real(wp), allocatable :: f_bin(:)
      real(wp), allocatable :: Db_bin(:)
      real(wp), allocatable :: Vb_bin(:)
      integer :: ibin
      real(wp) :: weight
      real(wp), allocatable :: sum_weights(:)
      real(wp), allocatable :: sum_f_weighted(:)
      real(wp), allocatable :: sum_Db_weighted(:)
      integer, allocatable :: nk(:) ! # of f_i in bin k

      fcut_ = fcut_mod
      if (present(fcut)) fcut_ = fcut

      nlog_high_ = nlog_high_mod
      if (present(nlog_high)) nlog_high_ = nlog_high

      allocate(sum_weights(nlog_high_))
      allocate(sum_f_weighted(nlog_high_))
      allocate(sum_Db_weighted(nlog_high_))

      ! Count how many original frequencies satisfy f <= fcut
      nlow = 0
      do i = 1, nd%num - 1
        if (nd%f(i) > fcut_) exit
      end do
      nlow = i - 1
      if (nlow <= 0) error stop  &
        'fcut is below frequency range of data.'
      if (nd%f(nlow) > fcut_) error stop  &
        'Error in determining number of unbinned frequencies (a).'
      if (nd%f(nlow + 1) < fcut_) error stop  &
        'Error in determining number of unbinned frequencies (b).'

      ! Total number of bins
      nbin = nlow + nlog_high_
      allocate(f_bin(nbin))
      allocate(Db_bin(nbin))
      allocate(Vb_bin(nbin))
      allocate(nk(nbin))

      ! dlogf in binned frequency range above fcut
      logf_min_high = log(fcut_)
      logf_max = log(nd%f(nd%num))
      dlogf = (logf_max - logf_min_high) / real(nlog_high_, wp)

      ! Below fcut no binning -> just take over original frequencies
      !                          and simulated data
      f_bin(1 : nlow) = nd%f(1 : nlow)
      Db_bin(1 : nlow) = nd%Dbar(1 : nlow)
      Vb_bin(1 : nlow) = nd%variance(1 : nlow)
      nk(1 : nlow) = 1

      ! Above fcut binning -> average frequenvy per bin i with
      !                       weights w_i = 1 / variance_i
      sum_weights = 0.0e0_wp
      sum_f_weighted = 0.0e0_wp
      sum_Db_weighted = 0.0e0_wp
      nk(nlow + 1 : nbin) = 0
      do i = nlow + 1, nd%num
        if (nd%variance(i) <= 0.0e0_wp) cycle
        ibin = int((log(nd%f(i)) - logf_min_high) / dlogf) + 1
        ibin = max(1, min(nlog_high_, ibin))
        weight = 1.0e0_wp / nd%variance(i)
        sum_weights(ibin) = sum_weights(ibin) + weight
        sum_f_weighted(ibin) = sum_f_weighted(ibin) + weight * nd%f(i)
        sum_Db_weighted(ibin) = sum_Db_weighted(ibin) + weight * nd%Dbar(i)
        ibin = ibin + nlow
        if ((ibin <= nlow) .or. (ibin > nbin)) error stop 'Could not assign f_i to bin.'
        nk(ibin) = nk(ibin) + 1
      end do
      do i = 1, nlog_high_
        f_bin(nlow + i) = sum_f_weighted(i) / sum_weights(i)
        Db_bin(nlow + i) = sum_Db_weighted(i) / sum_weights(i)
        Vb_bin(nlow + i) = 1.0e0_wp / sum_weights(i)
      end do

      if (sum(nk) /= nd%num) error stop 'Not all f_i were binned.'

      nd%num_binned = nbin
      nd%f_binned = f_bin
      nd%Dbar_binned = Db_bin
      nd%variance_binned = Vb_bin
      nd%nk_binned = nk

    end subroutine create_binned_spectrum

end module gwlisa__simulation_binning
