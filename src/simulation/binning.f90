module gwlisa__simulation_binning
  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_data, only : noisy_data
  implicit none

  private

  real(wp), parameter :: fcut_mod = 1.0e-4_wp
  integer, parameter :: nlog_high_mod  = 500

  public :: create_binned_spectrum

  contains

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

      fcut_ = fcut_mod
      if (present(fcut)) fcut_ = fcut

      nlog_high_ = nlog_high_mod
      if (present(nlog_high)) nlog_high_ = nlog_higH

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

      ! dlogf in binned frequency range above fcut
      logf_min_high = log(fcut_)
      logf_max = log(nd%f(nd%num))
      dlogf = (logf_max - logf_min_high) / real(nlog_high_, wp)

      ! Below fcut no binning -> just take over original frequencies
      !                          and simulated data
      f_bin(1 : nlow) = nd%f(1 : nlow)
      Db_bin(1 : nlow) = nd%Dbar(1 : nlow)
      Vb_bin(1 : nlow) = nd%variance(1 : nlow)

      ! Above fcut binning -> average frequenvy per bin i with
      !                       weights w_i = 1 / variance_i
      sum_weights = 0.0e0_wp
      sum_f_weighted = 0.0e0_wp
      sum_Db_weighted = 0.0e0_wp
      do i = nlow + 1, nd%num
        if (nd%variance(i) <= 0.0e0_wp) cycle
        ibin = int((log(nd%f(i)) - logf_min_high) / dlogf) + 1
        ibin = max(1, min(nlog_high_, ibin))
        weight = 1.0e0_wp / nd%variance(i)
        sum_weights(ibin) = sum_weights(ibin) + weight
        sum_f_weighted(ibin) = sum_f_weighted(ibin) + weight * nd%f(i)
        sum_Db_weighted(ibin) = sum_Db_weighted(ibin) + weight * nd%Dbar(i)
      end do
      do i = 1, nlog_high_
        f_bin(nlow + i) = sum_f_weighted(i) / sum_weights(i)
        Db_bin(nlow + i) = sum_Db_weighted(i) / sum_weights(i)
        Vb_bin(nlow + i) = 1.0e0_wp / sum_weights(i)
      end do

      nd%num_binned = nbin
      nd%f_binned = f_bin
      nd%Dbar_binned = Db_bin
      nd%variance_binned = Vb_bin

    end subroutine create_binned_spectrum

end module gwlisa__simulation_binning
