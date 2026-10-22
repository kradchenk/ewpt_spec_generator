program gwlisa__app_calibration_runner
  use gwlisa__util_kinds, only : wp
  use gwlisa__simulation_binning, only : likelihood_calibration
  use csv_module, only : csv_file
  implicit none

  integer, parameter :: R = 10000
  real(wp), allocatable :: bk(:)
  real(wp), allocatable :: Nk(:)
  real(wp), allocatable :: Vark(:)
  integer :: nbin
  integer, allocatable :: ink(:)

  call likelihood_calibration(bk, Nk, Vark=Vark, nbin=nbin, ink=ink, R=R)

  call store(nbin, bk, Nk, Vark, ink)

  contains

    subroutine store(nbin, bk, Nk, Vk, ink)
      integer, intent(in) :: nbin
      real(wp), intent(in) :: bk(nbin)
      real(wp), intent(in) :: Nk(nbin)
      real(wp), intent(in) :: Vk(nbin)
      integer, intent(in) :: ink(nbin)

      type(csv_file) :: f
      logical :: status_ok
      integer :: i

      call f%initialize(verbose=.true.)
      call f%open("data/calibration/data.csv", n_cols=4, status_ok=status_ok)
      call f%add(["bk", "Nk", "Vk", "nk"])
      call f%next_row()
      do i = 1, nbin
        call f%add(bk(i), real_fmt="(es15.5)")
        call f%add(Nk(i), real_fmt="(es15.5)")
        call f%add(Vk(i), real_fmt="(es15.5)")
        call f%add(ink(i))
        call f%next_row()
      end do
      call f%close(status_ok)

    end subroutine

end program
