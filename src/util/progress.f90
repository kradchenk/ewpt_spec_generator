module gwlisa__util_progress
  use gwlisa__util_kinds, only : wp
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none

  private

  public :: progress_bar

  contains

    subroutine progress_bar(i, n, t0, rate)
      integer, intent(in) :: i, n
      integer, intent(in) :: t0, rate

      integer, parameter :: width = 40
      integer :: t
      integer :: filled
      real :: frac, elapsed, eta
      character(len=width) :: bar

      call system_clock(t)
      frac    = real(i) / real(n)
      filled  = int(frac * width)
      bar     = repeat('#', filled) // repeat('-', width - filled)
      elapsed = real(t - t0) / real(rate)
      eta     = elapsed * (1.0 - frac) / frac

      write(output_unit, '(a,"[",a,"] ",i0,"/",i0," (",f5.1,"%)  elapsed ",a,"  ETA ",a)', &
            advance='no') achar(13), bar, i, n, 100.0*frac, fmt_time(elapsed), fmt_time(eta)
      if (i == n) write(output_unit, *)   ! newline when finished
      flush(output_unit)

    end subroutine progress_bar

  function fmt_time(s) result(str)
    real, intent(in) :: s
    character(len=8) :: str

    integer :: h, m, sec

    sec = nint(s)
    h = sec / 3600
    m = mod(sec, 3600) / 60
    sec = mod(sec, 60)
    write(str, '(i2.2,":",i2.2,":",i2.2)') h, m, sec

  end function fmt_time

end module gwlisa__util_progress
