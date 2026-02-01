program gwlisa__generate_samples

  use gwlisa__util_kinds,        only : wp
  use gwlisa__simulation_data,  only : noisy_data
  use csv_module,               only : csv_file

  implicit none

  !------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------
  integer, parameter :: nbin    = 500
  integer, parameter :: ncurves = 5000

  ! Bounds for the ewpt parameters
  real(wp), parameter :: pmin(5) = [ &
       10.0_wp,  0.4_wp,  100.0_wp,  110.0_wp, 0.9_wp ]
  real(wp), parameter :: pmax(5) = [ &
       150.0_wp,  0.4_wp, 100.0_wp, 110.0_wp, 0.9_wp ]

  !------------------------------------------------------------
  ! Variables
  !------------------------------------------------------------
  
  type(noisy_data) :: nd
  type(csv_file)   :: f
  logical           :: status_ok

  integer :: i, j, k, ibin
  real(wp) :: pars(5), r

  real(wp) :: logf_min, logf_max, dlogf
  real(wp), allocatable :: f_bin(:), Db_bin(:)
  integer, allocatable :: count_bin(:)
  real(wp), allocatable :: spectrum(:,:)
  real(wp), allocatable :: pars_all(:,:)
  real :: t_start, t_end, t_elapsed, elapsed
  integer :: t0, t1, rate

  !------------------------------------------------------------
  ! Allocate arrays
  !------------------------------------------------------------
  allocate(f_bin(nbin))
  allocate(Db_bin(nbin))
  allocate(count_bin(nbin))
  allocate(spectrum(ncurves, nbin))
  allocate(pars_all(ncurves,5))

  call random_seed()
  !------------------------------------------------------------
  ! Generate noisy data
  !------------------------------------------------------------
  nd = noisy_data(  &
    15.0e0_wp, 3.0e0_wp,  &
    60.0e0_wp, 4.0e-1_wp, 100.0e0_wp,  &
    110.0e0_wp, 0.9e0_wp)

  logf_min = log(nd%f(1))
  logf_max = log(nd%f(nd%num-1))
  dlogf    = (logf_max - logf_min) / real(nbin, wp)

  do j = 1, nbin
    f_bin(j) = exp(logf_min + (real(j,wp) - 0.5_wp) * dlogf)
  end do

  !------------------------------------------------------------
  ! Loop over curves
  !------------------------------------------------------------
  call cpu_time(t_start)
  call system_clock(count_rate=rate)
  call system_clock(t0)
  do k = 1, ncurves

    !----------------------------------------------------------
    ! Draw random parameters
    !----------------------------------------------------------
    do i = 1, 5
      call random_number(r)
      pars(i) = pmin(i) + r * (pmax(i) - pmin(i))
      pars_all(k,i) = pars(i)
    end do

    !----------------------------------------------------------
    ! Generate noisy data with random parameters
    !----------------------------------------------------------
    nd = noisy_data( &
      15.0_wp, 3.0_wp, &
      pars(1), pars(2), pars(3), pars(4), pars(5) )

    Db_bin  = 0.0_wp
    count_bin = 0

    !----------------------------------------------------------
    ! Logarithmic binning
    !----------------------------------------------------------
    do i = 1, nd%num - 1
      ibin = int((log(nd%f(i)) - logf_min) / dlogf) + 1
      ibin = max(1, min(nbin, ibin))

      Db_bin(ibin) = Db_bin(ibin) + nd%Dbar(i)
      count_bin(ibin) = count_bin(ibin) + 1
    end do

    !----------------------------------------------------------
    ! Finalize bin averages
    !----------------------------------------------------------
    do j = 1, nbin
      if (count_bin(j) > 0) then
        Db_bin(j) = Db_bin(j) / count_bin(j)
      end if
    end do

    spectrum(k,:) = Db_bin(:)
    write(*,'(A,I0,A)') 'Generating curve ', k
  end do
  call cpu_time(t_end)
  call system_clock(t1)
  elapsed = real(t1 - t0) / real(rate)
  write(*,'(A,F8.2,A)') 'Elapsed clock time: ', elapsed, ' s'
  t_elapsed = t_end - t_start
  write(*,'(A,F6.2,A)') 'Elapsed CPU time: ', t_elapsed, ' s'
  !------------------------------------------------------------
  ! Write CSV
  !------------------------------------------------------------
  call f%initialize(verbose=.true.)
  call f%open( &
    "data/ewpt_general/signals_update.csv", &
    n_cols = 5 + nbin, status_ok=status_ok )

  ! Header
  call f%add(["p1","p2","p3","p4","p5"])
  do j = 1, nbin
    call f%add("S"//trim(adjustl(to_string(j))))
  end do
  call f%next_row()

  ! Data
  do k = 1, ncurves
    do i = 1, 5
      call f%add(pars_all(k,i), real_fmt="(es15.5)")
    end do
    do j = 1, nbin
      call f%add(spectrum(k,j), real_fmt="(es15.5)")
    end do
    call f%next_row()
  end do

  call f%close(status_ok)

  write(*,*) "Random signal dataset written to random_signals.csv"
block
  type(csv_file) :: ff
  logical :: ok2

  call ff%initialize(verbose=.true.)
  call ff%open( &
    "data/ewpt_general/frequencies.csv", &
    n_cols = 2, status_ok=ok2 )

  call ff%add(["bin","f  "])
  call ff%next_row()

  do j = 1, nbin
    call ff%add(j)
    call ff%add(f_bin(j), real_fmt="(es15.5)")
    call ff%next_row()
  end do

  call ff%close(ok2)
end block
contains

  function to_string(i) result(str)
    integer, intent(in) :: i
    character(len=16) :: str
    write(str,'(i0)') i
  end function to_string

end program gwlisa__generate_samples
