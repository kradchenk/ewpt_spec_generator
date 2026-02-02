program gwlisa__generate_samples_bin

  use gwlisa__util_kinds,        only : wp
  use gwlisa__simulation_data,   only : noisy_data
  use csv_module,               only : csv_file

  implicit none

  !------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------
  integer, parameter :: ncurves    = 1
  integer, parameter :: nlog_high  = 1000             ! number of log-spaced intervals above f_cut
  real(wp), parameter :: f_cut     = 1.0e-3_wp        ! keep all original points up to this frequency

  ! Bounds for the ewpt parameters
  real(wp), parameter :: pmin(5) = [ &
       150.0_wp,  0.4_wp,  100.0_wp,  110.0_wp, 0.9_wp ]
  real(wp), parameter :: pmax(5) = [ &
       150.0_wp,  0.4_wp, 100.0_wp, 110.0_wp, 0.9_wp ]

  !------------------------------------------------------------
  ! Variables
  !------------------------------------------------------------
  type(noisy_data) :: nd, nd0
  type(csv_file)   :: f
  logical          :: status_ok

  integer :: i, j, k, ibin, n_low, nbin
  real(wp) :: pars(5), r, w

  real(wp) :: logf_min_high, logf_max, dlogf_high
  real(wp), allocatable :: f_bin(:), Db_bin(:), Vb_bin(:)
  real(wp), allocatable :: sumw(:), sumwD(:), sumwf(:)
  real(wp), allocatable :: spectrum(:,:), var_spectrum(:,:)
  real(wp), allocatable :: pars_all(:,:)
  real :: t_start, t_end, t_elapsed, elapsed
  integer :: t0, t1, rate

  !------------------------------------------------------------
  ! Allocate arrays
  !------------------------------------------------------------
  allocate(f_bin(nbin))
  allocate(Db_bin(nbin))
  allocate(Vb_bin(nbin))

  allocate(sumw(nlog_high))
  allocate(sumwD(nlog_high))
  allocate(sumwf(nlog_high))

  allocate(spectrum(ncurves, nbin))
  allocate(var_spectrum(ncurves, nbin))
  allocate(pars_all(ncurves,5))
  
  call random_seed()

  !------------------------------------------------------------
  ! Generate noisy data
  !------------------------------------------------------------
  nd0 = noisy_data(  &
    15.0_wp, 3.0_wp,  &
    60.0_wp, 4.0e-1_wp, 100.0_wp,  &
    110.0_wp, 0.9_wp )

  ! Count how many original frequencies satisfy f <= f_cut
  n_low = 0
  do i = 1, nd0%num - 1
    if (nd0%f(i) <= f_cut) then
      n_low = i
    else
      exit
    end if
  end do
  if (n_low <= 0) stop

  ! Total number of points
  nbin = n_low + nlog_high

  logf_min_high = log(f_cut)
  logf_max      = log(nd0%f(nd0%num-1))
  dlogf_high    = (logf_max - logf_min_high) / real(nlog_high, wp)


  !------------------------------------------------------------
  ! New frequency grid f_bin(:)
  !  - low freq: original frequencies
  !  - high freq: weighted average frequency per log-bin
  !    with weights w_i = 1/sigma_i^2 = 1/variance_i
  !------------------------------------------------------------
  f_bin(1:n_low) = nd0%f(1:n_low)

  sumw  = 0.0_wp
  sumwf = 0.0_wp

  do i = n_low + 1, nd0%num - 1
    if (nd0%variance(i) <= 0.0_wp) cycle
    ibin = int((log(nd0%f(i)) - logf_min_high) / dlogf_high) + 1
    ibin = max(1, min(nlog_high, ibin))

    w = 1.0_wp / nd0%variance(i)
    sumw(ibin)  = sumw(ibin)  + w
    sumwf(ibin) = sumwf(ibin) + w * nd0%f(i)
  end do

  do j = 1, nlog_high
    f_bin(n_low + j) = sumwf(j) / sumw(j)
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


    Db_bin(1:n_low) = nd%Dbar(1:n_low)
    Vb_bin(1:n_low) = nd%variance(1:n_low)

    !----------------------------------------------------------
    ! High-frequency part: log-binning with inverse-variance weights
    !   Dbar_bin = (sum w_i Dbar_i) / (sum w_i)
    !   sigma_bin = (sum w_i)^(-1/2)  =>  variance_bin = 1 / (sum w_i)
    !----------------------------------------------------------
    sumw  = 0.0_wp
    sumwD = 0.0_wp

    do i = n_low + 1, nd%num - 1
      if (nd%variance(i) <= 0.0_wp) cycle
      ibin = int((log(nd%f(i)) - logf_min_high) / dlogf_high) + 1
      ibin = max(1, min(nlog_high, ibin))

      w = 1.0_wp / nd%variance(i)
      sumw(ibin)  = sumw(ibin)  + w
      sumwD(ibin) = sumwD(ibin) + w * nd%Dbar(i)
    end do

    do j = 1, nlog_high
      Db_bin(n_low + j) = sumwD(j) / sumw(j)
      Vb_bin(n_low + j) = 1.0_wp / sumw(j)
    end do

    spectrum(k,:)     = Db_bin(:)
    var_spectrum(k,:) = Vb_bin(:)

    if (mod(k,100) == 0) write(*,'(A,I0,A,I0)') 'Generating curve ', k, ' / ', ncurves
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
    "data/samples_Tfree_bin.csv", &
    n_cols = 5 + 2*nbin, status_ok=status_ok )

  ! Header
  call f%add(["p1","p2","p3","p4","p5"])
  do j = 1, nbin
    call f%add("S"//trim(adjustl(to_string(j))))
  end do
  do j = 1, nbin
    call f%add("V"//trim(adjustl(to_string(j))))
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
    do j = 1, nbin
      call f%add(var_spectrum(k,j), real_fmt="(es15.5)")
    end do
    call f%next_row()
  end do

  call f%close(status_ok)
  write(*,*) "Random signal dataset written to samples_Tfree_bin.csv"

  block
    type(csv_file) :: ff
    logical :: ok2

    call ff%initialize(verbose=.true.)
    call ff%open( &
      "data/frequencies_bin.csv", &
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


end program gwlisa__generate_samples_bin