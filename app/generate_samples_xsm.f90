program generate_samples_xsm

  use gwlisa__util_kinds,        only : wp
  use gwlisa__simulation_data,   only : noisy_data
  use csv_module,                only : csv_file

  implicit none

  integer, parameter :: nlog_high = 500
  integer, parameter :: max_rows = 50000
  real(wp), parameter :: f_cut = 1.0e-4_wp

  type(noisy_data) :: nd, nd0
  type(csv_file)   :: f
  logical          :: status_ok

  integer :: i, j, k, ibin
  integer :: n_low, nbin
  integer :: row, n_rows
  integer :: ios, unit_in

  real(wp) :: pars(5), w
  real(wp) :: c1, c2, c3
  real(wp) :: logf_min_high, logf_max, dlogf_high

  real(wp), allocatable :: f_bin(:), Db_bin(:), Vb_bin(:)
  real(wp), allocatable :: sumw(:), sumwD(:), sumwf(:)
  real(wp), allocatable :: spectrum(:,:), var_spectrum(:,:)
  real(wp), allocatable :: pars_all(:,:)
  real(wp), allocatable :: pars_input(:,:)
  real(wp) :: logOmegaGal, logOmegaExt
  real(wp), allocatable :: fg_all(:,:)

  real :: t_start, t_end, t_elapsed, elapsed
  integer :: t0, t1, rate

  integer :: ncurves
  character(len=512) :: line

  !------------------------------------------------------------
  ! Generate reference data
  !------------------------------------------------------------
  nd0 = noisy_data( &
    15.0_wp, 3.0_wp, &
    60.0_wp, 4.0e-1_wp, 100.0_wp, &
    110.0_wp, 0.9_wp, &
    logOmegaGal=-7.84_wp, &
    logOmegaExt=-12.38_wp )

  !------------------------------------------------------------
  ! Count low-frequency points
  !------------------------------------------------------------
  n_low = 0
  do i = 1, nd0%num - 1
    if (nd0%f(i) <= f_cut) then
      n_low = i
    else
      exit
    end if
  end do

  if (n_low <= 0) stop "ERROR: no low-frequency points found."

  nbin = n_low + nlog_high

  logf_min_high = log(f_cut)
  logf_max      = log(nd0%f(nd0%num-1))
  dlogf_high    = (logf_max - logf_min_high) / real(nlog_high, wp)

  !------------------------------------------------------------
  ! Allocate arrays
  !------------------------------------------------------------
  allocate(f_bin(nbin))
  allocate(Db_bin(nbin))
  allocate(Vb_bin(nbin))

  allocate(sumw(nlog_high))
  allocate(sumwD(nlog_high))
  allocate(sumwf(nlog_high))

  !------------------------------------------------------------
  ! Build frequency grid
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
    if (sumw(j) > 0.0_wp) then
      f_bin(n_low + j) = sumwf(j) / sumw(j)
    else
      f_bin(n_low + j) = exp(logf_min_high + (real(j,wp)-0.5_wp)*dlogf_high)
    end if
  end do

  !------------------------------------------------------------
  ! Read input CSV
  !
  ! combined.csv has a header.
  !
  ! Column mapping:
  !   pars(1) = column 5
  !   pars(2) = column 6
  !   pars(3) = column 7
  !   pars(4) = 110
  !   pars(5) = 0.8
  !------------------------------------------------------------
  open(newunit=unit_in, file="cxsm_uniform_50k_params.csv", status="old", action="read", iostat=ios)

  if (ios /= 0) stop "ERROR: could not open combined.csv."

  ! Skip header
  read(unit_in,'(A)',iostat=ios) line
  if (ios /= 0) stop "ERROR: could not read header from combined.csv."

  n_rows = 0

  do while (n_rows < max_rows)

    read(unit_in,'(A)',iostat=ios) line
    if (ios /= 0) exit

    if (len_trim(line) == 0) cycle

    read(line,*,iostat=ios) &
      c1, c2, c3, &
      pars(1), pars(2), pars(3)

    if (ios == 0) then
      n_rows = n_rows + 1
    end if

  end do

  if (n_rows <= 0) stop "ERROR: no valid rows found in combined.csv."

  allocate(pars_input(n_rows,5))

  rewind(unit_in)

  ! Skip header again
  read(unit_in,'(A)',iostat=ios) line

  row = 0

  do while (row < n_rows)

    read(unit_in,'(A)',iostat=ios) line
    if (ios /= 0) exit

    if (len_trim(line) == 0) cycle

    read(line,*,iostat=ios) &
      c1, c2, c3, &
      pars(1), pars(2), pars(3)

    if (ios == 0) then
      row = row + 1

      pars(4) = 110.0_wp
      pars(5) = 0.8_wp

      pars_input(row,1:5) = pars(1:5)
    end if

  end do

  close(unit_in)

  write(*,'(A,I0)') "Input rows read: ", n_rows

  !------------------------------------------------------------
  ! Total curves
  !------------------------------------------------------------
  ncurves = n_rows

  allocate(spectrum(ncurves, nbin))
  allocate(var_spectrum(ncurves, nbin))
  allocate(pars_all(ncurves,5))
  allocate(fg_all(ncurves,2))

  !------------------------------------------------------------
  ! Main loop
  !------------------------------------------------------------
  call cpu_time(t_start)
  call system_clock(count_rate=rate)
  call system_clock(t0)

  do k = 1, ncurves

    pars(:) = pars_input(k,:)

    pars_all(k,:) = pars(:)

    logOmegaGal = -7.84_wp  !+ 0.21_wp * randn()
    logOmegaExt = -12.38_wp !+ 0.17_wp * randn()

    fg_all(k,1) = logOmegaGal
    fg_all(k,2) = logOmegaExt

    nd = noisy_data( &
      15.0_wp, 3.0_wp, &
      pars(1), pars(2), pars(3), pars(4), pars(5), &
      logOmegaGal=logOmegaGal, &
      logOmegaExt=logOmegaExt )

    Db_bin(1:n_low) = nd%Dbar(1:n_low)
    Vb_bin(1:n_low) = nd%variance(1:n_low)

    sumw  = 0.0_wp
    sumwD = 0.0_wp

    do j = n_low + 1, nd%num - 1
      if (nd%variance(j) <= 0.0_wp) cycle

      ibin = int((log(nd%f(j)) - logf_min_high) / dlogf_high) + 1
      ibin = max(1, min(nlog_high, ibin))

      w = 1.0_wp / nd%variance(j)
      sumw(ibin)  = sumw(ibin)  + w
      sumwD(ibin) = sumwD(ibin) + w * nd%Dbar(j)
    end do

    do j = 1, nlog_high
      if (sumw(j) > 0.0_wp) then
        Db_bin(n_low + j) = sumwD(j) / sumw(j)
        Vb_bin(n_low + j) = 1.0_wp / sumw(j)
      else
        Db_bin(n_low + j) = 0.0_wp
        Vb_bin(n_low + j) = huge(1.0_wp)
      end if
    end do

    spectrum(k,:)     = Db_bin(:)
    var_spectrum(k,:) = Vb_bin(:)

    if (mod(k,1000) == 0) then
      write(*,'(A,I0,A,I0)') 'Curve ', k, ' / ', ncurves
    end if

  end do

  call cpu_time(t_end)
  call system_clock(t1)

  elapsed = real(t1 - t0) / real(rate)
  write(*,'(A,F8.2,A)') 'Elapsed clock time: ', elapsed, ' s'

  t_elapsed = t_end - t_start
  write(*,'(A,F6.2,A)') 'Elapsed CPU time: ', t_elapsed, ' s'

  !------------------------------------------------------------
  ! Write output CSV: samples_xsm.csv
  !------------------------------------------------------------
  call f%initialize(verbose=.true.)

  call f%open("data/samples_xsm_uniform.csv", &
    n_cols = 7 + 2*nbin, status_ok=status_ok)

  if (.not. status_ok) stop "ERROR: could not open data/samples_xsm.csv."

  call f%add("p1")
  call f%add("p2")
  call f%add("p3")
  call f%add("p4")
  call f%add("p5")

  call f%add("logOmegaGal")
  call f%add("logOmegaExt")
  do j = 1, nbin
    call f%add("S"//trim(adjustl(to_string(j))))
  end do

  do j = 1, nbin
    call f%add("V"//trim(adjustl(to_string(j))))
  end do

  call f%next_row()

  do k = 1, ncurves

    do i = 1, 5
      call f%add(pars_all(k,i), real_fmt="(es15.5)")
    end do
    call f%add(fg_all(k,1), real_fmt="(es15.5)")
    call f%add(fg_all(k,2), real_fmt="(es15.5)")

    do j = 1, nbin
      call f%add(spectrum(k,j), real_fmt="(es15.5)")
    end do

    do j = 1, nbin
      call f%add(var_spectrum(k,j), real_fmt="(es15.5)")
    end do

    call f%next_row()

  end do

  call f%close(status_ok)

  if (.not. status_ok) stop "ERROR: problem closing output CSV."

  !------------------------------------------------------------
  ! Write frequency CSV: freq_xsm.csv
  !------------------------------------------------------------
  block
    type(csv_file) :: ff
    logical :: ok2

    call ff%initialize(verbose=.true.)

    call ff%open( &
      "data/freq_xsm_uniform.csv", &
      n_cols = 2, status_ok=ok2 )

    if (.not. ok2) stop "ERROR: could not open data/freq_xsm.csv."

    call ff%add("bin")
    call ff%add("f")
    call ff%next_row()

    do j = 1, nbin
      call ff%add(j)
      call ff%add(f_bin(j), real_fmt="(es15.5)")
      call ff%next_row()
    end do

    call ff%close(ok2)

    if (.not. ok2) stop "ERROR: problem closing frequency CSV."
  end block

contains
  function randn() result(z)

    real(wp) :: z
    real(wp) :: u1, u2
    real(wp), parameter :: pi = acos(-1.0_wp)

    call random_number(u1)
    call random_number(u2)

    u1 = max(u1, tiny(1.0_wp))

    z = sqrt(-2.0_wp * log(u1)) * cos(2.0_wp * pi * u2)

  end function randn
  function to_string(i) result(str)
    integer, intent(in) :: i
    character(len=16) :: str
    write(str,'(i0)') i
  end function to_string

end program generate_samples_xsm
