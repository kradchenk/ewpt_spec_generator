module gwlisa__simulation_data

  use gwlisa__util_kinds, only : wp
  use gwlisa__signals_soundwave, only : spectrum_sw
  use gwlisa__signals_turbulance, only : spectrum_turb
  use gwlisa__signals_collision, only : spectrum_coll
  use gwlisa__sensitivity_curves, only : sens_curves
  use gwlisa__util_constants, only : cs_default
  use evortran__prng_rand, only : initialize_rands
  use evortran__prng_rand, only : are_rands_initializaed
  use evortran__prng_rand, only : Gi => randnormal
  use omp_lib

  implicit none

  private

  real(wp), parameter :: fmin = 3.0e-5_wp
  real(wp), parameter :: fmax = 5.0e-1_wp
  real(wp), parameter :: fdel = 1.0e-6_wp
  integer, parameter :: fnum = floor((fmax - fmin) / fdel + 1)
  integer, parameter :: num_chunks = 94

  ! Source: [1906.09244]
  type, public :: noisy_data
    type(spectrum_sw), public :: sw
    type(spectrum_turb), public :: turb
    type(spectrum_coll), public :: coll
    type(sens_curves), public :: sens
    real(wp), public, allocatable :: f(:) ! fi
    real(wp), public, allocatable :: S(:,:) ! Si
    real(wp), public, allocatable :: N(:,:) ! Ni
    real(wp), public, allocatable :: D(:,:) ! Di
    real(wp), public, allocatable :: Dbar(:) ! Dbari
    real(wp), public, allocatable :: variance(:) ! Dbari
    integer, public :: num
    logical, private :: incl_coll
    ! Foreground amplitudes
    real(wp), private :: logOmegaGal
    real(wp), private :: logOmegaExt
  contains
    procedure, public :: OmegaGW_hsq
    procedure, public :: OmegaS_hsq
    procedure, private :: OmegaGal_hsq
    procedure, private :: OmegaExt_hsq
    procedure, private :: construct_signal
  end type noisy_data

  interface noisy_data
    procedure create_noisy_data
  end interface noisy_data

contains

  function create_noisy_data(  &
    P, A,  &
    Tx, alpha, betaH, gx, vw, ep, cs, include_coll, &
    logOmegaGal, logOmegaExt) result(this)

    real(wp), intent(in) :: P
    real(wp), intent(in) :: A
    real(wp), intent(in) :: Tx
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: betaH
    real(wp), intent(in) :: gx
    real(wp), intent(in) :: vw
    real(wp), intent(in), optional :: ep
    real(wp), intent(in), optional :: cs
    logical, intent(in), optional :: include_coll
    real(wp), intent(in), optional :: logOmegaGal
    real(wp), intent(in), optional :: logOmegaExt
    type(noisy_data) :: this

    real(wp) :: css
    real(wp) :: epss

    if (present(logOmegaGal)) then
      this%logOmegaGal = logOmegaGal
    else
      this%logOmegaGal = -7.84_wp
    end if

    if (present(logOmegaExt)) then
      this%logOmegaExt = logOmegaExt
    else
      this%logOmegaExt = -12.38_wp
    end if

    if (present(cs)) then
      css = cs
    else
      css = cs_default
    end if

    if (present(ep)) then
      epss = ep
    else
      epss = 0.5_wp
    end if

    if (present(include_coll)) then
      this%incl_coll = include_coll
    else
      this%incl_coll = .false.
    end if

    this%sw = spectrum_sw(Tx, alpha, betaH, gx, vw, css)
    this%turb = spectrum_turb(Tx, alpha, betaH, gx, vw, epss, css)
    if (this%incl_coll) then
      this%coll = spectrum_coll(Tx, alpha, betaH, gx)
    end if
    this%sens = sens_curves(P, A)

    allocate(this%f(fnum))
    allocate(this%S(fnum, num_chunks))
    allocate(this%N(fnum, num_chunks))
    allocate(this%D(fnum, num_chunks))
    allocate(this%Dbar(fnum))
    allocate(this%variance(fnum))

    this%f = construct_frequencies()
    this%num = fnum

    if (.not. are_rands_initializaed) then
      call initialize_rands(mode='twister')
    end if

    call this%construct_signal()

  end function create_noisy_data

  function construct_frequencies() result(f)

    real(wp) :: f(fnum)

    integer :: i

    f(1) = fmin
    do i = 2, fnum
      f(i) = f(i - 1) + fdel
    end do

  end function construct_frequencies

  function OmegaGW_hsq(this, f0) result(y)

      class(noisy_data), intent(inout) :: this
      real(wp), intent(in) :: f0
      real(wp) :: y

      y = this%sw%OmegaGW_hsq(f0) + this%turb%OmegaGW_hsq(f0)
      if (this%incl_coll) then
        y = y + this%coll%OmegaGW_hsq(f0)
      end if

  end function OmegaGW_hsq

  function OmegaS_hsq(this, f0) result(y)

      class(noisy_data), intent(inout) :: this
      real(wp), intent(in) :: f0
      real(wp) :: y

      y = this%sens%OmegaS_hsq(f0)

  end function OmegaS_hsq

  function OmegaExt_hsq(this, f0) result(y)

      class(noisy_data), intent(in) :: this
      real(wp), intent(in) :: f0
      real(wp) :: y
      real(wp) :: OmegaExt

      OmegaExt = 10.0_wp**this%logOmegaExt

      y = OmegaExt * (f0 / 1.0e-3_wp)**(2.0_wp / 3.0_wp)

  end function OmegaExt_hsq
  function OmegaGal_hsq(this, f0) result(y)

      class(noisy_data), intent(in) :: this
      real(wp), intent(in) :: f0
      real(wp) :: y

      real(wp), parameter :: Tobs = 4.0_wp

      real(wp), parameter :: a1 = -0.15_wp
      real(wp), parameter :: b1 = -2.72_wp

      real(wp), parameter :: ak = -0.37_wp
      real(wp), parameter :: bk = -2.49_wp

      real(wp), parameter :: nu = 1.56_wp
      real(wp), parameter :: f2 = 6.7e-4_wp

      real(wp) :: f1
      real(wp) :: fknee
      real(wp) :: OmegaGal

      f1 = 10.0_wp**(a1 * log10(Tobs) + b1)
      fknee = 10.0_wp**(ak * log10(Tobs) + bk)

      OmegaGal = 10.0_wp**this%logOmegaGal

      y = 0.5_wp * f0**3                              &
	* (f0 / 1.0_wp)**(-7.0_wp / 3.0_wp)         &
	* (1.0_wp + tanh((fknee - f0) / f2))        &
	* exp(-(f0 / f1)**nu)                       &
	* OmegaGal

  end function OmegaGal_hsq

  subroutine construct_signal(this)

    class(noisy_data), intent(inout) :: this

    integer :: i
    integer :: j
    real(wp) :: f0
    real(wp) :: xGW
    real(wp) :: xGal
    real(wp) :: xExt
    real(wp) :: xS

    real(wp) :: SGW
    real(wp) :: SGal
    real(wp) :: SExt
    real(wp) :: SN

    !$omp parallel do  &
    !$omp IF(.NOT. omp_in_parallel())  &
    !$omp default(none)  &
    !$omp private(i, j, f0, xGW, xGal, xExt, xS, SGW, SGal, SExt, SN) &
    !$omp shared(this)
    do i = 1, this%num
      f0 = this%f(i)
      xGW = sqrt(this%OmegaGW_hsq(f0))
      xGal = sqrt(this%OmegaGal_hsq(f0))
      xExt = sqrt(this%OmegaExt_hsq(f0))
      xS = sqrt(this%OmegaS_hsq(f0))
      do j = 1, num_chunks
        SGW = abs(Gi(0.0_wp, xGW) + (0,1) * Gi(0.0_wp, xGW))**2 / 2.0_wp
        SGal = abs(Gi(0.0_wp, xGal) + (0,1) * Gi(0.0_wp, xGal))**2 / 2.0_wp
        SExt = abs(Gi(0.0_wp, xExt) + (0,1) * Gi(0.0_wp, xExt))**2 / 2.0_wp
        SN = abs(Gi(0.0_wp, xS) + (0,1) * Gi(0.0_wp, xS))**2 / 2.0_wp

        this%S(i,j) = SGW + SGal + SExt
        this%N(i,j) = SN

        this%D(i,j) = this%S(i,j) + this%N(i,j)
      end do
    end do
    !$omp end parallel do

    !$omp parallel do  &
    !$omp IF(.NOT. omp_in_parallel())  &
    !$omp default(none)  &
    !$omp private(i)  &
    !$omp shared(this)
    do i = 1, this%num
      this%Dbar(i) = sum(this%D(i,:)) / real(num_chunks, wp)
    end do
    !$omp end parallel do

    !$omp parallel do  &
    !$omp IF(.NOT. omp_in_parallel())  &
    !$omp default(none)  &
    !$omp private(i)  &
    !$omp shared(this)
    do i = 1, this%num
      this%variance(i) = sum((this%D(i,:) - this%Dbar(i)) ** 2) / real(num_chunks, wp)
    end do
    !$omp end parallel do

  end subroutine construct_signal

end module gwlisa__simulation_data
