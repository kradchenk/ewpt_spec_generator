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
  contains
    procedure, public :: OmegaGW_hsq
    procedure, public :: OmegaS_hsq
    procedure, private :: construct_signal
  end type noisy_data

  interface noisy_data
    procedure create_noisy_data
  end interface noisy_data

contains

  function create_noisy_data(  &
    P, A,  &
    Tx, alpha, betaH, gx, vw, ep, cs, include_coll) result(this)

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
    type(noisy_data) :: this

    real(wp) :: css
    real(wp) :: epss


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

  subroutine construct_signal(this)

    class(noisy_data), intent(inout) :: this

    integer :: i
    integer :: j
    real(wp) :: f0
    real(wp) :: xGW
    real(wp) :: xS

    !$omp parallel do  &
    !$omp IF(.NOT. omp_in_parallel())  &
    !$omp default(none)  &
    !$omp private(i, j, xGW, xS, f0)  &
    !$omp shared(this)
    do i = 1, this%num
      f0 = this%f(i)
      xGW = sqrt(this%OmegaGW_hsq(f0))
      xS = sqrt(this%OmegaS_hsq(f0))
      do j = 1, num_chunks
        this%S(i, j) = abs(Gi(0.0e0_wp, xGW) +  &
          (0, 1) * Gi(0.0e0_wp, xGW)) ** 2 / 2.0e0_wp
        this%N(i, j) = abs(Gi(0.0e0_wp, xS) +  &
          (0, 1) * Gi(0.0e0_wp, xS)) ** 2 / 2.0e0_wp
        this%D(i, j) = this%S(i, j) + this%N(i, j)
      end do
    end do
    !$omp end parallel do

    !$omp parallel do  &
    !$omp IF(.NOT. omp_in_parallel())  &
    !$omp default(none)  &
    !$omp private(i)  &
    !$omp shared(this)
    do i = 1, this%num
      this%Dbar(i) = sum(this%D(i,:)) / num_chunks
    end do
    !$omp end parallel do

    !$omp parallel do  &
    !$omp IF(.NOT. omp_in_parallel())  &
    !$omp default(none)  &
    !$omp private(i)  &
    !$omp shared(this)
    do i = 1, this%num
      this%variance(i) = sum((this%D(i,:) - this%Dbar(i)) ** 2) / num_chunks
    end do
    !$omp end parallel do

  end subroutine construct_signal

end module gwlisa__simulation_data
