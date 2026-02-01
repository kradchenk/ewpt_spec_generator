module gwlisa__signals_turbulance

  use gwlisa__util_kinds, only : wp
  use gwlisa__util_constants, only : pi
  use gwlisa__util_constants, only : cs_default
  use gwlisa__util_constants, only : ZP
  use gwlisa__util_constants, only : adiab_ratio
  use gwlisa__util_constants, only : h
  use gwlisa__signals_util, only : calc_Hx
  use gwlisa__signals_util, only : calc_kappa
  use gwlisa__signals_util, only : calc_HxRx

  implicit none

  private

  real(wp), parameter :: Amhd = 0.00437_wp
  real(wp), parameter :: a1  = 4.0_wp
  real(wp), parameter :: a2  = 2.15_wp
  real(wp), parameter :: n1  = 3.0_wp
  real(wp), parameter :: n2  = 1.0_wp
  real(wp), parameter :: n3  = -8.0/3.0_wp
  real(wp), parameter :: N  = 2.0_wp
  
  ! Source [2403.03723]
  type, public :: spectrum_turb
    real(wp) :: Tx
    real(wp) :: alpha
    real(wp) :: betaH
    real(wp) :: gx
    real(wp) :: vw
    real(wp) :: eps
    real(wp) :: cs
    real(wp) :: Hx0
    real(wp) :: HxRx
    real(wp) :: kappa
    real(wp) :: K
    real(wp) :: FGW0_hsq
    real(wp) :: FGW0
    real(wp) :: Omegas
    real(wp) :: f1
    real(wp) :: f2
    real(wp) :: Omega2_hsq
    real(wp) :: Omega2

  contains
    procedure, private :: calc_Hx0
    procedure, private :: calc_K
    procedure, private :: calc_FGW0_hsq
    procedure, private :: calc_f1
    procedure, private :: calc_f2
    procedure, private :: S2
    procedure, private :: calc_Omegas
    procedure, private :: calc_Omega2_hsq
    procedure, public :: omegaGW
    procedure, public :: omegaGW_hsq
  end type spectrum_turb

  interface spectrum_turb
    procedure create_spectrum_turb
  end interface spectrum_turb

contains

function create_spectrum_turb(  &
  Tx, alpha, betaH, gx, vw, eps, cs) result(this)

  real(wp), intent(in) :: Tx
  real(wp), intent(in) :: alpha
  real(wp), intent(in) :: betaH
  real(wp), intent(in) :: gx
  real(wp), intent(in) :: vw
  real(wp), intent(in), optional :: eps
  real(wp), intent(in), optional :: cs
  type(spectrum_turb) :: this

  this%Tx    = Tx
  this%alpha = alpha
  this%betaH = betaH
  this%gx    = gx
  this%vw    = vw

  this%Hx0 = this%calc_Hx0()
  this%kappa = calc_kappa(alpha, vw, this%cs)
  this%K = this%calc_K()
  this%HxRx = calc_HxRx(betaH, vw, this%cs)
  this%FGW0_hsq = this%calc_FGW0_hsq()
  this%FGW0 = this%FGW0_hsq / (h**2)
  this%Omegas = this%calc_Omegas()
  this%f1 = this%calc_f1()
  this%f2 = this%calc_f2()
  this%Omega2_hsq = this%calc_Omega2_hsq()
  this%Omega2 = this%Omega2_hsq / (h**2)

  end function create_spectrum_turb
  
  function calc_Hx0(this) result(Hx0)
    class(spectrum_turb), intent(in) :: this
    real(wp) :: Hx0
    Hx0 = 1.65e-5_wp * (this%gx / 100.0_wp)**(1.0_wp/6.0_wp) * (this%Tx / 100.0_wp)
  end function calc_Hx0

  function calc_K(this) result(K)
    class(spectrum_turb), intent(in) :: this
    real(wp) :: K
    K = 0.6_wp * this%kappa * this%alpha / (1.0_wp + this%alpha)
  end function calc_K

  function calc_FGW0_hsq(this) result(FGW0_hsq)
    class(spectrum_turb), intent(in) :: this
    real(wp) :: FGW0_hsq
    FGW0_hsq = 1.64e-5_wp * (100.0_wp / this%gx)**(1.0_wp/3.0_wp)
  end function calc_FGW0_hsq
 
   function calc_Omegas(this) result(Omegas)
    class(spectrum_turb), intent(in) :: this
    real(wp) :: Omegas
    Omegas = this%eps * this%K
  end function calc_Omegas
  
  function calc_f1(this) result(f1)
    class(spectrum_turb), intent(in) :: this
    real(wp) :: f1
    f1 = sqrt(3.0_wp*this%Omegas) / (2.0_wp*N) * (this%Hx0 / this%HxRx)
  end function calc_f1

  function calc_f2(this) result(f2)
    class(spectrum_turb), intent(in) :: this
    real(wp) :: f2
    f2 = 2.2_wp * (this%Hx0 / this%HxRx)
  end function calc_f2

  function calc_Omega2_hsq(this) result(Omega2_hsq)
    class(spectrum_turb), intent(in) :: this
    real(wp) :: Omega2_hsq
    ! Eq 2.23
    Omega2_hsq = this%FGW0_hsq * Amhd * (this%Omegas**2) * (this%HxRx**2)
  end function calc_Omega2_hsq

  function S2(this, f0) result(y)
    class(spectrum_turb), intent(in) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y
    real(wp) :: ff1, ff2, f1f2

    ! Eq 2.24
    ff1  = f0 / this%f1
    ff2  = f0 / this%f2
    f1f2 = this%f1 / this%f2

    y = (2.0_wp)**(11.0_wp/(3.0_wp*a2)) * f1f2 * (ff1**3.0_wp) * &
        (1.0_wp + ff1**a1)**(-2.0_wp/a1) * &
        (1.0_wp + ff2**a2)**(-11.0_wp/(3.0_wp*a2))
  end function S2

  function omegaGW(this, f0) result(y)
    class(spectrum_turb), intent(in) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y
    y = this%Omega2 * this%S2(f0)
  end function omegaGW

  function omegaGW_hsq(this, f0) result(y)
    class(spectrum_turb), intent(in) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y
    y = this%Omega2_hsq * this%S2(f0)
  end function omegaGW_hsq

end module gwlisa__signals_turbulance