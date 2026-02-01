module gwlisa__signals_soundwave

  use gwlisa__util_kinds, only : wp
  use gwlisa__util_constants, only : cs_default
  use gwlisa__util_constants, only : adiab_ratio
  use gwlisa__util_constants, only : pi
  use gwlisa__util_constants, only : ZP
  use gwlisa__util_constants, only : h
  use gwlisa__signals_util, only : calc_kappa
  use gwlisa__signals_util, only : calc_Ubarf
  use gwlisa__signals_util, only : calc_HxRx
  use gwlisa__signals_util, only : calc_HxtauSH
  use gwlisa__signals_util, only : calc_HxtauSW
  use gwlisa__signals_util, only : calc_gxFac

  implicit none

  private

  real(wp), parameter :: Asw = 0.11_wp
  real(wp), parameter :: a1  = 2.0_wp
  real(wp), parameter :: a2  = 4.0_wp
  real(wp), parameter :: n1  = 3.0_wp
  real(wp), parameter :: n2  = 1.0_wp
  real(wp), parameter :: n3  = -3.0_wp


  ! Source: [2403.03723]
  type, public :: spectrum_sw
    real(wp) :: Tx
    real(wp) :: alpha
    real(wp) :: betaH
    real(wp) :: gx
    real(wp) :: vw
    real(wp) :: cs
    real(wp) :: Hx0
    real(wp) :: kappa
    real(wp) :: K
    real(wp) :: Ubarf
    real(wp) :: HxRx
    real(wp) :: HxtauSH
    real(wp) :: HxtauSW
    real(wp) :: f1
    real(wp) :: f2
    real(wp) :: FGW0_hsq
    real(wp) :: FGW0
    real(wp) :: Omegaint_hsq
    real(wp) :: Omegaint
    real(wp) :: Omega2
    real(wp) :: Omega2_hsq

  contains
    procedure, private :: calc_Hx0
    procedure, private :: calc_K
    procedure, private :: calc_FGW0_hsq
    procedure, private :: calc_f1
    procedure, private :: calc_f2
    procedure, private :: S2
    procedure, private :: calc_Omegaint_hsq
    procedure, private :: calc_Omega2
    procedure, public :: omegaGW
    procedure, public :: omegaGW_hsq
  end type spectrum_sw

  interface spectrum_sw
    procedure create_spectrum_sw
  end interface spectrum_sw

contains

  function create_spectrum_sw(  &
    Tx, alpha, betaH, gx, vw, cs) result(this)

    real(wp), intent(in) :: Tx
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: betaH
    real(wp), intent(in) :: gx
    real(wp), intent(in) :: vw
    real(wp), intent(in), optional :: cs
    type(spectrum_sw) :: this

    this%Tx = Tx
    this%alpha = alpha
    this%betaH = betaH
    this%gx = gx
    this%vw = vw
    if (present(cs)) then
      this%cs = cs
    else
      this%cs = cs_default
    end if
    this%Hx0 = this%calc_Hx0()
    this%HxRx = calc_HxRx(betaH, vw, this%cs)
    this%kappa = calc_kappa(alpha, vw, this%cs)
    this%K = this%calc_K()
    this%Ubarf = calc_Ubarf(this%kappa, alpha)
    this%HxtauSH = calc_HxtauSH(this%HxRx, this%Ubarf)
    this%HxtauSW = calc_HxtauSW(this%HxtauSH)
    this%FGW0_hsq = this%calc_FGW0_hsq()
    this%FGW0 = this%FGW0_hsq / (h ** 2)
    this%f1 = this%calc_f1()
    this%f2 = this%calc_f2()
    this%Omegaint_hsq = this%calc_Omegaint_hsq()
    this%Omegaint = this%Omegaint_hsq / (h ** 2)
    this%Omega2 = this%calc_Omega2()
    this%Omega2_hsq = this%Omega2 * (h ** 2) 
  end function create_spectrum_sw

  function calc_K(this) result(K)

    class(spectrum_sw), intent(in) :: this
    real(wp) :: K

    K = 0.6_wp * this%kappa * this%alpha / (1.0_wp + this%alpha)

  end function calc_K

  function calc_Hx0(this) result(Hx0)

    class(spectrum_sw), intent(in) :: this
    real(wp) :: Hx0

    Hx0 = 1.65e-5 * (this%gx / 100.0e0_wp) **  &
      (1.0e0_wp / 6.0e0_wp) *  &
      (this%Tx / 100.0e0_wp)

  end function calc_Hx0
  
  function calc_FGW0_hsq(this) result(FGW0_hsq)

    class(spectrum_sw), intent(in) :: this
    real(wp) :: FGW0_hsq

    FGW0_hsq = 1.64e-5_wp * (100.0e0_wp / this%gx) **  &
      (1.0e0_wp / 3.0e0_wp)

  end function calc_FGW0_hsq
  
  function calc_f1(this) result(f1)

    class(spectrum_sw), intent(inout) :: this
    real(wp) :: f1

    f1 = 0.2_wp * this%Hx0 / this%HxRx

  end function calc_f1
  
  function calc_f2(this) result(f2)

    class(spectrum_sw), intent(inout) :: this
    real(wp) :: f2
    real(wp) :: xi_shell, Deltaw, denom

    xi_shell = abs(this%vw - this%cs)
    denom    = max(this%vw, this%cs)
    Deltaw   = xi_shell / denom
    
    f2 =(0.5_wp * this%Hx0 / this%HxRx) / Deltaw

  end function calc_f2
  
  function calc_Omegaint_hsq(this) result(Omegaint_hsq)

    class(spectrum_sw), intent(inout) :: this
    real(wp) :: Omegaint_hsq

    Omegaint_hsq = this%FGW0_hsq * Asw * (this%K**2) *  &
                   this%HxtauSW * this%HxRx
  end function calc_Omegaint_hsq

    function calc_Omega2(this) result(Omega2)
    class(spectrum_sw), intent(inout) :: this
    real(wp) :: Omega2
    real(wp) :: r

    r  = this%f2 / this%f1

    Omega2 = this%Omegaint * (1.0_wp/pi) *  &
             (sqrt(2.0_wp) + (2.0_wp*r)/(1.0_wp + r**2))
  end function calc_Omega2

  function S2(this, f0) result(y)
    class(spectrum_sw), intent(inout) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y
    real(wp) :: ff1, ff2, f2f1
    real(wp) :: num, den

    ff1  = f0 / this%f1
    ff2  = f0 / this%f2
    f2f1 = this%f2 / this%f1

    num = (ff1**n1) * (1.0_wp + ff1**a1)**((n2 - n1)/a1) *  &
          (1.0_wp + ff2**a2)**((n3 - n2)/a2)

    den = (f2f1**n1) * (1.0_wp + f2f1**a1)**((n2 - n1)/a1) *  &
          (1.0_wp + 1.0_wp**a2)**((n3 - n2)/a2)

    y = num / den
  end function S2

  function OmegaGW(this, f0) result(y)
    class(spectrum_sw), intent(inout) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y

    y = this%Omega2 * this%S2(f0)
  end function OmegaGW

  function OmegaGW_hsq(this, f0) result(y)
    class(spectrum_sw), intent(inout) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y

    y = this%Omega2_hsq * this%S2(f0)
  end function OmegaGW_hsq

end module gwlisa__signals_soundwave