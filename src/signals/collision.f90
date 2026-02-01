module gwlisa__signals_collision

  use gwlisa__util_kinds, only : wp
  use gwlisa__util_constants, only : h

  implicit none

  private

  real(wp), private, parameter :: Astr = 0.05e0_wp
  real(wp), private, parameter :: a1 = 1.2e0_wp
  real(wp), private, parameter :: n1 = 2.4e0_wp
  real(wp), private, parameter :: n2 = -2.4e0_wp

  ! Source: [2403.03723]
  type, public :: spectrum_coll
    real(wp) :: Tx
    real(wp) :: alpha
    real(wp) :: betaH
    real(wp) :: gx
    real(wp) :: Hx0
    real(wp) :: fp
    real(wp) :: fb
    real(wp) :: Ktilde
    real(wp) :: Ast
    real(wp) :: FGW0_hsq
    real(wp) :: FGW0
    real(wp) :: Omegap_hsq
    real(wp) :: Omegap
    real(wp) :: Omegab_hsq
    real(wp) :: Omegab
  contains
    procedure, private :: calc_Hx0
    procedure, private :: calc_fp
    procedure, private :: calc_Ktilde
    procedure, private :: calc_FGW0_hsq
    procedure, private :: calc_Omegap_hsq
    procedure, private :: C
    procedure, public :: omegaGW
    procedure, public :: omegaGW_hsq
  end type spectrum_coll

  interface spectrum_coll
    procedure create_spectrum_coll
  end interface spectrum_coll

contains

  function create_spectrum_coll(  &
    Tx, alpha, betaH, gx) result(this)

    real(wp), intent(in) :: Tx
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: betaH
    real(wp), intent(in) :: gx
    type(spectrum_coll) :: this

    this%Tx = Tx
    this%alpha = alpha
    this%betaH = betaH
    this%gx = gx
    this%Hx0 = this%calc_Hx0()
    this%fp = this%calc_fp()
    this%fb = this%fp
    this%Ktilde = this%calc_Ktilde()
    this%FGW0_hsq = this%calc_FGW0_hsq()
    this%FGW0 = this%FGW0_hsq / (h ** 2)
    this%Omegap_hsq = this%calc_Omegap_hsq()
    this%Omegap = this%Omegap_hsq / (h ** 2)
    this%Omegab_hsq = this%Omegap_hsq
    this%Omegab = this%Omegap

  end function create_spectrum_coll

  function calc_Hx0(this) result(Hx0)

    class(spectrum_coll), intent(in) :: this
    real(wp) :: Hx0

    Hx0 = 1.65e-5 * (this%gx / 100.0e0_wp) **  &
      (1.0e0_wp / 6.0e0_wp) *  &
      (this%Tx / 100.0e0_wp)

  end function calc_Hx0

  function calc_fp(this) result(fp)

    class(spectrum_coll), intent(in) :: this
    real(wp) :: fp

    fp = 0.11e0_wp * this%Hx0 * this%betaH

  end function calc_fp

  function calc_Ktilde(this) result(Ktilde)

    class(spectrum_coll), intent(in) :: this
    real(wp) :: Ktilde

    Ktilde = this%alpha / (1.0e0_wp + this%alpha)

  end function calc_Ktilde

  function calc_FGW0_hsq(this) result(FGW0_hsq)

    class(spectrum_coll), intent(in) :: this
    real(wp) :: FGW0_hsq

    FGW0_hsq = 1.64e-5_wp * (100.0e0_wp / this%gx) **  &
      (1.0e0_wp / 3.0e0_wp)

  end function calc_FGW0_hsq

  function calc_Omegap_hsq(this) result(Omegap_hsq)

    class(spectrum_coll), intent(in) :: this
    real(wp) :: Omegap_hsq

    Omegap_hsq = this%FGW0_hsq * Astr * this%Ktilde ** 2 /  &
      (this%betaH ** 2)

  end function calc_Omegap_hsq

  function C(this, f0) result(y)

    class(spectrum_coll), intent(in) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y

    real(wp) :: f

    f = f0 / this%fb
    y = f ** n1 * (0.5e0_wp + 0.5e0_wp * f ** a1) ** ((n2 - n1) / a1)

  end function C

  function OmegaGW(this, f0) result(y)

    class(spectrum_coll), intent(in) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y

    y = this%Omegab * this%C(f0)

  end function OmegaGW

  function OmegaGW_hsq(this, f0) result(y)

    class(spectrum_coll), intent(in) :: this
    real(wp), intent(in) :: f0
    real(wp) :: y

    y = this%Omegab_hsq * this%C(f0)

  end function OmegaGW_hsq

end module
