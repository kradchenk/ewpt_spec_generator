module gwlisa__sensitivity_curves

  use gwlisa__util_kinds, only : wp
  use gwlisa__util_constants, only : pi
  use gwlisa__util_constants, only : c
  use gwlisa__util_constants, only : H0
  use gwlisa__util_constants, only : h

  implicit none

  private

  real(wp), public, parameter :: L = 2.5e9_wp

  ! Source: [1906.09244]
  type, public :: sens_curves
    real(wp) :: P
    real(wp) :: A
  contains
    procedure, public :: psd_oms  ! (2.1)
    procedure, public :: psd_acc  ! (2.1)
    procedure, public :: psd_nX   ! (2.2)
    procedure, public :: response ! (2.10)
    procedure, public :: Sn       ! (2.3)
    procedure, public :: OmegaS   ! (2.14)
    procedure, public :: OmegaS_hsq
  end type sens_curves

  interface sens_curves
    procedure create_sens_curves
  end interface sens_curves

contains

  function create_sens_curves(P, A) result(this)

    real(wp), intent(in) :: P
    real(wp), intent(in) :: A
    type(sens_curves) :: this

    this%A = A
    this%P = P

  end function create_sens_curves

  function psd_oms(this, f) result(y)

    class(sens_curves), intent(in) :: this
    real(wp), intent(in) :: f
    real(wp) :: y

    y = this%P ** 2 * 1.0e-12_wp ** 2 *  &
      (1.0e0_wp + (2.0e-6_wp / f) ** 4) *  &
      (2.0e0_wp * pi * f / c)**2

  end function psd_oms

  function psd_acc(this, f) result(y)

    class(sens_curves), intent(in) :: this
    real(wp), intent(in) :: f
    real(wp) :: y

    y = this%A ** 2 * 1.0e-15_wp ** 2 *  &
      (1.0e0_wp + (0.4e-3_wp / f) ** 2) *  &
      (1.0e0_wp + (f / 8.0e-3_wp) ** 4) *  &
      (1.0e0_wp / (2.0e0_wp * pi * f)) ** 4 *  &
      (2.0e0_wp * pi * f / c) ** 2

  end function psd_acc

  function psd_nX(this, f) result(y)

    class(sens_curves), intent(in) :: this
    real(wp), intent(in) :: f
    real(wp) :: y

    real(wp) :: poms
    real(wp) :: pacc

    real(wp) :: ff

    ff = 2.0e0_wp * pi * f * L / c

    poms = this%psd_oms(f)
    pacc = this%psd_acc(f)

    y = 16.0e0_wp * sin(ff) ** 2 *  &
      (poms + (3.0e0_wp + cos(2.0e0_wp * ff)) * pacc)

  end function psd_nX

  function response(this, f) result(y)

    class(sens_curves), intent(in) :: this
    real(wp), intent(in) :: f
    real(wp) :: y

    real(wp) :: ff

    ff = 2.0e0_wp * pi * f * L / c

    y = 16.0e0_wp * sin(ff) ** 2 * (3.0e0_wp / 10.0e0_wp) *  &
      ff ** 2 / (1.0e0_wp + 0.6e0_wp * ff**2)

  end function response

  function Sn(this, f) result(y)

    class(sens_curves), intent(in) :: this
    real(wp), intent(in) :: f
    real(wp) :: y

    y = this%psd_nX(f) / this%response(f)

  end function Sn

  function OmegaS(this, f) result(y)

    class(sens_curves), intent(in) :: this
    real(wp), intent(in) :: f
    real(wp) :: y

    y = 4.0e0_wp * pi ** 2 * f ** 3 * this%Sn(f) /  &
      (3.0e0_wp * H0 ** 2)

  end function OmegaS

  function OmegaS_hsq(this, f) result(y)

    class(sens_curves), intent(in) :: this
    real(wp), intent(in) :: f
    real(wp) :: y

    y = h**2 * this%OmegaS(f)

  end function OmegaS_hsq

end module gwlisa__sensitivity_curves
