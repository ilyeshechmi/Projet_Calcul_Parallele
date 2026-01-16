module precision_mod
  implicit none

  integer, parameter :: pr  = kind(1.0d0)

  real(pr), parameter :: tol = 1.0e-12_pr
  real(pr), parameter :: pi  = acos(-1.0_pr)

  real(pr), parameter :: gamma = 1.4_pr 
  real(pr), parameter :: rho0  = 0.2_pr  

  real(pr), parameter :: theta   = 6.4_pr * pi / 180._pr
  real(pr), parameter :: F       = 0.8476_pr
  real(pr), parameter :: lambda  = 3.0_pr
  real(pr), parameter :: epsilon = 6.07e-3_pr
  real(pr), parameter :: kappa   = 3.866_pr
  real(pr), parameter :: alpha   = 6.07e-3_pr
  real(pr), parameter :: beta    = 2.243e-7_pr
  real(pr), parameter :: Re      = 19.33_pr
  real(pr), parameter :: freq    = 1.5

end module precision_mod
