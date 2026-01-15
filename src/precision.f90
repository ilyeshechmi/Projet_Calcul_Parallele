module precision_mod
  implicit none

  integer, parameter :: pr  = kind(1.0d0)

  real(pr), parameter :: tol = 1.0e-12_pr
  real(pr), parameter :: pi  = acos(-1.0_pr)

  real(pr), parameter :: gamma = 1.4_pr   


  real(pr), parameter :: theta   = 0.1_pr      
  real(pr), parameter :: F  = 1.0_pr      
  real(pr), parameter :: lambda  = 1.0_pr      
  real(pr), parameter :: epsilon = 0.05_pr
  real(pr), parameter :: kappa   = 1.0_pr
  real(pr), parameter :: alpha   = 1.0e-2_pr
  real(pr), parameter :: beta    = 1.0e-3_pr
  real(pr), parameter :: Re      = 20.0_pr

end module precision_mod
