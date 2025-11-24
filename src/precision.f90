module precision_mod
    
  implicit none
  
  integer, parameter :: pr  = kind(1.0d0)

  real(pr), parameter :: tol = 1.0e-3_pr

  real(pr), parameter :: pi = acos(-1.0_pr)

  real(pr), parameter :: gamma = 1.4_pr
  
end module 
