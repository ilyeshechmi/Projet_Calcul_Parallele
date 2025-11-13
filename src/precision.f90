module precision
    
  implicit none
  
  integer, parameter :: pr  = kind(1.0d0)

  real(pr), parameter :: tol = 1.0e-3_pr

  real(pr), parameter :: pi = acos(-1.0_pr)

end module precision
