module initialiser_mod 
  use precision_mod
  implicit none

contains
  subroutine initialiser(cas_test ,i_CL , i_CI, CL_periodique )
    integer, intent(in) :: cas_test
    integer, intent(out) ::  i_CI, CL_periodique ,i_CL

    ! Initialisation 
    i_CI=0
    i_CL=0
    CL_periodique=0

    select case (cas_test)
      case (1) ! cas test 1 : gaussienne periodique
        i_CI = 2  ! CI : gaussienne
        CL_periodique = 1  ! CL periodique

      case (2) !cas test 2 : créneau 
        i_CI = 1  ! CI : nulle partout
        i_CL = 1  ! CL : =1 

      case (3) 
        i_CI = 3  ! CI : nulle partout
        i_CL = 2  ! CL sinusoidale
      ! ---------- Équations d'euler 
      case (4,5)  ! cas test 4 : équations d'Euler (tube de choc de Sod)
        i_CI = 0
        i_CL = 0
        CL_periodique = 0
  
      
      case default
        print *, "Erreur : Cas de test inconnu"
        stop
    end select
  end subroutine initialiser

  !-------------------------------------------------------
  ! CI du tube de choc de Sod
  ! U est (nx,3) : (rho, rho*u, E)
  !-------------------------------------------------------
  subroutine init_euler(x, U, L)
    use precision_mod
    implicit none
    real(pr), intent(in)  :: x(:), L
    real(pr), intent(out) :: U(:,:)         ! (nx,3)
    integer :: i, nx
    real(pr) :: rho, V, p, E

    nx = size(x)

    ! États initiales
    V = 1.0_pr; p = 2.0_pr

    do i = 1, nx
      rho = 1._pr +rho0*sin(2*pi*x(i))
      E      = p/(gamma - 1.0_pr) + 0.5*rho*V*V
      U(i,1) = rho
      U(i,2) = rho*V
      U(i,3) = E
    end do
  end subroutine init_euler 

  subroutine init_saint_venant(u)
  use precision_mod
  implicit none

  real(pr), intent(out) :: u(:,:)
  integer :: i, nx

  nx = size(u,1)

  do i = 1, nx
     u(i,1) = 1._pr   ! h
     u(i,2) = 1._pr   ! hU = h*U
     u(i,3) = 1._pr   ! hη
     u(i,4) = 0._pr   ! hw
     u(i,5) = 0._pr   ! p
  end do

end subroutine init_saint_venant





end module 