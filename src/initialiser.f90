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
      case (4)  ! cas test 4 : équations d'Euler (tube de choc de Sod)
        i_CI = 0
        i_CL = 0
        CL_periodique = 0
      
      case default
        print *, "Erreur : Cas de test inconnu"
        stop
    end select
  end subroutine initialiser


  subroutine init_euler_sod(x, U, gamma, L)
    use precision_mod
    implicit none
    real(pr), intent(in)  :: x(:), gamma, L
    real(pr), intent(out) :: U(:,:)         ! (nx,3)
    integer :: i, nx
    real(pr) :: rhoL, uL, pL, rhoR, uR, p_R, E

    nx = size(x)

    ! États gauche et droite classiques du test de Sod
    rhoL = 1.0
    uL   = 0.0
    pL   = 1.0

    rhoR = 0.125
    uR   = 0.0
    p_R   = 0.1

    do i = 1, nx
      if (x(i) < 0.5*L) then
        ! gauche
        E      = pL/(gamma - 1.0) + 0.5*rhoL*uL*uL
        U(i,1) = rhoL
        U(i,2) = rhoL*uL
        U(i,3) = E
      else
        ! droite
        E      = pR/(gamma - 1.0) + 0.5*rhoR*uR*uR
        U(i,1) = rhoR
        U(i,2) = rhoR*uR
        U(i,3) = E
      end if
    end do
  end subroutine init_euler_sod




end module 