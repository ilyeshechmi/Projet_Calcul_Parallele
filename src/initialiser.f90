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
          i_CI = 1  ! CI : nulle partout
          i_CL = 2  ! CL sinusoidale
        case default
          print *, "Erreur : Cas de test inconnu"
          stop
      end select
  end subroutine initialiser


  function C_init(i_CL, x) result(u0)
    implicit none
    integer,  intent(in) :: i_CL
    real(pr), intent(in) :: x
    real(pr)             :: u0

    select case(i_CL)

    case(1)
        ! --- Nulle ---
      u0= 0.0_pr

    case(2)
      ! --- Gaussienne ---
      u0 = exp( - ((x - 0.5_pr)/0.1_pr)**2 )

    case(3)
      ! --- Sinus ---
      u0 = sin(2.0_pr*pi*x)

    case default
      print *, "Erreur : Condition initiale inconnue"
      stop
    end select

  end function C_init


  function C_limite(i_CL, t) result(uL)
    implicit none
    integer,  intent(in) :: i_CL
    real(pr), intent(in) :: t
    real(pr)             :: uL

    select case(i_CL)

    case(1)
      ! --- Valeur constante ---
      uL = 1.0_pr

    case(2)
      ! --- Sinusoidale ---
      uL = sin(2.0_pr*pi*t)

    case default
      print *, "Erreur : Condition limite inconnue"
      stop  
    end select

  end function C_limite



end module 