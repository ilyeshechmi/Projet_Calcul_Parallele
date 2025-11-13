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


end module 