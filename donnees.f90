module donnees
  use precision
  implicit none

  type :: Parametres
    integer  :: nx         = 200
    real(pr) :: L          = 1.0_pr
    real(pr) :: a          = 1.0_pr
    real(pr) :: dt         = 2.0e-3_pr
    real(pr) :: T_final    = 0.5_pr
    integer  :: save_every = 50
    character(len=256) :: outfile = "results.dat"
  end type Parametres

contains

  subroutine lire_parametres(fname, params)
    character(len=*), intent(in)  :: fname
    type(Parametres),      intent(out)  :: params
    integer :: ios, u
    integer :: nx, save_every
    real(pr) :: L, a, dt, T_final
    character(len=256) :: outfile

    namelist /input/ nx, L, a, dt, T_final, outfile, save_every

    ! initialise avec les valeurs par défaut du type
    nx         = params%nx
    L          = params%L
    a          = params%a
    dt         = params%dt
    T_final    = params%T_final
    save_every = params%save_every
    outfile    = params%outfile

    open(newunit=u, file=fname, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(*,*) " Impossible d'ouvrir ", trim(fname), " — on garde les valeurs par défaut."
    else
      read(u, nml=input, iostat=ios)
      close(u)
      if (ios /= 0) then
        write(*,*) "ATTENTION: Lecture namelist échouée — on garde les valeurs par défaut."
      endif
    endif

    ! copie vers params
    params%nx         = nx
    params%L          = L
    params%a          = a
    params%dt         = dt
    params%T_final    = T_final
    params%save_every = save_every
    params%outfile    = outfile

  end subroutine lire_parametres


  subroutine ecrire(fichier_sortie, t, x, u)
  
    character(len=*), intent(in) :: fichier_sortie   
    real(pr),         intent(in) :: t           ! temps courant
    real(pr),         intent(in) :: x(:), u(:)  ! vecteurs (même taille)

    character(len=64)  :: t_str
    character(len=512) :: fname
    integer :: uo, ios, i

    ! met t en chaîne sans espaces à gauche
    write(t_str,'(F12.6)') t
    t_str = adjustl(t_str)

    ! nom de fichier: <base>_t=<t>.dat
    fname = trim(fichier_sortie)//'_t='//trim(t_str)//'.dat'

    open(newunit=uo, file=fname, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(*,*) "ERREUR: ouverture de ", trim(fname), " impossible."
        return
    end if

    ! Écriture
    write(uo,'(a,1x,a)') '!t=', trim(t_str)
    do i = 1, size(x)
        write(uo,'(2(ES20.10,1X))') x(i), u(i)
    end do
    close(uo)
    end subroutine ecrire


end module donnees
