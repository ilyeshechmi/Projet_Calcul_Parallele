program advection_rusanov
    
  use precision_mod
  use fonctions_mod
  use donnees_mod
  use schema_rusanov_mod
  implicit none

  type(Parametres) :: params
  integer  :: nx, nsteps, n, save_every, i ,cl_periodique ,cas_test 
  real(pr) :: L, a, dt, Tfinal, dx, t, cfl, uG
  real(pr), allocatable :: x(:), u(:,:),u_ex(:),xeff(:)
  real(pr) ::  errL2, errLinf
  character(len=256) :: fichier_param
  integer :: nargs, i_CI, i_CL , ios

  ! Récupérer le nom du fichier en ligne de commande
  nargs = command_argument_count()
  if (nargs < 1) then
    write(*,*) "Erreur: Veuillez fournir le fichier de paramètres en argument."
    stop
  end if
  call get_command_argument(1, fichier_param)

  ! --- Lecture paramètres ---
  call lire_parametres(fichier_param, params)
  cas_test   = params%cas_test
  nx         = params%nx
  L          = params%L
  a          = params%a
  Tfinal     = params%T_final
  save_every = params%save_every
  CFL        = params%CFL
  dx         = params%dx
  dt         = params%dt
  call initialiser(cas_test ,i_CL , i_CI, CL_periodique )

  ! --- Grille uniforme (centres de mailles) ---
  allocate(x(nx), u(nx,1),xeff(nx))
  x  = [ ( (i-0.5_pr)*dx, i=1,nx ) ]

  ! --- Condition initiale ---
  do i = 1, nx
    u(i,1) = C_init(i_CI ,x(i))
  end do
  

  ! --- Sortie initiale ---
  t = 0.0_pr
  call ecrire(trim(params%outfile), t, x, u(:,1))

  ! --- Boucle en temps ---

  nsteps = ceiling( Tfinal / dt )
  n=0
  do while (t< Tfinal)

    if (t + dt > Tfinal) then ! Cas ou t_final n'est pas un multiple de dt
      dt = Tfinal - t
    end if
    call avancer_Rusanov(u(:,1), a, dx, dt, cl_periodique, t, i_CL)
    t = t + dt
    n=n+1
    ! Écriure des résultats
    if (mod(n, save_every) == 0 .or. n == nsteps) then
      call ecrire(trim(params%outfile), t, x, u(:,1)) 
    end if
  end do

  write(*,*) "Terminé. Fichiers générés: "//trim(params%outfile)//"_t=<t>.dat"
  write(*,*) "============================================================"
  
  
  
  !------- Calcul de l'erreur -----------
  allocate(u_ex(nx))
  u_ex=U_exa(cas_test,x,t,a,L)
  
  ! Calcul des erreurs
  errL2   = erreur_L2(u_ex, u(:,1), dx)
  errLinf = erreur_Linf(u_ex, u(:,1))
  ! Affichage des erreurs
  write(*,'(A25, ES10.3)')  "  Erreur L2:",           errL2
  write(*,'(A25, ES10.3)')  "  Erreur Linf:",         errLinf
  ! Enregistrement des erreurs dans un fichier
  open(unit=16, file="erreurs.dat", status="unknown", position="append", iostat=ios)

  if (ios /= 0) then
      write(*,*) "Erreur : open(unit=16) a échoué"
      stop
  end if
  write(16,*)'zebi'
  write(16,'(3ES16.8)') dx, errL2, errLinf
  close(16)
  write(*,*) "Enregistrement de l'erreur dans erreurs.dat"
  
  deallocate(x, u,u_ex,xeff)
  
  
  write(*,*) "========================FIN du programme====================="


end program advection_rusanov
