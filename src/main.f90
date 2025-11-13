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
  integer :: nargs

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
  

  ! --- Grille uniforme (centres de mailles) ---
  allocate(x(nx), u(nx,1),xeff(nx))
  x  = [ ( (i-0.5_pr)*dx, i=1,nx ) ]

  ! --- Condition initiale ---
   u(:,1) = sin( 2.0_pr * pi * x / L )
  !u(:,1) = gaussienne(x)
  
  ! --- Condition à Gauche ---
  uG = 1._pr


  ! --- Sortie initiale ---
  t = 0.0_pr
  call ecrire(trim(params%outfile), t, x, u(:,1))

  ! --- Boucle en temps ---
  dt= CFL*dx/abs(a)

  nsteps = ceiling( Tfinal / dt )
  do n = 1, nsteps
    call avancer_Rusanov(u(:,1), a, dx, dt,uG,cl_periodique)
    t = t + dt
    ! Écriure des résultats
    if (mod(n, save_every) == 0 .or. n == nsteps) then
      call ecrire(trim(params%outfile), t, x, u(:,1)) 
    end if
  end do

  write(*,*) "Terminé. Fichiers générés: "//trim(params%outfile)//"_t=<t>.dat"
  write(*,*) "============================================================"

  ! Calcul de l'erreur
  allocate(u_ex(nx))

! ! Calcul des erreurs
! errL2   = erreur_L2(u_ex, u(:,1), dx)
! errLinf = erreur_Linf(u_ex, u(:,1))

! ! Ouverture du fichier en mode ajout
! open(unit=20, file=trim(params%outfile)//'_erreurs.dat', status='unknown', position='append', action='write')

! ! Écriture d'une ligne : dx, errL2, errLinf
! write(20, '(F12.6, 2X, F12.6, 2X, F12.6)') dx, errL2, errLinf

! close(20)

! write(*,*) 'Erreurs calculées et sauvegardées dans ', trim(params%outfile)//'_erreurs.dat'

  ! --- Sauvegarde des solutions finale et exacte ---
! open(unit=30, file=trim(params%outfile)//'_compare.dat', status='replace')
! write(30, '(A)') '# x   u_num   u_exact'
! do i = 1, nx
!     write(30, '(3E20.10)') x(i), u(i,1),! u_ex(i)
! end do
! close(30)


  deallocate(x, u,u_ex,xeff)




end program advection_rusanov
