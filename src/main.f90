program advection_rusanov
    
  use precision_mod
  use fonctions_mod
  use donnees_mod
  use schema_rusanov_mod
  use initialiser_mod
  use Euler_mod
  implicit none

  !---------------- Variables & paramètres ----------------
  type(Parametres) :: params
  integer  :: nx, nsteps, n, save_every
  integer  :: i, cl_periodique, cas_test
  integer  :: nargs, i_CI, i_CL, ios
  integer  :: nvar                    ! nb de variables : 1 (advection) ou 3 (Euler)

  real(pr) :: L, a, dt, Tfinal, dx, t, cfl
  real(pr) :: errL2, errLinf

  real(pr), allocatable :: x(:), u(:,:), u_ex(:,:), xeff(:)
  character(len=256) :: fichier_param

  !=======================================================
  ! 1) Lecture du fichier de paramètres
  !=======================================================

  nargs = command_argument_count()
  if (nargs < 1) then
    write(*,*) "Erreur: Veuillez fournir le fichier de paramètres en argument."
    stop
  end if
  call get_command_argument(1, fichier_param)

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

  !=======================================================
  ! 2) Initialisation des indicateurs de CI/CL
  !=======================================================

  call initialiser(cas_test, i_CL, i_CI, cl_periodique)

  ! Détermination du nombre de variables
  select case (cas_test)
  case (1, 2, 3)
     nvar  = 1       ! advection scalaire
  case (4)
     nvar  = 3       ! Euler 1D : (rho, rho*u, E)
  case default
     write(*,*) "Erreur: cas_test inconnu = ", cas_test
     stop
  end select

  !=======================================================
  ! 3) Grille spatiale et allocation de u
  !=======================================================

  allocate(x(nx), xeff(nx))
  x = [ ((i-0.5_pr)*dx, i=1,nx) ]

  
  allocate(u(nx, nvar))

  !=======================================================
  ! 4) Condition initiale
  !=======================================================

  if (cas_test /= 4) then
    !----------------------------------------------
    ! Cas 1–3 : advection scalaire
    ! u(i,1) = C_init(i_CI, x(i)write(*,*) "Pas d’erreur analytique )
    !----------------------------------------------
    do i = 1, nx
       u(i,1) = C_init(i_CI, x(i))
    end do

  else
    !----------------------------------------------
    ! Cas 4 : Euler – tube de choc de Sod
    ! u(:,1) = rho, u(:,2) = rho*u, u(:,3) = E
    !----------------------------------------------
    call init_euler_sod(x, u, L)
  end if

  !=======================================================
  ! 5) Sortie initiale
  !=======================================================

  t = 0.0_pr
  call ecrire(trim(params%outfile), t, x, u)


  !=======================================================
  ! 6) Boucle en temps
  !=======================================================

  nsteps = ceiling( Tfinal / dt )
  n      = 0

  do while (t < Tfinal)

    if (t + dt > Tfinal) then
      dt = Tfinal - t
    end if

    if (cas_test /= 4) then
       call avancer_Rusanov(u(:,1), a, dx, dt, cl_periodique, t, i_CL)
    else
       call avancer_rusanov_euler(u, dx, dt, gamma)
    end if

    t = t + dt
    n = n + 1

    if (mod(n, save_every) == 0 .or. n == nsteps) then  
          call ecrire(trim(params%outfile), t, x, u)
       
    end if

  end do

  write(*,*) "Terminé. Fichiers générés: "//trim(params%outfile)//"_t=<t>.dat"
  write(*,*) "============================================================"

  !=======================================================
  ! 7) Erreur (advection uniquement)
  !=======================================================

  if (cas_test /= 4) then
     allocate(u_ex(nx,1))
     u_ex(:,1) = U_exa(cas_test, x, t, a, L)
   
     errL2   = erreur_L2(u_ex, u, dx)
     errLinf = erreur_Linf(u_ex, u)

     write(*,'(A25, ES10.3)')  "  Erreur L2:",   errL2
     write(*,'(A25, ES10.3)')  "  Erreur Linf:", errLinf

     open(unit=16, file="erreurs.dat", status="unknown", position="append", iostat=ios)
     write(16,'(3ES16.8)') dx, errL2, errLinf
     close(16)

     deallocate(u_ex)
  else
   allocate(u_ex(nx,3))
   call sod_solution(u_ex, x, nx, t)
   errL2   = erreur_L2(u_ex, u, dx)
   errLinf = erreur_Linf(u_ex, u)

   do i = 1,3
      print *, "uex(:,",i," )=", u_ex(:,i)
      print *, "-------------------------------------"
      print *, "u(:,",i," ) =", u(:,i)
      print *, "-------------------------------------"   
      print *, "uex(:,",i," ) - u(:,",i," ) =", u_ex(:,i) - u(:,i)
      print *, "=====================================" 
   end do
     write(*,'(A25, ES10.3)')  "  Erreur L2:",   errL2
     write(*,'(A25, ES10.3)')  "  Erreur Linf:", errLinf

     open(unit=16, file="erreurs.dat", status="unknown", position="append", iostat=ios)
     write(16,'(3ES16.8)') dx, errL2, errLinf
     close(16)

     deallocate(u_ex)
  end if

  !=======================================================
  ! 8) Libération mémoire
  !=======================================================

  if (allocated(u))    deallocate(u)
  if (allocated(x))    deallocate(x)
  if (allocated(xeff)) deallocate(xeff)

  write(*,*) "========================FIN du programme====================="

end program advection_rusanov
