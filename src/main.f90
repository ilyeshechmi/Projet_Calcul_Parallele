program advection_rusanov
    
  use precision_mod
  use fonctions_mod
  use donnees_mod
  use schema_rusanov_mod
  use initialiser_mod
  use Euler_mod
  use exact_sod_mod
  implicit none

  !---------------- Variables & paramètres ----------------
  type(Parametres) :: params
  integer  :: nx, nsteps, n, save_every
  integer  :: i, cl_periodique, cas_test,i_schema
  integer  :: nargs, i_CI, i_CL, ios
  integer  :: nvar                    ! nb de variables : 1 (advection) ou 3 (Euler)

  real(pr) :: L, dt, Tfinal, dx, t, cfl,rho0
  real(pr) :: rho, V, p,E
  real(pr) :: errL2, errLinf,errL1

  real(pr), allocatable :: x(:), u(:,:), u_ex(:,:), xeff(:),a(:)
  character(len=256) :: fichier_param
  integer, parameter :: Nf_ref = 20480   ! maillage de référence FIXE
  integer :: nx_f, r_ref, k, m, i0
  real(pr) :: dx_f, dt_f, t_f
  real(pr), allocatable :: x_f(:), u_f(:,:), u_ref(:,:)


  !=======================================================
  ! 1) Lecture du fichier de paramètres
  !=======================================================
   allocate(a(1))
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
  a(:)          = params%a
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
  
  i_schema= 2! 0: rusanov , 1: muscl , 2 : muscl hancock
  
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
    rho0 =0.2_pr
    call init_euler(x, u, L,rho0)
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
    if ( cas_test ==4  ) then
      dt = dt_CFL_euler(u, dx, CFL)
    end if
    
    if (t + dt > Tfinal) then
      dt = Tfinal - t
    end if
    
    call avancer_Rusanov(u, x,a, dx, dt, cl_periodique, t, i_CL,i_schema,cas_test)

    t = t + dt
    n = n + 1

    if (mod(n, save_every) == 0 .or. n == nsteps) then  
      call ecrire(trim(params%outfile), t, x, u)
    end if

  end do

  write(*,*) "Terminé. Fichiers générés: "//trim(params%outfile)//"_t=<t>.dat"
  write(*,*) "============================================================"

  !=======================================================
  ! 7) Erreur 
  !=======================================================
!   !=======================================================
!   ! Référence numérique sur maillage très fin (FIXE)
!   !=======================================================
!   dx_f = L / real(Nf_ref, pr)

!   allocate(x_f(Nf_ref), u_f(Nf_ref,nvar))
!   do i = 1, Nf_ref
!     x_f(i) = (i-0.5_pr) * dx_f
!   end do

!   ! CI fine
!   if (cas_test /= 4) then
!     do i = 1, Nf_ref
!       u_f(i,1) = C_init(i_CI, x_f(i))
!     end do
!   else
!     call init_euler(x_f, u_f, L)
!   end if

!   t_f  = 0._pr


!   do while (t_f < Tfinal)
!     if ( cas_test == 4 ) then
!       dt_f =  dt_CFL_euler(u_f, dx_f, CFL)
!     end if 

!     if (t_f + dt_f > Tfinal) then 
!       dt_f = Tfinal - t_f
!     end if 

!     call avancer_Rusanov(u_f, x_f, a, dx_f, dt_f, cl_periodique, t_f, i_CL, i_schema, cas_test)

!     t_f = t_f + dt_f
!   end do

! !=======================================================
! ! Restriction de la référence fine vers la grille grossière
! !=======================================================

!   if (mod(Nf_ref, nx) /= 0) then
!     write(*,*) "ERREUR: Nf_ref doit etre multiple de nx"
!     stop
!   end if

!   r_ref = Nf_ref / nx
!   allocate(u_ref(nx,nvar))
!   u_ref(:,:) = 0._pr

!   do i = 1, nx
!     i0 = (i-1)*r_ref ! indice de départ dans la grille fine
!     do k = 1, r_ref
!       do m = 1, nvar
!         u_ref(i,m) = u_ref(i,m) + u_f(i0+k,m) 
!       end do
!     end do
!     u_ref(i,:) = u_ref(i,:) / real(r_ref,pr) ! moyenne des cellules fines contenue dans la cellule grossière
!   end do

  if (cas_test /= 4) then
    allocate(u_ex(nx,1))
    u_ex(:,1) = U_exa(cas_test, x, t, a(1), L)
    errL2   = erreur_L2(u_ex, u, dx)
    errLinf = erreur_Linf(u_ex, u)
    errL1   = erreur_L1(u_ex, u, dx)
    
   ! errL2   = erreur_L2(u_ref, u, dx)
   ! errLinf = erreur_Linf(u_ref, u)
   ! errL1   = erreur_L1(u_ref, u, dx)
     
     write(*,'(A25, ES10.3)')  "  Erreur L1:", errL1
     write(*,'(A25, ES10.3)')  "  Erreur L2:",   errL2
     write(*,'(A25, ES10.3)')  "  Erreur Linf:", errLinf

     open(unit=16, file="erreurs.dat", status="unknown", position="append", iostat=ios)
     write(16,'(4ES16.8)') dx, errL1,errL2, errLinf
     close(16)

     deallocate(u_ex)
    
  else
    allocate(u_ex(nx,3))
    if (Tfinal==1.0_pr) then
      call init_euler(x, u_ex, L,rho0)
    else
      V = 1.0_pr; p = 2.0_pr
      
      do i = 1, nx
        rho = 1._pr +rho0*sin(2*pi*(x(i)-Tfinal))
        E      = p/(gamma - 1.0_pr) + 0.5*rho*V*V
        u_ex(i,1) = rho
        u_ex(i,2) = rho*V
        u_ex(i,3) = E
      end do  
    end if
  

    
    errL2   = erreur_L2(u_ex, u, dx)
    errLinf = erreur_Linf(u_ex, u)
    errL1   = erreur_L1(u_ex, u, dx)

  !  errL2   = erreur_L2(u_ref, u, dx)
  !  errLinf = erreur_Linf(u_ref, u)
  !  errL1   = erreur_L1(u_ref, u, dx)

    write(*,'(A25, ES10.3)')  "  Erreur L1:", errL1
    write(*,'(A25, ES10.3)')  "  Erreur L2:",   errL2
    write(*,'(A25, ES10.3)')  "  Erreur Linf:", errLinf

    open(unit=16, file="erreurs.dat", status="unknown", position="append", iostat=ios)
    write(16,'(4ES16.8)') dx, errL1,errL2, errLinf
    close(16)

    deallocate(u_ex)
  
  end if

  !=======================================================
  ! 8) Libération mémoire
  !=======================================================

  if (allocated(u))    deallocate(u)
  if (allocated(x))    deallocate(x)
  if (allocated(xeff)) deallocate(xeff)
  if (allocated(x_f)) deallocate(x_f, u_f, u_ref)

  write(*,*) "========================FIN du programme====================="

end program advection_rusanov
