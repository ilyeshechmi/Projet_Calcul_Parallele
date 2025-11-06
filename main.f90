program advection_rusanov
    
  use precision 
  use fonctions
  use donnees
  use schema_rusanov
  implicit none

  type(Parametres) :: params
  integer  :: nx, nsteps, n, save_every, i
  real(pr) :: L, a, dt, Tfinal, dx, t, cfl, uG
  real(pr), allocatable :: x(:), u(:)

  ! --- Lecture paramètres ---
  call lire_parametres("parametres.txt", params)
  nx         = params%nx
  L          = params%L
  a          = params%a
  dt         = params%dt
  Tfinal     = params%T_final
  save_every = params%save_every
  CFL        = params%CFL
  dx         = params%dx

  ! --- Grille uniforme (centres de mailles) ---
  allocate(x(nx), u(nx))
  x  = [ ( (i-0.5_pr)*dx, i=1,nx ) ]

  ! --- Condition initiale ---
  ! u = sin( 2.0_pr * pi * x / L )
  u = gaussienne(x)
  
  ! --- Condition à Gauche ---
  uG = 1._pr


  ! --- Sortie initiale ---
  t = 0.0_pr
  call ecrire(trim(params%outfile), t, x, u)

  ! --- Boucle en temps ---
  nsteps = ceiling( Tfinal / dt )
  do n = 1, nsteps
    call avancer_Rusanov(u, a, dx, dt,uG)
    t = t + dt
    ! Écriure des résultats
    if (mod(n, save_every) == 0 .or. n == nsteps) then
      call ecrire(trim(params%outfile), t, x, u) 
    end if
  end do

  write(*,*) "Terminé. Fichiers générés: "//trim(params%outfile)//"_t=<t>.dat"
  write(*,*) "============================================================"

  deallocate(x, u)
end program advection_rusanov
