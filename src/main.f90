program advection_rusanov
    
  use precision 
  use fonctions
  use donnees
  use schema_rusanov
  implicit none

  type(Parametres) :: params
  integer  :: nx, nsteps, n, save_every, i ,cl_periodique
  real(pr) :: L, a, dt, Tfinal, dx, t, cfl, uG
  real(pr), allocatable :: x(:), u(:,:),u_ex(:)
  real(pr) ::  errL2, errLinf
  
  
  ! --- Lecture paramètres ---
  call lire_parametres("parametres.txt", params)
  nx         = params%nx
  L          = params%L
  a          = params%a
  Tfinal     = params%T_final
  save_every = params%save_every
  CFL        = params%CFL
  dx         = params%dx
  cl_periodique = params%cl_periodique
  

  ! --- Grille uniforme (centres de mailles) ---
  allocate(x(nx), u(1,nx))
  x  = [ ( (i-0.5_pr)*dx, i=1,nx ) ]

  ! --- Condition initiale ---
  ! u = sin( 2.0_pr * pi * x / L )
  u(1,:) = gaussienne(x)
  
  ! --- Condition à Gauche ---
  uG = 0._pr


  ! --- Sortie initiale ---
  t = 0.0_pr
  call ecrire(trim(params%outfile), t, x, u(1,:))

  ! --- Boucle en temps ---
  dt= CFL*dx/abs(a)

  nsteps = ceiling( Tfinal / dt )
  do n = 1, nsteps
    call avancer_Rusanov(u(1,:), a, dx, dt,uG,cl_periodique)
    t = t + dt
    ! Écriure des résultats
    if (mod(n, save_every) == 0 .or. n == nsteps) then
      call ecrire(trim(params%outfile), t, x, u(1,:)) 
    end if
  end do

  write(*,*) "Terminé. Fichiers générés: "//trim(params%outfile)//"_t=<t>.dat"
  write(*,*) "============================================================"

  ! Calcul de l'erreur
  allocate(u_ex(nx))
  do i = 1, nx
    if (x(i) - a*t > 0._pr) then
      u_ex(i) = exp(-(x(i) - a*t-5._pr)**2)
    else
      u_ex(i) = uG
    end if
  end do

  errL2   = erreur_L2(u_ex, u(1,:), dx)
  errLinf = erreur_Linf(u_ex, u(1,:))

  open(unit=20, file=trim(params%outfile)//'_erreurs.dat', status='replace', action='write')
  write(20, '(A, F12.6)') 'Erreur L2   = ', errL2
  write(20, '(A, F12.6)') 'Erreur Linf = ', errLinf
  close(20)
  write(*,*) 'Erreurs calculées et sauvegardées dans ', trim(params%outfile)//'_erreurs.dat'

  ! --- Sauvegarde des solutions finale et exacte ---
open(unit=30, file=trim(params%outfile)//'_compare.dat', status='replace')
write(30, '(A)') '# x   u_num   u_exact'
do i = 1, nx
    write(30, '(3E20.10)') x(i), u(1,i), u_ex(i)
end do
close(30)

! --- Génération automatique du script gnuplot ---
open(unit=31, file='plot_compare.gnu', status='replace')
write(31, '(A)') "set term pngcairo size 1000,600"
write(31, '(A)') "set output '"//trim(params%outfile)//"_compare.png'"
write(31, '(A)') "set xlabel 'x'"
write(31, '(A)') "set ylabel 'u(x,t)'"
write(31, '(A)') "set title 'Comparaison Rusanov vs solution exacte à t="//trim(adjustl(to_string(t)))//"'"
write(31, '(A)') "set grid"
write(31, '(A)') "plot '"//trim(params%outfile)//"_compare.dat' using 1:2 with lines lw 2 lc 'blue' title 'u numérique', \"
write(31, '(A)') "     '"//trim(params%outfile)//"_compare.dat' using 1:3 with lines lw 2 lc 'red' title 'u exacte'"
close(31)

! --- Lancement de gnuplot ---
call execute_command_line("gnuplot plot_compare.gnu")

write(*,*) "Graphique généré : ", trim(params%outfile)//"_compare.png"



  deallocate(x, u,u_ex)

contains
  function to_string(r) result(str)
    real(pr), intent(in) :: r
    character(len=50) :: str
    write(str, '(F8.4)') r
  end function to_string


end program advection_rusanov
