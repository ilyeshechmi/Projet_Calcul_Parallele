module fonctions_mod
    use precision_mod
    use initialiser_mod
    implicit none 
contains 

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
      u0 = exp( - (x - 5._pr)**2 )

    case(3)
       u0 =1._pr
    case(4)
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
        uL = sin(2.0_pr*pi*t)+1.0_pr

      case default
        print *, "Erreur : Condition limite inconnue"
        stop  
    end select

  end function C_limite


    
function erreur_L2(u_exact, u_num, dx) result(errL2)
    real(pr), intent(in) :: u_exact(:,:), u_num(:,:)
    real(pr), intent(in) :: dx
    real(pr) :: errL2
    real(pr) :: num_sum, den_sum
    integer  :: i, m, nx,dim

    nx = size(u_exact,1)
    dim = size(u_exact,2)
    num_sum = 0._pr
    den_sum = 0._pr

    do m = 1, dim
        do i = 1, nx
            num_sum = num_sum + dx * (u_exact(i,m) - u_num(i,m))**2
            den_sum = den_sum + dx * (u_exact(i,m))**2
        end do
    end do

    errL2 = sqrt(num_sum / den_sum)
end function erreur_L2




function erreur_Linf(u_exact, u_num) result(errLinf)
    real(pr), intent(in) :: u_exact(:,:), u_num(:,:)
    real(pr) :: errLinf
    real(pr) :: max_err, max_u

    max_err = maxval( abs(u_exact - u_num) )
    max_u   = maxval( abs(u_exact) )

    errLinf = max_err / max_u
end function erreur_Linf


function U_exa(cas_test, x, t, a, L) result(u_ex)

    integer,  intent(in) :: cas_test
    real(pr), intent(in) :: x(:), t, a, L
    real(pr), allocatable :: u_ex(:), xeff(:)
    integer :: i, nx
    integer :: i_CI, i_CL, CL_periodique

    call initialiser(cas_test, i_CL, i_CI, CL_periodique)

    nx = size(x)
    allocate(u_ex(nx), xeff(nx))
    do i = 1, nx
        if (CL_periodique == 1) then
            ! --- Condition périodique ---
            xeff(i) = modulo(x(i) - a*t, L)
            u_ex(i) = C_init(i_CI, xeff(i))

        else
            ! --- Condition non périodique (advection vers la droite) ---
            if (x(i) - a*t > 0._pr) then
                u_ex(i) = C_init(i_CI, x(i) - a*t)
            else
                u_ex(i) = C_limite(i_CL, t)
            end if

        end if
    end do

    deallocate(xeff)

end function U_exa


subroutine sod_solution(U_exa, x, nx, t)
  implicit none
  integer, intent(in) :: nx
  real(pr), intent(in) :: x(nx), t
  real(pr), intent(out) :: U_exa(nx,3)

  ! États gauche et droite
  real(pr) :: rhoL, uL, pL
  real(pr) :: rhoR0, uR, pR0    

  ! États étoile
  real(pr) :: pstar, ustar

  ! Vitesses du son et/ou dans la région étoile
  real(pr) :: cL, cR, cLstar

  ! Vitesses caractéristiques des ondes
  real(pr) :: shl, stl, stc, shR

  ! État point-échantillon
  real(pr) :: rho, u, p
  integer :: i

  rhoL = 1._pr
  uL   = 0._pr
  pL   = 1._pr

  rhoR0 = 0.125_pr
  uR    = 0._pr
  pR0   = 0.1_pr

  ! Résolution du problème de Riemann (Newton pour p*)
  call solve_star_region(pstar, ustar, rhoL, uL, pL, rhoR0, uR, pR0)

  ! Vitesses du son
  cL = sqrt(gamma * pL / rhoL)
  cR = sqrt(gamma * pR0 / rhoR0)

  ! c* dans région étoile gauche
  cLstar = cL * (pstar/pL)**((gamma - 1._pr) / (2._pr * gamma))

  ! Vitesses des ondes
  shl = uL - cL                  ! tête de la raréfaction gauche
  stl = ustar - cLstar           ! queue de la raréfaction gauche
  stc = ustar                    ! discontinuité de contact
  shR = uR + cR * sqrt( ((gamma+1._pr)/(2._pr*gamma))*(pstar/pR0) &
                        + (gamma-1._pr)/(2._pr*gamma) )  ! choc droit


  ! Boucle sur les points du maillage
  do i = 1, nx
     call sample_state(rho, u, p, x(i)/t, &
                       rhoL, uL, pL, rhoR0, uR, pR0, &
                       pstar, ustar, cL, cR, cLstar, &
                       shl, stl, stc, shR)

     ! Conversion variables primitives en variables conservatives
     U_exa(i,1) = rho
     U_exa(i,2) = rho*u
     U_exa(i,3) = p/(gamma-1._pr) + 0.5_pr*rho*u*u
  end do

end subroutine sod_solution

!  Résolution du Riemann exact (Newton) pour p* et u*
subroutine solve_star_region(pstar, ustar, rhoL,uL,pL, rhoR,uR,pR0)
  implicit none
  real(pr), intent(out) :: pstar, ustar
  real(pr), intent(in)  :: rhoL,uL,pL, rhoR,uR,pR0

  real(pr) :: p_old, fL, fR, dfL, dfR, f, df
  integer :: iter


  p_old = 0.5_pr*(pL + pR0)

  do iter = 1, 40
     call pressure_functions(fL,dfL, p_old, rhoL, pL)
     call pressure_functions(fR,dfR, p_old, rhoR, pR0)

     f  = fL + fR + (uR - uL)
     df = dfL + dfR

     pstar = p_old - f/df

     ! Convergence (seuil uniformisé en pr)
     if (abs((pstar - p_old)/pstar) < 1e-8_pr) exit

     if (pstar < 0._pr) pstar = 1e-6_pr
     p_old = pstar
  end do

  ! État étoile (écoulement)
  ustar = 0.5_pr*(uL + uR + fR - fL)

end subroutine solve_star_region


! Fonction f(p) pour choc ou raréfaction

subroutine pressure_functions(f, df, p, rho0, p0)
  implicit none
  real(pr), intent(out) :: f, df
  real(pr), intent(in)  :: p, rho0, p0

  real(pr) :: A, B, c0

  if (p <= p0) then
     ! Raréfaction
     c0 = sqrt(gamma*p0/rho0)
     f  = (2._pr*c0/(gamma-1._pr)) * ( (p/p0)**((gamma-1._pr)/(2._pr*gamma)) - 1._pr )
     df = (1._pr/(rho0*c0)) * (p/p0)**(-(gamma+1._pr)/(2._pr*gamma))
  else
     ! Choc
     A = 2._pr / ((gamma+1._pr)*rho0)
     B = (gamma-1._pr)/(gamma+1._pr) * p0
     f  = (p - p0) * sqrt(A/(p + B))
     df = sqrt(A/(p + B)) * (1._pr - 0.5_pr*(p - p0)/(p + B))
  end if

end subroutine pressure_functions


!  Échantillonnage de la solution exacte au point xi = x/t
subroutine sample_state(rho,u,p, xi, rhoL,uL,pL, rhoR,uR,pR0, &
                        pstar,ustar, cL,cR, cLstar, shl,stl,stc,shR)
  implicit none

  real(pr), intent(out) :: rho, u, p
  real(pr), intent(in)  :: xi
  real(pr), intent(in)  :: rhoL,uL,pL, rhoR,uR,pR0
  real(pr), intent(in)  :: pstar,ustar
  real(pr), intent(in)  :: cL,cR,cLstar
  real(pr), intent(in)  :: shl,stl,stc,shR

  real(pr) :: fact

  if (xi < ustar) then
     ! ------------------ CÔTÉ GAUCHE ------------------
     if (pstar > pL) then
        ! ---- Choc gauche ----
        if (xi < shl) then
           rho = rhoL ; u = uL ; p = pL
        else
           rho = rhoL*((pstar/pL + (gamma-1._pr)/(gamma+1._pr)) / &
                      ((gamma-1._pr)/(gamma+1._pr)*pstar/pL + 1._pr))
           u = ustar ; p = pstar
        end if

     else
        ! ---- Raréfaction gauche ----
        if (xi < shl) then
           rho = rhoL ; u = uL ; p = pL
        elseif (xi > stl) then
           rho = rhoL*(pstar/pL)**(1._pr/gamma)
           u = ustar ; p = pstar
        else
           u    = (2._pr/(gamma+1._pr))*(cL + xi)
           fact = 1._pr - (gamma-1._pr)/2._pr * u/cL
           rho  = rhoL * fact**(2._pr/(gamma-1._pr))
           p    = pL   * fact**(2._pr*gamma/(gamma-1._pr))
        end if
     end if

  else
     ! ------------------ CÔTÉ DROIT ------------------
     if (pstar > pR0) then
        ! ---- Choc droit ----
        if (xi > shR) then
           rho = rhoR ; u = uR ; p = pR0
        else
           rho = rhoR*((pstar/pR0 + (gamma-1._pr)/(gamma+1._pr)) / &
                      ((gamma-1._pr)/(gamma+1._pr)*pstar/pR + 1._pr))
           u = ustar ; p = pstar
        end if

     else
        ! ---- Raréfaction droite ----
        if (xi > uR + cR) then
           rho = rhoR ; u = uR ; p = pR0
        elseif (xi < stc) then
           rho = rhoR*(pstar/pR0)**(1._pr/gamma)
           u = ustar ; p = pstar
        else
           u    = (2._pr/(gamma+1._pr))*(-cR + xi)
           fact = 1._pr + (gamma-1._pr)/2._pr * u/cR
           rho  = rhoR * fact**(2._pr/(gamma-1._pr))
           p    = pR   * fact**(2._pr*gamma/(gamma-1._pr))
        end if
     end if
  end if

end subroutine sample_state

end module 