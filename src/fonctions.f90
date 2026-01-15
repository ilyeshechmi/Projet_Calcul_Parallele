module fonctions_mod
  use precision_mod
  use initialiser_mod
  implicit none
contains

  !========================
  ! Conditions initiales
  !========================
  function C_init(i_CL, x) result(u0)
    integer,  intent(in) :: i_CL
    real(pr), intent(in) :: x
    real(pr)             :: u0

    select case(i_CL)
    case(1)
      u0 = 0.0_pr
    case(2)
      u0 = exp( - (x - 5._pr)**2 )
    case(3)
      u0 = 1._pr
    case(4)
      u0 = sin(2.0_pr*pi*x)
    case default
      print *, "Erreur : Condition initiale inconnue"
      stop
    end select
  end function C_init

  !========================
  ! Conditions aux limites
  !========================
  function C_limite(i_CL, t) result(uL)
    integer,  intent(in) :: i_CL
    real(pr), intent(in) :: t
    real(pr)             :: uL

    select case(i_CL)
    case(1)
      uL = 1.0_pr
    case(2)
      uL = sin(2.0_pr*pi*t) + 1.0_pr
    case default
      print *, "Erreur : Condition limite inconnue"
      stop
    end select
  end function C_limite

  !========================
  ! Solution exacte advection
  !========================
  
  function U_exa(cas_test, x, t, a, L) result(u_ex)
    integer,  intent(in) :: cas_test
    real(pr), intent(in) :: x(:), t, a, L
    real(pr), allocatable :: u_ex(:)
    real(pr), allocatable :: xeff(:)
    integer :: i, nx
    integer :: i_CI, i_CL, CL_periodique
    real(pr) :: xb, tb

    call initialiser(cas_test, i_CL, i_CI, CL_periodique)

    nx = size(x)
    allocate(u_ex(nx), xeff(nx))

    if (CL_periodique == 1) then
      do i = 1, nx
        xeff(i) = modulo(x(i) - a*t, L)
        u_ex(i) = C_init(i_CI, xeff(i))
      end do
      deallocate(xeff)
      return
    end if

    ! --------- NON PERIODIQUE ----------
    if (a > 0._pr) then
      ! inflow à gauche : x=0
      do i = 1, nx
        xb = x(i) - a*t   ! pied de caractéristique à t=0
        if (xb >= 0._pr) then
          u_ex(i) = C_init(i_CI, xb)
        else
          tb = t - x(i)/a          ! temps d'entrée sur le bord gauche
          u_ex(i) = C_limite(i_CL, tb)
        end if
      end do

    else if (a < 0._pr) then
      ! inflow à droite : x=L
      do i = 1, nx
        xb = x(i) - a*t   ! = x + |a| t
        if (xb <= L) then
          u_ex(i) = C_init(i_CI, xb)
        else
          tb = t - (L - x(i))/(-a) ! temps d'entrée sur le bord droit
          u_ex(i) = C_limite(i_CL, tb)
        end if
      end do

    else
      ! a = 0 : pas d'advection
      do i = 1, nx
        u_ex(i) = C_init(i_CI, x(i))
      end do
    end if

    deallocate(xeff)
  end function U_exa


  !========================
  ! Erreurs (matriciel nx x dim)
  !========================
  function erreur_L1(u_exact, u_num, dx) result(errL1)
    real(pr), intent(in) :: u_exact(:,:), u_num(:,:)
    real(pr), intent(in) :: dx
    real(pr) :: errL1
    integer  :: i, m, nx, dim
    real(pr) :: num_sum, den_sum

    nx  = size(u_exact,1)
    dim = size(u_exact,2)

    num_sum = 0._pr
    den_sum = 0._pr

    do i = 1, nx
      do m = 1, dim
        num_sum = num_sum + dx * abs(u_exact(i,m) - u_num(i,m))
        den_sum = den_sum + dx * abs(u_exact(i,m))
      end do
    end do

    if (den_sum <= 0._pr) then
      errL1 = num_sum   
    else
      errL1 = num_sum / den_sum
    end if
  end function erreur_L1

  function erreur_L2(u_exact, u_num, dx) result(errL2)
    real(pr), intent(in) :: u_exact(:,:), u_num(:,:)
    real(pr), intent(in) :: dx
    real(pr) :: errL2
    integer  :: i, m, nx, dim
    real(pr) :: num_sum, den_sum

    nx  = size(u_exact,1)
    dim = size(u_exact,2)

    num_sum = 0._pr
    den_sum = 0._pr

    do i = 1, nx
      do m = 1, dim
        num_sum = num_sum + dx * (u_exact(i,m) - u_num(i,m))**2
        den_sum = den_sum + dx * (u_exact(i,m))**2
      end do
    end do

    if (den_sum <= 0._pr) then
      errL2 = sqrt(num_sum)   
    else
      errL2 = sqrt(num_sum / den_sum)
    end if
  end function erreur_L2

  function erreur_Linf(u_exact, u_num) result(errLinf)
    real(pr), intent(in) :: u_exact(:,:), u_num(:,:)
    real(pr) :: errLinf
    real(pr) :: max_err, max_ref

    max_err = maxval(abs(u_exact - u_num))
    max_ref = maxval(abs(u_exact))

    if (max_ref <= 0._pr) then
      errLinf = max_err
    else
      errLinf = max_err / max_ref
    end if
  end function erreur_Linf

  function source_term(U, x, t, cas_test) result(S)
  real(pr), intent(in) :: U(:), x, t
  integer,  intent(in) :: cas_test
  real(pr) :: S(size(U))
  real(pr) :: h, hU, heta, hw
  real(pr) :: vel, eta, w

  ! Valeur par défaut
  S(:) = 0._pr

  select case (cas_test)

  case (1,2,3)
    S(:) = 0._pr
  case (4)
    S(:) = 0._pr
  case (5)
    h    = U(1)
    hU   = U(2)
    heta = U(3)
    hw   = U(4)

    vel = hU / h
    eta = heta / h
    w   = hw / h

    S(1) = 0._pr
    S(2) = (lambda*h)/(epsilon*Re) &
         - (3._pr*vel)/(h*epsilon*Re)
    S(3) = hw
    S(4) = (1._pr/(alpha*beta))*(1._pr - eta/h) &
         - (9._pr*epsilon/(2._pr*beta*eta*Re))*w
    S(5) = 0._pr

  end select

end function source_term


function flux_advection(a, u) result(FU)
  real(pr), intent(in) :: a(:), u(:)
  real(pr) :: FU(size(u))

  FU(:) = a(:) * u(:)
end function flux_advection


function flux_saint_venant(U) result(Fu)
  real(pr), intent(in) :: U(5)
  real(pr) :: Fu(5)
  real(pr) :: h, hU, h_eta, h_w, p
  real(pr) :: Uv, eta
  real(pr) :: Press


  h     = U(1)
  hU    = U(2)
  h_eta = U(3)
  h_w   = U(4)
  p     = U(5)

  Uv  = hU / h
  eta = h_eta / h

  Press = (2._pr*lambda**2/225._pr)*h**5 + (cos(theta)/(2._pr*F**2))*h**2


  Fu(1) = hU

  Fu(2) = hU*Uv + Press &
         + (eta/alpha)*(1._pr - eta/h) &
         + (epsilon**2*kappa**2/(2._pr*F**2))*p**2

  Fu(3) = h_eta * Uv

  Fu(4) = h_w * Uv - (epsilon**2*kappa/(beta*F**2))*p

  Fu  (5) = p*Uv - h_w/h

end function flux_saint_venant



end module fonctions_mod
