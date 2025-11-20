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
    real(pr), intent(in) :: u_exact(:), u_num(:)
    real(pr), intent(in) :: dx
    real(pr) :: errL2
    real(pr) :: sum, exact
    integer  :: i 

    sum   = 0._pr
    exact = 0._pr
    do i = 1,size(u_exact)
        sum   = sum   + dx*(u_exact(i)-u_num(i))**2
        exact = exact + dx*(u_exact(i))**2
    end do
    
    errL2 = sqrt(sum / exact)
end function erreur_L2



function erreur_Linf(u_exact, u_num) result(errLinf)
    real(pr), intent(in) :: u_exact(:), u_num(:)
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

end module 