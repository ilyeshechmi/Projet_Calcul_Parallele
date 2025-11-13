module fonctions_mod
    use precision_mod
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
      u0 = exp( - ((x - 0.5_pr)/0.1_pr)**2 )

    case(3)
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
      uL = sin(2.0_pr*pi*t)

    case default
      print *, "Erreur : Condition limite inconnue"
      stop  
    end select

  end function C_limite


    
 function erreur_L2(u_exact, u_num, dx) result(errL2)
    implicit none
    real(pr), intent(in) :: u_exact(:), u_num(:)
    real(pr), intent(in) :: dx
    real(pr) :: errL2
    real(pr) :: sum,exact
    integer  :: i 

    sum=0._pr
    exact=0._pr
    do i = 1,size(u_exact)
        sum = sum+ dx*(u_exact(i)-u_num(i))**2
        exact = exact+ dx*(u_exact(i))**2
    end do
    
    errL2=sqrt(sum/exact)

  end function erreur_L2


  function erreur_Linf(u_exact, u_num) result(errLinf)
    implicit none
    real(pr), intent(in) :: u_exact(:), u_num(:)
    real(pr) :: errLinf
    real(pr) :: max_err, max_u

    max_err = maxval( abs(u_exact - u_num) )
    max_u   = maxval( abs(u_exact) )

    errLinf = max_err / max_u

  end function erreur_Linf

  function f_1(a,U) result(R)
    implicit none
    real(pr),intent(in)  :: U(:,:)
    real(pr),intent(in)  :: a(:)
    real(pr), allocatable :: R(:,:)
    integer :: i 
    allocate(R(size(a(:)),size(U(1,:))))

    do i = 1,size(a(:))
      R(i,:)=a(i)*U(i,:)
    end do

  end function f_1

  ! function U_exa() result(u_ex)


  !    allocate(u_ex(nx))
  ! do i = 1, nx
  !   if ( cl_periodique ==1 ) then
  !     xeff(i)=modulo(x(i)-a*t,10._pr)
  !     !u_ex(i) = exp(-(xeff(i)-5._pr)**2)
  !     u_ex(i) = sin(2._pr *pi*xeff(i)/L)
  !   else
  !     if (x(i) - a*t > 0._pr) then
  !       !u_ex(i) = exp(-(x(i) - a*t-5._pr)**2)
  !       u_ex(i) = sin(2._pr *pi*(x(i)-a*t)/L)
  !     else
  !       u_ex(i) = uG
  !     end if
  !   end if
   
end module 