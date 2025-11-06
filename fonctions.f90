module fonctions
    use precision
    implicit none 
contains 

    function gaussienne(x) result(res)
        real(pr), intent(in) :: x(:)
        real(pr) :: res(size(x,1))

        res = exp(-(x-5._pr)**2)
    end function gaussienne
    
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
    real(8), intent(in) :: u_exact(:), u_num(:)
    real(8) :: errLinf
    real(8) :: max_err, max_u

    max_err = maxval( abs(u_exact - u_num) )
    max_u   = maxval( abs(u_exact) )

    errLinf = max_err / max_u

  end function erreur_Linf

end module 