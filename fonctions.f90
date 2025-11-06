module fonctions
    use precision
    implicit none 
contains 

    function gaussienne(x) result(res)
        real(pr), intent(in) :: x(:)
        real(pr) :: res(size(x,1))

        res = exp(-(x-5._pr)**2)
    end function gaussienne

end module 