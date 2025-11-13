module schema_rusanov
  use precision
  implicit none
contains

  function flux_rusanov(a, ul, ur) result(F)
    real(pr), intent(in) :: a, ul, ur
    real(pr) :: F, smax
    smax = abs(a)
    F = 0.5_pr*(a*ul + a*ur) - 0.5_pr*smax*(ur - ul)
  end function flux_rusanov

  ! Avance d'un pas

  subroutine avancer_Rusanov(u, a, dx, dt, uG, cl_periodique)
    real(pr), intent(inout) :: u(:)
    real(pr), intent(in)    :: a, dx, dt, uG
    integer,  intent(in)    :: cl_periodique   ! condition aux limites 0=non périodique, 1=périodique

    integer  :: nx, i
    real(pr) :: fl, fr
    real(pr), allocatable :: unp1(:)

    nx = size(u)
    allocate(unp1(nx))

    !Cas periodique 
    if ( cl_periodique ==1 ) then
     
      fl = flux_rusanov(a, u(nx), u(1)) !ghost gauche = u(nx)
      fr = flux_rusanov(a, u(1), u(2))
      unp1(1) = u(1) - (dt/dx)*( fr - fl )
      
      do i  = 2, nx-1
        
        fl=fr !Flux gauche Dejà calculé
        fr=flux_rusanov(a, u(i), u(i+1)) ! Flux droit
        unp1(i) = u(i) - (dt/dx)*( fr - fl )
      end do

      fl=fr 
      fr = flux_rusanov(a, u(nx), u(1)) !ghost droit = u(1)
      unp1(nx)=u(nx) -(dt/dx)*(fr-fl)

    else
      !Cas non périodique

      fl = flux_rusanov(a, uG, u(1)) ! ghost gauche = uG
      fr = flux_rusanov(a, u(1), u(2))
      unp1(1) = u(1) - (dt/dx)*( fr - fl )
      
      do i  = 2, nx-1
        fl=fr
        fr=flux_rusanov(a, u(i), u(i+1))
        unp1(i) = u(i) - (dt/dx)*( fr - fl )
      end do

      fl=fr
      fr = flux_rusanov(a, u(nx), u(nx)) !ghost droit = u(nx)
      unp1(nx)=u(nx) -(dt/dx)*(fr-fl)

    end if

    u = unp1
    deallocate(unp1)
  end subroutine avancer_Rusanov

end module schema_rusanov
