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

    do i = 1, nx

      ! ---- Flux gauche F_{i-1/2} ----
      if (cl_periodique == 1) then
        ! Périodique
        if (i == 1) then
          fl = flux_rusanov(a, u(nx), u(1))
        else
          fl = flux_rusanov(a, u(i-1), u(i))
        end if
      else
        ! CL inflow à gauche
        if (i == 1) then
          fl = flux_rusanov(a, uG, u(1))   ! ghost gauche = uG
        else
          fl = flux_rusanov(a, u(i-1), u(i))
        end if
      end if

      ! ---- Flux droit F_{i+1/2} ----
      if (cl_periodique == 1) then
        ! Périodique
        if (i == nx) then
          fr = flux_rusanov(a, u(nx), u(1))   ! ghost droit = u(1)
        else
          fr = flux_rusanov(a, u(i), u(i+1))
        end if
      else
        ! Outflow libre à droite
        if (i == nx) then
          fr = flux_rusanov(a, u(nx), u(nx))  ! ghost droit = u(nx)
        else
          fr = flux_rusanov(a, u(i), u(i+1))
        end if
      end if

      unp1(i) = u(i) - (dt/dx)*( fr - fl )
    end do

    u = unp1
    deallocate(unp1)
  end subroutine avancer_Rusanov

end module schema_rusanov
