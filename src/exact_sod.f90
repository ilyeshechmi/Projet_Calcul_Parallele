module exact_sod_mod
    use precision_mod
    implicit none
contains
 !==========================================================
  ! Solution exacte Sod (propre, gère t=0)
  !==========================================================
  subroutine sod_solution(U_ex, x, nx, t)
    integer,  intent(in)  :: nx
    real(pr), intent(in)  :: x(nx), t
    real(pr), intent(out) :: U_ex(nx,3)

    real(pr), parameter :: rhoL0 = 1.0_pr,   uL0 = 0.0_pr, pL0 = 1.0_pr
    real(pr), parameter :: rhoR0 = 0.125_pr, uR0 = 0.0_pr, pR0 = 0.1_pr

    real(pr) :: x0, s, rho, uvel, p
    real(pr) :: pstar, ustar
    integer  :: i

    x0 = 0.5_pr*(x(1) + x(nx))

    ! --- t = 0 : exact = CI ---
    if (t <= 0._pr) then
      do i = 1, nx
        if (x(i) < x0) then
          rho = rhoL0; uvel = uL0; p = pL0
        else
          rho = rhoR0; uvel = uR0; p = pR0
        end if
        U_ex(i,1) = rho
        U_ex(i,2) = rho*uvel
        U_ex(i,3) = p/(gamma-1._pr) + 0.5_pr*rho*uvel*uvel
      end do
      return
    end if

    call solve_star_region(pstar, ustar, rhoL0,uL0,pL0, rhoR0,uR0,pR0)

    do i = 1, nx
      s = (x(i) - x0)/t
      call sample_state(rho, uvel, p, s, rhoL0,uL0,pL0, rhoR0,uR0,pR0, pstar, ustar)

      U_ex(i,1) = rho
      U_ex(i,2) = rho*uvel
      U_ex(i,3) = p/(gamma-1._pr) + 0.5_pr*rho*uvel*uvel
    end do
  end subroutine sod_solution

  !========================
  ! Newton: pstar, ustar
  !========================
  subroutine solve_star_region(pstar, ustar, rhoL,uL,pL, rhoR,uR,pR0)
    real(pr), intent(out) :: pstar, ustar
    real(pr), intent(in)  :: rhoL,uL,pL, rhoR,uR,pR0
    real(pr) :: p_old, fL, fR, dfL, dfR, f, df
    integer  :: iter

    p_old = 0.5_pr*(pL + pR0)

    do iter = 1, 50
      call pressure_functions(fL, dfL, p_old, rhoL, pL)
      call pressure_functions(fR, dfR, p_old, rhoR, pR0)

      f  = fL + fR + (uR - uL)
      df = dfL + dfR

      pstar = p_old - f/df
      if (pstar < 1e-12_pr) pstar = 1e-12_pr

      if (abs(pstar - p_old)/pstar < 1e-10_pr) exit
      p_old = pstar
    end do

    ustar = 0.5_pr*(uL + uR + fR - fL)
  end subroutine solve_star_region

  subroutine pressure_functions(f, df, p, rho0, p0)
    real(pr), intent(out) :: f, df
    real(pr), intent(in)  :: p, rho0, p0
    real(pr) :: A, B, c0, pratio

    if (p <= p0) then
      c0     = sqrt(gamma*p0/rho0)
      pratio = p/p0
      f  = (2._pr*c0/(gamma-1._pr))*(pratio**((gamma-1._pr)/(2._pr*gamma)) - 1._pr)
      df = (1._pr/(rho0*c0))*pratio**(-(gamma+1._pr)/(2._pr*gamma))
    else
      A  = 2._pr/((gamma+1._pr)*rho0)
      B  = (gamma-1._pr)/(gamma+1._pr)*p0
      f  = (p-p0)*sqrt(A/(p+B))
      df = sqrt(A/(p+B))*(1._pr - 0.5_pr*(p-p0)/(p+B))
    end if
  end subroutine pressure_functions

  !========================
  ! Sample exact Riemann state at s = (x-x0)/t
  ! (version compacte, suffisante pour Sod standard)
  !========================
  subroutine sample_state(rho,uvel,p, s, rhoL,uL,pL, rhoR,uR,pR0, pstar, ustar)
    real(pr), intent(out) :: rho,uvel,p
    real(pr), intent(in)  :: s
    real(pr), intent(in)  :: rhoL,uL,pL, rhoR,uR,pR0, pstar, ustar
    real(pr) :: cL, cR, cLstar, SHL, STL, SL, SR, qL, qR, fact

    cL = sqrt(gamma*pL/rhoL)
    cR = sqrt(gamma*pR0/rhoR)

    if (s <= ustar) then
      ! ---- Left side
      if (pstar > pL) then
        qL = pstar/pL
        SL = uL - cL*sqrt((gamma+1._pr)/(2._pr*gamma)*qL + (gamma-1._pr)/(2._pr*gamma))
        if (s <= SL) then
          rho=rhoL; uvel=uL; p=pL
        else
          rho = rhoL * ( (qL + (gamma-1._pr)/(gamma+1._pr)) / ((gamma-1._pr)/(gamma+1._pr)*qL + 1._pr) )
          uvel=ustar; p=pstar
        end if
      else
        SHL    = uL - cL
        cLstar = cL*(pstar/pL)**((gamma-1._pr)/(2._pr*gamma))
        STL    = ustar - cLstar

        if (s <= SHL) then
          rho=rhoL; uvel=uL; p=pL
        elseif (s >= STL) then
          rho = rhoL*(pstar/pL)**(1._pr/gamma)
          uvel=ustar; p=pstar
        else
          uvel = (2._pr/(gamma+1._pr))*(cL + 0.5_pr*(gamma-1._pr)*uL + s)
          fact = 1._pr - (gamma-1._pr)/(2._pr*cL)*(uvel - uL)
          rho  = rhoL*fact**(2._pr/(gamma-1._pr))
          p    = pL *fact**(2._pr*gamma/(gamma-1._pr))
        end if
      end if
    else
      ! ---- Right side (Sod standard => shock)
      if (pstar > pR0) then
        qR = pstar/pR0
        SR = uR + cR*sqrt((gamma+1._pr)/(2._pr*gamma)*qR + (gamma-1._pr)/(2._pr*gamma))
        if (s >= SR) then
          rho=rhoR; uvel=uR; p=pR0
        else
          rho = rhoR * ( (qR + (gamma-1._pr)/(gamma+1._pr)) / ((gamma-1._pr)/(gamma+1._pr)*qR + 1._pr) )
          uvel=ustar; p=pstar
        end if
      else
        ! Rarefaction right (non standard Sod, but kept)
        rho = rhoR*(pstar/pR0)**(1._pr/gamma)
        uvel=ustar; p=pstar
      end if
    end if
  end subroutine sample_state




end module exact_sod_mod