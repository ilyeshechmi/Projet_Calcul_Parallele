module Euler_mod
  use precision_mod
  implicit none
contains

  !========================================================
  ! Conversion conservatif -> primitif
  ! Uc = [rho, rho*vel, E]
  !========================================================
  subroutine cons_to_prim(Uc, rho, vel, pres)
    real(pr), intent(in)  :: Uc(3)
    real(pr), intent(out) :: rho, vel, pres
    real(pr) :: E, kinetic

    rho = Uc(1)
    vel = Uc(2) / rho
    E   = Uc(3)

    kinetic = 0.5_pr * rho * vel * vel
    pres    = (gamma - 1._pr) * (E - kinetic)
  end subroutine cons_to_prim


  !========================================================
  ! Flux physique Euler F(U)
  !========================================================
  function flux_euler(Uc) result(F)
    real(pr), intent(in) :: Uc(3)
    real(pr) :: F(3)
    real(pr) :: rho, vel, pres

    call cons_to_prim(Uc, rho, vel, pres)

    F(1) = rho * vel
    F(2) = rho * vel * vel + pres
    F(3) = vel * (Uc(3) + pres)
  end function flux_euler


  !========================================================
  ! Vitesse max de propagation |u| + c
  !========================================================
  function max_wave_speed_euler(UcL, UcR) result(smax)
    real(pr), intent(in) :: UcL(3), UcR(3)
    real(pr) :: smax
    real(pr) :: rhoL, velL, presL
    real(pr) :: rhoR, velR, presR
    real(pr) :: cL, cR

    call cons_to_prim(UcL, rhoL, velL, presL)
    call cons_to_prim(UcR, rhoR, velR, presR)

    cL = sqrt(gamma * presL / rhoL)
    cR = sqrt(gamma * presR / rhoR)

    smax = max( abs(velL) + cL, abs(velR) + cR )
  end function max_wave_speed_euler


  !========================================================
  ! Flux numérique de Rusanov pour Euler
  !========================================================
  function flux_rusanov_euler(UcL, UcR) result(Fnum)
    real(pr), intent(in) :: UcL(3), UcR(3)
    real(pr) :: Fnum(3)
    real(pr) :: FL(3), FR(3)
    real(pr) :: smax

    FL = flux_euler(UcL)
    FR = flux_euler(UcR)

    smax = max_wave_speed_euler(UcL, UcR)

    Fnum = 0.5_pr * (FL + FR) - 0.5_pr * smax * (UcR - UcL)
  end function flux_rusanov_euler


  !========================================================
  ! CFL Euler : dt = CFL * dx / max(|u|+c)
  !========================================================
  function dt_CFL_euler(U, dx, CFL) result(dt)
    real(pr), intent(in) :: U(:,:), dx, CFL
    real(pr) :: dt
    integer :: i, nx
    real(pr) :: rho, vel, pres, c
    real(pr) :: smax
    real(pr) :: Uc(3)

    nx   = size(U,1)
    smax = 0._pr

    do i = 1, nx
       Uc = U(i,:)
       call cons_to_prim(Uc, rho, vel, pres)

       if (rho <= 0._pr .or. pres <= 0._pr) then
          write(*,*) "ERREUR: état non physique, i=", i
          write(*,*) "U =", Uc
          stop
       end if

       c    = sqrt(gamma * pres / rho)
       smax = max(smax, abs(vel) + c)
    end do

    dt = CFL * dx / smax
  end function dt_CFL_euler

end module Euler_mod
