module Euler_mod
  use precision_mod
  implicit none

contains

  !-------------------------------------------------------
  ! Conservé -> primitif : (rho, m, E) -> (rho, u, p)
  ! Uv(1) = rho, Uv(2) = rho*u, Uv(3) = E
  !-------------------------------------------------------
  subroutine cons_to_prim(Uv, rho, uvel, p, gamma)
    real(pr), intent(in)  :: Uv(3)
    real(pr), intent(in)  :: gamma
    real(pr), intent(out) :: rho, uvel, p
    real(pr) :: E, kinetic

    rho   = Uv(1)
    uvel  = Uv(2) / rho
    E     = Uv(3)
    kinetic = 0.5 * rho * uvel*uvel
    p     = (gamma - 1.0) * (E - kinetic)
  end subroutine cons_to_prim

  !-------------------------------------------------------
  ! Flux physique F(U)
  !-------------------------------------------------------
  subroutine flux_euler(Uv, F, gamma)
    real(pr), intent(in)  :: Uv(3)
    real(pr), intent(in)  :: gamma
    real(pr), intent(out) :: F(3)
    real(pr) :: rho, uvel, p, E

    call cons_to_prim(Uv, rho, uvel, p, gamma)
    E = Uv(3)

    F(1) = rho*uvel
    F(2) = rho*uvel*uvel + p
    F(3) = uvel*(E + p)
  end subroutine flux_euler

  !-------------------------------------------------------
  ! Vitesse max de propagation : |u| + c
  !-------------------------------------------------------
  function max_wave_speed(ULv, URv, gamma) result(smax)
    real(pr), intent(in) :: ULv(3), URv(3)
    real(pr), intent(in) :: gamma
    real(pr) :: smax
    real(pr) :: rhoL, uL, pL, rhoR, uR, pR, cL, cR

    call cons_to_prim(ULv, rhoL, uL, pL, gamma)
    call cons_to_prim(URv, rhoR, uR, pR, gamma)

    cL = sqrt(gamma * pL / rhoL)
    cR = sqrt(gamma * pR / rhoR)

    smax = max( abs(uL) + cL, abs(uR) + cR )
  end function max_wave_speed

  !-------------------------------------------------------
  ! Flux numérique de Rusanov pour Euler
  !-------------------------------------------------------
  subroutine flux_rusanov_euler(ULv, URv, Fnum, gamma)
    real(pr), intent(in)  :: ULv(3), URv(3)
    real(pr), intent(in)  :: gamma
    real(pr), intent(out) :: Fnum(3)
    real(pr) :: FL(3), FR(3), smax

    call flux_euler(ULv, FL, gamma)
    call flux_euler(URv, FR, gamma)
    smax = max_wave_speed(ULv, URv, gamma)

    Fnum(:) = 0.5*(FL(:) + FR(:)) - 0.5*smax*(URv(:) - ULv(:))
  end subroutine flux_rusanov_euler

  !-------------------------------------------------------
  ! Avancement en temps Rusanov pour Euler
  ! U est de taille (nx,3) : U(i,1)=rho, U(i,2)=rho*u, U(i,3)=E
  !-------------------------------------------------------
  subroutine avancer_rusanov_euler(U, dx, dt, gamma)
    real(pr), intent(inout) :: U(:,:)     ! (nx,3)
    real(pr), intent(in)    :: dx, dt, gamma
    integer :: nx, i
    real(pr), allocatable :: Unp1(:,:), ULv(:), URv(:), Fg(:), Fd(:)

    nx = size(U, 1)   ! nb de points d'espace

    allocate(Unp1(nx,3), ULv(3), URv(3), Fg(3), Fd(3))

    ! --- CL simples : on gèle les bords (copy) ---
    do i = 1, nx
       ! flux à gauche
       if (i == 1) then
          ULv = U(1,:)
          URv = U(1,:)
       else
          ULv = U(i-1,:)
          URv = U(i,:)
       end if
       call flux_rusanov_euler(ULv, URv, Fg, gamma)

       ! flux à droite
       if (i == nx) then
          ULv = U(nx,:)
          URv = U(nx,:)
       else
          ULv = U(i,:)
          URv = U(i+1,:)
       end if
       call flux_rusanov_euler(ULv, URv, Fd, gamma)

       Unp1(i,:) = U(i,:) - (dt/dx)*(Fd(:) - Fg(:))
    end do

    U(:,:) = Unp1(:,:)

    deallocate(Unp1, ULv, URv, Fg, Fd)
  end subroutine avancer_rusanov_euler

  !-------------------------------------------------------
  ! CI du tube de choc de Sod
  ! U est (nx,3) : (rho, rho*u, E)
  !-------------------------------------------------------
  subroutine init_euler_sod(x, U, gamma, L)
    real(pr), intent(in)  :: x(:), gamma, L
    real(pr), intent(out) :: U(:,:)      ! (nx,3)
    integer :: i, nx
    real(pr) :: rhoL, uL, pL, rhoR, uR, pR, E

    nx = size(x)

    ! États gauche et droite Sod
    rhoL = 1.0
    uL   = 0.0
    pL   = 1.0

    rhoR = 0.125
    uR   = 0.0
    pR   = 0.1

    do i = 1, nx
      if (x(i) < 0.5*L) then
        ! gauche
        E      = pL/(gamma - 1.0) + 0.5*rhoL*uL*uL
        U(i,1) = rhoL
        U(i,2) = rhoL*uL
        U(i,3) = E
      else
        ! droite
        E      = pR/(gamma - 1.0) + 0.5*rhoR*uR*uR
        U(i,1) = rhoR
        U(i,2) = rhoR*uR
        U(i,3) = E
      end if
    end do
  end subroutine init_euler_sod

end module Euler_mod
