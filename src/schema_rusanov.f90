module schema_rusanov_mod
  use precision_mod
  implicit none
contains

  function flux_rusanov(a, ul, ur) result(F)
    real(pr), intent(in) :: a, ul, ur
    real(pr) :: F, smax
    smax = abs(a)
    F = 0.5_pr*(a*ul + a*ur) - 0.5_pr*smax*(ur - ul)
  end function flux_rusanov

  ! Avance d'un pas

  subroutine avancer_Rusanov(u, a, dx, dt, cl_periodique, t, i_CL)
  use precision_mod
  use fonctions_mod        ! pour flux_rusanov
  use initialiser_mod      ! pour C_limite
  implicit none

  real(pr), intent(inout) :: u(:)         ! solution au temps n -> mise à jour vers n+1
  real(pr), intent(in)    :: a, dx, dt
  integer,  intent(in)    :: cl_periodique
  real(pr), intent(in)    :: t
  integer,  intent(in)    :: i_CL

  integer :: nx, i
  real(pr), allocatable :: unp1(:)
  real(pr) :: fl, fr
  real(pr) :: uL, uR

  nx = size(u)
  allocate(unp1(nx))

  !=============================
  ! Cas périodique
  !=============================
  if (cl_periodique == 1) then

     ! --- Point i=1 (bord gauche) ---
     fl = flux_rusanov(a, u(nx), u(1))    ! flux entre i=nx et i=1
     fr = flux_rusanov(a, u(1),  u(2))
     unp1(1) = u(1) - (dt/dx)*(fr - fl)

     ! --- Points internes ---
     do i = 2, nx-1
        fl = flux_rusanov(a, u(i-1), u(i))
        fr = flux_rusanov(a, u(i),   u(i+1))
        unp1(i) = u(i) - (dt/dx)*(fr - fl)
     end do

     ! --- Point i=nx (bord droit) ---
     fl = flux_rusanov(a, u(nx-1), u(nx))
     fr = flux_rusanov(a, u(nx),   u(1))
     unp1(nx) = u(nx) - (dt/dx)*(fr - fl)

  else
  !=============================
  ! Cas NON périodique
  !=============================

     !-------------------------------------------
     ! 1) Condition limite AMONT selon le signe de a
     !-------------------------------------------
     if (a > 0._pr) then
        uL = C_limite(i_CL, t)      ! entrée à gauche : x=0
        uR = u(nx)                  ! sortie à droite : extrapolation
     else
        uL = u(1)                   ! sortie à gauche
        uR = C_limite(i_CL, t)      ! entrée à droite : x=L
     end if

     ! --- Point i=1 ---
     fl = flux_rusanov(a, uL, u(1))
     fr = flux_rusanov(a, u(1), u(2))
     unp1(1) = u(1) - (dt/dx)*(fr - fl)

     ! --- Points internes ---
     do i = 2, nx-1
        fl = flux_rusanov(a, u(i-1), u(i))
        fr = flux_rusanov(a, u(i),   u(i+1))
        unp1(i) = u(i) - (dt/dx)*(fr - fl)
     end do

     ! --- Point i=nx ---
     fl = flux_rusanov(a, u(nx-1), u(nx))
     fr = flux_rusanov(a, u(nx),   uR)
     unp1(nx) = u(nx) - (dt/dx)*(fr - fl)

  end if

  ! Mise à jour finale
  u = unp1
  deallocate(unp1)

end subroutine avancer_Rusanov






end module 
