module schema_rusanov_mod
   use precision_mod
   use fonctions_mod
   use Euler_mod
   use initialiser_mod    
   implicit none
contains

  !========================================================
  ! Flux de Rusanov pour l’advection linéaire scalaire
  !========================================================
  function flux_rusanov(a, ul, ur) result(F)
   real(pr), intent(in) :: a(:), ul(:), ur(:)
   real(pr) :: F(size(ul)), FL(Size(ul)), FR(Size(ul))
   
   FL = flux_advection(a, ul)
   FR = flux_advection(a, ur)

   F = 0.5_pr*(FL + Fr) - 0.5_pr*abs(a)*(ur - ul)
   
  end function flux_rusanov 

  !========================================================
  ! Flux numérique de Rusanov pour Euler
  !========================================================
  function flux_rusanov_euler(UcL, UcR) result(Fnum)
   real(pr), intent(in) :: UcL(:), UcR(:)
   real(pr) :: Fnum(Size(UcL))
   real(pr) :: FL(size(UcL)), FR(size(UcL))
   real(pr) :: smax

   FL = flux_euler(UcL)
   FR = flux_euler(UcR)

   smax = max_wave_speed_euler(UcL, UcR)

   Fnum = 0.5_pr * (FL + FR) - 0.5_pr * smax * (UcR - UcL)
  end function flux_rusanov_euler


  !========================================================
  ! Limiteur de pente : minmod(a, b)
  !========================================================
  function minmod(a, b) result(mm)
   real(pr), intent(in) :: a(:), b(:)
   real(pr) :: mm (size(a))
   integer :: i
   
   do i = 1, size(a)
      if (a(i)*b(i) <= 0._pr) then
         mm(i) = 0._pr
      else if (abs(a(i)) < abs(b(i))) then
         mm(i) = a(i)
      else
         mm(i) = b(i)
      end if
   end do 

  end function minmod

  !========================================================
  ! Avancement d’un pas de temps
  !   - muscl = 0 : schéma rusanov ordre 1
  !   - muscl = 1 : MUSCL + minmod
  !   - muscl = 2 : MUSCL-Hancock + minmod
  !   - BC périodiques ou inflow/outflow via C_limite
  !========================================================
  subroutine avancer_Rusanov(u,x, a, dx, dt, cl_periodique, t, i_CL, muscl,cas_test)
   real(pr), intent(inout) :: u(:,:)
   real(pr), intent(in)    :: a(:), dx, dt,x(:), t
   integer,  intent(in)    :: cl_periodique,cas_test,i_CL
   integer,  intent(in)    :: muscl     ! 0 : ordre 1, 1 : MUSCL
   
   real(pr), allocatable :: unp1(:,:)
   real(pr), allocatable :: pente(:,:)       ! pentes par cellule
   real(pr), allocatable :: uL(:,:), uR(:,:)   ! états aux interfaces
   real(pr), allocatable :: flux(:,:)        ! flux aux interfaces
   real(pr), allocatable :: uLc(:,:), uRc(:,:), uLh(:,:), uRh(:,:)
   real(pr), allocatable :: uH(:,:)
   integer :: nx, i,dim
   integer  :: il, ir                      ! indices cellule gauche/droite

   nx = size(u,1 )
   dim = size(u,2)

   allocate(unp1(nx,dim))
   allocate(uL(0:nx,dim), uR(0:nx,dim), flux(0:nx,dim))

   !=====================================================
   ! 1) Si MUSCL : calcul des pentes limitees
   !=====================================================
   if (muscl >= 1) then
      allocate(pente(nx,dim))

      if (cl_periodique == 1) then
         ! ----- CAS PERIODIQUE -----
         do i = 1, nx
           il = i - 1; if (il < 1)  il = nx
           ir = i + 1; if (ir > nx) ir = 1
           pente(i,:) = minmod(u(i,:) - u(il,:), u(ir,:) - u(i,:))
         end do
      else
         ! ----- CAS NON PERIODIQUE -----
         pente(1,:)  = 0._pr
         do i = 2, nx-1
           pente(i,:) = minmod(u(i,:) - u(i-1,:), u(i+1,:) - u(i,:))
         end do
         pente(nx,:) = 0._pr
      end if
   end if

   !=====================================================
   ! 2) Etats aux interfaces
   !=====================================================

   if (muscl == 0) then
      call build_interface_states(u, a, t, i_CL, 0, cas_test, cl_periodique, uL, uR)
   else
      ! MUSCL ou MUSCL–Hancock
      do i = 0, nx
         il = i;   if (il < 1) il = nx
         ir = i+1; if (ir > nx) ir = 1

         uL(i,:) = u(il,:) + 0.5_pr*pente(il,:)
         uR(i,:) = u(ir,:) - 0.5_pr*pente(ir,:)
      end do

      !-----------------------------------------------
      ! Hancock : prédiction à dt/2
      !-----------------------------------------------
      if (muscl == 2) then
         allocate(uLc(nx,dim), uRc(nx,dim), uLh(nx,dim), uRh(nx,dim))
         ! 1) reconstruction cellulaire au temps n
         do i = 1, nx
            uLc(i,:) = u(i,:) - 0.5_pr*pente(i,:)
            uRc(i,:) = u(i,:) + 0.5_pr*pente(i,:)
         end do

  ! 2) Hancock : prédiction à n+1/2 dans chaque cellule
         do i = 1, nx
            select case (cas_test)
               case (1,2,3)
               ! advection: F = a*u (a est un vecteur de taille dim ici)
                  uLh(i,:) = uLc(i,:) - 0.5_pr*(dt/dx) * ( flux_advection(a, uRc(i,:)) - flux_advection(a, uLc(i,:)) )
                  uRh(i,:) = uRc(i,:) - 0.5_pr*(dt/dx) * ( flux_advection(a, uRc(i,:)) - flux_advection(a, uLc(i,:)) )
               case (4)
               ! Euler: F = flux_euler(U)
                  uLh(i,:) = uLc(i,:) - 0.5_pr*(dt/dx) * ( flux_euler(uRc(i,:)) - flux_euler(uLc(i,:)) )
                  uRh(i,:) = uRc(i,:) - 0.5_pr*(dt/dx) * ( flux_euler(uRc(i,:)) - flux_euler(uLc(i,:)) )
            end select
         end do

  ! 3) construire les états aux interfaces à partir des états prédits
  !    interface i = (cellule i | cellule i+1) sur 0..nx
         if (cl_periodique == 1) then
            do i = 0, nx
               il = i;   if (il < 1) il = nx
               ir = i+1; if (ir > nx) ir = 1
               uL(i,:) = uRh(il,:)   ! côté droit prédît de la cellule gauche
               uR(i,:) = uLh(ir,:)   ! côté gauche prédît de la cellule droite
            end do
         else
            allocate(uH(nx,dim))
            uH(:,:) = 0.5_pr*(uLh(:,:) + uRh(:,:))
            call build_interface_states(uH, a, t+0.5_pr*dt, i_CL, 1, cas_test, cl_periodique, uL, uR, pente)
            deallocate(uH)
         end if
      deallocate(uLc, uRc, uLh, uRh)
      end if

   end if

   !=====================================================
   ! 3) Flux de Rusanov aux interfaces
   !=====================================================
   do i = 0, nx
     select case (cas_test)
     case (1,2,3)   
      flux(i,:) = flux_rusanov(a(:), uL(i,:), uR(i,:))
     case (4)
      flux(i,:) = flux_rusanov_euler(uL(i,:), uR(i,:))
      end select  
   end do

   !=====================================================
   ! 4) Mise à jour des cellules
   !    u^{n+1}_i = u^n_i - dt/dx ( F_{i+1/2} - F_{i-1/2} )
   !                      + terme source
   !      -> flux(i) = F_{i+1/2}
   !=====================================================
   do i = 1, nx
      unp1(i,:) = u(i,:) - (dt/dx) * ( flux(i,:) - flux(i-1,:) ) + dt * source_term(u(i,:), x(i), t, cas_test)
   end do

   !=====================================================
   ! 5) Mise à jour finale et nettoyage
   !=====================================================
   u = unp1

   deallocate(unp1, uL, uR, flux)
   if (allocated(pente)) deallocate(pente)

  end subroutine avancer_Rusanov
  !=====================================================
  ! Construction des états uL,uR aux interfaces (0..nx)
  !  - muscl = 0 : Rusanov ordre 1
  !  - muscl = 1 : reconstruction MUSCL via pente
  !  - périodique : wrap
  !  - non périodique :
  !      * Euler (cas_test=4) : extrapolation aux bords
  !      * advection : inflow/outflow via C_limite (a(1))
  !=====================================================
  subroutine build_interface_states(u, a, t, i_CL, muscl, cas_test, cl_periodique, uL, uR, pente)
   real(pr), intent(in)  :: u(:,:), a(:), t
   integer,  intent(in)  :: i_CL, muscl, cas_test, cl_periodique
   real(pr), intent(in), optional :: pente(:,:)
   real(pr), intent(out) :: uL(0:,:), uR(0:,:)
   
   integer :: nx, dim, i, il, ir
   logical :: use_muscl
   real(pr) :: uL_bc(size(u,2)), uR_bc(size(u,2))

   nx  = size(u,1)
   dim = size(u,2)
   use_muscl = (muscl == 1 .and. present(pente))

  !====================================================
  ! CAS PERIODIQUE
  !====================================================
  if (cl_periodique == 1) then
     do i = 0, nx
        il = i
        if (il < 1) il = nx

        ir = i + 1
        if (ir > nx) ir = 1

        if (use_muscl) then
           uL(i,:) = u(il,:) + 0.5_pr * pente(il,:)
           uR(i,:) = u(ir,:) - 0.5_pr * pente(ir,:)
        else
           uL(i,:) = u(il,:)
           uR(i,:) = u(ir,:)
        end if
     end do
     return
  end if

  !====================================================
  ! CAS NON PERIODIQUE : valeurs fantômes
  !====================================================
  if (cas_test == 4) then
     ! Euler : extrapolation aux deux bords
     uL_bc(:) = u(1,:)
     uR_bc(:) = u(nx,:)
  else
     ! Advection inflow / outflow (dim = 1 attendu)
     if (a(1) > 0._pr) then
        uL_bc(1) = C_limite(i_CL, t)
        uR_bc(1) = u(nx,1)
     else
        uL_bc(1) = u(1,1)
        uR_bc(1) = C_limite(i_CL, t)
     end if

     if (dim > 1) then
        uL_bc(2:dim) = u(1,2:dim)
        uR_bc(2:dim) = u(nx,2:dim)
     end if
  end if

  !====================================================
  ! Interfaces
  !====================================================

  ! i = 0 : (ghost gauche | cellule 1)
  uL(0,:) = uL_bc
  if (use_muscl) then
     uR(0,:) = u(1,:) - 0.5_pr * pente(1,:)
  else
     uR(0,:) = u(1,:)
  end if

  ! i = 1 .. nx-1 : (cellule i | cellule i+1)
  do i = 1, nx-1
     if (use_muscl) then
        uL(i,:) = u(i,:)   + 0.5_pr * pente(i,:)
        uR(i,:) = u(i+1,:) - 0.5_pr * pente(i+1,:)
     else
        uL(i,:) = u(i,:)
        uR(i,:) = u(i+1,:)
     end if
  end do

  ! i = nx : (cellule nx | ghost droite)
  if (use_muscl) then
     uL(nx,:) = u(nx,:) + 0.5_pr * pente(nx,:)
  else
     uL(nx,:) = u(nx,:)
  end if
  uR(nx,:) = uR_bc

end subroutine build_interface_states


end module schema_rusanov_mod

