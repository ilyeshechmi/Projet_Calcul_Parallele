module schema_rusanov_mod
  use precision_mod
   use fonctions_mod
   use Euler_mod
  implicit none
contains

  !========================================================
  ! Flux de Rusanov pour l’advection linéaire scalaire
  !========================================================
  function flux_rusanov(a, ul, ur) result(F)
    real(pr), intent(in) :: a(:), ul(:), ur(:)
    real(pr) :: F(size(ul))

    F(:) = 0.5_pr*(a(:)*ul(:) + a(:)*ur(:)) - 0.5_pr*abs(a(:))*(ur(:) - ul(:))
  end function flux_rusanov 

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
  !   - BC périodiques ou inflow/outflow via C_limite
  !========================================================
  subroutine avancer_Rusanov(u,x, a, dx, dt, cl_periodique, t, i_CL, muscl,cas_test)
    use initialiser_mod    
    use fonctions_mod   
    implicit none

    real(pr), intent(inout) :: u(:,:)
    real(pr), intent(in)    :: a(:), dx, dt,x(:)
    integer,  intent(in)    :: cl_periodique,cas_test
    real(pr), intent(in)    :: t
    integer,  intent(in)    :: i_CL
    integer,  intent(in)    :: muscl     ! 0 : ordre 1, 1 : MUSCL

    integer :: nx, i,dim
    real(pr), allocatable :: unp1(:,:)
    real(pr), allocatable :: pente(:,:)       ! pentes par cellule
    real(pr), allocatable :: uL(:,:), uR(:,:)   ! états aux interfaces
    real(pr), allocatable :: flux(:,:)        ! flux aux interfaces
    integer  :: il, ir                      ! indices cellule gauche/droite

    nx = size(u,1 )
    dim = size(u,2)

    allocate(unp1(nx,dim))
    allocate(uL(0:nx,dim), uR(0:nx,dim), flux(0:nx,dim))

    !=====================================================
    ! 1) Si MUSCL : calcul des pentes limitees
    !=====================================================
    if (muscl == 1) then
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
    ! 2) COnstruction des états uL, uR aux interfaces
    !=====================================================
    if (muscl == 1) then
       call build_interface_states(u, a, t, i_CL, muscl, cas_test, cl_periodique, uL, uR,pente)
    else
       ! pas de pente requis
       call build_interface_states(u, a, t, i_CL, muscl, cas_test, cl_periodique, uL, uR)
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
  use precision_mod
  use fonctions_mod
  implicit none

  real(pr), intent(in)  :: u(:,:), a(:), t
  integer,  intent(in)  :: i_CL, muscl, cas_test, cl_periodique
  real(pr), intent(in), optional :: pente(:,:)
  real(pr), intent(out) :: uL(0:,:), uR(0:,:)

  integer :: nx, dim, i, il, ir
  real(pr) :: uL_bc(size(u,2)), uR_bc(size(u,2))
  logical  :: use_muscl

  nx = size(u,1);  dim = size(u,2)
  use_muscl = (muscl == 1)

  !====================================================
  ! CAS PERIODIQUE
  !====================================================
  if (cl_periodique == 1) then
     do i = 0, nx
        il = merge(i, nx, i>=1);   ir = merge(i+1, 1, i<nx)
        uL(i,:) = u(il,:) + merge(0.5_pr*pente(il,:), 0._pr, use_muscl)
        uR(i,:) = u(ir,:) - merge(0.5_pr*pente(ir,:), 0._pr, use_muscl)
     end do
     return
  end if

  !====================================================
  ! CAS NON PERIODIQUE : valeurs fantômes
  !====================================================
  if (cas_test == 4) then
     ! Euler : extrapolation
     uL_bc = u(1,:)
     uR_bc = u(nx,:)
  else
     ! Advection inflow/outflow
     if (a(1) > 0._pr) then
        uL_bc(1) = C_limite(i_CL, t);  uR_bc(1) = u(nx,1)
     else
        uL_bc(1) = u(1,1);             uR_bc(1) = C_limite(i_CL, t)
     end if
     if (dim > 1) then
        uL_bc(2:) = u(1,2:);  uR_bc(2:) = u(nx,2:)
     end if
  end if

  !====================================================
  ! Interfaces
  !====================================================
  ! i = 0
  uL(0,:) = uL_bc
  uR(0,:) = u(1,:) - merge(0.5_pr*pente(1,:), 0._pr, use_muscl)

  ! i = 1 .. nx-1
  do i = 1, nx-1
     uL(i,:) = u(i,:)   + merge(0.5_pr*pente(i,:),   0._pr, use_muscl)
     uR(i,:) = u(i+1,:) - merge(0.5_pr*pente(i+1,:), 0._pr, use_muscl)
  end do

  ! i = nx
  uL(nx,:) = u(nx,:) + merge(0.5_pr*pente(nx,:), 0._pr, use_muscl)
  uR(nx,:) = uR_bc

end subroutine build_interface_states

end module schema_rusanov_mod

