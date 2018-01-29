!======================================================================
!! R.Nair NCAR/scd 10/02
!! Rewritten 06/2016 by Francois Hebert, SIParCS/Cornell
!! Dictionary for the indices/directions/words:
!! index 1 = -y = south
!! index 2 = +x = east
!! index 3 = +y = north
!! index 4 = -x = west
!! index 5 = -z = bottom
!! index 6 = +z = top
!! Implements the weak form DG method for evaluating the
!! semi-discrete RHS of the advection equation.
!! Computes the RHS of `dU/dt = rhs'
!======================================================================

Module Euler_flux3d_mod 
  use basic_mod
  use physical_const_mod
  use element_mod
  use dg3d_rhs_mod, only : Horizontal_strong_grads
  use testcases_nh3d_mod, only : Compute_Pressure,Compute_SoundSpeed,Pressure,SoundSpeed
  use cudafor

  implicit none
  private
  public :: Make_Euler3d_Flux , Compute_Ext_Vars,nh_flux_maker,Make_ext_normalflux,Apply_Noflux_TopbotBC,Make_Flux, nh_maker, Make_normalflux, Apply_TopbotBC
  public :: Recover_State_Vars, Store_Hydrostatic_edgevars

contains

!!============================================================================================

  Subroutine nh_flux_maker()
  implicit none

    real(kind=real_kind), dimension(neq) :: svec
    real(kind=real_kind) :: u,v,w,spd,prp,the,rho,fcori
    integer :: eq,e,i,j,k,ie,je,ke,nlx,nly,nlz

   nlx = size(vars,1)
   nly = size(vars,2)
   nlz = size(vars,3)

   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T

    ! vec(:,:,:,:) = vars%vec(:,:,:,:) 

        do ke = 1, nlz
        do je = 1, nly
        do ie = 1, nlx
       do k = 1, nx
        do j = 1, nx
        do i = 1, nx

     u = vars(ie,je,ke)%u(i,j,k) 
     v = vars(ie,je,ke)%v(i,j,k) 
     w = vars(ie,je,ke)%w(i,j,k) 

     rho = vars(ie,je,ke)%rho(i,j,k) 
     the = vars(ie,je,ke)%the(i,j,k) 
     prp = vars(ie,je,ke)%prp(i,j,k)  
     spd = vars(ie,je,ke)%spd(i,j,k) 

     !Full pressure in the discretization 
     !prs(:,:,:) =  vars%prb(:,:,:) 
     !Call Horizontal_strong_grads(prs,prs_x,prs_y) 

     fcori = vars(ie,je,ke)%fcori(i,j,k) 

  !Construct state vecror from the evolved vars 
              svec(1) = rho - vars(ie,je,ke)%rhb(i,j,k)
              svec(2) = rho * u
              svec(3) = rho * v
              svec(4) = rho * w
              svec(5) = rho * the - vars(ie,je,ke)%rtb(i,j,k)

  ! Flux x-direction 
    vars(ie,je,ke)%flx(i,j,k,1) = svec(2) 
    vars(ie,je,ke)%flx(i,j,k,2) = svec(2) * u + prp
    vars(ie,je,ke)%flx(i,j,k,3) = svec(2) * v
    vars(ie,je,ke)%flx(i,j,k,4) = svec(2) * w
    vars(ie,je,ke)%flx(i,j,k,5) = svec(2) * the

  ! Flux y-direction 
    vars(ie,je,ke)%fly(i,j,k,1) = svec(3) 
    vars(ie,je,ke)%fly(i,j,k,2) = svec(3) * u
    vars(ie,je,ke)%fly(i,j,k,3) = svec(3) * v + prp
    vars(ie,je,ke)%fly(i,j,k,4) = svec(3) * w
    vars(ie,je,ke)%fly(i,j,k,5) = svec(3) * the

  ! Flux z-direction 
    vars(ie,je,ke)%flz(i,j,k,1) = svec(4) 
    vars(ie,je,ke)%flz(i,j,k,2) = svec(4) * u
    vars(ie,je,ke)%flz(i,j,k,3) = svec(4) * v
    vars(ie,je,ke)%flz(i,j,k,4) = svec(4) * w + prp
    vars(ie,je,ke)%flz(i,j,k,5) = svec(4) * the

  ! Source term         
    vars(ie,je,ke)%src(i,j,k,1) = 0.0D0 
    vars(ie,je,ke)%src(i,j,k,2) =  fcori*v !- prs_x(:,:,:) 
    vars(ie,je,ke)%src(i,j,k,3) = -fcori*u !- prs_y(:,:,:) 
    vars(ie,je,ke)%src(i,j,k,4) = -svec(1) * grv 
    vars(ie,je,ke)%src(i,j,k,5) = 0.0D0

        end do
        end do
        end do

        end do
        end do
        end do

end subroutine nh_flux_maker

!!============================================================================================

 attributes(global) Subroutine nh_maker()
  implicit none

    real(kind=real_kind), dimension(neq) :: svec
    real(kind=real_kind) :: u,v,w,spd,prp,the,rho,fcori
    integer :: eq,e,i,j,k,ie,je,ke,nlx,nly,nlz

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)

   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T

    ! vec(:,:,:,:) = vars%vec(:,:,:,:) 

!        do ke = 1, nlz
!        do je = 1, nly
!        do ie = 1, nlx
!       do k = 1, nx
!        do j = 1, nx
!        do i = 1, nx
i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z

     u = vars_d(ie,je,ke)%u(i,j,k) 
     v = vars_d(ie,je,ke)%v(i,j,k) 
     w = vars_d(ie,je,ke)%w(i,j,k) 

     rho = vars_d(ie,je,ke)%rho(i,j,k) 
     the = vars_d(ie,je,ke)%the(i,j,k) 
     prp = vars_d(ie,je,ke)%prp(i,j,k)  
     spd = vars_d(ie,je,ke)%spd(i,j,k) 

     !Full pressure in the discretization 
     !prs(:,:,:) =  vars%prb(:,:,:) 
     !Call Horizontal_strong_grads(prs,prs_x,prs_y) 

     fcori = vars_d(ie,je,ke)%fcori(i,j,k) 

  !Construct state vecror from the evolved vars 
              svec(1) = rho - vars_d(ie,je,ke)%rhb(i,j,k)
              svec(2) = rho * u
              svec(3) = rho * v
              svec(4) = rho * w
              svec(5) = rho * the - vars_d(ie,je,ke)%rtb(i,j,k)

  ! Flux x-direction 
    vars_d(ie,je,ke)%flx(i,j,k,1) = svec(2) 
    vars_d(ie,je,ke)%flx(i,j,k,2) = svec(2) * u + prp
    vars_d(ie,je,ke)%flx(i,j,k,3) = svec(2) * v
    vars_d(ie,je,ke)%flx(i,j,k,4) = svec(2) * w
    vars_d(ie,je,ke)%flx(i,j,k,5) = svec(2) * the

  ! Flux y-direction 
    vars_d(ie,je,ke)%fly(i,j,k,1) = svec(3) 
    vars_d(ie,je,ke)%fly(i,j,k,2) = svec(3) * u
    vars_d(ie,je,ke)%fly(i,j,k,3) = svec(3) * v + prp
    vars_d(ie,je,ke)%fly(i,j,k,4) = svec(3) * w
    vars_d(ie,je,ke)%fly(i,j,k,5) = svec(3) * the

  ! Flux z-direction 
    vars_d(ie,je,ke)%flz(i,j,k,1) = svec(4) 
    vars_d(ie,je,ke)%flz(i,j,k,2) = svec(4) * u
    vars_d(ie,je,ke)%flz(i,j,k,3) = svec(4) * v
    vars_d(ie,je,ke)%flz(i,j,k,4) = svec(4) * w + prp
    vars_d(ie,je,ke)%flz(i,j,k,5) = svec(4) * the

  ! Source term         
    vars_d(ie,je,ke)%src(i,j,k,1) = 0.0D0 
    vars_d(ie,je,ke)%src(i,j,k,2) =  fcori*v !- prs_x(:,:,:) 
    vars_d(ie,je,ke)%src(i,j,k,3) = -fcori*u !- prs_y(:,:,:) 
    vars_d(ie,je,ke)%src(i,j,k,4) = -svec(1) * grv 
    vars_d(ie,je,ke)%src(i,j,k,5) = 0.0D0

!        end do
!        end do
!        end do
!
!        end do
!        end do
!        end do

end subroutine nh_maker

!!============================================================================================
  Subroutine Make_Euler3d_Flux()
  implicit none
!    type(element_t), intent(inout) :: vars
!    type(edge_t), intent(inout), dimension(6) :: edges
!    integer, intent(in) :: ke, ie,je  

    integer :: eq, e, i,j,k,ie,je,ke,nlx,nly,nlz 

   nlx = size(vars,1)
   nly = size(vars,2)
   nlz = size(vars,3)

        do ke = 1, nlz
        do je = 1, nly
        do ie = 1, nlx
        do j = 1, nx
        do i = 1, nx

    edgevars(ie,je,ke)%edges(1)%spd_int(i,j) = vars(ie,je,ke)%spd(i,1,j)
    edgevars(ie,je,ke)%edges(2)%spd_int(i,j) = vars(ie,je,ke)%spd(nx,i,j)
    edgevars(ie,je,ke)%edges(3)%spd_int(i,j) = vars(ie,je,ke)%spd(i,nx,j) 
    edgevars(ie,je,ke)%edges(4)%spd_int(i,j) = vars(ie,je,ke)%spd(1,i,j)
    edgevars(ie,je,ke)%edges(5)%spd_int(i,j) = vars(ie,je,ke)%spd(i,j,1)
    edgevars(ie,je,ke)%edges(6)%spd_int(i,j) = vars(ie,je,ke)%spd(i,j,nx)


   do eq = 1, neq 
    edgevars(ie,je,ke)%edges(1)%normflux_int(i,j,eq) = vars(ie,je,ke)%fly(i,1,j,eq)
    edgevars(ie,je,ke)%edges(2)%normflux_int(i,j,eq) = vars(ie,je,ke)%flx(nx,i,j,eq)
    edgevars(ie,je,ke)%edges(3)%normflux_int(i,j,eq) = vars(ie,je,ke)%fly(i,nx,j,eq)
    edgevars(ie,je,ke)%edges(4)%normflux_int(i,j,eq) = vars(ie,je,ke)%flx(1,i,j,eq)
    edgevars(ie,je,ke)%edges(5)%normflux_int(i,j,eq) = vars(ie,je,ke)%flz(i,j,1,eq)
    edgevars(ie,je,ke)%edges(6)%normflux_int(i,j,eq) = vars(ie,je,ke)%flz(i,j,nx,eq)
   end do 

        end do
        end do
        end do
        end do
        end do

!   Call Make_ext_normalflux()

!   Call Apply_Noflux_TopbotBC() 

  End Subroutine Make_Euler3d_Flux

!!============================================================================================
 attributes(global) Subroutine Make_Flux()
  implicit none
!    type(element_t), intent(inout) :: vars
!    type(edge_t), intent(inout), dimension(6) :: edges
!    integer, intent(in) :: ke, ie,je  

    integer :: eq, e, i,j,k,ie,je,ke,nlx,nly,nlz 

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)

!        do ke = 1, nlz
!        do je = 1, nly
!        do ie = 1, nlx
!        do j = 1, nx
!        do i = 1, nx
i = threadIdx%x
j = threadIdx%y
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z

    edgevars_d(ie,je,ke)%edges(1)%spd_int(i,j) = vars_d(ie,je,ke)%spd(i,1,j)
    edgevars_d(ie,je,ke)%edges(2)%spd_int(i,j) = vars_d(ie,je,ke)%spd(nx,i,j)
    edgevars_d(ie,je,ke)%edges(3)%spd_int(i,j) = vars_d(ie,je,ke)%spd(i,nx,j) 
    edgevars_d(ie,je,ke)%edges(4)%spd_int(i,j) = vars_d(ie,je,ke)%spd(1,i,j)
    edgevars_d(ie,je,ke)%edges(5)%spd_int(i,j) = vars_d(ie,je,ke)%spd(i,j,1)
    edgevars_d(ie,je,ke)%edges(6)%spd_int(i,j) = vars_d(ie,je,ke)%spd(i,j,nx)


   do eq = 1, neq 
    edgevars_D(ie,je,ke)%edges(1)%normflux_int(i,j,eq) = vars_d(ie,je,ke)%fly(i,1,j,eq)
    edgevars_d(ie,je,ke)%edges(2)%normflux_int(i,j,eq) = vars_d(ie,je,ke)%flx(nx,i,j,eq)
    edgevars_d(ie,je,ke)%edges(3)%normflux_int(i,j,eq) = vars_d(ie,je,ke)%fly(i,nx,j,eq)
    edgevars_d(ie,je,ke)%edges(4)%normflux_int(i,j,eq) = vars_d(ie,je,ke)%flx(1,i,j,eq)
    edgevars_d(ie,je,ke)%edges(5)%normflux_int(i,j,eq) = vars_d(ie,je,ke)%flz(i,j,1,eq)
    edgevars_d(ie,je,ke)%edges(6)%normflux_int(i,j,eq) = vars_d(ie,je,ke)%flz(i,j,nx,eq)
   end do 

!        end do
!        end do
!        end do
!        end do
!        end do

!   Call Make_ext_normalflux()

!   Call Apply_Noflux_TopbotBC() 

  End Subroutine Make_Flux

!!============================================================================================
  Subroutine Apply_Noflux_TopbotBC()
  implicit none
    integer :: ie,je, ke ,nlx,nly,nlz,i,j

   nlx = size(vars,1)
   nly = size(vars,2)
   nlz = size(vars,3)

        do ke = 1, nlz
        do je = 1, nly
        do ie = 1, nlx
        do j = 1, nx
        do i = 1, nx

     if (ke == 1) then 
          edgevars(ie,je,ke)%edges(5)%vec_ext(i,j,:) =  edgevars(ie,je,ke)%edges(5)%vec_int(i,j,:)
          edgevars(ie,je,ke)%edges(5)%vec_ext(i,j,4) = -edgevars(ie,je,ke)%edges(5)%vec_int(i,j,4)

          edgevars(ie,je,ke)%edges(5)%normflux_ext(i,j,:) = -edgevars(ie,je,ke)%edges(5)%normflux_int(i,j,:)
          edgevars(ie,je,ke)%edges(5)%normflux_ext(i,j,4) =  edgevars(ie,je,ke)%edges(5)%normflux_int(i,j,4)
          edgevars(ie,je,ke)%edges(5)%spd_ext(i,j) =  edgevars(ie,je,ke)%edges(5)%spd_int(i,j)
     endif 

     if (ke == nez) then 
          edgevars(ie,je,ke)%edges(6)%vec_ext(i,j,:) =  edgevars(ie,je,ke)%edges(6)%vec_int(i,j,:)
          edgevars(ie,je,ke)%edges(6)%vec_ext(i,j,4) = -edgevars(ie,je,ke)%edges(6)%vec_int(i,j,4)

          edgevars(ie,je,ke)%edges(6)%normflux_ext(i,j,:) = -edgevars(ie,je,ke)%edges(6)%normflux_int(i,j,:)
          edgevars(ie,je,ke)%edges(6)%normflux_ext(i,j,4) =  edgevars(ie,je,ke)%edges(6)%normflux_int(i,j,4)
          edgevars(ie,je,ke)%edges(6)%spd_ext(i,j) = edgevars(ie,je,ke)%edges(6)%spd_int(i,j)
     endif 

end do
end do
end do
end do
end do

!    if (ie == nex) then !east  
!         edge(2)%vec_ext(:,:,:) =  edge(2)%vec_int(:,:,:)
!         edge(2)%vec_ext(:,:,2) = -edge(2)%vec_int(:,:,2)

!         edge(2)%normflux_ext(:,:,:) = -edge(2)%normflux_int(:,:,:)
!         edge(2)%normflux_ext(:,:,2) =  edge(2)%normflux_int(:,:,2)
!         edge(2)%spd_ext(:,:) = edge(2)%spd_int(:,:)
!    endif
 
!    if (je == ney) then !north 
!         edge(3)%vec_ext(:,:,:) =  edge(3)%vec_int(:,:,:)
!         edge(3)%vec_ext(:,:,3) = -edge(3)%vec_int(:,:,3)

!         edge(3)%normflux_ext(:,:,:) = -edge(3)%normflux_int(:,:,:)
!         edge(3)%normflux_ext(:,:,3) =  edge(3)%normflux_int(:,:,3)
!         edge(3)%spd_ext(:,:) = edge(3)%spd_int(:,:)
!    endif
 
  End Subroutine Apply_Noflux_TopbotBC

!!============================================================================================
 attributes(global) Subroutine Apply_TopbotBC()
  implicit none
    integer :: ie,je, ke ,nlx,nly,nlz,i,j

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)

!        do ke = 1, nlz
!        do je = 1, nly
!        do ie = 1, nlx
!        do j = 1, nx
!        do i = 1, nx
i = threadIdx%x
j = threadIdx%y
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z

     if (ke == 1) then 
          edgevars_d(ie,je,ke)%edges(5)%vec_ext(i,j,:) =  edgevars_d(ie,je,ke)%edges(5)%vec_int(i,j,:)
          edgevars_d(ie,je,ke)%edges(5)%vec_ext(i,j,4) = -edgevars_d(ie,je,ke)%edges(5)%vec_int(i,j,4)

          edgevars_d(ie,je,ke)%edges(5)%normflux_ext(i,j,:) = -edgevars_d(ie,je,ke)%edges(5)%normflux_int(i,j,:)
          edgevars_d(ie,je,ke)%edges(5)%normflux_ext(i,j,4) =  edgevars_d(ie,je,ke)%edges(5)%normflux_int(i,j,4)
          edgevars_d(ie,je,ke)%edges(5)%spd_ext(i,j) =  edgevars_d(ie,je,ke)%edges(5)%spd_int(i,j)
     endif 

     if (ke == nez) then 
          edgevars_d(ie,je,ke)%edges(6)%vec_ext(i,j,:) =  edgevars_d(ie,je,ke)%edges(6)%vec_int(i,j,:)
          edgevars_d(ie,je,ke)%edges(6)%vec_ext(i,j,4) = -edgevars_d(ie,je,ke)%edges(6)%vec_int(i,j,4)

          edgevars_d(ie,je,ke)%edges(6)%normflux_ext(i,j,:) = -edgevars_d(ie,je,ke)%edges(6)%normflux_int(i,j,:)
          edgevars_d(ie,je,ke)%edges(6)%normflux_ext(i,j,4) =  edgevars_d(ie,je,ke)%edges(6)%normflux_int(i,j,4)
          edgevars_d(ie,je,ke)%edges(6)%spd_ext(i,j) = edgevars_d(ie,je,ke)%edges(6)%spd_int(i,j)
     endif 

!end do
!end do
!end do
!end do
!end do

!    if (ie == nex) then !east  
!         edge(2)%vec_ext(:,:,:) =  edge(2)%vec_int(:,:,:)
!         edge(2)%vec_ext(:,:,2) = -edge(2)%vec_int(:,:,2)

!         edge(2)%normflux_ext(:,:,:) = -edge(2)%normflux_int(:,:,:)
!         edge(2)%normflux_ext(:,:,2) =  edge(2)%normflux_int(:,:,2)
!         edge(2)%spd_ext(:,:) = edge(2)%spd_int(:,:)
!    endif
 
!    if (je == ney) then !north 
!         edge(3)%vec_ext(:,:,:) =  edge(3)%vec_int(:,:,:)
!         edge(3)%vec_ext(:,:,3) = -edge(3)%vec_int(:,:,3)

!         edge(3)%normflux_ext(:,:,:) = -edge(3)%normflux_int(:,:,:)
!         edge(3)%normflux_ext(:,:,3) =  edge(3)%normflux_int(:,:,3)
!         edge(3)%spd_ext(:,:) = edge(3)%spd_int(:,:)
!    endif
 
  End Subroutine Apply_TopbotBC

!! ============================================================================================
  Subroutine Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)
  implicit none
    integer, intent(in) :: face
    real(kind=real_kind), intent(in)  :: vec_ext(neq),rhb,rtb,prb

    real(kind=real_kind), intent(out) :: u_ext, v_ext,w_ext, &
                                                          the_ext,prp_ext,spd_ext
    real(kind=real_kind) :: rho_ext, prs_ext

    ! NH variables along exterior boundaries 
                    rho_ext = vec_ext(1)   + rhb
                      u_ext = vec_ext(2) / rho_ext
                      v_ext = vec_ext(3) / rho_ext
                      w_ext = vec_ext(4) / rho_ext
                    the_ext = (vec_ext(5) + rtb)/rho_ext

!       Call  Compute_Ext_Vars(rho_ext,the_ext, prs_ext,spd_ext)
prs_ext = Compute_Pressure(rho_ext,the_ext)
spd_ext = Compute_SoundSpeed(prs_ext,rho_ext)

                    prp_ext = prs_ext - prb

  End  Subroutine  Extract_ext_Vars

!! ============================================================================================
 attributes(device) Subroutine Extract_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)
  implicit none
    integer, intent(in) :: face
    real(kind=real_kind), intent(in)  :: vec_ext(neq),rhb,rtb,prb

    real(kind=real_kind), intent(out) :: u_ext, v_ext,w_ext, &
                                                          the_ext,prp_ext,spd_ext
    real(kind=real_kind) :: rho_ext, prs_ext

    ! NH variables along exterior boundaries 
                    rho_ext = vec_ext(1)   + rhb
                      u_ext = vec_ext(2) / rho_ext
                      v_ext = vec_ext(3) / rho_ext
                      w_ext = vec_ext(4) / rho_ext
                    the_ext = (vec_ext(5) + rtb)/rho_ext

!       Call  Compute_Ext_Vars(rho_ext,the_ext, prs_ext,spd_ext)
!prs_ext = Compute_Pressure(rho_ext,the_ext)
!spd_ext = Compute_SoundSpeed(prs_ext,rho_ext)
call Pressure(rho_ext,the_ext,prs_ext)
call SoundSpeed(prs_ext,rho_ext,spd_ext)

                    prp_ext = prs_ext - prb

  End  Subroutine  Extract_Vars

!! ============================================================================================

   Subroutine  Make_ext_normalflux()
    implicit none

    real(kind=real_kind) :: vec_ext(neq)
    real(kind=real_kind) :: u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb
    integer :: face,i,j,ie,je,ke,nlx,nly,nlz

   nlx = size(vars,1)
   nly = size(vars,2)
   nlz = size(vars,3)

do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx
do j = 1, nx
do i = 1, nx
    !South face (y-flux)  
    !------------------
    do  face = 1,3, 2 
        vec_ext(:) = edgevars(ie,je,ke)%edges(face)%vec_ext(i,j,:)
        rhb = edgevars(ie,je,ke)%edges(face)%rhb(i,j)
        rtb = edgevars(ie,je,ke)%edges(face)%rtb(i,j)
        prb = edgevars(ie,je,ke)%edges(face)%prb(i,j)
     Call   Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,1) =  vec_ext(3)    !rho*v 
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,2) =  vec_ext(3) * u_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,3) =  vec_ext(3) * v_ext + prp_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,4) =  vec_ext(3) * w_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,5) =  vec_ext(3) * the_ext
     edgevars(ie,je,ke)%edges(face)%spd_ext(i,j) =  spd_ext
    end do 

    !North (y-flux)  
    !------------------
   !  face = 3 
   !     vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
   !  Call   Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)
!
!     edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,3)
!     edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,3) * u_ext(:,:)
!     edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,3) * v_ext(:,:) + prp_ext(:,:)
!     edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,3) * w_ext(:,:)
!     edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,3) * the_ext(:,:)
!     edge(face)%spd_ext(:,:) =  spd_ext(:,:)

    !East  face (x-flux)  
    !------------------
    do  face = 2,4, 2 
        vec_ext(:) = edgevars(ie,je,ke)%edges(face)%vec_ext(i,j,:)
        rhb = edgevars(ie,je,ke)%edges(face)%rhb(i,j)
        rtb = edgevars(ie,je,ke)%edges(face)%rtb(i,j)
        prb = edgevars(ie,je,ke)%edges(face)%prb(i,j)
     Call   Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,1) =  vec_ext(2)    !rho*u 
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,2) =  vec_ext(2) * u_ext + prp_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,3) =  vec_ext(2) * v_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,4) =  vec_ext(2) * w_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,5) =  vec_ext(2) * the_ext
     edgevars(ie,je,ke)%edges(face)%spd_ext(i,j) =  spd_ext
    end do 

!   !West (x-flux)  
!   !------------------
!    face = 4
!       vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
!    Call   Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

!    edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,2)
!    edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,2) * u_ext(:,:) + prp_ext(:,:)
!    edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,2) * v_ext(:,:)
!    edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,2) * w_ext(:,:)
!    edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,2) * the_ext(:,:)
!    edge(face)%spd_ext(:,:) =  spd_ext(:,:)

    !Bottom (z-flux)  
    !------------------
    do face = 5, 6 
        vec_ext(:) = edgevars(ie,je,ke)%edges(face)%vec_ext(i,j,:)
        rhb = edgevars(ie,je,ke)%edges(face)%rhb(i,j)
        rtb = edgevars(ie,je,ke)%edges(face)%rtb(i,j)
        prb = edgevars(ie,je,ke)%edges(face)%prb(i,j)
     Call   Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,1) =  vec_ext(4)    !rho*w 
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,2) =  vec_ext(4) * u_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,3) =  vec_ext(4) * v_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,4) =  vec_ext(4) * w_ext + prp_ext
     edgevars(ie,je,ke)%edges(face)%normflux_ext(i,j,5) =  vec_ext(4) * the_ext
     edgevars(ie,je,ke)%edges(face)%spd_ext(i,j) =  spd_ext
    end do 
    
    !Top (z-flux)  
    !------------------
!    face = 6
!       vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
!    Call   Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

!    edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,4)
!    edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,4) * u_ext(:,:) 
!    edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,4) * v_ext(:,:)
!    edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,4) * w_ext(:,:) + prp_ext(:,:)
!    edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,4) * the_ext(:,:)
!    edge(face)%spd_ext(:,:) =  spd_ext(:,:)
!   
end do
end do
end do
end do
end do

 End Subroutine  Make_ext_normalflux

!! ============================================================================================

  attributes(global) Subroutine Make_normalflux()
    implicit none

    real(kind=real_kind) :: vec_ext(neq)
    real(kind=real_kind) :: u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb
    integer :: face,i,j,ie,je,ke,nlx,nly,nlz

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)

!do ke = 1, nlz
!do je = 1, nly
!do ie = 1, nlx
!do j = 1, nx
!do i = 1, nx
i = threadIdx%x
j = threadIdx%y
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z

    !South face (y-flux)  
    !------------------
    do  face = 1,3, 2 
        vec_ext(:) = edgevars_d(ie,je,ke)%edges(face)%vec_ext(i,j,:)
        rhb = edgevars_d(ie,je,ke)%edges(face)%rhb(i,j)
        rtb = edgevars_d(ie,je,ke)%edges(face)%rtb(i,j)
        prb = edgevars_d(ie,je,ke)%edges(face)%prb(i,j)
     Call   Extract_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,1) =  vec_ext(3)    !rho*v 
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,2) =  vec_ext(3) * u_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,3) =  vec_ext(3) * v_ext + prp_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,4) =  vec_ext(3) * w_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,5) =  vec_ext(3) * the_ext
     edgevars_d(ie,je,ke)%edges(face)%spd_ext(i,j) =  spd_ext
    end do 

    !North (y-flux)  
    !------------------
   !  face = 3 
   !     vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
   !  Call   Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)
!
!     edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,3)
!     edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,3) * u_ext(:,:)
!     edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,3) * v_ext(:,:) + prp_ext(:,:)
!     edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,3) * w_ext(:,:)
!     edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,3) * the_ext(:,:)
!     edge(face)%spd_ext(:,:) =  spd_ext(:,:)

    !East  face (x-flux)  
    !------------------
    do  face = 2,4, 2 
        vec_ext(:) = edgevars_d(ie,je,ke)%edges(face)%vec_ext(i,j,:)
        rhb = edgevars_d(ie,je,ke)%edges(face)%rhb(i,j)
        rtb = edgevars_d(ie,je,ke)%edges(face)%rtb(i,j)
        prb = edgevars_d(ie,je,ke)%edges(face)%prb(i,j)
     Call   Extract_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,1) =  vec_ext(2)    !rho*u 
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,2) =  vec_ext(2) * u_ext + prp_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,3) =  vec_ext(2) * v_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,4) =  vec_ext(2) * w_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,5) =  vec_ext(2) * the_ext
     edgevars_d(ie,je,ke)%edges(face)%spd_ext(i,j) =  spd_ext
    end do 

!   !West (x-flux)  
!   !------------------
!    face = 4
!       vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
!    Call   Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

!    edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,2)
!    edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,2) * u_ext(:,:) + prp_ext(:,:)
!    edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,2) * v_ext(:,:)
!    edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,2) * w_ext(:,:)
!    edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,2) * the_ext(:,:)
!    edge(face)%spd_ext(:,:) =  spd_ext(:,:)

    !Bottom (z-flux)  
    !------------------
    do face = 5, 6 
        vec_ext(:) = edgevars_d(ie,je,ke)%edges(face)%vec_ext(i,j,:)
        rhb = edgevars_d(ie,je,ke)%edges(face)%rhb(i,j)
        rtb = edgevars_d(ie,je,ke)%edges(face)%rtb(i,j)
        prb = edgevars_d(ie,je,ke)%edges(face)%prb(i,j)
     Call   Extract_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,1) =  vec_ext(4)    !rho*w 
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,2) =  vec_ext(4) * u_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,3) =  vec_ext(4) * v_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,4) =  vec_ext(4) * w_ext + prp_ext
     edgevars_d(ie,je,ke)%edges(face)%normflux_ext(i,j,5) =  vec_ext(4) * the_ext
     edgevars_d(ie,je,ke)%edges(face)%spd_ext(i,j) =  spd_ext
    end do 
    
    !Top (z-flux)  
    !------------------
!    face = 6
!       vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
!    Call   Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

!    edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,4)
!    edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,4) * u_ext(:,:) 
!    edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,4) * v_ext(:,:)
!    edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,4) * w_ext(:,:) + prp_ext(:,:)
!    edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,4) * the_ext(:,:)
!    edge(face)%spd_ext(:,:) =  spd_ext(:,:)
!   
!end do
!end do
!end do
!end do
!end do

 End Subroutine  Make_normalflux

!=======================================================================================================!
   Subroutine Compute_Ext_Vars(rho_ext,the_ext, prs_ext,spd_ext) 
    Implicit None
 
    real(kind=real_kind), dimension(nx,nx), intent(in) ::   rho_ext,the_ext
    real(kind=real_kind), dimension(nx,nx), intent(out) ::  prs_ext,spd_ext
    real(kind=real_kind) ::  rho,the,prs 
    integer :: i,j 

            do j = 1, nx
              do i = 1, nx
               rho = rho_ext(i,j) 
               the = the_ext(i,j) 

               prs = Compute_Pressure(rho,the)
               prs_ext(i,j) = prs

               spd_ext(i,j) = Compute_SoundSpeed(prs,rho) 
                
              end do 
            end do 

   End Subroutine Compute_Ext_Vars 
!=======================================================================================================!
   Subroutine Compute_int_3dvars(rho,the,prs,spd)
    Implicit None

    real(kind=real_kind), dimension(nx,nx,nx), intent(in) ::   rho,the
    real(kind=real_kind), dimension(nx,nx,nx), intent(out) ::  prs,spd
    real(kind=real_kind) ::  ro,pr
    integer :: i,j,k
           do k = 1, nx
            do j = 1, nx
              do i = 1, nx
               ro = rho(i,j,k)
!               th = the(i,j,k)

               pr = Compute_Pressure(ro,the(i,j,k))
               prs(i,j,k) = pr

               spd(i,j,k) = Compute_SoundSpeed(pr,ro)

              end do
            end do
           end do

   End Subroutine Compute_int_3dvars 

!=======================================================================================================!
   Subroutine Volume_to_Surf(psi,psi_edg)
    Implicit None
 
    real(kind=real_kind), intent(in) ::  psi(nx,nx,nx)
    real(kind=real_kind), intent(out) :: psi_edg(6,nx,nx) 

          psi_edg(1,:,:) = psi(:,1,:)
          psi_edg(2,:,:) = psi(nx,:,:)
          psi_edg(3,:,:) = psi(:,nx,:)
          psi_edg(4,:,:) = psi(1,:,:)
          psi_edg(5,:,:) = psi(:,:,1)
          psi_edg(6,:,:) = psi(:,:,nx)

   End Subroutine Volume_to_Surf

!!==================================================================================
   Subroutine Store_Hydrostatic_edgevars()
!    type(element_t), intent(in), dimension(:,:,:)  :: vars
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars

!    real(kind=real_kind), dimension(nx,nx,nx) ::   rhb,prb,rtb
    real(kind=real_kind), dimension(6,nx,nx) ::   rhb_edg,prb_edg,rtb_edg 
    integer :: nlx, nly, nlz
    integer :: ii, jj, kk, e, i, j

   nlx = size(vars,1)
   nly = size(vars,2)
   nlz = size(vars,3)


 !if (itn == 1) then    ! time frozen hydro. balanced variables 
 ! print*, 'Hydro/edge initialization ' 
   do kk = 1, nlz
     do jj = 1, nly
       do ii = 1, nlx

!          Call Volume_to_Surf(vars(ii,jj,kk)%rhb(:,:,:),rhb_edg) 
!          Call Volume_to_Surf(vars(ii,jj,kk)%rtb(:,:,:),rtb_edg) 
!          Call Volume_to_Surf(vars(ii,jj,kk)%prb(:,:,:),prb_edg) 

        do j = 1, nx
        do i = 1, nx
!         do e = 1, 6

          edgevars(ii,jj,kk)%edges(1)%rhb(i,j) = vars(ii,jj,kk)%rhb(i,1,j) 
          edgevars(ii,jj,kk)%edges(2)%rhb(i,j) = vars(ii,jj,kk)%rhb(nx,i,j) 
          edgevars(ii,jj,kk)%edges(3)%rhb(i,j) = vars(ii,jj,kk)%rhb(i,nx,j) 
          edgevars(ii,jj,kk)%edges(4)%rhb(i,j) = vars(ii,jj,kk)%rhb(1,i,j) 
          edgevars(ii,jj,kk)%edges(5)%rhb(i,j) = vars(ii,jj,kk)%rhb(i,j,1) 
          edgevars(ii,jj,kk)%edges(6)%rhb(i,j) = vars(ii,jj,kk)%rhb(i,j,nx)
 
          edgevars(ii,jj,kk)%edges(1)%rtb(i,j) = vars(ii,jj,kk)%rtb(i,1,j) 
          edgevars(ii,jj,kk)%edges(2)%rtb(i,j) = vars(ii,jj,kk)%rtb(nx,i,j) 
          edgevars(ii,jj,kk)%edges(3)%rtb(i,j) = vars(ii,jj,kk)%rtb(i,nx,j) 
          edgevars(ii,jj,kk)%edges(4)%rtb(i,j) = vars(ii,jj,kk)%rtb(1,i,j) 
          edgevars(ii,jj,kk)%edges(5)%rtb(i,j) = vars(ii,jj,kk)%rtb(i,j,1) 
          edgevars(ii,jj,kk)%edges(6)%rtb(i,j) = vars(ii,jj,kk)%rtb(i,j,nx)
 
          edgevars(ii,jj,kk)%edges(1)%prb(i,j) = vars(ii,jj,kk)%prb(i,1,j) 
          edgevars(ii,jj,kk)%edges(2)%prb(i,j) = vars(ii,jj,kk)%prb(nx,i,j) 
          edgevars(ii,jj,kk)%edges(3)%prb(i,j) = vars(ii,jj,kk)%prb(i,nx,j) 
          edgevars(ii,jj,kk)%edges(4)%prb(i,j) = vars(ii,jj,kk)%prb(1,i,j) 
          edgevars(ii,jj,kk)%edges(5)%prb(i,j) = vars(ii,jj,kk)%prb(i,j,1) 
          edgevars(ii,jj,kk)%edges(6)%prb(i,j) = vars(ii,jj,kk)%prb(i,j,nx) 
!         end do 
        end do
        end do
    end do 
    end do 
    end do 



 ! endif 

   End Subroutine Store_Hydrostatic_edgevars

!!==================================================================================
  attributes(global) Subroutine store()
!    type(element_t), intent(in), dimension(:,:,:)  :: vars
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars

!    real(kind=real_kind), dimension(nx,nx,nx) ::   rhb,prb,rtb
    real(kind=real_kind), dimension(6,nx,nx) ::   rhb_edg,prb_edg,rtb_edg 
    integer :: nlx, nly, nlz
    integer :: ii, jj, kk, e, i, j

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)


 !if (itn == 1) then    ! time frozen hydro. balanced variables 
 ! print*, 'Hydro/edge initialization ' 
!   do kk = 1, nlz
!     do jj = 1, nly
!       do ii = 1, nlx

!        do j = 1, nx
!        do i = 1, nx
i = threadIdx%x
j = threadIdx%y
ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z

          edgevars_d(ii,jj,kk)%edges(1)%rhb(i,j) = vars_d(ii,jj,kk)%rhb(i,1,j) 
          edgevars_d(ii,jj,kk)%edges(2)%rhb(i,j) = vars_d(ii,jj,kk)%rhb(nx,i,j) 
          edgevars_d(ii,jj,kk)%edges(3)%rhb(i,j) = vars_d(ii,jj,kk)%rhb(i,nx,j) 
          edgevars_d(ii,jj,kk)%edges(4)%rhb(i,j) = vars_d(ii,jj,kk)%rhb(1,i,j) 
          edgevars_d(ii,jj,kk)%edges(5)%rhb(i,j) = vars_d(ii,jj,kk)%rhb(i,j,1) 
          edgevars_d(ii,jj,kk)%edges(6)%rhb(i,j) = vars_d(ii,jj,kk)%rhb(i,j,nx)
 
          edgevars_d(ii,jj,kk)%edges(1)%rtb(i,j) = vars_d(ii,jj,kk)%rtb(i,1,j) 
          edgevars_d(ii,jj,kk)%edges(2)%rtb(i,j) = vars_d(ii,jj,kk)%rtb(nx,i,j) 
          edgevars_d(ii,jj,kk)%edges(3)%rtb(i,j) = vars_d(ii,jj,kk)%rtb(i,nx,j) 
          edgevars_d(ii,jj,kk)%edges(4)%rtb(i,j) = vars_d(ii,jj,kk)%rtb(1,i,j) 
          edgevars_d(ii,jj,kk)%edges(5)%rtb(i,j) = vars_d(ii,jj,kk)%rtb(i,j,1) 
          edgevars_d(ii,jj,kk)%edges(6)%rtb(i,j) = vars_d(ii,jj,kk)%rtb(i,j,nx)
 
          edgevars_d(ii,jj,kk)%edges(1)%prb(i,j) = vars_d(ii,jj,kk)%prb(i,1,j) 
          edgevars_d(ii,jj,kk)%edges(2)%prb(i,j) = vars_d(ii,jj,kk)%prb(nx,i,j) 
          edgevars_d(ii,jj,kk)%edges(3)%prb(i,j) = vars_d(ii,jj,kk)%prb(i,nx,j) 
          edgevars_d(ii,jj,kk)%edges(4)%prb(i,j) = vars_d(ii,jj,kk)%prb(1,i,j) 
          edgevars_d(ii,jj,kk)%edges(5)%prb(i,j) = vars_d(ii,jj,kk)%prb(i,j,1) 
          edgevars_d(ii,jj,kk)%edges(6)%prb(i,j) = vars_d(ii,jj,kk)%prb(i,j,nx) 
!        end do
!        end do
!    end do 
!    end do 
!    end do 
 

   End Subroutine store

!!==================================================================================
!! Compute NH vars from the evolved state vector 

   Subroutine Recover_State_Vars()
!    type(element_t), intent(inout), dimension(:,:,:)  :: vars
!    type(rkvec_t), intent(in), dimension(:,:,:)  :: rk_temp   
!    real(kind=real_kind), dimension(nx,nx,nx,neq) :: vec
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,the
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,prp,the,rho, rhb,thb,rtb
   integer :: nlx, nly, nlz
   integer :: ii, jj, kk, eq, e, i,j,k 
real(kind=real_kind) :: tmp, rho, prs, spd, the

   nlx = size(vars,1)
   nly = size(vars,2)
   nlz = size(vars,3)


   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T
   ! "communicaed" vars across processors:  [rho,u,v,w,the] 

   do kk = 1, nlz
     do jj = 1, nly
       do ii = 1, nlx
        do k = 1, nx
        do j = 1, nx
        do i = 1, nx
    ! Predicted state vector by RK  
!        vec(:,:,:,:) = rk_stage(ii,jj,kk)%svec(:,:,:,:)

    !Time-frozen variable from strorage 
!        rhb(:,:,:) = vars(ii,jj,kk)%rhb(:,:,:)
!        thb(:,:,:) = vars(ii,jj,kk)%thb(:,:,:) 
!        rtb(:,:,:) = vars(ii,jj,kk)%rtb(:,:,:) 

    !Update NH variables which evolve in time 
            rho = rk_stage(ii,jj,kk)%svec(i,j,k,1)   + vars(ii,jj,kk)%rhb(i,j,k)
!            tmp = 1/rho
       vars(ii,jj,kk)%rho(i,j,k) = rho
       vars(ii,jj,kk)%u(i,j,k) = rk_stage(ii,jj,kk)%svec(i,j,k,2) / rho
       vars(ii,jj,kk)%v(i,j,k) = rk_stage(ii,jj,kk)%svec(i,j,k,3) / rho
       vars(ii,jj,kk)%w(i,j,k) = rk_stage(ii,jj,kk)%svec(i,j,k,4) /rho 

            the = (rk_stage(ii,jj,kk)%svec(i,j,k,5) + vars(ii,jj,kk)%rtb(i,j,k)) / rho
       vars(ii,jj,kk)%the(i,j,k) = the
       vars(ii,jj,kk)%thp(i,j,k) = the - vars(ii,jj,kk)%thb(i,j,k)

!           Call Compute_int_3dvars(rho,the,prs,spd)
        prs = Compute_Pressure(rho, the)

       vars(ii,jj,kk)%prs(i,j,k) = prs
       vars(ii,jj,kk)%spd(i,j,k) = Compute_SoundSpeed(prs, the)
       vars(ii,jj,kk)%prp(i,j,k) = prs - vars(ii,jj,kk)%prb(i,j,k)

    !    if ((ii==8).and.(jj==8).and.(kk==2)) then
    !     print*,  "max speed", minval(spd), maxval(spd) 
    !    endif


        end do
        end do
        end do
       end do 
     end do 
   end do 

  End Subroutine Recover_State_Vars

!!==================================================================================
!! Compute NH vars from the evolved state vector 

 attributes(global) Subroutine recover()
!    type(element_t), intent(inout), dimension(:,:,:)  :: vars
!    type(rkvec_t), intent(in), dimension(:,:,:)  :: rk_temp   
!    real(kind=real_kind), dimension(nx,nx,nx,neq) :: vec
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,the
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,prp,the,rho, rhb,thb,rtb
!   integer :: nlx, nly, nlz
   integer :: ii, jj, kk, eq, e, i,j,k 
real(kind=real_kind) :: tmp, rho, prs, spd, the

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)


   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T
   ! "communicaed" vars across processors:  [rho,u,v,w,the] 
i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
!   do kk = 1, nlz
!     do jj = 1, nly
!       do ii = 1, nlx
!        do k = 1, nx
!        do j = 1, nx
!        do i = 1, nx
    ! Predicted state vector by RK  
!        vec(:,:,:,:) = rk_stage(ii,jj,kk)%svec(:,:,:,:)

    !Time-frozen variable from strorage 
!        rhb(:,:,:) = vars(ii,jj,kk)%rhb(:,:,:)
!        thb(:,:,:) = vars(ii,jj,kk)%thb(:,:,:) 
!        rtb(:,:,:) = vars(ii,jj,kk)%rtb(:,:,:) 

    !Update NH variables which evolve in time 
            rho = rk_stage_d(ii,jj,kk)%svec(i,j,k,1)   + vars_d(ii,jj,kk)%rhb(i,j,k)
!            tmp = 1/rho
       vars_d(ii,jj,kk)%rho(i,j,k) = rho
       vars_d(ii,jj,kk)%u(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,2) / rho
       vars_d(ii,jj,kk)%v(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,3) / rho
       vars_d(ii,jj,kk)%w(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,4) /rho 

            the = (rk_stage_d(ii,jj,kk)%svec(i,j,k,5) + vars_d(ii,jj,kk)%rtb(i,j,k)) / rho
       vars_d(ii,jj,kk)%the(i,j,k) = the
       vars_d(ii,jj,kk)%thp(i,j,k) = the - vars_d(ii,jj,kk)%thb(i,j,k)

!           Call Compute_int_3dvars(rho,the,prs,spd)
!        prs = Compute_Pressure(rho, the)
call Pressure(rho,the,prs)
call SoundSpeed(prs,rho,spd)

       vars_d(ii,jj,kk)%prs(i,j,k) = prs
       vars_d(ii,jj,kk)%spd(i,j,k) = spd
       vars_d(ii,jj,kk)%prp(i,j,k) = prs - vars_d(ii,jj,kk)%prb(i,j,k)

    !    if ((ii==8).and.(jj==8).and.(kk==2)) then
    !     print*,  "max speed", minval(spd), maxval(spd) 
    !    endif


!        end do
!        end do
!        end do
!       end do 
!     end do 
!   end do 

  End Subroutine recover
!!==================================================================================
!! Compute NH vars from the evolved state vector 

 attributes(global) Subroutine recover1()
!    type(element_t), intent(inout), dimension(:,:,:)  :: vars
!    type(rkvec_t), intent(in), dimension(:,:,:)  :: rk_temp   
!    real(kind=real_kind), dimension(nx,nx,nx,neq) :: vec
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,the
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,prp,the,rho, rhb,thb,rtb
!   integer :: nlx, nly, nlz
   integer :: ii, jj, kk, eq, e, i,j,k 
real(kind=real_kind) :: tmp, rho, prs, spd, the

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)


   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T
   ! "communicaed" vars across processors:  [rho,u,v,w,the] 
i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
!   do kk = 1, nlz
!     do jj = 1, nly
!       do ii = 1, nlx
!        do k = 1, nx
!        do j = 1, nx
!        do i = 1, nx
    ! Predicted state vector by RK  
!        vec(:,:,:,:) = rk_stage(ii,jj,kk)%svec(:,:,:,:)

    !Time-frozen variable from strorage 
!        rhb(:,:,:) = vars(ii,jj,kk)%rhb(:,:,:)
!        thb(:,:,:) = vars(ii,jj,kk)%thb(:,:,:) 
!        rtb(:,:,:) = vars(ii,jj,kk)%rtb(:,:,:) 

    !Update NH variables which evolve in time 
            rho = rk_stage_d(ii,jj,kk)%svec(i,j,k,1)   + vars_d(ii,jj,kk)%rhb(i,j,k)
!            tmp = 1/rho
       vars_d(ii,jj,kk)%rho(i,j,k) = rho
       vars_d(ii,jj,kk)%u(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,2) / rho
       vars_d(ii,jj,kk)%v(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,3) / rho
       vars_d(ii,jj,kk)%w(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,4) /rho 

            the = (rk_stage_d(ii,jj,kk)%svec(i,j,k,5) + vars_d(ii,jj,kk)%rtb(i,j,k)) / rho
       vars_d(ii,jj,kk)%the(i,j,k) = the
       vars_d(ii,jj,kk)%thp(i,j,k) = the - vars_d(ii,jj,kk)%thb(i,j,k)

!           Call Compute_int_3dvars(rho,the,prs,spd)
!        prs = Compute_Pressure(rho, the)
call Pressure(rho,the,prs)
call SoundSpeed(prs,rho,spd)

       vars_d(ii,jj,kk)%prs(i,j,k) = prs
       vars_d(ii,jj,kk)%spd(i,j,k) = spd
       vars_d(ii,jj,kk)%prp(i,j,k) = prs - vars_d(ii,jj,kk)%prb(i,j,k)

    !    if ((ii==8).and.(jj==8).and.(kk==2)) then
    !     print*,  "max speed", minval(spd), maxval(spd) 
    !    endif


!        end do
!        end do
!        end do
!       end do 
!     end do 
!   end do 

  End Subroutine recover1
!!==================================================================================
!! Compute NH vars from the evolved state vector 

 attributes(global) Subroutine recover2()
!    type(element_t), intent(inout), dimension(:,:,:)  :: vars
!    type(rkvec_t), intent(in), dimension(:,:,:)  :: rk_temp   
!    real(kind=real_kind), dimension(nx,nx,nx,neq) :: vec
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,the
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,prp,the,rho, rhb,thb,rtb
!   integer :: nlx, nly, nlz
   integer :: ii, jj, kk, eq, e, i,j,k 
real(kind=real_kind) :: tmp, rho, prs, spd, the

!   nlx = size(vars,1)
!   nly = size(vars,2)
!   nlz = size(vars,3)


   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T
   ! "communicaed" vars across processors:  [rho,u,v,w,the] 
i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
!   do kk = 1, nlz
!     do jj = 1, nly
!       do ii = 1, nlx
!        do k = 1, nx
!        do j = 1, nx
!        do i = 1, nx
    ! Predicted state vector by RK  
!        vec(:,:,:,:) = rk_stage(ii,jj,kk)%svec(:,:,:,:)

    !Time-frozen variable from strorage 
!        rhb(:,:,:) = vars(ii,jj,kk)%rhb(:,:,:)
!        thb(:,:,:) = vars(ii,jj,kk)%thb(:,:,:) 
!        rtb(:,:,:) = vars(ii,jj,kk)%rtb(:,:,:) 

    !Update NH variables which evolve in time 
            rho = rk_stage_d(ii,jj,kk)%svec(i,j,k,1)   + vars_d(ii,jj,kk)%rhb(i,j,k)
!            tmp = 1/rho
       vars_d(ii,jj,kk)%rho(i,j,k) = rho
       vars_d(ii,jj,kk)%u(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,2) / rho
       vars_d(ii,jj,kk)%v(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,3) / rho
       vars_d(ii,jj,kk)%w(i,j,k) = rk_stage_d(ii,jj,kk)%svec(i,j,k,4) /rho 

            the = (rk_stage_d(ii,jj,kk)%svec(i,j,k,5) + vars_d(ii,jj,kk)%rtb(i,j,k)) / rho
       vars_d(ii,jj,kk)%the(i,j,k) = the
       vars_d(ii,jj,kk)%thp(i,j,k) = the - vars_d(ii,jj,kk)%thb(i,j,k)

!           Call Compute_int_3dvars(rho,the,prs,spd)
!        prs = Compute_Pressure(rho, the)
call Pressure(rho,the,prs)
call SoundSpeed(prs,rho,spd)

       vars_d(ii,jj,kk)%prs(i,j,k) = prs
       vars_d(ii,jj,kk)%spd(i,j,k) = spd
       vars_d(ii,jj,kk)%prp(i,j,k) = prs - vars_d(ii,jj,kk)%prb(i,j,k)

    !    if ((ii==8).and.(jj==8).and.(kk==2)) then
    !     print*,  "max speed", minval(spd), maxval(spd) 
    !    endif


!        end do
!        end do
!        end do
!       end do 
!     end do 
!   end do 

  End Subroutine recover2

!-------------

End Module Euler_flux3d_mod 

