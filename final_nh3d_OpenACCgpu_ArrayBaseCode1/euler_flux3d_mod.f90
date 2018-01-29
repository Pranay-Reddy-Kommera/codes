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
  use testcases_nh3d_mod, only : Compute_Pressure, Compute_SoundSpeed
  use openacc

  implicit none
  private
  public :: Make_Euler3d_Flux , Compute_Ext_Vars, nh_flux_maker 
  public :: Recover_State_Vars, Store_Hydrostatic_edgevars

contains

!!============================================================================================
   
   Subroutine nh_flux_maker()
   implicit none

    real(kind=real_kind), dimension(neq) :: svec
!    real(kind=real_kind) :: fcori
!    real(kind=real_kind), dimension(nx,nx,nx) :: fcori
    real(kind=real_kind) :: u,v,w,spd,prp,the,rho,fcori
!    real(kind=real_kind), dimension(nx,nx,nx) :: u,v,w,spd,prp,the,rho
!    real(kind=real_kind), dimension(6,nx,nx) :: spd_int,prp_ext
!    real(kind=real_kind), dimension(nx,nx) ::  prp_e,the_e,rho_e,prs_e,spd_e,
!    rtp_e 
    integer :: eq, e, i,j,k,ie,je,ke,nlx,nly,nlz

   nlx = size(vars_rho,4)
   nly = size(vars_rho,5)
   nlz = size(vars_rho,6)

   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T
    ! vec(:,:,:,:) = vars%vec(:,:,:,:) 

!$acc parallel num_workers(8) vector_length(16)
!$acc loop gang worker vector collapse(6) private(u,v,w,rho,the,prp,spd,fcori,svec)
         do ke = 1, nlz
        do je = 1, nly
        do ie = 1, nlx
       do k = 1, nx
        do j = 1, nx
        do i = 1, nx

     u = vars_u(i,j,k,ie,je,ke) 
     v = vars_v(i,j,k,ie,je,ke) 
     w = vars_w(i,j,k,ie,je,ke) 

     rho = vars_rho(i,j,k,ie,je,ke) 
     the = vars_the(i,j,k,ie,je,ke) 
     prp = vars_prp(i,j,k,ie,je,ke)  
     spd = vars_spd(i,j,k,ie,je,ke) 

     !Full pressure in the discretization 
     !prs(:,:,:) =  vars%prb(:,:,:) 
     !Call Horizontal_strong_grads(prs,prs_x,prs_y) 

     fcori = vars_fcori(i,j,k,ie,je,ke) 

!Construct state vecror from the evolved vars 
              svec(1) = rho - vars_rhb(i,j,k,ie,je,ke)
              svec(2) = rho * u
              svec(3) = rho * v
              svec(4) = rho * w
              svec(5) = rho * the - vars_rtb(i,j,k,ie,je,ke)

  ! Flux x-direction 
    vars_flx(i,j,k,1,ie,je,ke) = svec(2)
    vars_flx(i,j,k,2,ie,je,ke) = svec(2) * u + prp
    vars_flx(i,j,k,3,ie,je,ke) = svec(2) * v
    vars_flx(i,j,k,4,ie,je,ke) = svec(2) * w
    vars_flx(i,j,k,5,ie,je,ke) = svec(2) * the

  ! Flux y-direction 
    vars_fly(i,j,k,1,ie,je,ke) = svec(3)
    vars_fly(i,j,k,2,ie,je,ke) = svec(3) * u
    vars_fly(i,j,k,3,ie,je,ke) = svec(3) * v + prp
    vars_fly(i,j,k,4,ie,je,ke) = svec(3) * w
    vars_fly(i,j,k,5,ie,je,ke) = svec(3) * the

  ! Flux z-direction 
    vars_flz(i,j,k,1,ie,je,ke) = svec(4)
    vars_flz(i,j,k,2,ie,je,ke) = svec(4) * u
    vars_flz(i,j,k,3,ie,je,ke) = svec(4) * v
    vars_flz(i,j,k,4,ie,je,ke) = svec(4) * w + prp
    vars_flz(i,j,k,5,ie,je,ke) = svec(4) * the

  ! Source term         
    vars_src(i,j,k,1,ie,je,ke) = 0.0D0
    vars_src(i,j,k,2,ie,je,ke) =  fcori*v !- prs_x(:,:,:) 
    vars_src(i,j,k,3,ie,je,ke) = -fcori*u !- prs_y(:,:,:) 
    vars_src(i,j,k,4,ie,je,ke) = -svec(1) * grv
    vars_src(i,j,k,5,ie,je,ke) = 0.0D0

        end do
        end do
        end do

        end do
        end do
        end do
!$acc end parallel

end subroutine nh_flux_maker

!!============================================================================================
  Subroutine Make_Euler3d_Flux()
  implicit none
!    type(element_t), intent(inout) :: vars
!    type(edge_t), intent(inout), dimension(6) :: edges
!    integer, intent(in) :: ke, ie,je  

    integer :: eq, e, i,j,k,ie,je,ke,nlx,nly,nlz

   nlx = size(vars_rho,4)
   nly = size(vars_rho,5)
   nlz = size(vars_rho,6)

!$acc parallel num_workers(8) vector_length(16)
!$acc loop gang worker vector collapse(5) private(eq)
        do ke = 1, nlz
        do je = 1, nly
        do ie = 1, nlx
        do j = 1, nx
        do i = 1, nx

    edgevars_edges_spd_int(i,j,1,ie,je,ke) = vars_spd(i,1,j,ie,je,ke)
    edgevars_edges_spd_int(i,j,2,ie,je,ke) = vars_spd(nx,i,j,ie,je,ke)
    edgevars_edges_spd_int(i,j,3,ie,je,ke) = vars_spd(i,nx,j,ie,je,ke)
    edgevars_edges_spd_int(i,j,4,ie,je,ke) = vars_spd(1,i,j,ie,je,ke)
    edgevars_edges_spd_int(i,j,5,ie,je,ke) = vars_spd(i,j,1,ie,je,ke)
    edgevars_edges_spd_int(i,j,6,ie,je,ke) = vars_spd(i,j,nx,ie,je,ke)


   do eq = 1, neq
    edgevars_edges_normflux_int(i,j,eq,1,ie,je,ke) = vars_fly(i,1,j,eq,ie,je,ke)
    edgevars_edges_normflux_int(i,j,eq,2,ie,je,ke) = vars_flx(nx,i,j,eq,ie,je,ke)
    edgevars_edges_normflux_int(i,j,eq,3,ie,je,ke) = vars_fly(i,nx,j,eq,ie,je,ke)
    edgevars_edges_normflux_int(i,j,eq,4,ie,je,ke) = vars_flx(1,i,j,eq,ie,je,ke)
    edgevars_edges_normflux_int(i,j,eq,5,ie,je,ke) = vars_flz(i,j,1,eq,ie,je,ke)
    edgevars_edges_normflux_int(i,j,eq,6,ie,je,ke) = vars_flz(i,j,nx,eq,ie,je,ke)
   end do

        end do
        end do
        end do
        end do
        end do
!$acc end parallel

   Call Make_ext_normalflux()

   Call Apply_Noflux_TopbotBC()

  End Subroutine Make_Euler3d_Flux

!!============================================================================================

Subroutine Apply_Noflux_TopbotBC()
  implicit none

    integer :: ie,je, ke,nlx,nly,nlz,i,j

   nlx = size(vars_rho,4)
   nly = size(vars_rho,5)
   nlz = size(vars_rho,6)

!$acc parallel num_workers(8) vector_length(16)
!$acc loop gang worker vector collapse(5)
        do ke = 1, nlz
        do je = 1, nly
        do ie = 1, nlx
        do j = 1, nx
        do i = 1, nx

     if (ke == 1) then
          edgevars_edges_vec_ext(i,j,:,5,ie,je,ke) = edgevars_edges_vec_int(i,j,:,5,ie,je,ke)
          edgevars_edges_vec_ext(i,j,4,5,ie,je,ke) = -edgevars_edges_vec_int(i,j,4,5,ie,je,ke)

          edgevars_edges_normflux_ext(i,j,:,5,ie,je,ke) = -edgevars_edges_normflux_int(i,j,:,5,ie,je,ke)
          edgevars_edges_normflux_ext(i,j,4,5,ie,je,ke) =  edgevars_edges_normflux_int(i,j,4,5,ie,je,ke)
          edgevars_edges_spd_ext(i,j,5,ie,je,ke) =  edgevars_edges_spd_int(i,j,5,ie,je,ke)
     endif

     if (ke == nez) then
          edgevars_edges_vec_ext(i,j,:,6,ie,je,ke) = edgevars_edges_vec_int(i,j,:,6,ie,je,ke)
          edgevars_edges_vec_ext(i,j,4,6,ie,je,ke) = -edgevars_edges_vec_int(i,j,4,6,ie,je,ke)

          edgevars_edges_normflux_ext(i,j,:,6,ie,je,ke) = -edgevars_edges_normflux_int(i,j,:,6,ie,je,ke)
          edgevars_edges_normflux_ext(i,j,4,6,ie,je,ke) =  edgevars_edges_normflux_int(i,j,4,6,ie,je,ke)
          edgevars_edges_spd_ext(i,j,6,ie,je,ke) = edgevars_edges_spd_int(i,j,6,ie,je,ke)
     endif

end do
end do
end do
end do
end do
!$acc end parallel

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

!! ============================================================================================

Subroutine Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)
!$acc routine
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

Subroutine  Make_ext_normalflux()
    implicit none

    real(kind=real_kind) :: vec_ext(neq)
    real(kind=real_kind) :: u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb
    integer :: face,i,j,ie,je,ke,nlx,nly,nlz

   nlx = size(vars_rho,4)
   nly = size(vars_rho,5)
   nlz = size(vars_rho,6)

!$acc parallel num_workers(16) vector_length(16)
!$acc loop gang worker vector collapse(5) private(face,vec_ext,rhb,rtb,prb,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)
do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx
do j = 1, nx
do i = 1, nx

    !South face (y-flux)  
    !------------------
    do  face = 1,3, 2
        vec_ext(:) = edgevars_edges_vec_ext(i,j,:,face,ie,je,ke)
        rhb = edgevars_edges_rhb(i,j,face,ie,je,ke)
        rtb = edgevars_edges_rtb(i,j,face,ie,je,ke)
        prb = edgevars_edges_prb(i,j,face,ie,je,ke)
     Call Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars_edges_normflux_ext(i,j,1,face,ie,je,ke) =  vec_ext(3)    !rho*v 
     edgevars_edges_normflux_ext(i,j,2,face,ie,je,ke) =  vec_ext(3) * u_ext
     edgevars_edges_normflux_ext(i,j,3,face,ie,je,ke) =  vec_ext(3) * v_ext + prp_ext
     edgevars_edges_normflux_ext(i,j,4,face,ie,je,ke) =  vec_ext(3) * w_ext
     edgevars_edges_normflux_ext(i,j,5,face,ie,je,ke) =  vec_ext(3) * the_ext
     edgevars_edges_spd_ext(i,j,face,ie,je,ke) =  spd_ext
    end do

    !North (y-flux)  
    !------------------
   !  face = 3 
   !     vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
   !  Call
   !  Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)
!
!     edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,3)
!     edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,3) * u_ext(:,:)
!     edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,3) * v_ext(:,:) +
!     prp_ext(:,:)
!     edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,3) * w_ext(:,:)
!     edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,3) * the_ext(:,:)
!     edge(face)%spd_ext(:,:) =  spd_ext(:,:)

!East  face (x-flux)  
    !------------------
    do  face = 2,4, 2
        vec_ext(:) = edgevars_edges_vec_ext(i,j,:,face,ie,je,ke)
        rhb = edgevars_edges_rhb(i,j,face,ie,je,ke)
        rtb = edgevars_edges_rtb(i,j,face,ie,je,ke)
        prb = edgevars_edges_prb(i,j,face,ie,je,ke)
     Call Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars_edges_normflux_ext(i,j,1,face,ie,je,ke) =  vec_ext(2)    !rho*u 
     edgevars_edges_normflux_ext(i,j,2,face,ie,je,ke) =  vec_ext(2) * u_ext + prp_ext
     edgevars_edges_normflux_ext(i,j,3,face,ie,je,ke) =  vec_ext(2) * v_ext
     edgevars_edges_normflux_ext(i,j,4,face,ie,je,ke) =  vec_ext(2) * w_ext
     edgevars_edges_normflux_ext(i,j,5,face,ie,je,ke) =  vec_ext(2) * the_ext
     edgevars_edges_spd_ext(i,j,face,ie,je,ke) =  spd_ext
    end do

!   !West (x-flux)  
!   !------------------
!    face = 4
!       vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
!    Call
!    Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

!    edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,2)
!    edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,2) * u_ext(:,:) +
!    prp_ext(:,:)
!    edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,2) * v_ext(:,:)
!    edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,2) * w_ext(:,:)
!    edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,2) * the_ext(:,:)
!    edge(face)%spd_ext(:,:) =  spd_ext(:,:)

    !Bottom (z-flux)  
    !------------------
    do face = 5, 6
        vec_ext(:) = edgevars_edges_vec_ext(i,j,:,face,ie,je,ke)
        rhb = edgevars_edges_rhb(i,j,face,ie,je,ke)
        rtb = edgevars_edges_rtb(i,j,face,ie,je,ke)
        prb = edgevars_edges_prb(i,j,face,ie,je,ke)
     Call Extract_ext_Vars(face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext,rhb,rtb,prb)

     edgevars_edges_normflux_ext(i,j,1,face,ie,je,ke) =  vec_ext(4)    !rho*w 
     edgevars_edges_normflux_ext(i,j,2,face,ie,je,ke) =  vec_ext(4) * u_ext
     edgevars_edges_normflux_ext(i,j,3,face,ie,je,ke) =  vec_ext(4) * v_ext
     edgevars_edges_normflux_ext(i,j,4,face,ie,je,ke) =  vec_ext(4) * w_ext + prp_ext
     edgevars_edges_normflux_ext(i,j,5,face,ie,je,ke) =  vec_ext(4) * the_ext
     edgevars_edges_spd_ext(i,j,face,ie,je,ke) =  spd_ext
    end do

    !Top (z-flux)  
    !------------------
!    face = 6
!       vec_ext(:,:,:) = edge(face)%vec_ext(:,:,:)
!    Call
!    Extract_ext_Vars(edge,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

!    edge(face)%normflux_ext(:,:,1) =  vec_ext(:,:,4)
!    edge(face)%normflux_ext(:,:,2) =  vec_ext(:,:,4) * u_ext(:,:) 
!    edge(face)%normflux_ext(:,:,3) =  vec_ext(:,:,4) * v_ext(:,:)
!    edge(face)%normflux_ext(:,:,4) =  vec_ext(:,:,4) * w_ext(:,:) +
!    prp_ext(:,:)
!    edge(face)%normflux_ext(:,:,5) =  vec_ext(:,:,4) * the_ext(:,:)
!    edge(face)%spd_ext(:,:) =  spd_ext(:,:)
!   
end do
end do
end do
end do
end do
!$acc end parallel

 End Subroutine  Make_ext_normalflux

!=======================================================================================================!
   Subroutine Compute_Ext_Vars(rho_ext,the_ext, prs_ext,spd_ext) 
!$acc routine
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
!$acc routine
    Implicit None

    real(kind=real_kind), dimension(nx,nx,nx), intent(in) ::   rho,the
    real(kind=real_kind), dimension(nx,nx,nx), intent(out) ::  prs,spd
    real(kind=real_kind) ::  ro,th,pr
    integer :: i,j,k

           do k = 1, nx
            do j = 1, nx
              do i = 1, nx
               ro = rho(i,j,k)
               th = the(i,j,k)

               pr = Compute_Pressure(ro,th)
               prs(i,j,k) = pr

               spd(i,j,k) = Compute_SoundSpeed(pr,ro)

              end do
            end do
           end do

   End Subroutine Compute_int_3dvars 

!=======================================================================================================!
   Subroutine Volume_to_Surf(psi,psi_edg)
!$acc routine
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

   nlx = size(vars_rho,4)
   nly = size(vars_rho,5)
   nlz = size(vars_rho,6)


 !if (itn == 1) then    ! time frozen hydro. balanced variables 
 ! print*, 'Hydro/edge initialization ' 
!$acc parallel num_workers(4) vector_length(16)
!$acc loop gang worker vector collapse(5)
   do kk = 1, nlz
     do jj = 1, nly
       do ii = 1, nlx
!          rhb(:,:,:) = vars_rhb(:,:,:,ii,jj,kk) 
!          rtb(:,:,:) = vars_rtb(:,:,:,ii,jj,kk) 
!          prb(:,:,:) = vars_prb(:,:,:,ii,jj,kk) 

!          Call Volume_to_Surf(rhb,rhb_edg) 
!          Call Volume_to_Surf(rtb,rtb_edg) 
!          Call Volume_to_Surf(prb,prb_edg) 

        do j = 1, nx
        do i = 1, nx
!         do e = 1, 6
          edgevars_edges_rhb(i,j,1,ii,jj,kk) = vars_rhb(i,1,j,ii,jj,kk) 
          edgevars_edges_rhb(i,j,2,ii,jj,kk) = vars_rhb(nx,i,j,ii,jj,kk) 
          edgevars_edges_rhb(i,j,3,ii,jj,kk) = vars_rhb(i,nx,j,ii,jj,kk) 
          edgevars_edges_rhb(i,j,4,ii,jj,kk) = vars_rhb(1,i,j,ii,jj,kk) 
          edgevars_edges_rhb(i,j,5,ii,jj,kk) = vars_rhb(i,j,1,ii,jj,kk) 
          edgevars_edges_rhb(i,j,6,ii,jj,kk) = vars_rhb(i,j,nx,ii,jj,kk) 

          edgevars_edges_rtb(i,j,1,ii,jj,kk) = vars_rtb(i,1,j,ii,jj,kk) 
          edgevars_edges_rtb(i,j,2,ii,jj,kk) = vars_rtb(nx,i,j,ii,jj,kk) 
          edgevars_edges_rtb(i,j,3,ii,jj,kk) = vars_rtb(i,nx,j,ii,jj,kk) 
          edgevars_edges_rtb(i,j,4,ii,jj,kk) = vars_rtb(1,i,j,ii,jj,kk) 
          edgevars_edges_rtb(i,j,5,ii,jj,kk) = vars_rtb(i,j,1,ii,jj,kk) 
          edgevars_edges_rtb(i,j,6,ii,jj,kk) = vars_rtb(i,j,nx,ii,jj,kk) 

          edgevars_edges_prb(i,j,1,ii,jj,kk) = vars_prb(i,1,j,ii,jj,kk) 
          edgevars_edges_prb(i,j,2,ii,jj,kk) = vars_prb(nx,i,j,ii,jj,kk) 
          edgevars_edges_prb(i,j,3,ii,jj,kk) = vars_prb(i,nx,j,ii,jj,kk) 
          edgevars_edges_prb(i,j,4,ii,jj,kk) = vars_prb(1,i,j,ii,jj,kk) 
          edgevars_edges_prb(i,j,5,ii,jj,kk) = vars_prb(i,j,1,ii,jj,kk) 
          edgevars_edges_prb(i,j,6,ii,jj,kk) = vars_prb(i,j,nx,ii,jj,kk) 
!         end do 
        end do
        end do
    end do 
    end do 
    end do 
!$acc end parallel

 ! endif 

   End Subroutine Store_Hydrostatic_edgevars

!!==================================================================================
!! Compute NH vars from the evolved state vector 

   Subroutine Recover_State_Vars()
!    type(element_t), intent(inout), dimension(:,:,:)  :: vars
!    type(rkvec_t), intent(in), dimension(:,:,:)  :: rk_temp   
!    real(kind=real_kind), dimension(nx,nx,nx,neq) :: vec
!    real(kind=real_kind), dimension(nx,nx,nx) ::  prs,spd,prp,the,rho, rhb,thb,rtb
   integer :: nlx, nly, nlz
   integer :: ii, jj, kk, eq, e, i,j,k 
   real(kind=real_kind) :: tmp, rho, prs, spd, the

   nlx = size(vars_rho,4)
   nly = size(vars_rho,5)
   nlz = size(vars_rho,6)


   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T
   ! "communicaed" vars across processors:  [rho,u,v,w,the] 

!$acc parallel num_workers(16) vector_length(32)
!$acc loop gang worker vector collapse(6) private(prs, rho, the)
   do kk = 1, nlz
     do jj = 1, nly
       do ii = 1, nlx
        do k = 1, nx
        do j = 1, nx
        do i = 1, nx

    ! Predicted state vector by RK  
!        vec(:,:,:,:) = rk_stage_svec(:,:,:,:,ii,jj,kk)

    !Time-frozen variable from strorage 
!        rhb(:,:,:) = vars_rhb(:,:,:,ii,jj,kk)
!        thb(:,:,:) = vars_thb(:,:,:,ii,jj,kk) 
!        rtb(:,:,:) = vars_rtb(:,:,:,ii,jj,kk) 

    !Update NH variables which evolve in time 
            rho = rk_stage_svec(i,j,k,1,ii,jj,kk)   + vars_rhb(i,j,k,ii,jj,kk)
       vars_rho(i,j,k,ii,jj,kk) = rho
       vars_u(i,j,k,ii,jj,kk) = rk_stage_svec(i,j,k,2,ii,jj,kk) / rho
       vars_v(i,j,k,ii,jj,kk) = rk_stage_svec(i,j,k,3,ii,jj,kk) / rho
       vars_w(i,j,k,ii,jj,kk) = rk_stage_svec(i,j,k,4,ii,jj,kk) / rho

            the = (rk_stage_svec(i,j,k,5,ii,jj,kk) + vars_rtb(i,j,k,ii,jj,kk)) / rho
       vars_the(i,j,k,ii,jj,kk) = the
       vars_thp(i,j,k,ii,jj,kk) = the - vars_thb(i,j,k,ii,jj,kk)

!           Call Compute_int_3dvars(rho,the,prs,spd)
        prs = Compute_Pressure(rho, the)

       vars_prs(i,j,k,ii,jj,kk) = prs
       vars_spd(i,j,k,ii,jj,kk) = Compute_SoundSpeed(prs, the)
       vars_prp(i,j,k,ii,jj,kk) = prs - vars_prb(i,j,k,ii,jj,kk)

    !    if ((ii==8).and.(jj==8).and.(kk==2)) then
    !     print*,  "max speed", minval(spd), maxval(spd) 
    !    endif

        end do
        end do  
        end do
       end do 
     end do 
   end do 
!$acc end parallel

  End Subroutine Recover_State_Vars

!-------------

End Module Euler_flux3d_mod 

