!!-------------------------------------------------------------------
!! Ram Nair, 07/2016, 11/2016 
!! Implements the weak form DG method for evaluating the
!! semi-discrete RHS of the advection equation.
!!-------------------------------------------------------------------

module vertical_split_mod
  use basic_mod
  use physical_const_mod, only: grv 
  use element_mod
  use gauss_quadrature_mod
  use dg3d_rhs_mod, only : Horizontal_strong_grads
  use euler_flux3d_mod, only: Compute_Ext_Vars 
  use openacc

  implicit none
  private
   public :: Compute_SplitDG_rhs, Compute_Column_rhs
  private :: Vertical_DGEuler_solver, Vertical_LXF_Flux
  private :: DGNH_Horizontal_rhs, Extract_ext_Vars, Compute_2DFluxJacobian
  private :: LxF_NH2d_Flux

  Integer, private, parameter :: ibot=1, itop=2

contains

!!===========================================================================
  Subroutine Compute_Column_rhs()

!    type(element_t), dimension(:,:,:), intent(inout) :: vars
!    type(column_t) :: velem
    integer :: i, j, k, ie, je, ke
    integer :: nlx, nly, nlz

    nlx = size(vars_rho,4)
    nly = size(vars_rho,5)
    nlz = size(vars_rho,6)

!$acc parallel num_workers(32) vector_length(32)
!$acc loop gang worker vector collapse(4)
    do je = 1, nly
    do ie = 1, nlx

        do j = 1, nx
        do i = 1, nx

            do ke = 1, nlz
              do k = 1, nx
                velem_u(k,ke)  = vars_u(i,j,k,ie,je,ke)
                velem_v(k,ke)  = vars_v(i,j,k,ie,je,ke)
                velem_w(k,ke)  = vars_w(i,j,k,ie,je,ke) 
                velem_rho(k,ke)= vars_rho(i,j,k,ie,je,ke)
                velem_the(k,ke)= vars_the(i,j,k,ie,je,ke)
                velem_prp(k,ke)= vars_prp(i,j,k,ie,je,ke)
                velem_spd(k,ke)= vars_spd(i,j,k,ie,je,ke)
                velem_vec(k,ke,:)= vars_vec(i,j,k,:,ie,je,ke)
              end do
            enddo

         !DG solver for 1d system of Euler equations in the vertical 
            Call Vertical_DGEuler_solver()

            do ke = 1, nlz
              do k = 1, nx
                vars_zsrc(i,j,k,:,ie,je,ke) = velem_rhs(k,ke,:) 
              end do
            enddo

        end do
        end do

    end do
    end do
!$acc end parallel

  End Subroutine Compute_Column_rhs

 !!===========================================================================
  Subroutine Vertical_DGEuler_solver()
!$acc routine
!    type(column_t), intent(inout)  :: velem

    real (kind=real_kind), dimension(nx,nez,neq) :: flz, svec, nh_src, num_flux  
    real (kind=real_kind), dimension(2,0:nez+1) :: spd_ext,w_ext,spd_int,w_int
    real (kind=real_kind), dimension(2,0:nez+1,neq) :: svec_int, svec_ext, flz_ext, flz_int
    real (kind=real_kind), dimension(2,nez) :: fj_max, fj_int, fj_ext 
    real (kind=real_kind), dimension(nx) :: flux, grad, vimx
    real (kind=real_kind), dimension(nx,nx) :: vder
    real (kind=real_kind) :: rhow, sm
    integer :: i, k, eq, ke

    !! Note:  ibot ==> left edge & itop ==> right edge of 1d vertical  element in a FV sense


!vertical Euler fluxes 
      w_int(:,:) = 0.0D0 
    spd_int(:,:) = 0.0D0 

    do k=1,nx
      vimx(k) = 2.0D0 / (delz*gllw(k))

      ! derivative matrix with gauss-wts
      do i=1, nx
        vder(i,k) =  der(i,k) * gllw(i)
      end do
    end do

     do ke=1, nez 
     do  k=1, nx
           rhow =  velem_rho(k,ke) * velem_w(k,ke)
       flz(k,ke,1) = rhow 
       flz(k,ke,2) = rhow * velem_u(k,ke) 
       flz(k,ke,3) = rhow * velem_v(k,ke) 
       flz(k,ke,4) = rhow * velem_w(k,ke)  + velem_prp(k,ke) 
       flz(k,ke,5) = rhow * velem_the(k,ke)

      svec(k,ke,:) = velem_vec(k,ke,:)
      nh_src(k,ke,:) = 0.0D0
     enddo
      spd_int(itop,ke) = velem_spd(nx,ke)
      spd_int(ibot,ke) = velem_spd(1,ke)
       w_int(itop,ke) = velem_w(nx,ke)
       w_int(ibot,ke) = velem_w(1,ke)
     enddo

    nh_src(:,:,4) = - velem_vec(:,:,1) * grv 

   do eq = 1, neq 
     do ke=1, nez 
        svec_int(ibot,ke,eq) = svec(1,ke,eq)
        svec_int(itop,ke,eq) = svec(nx,ke,eq)
         flz_int(ibot,ke,eq) = flz(1,ke,eq)
         flz_int(itop,ke,eq) = flz(nx,ke,eq)
     end do
    end do

    do ke=1,nez
     svec_ext(ibot,ke,:) = svec_int(itop,ke-1,:)
     svec_ext(itop,ke,:) = svec_int(ibot,ke+1,:)
      flz_ext(ibot,ke,:) = flz_int(itop,ke-1,:)
      flz_ext(itop,ke,:) = flz_int(ibot,ke+1,:)

      spd_ext(ibot,ke) = spd_int(itop,ke-1)
      spd_ext(itop,ke) = spd_int(ibot,ke+1)
       w_ext(ibot,ke) =  w_int(itop,ke-1)
       w_ext(itop,ke) =  w_int(ibot,ke+1)
    enddo

 ! No flux BC at top/bottom boundaries 
        svec_ext(itop,nez,:) =  svec_int(itop,nez,:)
        svec_ext(itop,nez,4) = -svec_int(itop,nez,4)

         flz_ext(itop,nez,:) = -flz_int(itop,nez,:)
         flz_ext(itop,nez,4) =  flz_int(itop,nez,4)
         spd_ext(itop,nez) =  spd_int(itop,nez)

        svec_ext(ibot,1,:) =  svec_int(ibot,1,:)
        svec_ext(ibot,1,4) = -svec_int(ibot,1,4)

         flz_ext(ibot,1,:) = -flz_int(ibot,1,:)
         flz_ext(ibot,1,4) =  flz_int(ibot,1,4)
         spd_ext(ibot,1) =  spd_int(ibot,1)

   ! fluxJacobian for LxF 

    do ke=1,nez
     fj_int(:,ke) = spd_int(:,ke)  + abs(w_int(:,ke))
     fj_ext(:,ke) = spd_ext(:,ke)  + abs(w_ext(:,ke)) 
     fj_max(itop,ke) = max(fj_int(itop,ke), fj_ext(itop,ke))
     fj_max(ibot,ke) = max(fj_int(ibot,ke), fj_ext(ibot,ke))
    end do

    num_flux(:,:,:) = Vertical_LXF_Flux(fj_max,svec_int,svec_ext,flz_int,flz_ext)

    ! DG discretization

    do eq = 1, neq 
    do ke=1,nez
      flux(:) = flz(:,ke,eq)
     !Gradient weak (DG)
      do k=1,nx
        sm = 0.0D0
        do i=1,nx
          sm = sm + flux(i) * vder(i,k)
        enddo
        grad(k) = sm
      enddo

      do k=1,nx
        velem_rhs(k,ke,eq) = (grad(k) -  num_flux(k,ke,eq)) * vimx(k) + nh_src(k,ke,eq) 
      enddo

    enddo
    enddo


  End Subroutine Vertical_DGEuler_solver

!=======================================================================================!
   Function Vertical_LXF_Flux(fj_max,vec_int,vec_ext,flz_int,flz_ext) result(num_flux)
!$acc routine
    Implicit None
    real (kind=real_kind), dimension(2,0:nez+1,5) :: vec_int, vec_ext, flz_int, flz_ext
    real (kind=real_kind), intent(in) :: fj_max(2,nez)
    real (kind=real_kind), dimension(nx) :: c_itop, c_ibot
    real (kind=real_kind) :: flz_up,flz_dn, vec_up,vec_dn , lxf_ibot, lxf_itop
    real (kind=real_kind) :: num_flux(nx,nez,neq)
    integer :: k,eq,ke

      c_itop(:) = 0.0D0
      c_ibot(:) = 0.0D0

      c_itop(nx) = 1.0D0
      c_ibot(1) = 1.0D0


      do eq=1, neq 
       do ke=1, nez
           flz_up = flz_int(ibot,ke,eq)
           flz_dn = flz_ext(ibot,ke,eq)
           vec_up = vec_int(ibot,ke,eq)
           vec_dn = vec_ext(ibot,ke,eq)
           lxf_ibot = 0.5D0 * ((flz_up + flz_dn) - fj_max(ibot,ke) *(vec_up - vec_dn))

           flz_up = flz_ext(itop,ke,eq)
           flz_dn = flz_int(itop,ke,eq)
           vec_up = vec_ext(itop,ke,eq)
           vec_dn = vec_int(itop,ke,eq)
           lxf_itop = 0.5D0 * ((flz_up + flz_dn) - fj_max(itop,ke) *(vec_up - vec_dn))

         do k = 1, nx 
          !num_flux(k,ke,eq) =   lxf_itop - lxf_ibot
           num_flux(k,ke,eq) =  c_itop(k) * lxf_itop - c_ibot(k) * lxf_ibot
         end do
       end do
      end do

   End Function Vertical_LXF_Flux
                                        
!==============================================================================!
  Subroutine Compute_SplitDG_rhs(ie,je,ke,rhs_h)

!    type(element_t), intent(inout) :: vars
!    type(edge_t), dimension(6), intent(inout) :: edgevars
integer, intent(in) :: ie, je, ke
    real(kind=real_kind), intent(out), dimension(nx,nx,nx,neq) :: rhs_h 

    real(kind=real_kind), dimension(nx,nx,nx,neq) :: zsrc
    integer :: i, j, col

    do col = 1, nx
      do j = 1, nx
        do i = 1, nx
          zsrc(i,j,col,:) = vars_zsrc(i,j,col,:,ie,je,ke)
        end do
      end do
    end do

    Call DG2d_Horizontal_solver(ie,je,ke,zsrc,rhs_h)

  End Subroutine Compute_SplitDG_rhs

 !!===========================================================================
  Subroutine DG2d_Horizontal_solver(ie,je,ke,z_src,rhs_h)

!    type(element_t), intent(inout) :: vars
!    type(edge_t), intent(inout), dimension(6) :: edges
integer, intent(in) :: ie, je, ke
    real(kind=real_kind), intent(in), dimension(nx,nx,nx,neq) :: z_src
    real(kind=real_kind), intent(out), dimension(nx,nx,nx,neq) :: rhs_h

    real(kind=real_kind), dimension(nx,nx,nx,neq) :: svec 
    real(kind=real_kind), dimension(nx,nx,nx) :: u,v,w, prs,spd,prp,the,rho,rtb , rtp
    real(kind=real_kind), dimension(nx,nx,nx) :: fcori,prs_x,prs_y
    real(kind=real_kind), dimension(nx,nx) ::  prp_e,the_e,rho_e,prs_e,spd_e, rtp_e
    real(kind=real_kind), dimension(nx,nx,neq)  :: vec_ext
    real(kind=real_kind), dimension(nx,nx) :: u_ext, v_ext,w_ext, the_ext,prp_ext,spd_ext


    real(kind=real_kind), dimension(nx,nx,nx) :: flx, fly, grads
    real(kind=real_kind), dimension(4,nx,nx) :: numflux
    real(kind=real_kind) :: dflx, dfly , dx2, dy2, maxspd
    Integer :: col, face, i,j, n, eq, e



   ! State vector: vec() => [rho', rho*u, rho*v, rho*w, (rho*theta)']^T

     u(:,:,:) = vars_u(:,:,:,ie,je,ke)
     v(:,:,:) = vars_v(:,:,:,ie,je,ke)
     w(:,:,:) = vars_w(:,:,:,ie,je,ke)

     rho(:,:,:) = vars_rho(:,:,:,ie,je,ke)
     the(:,:,:) = vars_the(:,:,:,ie,je,ke)
     prp(:,:,:) = vars_prp(:,:,:,ie,je,ke)
     spd(:,:,:) = vars_spd(:,:,:,ie,je,ke)

   !full pressure in the disretization 
   ! prs(:,:,:) = vars%prb(:,:,:)

   !Call Horizontal_strong_grads(prs,prs_x,prs_y)

     fcori(:,:,:) = vars_fcori(:,:,:,ie,je,ke)

  !Construct state vecror from the evolved vars 
              svec(:,:,:,1) = rho(:,:,:) - vars_rhb(:,:,:,ie,je,ke)
              svec(:,:,:,2) = rho(:,:,:) * u(:,:,:)
              svec(:,:,:,3) = rho(:,:,:) * v(:,:,:)
              svec(:,:,:,4) = rho(:,:,:) * w(:,:,:)
              svec(:,:,:,5) = rho(:,:,:) * the(:,:,:) - vars_rtb(:,:,:,ie,je,ke)

  ! Flux x-direction 
    vars_flx(:,:,:,1,ie,je,ke) = svec(:,:,:,2)
    vars_flx(:,:,:,2,ie,je,ke) = svec(:,:,:,2) * u(:,:,:) + prp(:,:,:)
    vars_flx(:,:,:,3,ie,je,ke) = svec(:,:,:,2) * v(:,:,:)
    vars_flx(:,:,:,4,ie,je,ke) = svec(:,:,:,2) * w(:,:,:)
    vars_flx(:,:,:,5,ie,je,ke) = svec(:,:,:,2) * the(:,:,:)

  ! Flux y-direction 
    vars_fly(:,:,:,1,ie,je,ke) = svec(:,:,:,3)
    vars_fly(:,:,:,2,ie,je,ke) = svec(:,:,:,3) * u(:,:,:)
    vars_fly(:,:,:,3,ie,je,ke) = svec(:,:,:,3) * v(:,:,:) + prp(:,:,:)
    vars_fly(:,:,:,4,ie,je,ke) = svec(:,:,:,3) * w(:,:,:)
    vars_fly(:,:,:,5,ie,je,ke) = svec(:,:,:,3) * the(:,:,:)

 ! Source term (Horizontal)         
    vars_src(:,:,:,1,ie,je,ke) = 0.0D0
    vars_src(:,:,:,2,ie,je,ke) =  fcori(:,:,:)*v(:,:,:) !- prs_x(:,:,:)
    vars_src(:,:,:,3,ie,je,ke) = -fcori(:,:,:)*u(:,:,:) !- prs_y(:,:,:)
    vars_src(:,:,:,4,ie,je,ke) = 0.0D0 ! computed by the vertical dicretization 
    vars_src(:,:,:,5,ie,je,ke) = 0.0D0

    edgevars_edges_spd_int(:,:,1,ie,je,ke) = spd(:,1, :)
    edgevars_edges_spd_int(:,:,2,ie,je,ke) = spd(nx,:,:)
    edgevars_edges_spd_int(:,:,3,ie,je,ke) = spd(:,nx,:)
    edgevars_edges_spd_int(:,:,4,ie,je,ke) = spd(1,:, :)


   do eq = 1, neq
    edgevars_edges_normflux_int(:,:,eq,1,ie,je,ke) = vars_fly(:,1, :,eq,ie,je,ke)
    edgevars_edges_normflux_int(:,:,eq,2,ie,je,ke) = vars_flx(nx,:,:,eq,ie,je,ke)
    edgevars_edges_normflux_int(:,:,eq,3,ie,je,ke) = vars_fly(:,nx,:,eq,ie,je,ke)
    edgevars_edges_normflux_int(:,:,eq,4,ie,je,ke) = vars_flx(1,:, :,eq,ie,je,ke)
   end do


  !South & North  face (y-flux)  
    !------------------
    do  face = 1,3, 2
        vec_ext(:,:,:) = edgevars_edges_vec_ext(:,:,:,face,ie,je,ke)
     Call   Extract_ext_Vars(ie,je,ke,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

     edgevars_edges_normflux_ext(:,:,1,face,ie,je,ke) =  vec_ext(:,:,3)    !rho*v 
     edgevars_edges_normflux_ext(:,:,2,face,ie,je,ke) =  vec_ext(:,:,3) * u_ext(:,:)
     edgevars_edges_normflux_ext(:,:,3,face,ie,je,ke) =  vec_ext(:,:,3) * v_ext(:,:) + prp_ext(:,:)
     edgevars_edges_normflux_ext(:,:,4,face,ie,je,ke) =  vec_ext(:,:,3) * w_ext(:,:)
     edgevars_edges_normflux_ext(:,:,5,face,ie,je,ke) =  vec_ext(:,:,3) * the_ext(:,:)
     edgevars_edges_spd_ext(:,:,face,ie,je,ke) =  spd_ext(:,:)
    end do


    !East/West  face (x-flux)  
    !------------------
    do  face = 2,4, 2
        vec_ext(:,:,:) = edgevars_edges_vec_ext(:,:,:,face,ie,je,ke)
     Call   Extract_ext_Vars(ie,je,ke,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)

     edgevars_edges_normflux_ext(:,:,1,face,ie,je,ke) =  vec_ext(:,:,2)    !rho*u 
     edgevars_edges_normflux_ext(:,:,2,face,ie,je,ke) =  vec_ext(:,:,2) * u_ext(:,:) + prp_ext(:,:)
     edgevars_edges_normflux_ext(:,:,3,face,ie,je,ke) =  vec_ext(:,:,2) * v_ext(:,:)
     edgevars_edges_normflux_ext(:,:,4,face,ie,je,ke) =  vec_ext(:,:,2) * w_ext(:,:)
     edgevars_edges_normflux_ext(:,:,5,face,ie,je,ke) =  vec_ext(:,:,2) * the_ext(:,:)
     edgevars_edges_spd_ext(:,:,face,ie,je,ke) =  spd_ext(:,:)
    end do

          Call  DGNH_Horizontal_rhs(ie,je,ke,z_src,rhs_h)  

  End Subroutine DG2d_Horizontal_solver

!!=============================================================================================
  Subroutine  Extract_ext_Vars(ie,je,ke,face,vec_ext,u_ext,v_ext,w_ext,the_ext,prp_ext,spd_ext)
  implicit none
!    type(edge_t), intent(in), dimension(6) :: edge
    integer, intent(in) :: face, ie, je, ke
    real(kind=real_kind), intent(in)  :: vec_ext(nx,nx,neq)

    real(kind=real_kind), intent(out),dimension(nx,nx) :: u_ext, v_ext,w_ext, &
                                                          the_ext,prp_ext,spd_ext
    real(kind=real_kind) :: rho_ext(nx,nx), prs_ext(nx,nx)

    ! NH variables along exterior boundaries for "nx" vertical column  
                    rho_ext(:,:) = vec_ext(:,:,1)   + edgevars_edges_rhb(:,:,face,ie,je,ke)
                      u_ext(:,:) = vec_ext(:,:,2) / rho_ext(:,:)
                      v_ext(:,:) = vec_ext(:,:,3) / rho_ext(:,:)
                      w_ext(:,:) = vec_ext(:,:,4) / rho_ext(:,:)
                    the_ext(:,:) = (vec_ext(:,:,5) + edgevars_edges_rtb(:,:,face,ie,je,ke))/rho_ext(:,:)

       Call  Compute_Ext_Vars(rho_ext,the_ext, prs_ext,spd_ext)

                    prp_ext(:,:) = prs_ext(:,:) - edgevars_edges_prb(:,:,face,ie,je,ke)

  End  Subroutine  Extract_ext_Vars

!!==================================================================================
  Subroutine DGNH_Horizontal_rhs(ie,je,ke, z_src, rhs_h)
   Implicit none
!    type(element_t), intent(in) :: vars
!    type(edge_t), intent(inout), dimension(6) :: edgevars
integer, intent(in) :: ie, je, ke
    real(kind=real_kind), intent(in), dimension(nx,nx,nx,neq) :: z_src 
    real(kind=real_kind), intent(out), dimension(nx,nx,nx,neq) :: rhs_h

    real(kind=real_kind), dimension(nx,nx,nx) :: flx, fly, divf, rhs ,src, h_src 
    real(kind=real_kind), dimension(4,nx,nx) :: numflux, nvec_int, nvec_ext
    real(kind=real_kind) :: dflx,dfly, dx2,dy2, fjmax(4)
    integer :: col, i,j, n, eq

    dx2 = 2.0D0 / delx
    dy2 = 2.0D0 / dely

    fjmax(:) = Compute_2DFluxJacobian(ie,je,ke)

    do eq = 1, neq

    Call  LxF_NH2d_Flux(ie,je,ke,eq,fjmax,numflux)

      flx(:,:,:) = vars_flx(:,:,:,eq,ie,je,ke)
      fly(:,:,:) = vars_fly(:,:,:,eq,ie,je,ke)
      h_src(:,:,:) = vars_src(:,:,:,eq,ie,je,ke)

    ! Horizontal divergence 
    do col = 1, nx     !vertical colum
      do j = 1, nx
      do i = 1, nx

          dflx = 0.0D0
          do n = 1, nx
            dflx = dflx + flx(n,j,col) * weakder(n,i)
          end do

          dfly = 0.0D0
          do n = 1, nx
            dfly = dfly + fly(i,n,col) * weakder(n,j)
          end do

          divf(i,j,col) = dflx * dx2  + dfly * dy2
      end do
      end do
    end do


    ! RHS = divF + \sum bdryF
    ! start with volume divF term, then add each bdry term one-by-one
    rhs(:,:,:) = divf(:,:,:)

 !! Numerical flux follows  Rt. minus Lt. wall,  in a FV sense 

    ! "upper" edges -- east, north
    rhs(nx,:,:) = rhs(nx,:,:) - numflux(2,:,:) * 2.0D0/(delx*gllw(nx))
    rhs(:,nx,:) = rhs(:,nx,:) - numflux(3,:,:) * 2.0D0/(dely*gllw(nx))

    ! "lower" edges -- west, south
    rhs(1,:,:) = rhs(1,:,:) - numflux(4,:,:) * 2.0D0/(delx*gllw(1))
    rhs(:,1,:) = rhs(:,1,:) - numflux(1,:,:) * 2.0D0/(dely*gllw(1))


   ! combining horizontal and vertical contributions 
    rhs_h(:,:,:,eq) = rhs(:,:,:) + h_src(:,:,:) + z_src(:,:,:,eq)
   enddo

  End subroutine DGNH_Horizontal_rhs

!!==================================================================================
   Function  Compute_2DFluxJacobian(ie,je,ke) result(fj_max)
   Implicit none
integer, intent(in) :: ie, je, ke
    Integer, parameter :: isouth=1,inorth=3,iwest=4,ieast=2

!    type(edge_t), intent(in), dimension(6) :: edgevars
    real(kind=real_kind) :: fj_max(4)

    real(kind=real_kind), dimension(4,nx,nx) :: nvec_int, nvec_ext
    real(kind=real_kind), dimension(nx,nx) :: val_int, val_ext, ev(4)
    integer :: i,j, e

 ! Collect Normal component of velocity on each face  for "nx" vertical points 
      nvec_int(isouth,:,:) =  edgevars_edges_v_int(:,:,1,ie,je,ke)
      nvec_int(ieast ,:,:) =  edgevars_edges_u_int(:,:,2,ie,je,ke)
      nvec_int(inorth,:,:) =  edgevars_edges_v_int(:,:,3,ie,je,ke)
      nvec_int(iwest ,:,:) =  edgevars_edges_u_int(:,:,4,ie,je,ke)

      nvec_ext(isouth,:,:) =  edgevars_edges_v_ext(:,:,1,ie,je,ke)
      nvec_ext(ieast ,:,:) =  edgevars_edges_u_ext(:,:,2,ie,je,ke)
      nvec_ext(inorth,:,:) =  edgevars_edges_v_ext(:,:,3,ie,je,ke)
      nvec_ext(iwest ,:,:) =  edgevars_edges_u_ext(:,:,4,ie,je,ke)

    do e = 1, 4
        do j = 1, nx
        do i = 1, nx
           val_int(i,j) =  abs(nvec_int(e,i,j)) + edgevars_edges_spd_int(i,j,e,ie,je,ke)
           val_ext(i,j) =  abs(nvec_ext(e,i,j)) + edgevars_edges_spd_ext(i,j,e,ie,je,ke)
        end do
        end do
        fj_max(e)  = max(maxval(val_int),maxval(val_ext))
    end do

   End Function  Compute_2DFluxJacobian 

!!==================================================================================
   Subroutine LxF_NH2d_Flux(ie,je,ke,eqn,fj_max,num_flux)
   Implicit none
    integer, intent(in) :: eqn,ie,je,ke
!    type(edge_t), intent(in), dimension(6) :: edgevars
    real(kind=real_kind), intent(in)  ::  fj_max(4)
    real(kind=real_kind), intent(out), dimension(4,nx,nx) :: num_flux

    real(kind=real_kind), dimension(nx,nx) :: psi_int,psi_ext,flx_int,flx_ext
    real(kind=real_kind) ::  sn(4)
    integer :: e

 ! normal vectors on horizontal sides 
   sn(1) = -1.0D0
   sn(2) = 1.0D0
   sn(3) = 1.0D0
   sn(4) = -1.0D0

 ! Lax-Fred Numerical fluxes on each face for every "eqn"  

    do e = 1, 4
       psi_int(:,:) = edgevars_edges_vec_int(:,:,eqn,e,ie,je,ke)
       psi_ext(:,:) = edgevars_edges_vec_ext(:,:,eqn,e,ie,je,ke)
       flx_int(:,:) = edgevars_edges_normflux_int(:,:,eqn,e,ie,je,ke) !* sn(e)
       flx_ext(:,:) = edgevars_edges_normflux_ext(:,:,eqn,e,ie,je,ke) !* sn(e)
      num_flux(e,:,:) =  ((flx_ext(:,:) + flx_int(:,:))*sn(e)  -  &
                           fj_max(e) * (psi_ext(:,:) - psi_int(:,:))) * 0.5D0
    end do

  End subroutine LxF_NH2d_Flux


!!=========================================================================================================

End module vertical_split_mod
