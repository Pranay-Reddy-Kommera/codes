!-------------------------------------------------------------------
! R.Nair NCAR/scd 10/02
! Rewritten 06/2016 by Francois Hebert, SIParCS/Cornell
!
! Implements the weak form DG method for evaluating the
! semi-discrete RHS of the advection equation.
! Computes the RHS of `dU/dt = rhs'
!-------------------------------------------------------------------

module dg3d_rhs_mod
  use basic_mod
  use element_mod
  use gauss_quadrature_mod
  use openacc

  implicit none
  private :: LxF_NH3d_Flux,Compute_MaxFluxJacobian 
  public :: compute_dgnh_rhs, Horizontal_strong_grads

contains

!!==================================================================================

Subroutine compute_dgnh_rhs(nlx,nly,nlz)
   Implicit none

    integer, intent(in) :: nlx, nly, nlz
    integer :: eqn,ie,je,ke,i,j,k

    Call  Compute_MaxFluxJacobian(nlx,nly,nlz)

!    do eqn = 1, neq

    Call  LxF_NH3d_Flux(nlx,nly,nlz)

!      flx(:,:,:) = vars_flx(ie,je,ke,:,:,:,eqn) 
!      fly(:,:,:) = vars_fly(ie,je,ke,:,:,:,eqn)
!      flz(:,:,:) = vars_flz(ie,je,ke,:,:,:,eqn)

!      src(:,:,:) = vars_src(ie,je,ke,:,:,:,eqn) 

    Call flux_divergence(nlx,nly,nlz)


!$acc parallel num_workers(8) vector_length(16)
!$acc loop gang worker vector collapse(4)
do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx
do eqn = 1, neq
    ! RHS = divF + \sum bdryF
    ! start with volume divF term, then add each bdry term one-by-one
!    rhs(:,:,:) = divf(:,:,:)

!! Numerical flux follows  Rt. minus Lt. wall,  in a FV sense 

    ! "upper" edges -- east, north, top
    divf(nx,:,:,eqn,ie,je,ke) = divf(nx,:,:,eqn,ie,je,ke) - numflux(:,:,eqn,2,ie,je,ke) * 2.0D0/(delx*gllw(nx))
    divf(:,nx,:,eqn,ie,je,ke) = divf(:,nx,:,eqn,ie,je,ke) - numflux(:,:,eqn,3,ie,je,ke) * 2.0D0/(dely*gllw(nx))
    divf(:,:,nx,eqn,ie,je,ke) = divf(:,:,nx,eqn,ie,je,ke) - numflux(:,:,eqn,6,ie,je,ke) * 2.0D0/(delz*gllw(nx))

    ! "lower" edges -- west, south, bottom
    divf(1,:,:,eqn,ie,je,ke) = divf(1,:,:,eqn,ie,je,ke) - numflux(:,:,eqn,4,ie,je,ke) * 2.0D0/(delx*gllw(1))
    divf(:,1,:,eqn,ie,je,ke) = divf(:,1,:,eqn,ie,je,ke) - numflux(:,:,eqn,1,ie,je,ke) * 2.0D0/(dely*gllw(1))
    divf(:,:,1,eqn,ie,je,ke) = divf(:,:,1,eqn,ie,je,ke) - numflux(:,:,eqn,5,ie,je,ke) * 2.0D0/(delz*gllw(1))


    elem_rhs(:,:,:,eqn,ie,je,ke) = divf(:,:,:,eqn,ie,je,ke) + vars_src(:,:,:,eqn,ie,je,ke)

   end do
   end do
   end do
   end do
!$acc end parallel

  end subroutine compute_dgnh_rhs

!!==================================================================================

Subroutine Compute_MaxFluxJacobian(nlx,nly,nlz)
   Implicit none
    Integer, parameter :: isouth=1,inorth=3,iwest=4,ieast=2,ibot=5,itop=6
integer, intent(in) :: nlx,nly,nlz

    real(kind=real_kind), dimension(nx,nx,6) :: nvec_int, nvec_ext
    real(kind=real_kind), dimension(nx,nx) :: val_int, val_ext, ev(6)
    integer :: i,j, e,k,l,ie,je,ke

!$acc parallel num_workers(8) vector_length(16)
!$acc loop gang worker vector collapse(3) private(nvec_int,nvec_ext,val_int,val_ext)
do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx

 ! Collect Normal component of velocity on each face of 3d element 
      nvec_int(:,:,isouth) =  edgevars_edges_v_int(:,:,1,ie,je,ke)
      nvec_int(:,:,ieast ) =  edgevars_edges_u_int(:,:,2,ie,je,ke)
      nvec_int(:,:,inorth) =  edgevars_edges_v_int(:,:,3,ie,je,ke)
      nvec_int(:,:,iwest ) =  edgevars_edges_u_int(:,:,4,ie,je,ke)
      nvec_int(:,:,ibot  ) =  edgevars_edges_w_int(:,:,5,ie,je,ke)
      nvec_int(:,:,itop  ) =  edgevars_edges_w_int(:,:,6,ie,je,ke)

      nvec_ext(:,:,isouth) =  edgevars_edges_v_ext(:,:,1,ie,je,ke)
      nvec_ext(:,:,ieast ) =  edgevars_edges_u_ext(:,:,2,ie,je,ke)
      nvec_ext(:,:,inorth) =  edgevars_edges_v_ext(:,:,3,ie,je,ke)
      nvec_ext(:,:,iwest ) =  edgevars_edges_u_ext(:,:,4,ie,je,ke)
      nvec_ext(:,:,ibot  ) =  edgevars_edges_w_ext(:,:,5,ie,je,ke)
      nvec_ext(:,:,itop  ) =  edgevars_edges_w_ext(:,:,6,ie,je,ke)


    do e = 1, 6
        do j = 1, nx
        do i = 1, nx
           val_int(i,j) =  abs(nvec_int(i,j,e)) + edgevars_edges_spd_int(i,j,e,ie,je,ke)
           val_ext(i,j) =  abs(nvec_ext(i,j,e)) + edgevars_edges_spd_ext(i,j,e,ie,je,ke)
        end do
        end do
        fjmax(e,ie,je,ke)  = max(maxval(val_int),maxval(val_ext))
    end do

end do
end do
end do
!$acc end parallel

   end Subroutine Compute_MaxFluxJacobian

!!==================================================================================

Subroutine LxF_NH3d_Flux(nlx,nly,nlz)
   Implicit none

integer, intent(in) :: nlx,nly,nlz
    real(kind=real_kind) :: psi_int,psi_ext,flx_int,flx_ext
    real(kind=real_kind) ::  sn(6)
    integer :: e,ie,je,ke,i,j,eqn

   sn(1) = -1.0D0
   sn(2) = 1.0D0
   sn(3) = 1.0D0
   sn(4) = -1.0D0
   sn(5) = -1.0D0
   sn(6) = 1.0D0

 ! Lax-Fred Numerical fluxes on each face for every "eqn"  

!$acc parallel num_workers(8) vector_length(16)
!$acc loop gang worker vector collapse(7) private(psi_int,psi_ext,flx_int,flx_ext)
do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx
    do e = 1, 6
    do eqn = 1, neq
    do j = 1, nx
    do i = 1, nx
       psi_int = edgevars_edges_vec_int(i,j,eqn,e,ie,je,ke)
       psi_ext = edgevars_edges_vec_ext(i,j,eqn,e,ie,je,ke)
       flx_int = edgevars_edges_normflux_int(i,j,eqn,e,ie,je,ke) !* sn(e)
       flx_ext = edgevars_edges_normflux_ext(i,j,eqn,e,ie,je,ke) !* sn(e)

      numflux(i,j,eqn,e,ie,je,ke) =  ((flx_ext + flx_int)*sn(e)  - fjmax(e,ie,je,ke) * (psi_ext - psi_int)) * 0.5D0

    end do
end do
end do
end do
end do
end do
end do
!$acc end parallel

  End subroutine LxF_NH3d_Flux

!!==================================================================================
  ! Compute the flux divergence according to the weak form of DG
  Subroutine Flux_Divergence(nlx,nly,nlz)
   Implicit none
integer, intent(in) :: nlx, nly, nlz
    real(kind=real_kind) :: dflx, dfly, dflz
    integer :: i,j,k,a,eqn,ie,je,ke

    ! Use pre-computed `weak form deriv matrix' weakder_ij = der_ij w_i / w_j
    ! to compute the derivative in each dimension.

!$acc parallel num_workers(8) vector_length(16)
!$acc loop gang worker vector collapse(7) private(a,dflx,dfly,dflz)
do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx

do eqn = 1, neq
    do k = 1, nx
      do j = 1, nx
        do i = 1, nx

          dflx = 0.0D0
          do a = 1, nx
            dflx = dflx + vars_flx(a,j,k,eqn,ie,je,ke) * weakder(a,i)
          end do

          dfly = 0.0D0
          do a = 1, nx
            dfly = dfly + vars_fly(i,a,k,eqn,ie,je,ke) * weakder(a,j)
          end do

          dflz = 0.0D0
          do a = 1, nx
            dflz = dflz + vars_flz(i,j,a,eqn,ie,je,ke) * weakder(a,k)
          end do

          ! add the 2/dx geometric factors
          !divf(i,j,k) = dflx * 2/delx + dfly * 2/dely + dflz * 2/delz
          divf(i,j,k,eqn,ie,je,ke) = (dflx /delx + dfly /dely + dflz /delz)*2.0D0
        end do
      end do
    end do
end do
end do
end do
end do
!$acc end parallel

  end subroutine Flux_Divergence

!!==================================================================================
  ! Compute the flux divergence according to the weak form of DG
  Subroutine Horizontal_strong_grads(psi,psi_x,psi_y)
   Implicit none
    real(kind=real_kind), intent(in), dimension(nx,nx,nx) :: psi              
    real(kind=real_kind), intent(out), dimension(nx,nx,nx) :: psi_x, psi_y
    real(kind=real_kind) :: dflx, dfly, dflz
    integer :: i,j,k,a

    ! Use pre-computed `weak form deriv matrix' weakder_ij = der_ij w_i / w_j
    ! to compute the derivative in each dimension.
    do k = 1, nx
      do j = 1, nx
        do i = 1, nx

          dflx = 0.0D0
          do a = 1, nx
            dflx = dflx + psi(a,j,k) * der(i,a)
          end do

          dfly = 0.0D0
          do a = 1, nx
            dfly = dfly + psi(i,a,k) * der(j,a)
          end do

         !dflz = 0.0D0
         !do a = 1, nx
         !  dflz = dflz + flz(i,j,a) * weakder(a,k)
         !end do

            psi_x(i,j,k) = dflx * (2.0D0/delx) 
            psi_y(i,j,k) = dfly * (2.0D0/dely) 
          ! add the 2/dx geometric factors
          !divf(i,j,k) = dflx * 2/delx + dfly * 2/dely + dflz * 2/delz
        end do
      end do
    end do

  end subroutine Horizontal_strong_grads

!!==================================================================================

end module dg3d_rhs_mod
