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
  public :: LxF_NH3d_Flux,Compute_MaxFluxJacobian,flux_divergence
  public :: compute_dgnh_rhs, Horizontal_strong_grads

contains

!!==================================================================================
  Subroutine compute_dgnh_rhs(nlx,nly,nlz)
   Implicit none 

    integer, intent(in) :: nlx, nly, nlz
    integer :: eqn,ie,je,ke,i,j,k

!    Call  Compute_MaxFluxJacobian(nlx,nly,nlz)

!!    do eqn = 1, neq 

!    Call  LxF_NH3d_Flux(nlx,nly,nlz)

!!      flx(:,:,:) = vars%flx(:,:,:,eqn) 
!!      fly(:,:,:) = vars%fly(:,:,:,eqn)
!!      flz(:,:,:) = vars%flz(:,:,:,eqn)

!!      src(:,:,:) = vars%src(:,:,:,eqn) 

!    Call flux_divergence(nlx,nly,nlz)

do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx
do eqn = 1, neq

    ! RHS = divF + \sum bdryF
    ! start with volume divF term, then add each bdry term one-by-one
!    rhs(:,:,:) = divf(:,:,:)   

 !! Numerical flux follows  Rt. minus Lt. wall,  in a FV sense 

    ! "upper" edges -- east, north, top
    elem(ie,je,ke)%divf(nx,:,:,eqn) = elem(ie,je,ke)%divf(nx,:,:,eqn) - elem(ie,je,ke)%numflux(:,:,eqn,2) * 2.0D0/(delx*gllw(nx))
    elem(ie,je,ke)%divf(:,nx,:,eqn) = elem(ie,je,ke)%divf(:,nx,:,eqn) - elem(ie,je,ke)%numflux(:,:,eqn,3) * 2.0D0/(dely*gllw(nx))
    elem(ie,je,ke)%divf(:,:,nx,eqn) = elem(ie,je,ke)%divf(:,:,nx,eqn) - elem(ie,je,ke)%numflux(:,:,eqn,6) * 2.0D0/(delz*gllw(nx))

    ! "lower" edges -- west, south, bottom
    elem(ie,je,ke)%divf(1,:,:,eqn) = elem(ie,je,ke)%divf(1,:,:,eqn) - elem(ie,je,ke)%numflux(:,:,eqn,4) * 2.0D0/(delx*gllw(1))
    elem(ie,je,ke)%divf(:,1,:,eqn) = elem(ie,je,ke)%divf(:,1,:,eqn) - elem(ie,je,ke)%numflux(:,:,eqn,1) * 2.0D0/(dely*gllw(1))
    elem(ie,je,ke)%divf(:,:,1,eqn) = elem(ie,je,ke)%divf(:,:,1,eqn) - elem(ie,je,ke)%numflux(:,:,eqn,5) * 2.0D0/(delz*gllw(1))


    elem(ie,je,ke)%rhs(:,:,:,eqn) = elem(ie,je,ke)%divf(:,:,:,eqn) + vars(ie,je,ke)%src(:,:,:,eqn)  

   end do
   end do
   end do
   end do

  end subroutine compute_dgnh_rhs

!!==================================================================================
 attributes(global) Subroutine compute_rhs()
   Implicit none 

!    integer, intent(in) :: nlx, nly, nlz
    integer :: eqn,ie,je,ke,i,j,k

!    Call  Compute_MaxFluxJacobian(nlx,nly,nlz)

!!    do eqn = 1, neq 

!    Call  LxF_NH3d_Flux(nlx,nly,nlz)

!!      flx(:,:,:) = vars%flx(:,:,:,eqn) 
!!      fly(:,:,:) = vars%fly(:,:,:,eqn)
!!      flz(:,:,:) = vars%flz(:,:,:,eqn)

!!      src(:,:,:) = vars%src(:,:,:,eqn) 

!    Call flux_divergence(nlx,nly,nlz)

!do ke = 1, nlz
!do je = 1, nly
!do ie = 1, nlx
!do eqn = 1, neq
i = threadIdx%x
j = threadIdx%y
eqn = threadIdx%z
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z

    ! RHS = divF + \sum bdryF
    ! start with volume divF term, then add each bdry term one-by-one
!    rhs(:,:,:) = divf(:,:,:)   

 !! Numerical flux follows  Rt. minus Lt. wall,  in a FV sense 

    ! "upper" edges -- east, north, top
    elem_d(ie,je,ke)%divf(nx,i,j,eqn) = elem_d(ie,je,ke)%divf(nx,i,j,eqn) - elem_d(ie,je,ke)%numflux(i,j,eqn,2) * 2.0D0/(delx_d*gllw_d(nx))
    elem_d(ie,je,ke)%divf(i,nx,j,eqn) = elem_d(ie,je,ke)%divf(i,nx,j,eqn) - elem_d(ie,je,ke)%numflux(i,j,eqn,3) * 2.0D0/(dely_d*gllw_d(nx))
    elem_d(ie,je,ke)%divf(i,j,nx,eqn) = elem_d(ie,je,ke)%divf(i,j,nx,eqn) - elem_d(ie,je,ke)%numflux(i,j,eqn,6) * 2.0D0/(delz_D*gllw_d(nx))

    ! "lower" edges -- west, south, bottom
    elem_d(ie,je,ke)%divf(1,i,j,eqn) = elem_d(ie,je,ke)%divf(1,i,j,eqn) - elem_d(ie,je,ke)%numflux(i,j,eqn,4) * 2.0D0/(delx_d*gllw_d(1))
    elem_d(ie,je,ke)%divf(i,1,j,eqn) = elem_d(ie,je,ke)%divf(i,1,j,eqn) - elem_d(ie,je,ke)%numflux(i,j,eqn,1) * 2.0D0/(dely_d*gllw_d(1))
    elem_d(ie,je,ke)%divf(i,j,1,eqn) = elem_d(ie,je,ke)%divf(i,j,1,eqn) - elem_d(ie,je,ke)%numflux(i,j,eqn,5) * 2.0D0/(delz_d*gllw_d(1))


    elem_d(ie,je,ke)%rhs(i,j,:,eqn) = elem_d(ie,je,ke)%divf(i,j,:,eqn) + vars_d(ie,je,ke)%src(i,j,:,eqn)  

!   end do
!   end do
!   end do
!   end do

  end subroutine compute_rhs

!!==================================================================================
   Subroutine Compute_MaxFluxJacobian(nlx,nly,nlz)
   Implicit none
    Integer, parameter :: isouth=1,inorth=3,iwest=4,ieast=2,ibot=5,itop=6 
 integer, intent(in) :: nlx,nly,nlz

    real(kind=real_kind), dimension(nx,nx,6) :: nvec_int, nvec_ext
    real(kind=real_kind), dimension(nx,nx) :: val_int, val_ext, ev(6)
    integer :: i,j, e,k,l,ie,je,ke

do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx

 ! Collect Normal component of velocity on each face of 3d element 
      nvec_int(:,:,isouth) =  edgevars(ie,je,ke)%edges(1)%v_int(:,:)
      nvec_int(:,:,ieast ) =  edgevars(ie,je,ke)%edges(2)%u_int(:,:)
      nvec_int(:,:,inorth) =  edgevars(ie,je,ke)%edges(3)%v_int(:,:)
      nvec_int(:,:,iwest ) =  edgevars(ie,je,ke)%edges(4)%u_int(:,:)
      nvec_int(:,:,ibot  ) =  edgevars(ie,je,ke)%edges(5)%w_int(:,:)
      nvec_int(:,:,itop  ) =  edgevars(ie,je,ke)%edges(6)%w_int(:,:)

      nvec_ext(:,:,isouth) =  edgevars(ie,je,ke)%edges(1)%v_ext(:,:)
      nvec_ext(:,:,ieast ) =  edgevars(ie,je,ke)%edges(2)%u_ext(:,:)
      nvec_ext(:,:,inorth) =  edgevars(ie,je,ke)%edges(3)%v_ext(:,:)
      nvec_ext(:,:,iwest ) =  edgevars(ie,je,ke)%edges(4)%u_ext(:,:)
      nvec_ext(:,:,ibot  ) =  edgevars(ie,je,ke)%edges(5)%w_ext(:,:)
      nvec_ext(:,:,itop  ) =  edgevars(ie,je,ke)%edges(6)%w_ext(:,:)


    do e = 1, 6
        do j = 1, nx 
        do i = 1, nx 
           val_int(i,j) =  abs(nvec_int(i,j,e)) + edgevars(ie,je,ke)%edges(e)%spd_int(i,j) 
           val_ext(i,j) =  abs(nvec_ext(i,j,e)) + edgevars(ie,je,ke)%edges(e)%spd_ext(i,j) 
        end do 
        end do 
        elem(ie,je,ke)%fjmax(e)  = max(maxval(val_int),maxval(val_ext)) 
    end do

end do
end do
end do

   end Subroutine Compute_MaxFluxJacobian 

!!==================================================================================
  attributes(global) Subroutine Compute_Jacobian()
   Implicit none
    Integer, parameter :: isouth=1,inorth=3,iwest=4,ieast=2,ibot=5,itop=6 
! integer, intent(in) :: nlx,nly,nlz

    real(kind=real_kind), dimension(nx,nx,6) :: nvec_int, nvec_ext
    real(kind=real_kind), dimension(nx,nx) :: val_int, val_ext, ev(6)
    integer :: i,j, e,k,l,ie,je,ke
    real(kind=real_kind) :: tmp

!do ke = 1, nlz
!do je = 1, nly
!do ie = 1, nlx
ie = threadIdx%x
je = blockIdx%x
ke = blockIdx%y
 ! Collect Normal component of velocity on each face of 3d element 
      nvec_int(:,:,isouth) =  edgevars_d(ie,je,ke)%edges(1)%v_int(:,:)
      nvec_int(:,:,ieast ) =  edgevars_d(ie,je,ke)%edges(2)%u_int(:,:)
      nvec_int(:,:,inorth) =  edgevars_d(ie,je,ke)%edges(3)%v_int(:,:)
      nvec_int(:,:,iwest ) =  edgevars_d(ie,je,ke)%edges(4)%u_int(:,:)
      nvec_int(:,:,ibot  ) =  edgevars_d(ie,je,ke)%edges(5)%w_int(:,:)
      nvec_int(:,:,itop  ) =  edgevars_d(ie,je,ke)%edges(6)%w_int(:,:)

      nvec_ext(:,:,isouth) =  edgevars_d(ie,je,ke)%edges(1)%v_ext(:,:)
      nvec_ext(:,:,ieast ) =  edgevars_d(ie,je,ke)%edges(2)%u_ext(:,:)
      nvec_ext(:,:,inorth) =  edgevars_d(ie,je,ke)%edges(3)%v_ext(:,:)
      nvec_ext(:,:,iwest ) =  edgevars_d(ie,je,ke)%edges(4)%u_ext(:,:)
      nvec_ext(:,:,ibot  ) =  edgevars_d(ie,je,ke)%edges(5)%w_ext(:,:)
      nvec_ext(:,:,itop  ) =  edgevars_d(ie,je,ke)%edges(6)%w_ext(:,:)


    do e = 1, 6
        do j = 1, nx 
        do i = 1, nx 
           val_int(i,j) =  abs(nvec_int(i,j,e)) + edgevars_d(ie,je,ke)%edges(e)%spd_int(i,j) 
           val_ext(i,j) =  abs(nvec_ext(i,j,e)) + edgevars_d(ie,je,ke)%edges(e)%spd_ext(i,j) 
        end do 
        end do 
!        elem_d(ie,je,ke)%fjmax(e)  = max(maxval(val_int),maxval(val_ext)) 
        tmp=val_int(1,1)
        do k = 1, nx
        do l = 1, nx
                tmp = max(tmp, val_int(l,k))
                tmp = max(tmp, val_ext(l,k))
        end do
        end do
        elem_d(ie,je,ke)%fjmax(e) = tmp

    end do

!end do
!end do
!end do

   end Subroutine Compute_Jacobian 

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

do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx
    do e = 1, 6
    do eqn = 1, neq
    do j = 1, nx
    do i = 1, nx

       psi_int = edgevars(ie,je,ke)%edges(e)%vec_int(i,j,eqn) 
       psi_ext = edgevars(ie,je,ke)%edges(e)%vec_ext(i,j,eqn) 
       flx_int = edgevars(ie,je,ke)%edges(e)%normflux_int(i,j,eqn) !* sn(e)
       flx_ext = edgevars(ie,je,ke)%edges(e)%normflux_ext(i,j,eqn) !* sn(e)

      elem(ie,je,ke)%numflux(i,j,eqn,e) =  ((flx_ext + flx_int)*sn(e)  - elem(ie,je,ke)%fjmax(e) * (psi_ext - psi_int)) * 0.5D0

    end do
end do
end do
end do
end do
end do
end do

  End subroutine LxF_NH3d_Flux

!!==================================================================================
  attributes(global) Subroutine LxF_Flux()
   Implicit none

!    integer, intent(in) :: nlx,nly,nlz
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
i = threadIdx%x
j = threadIdx%y
eqn = threadIdx%z
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z

!do ke = 1, nlz
!do je = 1, nly
!do ie = 1, nlx
    do e = 1, 6
!    do eqn = 1, neq
!    do j = 1, nx
!    do i = 1, nx

       psi_int = edgevars_d(ie,je,ke)%edges(e)%vec_int(i,j,eqn) 
       psi_ext = edgevars_d(ie,je,ke)%edges(e)%vec_ext(i,j,eqn) 
       flx_int = edgevars_d(ie,je,ke)%edges(e)%normflux_int(i,j,eqn) !* sn(e)
       flx_ext = edgevars_d(ie,je,ke)%edges(e)%normflux_ext(i,j,eqn) !* sn(e)

      elem_d(ie,je,ke)%numflux(i,j,eqn,e) =  ((flx_ext + flx_int)*sn(e)  - elem_d(ie,je,ke)%fjmax(e) * (psi_ext - psi_int)) * 0.5D0

    end do
!end do
!end do
!end do
!end do
!end do
!end do

  End subroutine LxF_Flux

!!==================================================================================
  ! Compute the flux divergence according to the weak form of DG
  Subroutine Flux_Divergence(nlx,nly,nlz)
   Implicit none
integer, intent(in) :: nlx,nly,nlz
    real(kind=real_kind) :: dflx, dfly, dflz
    integer :: i,j,k,a,eqn,ie,je,ke

    ! Use pre-computed `weak form deriv matrix' weakder_ij = der_ij w_i / w_j
    ! to compute the derivative in each dimension.

do ke = 1, nlz
do je = 1, nly
do ie = 1, nlx

do eqn = 1, neq
    do k = 1, nx
      do j = 1, nx
        do i = 1, nx

          dflx = 0.0D0
          do a = 1, nx
            dflx = dflx + vars(ie,je,ke)%flx(a,j,k,eqn) * weakder(a,i)
          end do

          dfly = 0.0D0
          do a = 1, nx
            dfly = dfly + vars(ie,je,ke)%fly(i,a,k,eqn) * weakder(a,j)
          end do

          dflz = 0.0D0
          do a = 1, nx
            dflz = dflz + vars(ie,je,ke)%flz(i,j,a,eqn) * weakder(a,k)
          end do

          ! add the 2/dx geometric factors
          !divf(i,j,k) = dflx * 2/delx + dfly * 2/dely + dflz * 2/delz
          elem(ie,je,ke)%divf(i,j,k,eqn) = (dflx /delx + dfly /dely + dflz /delz)*2.0D0
        end do
      end do
    end do
end do
end do
end do
end do

  end subroutine Flux_Divergence

!!==================================================================================
  ! Compute the flux divergence according to the weak form of DG
 attributes(global) Subroutine Flux_Div()
   Implicit none
!integer, intent(in) :: nlx,nly,nlz
    real(kind=real_kind) :: dflx, dfly, dflz
    integer :: i,j,k,a,eqn,ie,je,ke

    ! Use pre-computed `weak form deriv matrix' weakder_ij = der_ij w_i / w_j
    ! to compute the derivative in each dimension.

i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z
!do ke = 1, nlz
!do je = 1, nly
!do ie = 1, nlx

do eqn = 1, neq
!    do k = 1, nx
!      do j = 1, nx
!        do i = 1, nx

          dflx = 0.0D0
          do a = 1, nx
            dflx = dflx + vars_d(ie,je,ke)%flx(a,j,k,eqn) * weakder_d(a,i)
          end do

          dfly = 0.0D0
          do a = 1, nx
            dfly = dfly + vars_d(ie,je,ke)%fly(i,a,k,eqn) * weakder_d(a,j)
          end do

          dflz = 0.0D0
          do a = 1, nx
            dflz = dflz + vars_d(ie,je,ke)%flz(i,j,a,eqn) * weakder_d(a,k)
          end do

          ! add the 2/dx geometric factors
          !divf(i,j,k) = dflx * 2/delx + dfly * 2/dely + dflz * 2/delz
          elem_d(ie,je,ke)%divf(i,j,k,eqn) = (dflx /delx_d + dfly /dely_d + dflz /delz_d)*2.0D0
!        end do
!      end do
!    end do
!end do
!end do
!end do
end do

  end subroutine Flux_Div

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
