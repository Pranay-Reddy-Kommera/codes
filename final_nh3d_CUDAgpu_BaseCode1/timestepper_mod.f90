!-------------------------------------------------------------------
! 11/04 R.Nair NCAR/SCD
! Rewritten 06/2016 by Francois Hebert, SIParCS/Cornell
!
! Third-order Runge-Kutta (RK3) time integrator
! Total variation diminishing (TVD) [or strong stability preserving]
!-------------------------------------------------------------------

module timestepper_mod
  use basic_mod
  use element_mod
  use prepare_edgevars_mod
  use euler_flux3d_mod, only: Make_Euler3d_Flux, Recover_State_Vars, Store_Hydrostatic_edgevars, nh_flux_maker,Make_ext_normalflux,Apply_Noflux_TopbotBC,Make_Flux,recover,store,nh_maker,Make_normalflux,Apply_TopbotBC
  use dg3d_rhs_mod 
  use vertical_split_mod, only: Compute_column_rhs,Compute_SplitDG_rhs
  use mpi_mod, only: send_data, recv_data, send,recv
  use cudafor

  implicit none
  private
  public :: ssp_rk3,hostToDevice,deviceToHost

  ! When running in parallel, the number of elements owned by this rank
  ! is not known a priori, so we must dynamiCally allocate the storage:
  !type(element_t), dimension(:,:,:), allocatable :: temp1, temp2, rk_stage 
! type(element_t), dimension(:,:,:), allocatable :: rk_stage 
!  type(surface_t), dimension(:,:,:), allocatable :: edgevars
!  type(rkvec_t), dimension(:,:,:), allocatable :: rk_stage

contains
!!==================================================================================
  Subroutine ssp_rk3(t, dt)
    Implicit none 
    real(kind=real_kind), intent(inout) :: t
    real(kind=real_kind), intent(in) :: dt
!    type(element_t), intent(inout), dimension(:,:,:) :: vars

    real(kind=real_kind), dimension(nx,nx,nx,neq) :: rhs, vec 
    real(kind=real_kind), dimension(nx,nx,nx) ::  the,rho, rhb,thb,rtb
    integer :: nlx, nly, nlz, eq 
    integer :: ie, je, ke, i,j,k,l,istat,ierr
    type(dim3) :: gridRec, blockRec, gridSto, blockSto, gridNes, blockNes, gridExt, blockExt, blockPac, gridPac1, gridPac2

    nlx = size(vars,1)
    nly = size(vars,2)
    nlz = size(vars,3)

!write (*,*) 'nlx is ',nlx,'nly is ',nly,'nlz is ',nlz

!    ! Allocate the temp storage, this only happens on the FIRST Call
!    if (.not. allocated(rk_stage)) allocate(rk_stage(nlx,nly,nlz))
!    if (.not. allocated(edgevars)) allocate(edgevars(nlx,nly,nlz))

blockExt = dim3(4,4,5)
gridExt = dim3(64,32,16)
blockPac = dim3(4,4,1)
gridPac1 = dim3(64,16,1)
gridPac2 = dim3(32,16,1)
blockSto = dim3(4,4,1)
gridSto = dim3(64,32,16)
blockNes = dim3(64,1,1)
gridNes = dim3(32,16,1)
blockRec = dim3(4,4,4)
gridRec = dim3(64,32,16)

    ! Communicate basic NH varibels vars%(rho,the,u,v,w), and copy the edge values 
!=====================================================
!call hostToDevice()
!     Call extract_edgevars()
        call extract_edges<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!=====================================================
if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)
!=====================================================
!!call hostToDevice()
!     Call pack_edgevars1()
!     Call pack_edgevars2()
        call pack_edges1<<<gridPac1, blockPac>>>()
        call pack_edges2<<<gridPac2, blockPac>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!======================================================
!!call hostToDevice()
!     Call communicate_edgevars_work()
        call communicate_edges_work<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!======================================================
!     Call send_data()
!     Call recv_data()
        Call send()
        Call recv()
!======================================================
!!call hostToDevice()
!     Call unpack_edgevars1()
!     Call unpack_edgevars2()
        call unpack_edges1<<<gridPac1, blockPac>>>()
        call unpack_edges2<<<gridPac2, blockPac>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================

          if (t == 0.0D0) then !one time initialization 
!========================================================
!call hostToDevice()
!            Call Store_Hydrostatic_edgevars()
        Call store<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
          endif 


    !----------------
    ! RK stage 1
    !----------------
!=======================================================
!call hostToDevice()
!        Call nh_flux_maker()
        call nh_maker<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!            Call Make_Euler3d_Flux()
        call Make_Flux<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!        Call Make_ext_normalflux()
        call Make_normalflux<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!        Call Apply_Noflux_TopbotBC()
        call Apply_TopbotBC<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call Compute_MaxFluxJacobian(nlx,nly,nlz)
        call Compute_Jacobian<<<gridNes, blockNes>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call LxF_NH3d_Flux(nlx,nly,nlz)
        Call LxF_Flux<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call flux_divergence(nlx,nly,nlz)
        Call flux_div<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!            Call Compute_dgnh_rhs(nlx,nly,nlz)
        Call Compute_rhs<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
    
!    if (use_dim_split_scheme) Call Compute_column_rhs()

!call hostToDevice()
!        Call nested1()
        Call nested1_cuda<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost

    !Compute/update  NH variables from the state vector
!=============================================
!call hostToDevice() 
!     Call Recover_State_Vars() 
        call recover<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=============================================

    !----------------
    ! RK stage 2
    !----------------

!! Communicate basic NH varibels vars%(rho,the,u,v,w), and copy the edge values 
!     Call extract_edgevars()
!     Call pack_edgevars1()
!     Call pack_edgevars2()
!     Call communicate_edgevars_work()
!     Call send_data()
!     Call recv_data()
!     Call unpack_edgevars1()
!     Call unpack_edgevars2()
    ! Communicate basic NH varibels vars%(rho,the,u,v,w), and copy the edge values 
!=====================================================
!call hostToDevice()
!     Call extract_edgevars()
        call extract_edges1<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!=====================================================
!if (.not. allocated(bufsendS)) call allocate_buffers1(nlx, nly, nlz)
!=====================================================
!!call hostToDevice()
!     Call pack_edgevars1()
!     Call pack_edgevars2()
        call pack_edges11<<<gridPac1, blockPac>>>()
        call pack_edges21<<<gridPac2, blockPac>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!======================================================
!!call hostToDevice()
!     Call communicate_edgevars_work()
        call communicate_edges_work1<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!======================================================
!     Call send_data()
!     Call recv_data()
        Call send()
        Call recv()
!======================================================
!!call hostToDevice()
!     Call unpack_edgevars1()
!     Call unpack_edgevars2()
        call unpack_edges11<<<gridPac1, blockPac>>>()
        call unpack_edges21<<<gridPac2, blockPac>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================

!        Call nh_flux_maker()
!            Call Make_Euler3d_Flux()
!        Call Make_ext_normalflux()
!        Call Apply_Noflux_TopbotBC()
!
!        Call Compute_MaxFluxJacobian(nlx,nly,nlz)
!        Call LxF_NH3d_Flux(nlx,nly,nlz)
!        Call flux_divergence(nlx,nly,nlz)
!            Call Compute_dgnh_rhs(nlx,nly,nlz)

!=======================================================
!call hostToDevice()
!        Call nh_flux_maker()
        call nh_maker<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!            Call Make_Euler3d_Flux()
        call Make_Flux<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!        Call Make_ext_normalflux()
        call Make_normalflux<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!        Call Apply_Noflux_TopbotBC()
        call Apply_TopbotBC<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call Compute_MaxFluxJacobian(nlx,nly,nlz)
        call Compute_Jacobian<<<gridNes, blockNes>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call LxF_NH3d_Flux(nlx,nly,nlz)
        Call LxF_Flux<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call flux_divergence(nlx,nly,nlz)
        Call flux_div<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!            Call Compute_dgnh_rhs(nlx,nly,nlz)
        Call Compute_rhs<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================

!!   if (use_dim_split_scheme) Call Compute_column_rhs()
!
!        Call nested2()
!
!   Call Recover_State_Vars() 
!    
!!    if (use_dim_split_scheme) Call Compute_column_rhs()

!call hostToDevice()
!        Call nested2()
        Call nested2_cuda<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost

    !Compute/update  NH variables from the state vector
!=============================================
!call hostToDevice() 
!     Call Recover_State_Vars() 
        call recover1<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=============================================

!   !----------------
!   ! RK stage 3
!   !----------------

!! Communicate basic NH varibels vars%(rho,the,u,v,w), and copy the edge values 
!     Call extract_edgevars()
!     Call pack_edgevars1()
!     Call pack_edgevars2()
!     Call communicate_edgevars_work()
!     Call send_data()
!     Call recv_data()
!     Call unpack_edgevars1()
!     Call unpack_edgevars2()

    ! Communicate basic NH varibels vars%(rho,the,u,v,w), and copy the edge values 
!=====================================================
!call hostToDevice()
!     Call extract_edgevars()
        call extract_edges2<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!=====================================================
!if (.not. allocated(bufsendS)) call allocate_buffers2(nlx, nly, nlz)
!=====================================================
!!call hostToDevice()
!     Call pack_edgevars1()
!     Call pack_edgevars2()
        call pack_edges12<<<gridPac1, blockPac>>>()
        call pack_edges22<<<gridPac2, blockPac>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!======================================================
!!call hostToDevice()
!     Call communicate_edgevars_work()
        call communicate_edges_work2<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!!call deviceToHost()
!======================================================
!     Call send_data()
!     Call recv_data()
        Call send()
        Call recv()
!======================================================
!!call hostToDevice()
!     Call unpack_edgevars1()
!     Call unpack_edgevars2()
        call unpack_edges12<<<gridPac1, blockPac>>>()
        call unpack_edges22<<<gridPac2, blockPac>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================

!        Call nh_flux_maker()
!            Call Make_Euler3d_Flux()
!        Call Make_ext_normalflux()
!        Call Apply_Noflux_TopbotBC()
!
!        Call Compute_MaxFluxJacobian(nlx,nly,nlz)
!        Call LxF_NH3d_Flux(nlx,nly,nlz)
!        Call flux_divergence(nlx,nly,nlz)
!            Call Compute_dgnh_rhs(nlx,nly,nlz)

!=======================================================
!call hostToDevice()
!        Call nh_flux_maker()
        call nh_maker<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!            Call Make_Euler3d_Flux()
        call Make_Flux<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!        Call Make_ext_normalflux()
        call Make_normalflux<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=======================================================
!call hostToDevice()
!        Call Apply_Noflux_TopbotBC()
        call Apply_TopbotBC<<<gridSto, blockSto>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call Compute_MaxFluxJacobian(nlx,nly,nlz)
        call Compute_Jacobian<<<gridNes, blockNes>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call LxF_NH3d_Flux(nlx,nly,nlz)
        Call LxF_Flux<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!        Call flux_divergence(nlx,nly,nlz)
        Call flux_div<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================
!call hostToDevice()
!            Call Compute_dgnh_rhs(nlx,nly,nlz)
        Call Compute_rhs<<<gridExt, blockExt>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost
!=======================================================

!!    if (use_dim_split_scheme) Call Compute_column_rhs()
!
!        Call nested3()
!
!    Call Recover_State_Vars() 
!    
!!    if (use_dim_split_scheme) Call Compute_column_rhs()

!call hostToDevice()
!        Call nested3()
        Call nested3_cuda<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost

    !Compute/update  NH variables from the state vector
!=============================================
!call hostToDevice() 
!     Call Recover_State_Vars() 
        call recover2<<<gridRec, blockRec>>>()
istat = cudaDeviceSynchronize()
!call deviceToHost()
!=============================================

    ! update global time
    t = t + dt

  end subroutine ssp_rk3 

!!==================================================================================
Subroutine hostToDevice()
implicit none

edgevars_d = edgevars
vars_d = vars
weakder_d = weakder
der_d = der
gllw_d = gllw
delx_d = delx
dely_d = dely
delz_d = delz
rk_stage_d = rk_stage
elem_d = elem

End Subroutine hostToDevice
!!==================================================================================

Subroutine deviceToHost()
implicit none

edgevars = edgevars_d
vars = vars_d
weakder = weakder_d
der = der_d
gllw = gllw_d
delx = delx_d
dely = dely_d
delz = delz_d
rk_stage = rk_stage_d
elem = elem_d

End Subroutine deviceToHost
!!==================================================================================

Subroutine nested1()
implicit none

integer :: ie,je,ke,i,j,k,eq,nlx,nly,nlz

    nlx = size(vars,1)
    nly = size(vars,2)
    nlz = size(vars,3)

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
do eq = 1, neq
do k = 1, nx
do j = 1, nx
do i = 1, nx

           rk_stage(ie,je,ke)%svec(i,j,k,eq) = vars(ie,je,ke)%vec(i,j,k,eq) + dt*elem(ie,je,ke)%rhs(i,j,k,eq)

        end do
      end do
    end do
end do
end do
end do
end do

End Subroutine nested1

!!==================================================================================
attributes(global) Subroutine nested1_cuda()
implicit none

integer :: ie,je,ke,i,j,k,eq,nlx,nly,nlz

!    nlx = size(vars,1)
!    nly = size(vars,2)
!    nlz = size(vars,3)

i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z
 !   do ke = 1, nlz
 !     do je = 1, nly
 !       do ie = 1, nlx
do eq = 1, neq
!do k = 1, nx
!do j = 1, nx
!do i = 1, nx

           rk_stage_d(ie,je,ke)%svec(i,j,k,eq) = vars_d(ie,je,ke)%vec(i,j,k,eq) + dt*elem_d(ie,je,ke)%rhs(i,j,k,eq)

!        end do
!      end do
!    end do
end do
!end do
!end do
!end do

End Subroutine nested1_cuda
!!==================================================================================

Subroutine nested2()
implicit none

integer :: ie,je,ke,i,j,k,eq,nlx,nly,nlz

    nlx = size(vars,1)
    nly = size(vars,2)
    nlz = size(vars,3)

do ke = 1, nlz
     do je = 1, nly
       do ie = 1, nlx
do eq = 1, neq
do k = 1, nx
do j = 1, nx
do i = 1, nx
          rk_stage(ie,je,ke)%svec(i,j,k,eq) = (rk_stage(ie,je,ke)%svec(i,j,k,eq)  +    & 
                         3.0D0*vars(ie,je,ke)%vec(i,j,k,eq) + dt*elem(ie,je,ke)%rhs(i,j,k,eq) ) *0.25D0

       end do
     end do
   end do
end do
end do
end do
end do

End Subroutine nested2

!!==================================================================================

attributes(global) Subroutine nested2_cuda()
implicit none

integer :: ie,je,ke,i,j,k,eq,nlx,nly,nlz

!    nlx = size(vars,1)
!    nly = size(vars,2)
!    nlz = size(vars,3)
i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z
!do ke = 1, nlz
!     do je = 1, nly
!       do ie = 1, nlx
do eq = 1, neq
!do k = 1, nx
!do j = 1, nx
!do i = 1, nx
          rk_stage_d(ie,je,ke)%svec(i,j,k,eq) = (rk_stage_d(ie,je,ke)%svec(i,j,k,eq)  +    & 
                         3.0D0*vars_d(ie,je,ke)%vec(i,j,k,eq) + dt*elem_d(ie,je,ke)%rhs(i,j,k,eq) ) *0.25D0

!       end do
!     end do
!   end do
!end do
!end do
!end do
end do

End Subroutine nested2_cuda
!!==================================================================================

Subroutine nested3()
implicit none

integer :: ie,je,ke,i,j,k,eq,nlx,nly,nlz

    nlx = size(vars,1)
    nly = size(vars,2)
    nlz = size(vars,3)

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
do eq = 1, neq
do k = 1, nx
do j = 1, nx
do i = 1, nx

          vars(ie,je,ke)%vec(i,j,k,eq) = (vars(ie,je,ke)%vec(i,j,k,eq) +       & 
                     2.0D0*(rk_stage(ie,je,ke)%svec(i,j,k,eq) + dt*elem(ie,je,ke)%rhs(i,j,k,eq)) ) /3.0D0

          rk_stage(ie,je,ke)%svec(i,j,k,eq) = vars(ie,je,ke)%vec(i,j,k,eq)

        end do
      end do
    end do
end do
end do
end do
end do

End Subroutine nested3

!!==================================================================================

attributes(global) Subroutine nested3_cuda()
implicit none

integer :: ie,je,ke,i,j,k,eq,nlx,nly,nlz

!    nlx = size(vars,1)
!    nly = size(vars,2)
!    nlz = size(vars,3)
i = threadIdx%x
j = threadIdx%y
k = threadIdx%z
ie = blockIdx%x
je = blockIdx%y
ke = blockIdx%z
!    do ke = 1, nlz
!      do je = 1, nly
!        do ie = 1, nlx
do eq = 1, neq
!do k = 1, nx
!do j = 1, nx
!do i = 1, nx

          vars_d(ie,je,ke)%vec(i,j,k,eq) = (vars_d(ie,je,ke)%vec(i,j,k,eq) +       & 
                     2.0D0*(rk_stage_d(ie,je,ke)%svec(i,j,k,eq) + dt*elem_d(ie,je,ke)%rhs(i,j,k,eq)) ) /3.0D0

          rk_stage_d(ie,je,ke)%svec(i,j,k,eq) = vars_d(ie,je,ke)%vec(i,j,k,eq)

!        end do
!      end do
!    end do
end do
!end do
!end do
!end do

End Subroutine nested3_cuda

!!==================================================================================
End module timestepper_mod

