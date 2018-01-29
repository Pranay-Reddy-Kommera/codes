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
  use euler_flux3d_mod, only: Make_Euler3d_Flux, Recover_State_Vars, Store_Hydrostatic_edgevars,nh_flux_maker
  use dg3d_rhs_mod, only: Compute_dgnh_rhs 
  use vertical_split_mod, only: Compute_column_rhs,Compute_SplitDG_rhs 
  use mpi_mod, only: send_data, recv_data

  implicit none
  private
  public :: ssp_rk3 

  ! When running in parallel, the number of elements owned by this rank
  ! is not known a priori, so we must dynamiCally allocate the storage:
  !type(element_t), dimension(:,:,:), allocatable :: temp1, temp2, rk_stage 
! type(element_t), dimension(:,:,:), allocatable :: rk_stage 
!  type(surface_t), dimension(:,:,:), allocatable :: edgevars
!  type(rkvec_t), dimension(:,:,:), allocatable :: rk_stage

contains
!!==================================================================================
!  Subroutine ssp_rk3(t, dt, vars)
  Subroutine ssp_rk3(t, dt)
    Implicit none 
    real(kind=real_kind), intent(inout) :: t
    real(kind=real_kind), intent(in) :: dt
!    type(element_t), intent(inout), dimension(:,:,:) :: vars

    real(kind=real_kind), dimension(nx,nx,nx,neq) :: rhs, vec 
    real(kind=real_kind), dimension(nx,nx,nx) ::  the,rho, rhb,thb,rtb
    integer :: nlx, nly, nlz, eq 
    integer :: ie, je, ke, i,j,k,l

    nlx = size(vars_rho,4)
    nly = size(vars_rho,5)
    nlz = size(vars_rho,6)

!write (*,*) 'nlx is ',nlx,'nly is ',nly,'nlz is ',nlz

!    ! Allocate the temp storage, this only happens on the FIRST Call
!    if (.not. allocated(rk_stage)) allocate(rk_stage(nlx,nly,nlz))
!    if (.not. allocated(edgevars)) allocate(edgevars(nlx,nly,nlz))

    ! Communicate basic NH varibels vars%(rho,the,u,v,w), and copy the edge values 
!    Call prepare_edgevars()
        Call extract_edgevars()
!       Call commuicate_edgevars()
        Call pack_edgevars()
        Call communicate_edgevars_work()
        Call send_data()
        Call recv_data()
        Call unpack_edgevars()

          if (t == 0.0D0) then !one time initialization 
            Call Store_Hydrostatic_edgevars()
          endif 

    !----------------
    ! RK stage 1
    !----------------
   Call nh_flux_maker()
            Call Make_Euler3d_Flux()
            Call Compute_dgnh_rhs(nlx,nly,nlz) 

!    if (use_dim_split_scheme) Call Compute_column_rhs()

!$acc parallel num_workers(16) vector_length(16)
!$acc loop gang worker vector collapse(7)

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
do eq = 1, neq
do k = 1, nx
do j = 1, nx
do i = 1, nx

           rk_stage_svec(i,j,k,eq,ie,je,ke) = vars_vec(i,j,k,eq,ie,je,ke) + dt*elem_rhs(i,j,k,eq,ie,je,ke)

        end do
      end do
    end do
end do
end do
end do
end do

!$acc end parallel

    !Compute/update  NH variables from the state vector 
     Call Recover_State_Vars() 

    !----------------
    ! RK stage 2
    !----------------

!   Call prepare_edgevars()
        Call extract_edgevars()
!       Call communicate_edgevars()
        Call pack_edgevars()
        Call communicate_edgevars_work()
        Call send_data()
        Call recv_data()
        Call unpack_edgevars()

        Call nh_flux_maker()
            Call Make_Euler3d_Flux()
            Call Compute_dgnh_rhs(nlx,nly,nlz)

!   if (use_dim_split_scheme) Call Compute_column_rhs()

!$acc parallel num_workers(16) vector_length(16)
!$acc loop gang worker vector collapse(7)
   do ke = 1, nlz
     do je = 1, nly
       do ie = 1, nlx
do eq = 1, neq
do k = 1, nx
do j = 1, nx
do i = 1, nx

          rk_stage_svec(i,j,k,eq,ie,je,ke) = (rk_stage_svec(i,j,k,eq,ie,je,ke)  +    & 
                         3.0D0*vars_vec(i,j,k,eq,ie,je,ke) + dt*elem_rhs(i,j,k,eq,ie,je,ke) ) *0.25D0         
       end do
     end do
   end do
end do
end do
end do
end do

!$acc end parallel

   Call Recover_State_Vars() 

!   !----------------
!   ! RK stage 3
!   !----------------

!   Call prepare_edgevars()
        Call extract_edgevars()
!       Call communicate_edgevars()
        Call pack_edgevars()
        Call communicate_edgevars_work()
        Call send_data()
        Call recv_data()
        Call unpack_edgevars()

        Call nh_flux_maker()
            Call Make_Euler3d_Flux()
            Call Compute_dgnh_rhs(nlx,nly,nlz)

!    if (use_dim_split_scheme) Call Compute_column_rhs()

!$acc parallel num_workers(16) vector_length(16)
!$acc loop gang worker vector collapse(7)
    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
do eq = 1, neq
do k = 1, nx
do j = 1, nx
do i = 1, nx

          vars_vec(i,j,k,eq,ie,je,ke) = (vars_vec(i,j,k,eq,ie,je,ke) +       & 
                     2.0D0*(rk_stage_svec(i,j,k,eq,ie,je,ke) + dt*elem_rhs(i,j,k,eq,ie,je,ke)) ) /3.0D0 
          rk_stage_svec(i,j,k,eq,ie,je,ke) = vars_vec(i,j,k,eq,ie,je,ke)

        end do
      end do
    end do
end do
end do
end do
end do
!$acc end parallel

    Call Recover_State_Vars() 

    ! update global time
    t = t + dt

  end subroutine ssp_rk3 

!!==================================================================================

End module timestepper_mod

