!-------------------------------------------------------------------
! Original code by R.Nair (date?)
! Changes 06/2016 by Francois Hebert, SIParCS/Cornell
!
! Some basic types to hold data associated with a DG element:
! * grid x,y,z
! * metric (for Gal-Chen & Somerville type transformation)
! * evolved variables in volume
! * evolved variables on boundary
!-------------------------------------------------------------------

module element_mod
  use basic_mod
  use openacc

  implicit none
  public

!  type :: element_t
!    sequence
!  !Basic thermodynamic variables 
!    real(kind=real_kind) :: rho(nx,nx,nx)
!    real(kind=real_kind) :: the(nx,nx,nx)
!    real(kind=real_kind) :: prs(nx,nx,nx)
!    real(kind=real_kind) :: spd(nx,nx,nx)
!    real(kind=real_kind) :: psi(nx,nx,nx)
!
!  !Velocity  variables 
!    real(kind=real_kind) :: u(nx,nx,nx)
!    real(kind=real_kind) :: v(nx,nx,nx)
!    real(kind=real_kind) :: w(nx,nx,nx)
!    real(kind=real_kind) :: wt(nx,nx,nx)
!
!   !Perturbed variables 
!    real(kind=real_kind) :: prp(nx,nx,nx)
!    real(kind=real_kind) :: rhp(nx,nx,nx)
!    real(kind=real_kind) :: thp(nx,nx,nx)
! 
!   !State and flux vectors for "neq" system 
!    real(kind=real_kind) :: vec(nx,nx,nx,neq)
!    real(kind=real_kind) :: flx(nx,nx,nx,neq)
!    real(kind=real_kind) :: fly(nx,nx,nx,neq)
!    real(kind=real_kind) :: flz(nx,nx,nx,neq)
!    real(kind=real_kind) :: src(nx,nx,nx,neq)
!
!   !source term from vertical 1D solver 
!    real(kind=real_kind) :: zsrc(nx,nx,nx,neq)
!
!   !Time frozen variables 
!    real(kind=real_kind) :: rhb(nx,nx,nx)
!    real(kind=real_kind) :: thb(nx,nx,nx)
!    real(kind=real_kind) :: rtb(nx,nx,nx)
!    real(kind=real_kind) :: prb(nx,nx,nx)
!    real(kind=real_kind) :: fcori(nx,nx,nx)
!  end type

!!!!!!!!!!!!!!!!!!!    Allocatable Arrays
!vars variable
real(kind=real_kind), allocatable :: vars_rho(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_the(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_prs(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_spd(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_psi(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_u(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_v(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_w(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_wt(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_prp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_rhp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_thp(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_vec(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_flx(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_fly(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_flz(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_src(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_zsrc(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_rhb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_thb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_rtb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_prb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_fcori(:,:,:,:,:,:)

!$acc declare create(vars_rho, vars_the, vars_prs, vars_spd, vars_psi, &
!$acc vars_u, vars_v, vars_w, vars_wt, vars_prp, vars_rhp, vars_thp, &
!$acc vars_vec, vars_flx, vars_fly, vars_flz, vars_src, vars_zsrc, &
!$acc vars_rhb, vars_thb, vars_rtb, vars_prb, vars_fcori)

!vars_init
real(kind=real_kind), allocatable :: vars_init_rho(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_the(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_prs(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_spd(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_psi(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_init_u(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_v(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_w(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_wt(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_init_prp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_rhp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_thp(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_init_vec(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_flx(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_fly(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_flz(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_src(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_init_zsrc(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_init_rhb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_thb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_rtb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_prb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_init_fcori(:,:,:,:,:,:)

!vars_phys
real(kind=real_kind), allocatable :: vars_phys_rho(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_the(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_prs(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_spd(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_psi(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_phys_u(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_v(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_w(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_wt(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_phys_prp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_rhp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_thp(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_phys_vec(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_flx(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_fly(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_flz(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_src(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_phys_zsrc(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: vars_phys_rhb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_thb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_rtb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_prb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: vars_phys_fcori(:,:,:,:,:,:)

!global_vars variable
real(kind=real_kind), allocatable :: global_vars_rho(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_the(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_prs(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_spd(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_psi(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_u(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_v(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_w(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_wt(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_prp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_rhp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_thp(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_vec(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_flx(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_fly(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_flz(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_src(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_zsrc(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_rhb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_thb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_rtb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_prb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_fcori(:,:,:,:,:,:)

!vars_init
real(kind=real_kind), allocatable :: global_vars_init_rho(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_the(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_prs(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_spd(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_psi(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_init_u(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_v(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_w(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_wt(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_init_prp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_rhp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_thp(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_init_vec(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_flx(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_fly(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_flz(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_src(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_init_zsrc(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_init_rhb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_thb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_rtb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_prb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_init_fcori(:,:,:,:,:,:)

!vars_phys
real(kind=real_kind), allocatable :: global_vars_phys_rho(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_the(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_prs(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_spd(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_psi(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_phys_u(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_v(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_w(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_wt(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_phys_prp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_rhp(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_thp(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_phys_vec(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_flx(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_fly(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_flz(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_src(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_phys_zsrc(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: global_vars_phys_rhb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_thb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_rtb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_prb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: global_vars_phys_fcori(:,:,:,:,:,:)

!  ! metric data for Gal-Chen Somerville type coordinate transform
!  type :: metric_t
!    sequence
!    real(kind=real_kind) :: mtn(nx,nx)     ! ground height as function of x,y
!    real(kind=real_kind) :: sg(nx,nx)      ! sqrt(det(g)) as function of x,y
!    real(kind=real_kind) :: sg13(nx,nx,nx) ! sqrt(det(g)) * dzeta/dx -- for computing wt
!    real(kind=real_kind) :: sg23(nx,nx,nx) ! sqrt(det(g)) * dzeta/dy -- for computing wt
!    real(kind=real_kind) :: ght(nx,nx,nx)  ! zeta(x,y,z)
!  end type

!!!!!!!!!!!!!!!   Allocatable Arrays
real(kind=real_kind), allocatable :: metric_mtn(:,:,:,:,:)
real(kind=real_kind), allocatable :: metric_sg(:,:,:,:,:)
real(kind=real_kind), allocatable :: metric_sg13(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: metric_sg23(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: metric_ght(:,:,:,:,:,:)

!$acc declare create(metric_mtn, metric_sg, metric_sg13, metric_sg23, metric_ght)

!  type :: grid_t
!    sequence
!    real(kind=real_kind) :: x(nx,nx,nx)
!    real(kind=real_kind) :: y(nx,nx,nx)
!    real(kind=real_kind) :: z(nx,nx,nx)
!  end type

!!!!!!!!!!!!!!!!!   Allocatable Arrays
real(kind=real_kind), allocatable :: grid_x(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: grid_y(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: grid_z(:,:,:,:,:,:)

!$acc declare create(grid_x, grid_y, grid_z)

!  type :: rkvec_t
!    real(kind=real_kind) :: svec(nx,nx,nx,neq)
!    real(kind=real_kind) :: svec0(nx,nx,nx,neq)
!  end type

!!!!!!!!!!!!!!!!!!  Allocatable Arrays
real(kind=real_kind), allocatable :: rk_stage_svec(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: rk_stage_svec0(:,:,:,:,:,:,:)

!$acc declare create(rk_stage_svec, rk_stage_svec0)

!  type :: edge_t
!    sequence
!    real(kind=real_kind) :: rho_int(nx,nx)
!    real(kind=real_kind) :: rho_ext(nx,nx)
!    real(kind=real_kind) :: the_int(nx,nx)
!    real(kind=real_kind) :: the_ext(nx,nx)
!    real(kind=real_kind) :: prp_int(nx,nx)
!    real(kind=real_kind) :: prp_ext(nx,nx)
!    real(kind=real_kind) :: spd_int(nx,nx)
!    real(kind=real_kind) :: spd_ext(nx,nx)
!
!    real(kind=real_kind) :: u_int(nx,nx)
!    real(kind=real_kind) :: u_ext(nx,nx)
!    real(kind=real_kind) :: v_int(nx,nx)
!    real(kind=real_kind) :: v_ext(nx,nx)
!    real(kind=real_kind) :: w_int(nx,nx)
!    real(kind=real_kind) :: w_ext(nx,nx)
!
!  !Surface fluxes 
!    real(kind=real_kind) :: vec_int(nx,nx,neq)
!    real(kind=real_kind) :: vec_ext(nx,nx,neq)
!    real(kind=real_kind) :: normflux_int(nx,nx,neq)
!    real(kind=real_kind) :: normflux_ext(nx,nx,neq)
!
!   !Time frozen variables for edge surfaces  
!    real(kind=real_kind) :: rtb(nx,nx)
!    real(kind=real_kind) :: rhb(nx,nx)
!    real(kind=real_kind) :: thb(nx,nx)
!    real(kind=real_kind) :: prb(nx,nx)
!  end type
!
!  type :: surface_t
!    sequence
!    type(edge_t) :: edges(6)
!  end type

!!!!!!!!!!!!!!  Allocatable Arrays
real(kind=real_kind), allocatable :: edgevars_edges_rho_int(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_rho_ext(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_the_int(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_the_ext(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_prp_int(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_prp_ext(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_spd_int(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_spd_ext(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: edgevars_edges_u_int(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_u_ext(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_v_int(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_v_ext(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_w_int(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_w_ext(:,:,:,:,:,:)

real(kind=real_kind), allocatable :: edgevars_edges_vec_int(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_vec_ext(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_normflux_int(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_normflux_ext(:,:,:,:,:,:,:)

real(kind=real_kind), allocatable :: edgevars_edges_rtb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_rhb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_thb(:,:,:,:,:,:)
real(kind=real_kind), allocatable :: edgevars_edges_prb(:,:,:,:,:,:)

!$acc declare create(edgevars_edges_rho_int, edgevars_edges_rho_ext, edgevars_edges_the_int, &
!$acc edgevars_edges_the_ext, edgevars_edges_prp_int, edgevars_edges_prp_ext, &
!$acc edgevars_edges_spd_int, edgevars_edges_spd_ext, edgevars_edges_u_int, &
!$acc edgevars_edges_u_ext, edgevars_edges_v_int, edgevars_edges_v_ext, &
!$acc edgevars_edges_w_int, edgevars_edges_w_ext, edgevars_edges_vec_int, &
!$acc edgevars_edges_vec_ext, edgevars_edges_normflux_int, edgevars_edges_normflux_ext, &
!$acc edgevars_edges_rtb, edgevars_edges_rhb, edgevars_edges_thb, &
!$acc edgevars_edges_prb)

  ! Used in the vertical-split DG scheme
!  type :: column_t
!    sequence
!!Basic thermodynamic variables 
!    real(kind=real_kind) :: rho(nx,nez)
!    real(kind=real_kind) :: the(nx,nez)
!    real(kind=real_kind) :: prs(nx,nez)
!    real(kind=real_kind) :: spd(nx,nez)
!    real(kind=real_kind) :: psi(nx,nez)
!
!  !Velocity  variables 
!    real(kind=real_kind) :: u(nx,nez)
!    real(kind=real_kind) :: v(nx,nez)
!    real(kind=real_kind) :: w(nx,nez)
!
!   !Perturbed variables 
!    real(kind=real_kind) :: prp(nx,nez)
!    real(kind=real_kind) :: rhp(nx,nez)
!    real(kind=real_kind) :: thp(nx,nez)
!
!   !State and flux vectors for "neq" system 
!    real(kind=real_kind) :: vec(nx,nez,neq)
!    real(kind=real_kind) :: flx(nx,nez,neq)
!    real(kind=real_kind) :: fly(nx,nez,neq)
!    real(kind=real_kind) :: flz(nx,nez,neq)
!    real(kind=real_kind) :: rhs(nx,nez,neq)
!
!  end type

real(kind=real_kind), allocatable :: velem_rho(:,:)
real(kind=real_kind), allocatable :: velem_the(:,:)
real(kind=real_kind), allocatable :: velem_prs(:,:)
real(kind=real_kind), allocatable :: velem_spd(:,:)
real(kind=real_kind), allocatable :: velem_psi(:,:)

real(kind=real_kind), allocatable :: velem_u(:,:)
real(kind=real_kind), allocatable :: velem_v(:,:)
real(kind=real_kind), allocatable :: velem_w(:,:)

real(kind=real_kind), allocatable :: velem_prp(:,:)
real(kind=real_kind), allocatable :: velem_rhp(:,:)
real(kind=real_kind), allocatable :: velem_thp(:,:)

real(kind=real_kind), allocatable :: velem_vec(:,:,:)
real(kind=real_kind), allocatable :: velem_flx(:,:,:)
real(kind=real_kind), allocatable :: velem_fly(:,:,:)
real(kind=real_kind), allocatable :: velem_flz(:,:,:)
real(kind=real_kind), allocatable :: velem_rhs(:,:,:)
real(kind=real_kind), allocatable :: elem_rhs(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: fjmax(:,:,:,:)
real(kind=real_kind), allocatable :: numflux(:,:,:,:,:,:,:)
real(kind=real_kind), allocatable :: divf(:,:,:,:,:,:,:)

!$acc declare create(velem_rho, velem_the, velem_prs, velem_spd, &
!$acc velem_psi, velem_u, velem_v, velem_w, velem_prp, velem_rhp, &
!$acc velem_thp, velem_vec, velem_flx, velem_fly, velem_flz, velem_rhs,&
!$acc elem_rhs,fjmax,numflux,divf)

contains

subroutine allocateData(x, y, z)
implicit none
integer, intent(in) :: x, y, z

!--------------- Allocate data
!vars
allocate(vars_rho(nx,nx,nx,x,y,z))
allocate(vars_the(nx,nx,nx,x,y,z))
allocate(vars_prs(nx,nx,nx,x,y,z))
allocate(vars_spd(nx,nx,nx,x,y,z))
allocate(vars_psi(nx,nx,nx,x,y,z))

allocate(vars_u(nx,nx,nx,x,y,z))
allocate(vars_v(nx,nx,nx,x,y,z))
allocate(vars_w(nx,nx,nx,x,y,z))
allocate(vars_wt(nx,nx,nx,x,y,z))

allocate(vars_prp(nx,nx,nx,x,y,z))
allocate(vars_rhp(nx,nx,nx,x,y,z))
allocate(vars_thp(nx,nx,nx,x,y,z))

allocate(vars_vec(nx,nx,nx,neq,x,y,z))
allocate(vars_flx(nx,nx,nx,neq,x,y,z))
allocate(vars_fly(nx,nx,nx,neq,x,y,z))
allocate(vars_flz(nx,nx,nx,neq,x,y,z))
allocate(vars_src(nx,nx,nx,neq,x,y,z))
allocate(elem_rhs(nx,nx,nx,neq,x,y,z))
allocate(fjmax(6,x,y,z))
allocate(numflux(nx,nx,neq,6,x,y,z))
allocate(divf(nx,nx,nx,neq,x,y,z))

allocate(vars_zsrc(nx,nx,nx,neq,x,y,z))

allocate(vars_rhb(nx,nx,nx,x,y,z))
allocate(vars_thb(nx,nx,nx,x,y,z))
allocate(vars_rtb(nx,nx,nx,x,y,z))
allocate(vars_prb(nx,nx,nx,x,y,z))
allocate(vars_fcori(nx,nx,nx,x,y,z))

!vars_init
allocate(vars_init_rho(nx,nx,nx,x,y,z))
allocate(vars_init_the(nx,nx,nx,x,y,z))
allocate(vars_init_prs(nx,nx,nx,x,y,z))
allocate(vars_init_spd(nx,nx,nx,x,y,z))
allocate(vars_init_psi(nx,nx,nx,x,y,z))

allocate(vars_init_u(nx,nx,nx,x,y,z))
allocate(vars_init_v(nx,nx,nx,x,y,z))
allocate(vars_init_w(nx,nx,nx,x,y,z))
allocate(vars_init_wt(nx,nx,nx,x,y,z))

allocate(vars_init_prp(nx,nx,nx,x,y,z))
allocate(vars_init_rhp(nx,nx,nx,x,y,z))
allocate(vars_init_thp(nx,nx,nx,x,y,z))

allocate(vars_init_vec(nx,nx,nx,neq,x,y,z))
allocate(vars_init_flx(nx,nx,nx,neq,x,y,z))
allocate(vars_init_fly(nx,nx,nx,neq,x,y,z))
allocate(vars_init_flz(nx,nx,nx,neq,x,y,z))
allocate(vars_init_src(nx,nx,nx,neq,x,y,z))

allocate(vars_init_zsrc(nx,nx,nx,neq,x,y,z))

allocate(vars_init_rhb(nx,nx,nx,x,y,z))
allocate(vars_init_thb(nx,nx,nx,x,y,z))
allocate(vars_init_rtb(nx,nx,nx,x,y,z))
allocate(vars_init_prb(nx,nx,nx,x,y,z))
allocate(vars_init_fcori(nx,nx,nx,x,y,z))

!vars_phys
allocate(vars_phys_rho(nx,nx,nx,x,y,z))
allocate(vars_phys_the(nx,nx,nx,x,y,z))
allocate(vars_phys_prs(nx,nx,nx,x,y,z))
allocate(vars_phys_spd(nx,nx,nx,x,y,z))
allocate(vars_phys_psi(nx,nx,nx,x,y,z))

allocate(vars_phys_u(nx,nx,nx,x,y,z))
allocate(vars_phys_v(nx,nx,nx,x,y,z))
allocate(vars_phys_w(nx,nx,nx,x,y,z))
allocate(vars_phys_wt(nx,nx,nx,x,y,z))

allocate(vars_phys_prp(nx,nx,nx,x,y,z))
allocate(vars_phys_rhp(nx,nx,nx,x,y,z))
allocate(vars_phys_thp(nx,nx,nx,x,y,z))

allocate(vars_phys_vec(nx,nx,nx,neq,x,y,z))
allocate(vars_phys_flx(nx,nx,nx,neq,x,y,z))
allocate(vars_phys_fly(nx,nx,nx,neq,x,y,z))
allocate(vars_phys_flz(nx,nx,nx,neq,x,y,z))
allocate(vars_phys_src(nx,nx,nx,neq,x,y,z))

allocate(vars_phys_zsrc(nx,nx,nx,neq,x,y,z))

allocate(vars_phys_rhb(nx,nx,nx,x,y,z))
allocate(vars_phys_thb(nx,nx,nx,x,y,z))
allocate(vars_phys_rtb(nx,nx,nx,x,y,z))
allocate(vars_phys_prb(nx,nx,nx,x,y,z))
allocate(vars_phys_fcori(nx,nx,nx,x,y,z))

!grid
allocate(grid_x(nx,nx,nx,x,y,z))
allocate(grid_y(nx,nx,nx,x,y,z))
allocate(grid_z(nx,nx,nx,x,y,z))

!rk_stage
allocate(rk_stage_svec(nx,nx,nx,neq,x,y,z))
allocate(rk_stage_svec0(nx,nx,nx,neq,x,y,z))

!edgevars
allocate(edgevars_edges_rho_int(nx,nx,6,x,y,z))
allocate(edgevars_edges_rho_ext(nx,nx,6,x,y,z))
allocate(edgevars_edges_the_int(nx,nx,6,x,y,z))
allocate(edgevars_edges_the_ext(nx,nx,6,x,y,z))
allocate(edgevars_edges_prp_int(nx,nx,6,x,y,z))
allocate(edgevars_edges_prp_ext(nx,nx,6,x,y,z))
allocate(edgevars_edges_spd_int(nx,nx,6,x,y,z))
allocate(edgevars_edges_spd_ext(nx,nx,6,x,y,z))

allocate(edgevars_edges_u_int(nx,nx,6,x,y,z))
allocate(edgevars_edges_u_ext(nx,nx,6,x,y,z))
allocate(edgevars_edges_v_int(nx,nx,6,x,y,z))
allocate(edgevars_edges_v_ext(nx,nx,6,x,y,z))
allocate(edgevars_edges_w_int(nx,nx,6,x,y,z))
allocate(edgevars_edges_w_ext(nx,nx,6,x,y,z))

allocate(edgevars_edges_vec_int(nx,nx,neq,6,x,y,z))
allocate(edgevars_edges_vec_ext(nx,nx,neq,6,x,y,z))
allocate(edgevars_edges_normflux_int(nx,nx,neq,6,x,y,z))
allocate(edgevars_edges_normflux_ext(nx,nx,neq,6,x,y,z))

allocate(edgevars_edges_rtb(nx,nx,6,x,y,z))
allocate(edgevars_edges_rhb(nx,nx,6,x,y,z))
allocate(edgevars_edges_thb(nx,nx,6,x,y,z))
allocate(edgevars_edges_prb(nx,nx,6,x,y,z))

!metric
allocate(metric_mtn(nx,nx,x,y,z))
allocate(metric_sg(nx,nx,x,y,z))
allocate(metric_sg13(nx,nx,nx,x,y,z))
allocate(metric_sg23(nx,nx,nx,x,y,z))
allocate(metric_ght(nx,nx,nx,x,y,z))

!velem
allocate(velem_rho(nx,nez))
allocate(velem_the(nx,nez))
allocate(velem_prs(nx,nez))
allocate(velem_spd(nx,nez))
allocate(velem_psi(nx,nez))

allocate(velem_u(nx,nez))
allocate(velem_v(nx,nez))
allocate(velem_w(nx,nez))

allocate(velem_prp(nx,nez))
allocate(velem_rhp(nx,nez))
allocate(velem_thp(nx,nez))

allocate(velem_vec(nx,nez,neq))
allocate(velem_flx(nx,nez,neq))
allocate(velem_fly(nx,nez,neq))
allocate(velem_flz(nx,nez,neq))
allocate(velem_rhs(nx,nez,neq))

end subroutine allocateData

subroutine allocateGlobalData(x,y,z)
implicit none
integer, intent(in) :: x, y, z

!global_vars
allocate(global_vars_rho(nx,nx,nx,x,y,z))
allocate(global_vars_the(nx,nx,nx,x,y,z))
allocate(global_vars_prs(nx,nx,nx,x,y,z))
allocate(global_vars_spd(nx,nx,nx,x,y,z))
allocate(global_vars_psi(nx,nx,nx,x,y,z))

allocate(global_vars_u(nx,nx,nx,x,y,z))
allocate(global_vars_v(nx,nx,nx,x,y,z))
allocate(global_vars_w(nx,nx,nx,x,y,z))
allocate(global_vars_wt(nx,nx,nx,x,y,z))

allocate(global_vars_prp(nx,nx,nx,x,y,z))
allocate(global_vars_rhp(nx,nx,nx,x,y,z))
allocate(global_vars_thp(nx,nx,nx,x,y,z))

allocate(global_vars_vec(nx,nx,nx,neq,x,y,z))
allocate(global_vars_flx(nx,nx,nx,neq,x,y,z))
allocate(global_vars_fly(nx,nx,nx,neq,x,y,z))
allocate(global_vars_flz(nx,nx,nx,neq,x,y,z))
allocate(global_vars_src(nx,nx,nx,neq,x,y,z))

allocate(global_vars_zsrc(nx,nx,nx,neq,x,y,z))

allocate(global_vars_rhb(nx,nx,nx,x,y,z))
allocate(global_vars_thb(nx,nx,nx,x,y,z))
allocate(global_vars_rtb(nx,nx,nx,x,y,z))
allocate(global_vars_prb(nx,nx,nx,x,y,z))
allocate(global_vars_fcori(nx,nx,nx,x,y,z))

!global_vars_init
allocate(global_vars_init_rho(nx,nx,nx,x,y,z))
allocate(global_vars_init_the(nx,nx,nx,x,y,z))
allocate(global_vars_init_prs(nx,nx,nx,x,y,z))
allocate(global_vars_init_spd(nx,nx,nx,x,y,z))
allocate(global_vars_init_psi(nx,nx,nx,x,y,z))

allocate(global_vars_init_u(nx,nx,nx,x,y,z))
allocate(global_vars_init_v(nx,nx,nx,x,y,z))
allocate(global_vars_init_w(nx,nx,nx,x,y,z))
allocate(global_vars_init_wt(nx,nx,nx,x,y,z))

allocate(global_vars_init_prp(nx,nx,nx,x,y,z))
allocate(global_vars_init_rhp(nx,nx,nx,x,y,z))
allocate(global_vars_init_thp(nx,nx,nx,x,y,z))

allocate(global_vars_init_vec(nx,nx,nx,neq,x,y,z))
allocate(global_vars_init_flx(nx,nx,nx,neq,x,y,z))
allocate(global_vars_init_fly(nx,nx,nx,neq,x,y,z))
allocate(global_vars_init_flz(nx,nx,nx,neq,x,y,z))
allocate(global_vars_init_src(nx,nx,nx,neq,x,y,z))

allocate(global_vars_init_zsrc(nx,nx,nx,neq,x,y,z))

allocate(global_vars_init_rhb(nx,nx,nx,x,y,z))
allocate(global_vars_init_thb(nx,nx,nx,x,y,z))
allocate(global_vars_init_rtb(nx,nx,nx,x,y,z))
allocate(global_vars_init_prb(nx,nx,nx,x,y,z))
allocate(global_vars_init_fcori(nx,nx,nx,x,y,z))

!gloabl_vars_phys
allocate(global_vars_phys_rho(nx,nx,nx,x,y,z))
allocate(global_vars_phys_the(nx,nx,nx,x,y,z))
allocate(global_vars_phys_prs(nx,nx,nx,x,y,z))
allocate(global_vars_phys_spd(nx,nx,nx,x,y,z))
allocate(global_vars_phys_psi(nx,nx,nx,x,y,z))

allocate(global_vars_phys_u(nx,nx,nx,x,y,z))
allocate(global_vars_phys_v(nx,nx,nx,x,y,z))
allocate(global_vars_phys_w(nx,nx,nx,x,y,z))
allocate(global_vars_phys_wt(nx,nx,nx,x,y,z))

allocate(global_vars_phys_prp(nx,nx,nx,x,y,z))
allocate(global_vars_phys_rhp(nx,nx,nx,x,y,z))
allocate(global_vars_phys_thp(nx,nx,nx,x,y,z))

allocate(global_vars_phys_vec(nx,nx,nx,neq,x,y,z))
allocate(global_vars_phys_flx(nx,nx,nx,neq,x,y,z))
allocate(global_vars_phys_fly(nx,nx,nx,neq,x,y,z))
allocate(global_vars_phys_flz(nx,nx,nx,neq,x,y,z))
allocate(global_vars_phys_src(nx,nx,nx,neq,x,y,z))

allocate(global_vars_phys_zsrc(nx,nx,nx,neq,x,y,z))

allocate(global_vars_phys_rhb(nx,nx,nx,x,y,z))
allocate(global_vars_phys_thb(nx,nx,nx,x,y,z))
allocate(global_vars_phys_rtb(nx,nx,nx,x,y,z))
allocate(global_vars_phys_prb(nx,nx,nx,x,y,z))
allocate(global_vars_phys_fcori(nx,nx,nx,x,y,z))

end subroutine allocateGlobalData

subroutine deallocateData()
implicit none

if (allocated(vars_rho)) then

!--------------- Deallocate Data
!vars
deallocate(vars_rho)
deallocate(vars_the)
deallocate(vars_prs)
deallocate(vars_spd)
deallocate(vars_psi)

deallocate(vars_u)
deallocate(vars_v)
deallocate(vars_w)
deallocate(vars_wt)

deallocate(vars_prp)
deallocate(vars_rhp)
deallocate(vars_thp)

deallocate(vars_vec)
deallocate(vars_flx)
deallocate(vars_fly)
deallocate(vars_flz)
deallocate(vars_src)

deallocate(vars_zsrc)

deallocate(vars_rhb)
deallocate(vars_thb)
deallocate(vars_rtb)
deallocate(vars_prb)
deallocate(vars_fcori)

!vars_init
deallocate(vars_init_rho)
deallocate(vars_init_the)
deallocate(vars_init_prs)
deallocate(vars_init_spd)
deallocate(vars_init_psi)

deallocate(vars_init_u)
deallocate(vars_init_v)
deallocate(vars_init_w)
deallocate(vars_init_wt)

deallocate(vars_init_prp)
deallocate(vars_init_rhp)
deallocate(vars_init_thp)

deallocate(vars_init_vec)
deallocate(vars_init_flx)
deallocate(vars_init_fly)
deallocate(vars_init_flz)
deallocate(vars_init_src)

deallocate(vars_init_zsrc)

deallocate(vars_init_rhb)
deallocate(vars_init_thb)
deallocate(vars_init_rtb)
deallocate(vars_init_prb)
deallocate(vars_init_fcori)

!vars
deallocate(vars_phys_rho)
deallocate(vars_phys_the)
deallocate(vars_phys_prs)
deallocate(vars_phys_spd)
deallocate(vars_phys_psi)

deallocate(vars_phys_u)
deallocate(vars_phys_v)
deallocate(vars_phys_w)
deallocate(vars_phys_wt)

deallocate(vars_phys_prp)
deallocate(vars_phys_rhp)
deallocate(vars_phys_thp)

deallocate(vars_phys_vec)
deallocate(vars_phys_flx)
deallocate(vars_phys_fly)
deallocate(vars_phys_flz)
deallocate(vars_phys_src)

deallocate(vars_phys_zsrc)

deallocate(vars_phys_rhb)
deallocate(vars_phys_thb)
deallocate(vars_phys_rtb)
deallocate(vars_phys_prb)
deallocate(vars_phys_fcori)

!grid
deallocate(grid_x)
deallocate(grid_y)
deallocate(grid_z)

!rk_stage
deallocate(rk_stage_svec)
deallocate(rk_stage_svec0)

!edgevars
deallocate(edgevars_edges_rho_int)
deallocate(edgevars_edges_rho_ext)
deallocate(edgevars_edges_the_int)
deallocate(edgevars_edges_the_ext)
deallocate(edgevars_edges_prp_int)
deallocate(edgevars_edges_prp_ext)
deallocate(edgevars_edges_spd_int)
deallocate(edgevars_edges_spd_ext)

deallocate(edgevars_edges_u_int)
deallocate(edgevars_edges_u_ext)
deallocate(edgevars_edges_v_int)
deallocate(edgevars_edges_v_ext)
deallocate(edgevars_edges_w_int)
deallocate(edgevars_edges_w_ext)

deallocate(edgevars_edges_vec_int)
deallocate(edgevars_edges_vec_ext)
deallocate(edgevars_edges_normflux_int)
deallocate(edgevars_edges_normflux_ext)

deallocate(edgevars_edges_rtb)
deallocate(edgevars_edges_rhb)
deallocate(edgevars_edges_thb)
deallocate(edgevars_edges_prb)

!metric
deallocate(metric_mtn)
deallocate(metric_sg)
deallocate(metric_sg13)
deallocate(metric_sg23)
deallocate(metric_ght)

endif

if (allocated(global_vars_rho)) then

!global_vars
deallocate(global_vars_rho)
deallocate(global_vars_the)
deallocate(global_vars_prs)
deallocate(global_vars_spd)
deallocate(global_vars_psi)

deallocate(global_vars_u)
deallocate(global_vars_v)
deallocate(global_vars_w)
deallocate(global_vars_wt)

deallocate(global_vars_prp)
deallocate(global_vars_rhp)
deallocate(global_vars_thp)

deallocate(global_vars_vec)
deallocate(global_vars_flx)
deallocate(global_vars_fly)
deallocate(global_vars_flz)
deallocate(global_vars_src)

deallocate(global_vars_zsrc)

deallocate(global_vars_rhb)
deallocate(global_vars_thb)
deallocate(global_vars_rtb)
deallocate(global_vars_prb)
deallocate(global_vars_fcori)

!global_vars_init
deallocate(global_vars_init_rho)
deallocate(global_vars_init_the)
deallocate(global_vars_init_prs)
deallocate(global_vars_init_spd)
deallocate(global_vars_init_psi)

deallocate(global_vars_init_u)
deallocate(global_vars_init_v)
deallocate(global_vars_init_w)
deallocate(global_vars_init_wt)

deallocate(global_vars_init_prp)
deallocate(global_vars_init_rhp)
deallocate(global_vars_init_thp)

deallocate(global_vars_init_vec)
deallocate(global_vars_init_flx)
deallocate(global_vars_init_fly)
deallocate(global_vars_init_flz)
deallocate(global_vars_init_src)

deallocate(global_vars_init_zsrc)

deallocate(global_vars_init_rhb)
deallocate(global_vars_init_thb)
deallocate(global_vars_init_rtb)
deallocate(global_vars_init_prb)
deallocate(global_vars_init_fcori)

!gloabl_vars_phys
deallocate(global_vars_phys_rho)
deallocate(global_vars_phys_the)
deallocate(global_vars_phys_prs)
deallocate(global_vars_phys_spd)
deallocate(global_vars_phys_psi)

deallocate(global_vars_phys_u)
deallocate(global_vars_phys_v)
deallocate(global_vars_phys_w)
deallocate(global_vars_phys_wt)

deallocate(global_vars_phys_prp)
deallocate(global_vars_phys_rhp)
deallocate(global_vars_phys_thp)

deallocate(global_vars_phys_vec)
deallocate(global_vars_phys_flx)
deallocate(global_vars_phys_fly)
deallocate(global_vars_phys_flz)
deallocate(global_vars_phys_src)

deallocate(global_vars_phys_zsrc)

deallocate(global_vars_phys_rhb)
deallocate(global_vars_phys_thb)
deallocate(global_vars_phys_rtb)
deallocate(global_vars_phys_prb)
deallocate(global_vars_phys_fcori)

end if

end subroutine deallocateData

subroutine varsToVarsPhys()
implicit none

vars_phys_rho = vars_rho
vars_phys_the = vars_the
vars_phys_prs = vars_prs
vars_phys_spd = vars_spd
vars_phys_psi = vars_psi

vars_phys_u = vars_u
vars_phys_v = vars_v
vars_phys_w = vars_w
vars_phys_wt = vars_wt

vars_phys_prp = vars_prp
vars_phys_rhp = vars_rhp
vars_phys_thp = vars_thp

vars_phys_vec = vars_vec
vars_phys_flx = vars_flx
vars_phys_fly = vars_fly
vars_phys_flz = vars_flz
vars_phys_src = vars_src

vars_phys_zsrc = vars_zsrc

vars_phys_rhb = vars_rhb
vars_phys_thb = vars_thb
vars_phys_rtb = vars_rtb
vars_phys_prb = vars_prb
vars_phys_fcori = vars_fcori

end subroutine varsToVarsPhys

subroutine varsInitToVars()
implicit none

vars_rho = vars_init_rho
vars_the = vars_init_the
vars_prs = vars_init_prs
vars_spd = vars_init_spd
vars_psi = vars_init_psi

vars_u = vars_init_u
vars_v = vars_init_v
vars_w = vars_init_w
vars_wt = vars_init_wt

vars_prp = vars_init_prp
vars_rhp = vars_init_rhp
vars_thp = vars_init_thp

vars_vec = vars_init_vec
vars_flx = vars_init_flx
vars_fly = vars_init_fly
vars_flz = vars_init_flz
vars_src = vars_init_src

vars_zsrc = vars_init_zsrc

vars_rhb = vars_init_rhb
vars_thb = vars_init_thb
vars_rtb = vars_init_rtb
vars_prb = vars_init_prb
vars_fcori = vars_init_fcori


end subroutine varsInitToVars

end module element_mod
