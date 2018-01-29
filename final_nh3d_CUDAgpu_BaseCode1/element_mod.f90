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
  use cudafor

  implicit none
  public

  type :: element_t
    sequence
  !Basic thermodynamic variables 
    real(kind=real_kind) :: rho(nx,nx,nx)
    real(kind=real_kind) :: the(nx,nx,nx)
    real(kind=real_kind) :: prs(nx,nx,nx)
    real(kind=real_kind) :: spd(nx,nx,nx)
    real(kind=real_kind) :: psi(nx,nx,nx)

  !Velocity  variables 
    real(kind=real_kind) :: u(nx,nx,nx)
    real(kind=real_kind) :: v(nx,nx,nx)
    real(kind=real_kind) :: w(nx,nx,nx)
    real(kind=real_kind) :: wt(nx,nx,nx)

   !Perturbed variables 
    real(kind=real_kind) :: prp(nx,nx,nx)
    real(kind=real_kind) :: rhp(nx,nx,nx)
    real(kind=real_kind) :: thp(nx,nx,nx)
 
   !State and flux vectors for "neq" system 
    real(kind=real_kind) :: vec(nx,nx,nx,neq)
    real(kind=real_kind) :: flx(nx,nx,nx,neq)
    real(kind=real_kind) :: fly(nx,nx,nx,neq)
    real(kind=real_kind) :: flz(nx,nx,nx,neq)
    real(kind=real_kind) :: src(nx,nx,nx,neq)

   !source term from vertical 1D solver 
    real(kind=real_kind) :: zsrc(nx,nx,nx,neq)

   !Time frozen variables 
    real(kind=real_kind) :: rhb(nx,nx,nx)
    real(kind=real_kind) :: thb(nx,nx,nx)
    real(kind=real_kind) :: rtb(nx,nx,nx)
    real(kind=real_kind) :: prb(nx,nx,nx)
    real(kind=real_kind) :: fcori(nx,nx,nx)
  end type

  ! metric data for Gal-Chen Somerville type coordinate transform
  type :: metric_t
    sequence
    real(kind=real_kind) :: mtn(nx,nx)     ! ground height as function of x,y
    real(kind=real_kind) :: sg(nx,nx)      ! sqrt(det(g)) as function of x,y
    real(kind=real_kind) :: sg13(nx,nx,nx) ! sqrt(det(g)) * dzeta/dx -- for computing wt
    real(kind=real_kind) :: sg23(nx,nx,nx) ! sqrt(det(g)) * dzeta/dy -- for computing wt
    real(kind=real_kind) :: ght(nx,nx,nx)  ! zeta(x,y,z)
  end type

  type :: grid_t
    sequence
    real(kind=real_kind) :: x(nx,nx,nx)
    real(kind=real_kind) :: y(nx,nx,nx)
    real(kind=real_kind) :: z(nx,nx,nx)
  end type

  type :: rkvec_t
    real(kind=real_kind) :: svec(nx,nx,nx,neq)
    real(kind=real_kind) :: svec0(nx,nx,nx,neq)
  end type

  type :: edge_t
    sequence
    real(kind=real_kind) :: rho_int(nx,nx)
    real(kind=real_kind) :: rho_ext(nx,nx)
    real(kind=real_kind) :: the_int(nx,nx)
    real(kind=real_kind) :: the_ext(nx,nx)
    real(kind=real_kind) :: prp_int(nx,nx)
    real(kind=real_kind) :: prp_ext(nx,nx)
    real(kind=real_kind) :: spd_int(nx,nx)
    real(kind=real_kind) :: spd_ext(nx,nx)

    real(kind=real_kind) :: u_int(nx,nx)
    real(kind=real_kind) :: u_ext(nx,nx)
    real(kind=real_kind) :: v_int(nx,nx)
    real(kind=real_kind) :: v_ext(nx,nx)
    real(kind=real_kind) :: w_int(nx,nx)
    real(kind=real_kind) :: w_ext(nx,nx)

  !Surface fluxes 
    real(kind=real_kind) :: vec_int(nx,nx,neq)
    real(kind=real_kind) :: vec_ext(nx,nx,neq)
    real(kind=real_kind) :: normflux_int(nx,nx,neq)
    real(kind=real_kind) :: normflux_ext(nx,nx,neq)

   !Time frozen variables for edge surfaces  
    real(kind=real_kind) :: rtb(nx,nx)
    real(kind=real_kind) :: rhb(nx,nx)
    real(kind=real_kind) :: thb(nx,nx)
    real(kind=real_kind) :: prb(nx,nx)
  end type

  type :: surface_t
    sequence
    type(edge_t) :: edges(6)
  end type

  ! Used in the vertical-split DG scheme
  type :: column_t
    sequence
!Basic thermodynamic variables 
    real(kind=real_kind) :: rho(nx,nez)
    real(kind=real_kind) :: the(nx,nez)
    real(kind=real_kind) :: prs(nx,nez)
    real(kind=real_kind) :: spd(nx,nez)
    real(kind=real_kind) :: psi(nx,nez)

  !Velocity  variables 
    real(kind=real_kind) :: u(nx,nez)
    real(kind=real_kind) :: v(nx,nez)
    real(kind=real_kind) :: w(nx,nez)

   !Perturbed variables 
    real(kind=real_kind) :: prp(nx,nez)
    real(kind=real_kind) :: rhp(nx,nez)
    real(kind=real_kind) :: thp(nx,nez)

   !State and flux vectors for "neq" system 
    real(kind=real_kind) :: vec(nx,nez,neq)
    real(kind=real_kind) :: flx(nx,nez,neq)
    real(kind=real_kind) :: fly(nx,nez,neq)
    real(kind=real_kind) :: flz(nx,nez,neq)
    real(kind=real_kind) :: rhs(nx,nez,neq)

  end type

   ! elems
  type :: elem_s
    sequence
    real(kind=real_kind) :: rhs(nx,nx,nx,neq)
    real(kind=real_kind) :: fjmax(6)
    real(kind=real_kind) :: numflux(nx,nx,neq,6)
    real(kind=real_kind) :: divf(nx,nx,nx,neq)

  end type


  ! allocatable because MPI partition is determined at runtime
  type(grid_t), dimension(:,:,:), allocatable :: grid
  type(element_t), dimension(:,:,:), allocatable :: vars
  type(element_t), device, dimension(:,:,:), allocatable :: vars_d
  type(element_t), dimension(:,:,:), allocatable :: vars_init
  type(element_t), dimension(:,:,:), allocatable :: vars_phys

  ! allocatable because we only want global data on master rank
  type(element_t), dimension(:,:,:), allocatable :: global_vars
  type(element_t), dimension(:,:,:), allocatable :: global_vars_init
  type(element_t), dimension(:,:,:), allocatable :: global_vars_phys

  type(metric_t), dimension(:,:,:), allocatable :: metric
  type(rkvec_t), dimension(:,:,:), allocatable :: rk_stage
  type(rkvec_t), device, dimension(:,:,:), allocatable :: rk_stage_d
  type(surface_t), dimension(:,:,:), allocatable :: edgevars
  type(surface_t), device, dimension(:,:,:), allocatable :: edgevars_d
  type(column_t) :: velem
  type(elem_s), dimension(:,:,:), allocatable :: elem
  type(elem_s), device, dimension(:,:,:), allocatable :: elem_d


contains

subroutine allocateData(x, y, z)
implicit none
integer, intent(in) :: x, y, z

!---------------- Allocate Data
allocate(vars(x, y, z))
allocate(vars_d(x, y, z))
allocate(vars_init(x, y, z))
allocate(vars_phys(x, y, z))
allocate(grid(x, y, z))
allocate(rk_stage(x, y, z))
allocate(rk_stage_d(x, y, z))
allocate(edgevars(x, y, z))
allocate(edgevars_d(x, y, z))
allocate(metric(x, y, z))
allocate(elem(x, y, z))
allocate(elem_d(x, y, z))

end subroutine allocateData

subroutine allocateGlobalData(x, y, z)
implicit none
integer, intent(in) :: x,y,z

allocate(global_vars(x,y,z))
allocate(global_vars_init(x,y,z))
allocate(global_vars_phys(x,y,z))

end subroutine allocateGlobalData

subroutine deallocateData()
implicit none

if (allocated(vars)) then
deallocate(vars)
deallocate(vars_d)
deallocate(vars_init)
deallocate(vars_phys)
deallocate(grid)
deallocate(rk_stage)
deallocate(rk_stage_d)
deallocate(edgevars)
deallocate(edgevars_d)
deallocate(metric)
deallocate(elem)
deallocate(elem_d)
endif

if (allocated(global_vars)) then
deallocate(global_vars)
deallocate(global_vars_init)
deallocate(global_vars_phys)
endif

end subroutine deallocateData

end module element_mod
