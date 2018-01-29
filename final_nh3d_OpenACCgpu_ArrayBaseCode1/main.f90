!-------------------------------------------------------------------
! 09/04; 06/10 Sciparcs R.Nair NCAR/IMAGe
! Rewritten 07/2016 by Francois Hebert, SIParCS/Cornell
!
! The main driver for a 3D advection test.
!   Domain: 3D cartesian tensor product
!   Numerical method: Disc. Galerkin, nodal, weak, GLL grid
!   Timestepper: RK3 TVD/SSP
!-------------------------------------------------------------------

program main
  use basic_mod
  use element_mod
  use mpi_mod
  use grid_setup_mod
  use testcases_nh3d_mod
  use timestepper_mod
  use gather_solution_mod
  use output_mod
  use openacc

  implicit none

! Commented as they are declared in element_mod.f90
!  ! allocatable because MPI partition is determined at runtime
!  type(grid_t), dimension(:,:,:), allocatable :: grid
!  type(element_t), dimension(:,:,:), allocatable :: vars
!  type(element_t), dimension(:,:,:), allocatable :: vars_init
!  type(element_t), dimension(:,:,:), allocatable :: vars_phys
!
!  ! allocatable because we only want global data on master rank
!  type(element_t), dimension(:,:,:), allocatable :: global_vars
!  type(element_t), dimension(:,:,:), allocatable :: global_vars_init
!  type(element_t), dimension(:,:,:), allocatable :: global_vars_phys

  real(kind=real_kind) :: t, timer_start, timer_end
  integer :: nlx, nly, nlz
  integer :: nsteps, itn, frameid, nfreq
integer*8 :: start_clock, stop_clock, rate_clock

  !-----------------------------------------------------------------
  ! setup work:
  ! init MPI, build grid and initial data, save initial data, etc
  !-----------------------------------------------------------------

! From mpi_mod.f90
  call init
! From mpi_mod.f90
  call compute_partitioning
! From mpi_mod.f90
  nlx = numcolx()   !# of elements in block-x for a given rank
  nly = numcoly()   !# of elements in block-y for a given rank
  nlz = nez         ! no partion in the vertical 

call allocateData(nlx, nly, nlz)
if ( rank_is_master ) then
call allocateGlobalData(nex, ney, nez)
end if


!  allocate(grid(nlx,nly,nlz))
! !allocate(metric(nlx,nly,nlz)) ! defined in mountain_grid_mod
!  allocate(vars(nlx,nly,nlz))
!  allocate(vars_init(nlx,nly,nlz))
!  allocate(vars_phys(nlx,nly,nlz))

  ! set up the grid and initial data
! From testcases_nh3d_mod.f90
  call load_params()
! From grid_setup_mod.f90
  call make_grid()
! From testcases_nh3d_mod.f90
  call initialize_testcase()
!  call initialize_testcase(grid, vars_init)
!  vars = vars_init
call varsInitToVars()

  if (rank_is_master)  then
! From testcases_nh3d_mod.f90
    !print*, "domain size: ", xmax, ymax, zmax
    Call print_vital_stat( ) 
  endif 

  ! save initial t=0 movie frame (if enabled)
  frameid = 0
  if (save_movie_frames) then
    call gather_solution_onto_master_rank()
!    call gather_solution_onto_master_rank(vars, global_vars_phys)
    if (rank_is_master) call save_ycut_frame(frameid)
!    if (rank_is_master) call save_ycut_frame(global_vars_phys, frameid)
    frameid = frameid+1
  end if


  !-----------------------------------------------------------------
  ! main evolution loop
  !-----------------------------------------------------------------

  call cpu_time(timer_start)
!call system_clock(start_clock, rate_clock)

 !if (rank_is_master) print*, 'writing initial data'
 ! call gather_solution_onto_master_rank(vars_init, global_vars_phys)
 ! if (rank_is_master) call save_global_state(global_vars_phys)

  ! evolve
  t = 0.0D0
! nsteps = int(endtime/dt)
! nsteps = 16
! print*, 'N-steps? '
! read*, nsteps

!For IGW with dt=0.2s, end_time=3000s
  nsteps = 15000   
!  nsteps = 200
  nfreq = 300    
!  nfreq = 1

!$acc update device(vars_rho, vars_the, vars_prs, vars_spd, vars_psi, &
!$acc vars_u, vars_v, vars_w, vars_wt, vars_prp, vars_rhp, vars_thp, &
!$acc vars_vec, vars_flx, vars_fly, vars_flz, vars_src, vars_zsrc, &
!$acc vars_rhb, vars_thb, vars_rtb, vars_prb, vars_fcori, &
!$acc metric_mtn, metric_sg, metric_sg13, metric_sg23, metric_ght, &
!$acc grid_x, grid_y, grid_z, rk_stage_svec, rk_stage_svec0, &
!$acc edgevars_edges_rho_int, edgevars_edges_rho_ext, edgevars_edges_the_int, &
!$acc edgevars_edges_the_ext, edgevars_edges_prp_int, edgevars_edges_prp_ext, &
!$acc edgevars_edges_spd_int, edgevars_edges_spd_ext, edgevars_edges_u_int, &
!$acc edgevars_edges_u_ext, edgevars_edges_v_int, edgevars_edges_v_ext, &
!$acc edgevars_edges_w_int, edgevars_edges_w_ext, edgevars_edges_vec_int, &
!$acc edgevars_edges_vec_ext, edgevars_edges_normflux_int, edgevars_edges_normflux_ext, &
!$acc edgevars_edges_rtb, edgevars_edges_rhb, edgevars_edges_thb, &
!$acc edgevars_edges_prb, velem_rho, velem_the, velem_prs, velem_spd, &
!$acc velem_psi, velem_u, velem_v, velem_w, velem_prp, velem_rhp, &
!$acc velem_thp, velem_vec, velem_flx, velem_fly, velem_flz, velem_rhs,elem_rhs)
!
!$acc update device(name_of_the_test, endtime, xmin, xmax, ymin, ymax, zmin, zmax, delx, dely, delz, gllp, gllw, der, weakder)
!$acc update device(bufsendS, bufrecvS, bufsendE, bufrecvE, bufsendN, bufrecvN, bufsendW, bufrecvW)

  do itn = 1, nsteps

!!$acc update device(vars_rho, vars_the, vars_prs, vars_spd, vars_psi, &
!!$acc vars_u, vars_v, vars_w, vars_wt, vars_prp, vars_rhp, vars_thp, &
!!$acc vars_vec, vars_flx, vars_fly, vars_flz, vars_src, vars_zsrc, &
!!$acc vars_rhb, vars_thb, vars_rtb, vars_prb, vars_fcori, &
!!$acc metric_mtn, metric_sg, metric_sg13, metric_sg23, metric_ght, &
!!$acc grid_x, grid_y, grid_z, rk_stage_svec, rk_stage_svec0, &
!!$acc edgevars_edges_rho_int, edgevars_edges_rho_ext, edgevars_edges_the_int, &
!!$acc edgevars_edges_the_ext, edgevars_edges_prp_int, edgevars_edges_prp_ext, &
!!$acc edgevars_edges_spd_int, edgevars_edges_spd_ext, edgevars_edges_u_int, &
!!$acc edgevars_edges_u_ext, edgevars_edges_v_int, edgevars_edges_v_ext, &
!!$acc edgevars_edges_w_int, edgevars_edges_w_ext, edgevars_edges_vec_int, &
!!$acc edgevars_edges_vec_ext, edgevars_edges_normflux_int, edgevars_edges_normflux_ext, &
!!$acc edgevars_edges_rtb, edgevars_edges_rhb, edgevars_edges_thb, &
!!$acc edgevars_edges_prb, velem_rho, velem_the, velem_prs, velem_spd, &
!!$acc velem_psi, velem_u, velem_v, velem_w, velem_prp, velem_rhp, &
!!$acc velem_thp, velem_vec, velem_flx, velem_fly, velem_flz, velem_rhs,elem_rhs)

!!$acc update device(name_of_the_test, endtime, xmin, xmax, ymin, ymax, zmin, zmax, delx, dely, delz, gllp, gllw, der, weakder)
!!$acc update device(bufsendS, bufrecvS, bufsendE, bufrecvE, bufsendN, bufrecvN, bufsendW, bufrecvW)

    ! take time step
     call ssp_rk3(t, dt)

!!$acc update self(vars_rho, vars_the, vars_prs, vars_spd, vars_psi, &
!!$acc vars_u, vars_v, vars_w, vars_wt, vars_prp, vars_rhp, vars_thp, &
!!$acc vars_vec, vars_flx, vars_fly, vars_flz, vars_src, vars_zsrc, &
!!$acc vars_rhb, vars_thb, vars_rtb, vars_prb, vars_fcori, &
!!$acc metric_mtn, metric_sg, metric_sg13, metric_sg23, metric_ght, &
!!$acc grid_x, grid_y, grid_z, rk_stage_svec, rk_stage_svec0, &
!!$acc edgevars_edges_rho_int, edgevars_edges_rho_ext, edgevars_edges_the_int, &
!!$acc edgevars_edges_the_ext, edgevars_edges_prp_int, edgevars_edges_prp_ext, &
!!$acc edgevars_edges_spd_int, edgevars_edges_spd_ext, edgevars_edges_u_int, &
!!$acc edgevars_edges_u_ext, edgevars_edges_v_int, edgevars_edges_v_ext, &
!!$acc edgevars_edges_w_int, edgevars_edges_w_ext, edgevars_edges_vec_int, &
!!$acc edgevars_edges_vec_ext, edgevars_edges_normflux_int, edgevars_edges_normflux_ext, &
!!$acc edgevars_edges_rtb, edgevars_edges_rhb, edgevars_edges_thb, &
!!$acc edgevars_edges_prb, velem_rho, velem_the, velem_prs, velem_spd, &
!!$acc velem_psi, velem_u, velem_v, velem_w, velem_prp, velem_rhp, &
!!$acc velem_thp, velem_vec, velem_flx, velem_fly, velem_flz, velem_rhs,elem_rhs)
!
!!$acc update self(name_of_the_test, endtime, xmin, xmax, ymin, ymax, zmin, zmax, delx, dely, delz, gllp, gllw, der, weakder)
!!$acc update self(bufsendS, bufrecvS, bufsendE, bufrecvE, bufsendN, bufrecvN, bufsendW, bufrecvW)

    ! save movie frame (if enabled)
  ! if (save_movie_frames .and. mod(itn, movie_frame_freq)==0) then
  !   call convert_vars_to_physical_frame(vars, vars_phys)
  !   call gather_solution_onto_master_rank(vars, global_vars_phys)
  !   if (rank_is_master) call save_ycut_frame(global_vars_phys, frameid)
  !   frameid = frameid+1
  ! end if

    if (mod(itn, nfreq)==0) then
      call gather_solution_onto_master_rank()
!      call gather_solution_onto_master_rank(vars, global_vars_phys)
      if (rank_is_master) call print_diagnostics(itn)
!      if (rank_is_master) call print_diagnostics(itn,global_vars_phys)
    end if

    ! print status message
    if (rank_is_master .and. mod(itn,nfreq) == 0) then
      print('(a,i0,a,f0.2)'), "    step = ", itn, "    simulation time = ", t
    end if
  end do

!$acc update self(vars_rho, vars_the, vars_prs, vars_spd, vars_psi, &
!$acc vars_u, vars_v, vars_w, vars_wt, vars_prp, vars_rhp, vars_thp, &
!$acc vars_vec, vars_flx, vars_fly, vars_flz, vars_src, vars_zsrc, &
!$acc vars_rhb, vars_thb, vars_rtb, vars_prb, vars_fcori, &
!$acc metric_mtn, metric_sg, metric_sg13, metric_sg23, metric_ght, &
!$acc grid_x, grid_y, grid_z, rk_stage_svec, rk_stage_svec0, &
!$acc edgevars_edges_rho_int, edgevars_edges_rho_ext, edgevars_edges_the_int, &
!$acc edgevars_edges_the_ext, edgevars_edges_prp_int, edgevars_edges_prp_ext, &
!$acc edgevars_edges_spd_int, edgevars_edges_spd_ext, edgevars_edges_u_int, &
!$acc edgevars_edges_u_ext, edgevars_edges_v_int, edgevars_edges_v_ext, &
!$acc edgevars_edges_w_int, edgevars_edges_w_ext, edgevars_edges_vec_int, &
!$acc edgevars_edges_vec_ext, edgevars_edges_normflux_int, edgevars_edges_normflux_ext, &
!$acc edgevars_edges_rtb, edgevars_edges_rhb, edgevars_edges_thb, &
!$acc edgevars_edges_prb, velem_rho, velem_the, velem_prs, velem_spd, &
!$acc velem_psi, velem_u, velem_v, velem_w, velem_prp, velem_rhp, &
!$acc velem_thp, velem_vec, velem_flx, velem_fly, velem_flz, velem_rhs,elem_rhs)

!$acc update self(name_of_the_test, endtime, xmin, xmax, ymin, ymax, zmin, zmax, delx, dely, delz, gllp, gllw, der, weakder)
!$acc update self(bufsendS, bufrecvS, bufsendE, bufrecvE, bufsendN, bufrecvN, bufsendW, bufrecvW)

  call cpu_time(timer_end)
!call system_clock(stop_clock, rate_clock)
  if (rank_is_master) then
!    print('(a,f0.1,a)'), '==> CPU time (in evolution loop) = ', (stop_clock - start_clock)/REAL(rate_clock), " sec"

    print('(a,f0.1,a)'), '==> CPU time (in evolution loop) = ', timer_end - timer_start, " sec"
  end if


  !-----------------------------------------------------------------
  ! post-run processing and data handling
  !-----------------------------------------------------------------

  ! output for visualization
  call gather_solution_onto_master_rank()
!  call gather_solution_onto_master_rank(vars, global_vars_phys)
  if (rank_is_master) call save_global_state()
!  if (rank_is_master) call save_global_state(global_vars_phys)

  ! error calculations
  call Gather_EulerSoln_onto_MasterRank(1)
!  call Gather_EulerSoln_onto_MasterRank(vars_init, global_vars_init)
  call Gather_EulerSoln_onto_MasterRank(2)
!  call Gather_EulerSoln_onto_MasterRank(vars, global_vars)
!! if (rank_is_master) call error_norms(global_vars_init, global_vars)
  if (rank_is_master) call error_norms(2)
!  if (rank_is_master) call error_norms(global_vars_init,global_vars,2)

call deallocateData()
  call finalize

end program main

