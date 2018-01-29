!-------------------------------------------------------------------
!   Numerical method: Disc. Galerkin, nodal, weak, GLL grid
!-------------------------------------------------------------------

program post_process
  use basic_mod
  use gauss_quadrature_mod
  use grids_mod    
  use interpol_mod 

  implicit none

  real(kind=real_kind), dimension(nx,nx,nx,nex,ney,nez)  :: data_read, tst_dat
  real(kind=real_kind) :: t, timer_start, timer_end, err3(3) 
  integer :: nlx, nly, nlz, g_type 

       type(grid_t),  dimension(nex,ney,nez)  :: grids
       Real(Kind=real_kind), dimension(nsp,nsp,nsp,nex,ney,nez)  :: src_fld, red_fld
       Real(Kind=real_kind), dimension(ntp,ntp,ntp,nex,ney,nez)  :: ref_fld, tgt_fld 
       Real(Kind=real_kind), dimension(nsp,nsp,nsp)  :: tst_fld
       Real(Kind=real_kind), dimension(ntp,ntp,ntp)  :: tgi_fld, exact_fld

  integer :: i,j,k, itn, ie,je,ke         

  !-----------------------------------------------------------------
  ! setup work:
  ! init MPI, build grid and initial data, save initial data, etc
  !-----------------------------------------------------------------

!  nlx = numcolx()   !# of elements in block-x for a given rank
!  nly = numcoly()   !# of elements in block-y for a given rank
!  nlz = nez         ! no partion in the vertical 


    Call quadrature_initialize( ) 


          g_type = 0 !source grid 
          Call Make_global_grid(grids,nsp,gll_s,g_type)
          Call Testcase_dat(grids,nsp,g_type,src_fld)

          g_type = 1 !target grid 
          Call Make_global_grid(grids,ntp,gll_t,g_type)
          Call Testcase_dat(grids,ntp,g_type,ref_fld)

          ! tst_fld(:,:,:) = src_fld(:,:,:,16,8,4)

          ! exact_fld(:,:,:) = ref_fld(:,:,:,16,8,4)

          ! tgi_fld(:,:,:) = Gll3d_2_Gll3d(tst_fld,nsp,ntp)


         open (unit=31,file='total3d.dat',form='unformatted') 
         write(31) ((( (((src_fld(i,j,k,ie,je,ke),i=1,nsp),j=1,nsp),k=1,nsp),ie=1,nex),je=1,ney),ke=1,nez)
         close(31) 
         
         open (unit=31,file='total3d.dat',form='unformatted') 
         read(31) ((( (((red_fld(i,j,k,ie,je,ke),i=1,nsp),j=1,nsp),k=1,nsp),ie=1,nex),je=1,ney),ke=1,nez)
         close(31) 
         

        do ke = 1, nez
        do je = 1, ney
        do ie = 1, nex
           
            tst_fld(:,:,:) = red_fld(:,:,:,ie,je,ke)
            tgi_fld(:,:,:) = Gll3d_2_Gll3d(tst_fld,nsp,ntp)
            tgt_fld(:,:,:,ie,je,ke) = tgi_fld(:,:,:) 

        enddo 
        enddo 
        enddo 

             Call Check_3derrors(ntp,gll_tw,ref_fld,tgt_fld,err3) 

        print*, " L-errors : ", err3(1), err3(2), err3(3) 

           !do k = 1,nsp 
           !  print*, k, tst_fld(k,2,2)          
           !enddo 
           !do k = 1,ntp
           !  print*, k, tgi_fld(1,k,3), exact_fld(1,k,3)
           !enddo


  !-----------------------------------------------------------------
  ! post-run processing and data handling
  !-----------------------------------------------------------------

! output for visualization

! error calculations
! call gather_solution_onto_master_rank(vars, global_vars)
! call gather_solution_onto_master_rank(vars_init, global_vars_init)
! if (rank_is_master) call error_norms(global_vars_init, global_vars)


end program post_process

