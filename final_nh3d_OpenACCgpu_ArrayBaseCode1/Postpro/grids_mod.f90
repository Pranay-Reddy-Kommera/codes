!!-------------------------------------------------------------------
!! Discontinuos Galerkin (DG), NH Grid initialization
!!-------------------------------------------------------------------
 Module grids_mod
     Use Basic_mod
       Use gauss_quadrature_mod
       Use Interpol_mod, only : Interpol_Matrix 

       Integer, public, parameter :: nsp = nx , ntp = nx - 2
       Public  :: quadrature_initilize, make_global_grid , Check_3derrors 
       public  :: make_basic_1dgrid, Testcase_dat, cosine_blob
       public  :: tst_function, Gll2Gll_1d,  Gll3d_2_Gll3d
       private  :: Element_int, Volume_int

  ! Local storage holding the 1D x/y/z grid locations.
  ! The 3D grids are constructed by the tensor product of these 1D arrays.
  !   global_x1d -> every grid point along this dimension
  !   viz_x1d    -> same, without dupes at interfaces; for visualization grid

  public
   !source pts 
   type :: grid_t
     sequence
     real(kind=real_kind) :: xs(nsp,nsp,nsp)
     real(kind=real_kind) :: ys(nsp,nsp,nsp)
     real(kind=real_kind) :: zs(nsp,nsp,nsp)
   !target pts 
     real(kind=real_kind) :: xt(ntp,ntp,ntp)
     real(kind=real_kind) :: yt(ntp,ntp,ntp)
     real(kind=real_kind) :: zt(ntp,ntp,ntp)
   end type

   real(kind=real_kind), public :: mat_ip(nsp,ntp) 
   real(Kind=real_kind), public, Dimension(nsp) :: gll_s, gll_sw
   real(Kind=real_kind), public, Dimension(ntp) :: gll_t, gll_tw


      Contains
!!======================================================================================
!! Initiltialize both GLL and GL quadrture with derivative matrices 
       Subroutine quadrature_initialize
       Implicit None
       Real(Kind=real_kind), Dimension(nx,nlg) :: hmat

       Real(Kind=real_kind), Dimension(nsp,nsp) :: dmat_s
       Real(Kind=real_kind), Dimension(ntp,ntp) :: dmat_t
      !Real(Kind=real_kind), Dimension(nsp) :: gll_s, gll_sw
      !Real(Kind=real_kind), Dimension(ntp) :: gll_t, gll_tw 
      !type(grid_t),  dimension(nex,ney,nez)  :: grids         
      !Real(Kind=real_kind), dimension(nsp,nsp,nsp,nex,ney,nez)  :: src_fld       
      !Real(Kind=real_kind), dimension(ntp,ntp,ntp,nex,ney,nez)  :: tar_fld       
      !Real(Kind=real_kind), dimension(nsp,nsp,nsp)  :: tst_fld       
      !Real(Kind=real_kind), dimension(ntp,ntp,ntp)  :: tgi_fld, exact_fld

       Real(Kind=real_kind) :: val, fs(nsp), fi(ntp), fa(ntp) 

       Integer :: i,j,k,l,  kj,n, g_type 

     ! Gauss-Lobatto points (1d) and Gaussin weights 

         Call Gauss_Lobatto(nsp,gllp,gllw,der_gll)
         Call  Legendre_Poly(nsp,gllp,pmx,dpm)
             gll_s(:) = gllp(:)
            gll_sw(:) = gllw(:) 

          ! do k = 1,nx
          !   print*, k, gllp(k),gllw(k) 
          ! enddo 
         Call Gauss_Lobatto(ntp,gll_t,gll_tw,dmat_t)

       !    do k = 1,nsp 
       !      fs(k) = tst_function(gllp(k))
       !    enddo 
       !    do k = 1,ntp 
       !      fa(k) = tst_function(gll_t(k))
       !    enddo 

    ! interpolation matrix 
         mat_ip(:,:) =  Interpol_Matrix(nsp,ntp)

         print*, 'Done quadrature initialization' 
         print*, 'Soure/Target polynomial order ' , nsp, ntp 
         print*, 'Number of elements in (x,y,z) ' , nex, ney, nez 
         

    ! test with a simple function 
          !    fi(:) = Gll2Gll_1d(fs,nsp,ntp) 
          ! do k = 1,ntp
          !   print*, k, fa(k), fi(k)               
          ! enddo 


     !    g_type = 0 !source grid 
     !    Call Make_global_grid(grids,nsp,gll_s,g_type) 
     !    Call Testcase_dat(grids,nsp,g_type,src_fld) 

     !    g_type = 1 !target grid 
     !    Call Make_global_grid(grids,ntp,gll_t,g_type) 
     !    Call Testcase_dat(grids,ntp,g_type,tar_fld) 

     !      tst_fld(:,:,:) = src_fld(:,:,:,16,8,4) 

     !      exact_fld(:,:,:) = tar_fld(:,:,:,16,8,4) 

     !      tgi_fld(:,:,:) = Gll3d_2_Gll3d(tst_fld,nsp,ntp) 

           !do k = 1,nsp 
           !  print*, k, tst_fld(k,2,2)          
           !enddo 
     !     do k = 1,ntp 
     !       print*, k, tgi_fld(1,k,3), exact_fld(1,k,3)           
     !     enddo 

       End Subroutine quadrature_initialize

!!======================================================================================
!   A test function for 1d interpolation 
     Function tst_function(xi) result(val)
      real(Kind=real_kind), intent(in) :: xi 
      real(Kind=real_kind) :: val 
          val = 100.0D0 * xi**3 - 10.0D0 * xi**2 + 7.0D0 * sin(xi) + 50.0D0
     End Function tst_function

!!======================================================================================
!! 3D GLL gridbox to gridbox interpolation 
     Function Gll3d_2_Gll3d(src,spt,tpt) result(tar)
      integer, intent(in) :: spt,tpt
      real(Kind=real_kind), intent(in) :: src(spt,spt,spt)  
      real(Kind=real_kind) :: tar(tpt,tpt,tpt)  , s1 
      real(Kind=real_kind) :: tmp1(tpt,spt,spt), tmp2(tpt,tpt,spt) 
      integer :: i,j,k
        
       do k = 1, spt 
       do j = 1, spt 
        tmp1(:,j,k) = Gll2Gll_1d(src(:,j,k),spt,tpt)
       enddo
       enddo

       do k = 1, spt 
       do i = 1, tpt 
        tmp2(i,:,k) = Gll2Gll_1d(tmp1(i,:,k),spt,tpt)
       enddo
       enddo

       do j = 1, tpt 
       do i = 1, tpt 
        tar(i,j,:) = Gll2Gll_1d(tmp2(i,j,:),spt,tpt)
       enddo
       enddo
     End Function Gll3d_2_Gll3d
!!======================================================================================
     Function Gll2Gll_1d(fs,spt,tpt) result(ft)
      integer, intent(in) :: spt,tpt
      real(Kind=real_kind), intent(in) :: fs(spt) 
      real(Kind=real_kind) :: ft(tpt)  , s1 
      integer :: i,k
        
       do k = 1, tpt 
        s1 = 0.0D0
         do i = 1, spt 
          s1 = s1 + mat_ip(i,k) * fs(i)
         enddo
        ft(k) = s1
       enddo

     End Function Gll2Gll_1d   
!!======================================================================================
 Subroutine Check_3derrors(ngll,gll_wt,f3d_ana,f3d_num,errs) 
  implicit none 
  integer, intent(in) :: ngll
  real(kind=real_kind), intent(in), dimension(ngll) :: gll_wt 
  real(kind=real_kind), intent(in), dimension(ngll,ngll,ngll,nex,ney,nez) :: f3d_ana, f3d_num 
  real(kind=real_kind), intent(out), dimension(3) :: errs   
  real(kind=real_kind), dimension(ngll,ngll,ngll,nex,ney,nez) :: dif, dummy 

    real(kind=real_kind) :: ref_l1, ref_l2, ref_li
    real(kind=real_kind) :: err_l1, err_l2, err_li
    integer :: ie, je, ke

    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex
          dif(:,:,:,ie,je,ke) = f3d_num(:,:,:,ie,je,ke) - f3d_ana(:,:,:,ie,je,ke)
        end do
      end do
    end do

!! sanity check for volume integrals 
!   dummy = 1.0D0 
!   print*, "normalized integral of 1 = ", Volume_int(ngll,gll_wt,dummy)

    ! norms: L1, L2, Linf
 !  ref_l1 = volume_int(ngll,gll_wt,f3d_ana)
 !  ref_l2 = sqrt(volume_int(ngll,gll_wt,f3d_ana*f3d_ana))
 !  ref_li = maxval(abs(f3d_ana))

    err_l1 = Volume_int(ngll,gll_wt,dif)
    err_l2 = sqrt(Volume_int(ngll,gll_wt,dif*dif))
    err_li = maxval(abs(dif))

    errs(1) = err_l1 
    errs(2) = err_l2 
    errs(3) = err_li 


 End Subroutine Check_3derrors
!!======================================================================================
  ! NORMALIZED volume integral over domain
  Function Volume_int(ngll,gll_wt,f3d) result(vol_int) 
  integer, intent(in) :: ngll
  real(kind=real_kind), intent(in), dimension(ngll) :: gll_wt 
    real(kind=real_kind) :: vol_int
    real(kind=real_kind), intent(in), dimension(ngll,ngll,ngll,nex,ney,nez) :: f3d 
    real(kind=real_kind) :: sms 
    integer :: ie, je, ke

    sms = 0.0D0
    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex
          sms = sms + Element_int(ngll,gll_wt,f3d(:,:,:,ie,je,ke))
        end do
      end do
    end do

    vol_int = sms / ((xmax-xmin) * (ymax-ymin) * (zmax-zmin))

  end Function volume_int


!!======================================================================================
  ! volume integral over one element
  Function Element_int(ngll,gll_wt,elm) result(val_int)
  integer, intent(in) :: ngll
  real(kind=real_kind), intent(in), dimension(ngll) :: gll_wt 
  real(kind=real_kind), intent(in), dimension(ngll,ngll,ngll) :: elm
  real(kind=real_kind) :: val_int 
    real(kind=real_kind) :: sms 
    integer :: i, j, k

    sms = 0.0D0
    do k = 1, ngll
      do j = 1, ngll
        do i = 1, ngll
          sms = sms + elm(i,j,k) * gll_wt(i) * gll_wt(j) * gll_wt(k)
        end do
      end do
    end do

    val_int  = sms * (delx * dely * delz)/8.0D0 
  end Function Element_int

!!======================================================================================
 Subroutine Make_global_grid(grid,ngll,gll_pt,g_type)
  implicit none 
  integer, intent(in) :: ngll, g_type  
  type(grid_t), intent(out), dimension(nex,ney,nez) :: grid
  real(kind=real_kind), intent(in), dimension(ngll) :: gll_pt 

  real(kind=real_kind), dimension(ngll,nex) :: global_x1d
  real(kind=real_kind), dimension(ngll,ney) :: global_y1d
  real(kind=real_kind), dimension(ngll,nez) :: global_z1d
  real(kind=real_kind), dimension((ngll-1)*nex+1) :: viz_x1d
  real(kind=real_kind), dimension((ngll-1)*ney+1) :: viz_y1d
  real(kind=real_kind), dimension((ngll-1)*nez+1) :: viz_z1d
    integer :: i, j, k, ie, je, ke, npx,npy,npz               

      npx = (ngll-1)*nex + 1 
      npy = (ngll-1)*ney + 1 
      npz = (ngll-1)*nez + 1 

  ! Domain specification for IGW test 
      xmin = 0.0D0
      xmax = 320000.0D0
      ymin = 0.0D0
      ymax = 160000.0D0
      zmin = 0.0D0
      zmax = 10000.0D0

   Call make_basic_1dgrid(xmin, xmax, ngll,nex, npx, gll_pt, delx, global_x1d, viz_x1d)
   Call make_basic_1dgrid(ymin, ymax, ngll,ney, npy, gll_pt, dely, global_y1d, viz_y1d)
   Call make_basic_1dgrid(zmin, zmax, ngll,nez, npz, gll_pt, delz, global_z1d, viz_z1d)

    ! Take tensor product to obtain full grid
    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex
          do k = 1, ngll
            do j = 1, ngll
              do i = 1, ngll
               if (g_type == 0 ) then 
                grid(ie,je,ke)%xs(i,j,k)  = global_x1d(i,ie)
                grid(ie,je,ke)%ys(i,j,k)  = global_y1d(j,je)
                grid(ie,je,ke)%zs(i,j,k)  = global_z1d(k,ke)
               elseif (g_type == 1) then
                grid(ie,je,ke)%xt(i,j,k)  = global_x1d(i,ie)
                grid(ie,je,ke)%yt(i,j,k)  = global_y1d(j,je)
                grid(ie,je,ke)%zt(i,j,k)  = global_z1d(k,ke)
               endif 
              end do
            end do
          end do
        end do
      end do
    end do

 if (g_type == 0)    print*, "source grid generated" 
 if (g_type == 1)    print*, "target grid generated" 

  if (g_type == 0 ) then 
   3 format(e24.16) ! use 'e' rather than 'd' notation for python plotting
    open (unit = 25, file = 'grid_geometry.dat')
    write(25,*) ngll, nex, ney, nez
    close(25)
    open (unit = 35, file = 'grid_xs.dat')
    write(35,3) (viz_x1d(i),i=1,npx)
    close(35)
    open (unit = 35, file = 'grid_ys.dat')
    write(35,3) (viz_y1d(i),i=1,npy)
    close(35)
    open (unit = 35, file = 'grid_zs.dat')
    write(35,3) (viz_z1d(i),i=1,npz)
    close(35)
  endif 

  end subroutine Make_global_grid 

!!======================================================================================
  Subroutine make_basic_1dgrid(smin, smax, ngll,nels, nviz, gll_pt, dels, global_s1d, viz_s1d)
    Implicit none 
    real(kind=real_kind), intent(in) :: smin, smax
    integer, intent(in) :: ngll, nels, nviz
    real(kind=real_kind), intent(in), dimension(ngll) :: gll_pt
    real(kind=real_kind), intent(out) :: dels
    real(kind=real_kind), intent(out), dimension(ngll,nels) :: global_s1d
    real(kind=real_kind), intent(out), dimension(nviz) :: viz_s1d

    real(kind=real_kind), dimension(nels+1) :: sedges
    integer :: i, k, ki

    ! widths of elements
    dels = (smax-smin) / dble(nels) 

    ! boundaries of each 1d interval within [xmin, xmax]
    do k = 1, nels+1
      sedges(k) = smin + dels * (k-1)
    end do

    ! mapped values of GLL points
    do k = 1, nels
      do i = 1, ngll
        ! set the full points (with dupes) for evolution
        global_s1d(i,k) = (sedges(k) + sedges(k+1) + dels * gll_pt(i)) / 2.0D0
        ! set the reduced points (without dupes) for visualization
        if (i==ngll) exit ! skip the dupe
        ki = (k-1)*(ngll-1) + i
        viz_s1d(ki) = global_s1d(i,k)
      end do
    end do
    viz_s1d(nviz) = smax ! last point

  end subroutine make_basic_1dgrid

!!===============================================================================
 subroutine testcase_dat(grid,ngll,g_type,vars)
    implicit none 
    integer, intent(in) :: ngll,g_type 
    type(grid_t), intent(in), dimension(:,:,:) :: grid
    real(kind=real_kind), intent(out), dimension(ngll,ngll,ngll,nex,ney,nez) :: vars

    real(kind=real_kind), dimension(ngll,ngll,ngll) :: psi, xg,yg,zg 
    real(kind=real_kind), dimension((ngll-1)*nex+1,(ngll-1)*ney+1,(ngll-1)*nez+1) :: fglob

    integer :: ie, je, ke, iu,ju,ku 
    integer :: i, j, k, npx,npy,npz, midx,midy,midz                    

      npx = (ngll-1)*nex + 1
      npy = (ngll-1)*ney + 1
      npz = (ngll-1)*nez + 1

    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex

               if (g_type == 0 ) then 
                xg(:,:,:) = grid(ie,je,ke)%xs(:,:,:) 
                yg(:,:,:) = grid(ie,je,ke)%ys(:,:,:) 
                zg(:,:,:) = grid(ie,je,ke)%zs(:,:,:) 
               elseif (g_type == 1) then
                xg(:,:,:) = grid(ie,je,ke)%xt(:,:,:) 
                yg(:,:,:) = grid(ie,je,ke)%yt(:,:,:) 
                zg(:,:,:) = grid(ie,je,ke)%zt(:,:,:) 
               endif 
          vars(:,:,:,ie,je,ke) = Cosine_Blob(ngll,xg,yg,zg) 
        end do
      end do
    end do

    ! rewrite psi data as a simple 3D grid with no dupes
    do ke = 1, nez
      do k = 1, ngll-1
        ku = (ke-1)*(ngll-1) + k

        do je = 1, ney
          do j = 1, ngll-1
            ju = (je-1)*(ngll-1) + j

            do ie = 1, nex
              do i = 1, ngll-1
                iu = (ie-1)*(ngll-1) + i
                fglob(iu,ju,ku) = vars(i,j,k,ie,je,ke)
              end do
            end do

          end do
        end do

      end do
    end do

   ! doubly periodic boundary (extension)
    fglob(npx,:,:) = fglob(1,:,:)
    fglob(:,npy,:) = fglob(:,1,:)
    fglob(:,:,npz) = fglob(:,:,1)

    midx = int((npx+1)/2.0D0)
    midy = int((npy+1)/2.0D0)
    midz = int((npz+1)/2.0D0)

!! 2d slice and full 3D data without dupes 
   3 format(e24.16) 
    open (unit = 30, file = 'data_Xcut.dat')
    write(30,3) ((fglob(midx,i,j),i=1,npy),j=1,npz)
    close(30)
    open (unit = 30, file = 'data_Ycut.dat')
    write(30,3) ((fglob(i,midy,j),i=1,npx),j=1,npz)
    close(30)
    open (unit = 30, file = 'data_Zcut.dat')
    write(30,3) ((fglob(i,j,midz),i=1,npx),j=1,npy)
    close(30)
      open (unit = 30, file = 'data_full.dat')
      write(30,3) (((fglob(i,j,k),i=1,npx),j=1,npy),k=1,npz)
      close(30)


  end subroutine testcase_dat
!----------------------------------------------------------------------

 Function cosine_blob(ngll,xg,yg,zg)

    integer, intent(in) :: ngll 
    real(kind=real_kind), dimension(ngll,ngll,ngll), intent(in) :: xg,yg,zg
    real(kind=real_kind):: x0, y0, z0, a0, b0, c0

    real(kind=real_kind), dimension(ngll,ngll,ngll) :: cosine_blob
    real(kind=real_kind) :: dist

    integer :: i,j,k

     x0 = xmax/2.0D0
     y0 = ymax/2.0D0
     z0 = zmax/2.0D0
     a0 = 120000.0D0
     b0 = 60000.0D0
     c0 = 4000.0D0

    do k = 1, ngll
      do j = 1, ngll
        do i = 1, ngll

          dist = sqrt( ((xg(i,j,k) - x0)/a0)**2 &
                     + ((yg(i,j,k) - y0)/b0)**2 &
                     + ((zg(i,j,k) - z0)/c0)**2 )
          if (dist <= 1.0D0 ) then
            cosine_blob(i,j,k) = 1.0D0 * (cos(pi*0.5D0 * dist))**2
          else
            cosine_blob(i,j,k) = 0.0D0
          endif

        end do
      end do
    end do

  end Function cosine_blob


!!===============================================================================



       End Module grids_mod 

