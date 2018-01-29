!! R. Nair IMAGe, 01/2009
!! ------- Inspired from SFV work ---------
!  Interpolation using native high-order basis function on 
!  GLL grid to any given target grid (including low-order GLL)
!
   Module   Interpol_mod   
      Use Basic_mod
      Use gauss_quadrature_mod, only : gllp,pmx,hm_lg,ime_gl, ime_der_gl 


       Integer, Parameter :: nv = 3 , nel = nex, nvt = nel*(nv-1), neq = 1

       Real(Kind=real_kind), Dimension(nx,nv), Public :: hi       

       Public :: intialize_lg_intp 
       Public :: Interpol_Matrix
       Public :: Gll_2_Gll, GLL_grids    
       Public :: gl_pts_wts              
       Public :: edge_intparrays_int, edge_intparrays_ext, edge_internal_lg
       Public :: gll_2_lg1d, gll_2_lg2d 
 
       Public :: edge_simple_copy1, edge_simple_copy2, edge_simple_copy0
       Public :: edge_copy_xorz
       Public :: edge_glx_only, edge_glz_only

       Public :: edge_gl_elem, edge_gl_1d, edge_Dergl_1d
       Public :: interpol_matrix_gl2edge, Interpol_DerMatrix_gl2edge

       Private ::  Legend_Poly,  Legend_Deri, Interpol_tst
       !Private ::  Postprocess

  Contains
      
!!===============++++++++++++++++++++++++++++++=====================
    Function edge_glx_only(ff) result(f2gl)

      Real(Kind=real_kind), intent(in), Dimension(nx,nx) :: ff
      Real(Kind=real_kind), Dimension(nlg,2) :: f2gl
      Real(Kind=real_kind), Dimension(nlg) ::  fgl
      Real(Kind=real_kind) :: lr(2)

      Integer :: k,i,j

     ! East and West    (x-side)
      do k = 1, nx
        fgl(:) = ff(:,k)
        lr(:)  = edge_gl_1d(fgl)

        f2gl(k,1) = lr(1)
        f2gl(k,2) = lr(2)
      enddo
    End Function edge_glx_only
!!-----------------------------------------------------
    Function edge_glz_only(ff) result(f2gl)

      Real(Kind=real_kind), intent(in), Dimension(nx,nx) :: ff
      Real(Kind=real_kind), Dimension(nlg,2) :: f2gl
      Real(Kind=real_kind), Dimension(nlg) ::  fgl
      Real(Kind=real_kind) :: lr(2)
      Integer :: k,i,j

          ! South & North    (z-side)
      do k = 1, nx
        fgl(:) = ff(k,:)
        lr(:)  = edge_gl_1d(fgl)

        f2gl(k,1) = lr(1)
        f2gl(k,2) = lr(2)
      enddo

    End Function edge_glz_only
!!-----------------------------------------------------
    Function edge_gl_elem(ff) result(f4gl)

      Real(Kind=real_kind), intent(in), Dimension(nx,nx) :: ff
      Real(Kind=real_kind), Dimension(nx,4) :: f4gl
      Real(Kind=real_kind), Dimension(nx) ::  fgl
      Real(Kind=real_kind) :: lr(2)

      Integer :: k,i,j

     !! Source grid GLL target grid GL of same size 

     ! East and West 
      do k = 1, nx
        fgl(:) = ff(:,k)
        lr(:)  = edge_gl_1d(fgl)

        f4gl(k,4) = lr(1)
        f4gl(k,2) = lr(2)
      enddo

          ! South & North 
      do k = 1, nx
        fgl(:) = ff(k,:)
        lr(:)  = edge_gl_1d(fgl)

        f4gl(k,1) = lr(1)
        f4gl(k,3) = lr(2)
      enddo

    End Function edge_gl_elem
!!------------------------------------------------------
    Function edge_gl_1d(f1) result(gle)

      Real(Kind=real_kind), intent(in), Dimension(nx) :: f1
      Real(Kind=real_kind) ::  v1, v2, gle(2)
      Integer :: k


        v1 = 0.0
        v2 = 0.0
        do k=1, nx
         v1  = v1 + ime_gl(k,1)*f1(k)
         v2  = v2 + ime_gl(k,2)*f1(k)
        enddo

       !! [gl(-1), gl(+1)] 

        gle(1) = v1   !left edge 
        gle(2) = v2   !right edge 

    End Function edge_gl_1d
!!------------------------------------------------------
    Function edge_Dergl_1d(f1) result(d_gle)

      Real(Kind=real_kind), intent(in), Dimension(nx) :: f1
      Real(Kind=real_kind) ::  v1, v2, d_gle(2)
      Integer :: k


        v1 = 0.0
        v2 = 0.0
        do k=1, nx
         v1  = v1 + ime_der_gl(k,1)*f1(k)
         v2  = v2 + ime_der_gl(k,2)*f1(k)
        enddo

       !! [d_gl(-1), d_gl(+1)] 

        d_gle(1) = v1   !left edge 
        d_gle(2) = v2   !right edge 

    End Function edge_Dergl_1d

!-------------------------------------------
      Function    Interpol_Matrix_gl2edge(glg,ngp) result(hm)

      Implicit None
      Integer, Intent(in) :: ngp
      Real(Kind=real_kind),  Intent(in), Dimension(ngp) :: glg
      Real(Kind=real_kind),  Dimension(ngp,2) :: hm
      Real(Kind=real_kind) :: tgrid(2), pt, deno, term, sum

        Integer :: i,j,k,l, trivial,nd

             tgrid(1) = -1.0D0   ! Edge grid points 
             tgrid(2) =  1.0D0

      !Interpolation using basis function ("ngp" order)
      do k = 1, 2
            pt = tgrid(k)  !target grid 
           term =   Legend_Poly(ngp,pt)
            do j = 1, ngp
             deno = (pt - glg(j)) * Legend_Deri(ngp,glg(j))
             hm(j,k) = term / deno
            enddo
      enddo

      End Function    Interpol_Matrix_gl2edge
!-------------------------------------------
      Function    Interpol_DerMatrix_gl2edge(glg,ngp) result(der_hm)

      Implicit None
      Integer, Intent(in) :: ngp
      Real(Kind=real_kind),  Intent(in), Dimension(ngp) :: glg
      Real(Kind=real_kind),  Dimension(ngp,2) :: der_hm
      Real(Kind=real_kind) :: tgrid(2), pt, deno, term, sum

        Integer :: i,j,k,l, trivial,nd

             tgrid(1) = -1.0D0   ! Edge grid points 
             tgrid(2) =  1.0D0

     !Interpolation using derivative of basis function ("ngp" order)
      do k = 1, 2
            pt = tgrid(k)  !target grid 
           term =   Legend_Deri(ngp,pt)
            do j = 1, ngp
             deno = (pt - glg(j)) * Legend_Deri(ngp,glg(j))
             der_hm(j,k) = term / deno
            enddo
      enddo

      End Function    Interpol_DerMatrix_gl2edge

!!=========================================
      Function  edge_internal_lg(si) result(si4)

      Implicit None
      Real (Kind=real_kind), Intent(in), Dimension(nx,nx) :: si
      Real (Kind=real_kind), Dimension(nlg,4) :: si4
      Real (Kind=real_kind) :: evs(nx)
      Integer ::   i,j,k

           si4 = 0.0D0

         !Copy to west (internal) 
           evs(:) = si(1,:)
           si4(:,4) = gll_2_lg1d(evs)

         !East (internal) 
           evs(:) = si(nx,:)
           si4(:,2) = gll_2_lg1d(evs)

         !South (internal) 
           evs(:) = si(:,1)
           si4(:,1) = gll_2_lg1d(evs)

         !North(internal) 
           evs(:) = si(:,nx)
           si4(:,3) = gll_2_lg1d(evs)

      End Function edge_internal_lg
!!-----------------------------------------
      Function  edge_simple_copy0(si) result(si4)

      Implicit None
      Real (Kind=real_kind), Intent(in), Dimension(nx,nx) :: si
      Real (Kind=real_kind), Dimension(nlg,4) :: si4
      Real (Kind=real_kind) :: evs(nx)
      Integer ::   i,j,k

           si4 = 0.0D0

         !Copy to west (internal) 
           evs(:) = si(1,:)
           si4(:,4) = evs(:)         

         !East (internal) 
           evs(:) = si(nx,:)
           si4(:,2) = evs(:)         

         !South (internal) 
           evs(:) = si(:,1)
           si4(:,1) = evs(:)         

         !North(internal) 
           evs(:) = si(:,nx)
           si4(:,3) = evs(:)         

      End Function edge_simple_copy0
!!-----------------------------------------

     Function  edge_simple_copy2(si) result(si4)

      Implicit None
      Real (Kind=real_kind), Intent(in), Dimension(nx,nx,neq) :: si
      Real (Kind=real_kind), Dimension(nlg,4,neq) :: si4
      Integer ::   i,j,k

           si4 = 0.0D0

         do k = 1, neq 
         !Copy to west (internal) 
           si4(:,4,k) = si(1,:,k)

         !East (internal) 
           si4(:,2,k) = si(nx,:,k)  

         !South (internal) 
           si4(:,1,k) = si(:,1,k)         

         !North(internal) 
           si4(:,3,k) = si(:,nx,k)          
          enddo 

      End Function edge_simple_copy2
!----------------------------------------------
    Function  edge_copy_xorz(mark,si) result(si4)

      Implicit None
      Integer :: mark
      Real (Kind=real_kind), Intent(in), Dimension(nx,nx) :: si
      Real (Kind=real_kind), Dimension(nlg,2) :: si4
      Integer ::   i,j,k

           si4 = 0.0D0

        if  (mark == 1) then       !x-direction
         !Copy to west (internal) 
           si4(:,1) = si(1,:)

         !East (internal) 
           si4(:,2) = si(nx,:)
        elseif (mark == 2) then    !z-direction
         !South (internal) 
           si4(:,1) = si(:,1)

         !North(internal) 
           si4(:,2) = si(:,nx)
        endif

      End Function edge_copy_xorz    
!----------------------------------------------
    Function  edge_simple_copy1(mark,si) result(si4)

      Implicit None
      Integer :: mark
      Real (Kind=real_kind), Intent(in), Dimension(nx,nx,neq) :: si
      Real (Kind=real_kind), Dimension(nlg,2,neq) :: si4
      Integer ::   i,j,k

           si4 = 0.0D0

        if  (mark == 1) then 
         do k = 1, neq
         !Copy to west (internal) 
           si4(:,1,k) = si(1,:,k)

         !East (internal) 
           si4(:,2,k) = si(nx,:,k)
         enddo 
        elseif (mark == 2) then  
         do k = 1, neq

         !South (internal) 
           si4(:,1,k) = si(:,1,k)

         !North(internal) 
           si4(:,2,k) = si(:,nx,k)
          enddo
        endif 

      End Function edge_simple_copy1
!----------------------------------------------
      Subroutine  edge_intparrays_int(si,si4)

      Implicit None
      Real (Kind=real_kind), Intent(in), Dimension(nx,nx) :: si
      Real (Kind=real_kind), Intent(out), Dimension(4,nlg) :: si4
      Real (Kind=real_kind) :: evt(nlg), evs(nx)
      Integer ::   i,j,k

           si4 = 0.0D0

         !west (internal) 
           evs(:) = si(1,:)
           si4(4,:) = gll_2_lg1d(evs)

         !East (internal) 
           evs(:) = si(nx,:)
           si4(2,:) = gll_2_lg1d(evs)

         !South (internal) 
           evs(:) = si(:,1)
           si4(1,:) = gll_2_lg1d(evs)

         !North(internal) 
           evs(:) = si(:,nx)
           si4(3,:) = gll_2_lg1d(evs)

      End Subroutine edge_intparrays_int
!!-----------------------------------------
      Subroutine  edge_intparrays_ext(si4,si_edge)

      Implicit None
      Real (Kind=real_kind), Intent(in), Dimension(4,nx) :: si4
      Real (Kind=real_kind), Intent(out), Dimension(4,nlg) :: si_edge
      Real (Kind=real_kind) :: evt(nlg), evs(nx)
      Integer ::   i,j,k


         !South(external) 
           evs(:) = si4(1,:)
           si_edge(2,:) = gll_2_lg1d(evs)

         !East (external) 
           evs(:) = si4(2,:)
           si_edge(2,:) = gll_2_lg1d(evs)

         !North(external) 
           evs(:) = si4(3,:)
           si_edge(3,:) = gll_2_lg1d(evs)

         !West (external) 
           evs(:) = si4(4,:)
           si_edge(4,:) = gll_2_lg1d(evs)


      End Subroutine edge_intparrays_ext
!!-----------------------------------------

 Function  gll_2_lg1d(fs) result(ft)

      Implicit None
      Real(Kind=real_kind), Intent(in), Dimension(nx) :: fs  !Source data 
      Real(Kind=real_kind), Dimension(nlg) :: ft              !Target data 
      Real(Kind=real_kind) ::  sums

        Integer :: i,k

     ! Interpolating from GLL(nx) to GL(nlg)

       do k = 1, nlg
        sums = 0.0D0
         do i = 1, nx
          sums = sums + hm_lg(i,k) * fs(i)
         enddo
        ft(k) = sums
       enddo

      End Function  gll_2_lg1d

!!===============++++++++++++++++++++++++++++++=====================

    Function  gll_2_lg2d(fs) result(ft)

      Implicit None
      Real(Kind=real_kind), Intent(in), Dimension(nx,nx) :: fs  !Source data 
      Real(Kind=real_kind), Dimension(nlg,nlg) :: ft, dm_gl     !Target data 
      Real(Kind=real_kind), Dimension(nx) :: dmy           

        Integer :: k,l 
         !note: nx = nlg 

       do k = 1, nx 
         dmy(:) = fs(:,k)
         dm_gl(:,k) = gll_2_lg1d(dmy)
       enddo 

       do l = 1, nlg 
         dmy(:) = dm_gl(l,:)
         ft(l,:) = gll_2_lg1d(dmy)
       enddo 

      End Function  gll_2_lg2d
!!===============++++++++++++++++++++++++++++++=====================
      Subroutine initialize_lg_intp(hmat)
      Implicit None

      Real(Kind=real_kind), intent(out),   Dimension(nx,nlg) :: hmat
      Real(Kind=real_kind),  Dimension(nlg) ::  tgl, lg_pt,lg_wt 
      Real(Kind=real_kind) :: pt, epsil, deno, fact, term, sums

        Integer :: i,j,k,l, trivial,nd

     ! snx --> source Gll dimension
     ! tnv --> traget Gll dimension
     ! snx > tnv (Generally)

        nd = nx-1  !Degree of the  source Legend polynomial 

     ! lgp(:) = Gll_grids(nlg)   ! Target LG  grid  copy       

     ! Construction of interpolation matrix "hm(snx,tnv)"
     ! for 1D scan  [ U(x) = sums_k{ U(k) * hm(k,x) }]

          fact =  dble(nd*(nd+1))
       epsil = 1.0D-12

      hmat(:,:) = 0.0D0

      Call     gl_pts_wts(nlg,lg_pt,lg_wt)

      print*, 'Interploation matrix, LG pts & wts Initialized'

        tgl(:) = lg_pt(:) 
          
      do k = 1, nlg 
            pt = tgl(k)
            trivial = 0
          do i = 1, nx
            if (abs(pt-gllp(i)) < epsil ) then
            hmat(i,k) = 1.0D0
            trivial = 1     !avoiding trivial points 
            exit
            endif
          enddo

          if (trivial == 0 ) then
           term = (pt*pt - 1.0D0 ) *   Legend_Deri(nd,pt)
            do i = 1, nx
             deno = (pt - gllp(i)) * fact * pmx(nx,i)
             hmat(i,k) = term / deno
            enddo
          endif
      enddo

     ! do i= 1, nx 
     ! do k = 1, nlg 
     ! !print*, k, lgp(k), lgw(k) 
     !  print*, i,k, hmat(i,k)       
     !  enddo 
     !  enddo 

      End Subroutine   initialize_lg_intp 

!-------------------------------------------------------------------
       Subroutine interpol_tst(ff)

       Implicit None

       Real(Kind=real_kind), intent(in), Dimension(nx,nx,nel,nel) :: ff
       Real(Kind=real_kind), Dimension(nv,nv,nel,nel) :: flo_gll  
       Real(Kind=real_kind), Dimension(nx,nx) :: fn
       Real(Kind=real_kind), Dimension(nv,nv) :: fg 

       Real(Kind=real_kind) :: pt, epsil, deno, fact, term, sums 

        Integer :: i,j,k, ie,je

!-------------------------------------------
!        Interploation test
!-------------------------------------------

! Generate data-independent Interpolation matrix with Dim(nx,nv) 

      hi(:,:) = Interpol_Matrix(nx,nv)

! Interpolating from GLL(nx,nx) to GLL(nv,nv)

      do je = 1, nel 
        do ie = 1, nel 
           fn(:,:) = ff(:,:,ie,je)

           flo_gll(:,:,ie,je) = Gll_2_Gll(fn,nx,nv)

        end do
      end do

!  Call  Postprocess(flo_gll)
          

       End  Subroutine Interpol_tst

!-------------------------------------------
!      Subroutine Postprocess(fin)

!      Implicit None

!      Real(Kind=real_kind), intent(in), Dimension(nv,nv,nel,nel) :: fin
!      Real(Kind=real_kind), Dimension(nvt+1,nvt+1) :: fout
!      Real(Kind=real_kind), Dimension(nvt+1) :: ax_reduce
!      Real(Kind=real_kind), Dimension(nv) :: gllg

!       Integer ::   i,j,k, ki, kj, k1,k2

!       !Elemnetal gridpoint values

!       do k2 = 1, nel
!       do j = 1, nv-1
!             kj = (k2-1)*(nv-1) + j

!        do k1 = 1, nel
!        do i = 1, nv-1
!             ki = (k1-1)*(nv-1) + i
!               fout(ki,kj) = fin(i,j,k1,k2)
!        enddo
!        enddo

!       enddo
!       enddo

 !Doubly  Periodic Boundary

!        do k2 = 1, nvt+1
!          fout(nvt+1,k2) = fout(1,k2)
!        enddo
!        do k1 = 1, nvt+1
!          fout(k1,nvt+1) = fout(k1,1)
!        enddo

 !New reduced physical grid "ax_reduce"

!     gllg(:) = Gll_grids(nv)

 ! Global grid of independent points (see grid_maker.f90 )

!         do k = 1, nel
!          do j = 1, nv-1
!             kj = (k-1)*(nv-1) + j
!               ax_reduce(kj) = (xgl(k) + xgl(k+1) +  del *gllg(j) )/2.0D0
!          enddo
!         enddo

!         ax_reduce(nvt+1) = xgl(nel+1)        !last point +1 (periodic)

!       ! do k = 1, nvt+1
!       !     print*, k, ax_reduce(k)
!       ! enddo

!       open (unit = 40, file = 'reduced_gll.dat')
!       open (unit = 42, file = 'reduced_geo.dat')
!       open (unit = 45, file = 'reduced_ax.dat')
!       write(42,*) nv-1, nel
!       write(40,3) ((fout(i,j),i=1,nvt+1),j=1,nvt+1)
!       write(45,3) (ax_reduce(i),i=1,nvt+1)
!        3 format(10f10.5)
!        close(30)
!        close(35)

!      End  Subroutine Postprocess
!-------------------------------------------
      Function  Gll_2_Gll(fs,snx,tnv) result(ft)

      Implicit None
      Integer, Intent(in) :: snx, tnv 
      Real(Kind=real_kind), Intent(in), Dimension(snx,snx) :: fs  !Source data 
      Real(Kind=real_kind), Dimension(tnv,tnv) :: ft              !Target data 
      Real(Kind=real_kind), Dimension(tnv,snx) :: fm              !Temporary   
      Real(Kind=real_kind), Dimension(tnv) ::  tgl 
      Real(Kind=real_kind) ::  sums 

        Integer :: i,j,k,l, trivial,nd


     ! Interpolating from GLL(snx,snx) to GLL(tnv,tnv)

      do j = 1, snx 
       do k = 1, tnv 
        sums = 0.0D0 
         do i = 1, snx 
          sums = sums + hi(i,k) * fs(i,j) 
         enddo 
        fm(k,j) = sums 
       enddo 
      enddo 
          
      do l = 1, tnv 
       do k = 1, tnv 
        sums = 0.0D0 
         do j = 1, snx 
          sums = sums + hi(j,k) * fm(l,j) 
         enddo 
        ft(l,k) = sums 
       enddo 
      enddo 
          
      End Function  Gll_2_Gll
!-------------------------------------------
      Function    Interpol_Matrix(snx,tnv) result(hm)

      Implicit None
      Integer, Intent(in) :: snx, tnv 
      Real(Kind=real_kind),  Dimension(snx,tnv) :: hm 
      Real(Kind=real_kind),  Dimension(tnv) ::  tgl 
      Real(Kind=real_kind) :: pt, epsil, deno, fact, term, sums 

        Integer :: i,j,k,l, trivial,nd

     ! snx --> source Gll dimension
     ! tnv --> traget Gll dimension
     ! snx > tnv (Generally)

        nd = snx-1  !Degree of the  source Legend polynomial 

      tgl(:) = Gll_grids(tnv)   ! Target Gll grid  copy       

     ! Construction of interpolation matrix "hm(snx,tnv)"
     ! for 1D scan  [ U(x) = Sum_k{ U(k) * hm(k,x) }]

          fact =  dble(nd*(nd+1))
       epsil = 1.0D-12

      hm(:,:) = 0.0D0 

      do k = 1, tnv
            pt = tgl(k) 
            trivial = 0 
          do i = 1, snx
            if (abs(pt-gllp(i)) < epsil ) then
            hm(i,k) = 1.0D0
            trivial = 1     !avoiding trivial points 
            exit   
            endif 
          enddo

          if (trivial == 0 ) then
           term = (pt*pt - 1.0D0 ) *   Legend_Deri(nd,pt)  
            do i = 1, snx
             deno = (pt - gllp(i)) * fact * pmx(snx,i)
             hm(i,k) = term / deno 
            enddo 
          endif 
      enddo 

      End Function    Interpol_Matrix
!-------------------------------------------
      Function    Legend_Poly(deg,x) result(pval)

      Implicit None
      Integer, Intent(in) :: deg 
      Real(Kind=real_kind), Intent(in)   :: x   
      Real(Kind=real_kind)   :: xx,  pval  

          xx = x*x 

         SelectCase(deg)
          Case(1) 
             pval= x            
          Case(2) 
             pval= 1.5D0 * xx - 0.5D0 
          Case(3) 
             pval= 0.5D0 * x *(5.0D0 * xx - 3.0D0)
          Case(4) 
             pval= 0.125D0 *( xx*(35.0D0 * xx - 30.0D0) + 3.0D0)
          Case(5) 
             pval= 0.125D0 * x*( xx*(63.0D0 * xx - 70.0D0) + 15.0D0)
          Case(6) 
             pval= 0.0625D0 * ( xx*( xx*(231.0D0 * xx - 315.0D0) + 105.0D0) - 5.0D0)
          Case(7) 
             pval= 0.0625D0 * x *( xx*( xx*(429.0D0 * xx - 693.0D0) + 315.0D0) - 35.0D0)
          Case(8) 
             pval= 0.0078125D0 * ( xx*( xx*(xx*(6435.0D0 * xx - 12012.0D0) + 6930.0D0) - 1260.0D0) + &
                                 35.0D0)
          Case(9) 
             pval= 0.0078125D0 * x *( xx*( xx*(xx*(12155.0D0 * xx - 25740.0D0) + 18018.0D0) - &
                                    4620.0D0) + 315.0D0)
         End Select
      End Function    Legend_Poly
!-------------------------------------------
      Function    Legend_Deri(deg,x) result(pval)

      Implicit None
      Integer, Intent(in) :: deg 
      Real(Kind=real_kind), Intent(in)   :: x   
      Real(Kind=real_kind)   :: xx,  pval  

          xx = x*x 

         SelectCase(deg)
          Case(1)           
             pval = 1.0D0 
          Case(2)           
             pval = 3.0D0 * x
          Case(3)
             pval = 1.5D0 *(5.0D0 * xx - 1.0D0)
          Case(4)  
             pval = 2.5D0 * x *(7.0D0 * xx - 3.0D0)
          Case(5)
             pval = 1.875D0 * (7.0D0 * xx* (3.0D0 * xx - 2.0D0) + 1.0D0)
          Case(6) 
             pval = 0.375D0 * x * (3.0D0 * xx* (77.0D0 * xx - 70.0D0) + 35.0D0)
          Case(7) 
             pval = 0.0625 * (xx*(xx*(3003.0D0 * xx - 3465.0D0) + 945.0D0) - 35.0D0)
          Case(8) 
             pval = 0.56250D0 *  x*( xx*(xx*(715.0D0 * xx - 1001.0D0)+ 385.0D0)- 35.0D0)
          Case(9) 
             pval = 0.3515625D0 * ( xx*( xx*(xx*(2431.0D0 * xx - 4004.0D0)+ 2002.0D0)- 308.0D0)+ 7.0D0)
         End Select
      End Function    Legend_Deri
!-------------------------------------------
      Function  Gll_grids(np) result(gllg)

      Implicit None
      Integer, Intent(in) :: np 
      Real(Kind=real_kind),  Dimension(np) :: gllg 
      Real(Kind=real_kind) :: gll3(3),gll4(4),gll5(5),gll6(6),gll7(7),gll8(8),gll9(9) , &
                           gll10(10),gll11(11) 

        if (np > 11 )  then
            print*, 'No GLL grid availabe' 
        endif 

          gll3(1) = -1.00000000000000D0     
          gll3(2) = 0.000000000000000D0     
          gll3(3) =  1.00000000000000D0 


           gll4(1) = -1.00000000000000D0     
           gll4(2) = -0.447213595499958D0     
           gll4(3) =  0.447213595499958D0     
           gll4(4) =  1.00000000000000D0    

           gll5(1) =   -1.00000000000000D0     
           gll5(2) =  -0.654653670707977D0     
           gll5(3) =   0.000000000000000D0
           gll5(4) =   0.654653670707977D0    
           gll5(5) =    1.00000000000000D0    


           gll6(1) =   -1.00000000000000D0     
           gll6(2) =  -0.765055323929465D0     
           gll6(3) =  -0.285231516480645D0     
           gll6(4) =   0.285231516480645D0     
           gll6(5) =   0.765055323929465D0     
           gll6(6) =    1.00000000000000D0   
      
 
           gll7(1) =   -1.00000000000000D0     
           gll7(2) =  -0.830223896278567D0     
           gll7(3) =  -0.468848793470714D0     
           gll7(4) =   0.000000000000000D0
           gll7(5) =   0.468848793470714D0     
           gll7(6) =   0.830223896278567D0     
           gll7(7) =   1.000000000000000D0 



           gll8(1) =  -1.000000000000000D0 
           gll8(2) =  -0.871740148509607D0     
           gll8(3) =  -0.591700181433142D0     
           gll8(4) =  -0.209299217902479D0     
           gll8(5) =   0.209299217902479D0     
           gll8(6) =   0.591700181433142D0     
           gll8(7) =   0.871740148509607D0     
           gll8(8) =   1.000000000000000D0   

           gll9(1) = -1.000000000000000D0     
           gll9(2) = -0.899757995411460D0     
           gll9(3) = -0.677186279510738D0     
           gll9(4) = -0.363117463826178D0     
           gll9(5) =  0.000000000000000D0
           gll9(6) =  0.363117463826178D0     
           gll9(7) =  0.677186279510738D0     
           gll9(8) =  0.899757995411460D0     
           gll9(9) =  1.000000000000000D0     
 
           gll10(1) = -1.00000000000000D0     
           gll10(2) = -0.919533908166459D0     
           gll10(3) = -0.738773865105505D0     
           gll10(4) = -0.477924949810444D0     
           gll10(5) = -0.165278957666387D0     
           gll10(6) =  0.165278957666387D0     
           gll10(7) =  0.477924949810444D0     
           gll10(8) =  0.738773865105505D0     
           gll10(9) =  0.919533908166459D0     
           gll10(10) =  1.00000000000000D0     

           gll11(1) = -1.00000000000000D0     
           gll11(2) = -0.934001430408059D0     
           gll11(3) = -0.784483473663144D0     
           gll11(4) = -0.565235326996205D0     
           gll11(5) = -0.295758135586939D0     
           gll11(6) =  0.000000000000000D0
           gll11(7) =  0.295758135586939D0     
           gll11(8) =  0.565235326996205D0     
           gll11(9) =  0.784483473663144D0     
          gll11(10) =  0.934001430408059D0     
          gll11(11) =  1.00000000000000D0     

 
         SelectCase(np)
          Case(3) 
             gllg(:) = gll3(:)
          Case(4) 
             gllg(:) = gll4(:)
          Case(5) 
             gllg(:) = gll5(:)
          Case(6) 
             gllg(:) = gll6(:)
          Case(7) 
             gllg(:) = gll7(:)
          Case(8) 
             gllg(:) = gll8(:)
          Case(9) 
             gllg(:) = gll9(:)
          Case(10) 
             gllg(:) = gll10(:)
          Case(11) 
             gllg(:) = gll11(:)
         End select

         
      End Function  Gll_grids
!============================================================
      Subroutine  gl_pts_wts(np,glpt,glwt)

      Implicit None
      Integer, Intent(in) :: np
      Real(Kind=real_kind),  Dimension(np), Intent(out)  :: glpt, glwt
      Real(Kind=real_kind) :: gl2(2),gl3(3),gl4(4),gl5(5),gl6(6),gl7(7), &
                           gw2(2),gw3(3),gw4(4),gw5(5),gw6(6),gw7(7)

        if (np > 7 )  then
            print*, 'No GL  grid availabe'
        endif
!============================================================
!                 GL(pts)                GL(wts) 
  !!GL-2 
  gl2(1) = -0.577350269189626D0
  gl2(2) =  0.577350269189626D0
  gw2(1) =  1.000000000000000D0
  gw2(2) =  1.000000000000000D0

  !!GL-3 
  gl3(1) = -0.774596669241483D0
  gl3(2) =  0.000000000000000D0
  gl3(3) =  0.774596669241483D0

  gw3(1) =  0.555555555555556D0
  gw3(2) =  0.888888888888889D0
  gw3(3) =  0.555555555555556D0

  !!GL-4 
   gl4(1) = -0.861136311594053D0
   gl4(2) = -0.339981043584856D0
   gl4(3) =  0.339981043584856D0
   gl4(4) =  0.861136311594053D0

   gw4(1) =  0.347854845137454D0
   gw4(2) =  0.652145154862546D0
   gw4(3) =  0.652145154862546D0
   gw4(4) =  0.347854845137454D0

  !!GL-5 
   gl5(1) = -0.906179845938664D0
   gl5(2) = -0.538469310105683D0
   gl5(3) =  0.000000000000000D0
   gl5(4) =  0.538469310105683D0
   gl5(5) =  0.906179845938664D0

   gw5(1) =  0.236926885056189D0
   gw5(2) =  0.478628670499366D0
   gw5(3) =  0.568888888888889D0
   gw5(4) =  0.478628670499366D0
   gw5(5) =  0.236926885056189D0
  !!GL-6 
  gl6(1) = -0.932469514203152D0
  gl6(2) = -0.661209386466264D0
  gl6(3) = -0.238619186083197D0
  gl6(4) =  0.238619186083197D0
  gl6(5) =  0.661209386466264D0
  gl6(6) =  0.932469514203152D0

  gw6(1) =  0.171324492379170D0
  gw6(2) =  0.360761573048139D0
  gw6(3) =  0.467913934572691D0
  gw6(4) =  0.467913934572691D0
  gw6(5) =  0.360761573048139D0
  gw6(6) =  0.171324492379170D0

  !!GL-7 
  gl7(1) = -0.949107912342758D0
  gl7(2) = -0.741531185599394D0
  gl7(3) = -0.405845151377397D0
  gl7(4) =  0.000000000000000D0
  gl7(5) =  0.405845151377397D0
  gl7(6) =  0.741531185599394D0
  gl7(7) =  0.949107912342758D0

  gw7(1) =  0.129484966168870D0
  gw7(2) =  0.279705391489277D0
  gw7(3) =  0.381830050505119D0
  gw7(4) =  0.417959183673469D0
  gw7(5) =  0.381830050505119D0
  gw7(6) =  0.279705391489277D0
  gw7(7) =  0.129484966168870D0

         SelectCase(np)
          Case(2)
             glpt(:) = gl2(:)
             glwt(:) = gw2(:)
          Case(3)
             glpt(:) = gl3(:)
             glwt(:) = gw3(:)
          Case(4)
             glpt(:) = gl4(:)
             glwt(:) = gw4(:)
          Case(5)
             glpt(:) = gl5(:)
             glwt(:) = gw5(:)
          Case(6)
             glpt(:) = gl6(:)
             glwt(:) = gw6(:)
          Case(7)
             glpt(:) = gl7(:)
             glwt(:) = gw7(:)
        end select
      End Subroutine  gl_pts_wts
                                                            
!--------------------------------------------------------

    End Module  Interpol_mod  
