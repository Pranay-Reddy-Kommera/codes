!-------------------------------------------------------------------
! Defines the test case to be solved:
! * sets the duration and spatial extent of the simulation
! * sets the initial data for the simulation
! * updates the advection velocities, if need be
!
! This module provides a single entry point into the various "work"
! routines for the different testcases. This is as close to OOP as I
! know how with Fortran...
!-------------------------------------------------------------------

module testcases_nh3d_mod
  use basic_mod
  use element_mod
  use physical_const_mod
  use cudafor

  implicit none
  private
  public :: load_params, initialize_testcase, print_vital_stat
  public :: Compute_Pressure, Compute_SoundSpeed,Pressure,SoundSpeed
  private :: nh3d_warm_bubble, nh3d_igw , Load_NHvars 
  private :: Iterative_z2eta, nh3d_steady_state

contains


!!==========================================================================================
  Subroutine print_vital_stat()
  real (kind=real_kind) :: dx, dy, dz, min_len 

        dx = delx/dble(nx-1)
        dy = dely/dble(nx-1)
        dz = delz/dble(nx-1)
       min_len = floor(min(dx,dy,dz)) 

      print*, " HOMAM_Cart: Vital Statistics ------------------------- "
      if (use_dim_split_scheme) then
            print*, "     Split 2D+1D Euler DG (nodal) Solver"
      else 
            print*, "     Fully 3D Euler DG (nodal) Solver"
      endif 

      print*, " Name of the Test               : ", name_of_the_test 
      print*, " Degree of the poly in DG-Space : ", nod              
      print*, " Number of Elements in (x,y,z)  : ", nex,ney,nez
      print*, " Element width del(x,y,z) in (m): ",  floor(delx),floor(dely), floor(delz)
      print*, " Grid Aspect Ratio, (dx,dy)/dz  : ", floor(delx/delz), floor(dely/delz)
      print*, " Global min grid-spacing in  (m): ",  min_len
      print('(a,f8.4)'),  "  Time step in seconds         : ",  dt         
      print*, " -------------------------------------------------- "

  End Subroutine print_vital_stat

!!==========================================================================================
  ! Wrapper for setting the duration and spatial extents in basic_mod with
  ! specific values defined for each testcase
  Subroutine load_params()
    select case(testcase)
    case (testcase_id_steady_state)
      name_of_the_test = "steady_state"
      xmin = 0.0D0 
      xmax = 40000.0D0 *1000.0D0 
      ymin = 0.0D0 
      ymax = 6000.0D0 *1000.0D0 
      zmin = 0.0D0
      zmax = 30.0D0 *1000.0D0 
    case (testcase_id_warm_bubble)
      name_of_the_test = "warm bubble"
      xmin = 0.0D0 
      xmax = 10000.0D0 *2.5D0
      ymin = 0.0D0 
      ymax = 10000.0D0 *2.5D0
      zmin = 0.0D0
      zmax = 10000.0D0 
    case (testcase_id_nh_igw)
      name_of_the_test = "nh_igw"
      xmin = 0.0D0 
      xmax = 320000.0D0 
      ymin = 0.0D0 
      ymax = 160000.0D0 
      zmin = 0.0D0
      zmax = 10000.0D0 
   !print*, "domain size: ", xmax, ymax, zmax
    case default
      print*, "unrecognized test case, can't load params!"
      stop
    end select
  End Subroutine load_params

!!==========================================================================================

  ! Wrapper for setting the initial values of psi,u,v,w for each testcase
  Subroutine initialize_testcase()
!    type(grid_t), intent(in), dimension(:,:,:) :: grid
!    type(element_t), intent(out), dimension(:,:,:) :: vars

    ! set metric values to identity
    ! testcases with deformed grids must then overwrite
  ! call initialize_identity_metric(grid, metric)

    select case(testcase)
    case (testcase_id_warm_bubble)
      call nh3d_warm_bubble()
    case (testcase_id_nh_igw)
      call nh3d_igw()
    case (testcase_id_steady_state)
      call nh3d_steady_state()
    case default
      print*, "unrecognized test case, can't compute initial data!"
      stop
    end select
  End Subroutine initialize_testcase

!!==========================================================================================
  Subroutine nh3d_igw()
   Implicit none 
!    type(grid_t), intent(in), dimension(:,:,:) :: grid
!    type(element_t), intent(out), dimension(:,:,:) :: vars

    real (kind=real_kind) :: prs,prb,spd,prp,the,thb,thp,rhb, pib
    real (kind=real_kind) :: rho,rhp,rtp,rtb, u,v,w , cori 
    real (kind=real_kind) :: dist, x,y,z, junk(nx,nx,nx)  
    integer :: nlx, nly, nlz, ie, je, ke, i,j,k

    real(kind=real_kind), parameter :: xc = 100000.0D0
    real(kind=real_kind), parameter :: yc = 80000.0D0   !ymax/2 
    real(kind=real_kind), parameter :: aa = 5000.0D0
    real(kind=real_kind), parameter :: zt = 10000.0D0    
    real(kind=real_kind), parameter :: t0 = 0.1D0


    if (size(grid) /= size(vars)) stop
    nlx = size(grid,1)
    nly = size(grid,2)
    nlz = size(grid,3)

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
          
          cori = 0.0D0

          do k = 1, nx 
          do j = 1, nx 
          do i = 1, nx 
              x = grid(ie,je,ke)%x(i,j,k) 
              y = grid(ie,je,ke)%y(i,j,k) 
              z = grid(ie,je,ke)%z(i,j,k) 

              u = 20.0D0 
              v = 0.0D0 
              w = 0.0D0 

              thb  = T_0 * exp(N_0*N_0 * z/grv)

              dist =   1.0D0 + ((x - xc)/aa)**2 + ((y - yc)/aa)**2 
              thp = t0  * (sin(pi*z /zt))/ dist 
              the = thb + thp 

                pib = 1.0D0 + (grv*grv / (C_p* T_0* N_0*N_0))*(T_0 - thb)/thb

              prb  = P_0 * pib**(Cp_Rd)  !Exner relation
              rhb = (prb / C_0 )**(1.0D0/Cp_Cv)  / thb   !Eqn of state 
              rho = (prb / C_0 )**(1.0D0/Cp_Cv)  / the 

              rhp =  rho - rhb
              rtb =  rhb * thb
              rtp =  rho * the - rtb

                prs = Compute_Pressure(rho,the)
                spd = Compute_SoundSpeed(prs,rho)
              cori = 0.0D0 

             Call Load_NHvars(ie,je,ke,i,j,k,u,v,w,the,thb,thp,prs,prb,spd,rho,rhb,rhp,rtb,rtp,cori)

            end do
            end do
            end do

          end do
        end do
      end do

  End Subroutine nh3d_igw

!!==========================================================================================
  Subroutine nh3d_steady_state()
   Implicit none 
!    type(grid_t), intent(in), dimension(:,:,:) :: grid
!    type(element_t), intent(out), dimension(:,:,:) :: vars

    real (kind=real_kind) :: prs,prb,spd,prp,the,thb,thp,rhb, pib
    real (kind=real_kind) :: rho,rhp,rtp,rtb, u,v,w, cori  
    real (kind=real_kind), dimension(nx,nx,nx)  ::  u_el,rho_el,the_el, fcori, prb_el,rhb_el, thb_el 
    integer :: nlx, nly, nlz, ie, je, ke, i,j,k

    if (size(grid) /= size(vars)) stop
    nlx = size(grid,1)
    nly = size(grid,2)
    nlz = size(grid,3)

       cori = 2.0D0 *omg_earth * sin(pi/4.0D0) 

   do ke = 1, nlz
     do je = 1, nly
       do ie = 1, nlx
         ! For each element 
        Call Iterative_z2eta(ie,je,ke,u_el,rho_el,the_el,rhb_el,thb_el,prb_el,fcori) 

        do k = 1, nx
        do j = 1, nx
        do i = 1, nx
            u = u_el(i,j,k)            
            v = 0.0D0
            w = 0.0D0

            thb  = thb_el(i,j,k)                   
            the  = the_el(i,j,k)                   
            thp = the- thb 

            rho = rho_el(i,j,k)          
            rhb = rhb_el(i,j,k) 

            rhp =  rho - rhb
            rtb =  rhb * thb
            rtp =  rho * the - rtb

            prb = prb_el(i,j,k) 

            prs = prb !Compute_Pressure(rho,the)
            spd = Compute_SoundSpeed(prs,rho)

           cori = fcori(i,j,k) 

           Call Load_NHvars(ie,je,ke,i,j,k,u,v,w,the,thb,thp,prs,prb,spd,rho,rhb,rhp,rtb,rtp,cori)

        end do
        end do
        end do

      end do 
    end do 
  end do 

  End Subroutine nh3d_steady_state

!!==========================================================================================
  Subroutine Iterative_z2eta(ie,je,ke, u,rho,the,rhb,thb,prb,fcori)
   Implicit none 
!    type(grid_t), intent(in):: grids
integer, intent(in) :: ie,je,ke
    real (kind=real_kind), intent(out), dimension(nx,nx,nx)   :: u, rho, the ,fcori, prb,rhb, thb 
    real (kind=real_kind) :: prs,trm1,trm2,trm3, Phi, Temp, phip 
    real (kind=real_kind) :: eta_tol,eta_tmp,eta_old,eta_log,eta_err
    real (kind=real_kind) :: x,y,z, t0,u0,y0,f0,b0,bb, gama, phi_bar,t_bar 
    real (kind=real_kind) :: F_eta, dF_eta, sin_t, eta(nx,nx,nx)  
    integer ::  it,itn, i,j,k

           bb = 5.0D0 !! different from U&J (mwr 2012) paper  
           y0 = ymax / 2.0D0 
           u0 = 35.0D0 
           f0 = 2.0D0 * omg_earth * sin(pi/4.0D0)
           b0 = 0.0D0  
           !b0 = 2.0D0 / r_earth * omg_earth * cos(pi/4.0D0) !for the beta plane 
         gama = 0.005 
           t0 = 288.0D0 
      eta_old = 1.0D-7 !! U&J choice --??  
      eta_tol = 1.0D-14
          itn = 25 !default iterations  

      eta_tmp = 0.007 !! U&J paper seems to be wrong, no convergence with their params 

  ! For a given element 
     do k = 1, nx 
      do j = 1, nx 
       do i = 1, nx 
         x = grid(ie,je,ke)%x(i,j,k) 
         y = grid(ie,je,ke)%y(i,j,k) 
         z = grid(ie,je,ke)%z(i,j,k) 

         eta_err = 1.0D0  
         it = 0
      do while (eta_err > eta_tol)   !!iteratively find "eta" for a given (x,y,z) 

         it = it + 1
         eta_old = eta_tmp 

         trm1 = eta_tmp**(R_d*gama/grv) 
         phi_bar = t0*grv/gama * (1.0D0 - trm1)  
           t_bar = t0*trm1  

         sin_t   = sin(2.0D0*pi*y/ymax) 
         eta_log = log(eta_tmp) 
         
         trm2 = (f0 - b0*y0)*(y - ymax/2.0D0 - ymax/(2.0D0*pi) * sin_t )
         trm3 = b0/2.0D0 * (y**2 -(ymax*y/pi) * sin_t - ymax*ymax/(2.0D0 * pi*pi)) * &
                           (cos(2.0D0*pi*y/ymax) + ymax*ymax/3.0D0) 

         phip = u0*0.5D0 * (trm2 + trm3) 
         Phi =  phi_bar + phip * eta_log * exp(-(eta_log/ bb)**2)
         Temp = t_bar + phip/R_d * (2.0D0 *(eta_log/bb)**2 -1.0D0) * exp(-eta_log/bb) 

         F_eta = -grv *z + Phi 
         dF_eta = -R_d/eta_tmp * Temp 

         eta_tmp = eta_tmp - F_eta / dF_eta  
         eta_err = abs(eta_old - eta_tmp) 

        !print*, it, eta_tmp

         if (it  > itn ) then
            print*, 'eta convergence failed '
            stop 
         endif 

      end do !while 

        !print*, i,j,k, eta_tmp 

       !! Compute u,rho,the with converged "eta" 

        fcori(i,j,k)  = 0.0D0 !f0   !Coriolis term in the source 

        u(i,j,k) = -u0* (sin(pi*y/ymax))**2 * eta_log * exp(-(eta_log/bb)**2) 
        prs = P_0 * eta_tmp 
        prb(i,j,k) = prs 
        rho(i,j,k) = prs/(R_d* Temp)
        rhb(i,j,k) = rho(i,j,k) !prs/(R_d* t_bar)

        the(i,j,k) = Temp * (P_0 / prs)**(Rd_Cp) 
        thb(i,j,k) = the(i,j,k) !t_bar * (P_0 / prs)**(Rd_Cp) 
       end do 
     end do 
   end do 

  End Subroutine Iterative_z2eta

!!==========================================================================================
  Subroutine nh3d_warm_bubble()
   Implicit none 
!    type(grid_t), intent(in), dimension(:,:,:) :: grid
!    type(element_t), intent(out), dimension(:,:,:) :: vars

    real (kind=real_kind) :: prs,prb,spd,prp,the,thb,thp,rhb, pib
    real (kind=real_kind) :: rho,rhp,rtp,rtb, u,v,w , cori 
    real (kind=real_kind) :: dist, x,y,z, junk(nx,nx,nx)  
    integer :: nlx, nly, nlz, ie, je, ke, i,j,k

    real(kind=real_kind), parameter :: xc = 5000.0D0*2.5D0
    real(kind=real_kind), parameter :: yc = 5000.0D0*2.5D0
    real(kind=real_kind), parameter :: zc = 2000.0D0
    real(kind=real_kind), parameter :: r0 = 2000.0D0
    real(kind=real_kind), parameter :: t0 = 2.0D0


    if (size(grid) /= size(vars)) stop
    nlx = size(grid,1)
    nly = size(grid,2)
    nlz = size(grid,3)

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
         
          cori = 0.0D0  

          do k = 1, nx 
          do j = 1, nx 
          do i = 1, nx 
             x = grid(ie,je,ke)%x(i,j,k) 
             y = grid(ie,je,ke)%y(i,j,k) 
             z = grid(ie,je,ke)%z(i,j,k) 

                 u = 0.0D0 
                 v = 0.0D0 
                 w = 0.0D0 

            dist =  sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2)
            thp = t0  * max(0.0D0,(1.0D0 - dist/r0))
            thb  = T_0      
                 
            the = thb + thp 
            pib = 1.0D0 - ( z* grv / (C_p* T_0))

            prb  = P_0 * pib**(Cp_Rd)  !Exner relation
            rhb = (prb / C_0 )**(1.0D0/Cp_Cv)  / T_0   !Eqn of state 

              rho =  rhb
              rhp =  rho - rhb

              rtb =  rhb * thb 
              rtp =  rho * the - rtb 
                
             prs = Compute_Pressure(rho,the)
             spd = Compute_SoundSpeed(prs,rho)

       Call Load_NHvars(ie,je,ke,i,j,k,u,v,w,the,thb,thp,prs,prb,spd,rho,rhb,rhp,rtb,rtp,cori)


          end do 
          end do 
          end do 

        end do
      end do
    end do

  end Subroutine nh3d_warm_bubble

!=======================================================================================================!
 Subroutine Load_NHvars(ie,je,ke,i,j,k,u,v,w,the,thb,thp,prs,prb,spd,rho,rhb,rhp,rtb,rtp,cori)
   Implicit none
!    type(element_t), intent(inout) :: var
    integer, intent(in)  :: i,j,k,ie,je,ke
    real(kind=real_kind), intent(in)  :: u,v,w,the,thb,thp,prs,prb,spd,rho,rhb,rhp,rtb,rtp ,cori 

             vars_init(ie,je,ke)%u(i,j,k)  = u
             vars_init(ie,je,ke)%v(i,j,k)  = v
             vars_init(ie,je,ke)%w(i,j,k)  = w

             vars_init(ie,je,ke)%thb(i,j,k)  = thb
             vars_init(ie,je,ke)%the(i,j,k)  = the
             vars_init(ie,je,ke)%thp(i,j,k)  = thp !display output           

             vars_init(ie,je,ke)%prb(i,j,k)  = prb
             vars_init(ie,je,ke)%prs(i,j,k)  = prs
             vars_init(ie,je,ke)%prp(i,j,k)  = prs  - prb
             vars_init(ie,je,ke)%spd(i,j,k)  = spd

             vars_init(ie,je,ke)%rhb(i,j,k)  = rhb
             vars_init(ie,je,ke)%rho(i,j,k)  = rho

             vars_init(ie,je,ke)%rtb(i,j,k)  = rtb
             vars_init(ie,je,ke)%rhp(i,j,k)  = rhp

             vars_init(ie,je,ke)%fcori(i,j,k)  = cori 

           !Store initial State vecor for Euler-3d system  
             vars_init(ie,je,ke)%vec(i,j,k,1)  = rhp
             vars_init(ie,je,ke)%vec(i,j,k,2)  = rho * u
             vars_init(ie,je,ke)%vec(i,j,k,3)  = rho * v
             vars_init(ie,je,ke)%vec(i,j,k,4)  = rho * w
             vars_init(ie,je,ke)%vec(i,j,k,5)  = rtp

 End Subroutine Load_NHvars

!=======================================================================================================!
  Function  Compute_SoundSpeed(prs,rho) result(speed)
    Implicit None
    real (kind=real_kind) :: rho, speed, prs

     speed =  sqrt(Cp_Cv * prs / rho)

   End Function Compute_SoundSpeed

!=======================================================================================================!
  Function  Compute_Pressure(rho,the) result(prs)
    Implicit None
    real (kind=real_kind) :: rho, the, prs
    real (kind=real_kind) :: tmp

     prs =  C_0 * (rho * the)**Cp_Cv
!     tmp = rho * the
!     prs =  C_0 * exp(Cp_Cv * log(tmp))
!      tmp = dlog(C_0) + (Cp_Cv * dlog(rho * the))
!      prs = dexp(tmp)

   End Function Compute_Pressure

!=======================================================================================================!

!=======================================================================================================!
!  attributes(device) Function  SoundSpeed(prs,rho) result(speed)
  attributes(device) Subroutine SoundSpeed(prs,rho,speed)
    Implicit None
    real (kind=real_kind), intent(in) :: rho, prs
    real (kind=real_kind), intent(out) :: speed

     speed =  sqrt(Cp_Cv * prs / rho)

!   End Function SoundSpeed
   End Subroutine SoundSpeed

!=======================================================================================================!
!  attributes(device) Function  Pressure(rho,the) result(prs)
  attributes(device) Subroutine Pressure(rho,the,prs)
    Implicit None
    real (kind=real_kind), intent(in) :: rho, the
    real (kind=real_kind), intent(out) :: prs
    real (kind=real_kind) :: tmp

!     prs =  C_0 * (rho * the)**Cp_Cv
      tmp = dlog(C_0) + (Cp_Cv * dlog(rho * the))
      prs = dexp(tmp)

   End Subroutine Pressure
!   End Function Pressure

!=======================================================================================================!


End module testcases_nh3d_mod
