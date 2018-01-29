!------------------------------------------------
! Gauss Qudarture (GLL-pts) and derivative matrix
! Legendre polynomials and derivatives
!------------------------------------------------
! R.Nair NCAR/scd 10/02
! Minor style changes 06/2016 by Francois Hebert, SIParCS/Cornell
!
! Gauss_Lobatto points [-1,1] and corresponding weights
! Reference book : Cantuo, Hussaini et al. Pg 524
!
! Note: xpt, dpt stores all the polynomial value up to degree "nod"
!   Roots of Lagendre polynomial via Jacobi polynomials
!   der ->  coefficients of derivative computation at nodes
!   weakder ->  (der_ij w_i/w_j) as shows up in weak-form DG
!------------------------------------------------

module gauss_quadrature_mod
  use basic_mod

  implicit none
  private
  public :: gauss_lobatto
  public :: legendre_poly
  public :: jacobi_pts

  real(kind=real_kind), public, dimension(nx) :: gllp, gllw
  real(kind=real_kind), public, dimension(nx,nx) :: der, weakder
  Real(Kind=real_kind), Public, Dimension(nx,nx) :: pmx,dpm,der_gll

      ! GLL or GL  depending on initialization 
       Real(Kind=real_kind), Public, Dimension(nx) :: gl, gw
       Real(Kind=real_kind), Public, Dimension(nx,nlg) :: hm_lg   !interplolation matrix for GLL-GL 

       Real(Kind=real_kind), Public, Dimension(nx,nx) :: mmx, mmxi

      ! GL specific 
       Real(Kind=real_kind), Public, Dimension(nx) :: lgp,lgw   !GL points & weights 
       Real(Kind=real_kind), Public, Dimension(nx,nx) :: der_gl !Derivative matrix  
       Real(Kind=real_kind), Public, Dimension(nx,2) :: ime_gl  !interpolation edge matrix 
       Real(Kind=real_kind), Public, Dimension(nx,2) :: ime_der_gl  !interpolation edge derivative matrix 
       Real(Kind=real_kind), Public, Dimension(nx,2) :: hed_gl  !interpolated  edge LG-basis function


contains
!!=================================================================================


       Subroutine  Gauss_Lobatto(ngp,gl,gw,der)

        Implicit None

        Integer, intent(in) :: ngp
        Real(Kind=real_kind), Intent(out), Dimension(ngp)  :: gl, gw
        Real(Kind=real_kind), Intent(out), Dimension(ngp,ngp)  :: der

        Real(Kind=real_kind), Dimension(0:ngp)  :: xj, glw, dxx, xjr
        Real(Kind=real_kind), Dimension(0:ngp)  :: xpt, dpt, xpt1,dpt1
        Real(Kind=real_kind), Dimension(0:ngp,0:ngp)  :: ddr

        Real(Kind=real_kind) :: alpha, beta
        Real(Kind=real_kind) :: rr, sm, eps, djp, xjp,dl
        Real(Kind=real_kind) :: dth,dd, aa,bb, cd,cs,ss,sd,tmp
        Real(Kind=real_kind) :: c0,c1,c2,c4

        Integer ::   k ,i,j ,itn , nh, npt, np1, deg


          itn = 35            !Max iteration number
          eps = 1.0D-18        !Tolerence

           c0 = 0.0D0 
           c1 = 1.0D0 
           c2 = 2.0D0
           c4 = 4.0D0 

         alpha = c0
          beta = c0

           npt = ngp
           np1 = ngp
           deg = ngp - 1


           xj(0) = c1               !known boundary values of GL pts [1,-1]
         xj(deg) = -c1

              nh = (deg+1)/2


       call Jacobi_pts(np1,npt,alpha,beta, c1,xpt, dpt)
       call Jacobi_pts(np1,npt,alpha,beta,-c1,xpt1,dpt1)

           dd =  xpt(deg)*xpt1(deg-1) - xpt1(deg)*xpt(deg-1)
           aa = (xpt1(np1)*xpt(deg-1) - xpt1(deg-1)*xpt(np1) )/dd
           bb = (xpt1(deg)*xpt(np1) - xpt1(np1)*xpt(deg) )/dd

       ! Initial guess

          dth = pi/dble(2*deg+1)
           cd = cos(c2*dth)
           sd = sin(c2*dth)
           cs = cos(dth)
           ss = sin(dth)

       do  k = 1, nh-1
               rr = cs
               dl = c1

               do i = 1, itn            !Newton-Raphson iterations

              call Jacobi_pts(np1,npt,alpha,beta,rr,xpt,dpt)

                xjp =  xpt(np1)+aa* xpt(deg)+bb* xpt(deg-1)
                djp =  dpt(np1)+aa* dpt(deg)+bb* dpt(deg-1)

                    sm = c0
                  do j = 0, k-1
                  sm = sm + c1 /(rr - xj(j))
                  enddo

                dl = -(xjp /(djp - xjp*sm))
               rr = rr + dl
                     ! print*, i
                  if (abs(dl)  < eps) Exit
              enddo

           xj(k) = rr

       !Initial guess for the next point

          tmp = cs*cd-ss*sd
          ss = cs*sd+ss*cd
          cs = tmp
       enddo

          do  k = 1, nh
            xj(deg-k) = -xj(k)    !roots by symmetry
          enddo
          if (mod(deg,2) == 0) xj(nh) = c0

      !Now Weight for the Gauss_Lobatto points

        do k = 0, deg
           rr = xj(k)
         call Jacobi_pts(deg,npt,alpha,beta,rr,xpt,dpt)
                   xjp = xpt(deg)
         glw(k) = c2 / (dble(deg*(deg+1)) * xjp*xjp)
         dxx(deg-k) = xjp
        enddo

        do j = 0, deg
          xjr(deg-j) = xj(j)      ![-1, 1] variation
          gl(deg+1-j) = xj(j)
          gw(deg+1-j) = glw(j)
        enddo

      !Derivative coefficients

         do i = 0, deg
          do j = 0, deg
             if ((i ==0).and.(j ==0)) then
               ddr(i,j) = - dble(deg*(deg +1)) /c4
              elseif ((i == deg).and.(j==deg)) then
               ddr(i,j) =   dble(deg*(deg +1)) /c4
              else

                if ( i /= j) then
                 ddr(i,j) = dxx(i)/ (dxx(j) *(xjr(i) - xjr(j)) )
                elseif ( (i == j).and.((i <= deg-1).and.(j <= deg-1)) ) then
                 ddr(i,j) = c0
                endif

            endif
              der(i+1,j+1) = ddr(i,j)
          enddo
         enddo

          sm = c0

      end    Subroutine gauss_lobatto

!!------------------------------------------------
 !R.Nair NCAR/scd 04/03
 !  Legendre ploynomials and derivatives
!!------------------------------------------------
       Subroutine  Legendre_Poly(ngp,gl,pmx,dpm)

        Implicit None

        Integer, intent(in) :: ngp
        Real(Kind=real_kind), Intent(in), Dimension(ngp) :: gl
        Real(Kind=real_kind), Intent(out), Dimension(ngp,ngp) :: pmx,dpm

        Real(Kind=real_kind), Dimension(0:ngp) :: px, dx, cm
        Real(Kind=real_kind), Dimension(ngp,ngp) :: dmy,chk, dsp
        Real(Kind=real_kind), Dimension(0:ngp-1,ngp) :: cx, qx
        Real(Kind=real_kind) :: alpha, beta, xi

        Integer ::   i,j,k,l, n1,n2


             alpha = 0.0D0  !Set for Lagendre poly
              beta = 0.0D0

               pmx = 0.0D0
               dpm = 0.0D0

 ! Compute Legendre polynomials of degree less than or equal to "nod"
 !   store it in "nod+1" positions starting from index 1, for each
 !   abcissa gl(i), i=1,..,nx
 !
        do i= 1, nx
              xi = gl(i)

           Call Jacobi_pts(ngp-1,ngp,alpha,beta,xi,px,dx)

           do l = 0, ngp-1
                cx(l,i) = px(l)    !Legendre Poly of order "l"
             pmx(l+1,i) = px(l)
             dpm(l+1,i) = dx(l)
           enddo
        enddo


      End    Subroutine Legendre_Poly

  !----------------------------------------------------
  ! Evaluates the Jacobi polynomials at "x" of the order
  ! up to "nox". Derivatives are also calculated,
  ! px(nox), dpx(nox) are the function  value and its drivative
  !----------------------------------------------------

  subroutine jacobi_pts(nox,npt,alpha,beta,x,px,dpx)

    integer, intent(in) :: nox,npt
    real(kind=real_kind), intent(in) :: alpha, beta, x
    real(kind=real_kind), intent(inout), dimension(0:npt) :: px, dpx

    real(kind=real_kind), dimension(0:npt) :: jpx, djpx
    real(kind=real_kind) :: c1,c2,c0, k1,dk, ab, ac1,ac2,ac3,acx
    integer :: k

    c0 = 0.0D0
    c1 = 1.0D0
    c2 = 2.0D0

    do k = 0, nox
      px(k) = c0
      dpx(k) = c0
    end do

    if (nox == 0) then
      px(nox) = c1
      dpx(nox) = c0

    else if (nox == 1) then
      px(0) = c1
      px(1) =  (c1 + alpha)*x
      dpx(0) = c0
      dpx(1) = c1 + alpha

    else                             ! general cases
      ab = alpha + beta

      jpx(0) = c1
      jpx(1) = (c1 + alpha)*x
      djpx(0) = c0
      djpx(1) = (c1 + alpha)

      ! recursive process
      do k = 1, nox-1
        dk = k
        k1 = k+1

        ac1 =  c2*k1* (dk+ ab +c1) * (c2*dk + ab)
        ac2 =  (c2*k1+ ab) * (c2*dk +ab + c1) * (c2*dk + ab)
        acx =  (c2*dk + ab + c1) * (alpha*alpha - beta*beta) + x *ac2
        ac3 =  c2*(dk + alpha) * (dk + beta) * (c2*k + ab + c2)

        jpx(k+1) = (acx*jpx(k) - ac3*jpx(k-1)) / ac1
        djpx(k+1) = (acx*djpx(k) + ac2*jpx(k) - ac3*djpx(k-1)) / ac1
      end do

    end if

    if (nox > 1) then
      do k = 0, nox          ! The last value at "nox" is the desired
        px(k) = jpx(k)       ! point value  (of degree "nox")
        dpx(k) = djpx(k)
      end do
    end if

  end subroutine jacobi_pts

end module gauss_quadrature_mod
