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

!  real(kind=real_kind), public, dimension(nx) :: gllp, gllw
!  real(kind=real_kind), public, dimension(nx,nx) :: der, weakder

contains

  subroutine gauss_lobatto(gllp,gllw,der)

    real(kind=real_kind), intent(out), dimension(nod+1) :: gllp, gllw
    real(kind=real_kind), intent(out), dimension(nod+1,nod+1) :: der

    real(kind=real_kind), dimension(0:nod+1) :: xj, glw, dxx, xjr
    real(kind=real_kind), dimension(0:nod+1) :: xpt, dpt, xpt1, dpt1
    real(kind=real_kind), dimension(0:nod+1,0:nod+1) :: ddr

    real(kind=real_kind) :: alpha, beta
    real(kind=real_kind) :: rr, sm, eps, djp, xjp, dl
    real(kind=real_kind) :: dth, dd, aa, bb, cd, cs, ss, sd, tmp
    real(kind=real_kind) :: c0,c1,c2,c4

    integer :: k,i,j,itn, nh, npt, np1

    itn = 35       ! max iteration number
    eps = 1.0D-18  ! tolerence

    c0 = 0.0D0
    c1 = 1.0D0
    c2 = 2.0D0
    c4 = 4.0D0

    alpha = c0
    beta = c0

    npt = nod + 1
    np1 = nod + 1

    xj(0) = c1 ! known boundary values of GL pts [1,-1]
    xj(nod) = -c1

    nh = (nod+1)/2

    call jacobi_pts(np1,npt,alpha,beta, c1,xpt, dpt)
    call jacobi_pts(np1,npt,alpha,beta,-c1,xpt1,dpt1)

    dd =  xpt(nod)*xpt1(nod-1) - xpt1(nod)*xpt(nod-1)
    aa = (xpt1(np1)*xpt(nod-1) - xpt1(nod-1)*xpt(np1) )/dd
    bb = (xpt1(nod)*xpt(np1) - xpt1(np1)*xpt(nod) )/dd

    ! initial guess

    dth = pi/(c2*nod+1)
    cd = cos(c2*dth)
    sd = sin(c2*dth)
    cs = cos(dth)
    ss = sin(dth)

    do  k = 1, nh-1
      rr = cs
      dl = c1

      do i = 1, itn ! Newton-Raphson iterations

        call Jacobi_pts(np1,npt,alpha,beta,rr,xpt,dpt)

        xjp =  xpt(np1)+aa* xpt(nod)+bb* xpt(nod-1)
        djp =  dpt(np1)+aa* dpt(nod)+bb* dpt(nod-1)

        sm = c0
        do j = 0, k-1
          sm = sm + c1 /(rr - xj(j))
        end do

        dl = -(xjp /(djp - xjp*sm))
        rr = rr + dl
        ! print*, i
        if (abs(dl)  < eps) Exit
      end do

      xj(k) = rr

      ! initial guess for the next point

      tmp = cs*cd-ss*sd
      ss = cs*sd+ss*cd
      cs = tmp
    end do

    do  k = 1, nh
      xj(nod-k) = -xj(k) ! roots by symmetry
    end do
    if (mod(nod,2) == 0) xj(nh) = c0

    ! now weights for the Gauss-Lobatto points

    do k = 0, nod
      rr = xj(k)
      call Jacobi_pts(nod,npt,alpha,beta,rr,xpt,dpt)
      xjp = xpt(nod)
      glw(k) = c2 / (nod*(nod+1) * xjp*xjp)
      dxx(nod-k) = xjp
    end do

    do j = 0, nod
      xjr(nod-j) = xj(j) ! [-1, 1] variation
      gllp(nod+1-j) = xj(j)
      gllw(nod+1-j) = glw(j)
    end do

    ! Derivative coefficients

    do i = 0, nod
      do j = 0, nod
        if ((i ==0).and.(j ==0)) then
          ddr(i,j) = - nod*(nod+1)/c4
        else if ((i == nod).and.(j==nod)) then
          ddr(i,j) =   nod*(nod+1)/c4
        else

          if ( i /= j) then
            ddr(i,j) = dxx(i)/ (dxx(j) *(xjr(i) - xjr(j)) )
          else if ( (i == j).and.((i <= nod-1).and.(j <= nod-1)) ) then
            ddr(i,j) = c0
          end if

        end if
        der(i+1,j+1) = ddr(i,j)
        weakder(i+1,j+1) = der(i+1,j+1) * gllw(i+1) / gllw(j+1)
      end do
    end do

    sm = c0

  end subroutine gauss_lobatto

  !------------------------------------------------
  ! R.Nair NCAR/scd 04/03
  ! Legendre ploynomials and derivatives
  !------------------------------------------------

  subroutine legendre_poly(gllp,pmx,dpm)

    real(kind=real_kind), intent(in), dimension(nx) :: gllp
    real(kind=real_kind), intent(out), dimension(nx,nx) :: pmx,dpm

    real(kind=real_kind), dimension(0:nx) :: px, dx
    real(kind=real_kind), dimension(0:nx-1,nx) :: cx
    real(kind=real_kind) :: alpha, beta, xi

    integer :: i,l

    alpha = 0 ! set for Lagendre poly
    beta = 0

    pmx = 0
    dpm = 0

    ! Compute Legendre polynomials of degree less than or equal to "nod"
    !   store it in "nod+1" positions starting from index 1, for each
    !   abcissa gllp(i), i=1,..,nx
    !
    do i= 1, nx
      xi = gllp(i)

      call jacobi_pts(nod,nx,alpha,beta,xi,px,dx)

      do l = 0, nod
        cx(l,i) = px(l)    ! Legendre Poly of order "l"
        pmx(l+1,i) = px(l)
        dpm(l+1,i) = dx(l)
      end do
    end do

  end subroutine legendre_poly

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
