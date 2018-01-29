!-------------------------------------------------------------------
! 04/04 R.Nair NCAR/SCD
! Rewritten 06/2016 by Francois Hebert, SIParCS/Cornell
!
! Build a simple cartesian grid of DG elements. 3D.
! Has functions to return several types of grids:
! a) the evolution grid for 'this' MPI process
! b) the 1D grids without duplicated points at element interfaces,
!    most useful for saving visualization data.
!
! Reminder on notation: the various constants used are,
!   nod = order of the polynomial
!   nx = nod+1 = number of GLL points
!   nex = number of elements in x
!   nux = number of unique x locations
!
! The physical domain spans [xmin,xmax] x [ymin,ymax] x [zmin,zmax]
! where the min/max variables are defined in basic_mod, and their
! values are set in each test case as appropriate.
!
!-------------------------------------------------------------------

module grid_setup_mod
  use basic_mod
  use gauss_quadrature_mod
  use element_mod

  implicit none
  private
  public :: make_grid, get_visualization_grid

  ! Local storage holding the 1D x/y/z grid locations.
  ! The 3D grids are constructed by the tensor product of these 1D arrays.
  !   global_x1d -> every grid point along this dimension
  !   viz_x1d    -> same, without dupes at interfaces; for visualization grid
  real(kind=real_kind), dimension(nx,nex) :: global_x1d
  real(kind=real_kind), dimension(nx,ney) :: global_y1d
  real(kind=real_kind), dimension(nx,nez) :: global_z1d
  real(kind=real_kind), dimension(nux) :: viz_x1d
  real(kind=real_kind), dimension(nuy) :: viz_y1d
  real(kind=real_kind), dimension(nuz) :: viz_z1d

  logical :: initialize_1d_grids_called = .false.

contains

  ! Computes the x,y,z coords of every point in this MPI process' block of elements.
  subroutine make_grid()
    use mpi_mod, only: offsetx, offsety

!    type(grid_t), intent(out), dimension(:,:,:) :: grid
    integer :: i, j, k, ie, je, ke, ie_global, je_global
    integer :: nlx, nly, nlz

    if (.not. initialize_1d_grids_called) call initialize_1d_grids()

    nlx = size(grid,1)
    nly = size(grid,2)
    nlz = nez

    ! Take tensor product to obtain full grid
    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
          ie_global = ie + offsetx()
          je_global = je + offsety()
          do k = 1, nx
            do j = 1, nx
              do i = 1, nx
                grid(ie,je,ke)%x(i,j,k) = global_x1d(i,ie_global)
                grid(ie,je,ke)%y(i,j,k) = global_y1d(j,je_global)
                grid(ie,je,ke)%z(i,j,k) = global_z1d(k,ke)
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine make_grid


!!====================================================

  ! Returns the de-duplicated 1D grid in each direction, for visualization.
  subroutine get_visualization_grid(ax, ay, az)
    real(kind=real_kind), intent(out), dimension(nux) :: ax
    real(kind=real_kind), intent(out), dimension(nuy) :: ay
    real(kind=real_kind), intent(out), dimension(nuz) :: az

    if (.not. initialize_1d_grids_called) call initialize_1d_grids()

    ax = viz_x1d
    ay = viz_y1d
    az = viz_z1d
  end subroutine get_visualization_grid

!!====================================================
  ! Wrapper to initialize the 1D coordinates in all three directions.
  ! Also initializes the GLL quantities: quadrature points, weights, derivs, etc
  subroutine initialize_1d_grids()
    integer :: k 
    if (initialize_1d_grids_called) return

    ! (1D) Gauss-Lobatto points and weights
    call gauss_lobatto(gllp, gllw, der)

!! For Boyd-Vandeeven filter 
       ! Call Legendre_Poly(nx,gllp,pmx,dpm)

       ! Call BV_filter_fn(bv_sigma)
       ! !precomputation of filter weights/coeff 
       !   do k = 1,nx
       !     pmx_bv(:,k) = pmx(:,k)*bv_sigma(:)
       !     pmx_gw(k,:) = pmx(k,:)*gllw(:)
       !    !print*, k, bv_sigma(k) 
       !   enddo


    ! Compute coordinates global_x1d,viz_x1d in each dimension
    call grid_1dwork(xmin, xmax, nex, nux, gllp, delx, global_x1d, viz_x1d)
    call grid_1dwork(ymin, ymax, ney, nuy, gllp, dely, global_y1d, viz_y1d)
    call grid_1dwork(zmin, zmax, nez, nuz, gllp, delz, global_z1d, viz_z1d)

    initialize_1d_grids_called = .true.
  end subroutine initialize_1d_grids


!!====================================================

  ! Initialize the coordinates for a 1D interval domain.
  ! This is the heavy lifting for the 3D case, which is then simply
  ! obtained as a tensor product of 1D grids.
  ! NOTE: Use a generic coordinate 's' in the variable names, because
  !       some compilers seem to poorly handle Fortran's scoping...
  !       exception: variable 'nx' which is used in all x/y/z directions.
  subroutine grid_1dwork(smin, smax, nes, nus, gllp, dels, global_s1d, viz_s1d)

    real(kind=real_kind), intent(in) :: smin, smax
    integer, intent(in) :: nes, nus
    real(kind=real_kind), intent(in), dimension(nx) :: gllp
    real(kind=real_kind), intent(out) :: dels
    real(kind=real_kind), intent(out), dimension(nx,nes) :: global_s1d
    real(kind=real_kind), intent(out), dimension(nus) :: viz_s1d

    real(kind=real_kind), dimension(nes+1) :: sedges
    integer :: i, k, ki

    ! widths of elements
    dels = (smax-smin) / nes

    ! boundaries of each 1d interval within [xmin, xmax]
    do k = 1, nes+1
      sedges(k) = smin + dels * (k-1)
    end do

    ! mapped values of GLL points
    do k = 1, nes
      do i = 1, nx
        ! set the full points (with dupes) for evolution
        global_s1d(i,k) = (sedges(k) + sedges(k+1) + dels * gllp(i)) / 2.0D0
        ! set the reduced points (without dupes) for visualization
        if (i==nx) exit ! skip the dupe
        ki = (k-1)*(nx-1) + i
        viz_s1d(ki) = global_s1d(i,k)
      end do
    end do
    viz_s1d(nus) = smax ! last point

  end subroutine grid_1dwork

end module grid_setup_mod

