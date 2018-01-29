!=======================================================================================!
!! Terrain Following Vertical coordinates (Gal-Chen & Sommerville, JCP 1975)
!! With option for Schar/Agnesi type topography
!! RDN, 07/12/2016
!! Changes 07/2016 by Francois Hebert, SIParCS/Cornell
!=======================================================================================!


Module mountain_grid_mod
  Use Basic_mod
  Use element_mod
  Use gauss_quadrature_mod
  Implicit None

  Private
  Public :: Compute_Zeta_Velocity, Compute_metric_terms
  Public :: Convert_vars_to_physical_frame
  Public :: Plot_mountain

  ! allocatable because MPI partition is determined at runtime
  type(metric_t), public, dimension(:,:,:), allocatable :: metric

Contains

  !!=========================================================================================
  !! Divide out the coordinate transformation Jacobian, to recover the
  ! physical-frame values. Useful to prepare values for output.
  Subroutine Convert_vars_to_physical_frame()
!    type(element_t), intent(in), dimension(:,:,:) :: vars
!    type(element_t), intent(out), dimension(:,:,:) :: vars_phys
    integer :: k, ie, je, ke, nlx, nly, nlz

    nlx = size(vars,1)
    nly = size(vars,2)
    nlz = size(vars,3)

    vars_phys = vars
    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
          do k = 1, nx
            vars_phys(ie,je,ke)%psi(:,:,k) = vars(ie,je,ke)%psi(:,:,k) / metric(ie,je,ke)%sg(:,:)
            vars_phys(ie,je,ke)%rho(:,:,k) = vars(ie,je,ke)%rho(:,:,k) / metric(ie,je,ke)%sg(:,:)
          end do
        end do
      end do
    end do
  End Subroutine Convert_vars_to_physical_frame


  !!=========================================================================================
  !! Compute values of metric terms over the grid (only within MPI partition)
  Subroutine Compute_metric_terms()
!    type(grid_t), intent(in), dimension(:,:,:) :: grids
    Real(Kind=real_kind) :: x, y, zeta
    Real(Kind=real_kind) :: tmp_mtn, tmp_ght, tmp_sg, tmp_sg13, tmp_sg23
    integer :: i,j,k, ie,je,ke, nlx,nly,nlz

    nlx = size(grid,1)
    nly = size(grid,2)
    nlz = size(grid,3)

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
          do k = 1, nx
            do j = 1, nx
              do i = 1, nx
                x = grid(ie,je,ke)%x(i,j,k)
                y = grid(ie,je,ke)%y(i,j,k)
                zeta = grid(ie,je,ke)%z(i,j,k)

                call Compute_metric_at_coords(x, y, zeta, &
                     tmp_mtn, tmp_ght, tmp_sg, tmp_sg13, tmp_sg23)

                ! mountain and coordinate heights
                metric(ie,je,ke)%mtn(i,j) = tmp_mtn
                metric(ie,je,ke)%ght(i,j,k) = tmp_ght

                ! metric terms
                metric(ie,je,ke)%sg(i,j) = tmp_sg
                metric(ie,je,ke)%sg13(i,j,k) = tmp_sg13
                metric(ie,je,ke)%sg23(i,j,k) = tmp_sg23
              end do
            end do
          end do
        end do
      end do
    end do

  End Subroutine Compute_metric_terms


  !!=========================================================================================
  !! Compute values of metric terms at a specific location (in computational frame)
  Subroutine Compute_metric_at_coords(x, y, zeta, mtn, ght, sg, sg13, sg23)
    Real(kind=real_kind), intent(in) :: x, y, zeta
    Real(kind=real_kind), intent(out) :: mtn, ght, sg, sg13, sg23

    Real(kind=real_kind) :: hxy, dh_x, dh_y, z_term

    select case(testcase)
    case (testcase_id_mountain_agnesi)
      Call Agnesi_Mountain_Surfdat(x,y,hxy,dh_x,dh_y)
    case(testcase_id_mountain_schar)
      Call Schar_Mountain_Surfdat(x,y,hxy,dh_x,dh_y)
    case default
      print*, "unrecognized mountain type"
      stop
    end select

    z_term = zeta/zmax - 1.0D0

    ! height of terrain at (x,y)
    mtn = hxy
    ! geometric height z for computational grid location (x,y,zeta)
    ght = mtn + zeta * (zmax - mtn)/zmax
    ! jacobian and metric terms
    sg = 1.0D0 - mtn/zmax
    sg13 = dh_x * z_term
    sg23 = dh_y * z_term

  End Subroutine Compute_metric_at_coords


  !!=========================================================================================
  Function Compute_Zeta_Velocity(ie, je, ke) result(wt)
    Implicit None
    Integer, intent(in) :: ie, je, ke
!    Type(element_t), intent(in) :: vars
    Real(Kind=real_kind), dimension(nx,nx,nx) :: wt

    Real(Kind=real_kind) :: sg,sg13,sg23
    Integer :: i,j,k

    do k = 1, nx
      do j = 1, nx
        do i = 1, nx
          sg13 = metric(ie,je,ke)%sg13(i,j,k)
          sg23 = metric(ie,je,ke)%sg23(i,j,k)
          sg = metric(ie,je,ke)%sg(i,j)

          ! Based on G-S transform
          wt(i,j,k) = (vars%w(i,j,k) + sg13 * vars%u(i,j,k) + sg23 * vars%v(i,j,k))/sg
        end do
      end do
    end do
  End Function Compute_Zeta_Velocity


  !!=========================================================================================
  !! Output mountain height and mapped vertical coordinates.
  Subroutine plot_mountain()
    Use mpi_mod, only: rank_is_master
    Use grid_setup_mod, only: get_visualization_grid

    Implicit none
    Real(Kind=real_kind), dimension(nux) :: viz_x
    Real(Kind=real_kind), dimension(nuy) :: viz_y
    Real(Kind=real_kind), dimension(nuz) :: viz_z
    Real(Kind=real_kind), dimension(nux,nuy) :: mtn
    Real(Kind=real_kind), dimension(nux,nuz) :: vslice
    Real(Kind=real_kind) :: x,y,zeta, hxy, dh_x, dh_y, ght, sg, sg13, sg23
    integer :: iu,ju,ku

    ! only write data on master rank!
    if (.not. rank_is_master) return

    call get_visualization_grid(viz_x, viz_y, viz_z)

    do ju = 1, nuy
      do iu = 1, nux
        x = viz_x(iu)
        y = viz_y(ju)
        select case(testcase)
        case (testcase_id_mountain_agnesi)
          Call Agnesi_Mountain_Surfdat(x,y,hxy,dh_x,dh_y)
        case(testcase_id_mountain_schar)
          Call Schar_Mountain_Surfdat(x,y,hxy,dh_x,dh_y)
        case default
          print*, "unrecognized mountain type"
          stop
        end select
        mtn(iu,ju) = hxy
      end do
    end do

    ! Display data for surface contouring
    3 format(e24.16)
    open (unit = 20, file = output_dir//'/mtn_surface.dat')
    write(20,3) ((mtn(iu,ju),iu=1,nux),ju=1,nuy)
    close(20)

    ju = int((nuy+1)/2.0D0)

    y = viz_y(ju)
    do ku = 1, nuz
      do iu = 1, nux
        x = viz_x(iu)
        zeta = viz_z(ku)
        call Compute_metric_at_coords(x, y, zeta, hxy, ght, sg, sg13, sg23)
        vslice(iu,ku) = ght
      end do
    end do

    open (unit = 22, file = output_dir//'/mtn_grid_xz.dat')
    write(22,3) ((vslice(iu,ku),iu=1,nux),ku=1,nuz)
    close(22)

  End Subroutine plot_mountain


  !!=========================================================================================
  Subroutine Agnesi_Mountain_Surfdat(x,y,hxy,dh_x,dh_y)
    Implicit None
    Real(Kind=real_kind), intent(in) :: x,y
    Real(Kind=real_kind), intent(out) :: hxy, dh_x, dh_y

    Real(Kind=real_kind) :: h0,a0,b0, x0,y0, deno, deno52

    h0 = 3000.0D0
    a0 = 12000.0D0
    b0 = 4000.0D0
    x0 = xmax/2.0D0
    y0 = ymax/2.0D0

    deno =  ((x-x0)/a0)**2 + ((y-y0)/b0)**2 + 1.0D0

    !mountain
    hxy = h0/(deno)**1.50D0

    deno52  =  1.0D0 / (deno)**2.5D0
    !derivative w.r.t x, y
    dh_x = -3.0D0 * (h0 /a0**2) * (x-x0) * deno52
    dh_y = -3.0D0 * (h0 /b0**2) * (y-y0) * deno52

  End Subroutine Agnesi_Mountain_Surfdat


  !!=========================================================================================
  !! Schar mountain profile h_s(x)
  Subroutine Schar_Mountain_Surfdat(x,y,mt,dh_x,dh_y)
    Implicit None
    Real(Kind=real_kind), intent(in) :: x,y
    Real(Kind=real_kind), intent(out) :: mt, dh_x, dh_y

    Real(Kind=real_kind) :: h0, a0,b0, x0,y0, r, dterm
    Real(Kind=real_kind), parameter :: real_eps = 1.0E-16

    h0 = 5000.0D0
    a0 = 10000.0D0
    b0 = 6000.0D0
    x0 = xmax/2.0D0
    y0 = ymax/2.0D0
    r = sqrt((x-x0)**2 + (y-y0)**2)
    mt = h0 * exp(-(r/a0)**2) * (cos(pi*r/b0))**2

    dterm = 2.0D0 * ( r/a0**2 + tan(pi*r/b0)*pi/b0 )

    if (r <= real_eps) then
      dh_x = 0.0D0
      dh_y = 0.0D0
    else
      dh_x = -mt * dterm * (x-x0)/r
      dh_y = -mt * dterm * (y-y0)/r
    endif

  End Subroutine Schar_Mountain_Surfdat

end module mountain_grid_mod

