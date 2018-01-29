!-------------------------------------------------------------------
! RD Nair (IMAGe/NCAR)  June 2010
! Changes 06/2016 by Francois Hebert, SIParCS/Cornell
!
! Defines an advection test case:
! Uniform advection of a cosine blob over a mountain.
!
! This problem has physical dimensions!
! * domain is 120km x 30km x 30km in size
! * the flow rate is 20m/s in the +x direction
!   domain is periodic => the cosine blob cycles through the box
!-------------------------------------------------------------------

module testcase_mountain_mod
  use basic_mod
  use element_mod
  use testcase_helper_mod
  use mountain_grid_mod
  use mpi_mod, only: rank_is_master

  implicit none
  private
  public :: tc_mountain_load_params
  public :: tc_mountain_initial_data
  public :: tc_mountain_update_velocities

  real(kind=real_kind), parameter :: tc_xmin = 0.0D0
  real(kind=real_kind), parameter :: tc_xmax = 120000.0D0 ! km
  real(kind=real_kind), parameter :: tc_ymin = 0.0D0
  real(kind=real_kind), parameter :: tc_ymax = 30000.0D0 ! km
  real(kind=real_kind), parameter :: tc_zmin = 0.0D0
  real(kind=real_kind), parameter :: tc_zmax = 30000.0D0 ! km

  real(kind=real_kind), parameter :: tc_endtime = 3000.0D0 ! seconds

contains

  subroutine tc_mountain_load_params()
    xmin = tc_xmin
    xmax = tc_xmax
    ymin = tc_ymin
    ymax = tc_ymax
    zmin = tc_zmin
    zmax = tc_zmax
    endtime = tc_endtime
  end subroutine tc_mountain_load_params


  subroutine tc_mountain_initial_data()

!    type(grid_t), intent(in), dimension(:,:,:) :: grid
!    type(metric_t), intent(in), dimension(:,:,:) :: metric
!    type(element_t), intent(out), dimension(:,:,:) :: vars

    real(kind=real_kind), dimension(nx,nx,nx) :: psi
    integer :: nlx, nly, nlz, ie, je, ke, k

    real(kind=real_kind), parameter :: x0 = tc_xmax/6.0D0
    real(kind=real_kind), parameter :: y0 = tc_ymax/2.0D0
    real(kind=real_kind), parameter :: z0 = tc_zmax/2.0D0
    real(kind=real_kind), parameter :: a0 = 16000.0D0
    real(kind=real_kind), parameter :: b0 = 8000.0D0
    real(kind=real_kind), parameter :: c0 = 8000.0D0
    real(kind=real_kind), parameter :: amp = 1.0D0

    if (size(grid) /= size(vars)) stop
    nlx = size(grid,1)
    nly = size(grid,2)
    nlz = size(grid,3)

    call compute_metric_terms(grid)

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx

          ! physical-frame velocity
          vars(ie,je,ke)%u = 20.0D0 ! m/s
          vars(ie,je,ke)%v = 0.0D0
          vars(ie,je,ke)%w = 0.0D0

          ! computational-frame velocity
          ! w => wt = d(zeta)/dt in computational grid
          vars(ie,je,ke)%wt = compute_zeta_velocity(ie,je,ke)

          ! We evolve the jacobian-weighted scalar in the computational frame
          psi = cosine_blob(ie,je,ke, x0, y0, z0, a0, b0, c0, amp)
          do k = 1, nx
            vars(ie,je,ke)%psi(:,:,k) =  psi(:,:,k) * metric(ie,je,ke)%sg(:,:)
          enddo

        end do
      end do
    end do

    if (rank_is_master) call plot_mountain()

  end subroutine tc_mountain_initial_data


  subroutine tc_mountain_update_velocities(vars_ref, vars)
    type(element_t), intent(in), dimension(:,:,:) :: vars_ref
    type(element_t), intent(inout), dimension(:,:,:) :: vars
    call copy_velocities(vars_ref, vars)
  end subroutine tc_mountain_update_velocities


end module testcase_mountain_mod
