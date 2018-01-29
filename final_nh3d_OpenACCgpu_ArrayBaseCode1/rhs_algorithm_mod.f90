!-------------------------------------------------------------------
! R.Nair NCAR/scd 10/02
! Rewritten 06/2016 by Francois Hebert, SIParCS/Cornell
!
! Implements the weak form DG method for evaluating the
! semi-discrete RHS of the advection equation.
!-------------------------------------------------------------------

module rhs_algorithm_mod
  use basic_mod
  use element_mod
  use gauss_quadrature_mod

  implicit none
  private
  public :: compute_rhs

contains

  ! Computes the RHS of `dU/dt = rhs'
  ! * Implements a weak form DG algorithm
  ! * Outsources the work of computing the RHS equations of the specific
  !   physical system to a subprogram.
  subroutine compute_rhs(ie,je,ke, rhs)

!    type(element_t), intent(in) :: vars
!    type(edge_t), intent(inout), dimension(6) :: edgevars
    real(kind=real_kind), intent(out), dimension(nx,nx,nx) :: rhs

    real(kind=real_kind), dimension(nx,nx,nx) :: flx, fly, flz, divf
    real(kind=real_kind), dimension(6,nx,nx) :: numflux

    Call  Extract_normalFlux(ie,je,ke,flx,fly,flz,numflux)

    !call advection_fluxes(vars, edgevars, flx, fly, flz)
    !call lf_numflux(edgevars, numflux)

    call flux_divergence(flx,fly,flz,divf)

    ! RHS = divF + \sum bdryF
    ! start with volume divF term, then add each bdry term one-by-one
    rhs(:,:,:) = divf(:,:,:)

    ! "lower" edges -- west, south, bottom
    ! NOTE: The numflux definition used includes the normal vector
    !       (because it's computed from 'outgoing' speed among u,v,w),
    !       so the expected (+) sign picks up a (-) on lower edges.
    rhs(1,:,:) = rhs(1,:,:) - numflux(4,:,:) * 2.0D0/(delx*gllw(1))
    rhs(:,1,:) = rhs(:,1,:) - numflux(1,:,:) * 2.0D0/(dely*gllw(1))
    rhs(:,:,1) = rhs(:,:,1) - numflux(5,:,:) * 2.0D0/(delz*gllw(1))

    ! "upper" edges -- east, north, top
    rhs(nx,:,:) = rhs(nx,:,:) - numflux(2,:,:) * 2.0D0/(delx*gllw(nx))
    rhs(:,nx,:) = rhs(:,nx,:) - numflux(3,:,:) * 2.0D0/(dely*gllw(nx))
    rhs(:,:,nx) = rhs(:,:,nx) - numflux(6,:,:) * 2.0D0/(delz*gllw(nx))

  end subroutine compute_rhs


  ! Compute the RHS terms for the advection equation.
  ! * Computes the flux in the volume and on the surface as F = psi*velocity
  ! * Currently no sources are used
  subroutine advection_fluxes(ie,je,ke, flx, fly, flz)

!    type(element_t), intent(in) :: vars
!    type(edge_t), intent(inout), dimension(6) :: edgevars
    real(kind=real_kind), intent(out), dimension(nx,nx,nx) :: flx, fly, flz
    integer :: e

    ! compute volume flux components
    ! NOTE: This code will need to be adjusted for deformed domains:
    !       it will need to include metric terms and couple u,v,w.
    flx(:,:,:) = vars_psi(:,:,:,ie,je,ke) * vars_u(:,:,:,ie,je,ke)
    fly(:,:,:) = vars_psi(:,:,:,ie,je,ke) * vars_v(:,:,:,ie,je,ke)
    flz(:,:,:) = vars_psi(:,:,:,ie,je,ke) * vars_wt(:,:,:,ie,je,ke) ! on zeta-grid, w -> wt

    ! compute surface flux components
    ! NOTE: see above, but changes to be made in extract_edgevars

    do e = 1, 6
      edgevars(e)%psiflux_int(:,:) = edgevars(e)%psi_int(:,:) * edgevars(e)%normalspd_int(:,:)
      edgevars(e)%psiflux_ext(:,:) =-edgevars(e)%psi_ext(:,:) * edgevars(e)%normalspd_ext(:,:)
    end do

  ! Call Extract_normalFlux(vars, edgevars) 

  end subroutine advection_fluxes

!-------------
   Subroutine Extract_normalFlux(vars, edgevars,flx_x,flx_y,flx_z,num_flux)
    Integer, parameter :: isouth=1,inorth=3,iwest=4,ieast=2,ibot=5,itop=6 
    type(element_t), intent(in) :: vars
    type(edge_t), intent(inout), dimension(6) :: edgevars
    real(kind=real_kind), intent(out), dimension(nx,nx,nx) :: flx_x, flx_y, flx_z
    real(kind=real_kind), intent(out), dimension(6,nx,nx) :: num_flux
    real(kind=real_kind), dimension(6,nx,nx) :: nvelo_int,nvelo_ext
    real(kind=real_kind), dimension(nx,nx) :: psi_int,psi_ext,flx_int,flx_ext,spd_int,spd_ext 
    real(kind=real_kind) :: maxspd
    integer :: e

  !Internal fluxes 
    flx_x(:,:,:) = vars%psi(:,:,:) * vars%u(:,:,:)
    flx_y(:,:,:) = vars%psi(:,:,:) * vars%v(:,:,:)
    flx_z(:,:,:) = vars%psi(:,:,:) * vars%wt(:,:,:) ! on zeta-grid, w -> wt

 ! Normal component of velocity for flux integrals 
      nvelo_int(isouth,:,:) = - edgevars(1)%v_int(:,:) 
      nvelo_int(ieast ,:,:) =   edgevars(2)%u_int(:,:) 
      nvelo_int(inorth,:,:) =   edgevars(3)%v_int(:,:) 
      nvelo_int(iwest ,:,:) = - edgevars(4)%u_int(:,:) 
      nvelo_int(ibot  ,:,:) = - edgevars(5)%w_int(:,:) 
      nvelo_int(itop  ,:,:) =   edgevars(6)%w_int(:,:) 

      nvelo_ext(isouth,:,:) = - edgevars(1)%v_ext(:,:) 
      nvelo_ext(ieast ,:,:) =   edgevars(2)%u_ext(:,:) 
      nvelo_ext(inorth,:,:) =   edgevars(3)%v_ext(:,:) 
      nvelo_ext(iwest ,:,:) = - edgevars(4)%u_ext(:,:) 
      nvelo_ext(ibot  ,:,:) = - edgevars(5)%w_ext(:,:) 
      nvelo_ext(itop  ,:,:) =   edgevars(6)%w_ext(:,:) 


 ! Lax-Fred Numerical fluxes 

    do e = 1, 6
!     edgevars(e)%psiflux_int(:,:) = edgevars(e)%psi_int(:,:) * nvelo_int(e,:,:)
!     edgevars(e)%psiflux_ext(:,:) = edgevars(e)%psi_ext(:,:) * nvelo_ext(e,:,:)
      psi_int(:,:) = edgevars(e)%psi_int(:,:) 
      psi_ext(:,:) = edgevars(e)%psi_ext(:,:) 
      spd_int(:,:) = nvelo_int(e,:,:) 
      spd_ext(:,:) = nvelo_ext(e,:,:) 
      flx_int(:,:) = psi_int(:,:)*spd_int(:,:) 
      flx_ext(:,:) = psi_ext(:,:)*spd_ext(:,:) 

        maxspd = max(maxval(abs(spd_int)), maxval(abs(spd_ext)))
        num_flux(e,:,:) = (flx_ext(:,:) + flx_int(:,:) - maxspd * (psi_ext(:,:) - psi_int(:,:))) / 2.0D0
    end do

  end subroutine Extract_normalFlux


!-------------

  ! Compute the flux divergence according to the weak form of DG
  subroutine flux_divergence(flx,fly,flz,divf)

    real(kind=real_kind), intent(in), dimension(nx,nx,nx) :: flx, fly, flz
    real(kind=real_kind), intent(out), dimension(nx,nx,nx) :: divf
    real(kind=real_kind) :: dflx, dfly, dflz
    integer :: i,j,k,a

    ! Use pre-computed `weak form deriv matrix' weakder_ij = der_ij w_i / w_j
    ! to compute the derivative in each dimension.
    do k = 1, nx
      do j = 1, nx
        do i = 1, nx

          dflx = 0.0D0
          do a = 1, nx
            dflx = dflx + flx(a,j,k) * weakder(a,i)
          end do

          dfly = 0.0D0
          do a = 1, nx
            dfly = dfly + fly(i,a,k) * weakder(a,j)
          end do

          dflz = 0.0D0
          do a = 1, nx
            dflz = dflz + flz(i,j,a) * weakder(a,k)
          end do

          ! add the 2/dx geometric factors
          divf(i,j,k) = dflx * 2/delx + dfly * 2/dely + dflz * 2/delz
        end do
      end do
    end do

  end subroutine flux_divergence


  ! Compute the Lax-Friedrichs numerical flux at an interface.
  ! * The maximum speed is computed as a max over the entire surface,
  !   and is therefore more dissipative than a pointwise maximum.
  ! * The `normalspd' variable holds the flow velocity normal to the
  !   surface, and is (+) for outgoing flow and (-) for ingoing flow.
  subroutine lf_numflux(edgevars, numflux)

    type(edge_t), intent(in), dimension(6) :: edgevars
    real(kind=real_kind), intent(out), dimension(6,nx,nx) :: numflux
    real(kind=real_kind) :: maxspd
    integer :: e

    do e = 1, 6
      associate(s_int => edgevars(e)%normalspd_int, &
                s_ext => edgevars(e)%normalspd_ext, &
                f_int => edgevars(e)%psiflux_int, &
                f_ext => edgevars(e)%psiflux_ext, &
                u_int => edgevars(e)%psi_int, &
                u_ext => edgevars(e)%psi_ext)
        maxspd = max(maxval(abs(s_int)), maxval(abs(s_ext)))
        numflux(e,:,:) = (f_ext(:,:) + f_int(:,:) - maxspd * (u_ext(:,:) - u_int(:,:))) / 2.0D0
      end associate
    end do

  end subroutine lf_numflux

end module rhs_algorithm_mod
