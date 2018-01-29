!-------------------------------------------------------------------
! 06/2016 Francois Hebert, SIParCS/Cornell
!
! Fills the interior/exterior edge data for each element.
!
! Interior data is directly copied from the volume data structure,
! because we use GLL points. Exterior data is communicated from the
! neighboring element, using MPI if necessary.
!
! Dictionary for the indices/directions/words:
! index 1 = -y = south
! index 2 = +x = east
! index 3 = +y = north
! index 4 = -x = west
! index 5 = -z = bottom
! index 6 = +z = top
!-------------------------------------------------------------------

module prepare_edgevars_mod
  use basic_mod
  use element_mod
  use mpi_mod, only: send_data, recv_data
  use openacc

  implicit none
  private
  public :: extract_edgevars, communicate_edgevars_work, pack_edgevars, unpack_edgevars

  ! allocatable because buffer length depends on the number of MPI
  ! partitions, which is determined at runtime.
!  real(kind=real_kind), dimension(:), allocatable :: bufsendS, bufrecvS
!  real(kind=real_kind), dimension(:), allocatable :: bufsendE, bufrecvE
!  real(kind=real_kind), dimension(:), allocatable :: bufsendN, bufrecvN
!  real(kind=real_kind), dimension(:), allocatable :: bufsendW, bufrecvW
!
!  integer, parameter :: number_of_send_recv_vars = 5

contains

!  ! Simple wrapper for the extraction and communication subroutines
!  subroutine prepare_edgevars()
!!    type(element_t), intent(inout), dimension(:,:,:) :: vars
!!    type(surface_t), intent(out), dimension(:,:,:) :: edgevars
!    call extract_edgevars()
!    call communicate_edgevars()
!
!  end subroutine prepare_edgevars

!!=========================================================================================
  ! Fill the interior edge data by extracting (i.e. copying) from the
  ! interior volume data. Quantities extracted are state vectors  (neq = 5)        

  subroutine extract_edgevars()
!    type(element_t), intent(in), dimension(:,:,:) :: vars
!    type(surface_t), intent(out), dimension(:,:,:) :: edgevars
    integer :: nlx, nly, nlz
    integer :: n, ii, jj, kk, i, j

    nlx = size(edgevars_edges_rho_int,4)
    nly = size(edgevars_edges_rho_int,5)
    nlz = size(edgevars_edges_rho_int,6)

!$acc parallel num_workers(5) vector_length(16)
!$acc loop gang worker vector collapse(6)
    do kk = 1, nlz
      do jj = 1, nly
        do ii = 1, nlx

         do n = 1, neq 
        do j = 1, nx
        do i = 1, nx
          edgevars_edges_vec_int(i,j,n,1,ii,jj,kk) = vars_vec(i,1,j, n,ii,jj,kk)
          edgevars_edges_vec_int(i,j,n,2,ii,jj,kk) = vars_vec(nx,i,j,n,ii,jj,kk)
          edgevars_edges_vec_int(i,j,n,3,ii,jj,kk) = vars_vec(i,nx,j,n,ii,jj,kk)
          edgevars_edges_vec_int(i,j,n,4,ii,jj,kk) = vars_vec(1,i,j, n,ii,jj,kk)
          edgevars_edges_vec_int(i,j,n,5,ii,jj,kk) = vars_vec(i,j,1, n,ii,jj,kk)
          edgevars_edges_vec_int(i,j,n,6,ii,jj,kk) = vars_vec(i,j,nx,n,ii,jj,kk)
         end do
        end do
        end do

        end do
      end do
    end do
!$acc end parallel

    if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)

  end subroutine extract_edgevars

!!!=========================================================================================
!  ! Fill exterior edge data by getting the neighbor's interior values.
!  ! If the neighboring element is on the same MPI rank, then this is a
!  ! simple copy. Otherwise, the data is sent/received using MPI.
!  subroutine communicate_edgevars()
!!    type(surface_t), intent(out), dimension(:,:,:) :: edgevars
!    integer :: nlx, nly, nlz
!    integer :: ii, jj, kk, ip, jp, kp, im, jm, km, n 
!
!    nlx = size(edgevars_edges_rho_int,4)
!    nly = size(edgevars_edges_rho_int,5)
!    nlz = size(edgevars_edges_rho_int,6)
!
!    if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)
!
!    ! Pack the interior edge data for those edges that interface with a
!    ! different MPI rank, and send the data off asynchronously.
!    call pack_edgevars(bufsendS, bufsendE, bufsendN, bufsendW)
!    call send_data(bufsendS, bufsendE, bufsendN, bufsendW)
!
!    ! This loops over every element on this MPI rank, and gets the external
!    ! edge vars from the neighboring elements. Periodicity is applied, so that
!    ! only local data is manipulated.
!    !  The external boundaries of this rank's elements will be overwritten
!    !       by the received MPI data => room for small optimization here.
!
!    do kk = 1, nlz
!      do jj = 1, nly
!        do ii = 1, nlx
!
!          ip = ii + 1
!          im = ii - 1
!          jp = jj + 1
!          jm = jj - 1
!          kp = kk + 1
!          km = kk - 1
!
!          ! apply periodicity in each dimension
!          if (ii == 1) im = nlx
!          if (ii == nlx) ip = 1
!          if (jj == 1) jm = nly
!          if (jj == nly) jp = 1
!          if (kk == 1) km =  kk !nlz
!          if (kk == nlz) kp = kk  !1
!
!         do n = 1, neq 
!          edgevars_edges_vec_ext(:,:,n,1,ii,jj,kk) = edgevars_edges_vec_int(:,:,n,3,ii,jm,kk)
!          edgevars_edges_vec_ext(:,:,n,2,ii,jj,kk) = edgevars_edges_vec_int(:,:,n,4,ip,jj,kk)
!          edgevars_edges_vec_ext(:,:,n,3,ii,jj,kk) = edgevars_edges_vec_int(:,:,n,1,ii,jp,kk)
!          edgevars_edges_vec_ext(:,:,n,4,ii,jj,kk) = edgevars_edges_vec_int(:,:,n,2,im,jj,kk)
!          edgevars_edges_vec_ext(:,:,n,5,ii,jj,kk) = edgevars_edges_vec_int(:,:,n,6,ii,jj,km)
!          edgevars_edges_vec_ext(:,:,n,6,ii,jj,kk) = edgevars_edges_vec_int(:,:,n,5,ii,jj,kp)
!         end do
!
!        end do
!      end do
!    end do
!!
!!   ! for no-flux BC, simply overwrite the periodic BC data with no-flux data
!!   if (use_noflux_bc_top_bottom) then
!
!!    do jj = 1, nly
!!      do ii = 1, nlx
!!        edgevars(ii,jj,1)%edges(5)%rho_ext(:,:) = edgevars(ii,jj,1)%edges(5)%rho_int(:,:)
!!        edgevars(ii,jj,nez)%edges(6)%rho_ext(:,:) = edgevars(ii,jj,nez)%edges(6)%rho_int(:,:)
!
!!        edgevars(ii,jj,1)%edges(5)%u_ext(:,:) = edgevars(ii,jj,1)%edges(5)%u_int(:,:)
!!        edgevars(ii,jj,nez)%edges(6)%u_ext(:,:) = edgevars(ii,jj,nez)%edges(6)%u_int(:,:)
!
!!        edgevars(ii,jj,1)%edges(5)%v_ext(:,:) = edgevars(ii,jj,1)%edges(5)%v_int(:,:)
!!        edgevars(ii,jj,nez)%edges(6)%v_ext(:,:) = edgevars(ii,jj,nez)%edges(6)%v_int(:,:)
!
!!        edgevars(ii,jj,1)%edges(5)%the_ext(:,:) = edgevars(ii,jj,1)%edges(5)%the_int(:,:)
!!        edgevars(ii,jj,nez)%edges(6)%the_ext(:,:) = edgevars(ii,jj,nez)%edges(6)%the_int(:,:)
!
!!        edgevars(ii,jj,1)%edges(5)%w_ext(:,:) =   edgevars(ii,jj,1)%edges(5)%w_int(:,:)
!!        edgevars(ii,jj,nez)%edges(6)%w_ext(:,:) =   edgevars(ii,jj,nez)%edges(6)%w_int(:,:) 
!!      end do
!!    end do
!
!  ! end if
!
!    ! Receive exterior data from edges that interface with a different MPI
!    ! rank, and unpack it into edge_vars.
!    call recv_data(bufrecvS, bufrecvE, bufrecvN, bufrecvW)
!    call unpack_edgevars(bufrecvS, bufrecvE, bufrecvN, bufrecvW)
!  end subroutine communicate_edgevars
!
!
!!=========================================================================================

subroutine communicate_edgevars_work()

    integer :: nlx, nly, nlz
    integer :: ii, jj, kk, ip, jp, kp, im, jm, km, n, i, j

    nlx = size(edgevars_edges_rho_int,4)
    nly = size(edgevars_edges_rho_int,5)
    nlz = size(edgevars_edges_rho_int,6)

!$acc parallel num_workers(5) vector_length(16)
!$acc loop gang collapse(3) private(ip, im, jp, jm, kp, km)
do kk = 1, nlz
      do jj = 1, nly
        do ii = 1, nlx

          ip = ii + 1
          im = ii - 1
          jp = jj + 1
          jm = jj - 1
          kp = kk + 1
          km = kk - 1

          ! apply periodicity in each dimension
          if (ii == 1) im = nlx
          if (ii == nlx) ip = 1
          if (jj == 1) jm = nly
          if (jj == nly) jp = 1
          if (kk == 1) km =  kk !nlz
          if (kk == nlz) kp = kk  !1

!$acc loop worker vector collapse(3)
         do n = 1, neq 
        do j = 1, nx
        do i = 1, nx
          edgevars_edges_vec_ext(i,j,n,1,ii,jj,kk) = edgevars_edges_vec_int(i,j,n,3,ii,jm,kk)
          edgevars_edges_vec_ext(i,j,n,2,ii,jj,kk) = edgevars_edges_vec_int(i,j,n,4,ip,jj,kk)
          edgevars_edges_vec_ext(i,j,n,3,ii,jj,kk) = edgevars_edges_vec_int(i,j,n,1,ii,jp,kk)
          edgevars_edges_vec_ext(i,j,n,4,ii,jj,kk) = edgevars_edges_vec_int(i,j,n,2,im,jj,kk)
          edgevars_edges_vec_ext(i,j,n,5,ii,jj,kk) = edgevars_edges_vec_int(i,j,n,6,ii,jj,km)
          edgevars_edges_vec_ext(i,j,n,6,ii,jj,kk) = edgevars_edges_vec_int(i,j,n,5,ii,jj,kp)
         end do
        end do
        end do

        end do
      end do
    end do
!$acc end parallel

end subroutine communicate_edgevars_work

!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

  subroutine pack_edgevars()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_edges_rho_int,4)
    nly = size(edgevars_edges_rho_int,5)
    nlz = size(edgevars_edges_rho_int,6)

!    idx = 1
!$acc parallel vector_length(16)
!$acc loop gang vector collapse(4) private(idx)
    do ke = 1, nlz
      do ie = 1, nlx
        do k = 1, nx
          do i = 1, nx
idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            bufsendS(idx) = edgevars_edges_vec_int(i,k,  1,1,ie,1,ke)
            bufsendN(idx) = edgevars_edges_vec_int(i,k,1,3,ie,nly,ke)
            idx = idx + 1

            bufsendS(idx) = edgevars_edges_vec_int(i,k,  2,1,ie,1,ke)
            bufsendN(idx) = edgevars_edges_vec_int(i,k,2,3,ie,nly,ke)
            idx = idx + 1

            bufsendS(idx) = edgevars_edges_vec_int(i,k,  3,1,ie,1,ke)
            bufsendN(idx) = edgevars_edges_vec_int(i,k,3,3,ie,nly,ke)
            idx = idx + 1

            bufsendS(idx) = edgevars_edges_vec_int(i,k,  4,1,ie,1,ke)
            bufsendN(idx) = edgevars_edges_vec_int(i,k,4,3,ie,nly,ke)
            idx = idx + 1

            bufsendS(idx) = edgevars_edges_vec_int(i,k,  5,1,ie,1,ke)
            bufsendN(idx) = edgevars_edges_vec_int(i,k,5,3,ie,nly,ke)
!            idx = idx + 1

          end do
        end do
      end do
    end do
!$acc end parallel

!    idx = 1
!$acc parallel vector_length(16)
!$acc loop gang vector collapse(4) private(idx)
    do ke = 1, nlz
      do je = 1, nly
        do k = 1, nx
          do j = 1, nx
idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            bufsendE(idx) = edgevars_edges_vec_int(j,k,1,2,nlx,je,ke)
            bufsendW(idx) = edgevars_edges_vec_int(j,k,  1,4,1,je,ke)
            idx = idx + 1

            bufsendE(idx) = edgevars_edges_vec_int(j,k,2,2,nlx,je,ke)
            bufsendW(idx) = edgevars_edges_vec_int(j,k,  2,4,1,je,ke)
            idx = idx + 1

            bufsendE(idx) = edgevars_edges_vec_int(j,k,3,2,nlx,je,ke)
            bufsendW(idx) = edgevars_edges_vec_int(j,k,  3,4,1,je,ke)
            idx = idx + 1

            bufsendE(idx) = edgevars_edges_vec_int(j,k,4,2,nlx,je,ke)
            bufsendW(idx) = edgevars_edges_vec_int(j,k,  4,4,1,je,ke)
            idx = idx + 1

            bufsendE(idx) = edgevars_edges_vec_int(j,k,5,2,nlx,je,ke)
            bufsendW(idx) = edgevars_edges_vec_int(j,k,  5,4,1,je,ke)
!            idx = idx + 1
          end do
        end do
      end do
    end do
!$acc end parallel

  end subroutine pack_edgevars


!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
  subroutine unpack_edgevars()
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_edges_rho_int,4)
    nly = size(edgevars_edges_rho_int,5)
    nlz = size(edgevars_edges_rho_int,6)

!    idx = 1
!$acc parallel vector_length(16)
!$acc loop gang vector collapse(4) private(idx)
    do ke = 1, nlz
      do ie = 1, nlx
        do k = 1, nx
          do i = 1, nx
idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            edgevars_edges_vec_ext(i,k,  1,1,ie,1,ke) = bufrecvS(idx)
            edgevars_edges_vec_ext(i,k,1,3,ie,nly,ke) = bufrecvN(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(i,k,  2,1,ie,1,ke) = bufrecvS(idx)
            edgevars_edges_vec_ext(i,k,2,3,ie,nly,ke) = bufrecvN(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(i,k,  3,1,ie,1,ke) = bufrecvS(idx)
            edgevars_edges_vec_ext(i,k,3,3,ie,nly,ke) = bufrecvN(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(i,k,  4,1,ie,1,ke) = bufrecvS(idx)
            edgevars_edges_vec_ext(i,k,4,3,ie,nly,ke) = bufrecvN(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(i,k,  5,1,ie,1,ke) = bufrecvS(idx)
            edgevars_edges_vec_ext(i,k,5,3,ie,nly,ke) = bufrecvN(idx)
!            idx = idx + 1

          end do
        end do
      end do
    end do
!$acc end parallel

!    idx = 1
!$acc parallel vector_length(16)
!$acc loop gang vector collapse(4) private(idx)
    do ke = 1, nlz
      do je = 1, nly
        do k = 1, nx
          do j = 1, nx
idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            edgevars_edges_vec_ext(j,k,1,2,nlx,je,ke) = bufrecvE(idx)
            edgevars_edges_vec_ext(j,k,  1,4,1,je,ke) = bufrecvW(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(j,k,2,2,nlx,je,ke) = bufrecvE(idx)
            edgevars_edges_vec_ext(j,k,  2,4,1,je,ke) = bufrecvW(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(j,k,3,2,nlx,je,ke) = bufrecvE(idx)
            edgevars_edges_vec_ext(j,k,  3,4,1,je,ke) = bufrecvW(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(j,k,4,2,nlx,je,ke) = bufrecvE(idx)
            edgevars_edges_vec_ext(j,k,  4,4,1,je,ke) = bufrecvW(idx)
            idx = idx + 1

            edgevars_edges_vec_ext(j,k,5,2,nlx,je,ke) = bufrecvE(idx)
            edgevars_edges_vec_ext(j,k,  5,4,1,je,ke) = bufrecvW(idx)
!            idx = idx + 1

          end do
        end do
      end do
    end do
!$acc end parallel

  end subroutine unpack_edgevars


!!=========================================================================================
  ! Allocate the send/recv buffers used in MPI communication
  subroutine allocate_buffers(nlx, nly, nlz)
    integer, intent(in) :: nlx, nly, nlz
    integer :: nSN, nEW
    nSN = number_of_send_recv_vars * nx*nx * nlx*nlz
    nEW = number_of_send_recv_vars * nx*nx * nly*nlz
    if (.not. allocated(bufsendS)) allocate(bufsendS(nSN))
    if (.not. allocated(bufrecvS)) allocate(bufrecvS(nSN))
    if (.not. allocated(bufsendE)) allocate(bufsendE(nEW))
    if (.not. allocated(bufrecvE)) allocate(bufrecvE(nEW))
    if (.not. allocated(bufsendN)) allocate(bufsendN(nSN))
    if (.not. allocated(bufrecvN)) allocate(bufrecvN(nSN))
    if (.not. allocated(bufsendW)) allocate(bufsendW(nEW))
    if (.not. allocated(bufrecvW)) allocate(bufrecvW(nEW))
  end subroutine allocate_buffers

end module prepare_edgevars_mod
