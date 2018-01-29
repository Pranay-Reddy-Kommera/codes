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
  use cudafor

  implicit none
  private
  public :: extract_edgevars, communicate_edgevars_work, pack_edgevars, unpack_edgevars,allocate_buffers,pack_edgevars1,pack_edgevars2,unpack_edgevars1,unpack_edgevars2

!  ! allocatable because buffer length depends on the number of MPI
!  ! partitions, which is determined at runtime.
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

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

    do kk = 1, nlz
      do jj = 1, nly
        do ii = 1, nlx

         do n = 1, neq 
        do j = 1, nx
        do i = 1, nx
          edgevars(ii,jj,kk)%edges(1)%vec_int(i,j,n) = vars(ii,jj,kk)%vec(i,1,j, n)
          edgevars(ii,jj,kk)%edges(2)%vec_int(i,j,n) = vars(ii,jj,kk)%vec(nx,i,j,n)
          edgevars(ii,jj,kk)%edges(3)%vec_int(i,j,n) = vars(ii,jj,kk)%vec(i,nx,j,n)
          edgevars(ii,jj,kk)%edges(4)%vec_int(i,j,n) = vars(ii,jj,kk)%vec(1,i,j, n)
          edgevars(ii,jj,kk)%edges(5)%vec_int(i,j,n) = vars(ii,jj,kk)%vec(i,j,1, n)
          edgevars(ii,jj,kk)%edges(6)%vec_int(i,j,n) = vars(ii,jj,kk)%vec(i,j,nx,n)
        end do
        end do
         end do

        end do
      end do
    end do

     if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)

  end subroutine extract_edgevars

!!=========================================================================================
  ! Fill the interior edge data by extracting (i.e. copying) from the
  ! interior volume data. Quantities extracted are state vectors  (neq = 5)        

 attributes(global) subroutine extract_edges()
!    type(element_t), intent(in), dimension(:,:,:) :: vars
!    type(surface_t), intent(out), dimension(:,:,:) :: edgevars
!    integer :: nlx, nly, nlz
    integer :: n, ii, jj, kk, i, j

!    nlx = size(edgevars,1)
!    nly = size(edgevars,2)
!    nlz = size(edgevars,3)

ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
i = threadIdx%x
j = threadIdx%y
n = threadIdx%z

!    do kk = 1, nlz
!      do jj = 1, nly
!        do ii = 1, nlx

!         do n = 1, neq 
!        do j = 1, nx
!        do i = 1, nx
          edgevars_d(ii,jj,kk)%edges(1)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,1,j, n)
          edgevars_d(ii,jj,kk)%edges(2)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(nx,i,j,n)
          edgevars_d(ii,jj,kk)%edges(3)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,nx,j,n)
          edgevars_d(ii,jj,kk)%edges(4)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(1,i,j, n)
          edgevars_d(ii,jj,kk)%edges(5)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,j,1, n)
          edgevars_d(ii,jj,kk)%edges(6)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,j,nx,n)
!        end do
!        end do
!         end do
!
!        end do
!      end do
!    end do

!     if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)

  end subroutine extract_edges
!!=========================================================================================
  ! Fill the interior edge data by extracting (i.e. copying) from the
  ! interior volume data. Quantities extracted are state vectors  (neq = 5)        

 attributes(global) subroutine extract_edges1()
!    type(element_t), intent(in), dimension(:,:,:) :: vars
!    type(surface_t), intent(out), dimension(:,:,:) :: edgevars
!    integer :: nlx, nly, nlz
    integer :: n, ii, jj, kk, i, j

!    nlx = size(edgevars,1)
!    nly = size(edgevars,2)
!    nlz = size(edgevars,3)

ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
i = threadIdx%x
j = threadIdx%y
n = threadIdx%z

!    do kk = 1, nlz
!      do jj = 1, nly
!        do ii = 1, nlx

!         do n = 1, neq 
!        do j = 1, nx
!        do i = 1, nx
          edgevars_d(ii,jj,kk)%edges(1)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,1,j, n)
          edgevars_d(ii,jj,kk)%edges(2)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(nx,i,j,n)
          edgevars_d(ii,jj,kk)%edges(3)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,nx,j,n)
          edgevars_d(ii,jj,kk)%edges(4)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(1,i,j, n)
          edgevars_d(ii,jj,kk)%edges(5)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,j,1, n)
          edgevars_d(ii,jj,kk)%edges(6)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,j,nx,n)
!        end do
!        end do
!         end do
!
!        end do
!      end do
!    end do

!     if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)

  end subroutine extract_edges1

!!=========================================================================================
  ! Fill the interior edge data by extracting (i.e. copying) from the
  ! interior volume data. Quantities extracted are state vectors  (neq = 5)        

 attributes(global) subroutine extract_edges2()
!    type(element_t), intent(in), dimension(:,:,:) :: vars
!    type(surface_t), intent(out), dimension(:,:,:) :: edgevars
!    integer :: nlx, nly, nlz
    integer :: n, ii, jj, kk, i, j

!    nlx = size(edgevars,1)
!    nly = size(edgevars,2)
!    nlz = size(edgevars,3)

ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
i = threadIdx%x
j = threadIdx%y
n = threadIdx%z

!    do kk = 1, nlz
!      do jj = 1, nly
!        do ii = 1, nlx

!         do n = 1, neq 
!        do j = 1, nx
!        do i = 1, nx
          edgevars_d(ii,jj,kk)%edges(1)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,1,j, n)
          edgevars_d(ii,jj,kk)%edges(2)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(nx,i,j,n)
          edgevars_d(ii,jj,kk)%edges(3)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,nx,j,n)
          edgevars_d(ii,jj,kk)%edges(4)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(1,i,j, n)
          edgevars_d(ii,jj,kk)%edges(5)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,j,1, n)
          edgevars_d(ii,jj,kk)%edges(6)%vec_int(i,j,n) = vars_d(ii,jj,kk)%vec(i,j,nx,n)
!        end do
!        end do
!         end do
!
!        end do
!      end do
!    end do

!     if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)

  end subroutine extract_edges2

!!!=========================================================================================
!  ! Fill exterior edge data by getting the neighbor's interior values.
!  ! If the neighboring element is on the same MPI rank, then this is a
!  ! simple copy. Otherwise, the data is sent/received using MPI.
!  subroutine communicate_edgevars()
!!    type(surface_t), intent(out), dimension(:,:,:) :: edgevars
!    integer :: nlx, nly, nlz
!    integer :: ii, jj, kk, ip, jp, kp, im, jm, km, n 
!
!    nlx = size(edgevars,1)
!    nly = size(edgevars,2)
!    nlz = size(edgevars,3)
!    if (.not. allocated(bufsendS)) call allocate_buffers(nlx, nly, nlz)
!
!    ! Pack the interior edge data for those edges that interface with a
!    ! different MPI rank, and send the data off asynchronously.
!    call pack_edgevars( bufsendS, bufsendE, bufsendN, bufsendW)
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
!          edgevars(ii,jj,kk)%edges(1)%vec_ext(:,:,n) = edgevars(ii,jm,kk)%edges(3)%vec_int(:,:,n)
!          edgevars(ii,jj,kk)%edges(2)%vec_ext(:,:,n) = edgevars(ip,jj,kk)%edges(4)%vec_int(:,:,n)
!          edgevars(ii,jj,kk)%edges(3)%vec_ext(:,:,n) = edgevars(ii,jp,kk)%edges(1)%vec_int(:,:,n)
!          edgevars(ii,jj,kk)%edges(4)%vec_ext(:,:,n) = edgevars(im,jj,kk)%edges(2)%vec_int(:,:,n)
!          edgevars(ii,jj,kk)%edges(5)%vec_ext(:,:,n) = edgevars(ii,jj,km)%edges(6)%vec_int(:,:,n)
!          edgevars(ii,jj,kk)%edges(6)%vec_ext(:,:,n) = edgevars(ii,jj,kp)%edges(5)%vec_int(:,:,n)
!         end do
!
!        end do
!      end do
!    end do
!
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
!    call unpack_edgevars( bufrecvS, bufrecvE, bufrecvN, bufrecvW)
!  end subroutine communicate_edgevars
!
!
!!!=========================================================================================

subroutine communicate_edgevars_work()

    integer :: nlx, nly, nlz
    integer :: ii, jj, kk, ip, jp, kp, im, jm, km, n, i, j

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

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

         do n = 1, neq
        do j = 1, nx
        do i = 1, nx
          edgevars(ii,jj,kk)%edges(1)%vec_ext(i,j,n) = edgevars(ii,jm,kk)%edges(3)%vec_int(i,j,n)
          edgevars(ii,jj,kk)%edges(2)%vec_ext(i,j,n) = edgevars(ip,jj,kk)%edges(4)%vec_int(i,j,n)
          edgevars(ii,jj,kk)%edges(3)%vec_ext(i,j,n) = edgevars(ii,jp,kk)%edges(1)%vec_int(i,j,n)
          edgevars(ii,jj,kk)%edges(4)%vec_ext(i,j,n) = edgevars(im,jj,kk)%edges(2)%vec_int(i,j,n)
          edgevars(ii,jj,kk)%edges(5)%vec_ext(i,j,n) = edgevars(ii,jj,km)%edges(6)%vec_int(i,j,n)
          edgevars(ii,jj,kk)%edges(6)%vec_ext(i,j,n) = edgevars(ii,jj,kp)%edges(5)%vec_int(i,j,n)
         end do
        end do
        end do

        end do
      end do
    end do


end subroutine communicate_edgevars_work

!!!=========================================================================================

attributes(global) subroutine communicate_edges_work()

    integer :: nlx, nly, nlz
    integer :: ii, jj, kk, ip, jp, kp, im, jm, km, n, i, j

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    do kk = 1, nlz
!      do jj = 1, nly
!        do ii = 1, nlx
ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
i = threadIdx%x
j = threadIdx%y
n = threadIdx%z

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

!         do n = 1, neq
!        do j = 1, nx
!        do i = 1, nx
          edgevars_d(ii,jj,kk)%edges(1)%vec_ext(i,j,n) = edgevars_d(ii,jm,kk)%edges(3)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(2)%vec_ext(i,j,n) = edgevars_d(ip,jj,kk)%edges(4)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(3)%vec_ext(i,j,n) = edgevars_d(ii,jp,kk)%edges(1)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(4)%vec_ext(i,j,n) = edgevars_d(im,jj,kk)%edges(2)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(5)%vec_ext(i,j,n) = edgevars_d(ii,jj,km)%edges(6)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(6)%vec_ext(i,j,n) = edgevars_d(ii,jj,kp)%edges(5)%vec_int(i,j,n)
!         end do
!        end do
!        end do

!        end do
!      end do
!    end do


end subroutine communicate_edges_work
!!!=========================================================================================

attributes(global) subroutine communicate_edges_work1()

    integer :: nlx, nly, nlz
    integer :: ii, jj, kk, ip, jp, kp, im, jm, km, n, i, j

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    do kk = 1, nlz
!      do jj = 1, nly
!        do ii = 1, nlx
ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
i = threadIdx%x
j = threadIdx%y
n = threadIdx%z

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

!         do n = 1, neq
!        do j = 1, nx
!        do i = 1, nx
          edgevars_d(ii,jj,kk)%edges(1)%vec_ext(i,j,n) = edgevars_d(ii,jm,kk)%edges(3)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(2)%vec_ext(i,j,n) = edgevars_d(ip,jj,kk)%edges(4)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(3)%vec_ext(i,j,n) = edgevars_d(ii,jp,kk)%edges(1)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(4)%vec_ext(i,j,n) = edgevars_d(im,jj,kk)%edges(2)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(5)%vec_ext(i,j,n) = edgevars_d(ii,jj,km)%edges(6)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(6)%vec_ext(i,j,n) = edgevars_d(ii,jj,kp)%edges(5)%vec_int(i,j,n)
!         end do
!        end do
!        end do

!        end do
!      end do
!    end do


end subroutine communicate_edges_work1
!!!=========================================================================================

attributes(global) subroutine communicate_edges_work2()

    integer :: nlx, nly, nlz
    integer :: ii, jj, kk, ip, jp, kp, im, jm, km, n, i, j

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    do kk = 1, nlz
!      do jj = 1, nly
!        do ii = 1, nlx
ii = blockIdx%x
jj = blockIdx%y
kk = blockIdx%z
i = threadIdx%x
j = threadIdx%y
n = threadIdx%z

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

!         do n = 1, neq
!        do j = 1, nx
!        do i = 1, nx
          edgevars_d(ii,jj,kk)%edges(1)%vec_ext(i,j,n) = edgevars_d(ii,jm,kk)%edges(3)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(2)%vec_ext(i,j,n) = edgevars_d(ip,jj,kk)%edges(4)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(3)%vec_ext(i,j,n) = edgevars_d(ii,jp,kk)%edges(1)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(4)%vec_ext(i,j,n) = edgevars_d(im,jj,kk)%edges(2)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(5)%vec_ext(i,j,n) = edgevars_d(ii,jj,km)%edges(6)%vec_int(i,j,n)
          edgevars_d(ii,jj,kk)%edges(6)%vec_ext(i,j,n) = edgevars_d(ii,jj,kp)%edges(5)%vec_int(i,j,n)
!         end do
!        end do
!        end do

!        end do
!      end do
!    end do


end subroutine communicate_edges_work2

!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

  subroutine pack_edgevars()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

!    idx = 1
    do ke = 1, nlz
      do ie = 1, nlx
        do k = 1, nx
          do i = 1, nx
idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  1)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,1)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  2)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,2)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  3)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,3)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  4)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,4)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  5)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,5)
!            idx = idx + 1

          end do
        end do
      end do
    end do

!    idx = 1
    do ke = 1, nlz
      do je = 1, nly
        do k = 1, nx
          do j = 1, nx
idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,1)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  1)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,2)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  2)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,3)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  3)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,4)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  4)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,5)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  5)
!            idx = idx + 1
          end do
        end do
      end do
    end do

  end subroutine pack_edgevars

!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

  subroutine pack_edgevars1()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

!    idx = 1
    do ke = 1, nlz
      do ie = 1, nlx
        do k = 1, nx
          do i = 1, nx
idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  1)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,1)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  2)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,2)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  3)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,3)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  4)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,4)
            idx = idx + 1

            bufsendS(idx) = edgevars(ie,1,ke)%edges(1)%vec_int(i,k,  5)
            bufsendN(idx) = edgevars(ie,nly,ke)%edges(3)%vec_int(i,k,5)
!            idx = idx + 1

          end do
        end do
      end do
    end do

end subroutine pack_edgevars1

!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

 attributes(global) subroutine pack_edges1()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do ie = 1, nlx
!        do k = 1, nx
!          do i = 1, nx
i = threadIdx%x
k = threadIdx%y
ie = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  1)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,1)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  2)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,2)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  3)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,3)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  4)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,4)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  5)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,5)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

end subroutine pack_edges1
!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

 attributes(global) subroutine pack_edges11()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do ie = 1, nlx
!        do k = 1, nx
!          do i = 1, nx
i = threadIdx%x
k = threadIdx%y
ie = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  1)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,1)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  2)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,2)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  3)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,3)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  4)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,4)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  5)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,5)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

end subroutine pack_edges11
!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

 attributes(global) subroutine pack_edges12()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do ie = 1, nlx
!        do k = 1, nx
!          do i = 1, nx
i = threadIdx%x
k = threadIdx%y
ie = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  1)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,1)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  2)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,2)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  3)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,3)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  4)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,4)
            idx = idx + 1

            bufsendS_d(idx) = edgevars_d(ie,1,ke)%edges(1)%vec_int(i,k,  5)
            bufsendN_d(idx) = edgevars_d(ie,nly,ke)%edges(3)%vec_int(i,k,5)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

end subroutine pack_edges12

!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

  subroutine pack_edgevars2()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

!    idx = 1
    do ke = 1, nlz
      do je = 1, nly
        do k = 1, nx
          do j = 1, nx
idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,1)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  1)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,2)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  2)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,3)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  3)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,4)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  4)
            idx = idx + 1

            bufsendE(idx) = edgevars(nlx,je,ke)%edges(2)%vec_int(j,k,5)
            bufsendW(idx) = edgevars(1,je,ke)%edges(4)%vec_int(j,k,  5)
!            idx = idx + 1
          end do
        end do
      end do
    end do

  end subroutine pack_edgevars2

!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

 attributes(global) subroutine pack_edges2()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do je = 1, nly
!        do k = 1, nx
!          do j = 1, nx
j = threadIdx%x
k = threadIdx%y
je = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,1)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  1)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,2)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  2)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,3)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  3)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,4)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  4)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,5)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  5)
!            idx = idx + 1
!          end do
!        end do
!      end do
!    end do

  end subroutine pack_edges2
!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

 attributes(global) subroutine pack_edges21()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do je = 1, nly
!        do k = 1, nx
!          do j = 1, nx
j = threadIdx%x
k = threadIdx%y
je = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,1)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  1)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,2)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  2)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,3)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  3)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,4)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  4)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,5)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  5)
!            idx = idx + 1
!          end do
!        end do
!      end do
!    end do

  end subroutine pack_edges21
!!=========================================================================================
  ! Pack the edge data from the rank's elements into 1-dim communication
  ! buffers, for sending using MPI.

 attributes(global) subroutine pack_edges22()
!    type(surface_t), intent(in), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(out), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(out), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do je = 1, nly
!        do k = 1, nx
!          do j = 1, nx
j = threadIdx%x
k = threadIdx%y
je = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,1)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  1)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,2)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  2)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,3)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  3)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,4)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  4)
            idx = idx + 1

            bufsendE_d(idx) = edgevars_d(nlx,je,ke)%edges(2)%vec_int(j,k,5)
            bufsendW_d(idx) = edgevars_d(1,je,ke)%edges(4)%vec_int(j,k,  5)
!            idx = idx + 1
!          end do
!        end do
!      end do
!    end do

  end subroutine pack_edges22

!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
  subroutine unpack_edgevars()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

!    idx = 1
    do ke = 1, nlz
      do ie = 1, nlx
        do k = 1, nx
          do i = 1, nx
idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  1) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,1) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  2) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,2) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  3) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,3) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  4) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,4) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  5) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,5) = bufrecvN(idx)
!            idx = idx + 1

          end do
        end do
      end do
    end do

!    idx = 1
    do ke = 1, nlz
      do je = 1, nly
        do k = 1, nx
          do j = 1, nx
idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,1) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  1) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,2) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  2) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,3) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  3) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,4) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  4) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,5) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  5) = bufrecvW(idx)
!            idx = idx + 1

          end do
        end do
      end do
    end do

  end subroutine unpack_edgevars

!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
  subroutine unpack_edgevars1()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

!    idx = 1
    do ke = 1, nlz
      do ie = 1, nlx
        do k = 1, nx
          do i = 1, nx
idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  1) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,1) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  2) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,2) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  3) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,3) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  4) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,4) = bufrecvN(idx)
            idx = idx + 1

            edgevars(ie,1,ke)%edges(1)%vec_ext(i,k,  5) = bufrecvS(idx)
            edgevars(ie,nly,ke)%edges(3)%vec_ext(i,k,5) = bufrecvN(idx)
!            idx = idx + 1

          end do
        end do
      end do
    end do

end subroutine unpack_edgevars1

!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
 attributes(global) subroutine unpack_edges1()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do ie = 1, nlx
!        do k = 1, nx
!          do i = 1, nx
i = threadIdx%x
k = threadIdx%y
ie = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  1) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,1) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  2) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,2) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  3) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,3) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  4) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,4) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  5) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,5) = bufrecvN_d(idx)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

end subroutine unpack_edges1
!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
 attributes(global) subroutine unpack_edges11()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do ie = 1, nlx
!        do k = 1, nx
!          do i = 1, nx
i = threadIdx%x
k = threadIdx%y
ie = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  1) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,1) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  2) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,2) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  3) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,3) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  4) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,4) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  5) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,5) = bufrecvN_d(idx)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

end subroutine unpack_edges11
!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
 attributes(global) subroutine unpack_edges12()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do ie = 1, nlx
!        do k = 1, nx
!          do i = 1, nx
i = threadIdx%x
k = threadIdx%y
ie = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nlx*nx*nx*5)+((ie-1)*nx*nx*5)+((k-1)*nx*5)+((i-1)*5)+1
            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  1) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,1) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  2) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,2) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  3) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,3) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  4) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,4) = bufrecvN_d(idx)
            idx = idx + 1

            edgevars_d(ie,1,ke)%edges(1)%vec_ext(i,k,  5) = bufrecvS_d(idx)
            edgevars_d(ie,nly,ke)%edges(3)%vec_ext(i,k,5) = bufrecvN_d(idx)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

end subroutine unpack_edges12

!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
  subroutine unpack_edgevars2()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars,1)
    nly = size(edgevars,2)
    nlz = size(edgevars,3)

!    idx = 1
    do ke = 1, nlz
      do je = 1, nly
        do k = 1, nx
          do j = 1, nx
idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,1) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  1) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,2) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  2) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,3) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  3) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,4) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  4) = bufrecvW(idx)
            idx = idx + 1

            edgevars(nlx,je,ke)%edges(2)%vec_ext(j,k,5) = bufrecvE(idx)
            edgevars(1,je,ke)%edges(4)%vec_ext(j,k,  5) = bufrecvW(idx)
!            idx = idx + 1

          end do
        end do
      end do
    end do

  end subroutine unpack_edgevars2

!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
 attributes(global) subroutine unpack_edges2()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do je = 1, nly
!        do k = 1, nx
!          do j = 1, nx
j = threadIdx%x
k = threadIdx%y
je = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,1) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  1) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,2) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  2) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,3) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  3) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,4) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  4) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,5) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  5) = bufrecvW_d(idx)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

  end subroutine unpack_edges2
!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
 attributes(global) subroutine unpack_edges21()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do je = 1, nly
!        do k = 1, nx
!          do j = 1, nx
j = threadIdx%x
k = threadIdx%y
je = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,1) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  1) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,2) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  2) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,3) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  3) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,4) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  4) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,5) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  5) = bufrecvW_d(idx)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

  end subroutine unpack_edges21
!!=========================================================================================
  ! Unpack the data from the 1-dim communication buffer into edge_vars.
  ! NOTE: The velocity is a VECTOR quantity, so `speed dot normal` has
  !       a sign flip across the interface.
 attributes(global) subroutine unpack_edges22()
!  subroutine unpack_edgevars( bufferS, bufferE, bufferN, bufferW)
!    type(surface_t), intent(inout), dimension(:,:,:) :: edgevars
!    real(kind=real_kind), intent(in), dimension(:) :: bufferS, bufferN
!    real(kind=real_kind), intent(in), dimension(:) :: bufferE, bufferW
    integer :: nlx, nly, nlz
    integer :: i, j, k, ie, je, ke, idx

    nlx = size(edgevars_d,1)
    nly = size(edgevars_d,2)
    nlz = size(edgevars_d,3)

!    idx = 1
!    do ke = 1, nlz
!      do je = 1, nly
!        do k = 1, nx
!          do j = 1, nx
j = threadIdx%x
k = threadIdx%y
je = blockIdx%x
ke = blockIdx%y

idx = ((ke-1)*nly*nx*nx*5)+((je-1)*nx*nx*5)+((k-1)*nx*5)+((j-1)*5)+1
            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,1) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  1) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,2) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  2) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,3) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  3) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,4) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  4) = bufrecvW_d(idx)
            idx = idx + 1

            edgevars_d(nlx,je,ke)%edges(2)%vec_ext(j,k,5) = bufrecvE_d(idx)
            edgevars_d(1,je,ke)%edges(4)%vec_ext(j,k,  5) = bufrecvW_d(idx)
!            idx = idx + 1

!          end do
!        end do
!      end do
!    end do

  end subroutine unpack_edges22

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
    if (.not. allocated(bufsendS_d)) allocate(bufsendS_d(nSN))
    if (.not. allocated(bufrecvS_d)) allocate(bufrecvS_d(nSN))
    if (.not. allocated(bufsendE_d)) allocate(bufsendE_d(nEW))
    if (.not. allocated(bufrecvE_d)) allocate(bufrecvE_d(nEW))
    if (.not. allocated(bufsendN_d)) allocate(bufsendN_d(nSN))
    if (.not. allocated(bufrecvN_d)) allocate(bufrecvN_d(nSN))
    if (.not. allocated(bufsendW_d)) allocate(bufsendW_d(nEW))
    if (.not. allocated(bufrecvW_d)) allocate(bufrecvW_d(nEW))
  end subroutine allocate_buffers
!!=========================================================================================
  ! Allocate the send/recv buffers used in MPI communication
  subroutine allocate_buffers1(nlx, nly, nlz)
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
    if (.not. allocated(bufsendS_d)) allocate(bufsendS_d(nSN))
    if (.not. allocated(bufrecvS_d)) allocate(bufrecvS_d(nSN))
    if (.not. allocated(bufsendE_d)) allocate(bufsendE_d(nEW))
    if (.not. allocated(bufrecvE_d)) allocate(bufrecvE_d(nEW))
    if (.not. allocated(bufsendN_d)) allocate(bufsendN_d(nSN))
    if (.not. allocated(bufrecvN_d)) allocate(bufrecvN_d(nSN))
    if (.not. allocated(bufsendW_d)) allocate(bufsendW_d(nEW))
    if (.not. allocated(bufrecvW_d)) allocate(bufrecvW_d(nEW))
  end subroutine allocate_buffers1
!!=========================================================================================
  ! Allocate the send/recv buffers used in MPI communication
  subroutine allocate_buffers2(nlx, nly, nlz)
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
    if (.not. allocated(bufsendS_d)) allocate(bufsendS_d(nSN))
    if (.not. allocated(bufrecvS_d)) allocate(bufrecvS_d(nSN))
    if (.not. allocated(bufsendE_d)) allocate(bufsendE_d(nEW))
    if (.not. allocated(bufrecvE_d)) allocate(bufrecvE_d(nEW))
    if (.not. allocated(bufsendN_d)) allocate(bufsendN_d(nSN))
    if (.not. allocated(bufrecvN_d)) allocate(bufrecvN_d(nSN))
    if (.not. allocated(bufsendW_d)) allocate(bufsendW_d(nEW))
    if (.not. allocated(bufrecvW_d)) allocate(bufrecvW_d(nEW))
  end subroutine allocate_buffers2

end module prepare_edgevars_mod
