!-------------------------------------------------------------------
! 06/2016 Francois Hebert, SIParCS/Cornell
! (inspired from Robert Kloefkorn's MPI module for the 2D NH Euler code)
!
! Wrappers around the MPI communication and environment routines.
!-------------------------------------------------------------------

module mpi_mod
  use basic_mod
! use mpi
  include 'mpif.h'

! implicit none
  private

  public :: init, finalize, mpirank, mpicomm, mpiprocs
  public :: compute_partitioning, offsetx, offsety, numcolx, numcoly
  public :: send_data, recv_data, gather_data_onto_master_rank, barrier

  ! whether to do I/O on this rank
  logical, public :: rank_is_master

  ! Data describing the MPI environment
  type :: environment_t
    integer :: rank, comm, procs
  end type
  type(environment_t) :: mpienv

  ! How the domain's elements are split among the MPI processes
  ! - the elements are grouped into (cube-like) blocks
  ! - each MPI process handles 1 block
  ! The master process needs to know everything, and other processes need
  ! to know some parts... so simply make each processor know all data.
  integer :: blocksx, blocksy
  type :: block_t
    ! number of columns (vertical stack of elements) in x,y that form this block
    integer :: numcolx, numcoly
    ! zero-based offset of the first column in x,y. global indices are:
    ! do ix = offsetx+1, offsetx+numcolx
    integer :: offsetx, offsety
    ! ranks of neighbors
    integer :: rankS, rankE, rankN, rankW
  end type
  type(block_t), dimension(:), allocatable :: partition


  ! These numerical values are arbitrary, but must be distinct:
  integer, parameter :: tagS = 11
  integer, parameter :: tagE = 12
  integer, parameter :: tagN = 13
  integer, parameter :: tagW = 14

contains

  ! Initialize the MPI environment
  subroutine init()
    integer :: ierr
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpienv%procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpienv%rank, ierr)
    mpienv%comm = MPI_COMM_WORLD
    rank_is_master = (mpienv%rank==master_rank)
    if (rank_is_master) then
      print('(a,i0)'), "==> initialized MPI with number of procs = ", mpienv%procs
    end if
  end subroutine init

  ! Clean up the MPI environment
  subroutine finalize()
    integer :: ierr
    call MPI_FINALIZE(ierr)
  end subroutine finalize


  pure integer function mpirank()
    mpirank = mpienv%rank
  end function mpirank
  pure integer function mpicomm()
    mpicomm = mpienv%comm
  end function mpicomm
  pure integer function mpiprocs()
    mpiprocs = mpienv%procs
  end function mpiprocs


  ! Compute the distribution of the domain into blocks across the MPI processes.
  subroutine compute_partitioning
    integer :: ix, iy, rank, extra

    ! Get number of blocks in x,y from the number of MPI processes
    ! TODO: is there a nicer way than a case statement? (algorithmic?)
    select case (mpienv%procs)
    case (1)
      blocksx = 1
      blocksy = 1
    case (2)
      blocksx = 2
      blocksy = 1
    case (3)
      blocksx = 3
      blocksy = 1
    case (4)
      blocksx = 2
      blocksy = 2
    case (6)
      blocksx = 3
      blocksy = 2
    case (8)
      blocksx = 4
      blocksy = 2
    case (12)
      blocksx = 4
      blocksy = 3
    case (16)
      blocksx = 4
      blocksy = 4
    case (32)
      blocksx = 8
      blocksy = 4
    case (48)
      blocksx = 8
      blocksy = 6
    case (64)
      blocksx = 8
      blocksy = 8
    case (128)
      blocksx = 16
      blocksy = 8
    case (256)
      blocksx = 16
      blocksy = 16
    case (512)
      blocksx = 32
      blocksy = 16
    case (1024)
      blocksx = 32
      blocksy = 32
    case default
      print*, "no provision for # procs = ", mpienv%procs
      stop
    end select

    ! Sanity checks on the partitioning
    if (mpienv%procs > nex*nez) then
      print*, "no provision for using more MPI procs than elements"
      stop
    else if (blocksx > nex) then
      print*, "too few elements in x, for this number of MPI procs"
      stop
    else if (blocksy > ney) then
      print*, "too few elements in y, for this number of MPI procs"
      stop
    else if (blocksx*blocksy /= mpienv%procs) then
      print*, "internal error in setting blocks in x,y"
      stop
    end if

    ! Find each rank/block's partition info
    if (.not. allocated(partition)) allocate(partition(0:mpienv%procs-1))
    do rank = 0, size(partition)-1
      ix = modulo(rank, blocksx)
      iy = rank / blocksx ! integer div intended

      associate(part => partition(rank))
        ! set up blocks in x
        part%numcolx = nex / blocksx ! integer div intentional
        extra = modulo(nex, blocksx)
        if (ix < extra) then
          part%numcolx = part%numcolx + 1
          part%offsetx = part%numcolx * ix
        else
          ! offsetx = extra * (numcolx+1) + (ix-extra) * numcolx
          !         = extra + numcolx * ix
          part%offsetx = extra + part%numcolx * ix
        end if

        ! blocks in y
        part%numcoly = ney / blocksy
        extra = modulo(ney, blocksy)
        if (iy < extra) then
          part%numcoly = part%numcoly + 1
          part%offsety = part%numcoly * iy
        else
          part%offsety = extra + part%numcoly * iy
        end if

        part%rankS = modulo(rank - blocksx, mpienv%procs)
        part%rankE = iy*blocksx + modulo(rank + 1, blocksx)
        part%rankN = modulo(rank + blocksx, mpienv%procs)
        part%rankW = iy*blocksx + modulo(rank - 1, blocksx)
      end associate
    end do
  end subroutine compute_partitioning


  pure integer function offsetx(rank)
    integer, optional, intent(in) :: rank
    if (present(rank)) then
      offsetx = partition(rank)%offsetx
    else
      offsetx = partition(mpienv%rank)%offsetx
    end if
  end function offsetx

  pure integer function offsety(rank)
    integer, optional, intent(in) :: rank
    if (present(rank)) then
      offsety = partition(rank)%offsety
    else
      offsety = partition(mpienv%rank)%offsety
    end if
  end function offsety

  pure integer function numcolx(rank)
    integer, optional, intent(in) :: rank
    if (present(rank)) then
      numcolx = partition(rank)%numcolx
    else
      numcolx = partition(mpienv%rank)%numcolx
    end if
  end function numcolx

  pure integer function numcoly(rank)
    integer, optional, intent(in) :: rank
    if (present(rank)) then
      numcoly = partition(rank)%numcoly
    else
      numcoly = partition(mpienv%rank)%numcoly
    end if
  end function numcoly


  ! Wrapper around MPI_SEND. Handles bookkeeping and groups sends to the four
  ! neighbors of each block (in the SENW directions)
  subroutine send_data()
!    real(kind=real_kind), intent(in), dimension(:) :: bufsendS, bufsendN
!    real(kind=real_kind), intent(in), dimension(:) :: bufsendE, bufsendW
    integer :: ierr, ireq
    integer :: rS, rE, rN, rW
    integer :: nSN, nEW

    rS = partition(mpienv%rank)%rankS
    rE = partition(mpienv%rank)%rankE
    rN = partition(mpienv%rank)%rankN
    rW = partition(mpienv%rank)%rankW

    nSN = size(bufsendS)
    nEW = size(bufsendE)

!$acc host_data use_device(bufsendS, bufsendE, bufsendN, bufsendW)
    call MPI_ISEND(bufsendS, nSN, MPI_DOUBLE_PRECISION, rS, tagS, mpienv%comm, ireq, ierr)
    call MPI_ISEND(bufsendE, nEW, MPI_DOUBLE_PRECISION, rE, tagE, mpienv%comm, ireq, ierr)
    call MPI_ISEND(bufsendN, nSN, MPI_DOUBLE_PRECISION, rN, tagN, mpienv%comm, ireq, ierr)
    call MPI_ISEND(bufsendW, nEW, MPI_DOUBLE_PRECISION, rW, tagW, mpienv%comm, ireq, ierr)
!$acc end host_data
  end subroutine send_data


  ! Wrapper around MPI_RECV. Handles bookkeeping and groups receives from the
  ! four neighbors of each block (in the SENW directions)
  subroutine recv_data()
!    real(kind=real_kind), intent(out), dimension(:) :: bufrecvS, bufrecvN
!    real(kind=real_kind), intent(out), dimension(:) :: bufrecvE, bufrecvW
    integer :: ierr, istat(MPI_STATUS_SIZE)
    integer :: rS, rE, rN, rW
    integer :: nSN, nEW

    rS = partition(mpienv%rank)%rankS
    rE = partition(mpienv%rank)%rankE
    rN = partition(mpienv%rank)%rankN
    rW = partition(mpienv%rank)%rankW

    nSN = size(bufrecvS)
    nEW = size(bufrecvE)

!$acc host_data use_device(bufrecvS, bufrecvE, bufrecvN, bufrecvW)
    call MPI_RECV(bufrecvS, nSN, MPI_DOUBLE_PRECISION, rS, tagN, mpienv%comm, istat, ierr)
    call MPI_RECV(bufrecvE, nEW, MPI_DOUBLE_PRECISION, rE, tagW, mpienv%comm, istat, ierr)
    call MPI_RECV(bufrecvN, nSN, MPI_DOUBLE_PRECISION, rN, tagS, mpienv%comm, istat, ierr)
    call MPI_RECV(bufrecvW, nEW, MPI_DOUBLE_PRECISION, rW, tagE, mpienv%comm, istat, ierr)
!$acc end host_data
  end subroutine recv_data


  ! Wrapper around MPI_GATHERV. Handles bookkeeping.
  subroutine gather_data_onto_master_rank(num_vars,buffer, global_buffer)
    Integer, intent(in) :: num_vars 
    real(kind=real_kind), intent(in), dimension(:) :: buffer
    real(kind=real_kind), intent(out), dimension(:) :: global_buffer

    logical :: has_extra_points
    integer, dimension(0:mpienv%procs-1) :: sizes, offsets
    integer :: i, ierr

    ! sanity check buffer size
    has_extra_points = modulo(size(global_buffer), nex*ney*nez*nx*nx*nx) /= 0
    if (rank_is_master .and. has_extra_points) then
      print*, "got bad global_buffer in gather_data"
      print*, "size(global_buffer) = ", size(global_buffer)
      print*, "expected size = ", num_vars*nx*nx*nx*nex*ney*nez
      stop
    end if

    do i = 0, size(sizes)-1
      sizes(i) = numcolx(i)*numcoly(i)*nez * nx*nx*nx * num_vars
    end do
    offsets(0) = 0
    do i = 1, size(sizes)-1
      offsets(i) = offsets(i-1) + sizes(i-1)
    end do

    call MPI_GATHERV(buffer, size(buffer), MPI_DOUBLE_PRECISION, &
                     global_buffer, sizes, offsets, MPI_DOUBLE_PRECISION, &
                     master_rank, mpienv%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      print*, "got error"
      stop
    end if
  end subroutine gather_data_onto_master_rank


  ! Trivial wrapper around MPI_BARRIER, useful when debugging MPI things.
  subroutine barrier()
    integer :: ierr
    call MPI_BARRIER(mpienv%comm, ierr)
  end subroutine barrier

end module mpi_mod

