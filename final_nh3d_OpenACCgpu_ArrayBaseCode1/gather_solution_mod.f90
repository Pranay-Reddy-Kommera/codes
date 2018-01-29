!-------------------------------------------------------------------
! 06/2016 Francois Hebert, SIParCS/Cornell
!
! Gather the solution onto a single MPI rank before output.
!
! Takes an array of grid_t, holding the vars on each MPI process'
! partition and returns a bigger array of grid_t on every element.
!
! Communicates psi,u,v,w
!
!-------------------------------------------------------------------

module gather_solution_mod
  use basic_mod
  use element_mod
  use mpi_mod

  implicit none
  private
  public :: gather_solution_onto_master_rank
  public :: Gather_EulerSoln_onto_MasterRank

  ! we do not know the size of this buffer at compile-time
  ! (number of elements on this MPI partition is determined at runtime)
  real(kind=real_kind), dimension(:), allocatable :: buffer

  ! we know size at compile-time, but we only want this buffer on one rank
  ! (size is nex*ney*nez -- the total number of elements)
  real(kind=real_kind), dimension(:), allocatable :: global_buffer

  integer, parameter :: number_of_gathered_vars = 5

contains
!!======================================================================
  ! gathers
  subroutine gather_solution_onto_master_rank()
!    type(element_t), intent(in), dimension(:,:,:) :: vars
    ! allocatable: we only want this on one rank
!    type(element_t), intent(out), dimension(:,:,:), allocatable :: global_vars
    integer :: nlx, nly, nlz

    nlx = size(vars_rho,4)
    nly = size(vars_rho,5)
    nlz = size(vars_rho,6)

    if (.not. allocated(buffer)) then
      allocate(buffer(number_of_gathered_vars*nlx*nly*nlz*nx*nx*nx))
    end if
    if (.not. allocated(global_buffer) .and. rank_is_master) then
      allocate(global_buffer(number_of_gathered_vars*nex*ney*nez*nx*nx*nx))
    end if

    ! all ranks pack data
    call pack_vars(buffer)

    ! all ranks send, one rank receives
    call gather_data_onto_master_rank(number_of_gathered_vars, buffer, global_buffer)

    ! one rank unpacks
    if (rank_is_master) then
!      if (.not. allocated(global_vars_rho)) call allocateGlobalData(nex, ney,nez)
      call unpack_vars(global_buffer)
    end if
  end subroutine gather_solution_onto_master_rank

!!======================================================================

  ! gathers (rho,the,u,v,w) 
  subroutine Gather_EulerSoln_onto_MasterRank(num)
!    type(element_t), intent(in), dimension(:,:,:) :: vars
    ! allocatable: we only want this on one rank
!    type(element_t), intent(out), dimension(:,:,:), allocatable :: global_vars
integer, intent(in) :: num
    integer :: nlx, nly, nlz

    nlx = size(vars_rho,4)
    nly = size(vars_rho,5)
    nlz = size(vars_rho,6)

    if (.not. allocated(buffer)) then
      allocate(buffer(number_of_gathered_vars*nlx*nly*nlz*nx*nx*nx))
    end if
    if (.not. allocated(global_buffer) .and. rank_is_master) then
      allocate(global_buffer(number_of_gathered_vars*nex*ney*nez*nx*nx*nx))
    end if

    ! all ranks pack data
    call pack_Euler_vars(buffer, num)

    ! all ranks send, one rank receives
    call gather_data_onto_master_rank(number_of_gathered_vars, buffer, global_buffer)

    ! one rank unpacks
    if (rank_is_master) then
!      if (.not. allocated(global_vars_rho)) call allocateGlobalData(nex, ney, nez)
      call unpack_Euler_vars(global_buffer,num)
    end if
  end subroutine Gather_EulerSoln_onto_MasterRank

!!======================================================================
  subroutine pack_Euler_vars(buffer,num)
!    type(element_t), intent(in), dimension(:,:,:) :: vars
    real(kind=real_kind), intent(out), dimension(:) :: buffer
integer, intent(in) :: num
    integer :: i, j, k, ie, je, ke, idx
    integer :: nlx, nly, nlz

    nlx = size(vars_rho,4)
    nly = size(vars_rho,5)
    nlz = size(vars_rho,6)

    idx = 1

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
          do k = 1, nx
            do j = 1, nx
              do i = 1, nx
if (num == 1) then
                buffer(idx) = vars_init_rho(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_init_the(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_init_u(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_init_v(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_init_w(i,j,k,ie,je,ke)
                idx = idx + 1
else if (num==2) then
                buffer(idx) = vars_rho(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_the(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_u(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_v(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_w(i,j,k,ie,je,ke)
                idx = idx + 1

end if
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine pack_Euler_vars

!!======================================================================
  subroutine unpack_Euler_vars(global_buffer, num)
!    type(element_t), intent(out), dimension(nex,ney,nez) :: global_vars
    real(kind=real_kind), intent(in), &
      dimension(number_of_gathered_vars*nex*ney*nez*nx*nx*nx) :: global_buffer
integer, intent(in) :: num
    integer :: i, j, k, ie, je, ke, iglobal, jglobal, idx
    integer :: rank, nlx, nly, nlz

    idx = 1

    ! data is concatenated rank-by-rank in global_buffer
    ! iterate over ranks to reconstruct the global_vars state
    do rank = 0, mpiprocs()-1
      nlx = numcolx(rank)
      nly = numcoly(rank)
      nlz = nez
      do ke = 1, nlz
        do je = 1, nly
          do ie = 1, nlx
            iglobal = ie + offsetx(rank)
            jglobal = je + offsety(rank)
            do k = 1, nx
              do j = 1, nx
                do i = 1, nx
if (num == 2) then
                  global_vars_rho(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_the(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_u(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_v(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_w(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
else if (num == 1) then
                  global_vars_init_rho(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_init_the(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_init_u(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_init_v(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_init_w(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1

end if
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine unpack_Euler_vars

!!======================================================================

  subroutine pack_vars(buffer)
!    type(element_t), intent(in), dimension(:,:,:) :: vars
    real(kind=real_kind), intent(out), dimension(:) :: buffer
    integer :: i, j, k, ie, je, ke, idx
    integer :: nlx, nly, nlz

    nlx = size(vars_rho,4)
    nly = size(vars_rho,5)
    nlz = size(vars_rho,6)

    idx = 1

    do ke = 1, nlz
      do je = 1, nly
        do ie = 1, nlx
          do k = 1, nx
            do j = 1, nx
              do i = 1, nx
                buffer(idx) = vars_thp(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_spd(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_u(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_v(i,j,k,ie,je,ke)
                idx = idx + 1
                buffer(idx) = vars_w(i,j,k,ie,je,ke)
                idx = idx + 1
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine pack_vars

!!======================================================================
  subroutine unpack_vars( global_buffer)
!    type(element_t), intent(out), dimension(nex,ney,nez) :: global_vars
    real(kind=real_kind), intent(in), &
      dimension(number_of_gathered_vars*nex*ney*nez*nx*nx*nx) :: global_buffer
    integer :: i, j, k, ie, je, ke, iglobal, jglobal, idx
    integer :: rank, nlx, nly, nlz

    idx = 1

    ! data is concatenated rank-by-rank in global_buffer
    ! iterate over ranks to reconstruct the global_vars state
    do rank = 0, mpiprocs()-1
      nlx = numcolx(rank)
      nly = numcoly(rank)
      nlz = nez
      do ke = 1, nlz
        do je = 1, nly
          do ie = 1, nlx
            iglobal = ie + offsetx(rank)
            jglobal = je + offsety(rank)
            do k = 1, nx
              do j = 1, nx
                do i = 1, nx
                  global_vars_phys_thp(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_phys_spd(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_phys_u(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_phys_v(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                  global_vars_phys_w(i,j,k,iglobal,jglobal,ke) = global_buffer(idx)
                  idx = idx + 1
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine unpack_vars

end module gather_solution_mod
