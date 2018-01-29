!-------------------------------------------------------------------
! R.Nair NCAR/scd 04/2004
! Rewritten 07/2016 by Francois Hebert, SIParCS/Cornell
!
! Various subroutines for saving simulation data to file.
!-------------------------------------------------------------------

module output_mod
  use basic_mod
  use element_mod

  implicit none
  private
  public :: save_global_state, save_ycut_frame, error_norms, print_diagnostics

  logical :: prepare_output_called = .false.

contains

  ! Some work must be done before 'real' data can be output. This work is
  ! gathered here, so it can be called by various 'real' output functions.
  subroutine prepare_output
    use grid_setup_mod, only: get_visualization_grid

    integer :: i
    real(kind=real_kind), dimension(nux) :: ax
    real(kind=real_kind), dimension(nuy) :: ay
    real(kind=real_kind), dimension(nuz) :: az

    if (prepare_output_called) return

    ! get the 1d grid in each x/y/z direction, with points de-duplicated
    call get_visualization_grid(ax,ay,az)

    ! NOTE: This command from Fortran 2008 is only supported on intel 16+.
    !       For running on Yellowstone with default modules, keep this
    !       commented out. For running on a system with gfortran or the
    !       latest intel compiler, can uncomment the line.. the code will
    !       then generate the output dir for you.
    !call execute_command_line("mkdir -p "//output_dir)

    3 format(e24.16) ! use 'e' rather than 'd' notation for python plotting
    open (unit = 25, file = output_dir//'/grid_geometry.dat')
    write(25,*) nx, nex, ney, nez
    close(25)
    open (unit = 35, file = output_dir//'/grid_xs.dat')
    write(35,3) (ax(i),i=1,nux)
    close(35)
    open (unit = 35, file = output_dir//'/grid_ys.dat')
    write(35,3) (ay(i),i=1,nuy)
    close(35)
    open (unit = 35, file = output_dir//'/grid_zs.dat')
    write(35,3) (az(i),i=1,nuz)
    close(35)

    prepare_output_called = .true.

  end subroutine prepare_output


  ! Save the global variable state to file.
  !   Duplicated gridpoints at element interfaces are not written
  !   to file twice. The lower (in x,y,z) value only is saved.
  !
  ! Creates files holding the grid information...
  ! * grid_geometry.dat -- nx,nex,ney,nez for NCL scripts
  ! * grid_(x,y,z)s.dat -- locations of GLL gridpoints in x,y,z
  !
  ! ...and files holding snapshots of psi...
  ! * data_full.dat -- (optional) full psi data
  ! * data_{X,Y,Z}cut.dat -- slice of the data at constant X,Y,Z
  !                          respectively. slice is at midpoint of
  !                          the domain in that coordinate.
  subroutine save_global_state()

!    type(element_t), intent(in), dimension(nex,ney,nez) :: vars
    real(kind=real_kind), dimension(nux,nuy,nuz) :: fglob
   real(kind=real_kind), dimension(nx,nx,nx,nex,ney,nez) :: dat3d      
    integer :: i, j, k, ie, je, ke, iu, ju, ku
    integer :: midx, midy, midz

    call prepare_output

    ! rewrite psi data as a simple 3D grid with no dupes
    do ke = 1, nez
      do k = 1, nx-1
        ku = (ke-1)*(nx-1) + k

        do je = 1, ney
          do j = 1, nx-1
            ju = (je-1)*(nx-1) + j

            do ie = 1, nex
              do i = 1, nx-1
                iu = (ie-1)*(nx-1) + i
                fglob(iu,ju,ku) = global_vars_phys(ie,je,ke)%thp(i,j,k) 
              end do
            end do

          end do
        end do

      end do
    end do


!! Global data on 3D domain 

    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex
          dat3d(:,:,:,ie,je,ke) = global_vars_phys(ie,je,ke)%thp(:,:,:) 
        end do
      end do
    end do


    ! doubly periodic boundary (extension)
    fglob(nux,:,:) = fglob(1,:,:)
    fglob(:,nuy,:) = fglob(:,1,:)
    fglob(:,:,nuz) = fglob(:,:,1)

    midx = int((nux+1)/2.0D0)
    midy = int((nuy+1)/2.0D0)
    midz = int((nuz+1)/2.0D0)

                print*, "Global data min/max",  minval(fglob), maxval(fglob) 

    3 format(e24.16) ! use 'e' rather than 'd' notation for python plotting
    open (unit = 30, file = output_dir//'/data_Xcut.dat')
    write(30,3) ((fglob(midx,i,j),i=1,nuy),j=1,nuz)
    close(30)
    open (unit = 30, file = output_dir//'/data_Ycut.dat')
    write(30,3) ((fglob(i,midy,j),i=1,nux),j=1,nuz)
    close(30)
    open (unit = 30, file = output_dir//'/data_Zcut.dat')
    write(30,3) ((fglob(i,j,midz),i=1,nux),j=1,nuy)
    close(30)
    if (output_full_volume_data) then
     !open (unit = 30, file = output_dir//'/data_full.dat')
     !write(30,3) (((fglob(i,j,k),i=1,nux),j=1,nuy),k=1,nuz)
     !close(30)
     print*, "writing full data" 
    !open (unit = 33, file = output_dir//'/nh3d_full.dat')
     open (unit = 33, file = output_dir//'/full_32x16x8nx6.dat')
     write(33,*) ((( (((dat3d(i,j,k,ie,je,ke),i=1,nx),j=1,nx),k=1,nx),ie=1,nex),je=1,ney),ke=1,nez)
     close(33) 

    end if

  end subroutine save_global_state

!!=====================================================================================
  Subroutine print_diagnostics(itn)

    integer, intent(in) :: itn 
!    type(element_t), intent(in), dimension(nex,ney,nez) :: vars
    real(kind=real_kind), dimension(nux,nuy,nuz) :: uf,vf,wf 
    integer :: i, j, k, ie, je, ke, iu, ju, ku
    integer :: midx, midy, midz

    ! rewrite psi data as a simple 3D grid with no dupes
    do ke = 1, nez
      do k = 1, nx-1
        ku = (ke-1)*(nx-1) + k

        do je = 1, ney
          do j = 1, nx-1
            ju = (je-1)*(nx-1) + j

            do ie = 1, nex
              do i = 1, nx-1
                iu = (ie-1)*(nx-1) + i
                uf(iu,ju,ku) = global_vars_phys(ie,je,ke)%u(i,j,k)
               !vf(iu,ju,ku) = vars(ie,je,ke)%spd(i,j,k)
                vf(iu,ju,ku) = global_vars_phys(ie,je,ke)%v(i,j,k)
                wf(iu,ju,ku) = global_vars_phys(ie,je,ke)%w(i,j,k) 
              end do
            end do

          end do
        end do

      end do
    end do

    ! doubly periodic boundary (extension)
    uf(nux,:,:) = uf(1,:,:)
    uf(:,nuy,:) = uf(:,1,:)
    uf(:,:,nuz) = uf(:,:,1)
    wf(nux,:,:) = wf(1,:,:)
    wf(:,nuy,:) = wf(:,1,:)
    wf(:,:,nuz) = wf(:,:,1)
    vf(nux,:,:) = vf(1,:,:)
    vf(:,nuy,:) = vf(:,1,:)
    vf(:,:,nuz) = vf(:,:,1)

     print*, "Iteration,  Time:", itn , itn*dt 
     print*, "Global u min/max:",  minval(uf), maxval(uf) 
     print*, "Global v min/max:",  minval(vf), maxval(vf) 
     print*, "Global w min/max:",  minval(wf), maxval(wf) 
     !print*, "Global speed  min/max",  minval(vf), maxval(vf) 
     print*, "        "

 End   Subroutine print_diagnostics

!!=====================================================================================
  ! Similar to save_global_state above, but only dumps the ycut component
  ! of the data, and appends a frame number to the data file. Intended to
  ! dump the data needed to generate a movie of the evolution.
  !
  ! (By using a separate subroutine specialized to operate only on the y-cut
  ! data, we can save computation and afford to dump more movie frames.)
  subroutine save_ycut_frame( frameid)

!    type(element_t), intent(in), dimension(nex,ney,nez) :: vars
    integer, intent(in) :: frameid
    real(kind=real_kind), dimension(nux,nuz) :: fcut
    integer :: i, j, k, ie, je, ke, iu, ju, ku
    character(len=16) :: buffer

    call prepare_output

    ! find the je,j pair that corresponds to ju
    ju = int((nuy+1)/2.0D0) ! middle of y dimension
    je = (ju-1)/(nx-1) + 1 ! integer div intentional
    j = ju - (je-1)*(nx-1)

    ! extract mid-y slice of psi data
    do ke = 1, nez
      do k = 1, nx-1
        ku = (ke-1)*(nx-1) + k

        do ie = 1, nex
          do i = 1, nx-1
            iu = (ie-1)*(nx-1) + i
            fcut(iu,ku) = global_vars_phys(ie,je,ke)%psi(i,j,k)
          end do
        end do

      end do
    end do

    ! doubly periodic boundary (extension)
    fcut(nux,:) = fcut(1,:)
    fcut(:,nuz) = fcut(:,1)

    3 format(e24.16) ! use 'e' rather than 'd' notation for python plotting
    write(buffer, "(I0.4)") frameid
    open (unit = 30, file = output_dir//'/frame_Ycut.'//trim(buffer)//'.dat')
    write(30,3) ((fcut(i,k),i=1,nux),k=1,nuz)
    close(30)

  end subroutine save_ycut_frame



  ! Error norms
  subroutine error_norms( vec_id)

!    type(element_t), intent(in), dimension(nex,ney,nez) ::  vars, vars_init 

    integer, intent(in) :: vec_id
    real(kind=real_kind), dimension(nx,nx,nx,nex,ney,nez) :: ref, err
    real(kind=real_kind) :: ref_l1, ref_l2, ref_li
    real(kind=real_kind) :: err_l1, err_l2, err_li
    integer :: ie, je, ke

    call prepare_output

   if (vec_id == 5) then
    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex
          ref(:,:,:,ie,je,ke) = global_vars_init(ie,je,ke)%rho(:,:,:) * global_vars_init(ie,je,ke)%the(:,:,:) 
          err(:,:,:,ie,je,ke) = global_vars(ie,je,ke)%rho(:,:,:) * global_vars(ie,je,ke)%the(:,:,:) - ref(:,:,:,ie,je,ke)
        end do
      end do
    end do
    print*, 'state_vec(5): rho*theta error' 
   endif 

   if (vec_id == 2) then
    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex
          ref(:,:,:,ie,je,ke) = global_vars_init(ie,je,ke)%u(:,:,:) * global_vars_init(ie,je,ke)%rho(:,:,:)
          err(:,:,:,ie,je,ke) = global_vars(ie,je,ke)%u(:,:,:)* global_vars(ie,je,ke)%rho(:,:,:) - ref(:,:,:,ie,je,ke)
        end do
      end do
    end do
    print*, 'state_vec(2): rho*u error'
   endif

    !sanity check for volume integrals 
    !print*, "normalized integral of 1.0 = ", volume_int(1+0*err)

    ! norms: L1, L2, Linf
    ref_l1 = volume_int(ref)
    ref_l2 = sqrt(volume_int(ref*ref))
    ref_li = maxval(abs(ref))
    err_l1 = volume_int(err)
    err_l2 = sqrt(volume_int(err*err))
    err_li = maxval(abs(err))

    ! write to file
    4 format(6e24.16) ! use 'e' rather than 'd' notation for python plotting
    open (unit = 25, file = output_dir//'/error_norms.dat')
    !write(25,*) '# ref_l1  ref_l2  ref_li  err_l1  err_l2  err_li'
    !write(25,4) ref_l1, ref_l2, ref_li, err_l1, err_l2, err_li
    write(25,*) '# L1_err, L2_err, Linf_err'
    write(25,4) err_l1/ref_l1, err_l2/ref_l2, err_li/ref_li 
    close(25)

    print*, ref_l1, ref_l2, ref_li, err_l1, err_l2, err_li
  end subroutine error_norms


  ! NORMALIZED volume integral over domain
  pure function volume_int(func)
    real(kind=real_kind) :: volume_int
    real(kind=real_kind), intent(in), dimension(nx,nx,nx,nex,ney,nez) :: func
    real(kind=real_kind) :: acc
    integer :: ie, je, ke

    acc = 0.0D0
    do ke = 1, nez
      do je = 1, ney
        do ie = 1, nex
          acc = acc + element_int(func(:,:,:,ie,je,ke))
        end do
      end do
    end do

    volume_int = acc / ((xmax-xmin) * (ymax-ymin) * (zmax-zmin))
  end function volume_int


  ! volume integral over one element
  pure function element_int(func)
!    use gauss_quadrature_mod, only: gllw
    real(kind=real_kind) :: element_int
    real(kind=real_kind), intent(in), dimension(nx,nx,nx) :: func
    real(kind=real_kind) :: acc
    integer :: i, j, k

    acc = 0.0D0
    do k = 1, nx
      do j = 1, nx
        do i = 1, nx
          acc = acc + func(i,j,k) * gllw(i) * gllw(j) * gllw(k)
        end do
      end do
    end do

    element_int = acc * (delx/2) * (dely/2) * (delz/2)
  end function element_int

end module output_mod
