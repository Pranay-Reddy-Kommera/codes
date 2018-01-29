!-------------------------------------------------------------------
! R.D. Nair IMAGe/NCAR 2010
! Changes 06/2016 by Francois Hebert, SIParCS/Cornell
!
! global constants and simulation parameters
!-------------------------------------------------------------------

module basic_mod
  use cudafor
  implicit none

  ! all defns in this module are meant to be globally accessible
  public

  ! ----------------------------------------------------------------
  ! define global constants

  integer, parameter :: real_kind = kind(0.0D0)
  real(kind=real_kind), parameter :: pi = 2*acos(0.0D0)
  real(kind=real_kind), parameter :: twopi = 2*pi

  integer, parameter :: testcase_id_uniflow_test = 1000
  integer, parameter :: testcase_id_warm_bubble  = 1001
  integer, parameter :: testcase_id_nh_igw       = 1002
  integer, parameter :: testcase_id_steady_state = 1003


  ! -- Choose Parameters -------------------------------------------
  ! control parameters for the simulation -- to be set by user

  logical, parameter :: use_dim_split_scheme = .false.
! logical, parameter :: use_dim_split_scheme = .true.
  logical, parameter :: use_noflux_bc_top_bottom = .true.

  ! resolution parameters
  integer, parameter :: nod = 3    ! degree of the polynomial 
  integer, parameter :: neq = 5    ! # of equations in the sys 

  ! nex/y/z is the number of elements in x/y/z

 !! Config for warm bubble 
 !integer, parameter :: testcase = testcase_id_warm_bubble  
 !integer, parameter :: nex = 16
 !integer, parameter :: ney = 16
 !integer, parameter :: nez = 16
 !real(kind=real_kind), parameter :: dt = 0.05D0

 !! Config for inertia-gravity wave 
   integer, parameter :: testcase = testcase_id_nh_igw            
   integer, parameter :: nex = 64
   integer, parameter :: ney = 64
   integer, parameter :: nez = 16 
   real(kind=real_kind), parameter :: dt = 0.10D0

 !! Config for steady-state flow       
 ! integer, parameter :: testcase = testcase_id_steady_state      
 ! integer, parameter :: nex = 48 
 ! integer, parameter :: ney = 12
 ! integer, parameter :: nez = 15
 ! real(kind=real_kind), parameter :: dt = 0.75D0

  character(len=20), public   :: name_of_the_test 

  ! data output parameters
  integer, parameter :: master_rank = 0 ! this rank will do I/O
  character(len=*), parameter :: output_dir = "data"
 !logical, parameter :: output_full_volume_data = .true.
  logical, parameter :: output_full_volume_data = .false.
  logical, parameter :: save_movie_frames = .false.
  integer, parameter :: movie_frame_freq = 100


  ! ----------------------------------------------------------------
  ! other parameters

  ! spatial and temporal domain of the simulation
  ! these are set by each testcase as appropriate
  real(kind=real_kind) :: endtime
  real(kind=real_kind) :: xmin, xmax, ymin, ymax, zmin, zmax

  ! size of each element, set in grid setup
  real(kind=real_kind) :: delx, dely, delz

  ! useful reparametrizations of nx, nex/y/z
  ! nx is the number of GLL points in each direction
  integer, parameter :: nx = nod+1
  ! nux/y/z is the number of UNIQUE x/y/z collocation points (GLL points overlap!)
  ! (this corresponds to 'nt+1' in Ram's 2D advection code)
  integer, parameter :: nux = nod*nex+1, nuy = nod*ney+1, nuz = nod*nez+1
  real(kind=real_kind), public, dimension(nx) :: gllp, gllw
  real(kind=real_kind), public, dimension(nx,nx) :: der, weakder

  real(kind=real_kind), device, public, dimension(nx) :: gllp_d, gllw_d
  real(kind=real_kind), device, public, dimension(nx,nx) :: der_d, weakder_d
  real(kind=real_kind), device :: delx_d, dely_d, delz_d

! allocatable because buffer length depends on the number of MPI
  ! partitions, which is determined at runtime.
  real(kind=real_kind), dimension(:), allocatable :: bufsendS, bufrecvS
  real(kind=real_kind), dimension(:), allocatable :: bufsendE, bufrecvE
  real(kind=real_kind), dimension(:), allocatable :: bufsendN, bufrecvN
  real(kind=real_kind), dimension(:), allocatable :: bufsendW, bufrecvW
  real(kind=real_kind), device, dimension(:), allocatable :: bufsendS_d, bufrecvS_d
  real(kind=real_kind), device, dimension(:), allocatable :: bufsendE_d, bufrecvE_d
  real(kind=real_kind), device, dimension(:), allocatable :: bufsendN_d, bufrecvN_d
  real(kind=real_kind), device, dimension(:), allocatable :: bufsendW_d, bufrecvW_d

integer, parameter :: number_of_send_recv_vars = 5


end module basic_mod
