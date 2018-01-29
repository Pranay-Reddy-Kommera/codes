!-------------------------------------------------------------------
! R.D. Nair IMAGe/NCAR 2016
! post processing 3D data 
! global constants and simulation parameters
!-------------------------------------------------------------------

module basic_mod
  implicit none

  ! all defns in this module are meant to be globally accessible
  public

  ! ----------------------------------------------------------------
  ! define global constants

  integer, parameter :: real_kind = kind(0.0D0)
  real(kind=real_kind), parameter :: pi = 2*acos(0.0D0)
  real(kind=real_kind), parameter :: twopi = 2*pi

  ! resolution parameters
  integer, parameter :: nod = 5    ! degree of the polynomial 
  integer, parameter :: nx = nod+1
  integer, Parameter :: nlg=nx  

  ! nex/y/z is the number of elements in x/y/z

   integer, parameter :: nex = 32
   integer, parameter :: ney = 16
   integer, parameter :: nez = 8 

  character(len=20), public   :: name_of_the_test 

  real(kind=real_kind) :: xmin, xmax, ymin, ymax, zmin, zmax

  ! size of each element, set in grid setup
  real(kind=real_kind) :: delx, dely, delz

  ! useful reparametrizations of nx, nex/y/z
  ! nx is the number of GLL points in each direction
  ! nux/y/z is the number of UNIQUE x/y/z collocation points (GLL points overlap!)
  ! (this corresponds to 'nt+1' in Ram's 2D advection code)
  integer, parameter :: nux = nod*nex+1, nuy = nod*ney+1, nuz = nod*nez+1

end module basic_mod
