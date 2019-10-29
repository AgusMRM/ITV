MODULE healpix_types
  ! This module sets the types used in the Fortran 90 modules
  ! of the HEALPIX distribution and follows the example of Numerical Recipes
  !
  ! Benjamin D. Wandelt October 1997
  ! Eric Hivon June 1998
  ! Eric Hivon Oct  2001, edited to be compatible with 'F' compiler
  ! Eric Hivon July 2002, addition of i8b, i2b, i1b
  !                       addition of max_i8b, max_i2b and max_i1b
  
  INTEGER, PARAMETER, public :: i8b = SELECTED_INT_KIND(18)
  INTEGER, PARAMETER, public :: i4b = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER, public :: i2b = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER, public :: i1b = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER, public :: sp  = SELECTED_REAL_KIND(5,30)
  INTEGER, PARAMETER, public :: dp  = SELECTED_REAL_KIND(12,200)
  INTEGER, PARAMETER, public :: lgt = KIND(.TRUE.)
  INTEGER, PARAMETER, public :: spc = KIND((1.0_sp, 1.0_sp))
  INTEGER, PARAMETER, public :: dpc = KIND((1.0_dp, 1.0_dp))
  !
  INTEGER(I8B),  PARAMETER, public :: max_i8b = HUGE(1_i8b)
  INTEGER,       PARAMETER, public :: max_i4b = HUGE(1_i4b)
  INTEGER,       PARAMETER, public :: max_i2b = HUGE(1_i2b)
  INTEGER,       PARAMETER, public :: max_i1b = HUGE(1_i1b)
  REAL(kind=sp), PARAMETER, public :: max_sp  = HUGE(1.0_sp)
  REAL(kind=dp), PARAMETER, public :: max_dp  = HUGE(1.0_dp)


  ! Numerical Constant (Double precision)
  REAL(kind=dp), PARAMETER, public :: QUARTPI=0.785398163397448309615660845819875721049_dp
  REAL(kind=dp), PARAMETER, public :: HALFPI= 1.570796326794896619231321691639751442099_dp
  REAL(kind=dp), PARAMETER, public :: PI    = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER, public :: TWOPI = 6.283185307179586476925286766559005768394_dp
  REAL(kind=dp), PARAMETER, public :: FOURPI=12.56637061435917295385057353311801153679_dp
  REAL(kind=dp), PARAMETER, public :: SQRT2 = 1.41421356237309504880168872420969807856967_dp
  REAL(kind=dp), PARAMETER, public :: EULER = 0.5772156649015328606065120900824024310422_dp
  REAL(kind=dp), PARAMETER, public :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_dp
  REAL(kind=dp), PARAMETER, public :: TWOTHIRD = 0.6666666666666666666666666666666666666666_dp

  ! Maximum length of filenames
  integer, parameter :: filenamelen = 1024


  ! ---- Normalisation and convention ----
  ! normalisation of spin weighted functions
  real(kind=dp), parameter, public ::  KvS = 1.0_dp ! 1.0 : CMBFAST (Healpix 1.2)
  ! sign of Q
  real(kind=dp), parameter, public :: sgQ = -1.0_dp ! -1 : CMBFAST (Healpix 1.2)
  ! sign of spin weighted function !
  real(kind=dp), parameter, public :: SW1 = -1.0_dp ! -1 : Healpix 1.2, bug correction

!  ! normalisation of spin weighted functions
!  real(kind=dp), parameter, public ::  KvS = 2.0_dp ! 2.0 : KKS  (Healpix 1.1)
!  ! sign of Q
!  real(kind=dp), parameter, public :: sgQ = +1.0_dp ! +1 : KKS (Healpix 1.1)
!  ! sign of spin weighted function !
!  real(kind=dp), parameter, public :: SW1 = +1.0_dp ! +1 : Healpix 1.1

  real(kind=dp), parameter, public :: iKvS = 1.0_dp / KvS  ! inverse of KvS

END MODULE healpix_types
