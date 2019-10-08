program
        use modulos
        use OMP_lib
        implicit none
        type trivect
              real(DP),dimension(3) :: v
        endtype trivect

  type(trivect),dimension(:),allocatable :: grid_pixels
  real(DP),dimension(3) :: v
  real(DP),parameter              :: rmin=1000.
  integer,parameter           :: nthread   = 32 
  integer,parameter           :: nside     = 60   ! QUE ES ?
  integer,parameter           :: nbindist  = 10
  integer,parameter           :: nsideheal = 2    ! QUE ES ?
  real,parameter              :: rmax_vu=10. , rmin_vu=0.05
  real :: xe,ye,ze,rv,ncen
  integer                     :: kk1,kk2,ipx,nwrites
  xe = 250
  ye = 250
  ze = 250
  rv = 9.75
  ncen = 1
        
  call reader()
  
  allocate(dens_h(nbindist))
  
  npix = 12*( nsideheal )**2 
  allocate(grid_pixels(npix))

  lrmax_vu=alog10(rmax_vu)
  lrmin_vu =alog10(rmin_vu)
  dbindist_vu = (lrmax_vu - lrmin_vu) / nbindist


  dens   = -999.999
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE()

!$OMP SHARED()

  allocate(   dens_pix(npix))
















endprogram
