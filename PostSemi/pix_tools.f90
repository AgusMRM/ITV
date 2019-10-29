module pix_tools
!==================================================================
!    Various subroutines for converting angles and unit vectors
!    to pixel numbers in both the NESTED and RING schemes,
!    as well as looking for neighbours and making pixel queries
!
!    last update : Oct 4, 2002, EH, computation of pixel vertex in 
!        pix2vec_ring and pix2vec_nest
!==================================================================

  USE healpix_types
  IMPLICIT none

  INTEGER(KIND=i4b), private, PARAMETER :: ns_max=8192 ! 2^13 : largest nside available

  !initialise array x2pix, y2pix and pix2x, pix2y used in several routines
  integer(KIND=i4b), private, save, dimension(128) :: x2pix=0,y2pix=0

  integer(KIND=i4b), private, save, dimension(0:1023) :: pix2x=0, pix2y=0

  private

  public :: remove_dipole, &
       & query_strip, &
       & query_polygon, &
       & query_triangle, &
       & query_disc, &
       & pix2ang_ring, pix2vec_ring, ang2pix_ring, vec2pix_ring, &
       & pix2ang_nest, pix2vec_nest, ang2pix_nest, vec2pix_nest, &
       & convert_nest2ring, convert_ring2nest, &
       & convert_inplace_real, convert_inplace_int, &
       & nest2ring, ring2nest, xy2pix_nest, pix2xy_nest, &
       & mk_pix2xy, mk_xy2pix, &
       & neighbours_nest, &
       & next_in_line_nest, &
       & ang2vec, vec2ang, &
       & npix2nside, nside2npix, &
       & surface_triangle, angdist, vect_prod

  public :: intrs_intrv, in_ring, ring_num ! arcane usage
  public :: getdisc_ring  ! obsolete
  public :: chpix_nest, introspixtion, spreadpix_nest, Int2Binary  ! ADDED
  

contains

  !==============================================================
  subroutine remove_dipole( nside, map, nmaps, ordering, degree, multipoles, cos_theta_cut, fmissval)
    !============================================================
    ! remove_dipole( nside, map, nmaps, ordering, degree, multipoles, cos_theta_cut, fmissval)
    !
    ! removes monopole (and dipole) from a map
    !
    ! Nside:     I4,       IN   : Healpix resolution parameter
    ! map:       R4, array,INOUT: Heapix map (see Notes below)
    ! ordering:  I4,       IN:   Healpix scheme 1:RING, 2: NESTED
    ! degree:    I4,       IN:   multipole to remove, 1: monopole, 2: monopole and dipole
    ! multipoles:R8, array,OUT:  value of monopole and dipole
    ! cos_theta_cut:R8,    IN:   value of symmetric cut off around equator
    ! fmissval:  R4, Option, IN: value used to flag bad pixel on input, default=-1.6375e30
    !
    ! note : if degree= 1, or 2, the map is modified on output
    !     * the monopole (and dipole) is/are removed
    !     * pixels within the symmetric cut parameterized 
    !       by cos_theta_cut are set to fmissval (or its default value)
    !  if degree = 0, nothing is done
    !  all other values of degree are invalid
    !
    ! v1.0, EH, Caltech, Jan-2002, based on homonyme IDL routine
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    implicit none
    integer(kind=i4b),                  intent(in)    :: nside
    real   (kind=sp), dimension(0:),    intent(inout) :: map
    integer(kind=i4b),                  intent(in)    :: nmaps
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP), dimension(0:),    intent(out)   :: multipoles
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=SP),                   intent(in), optional :: fmissval


    integer(kind=i4b)                 :: ipix, npix, n_mult
    integer(kind=i4b)                 :: i
    logical(lgt)                      :: dodipole
    real(kind=dp)                     :: flag, temp, wmin
    real(kind=dp), dimension(1:3)     :: vec
    real(kind=dp), dimension(0:3)     :: b, wdiag
    real(kind=dp), dimension(0:3,0:3) :: mat, umat, vmat
    real(kind=dp), dimension(0:3)     :: dmultipoles
    real(kind=SP)                     :: fmiss_effct
    character(len=*), parameter :: code = "REMOVE_DIPOLE"
    real(kind=sp),    parameter :: fbad_value = -1.6375e30_sp
    real(kind=dp)                     :: theta1, theta2
    integer(kind=i4b)                 :: ncpix, ncp
    integer(kind=i4b), dimension(:), allocatable :: cut_pixel
    !============================================================

    npix = nside2npix(nside)
    multipoles = 0.0_dp
    n_mult = (degree)**2
    fmiss_effct = fbad_value
    if (present(fmissval)) fmiss_effct = fmissval

    if (degree == 0) then
       print*," No monopole nor dipole removal"
       return
    elseif (degree == 1) then
       dodipole = .false.
    else if (degree == 2) then
       dodipole = .true.
    else
       print*,code//"> degree can only be "
       print*,"      1: monopole (l=0) removal or "
       print*,"      2: monopole and dipole (l=0,1) removal"
       print*,code//"> ABORT ! "
       stop
    endif

    !----------------------------------------------
    ! flag out pixels within symmetric galactic cut
    !----------------------------------------------
    if (cos_theta_cut > 0.0_dp) then
       ncpix = npix * cos_theta_cut * 1.1 + 10* nside
       allocate(cut_pixel(0:ncpix))
       theta1 = acos(cos_theta_cut)
       theta2 = PI - theta1
       call query_strip(nside, theta1, theta2, cut_pixel, ncp, nest=ordering-1)
!        map(cut_pixel(ncpix)) = fmiss_effct
       map(cut_pixel(0:ncp-1)) = fmiss_effct ! bug correction
       deallocate(cut_pixel)
    endif

    !----------------------------------------------
    ! generate least square linear system
    !----------------------------------------------
    mat = 0.0_dp
    b   = 0.0_dp
    do ipix = 0, npix-1

       ! flag = 1 for good values
       !      = 0 for bad values
       flag = 1.0_dp
       !     if ( abs(map(ipix) - fmiss_effct) < 1.e-5*fmiss_effct ) flag = 0.0_dp
!        if ( abs(map(ipix) - fmiss_effct) < 1.e-5*fmiss_effct ) goto 20
       if ( abs(map(ipix) - fmiss_effct) <= abs(1.e-5*fmiss_effct) ) goto 20

       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
       endif

       ! construct vector T*(1,x,y,z)
       temp = map(ipix) * flag
       b(0) = b(0) + temp
       if (dodipole) then
          ! computes dipole basis functions
          b(1:3) = b(1:3) + temp * vec(1:3)
       endif

       ! construct matrix (1,x,y,z)#(1,x,y,z)
       mat(0,0) = mat(0,0) + flag
       if (dodipole) then
          do i = 1, 3
             mat(i,0)   = mat(i,0)   + vec(i) * flag
             mat(i,1:3) = mat(i,1:3) + vec(i) * vec(1:3) * flag
          enddo
       endif

20     continue
    enddo

    ! first row = first column (and vice versa)
    mat(0,1:3) = mat(1:3,0)


    !----------------------------------------------
    ! solve system    mat . (mono, dip_x, dip_y, dip_z) = b
    !----------------------------------------------

    if (dodipole) then
       ! SVD decomposition
       umat = mat
       call dsvdcmp(umat, 4, 4, 4, 4, wdiag, vmat)
       ! thresholding
       wmin = maxval(wdiag)* 1.e-6_dp
       where (wdiag < wmin)
          wdiag = 0.0_dp
       end where
       ! back substitution
       call dsvbksb(umat, wdiag, vmat, 4, 4, 4, 4, b, dmultipoles)

    else
       dmultipoles(0) = b(0) / mat(0,0) ! average temperature
    endif

    !----------------------------------------------
    ! remove monopole and dipole
    !----------------------------------------------
    do ipix = 0, npix-1

!        if ( abs(map(ipix) - fmiss_effct) < 1.e-5*fmiss_effct ) goto 10
       if ( abs(map(ipix) - fmiss_effct) < abs(1.e-5*fmiss_effct) ) goto 10

       map(ipix) = map(ipix) - dmultipoles(0)
       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
          map(ipix) = map(ipix) - SUM( dmultipoles(1:3) * vec(1:3))
       endif

10     continue

    enddo

    multipoles(0:n_mult-1) = dmultipoles(0:n_mult-1)

    !

    return
  end subroutine remove_dipole

  !=======================================================================
  subroutine query_strip ( nside, theta1, theta2, listpix, nlist, nest)
  !=======================================================================
    ! query_strip ( nside, theta1, theta2, listpix, nlist, nest)
    !
    ! finds pixels having a colatitude (measured from North Pole): 
    !  theta1 < colatitude < theta2
    !  with 0 <= theta1 < theta2 <= Pi
    ! if theta2 < theta1 then pixels with
    ! 0 <= colatitude < theta2   or   theta1 < colatitude < Pi are returned
    ! 
    ! nside :          I4,       input,  resolution parameter
    ! theta1, theta2 : DP,       input,  bounds
    ! listpix :        I4 array, input,  list of pixels
    ! nlist :          I4,       output, number of pixels in list
    ! nest  :          I4 opt.,  input, =0 : RING scheme, =1 : NESTED scheme
    !
    ! v1.0, EH, Caltech, Jan-2002
    !=======================================================================
    integer(kind=I4B), intent(in)                    :: nside
    real(kind=DP),     intent(in)                    :: theta1, theta2
    integer(kind=I4B), intent(out), dimension(0:)    :: listpix
    integer(kind=I4B), intent(out)                   :: nlist
    integer(kind=I4B), intent(in),  optional         :: nest


    integer(kind=I4B)                  :: npix, nstrip
    integer(kind=I4B)                  :: iz, ip, is, irmin, irmax
    integer(kind=I4B)                  :: ilist, nir, nlost, list_size
    real(kind=DP)                      :: phi0, dphi
    real(kind=DP), dimension(1:4)      :: colrange
    character(len=*), parameter        :: code = "query_strip"
    integer(kind=I4B), dimension(:), allocatable :: listir

    !=======================================================================
    list_size = size(listpix)
    !     ---------- check inputs ----------------
    npix = nside2npix(nside)
    if (npix < 0) then
       print*,code//"> Nside should be a power of 2"
       print*,code//"> current value = ",nside
       STOP "> program abort "
    endif

    if    (theta1 < 0.0_dp .or. theta1 > PI .or. &
         & theta2 < 0.0_dp .or. theta2 > PI) then
       write(unit=*,fmt="(a)") code//"> the colatitudes are in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       print*,code//"> current value = ", theta1, theta2
       STOP "> program abort "
    endif

    if (theta1 <= theta2) then
       nstrip = 1
       colrange(1:2*nstrip) = (/ theta1, theta2 /)
    else
       nstrip = 2
       colrange(1:2*nstrip) = (/ 0.0_dp, theta2, theta1, PI/)
    endif
       

    !     ------------- loop on strips ---------------------
    ilist = -1
    allocate(listir(0:4*nside-1))
    do is =0, nstrip-1
       irmin = ring_num(nside, cos(colrange(2*is+1)) )
       irmax = ring_num(nside, cos(colrange(2*is+2)) )


       !     ------------- loop on ring number ---------------------
       do iz = irmin, irmax
          !        ------- finds pixels in the ring ---------
          phi0 = 0
          dphi = PI
          call in_ring(nside, iz, phi0, dphi, listir, nir, nest)

          nlost = ilist + nir + 1 - list_size
          if ( nlost > 0 ) then
             print*,code//"> listpix is too short, it will be truncated at ",nir
             print*,"                         pixels lost : ", nlost
             nir = nir - nlost
          endif
          do ip = 0, nir-1
             ilist = ilist + 1
             listpix(ilist) = listir(ip)
          enddo

       enddo
       
    enddo

    !     ------ total number of pixel in the strip --------
    nlist = ilist + 1

    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir) 

    return
  end subroutine query_strip
  !=======================================================================
  subroutine query_polygon ( nside, vlist, nv, listpix, nlist, nest, inclusive)
    !=======================================================================
    ! finds pixels that lie within a CONVEX polygon defined by its vertex on the sphere
    !
    ! nside             : IN
    ! vlist(1:3, 0:n-1) : IN, list of vertices
    ! nv                : IN, number of vertices to be used (nv <= n)
    ! listpix           : OUT
    ! nlist             : OUT
    ! nest              : IN, OPTIONAL
    ! inclusive         : IN, OPTIONAL
    ! 
    ! algorithm:
    !   the polygon is divided into triangles
    !   vertex 0 belongs to all the triangles
    !
    ! v1.0, EH, Caltech, Dec-2001
    !=======================================================================
    USE num_rec, ONLY : isort
    integer(kind=I4B), intent(in)                    :: nside
    real(kind=DP),     intent(in),  dimension(1:,1:) :: vlist
    integer(kind=I4B), intent(in)                    :: nv
    integer(kind=I4B), intent(out), dimension(0:)    :: listpix
    integer(kind=I4B), intent(out)                   :: nlist
    integer(kind=I4B), intent(in),  optional         :: nest

    real(kind=DP),     dimension(:), pointer         :: vp0, vp1, vp2
    real(kind=DP),     dimension(1:3)                :: vo
    real(kind=DP),     dimension(:,:), allocatable, target   :: vvlist
    real(kind=DP),     dimension(:),   allocatable   :: ss
    real(kind=DP)                                    :: surface, fsky, hand
    integer(kind=I4B)                                :: npix, n_in_trg, ilist, ntl
    integer(kind=I4B)                                :: n_remain
    integer(kind=I4B)                                :: i0, i1, i2, i, k
    integer(kind=I4B)                                :: np, nm, nlow
    integer(kind=I4B)                                :: list_size, nlost
    integer(kind=I4B), dimension(1:1)                :: ix
    integer(kind=I4B), dimension(:), allocatable     :: templist
    character(len=*), parameter :: code = "QUERY_POLYGON"
    integer(kind=i4b)                                :: inclusive

    !=======================================================================

    list_size = size(listpix)
    n_remain = nv
    ilist = -1
    listpix(:) = -1
    
    if (n_remain < 3) then
       print*,code//"> The polygon should have at least 3 vertices"
       stop
    endif

    if (size(vlist) < n_remain*3) then
       print*,code//"> ",size(vlist)/3," vertices are given"
       print*,code//"> expected : ",n_remain
       stop
    endif

    allocate(vvlist(1:3, 1:n_remain))
    vvlist = vlist

    !-----------------------------------------------------------------
    ! check that the polygon is convex or has only one concave vertex
    !-----------------------------------------------------------------
    if (n_remain == 3) goto 1000 ! a triangle is always convex

    allocate(ss(1:n_remain))
    do i1 = 0, n_remain-1
       i0 = modulo(i1-1, n_remain)  ! in [0, n_remain-1]
       i2 = modulo(i1+1, n_remain)  ! in [0, n_remain-1]

       vp0 => vvlist(:,i0+1)
       vp1 => vvlist(:,i1+1)
       vp2 => vvlist(:,i2+1)

       ! computes the handedness   (v0 x v2) . v1  for each vertex v1
       call vect_prod(vp0, vp2, vo)
       hand = dot_product(vo, vp1)
       ss(i1+1) = sign(1.0_dp, hand) ! either +1 or -1
    enddo
    np = count( ss > 0.0_dp) ! number of vertex with positive handedness
    nm = n_remain - np

    nlow = min(np,nm)
    if (nlow == 0) goto 1000
    if (nlow == 1) then
       ! only one concave vertex
       if (np == 1) then
          ix = maxloc(ss)
       else
          ix = minloc(ss)
       endif
       ! rotate pixel list to put that vertex in #0
       vvlist = cshift(vvlist, ix(1)-1, dim = 2)
    endif
    if (nlow > 1) then
       ! more than 1 concave vertex
       print*,"***************************************"
       print*,"The polygon is not convex"
       print*," and has more than one concave vertex"
       print*,"The result is unpredictable"
       print*,"***************************************"
    endif
    deallocate(ss)

1000 continue

    !--------------------------------------------
    ! fill the polygon, one triangle at a time
    !--------------------------------------------
    npix = nside2npix(nside)
    do
       ! build triangle from vertices #0, n-2, and n-1
       vp0 => vvlist(1:3,1)
       vp1 => vvlist(1:3,n_remain-1)
       vp2 => vvlist(1:3,n_remain  )

       ! computes its surface
       call surface_triangle( vp0, vp1, vp2, surface )
       fsky = surface/FOURPI
       n_in_trg = npix * (fsky * 1.4) + 12*nside
       n_in_trg = min(n_in_trg, npix)
       
       ! find pixels within triangle
       allocate(templist(0:n_in_trg-1))
       call query_triangle( nside, vp0, vp1, vp2, templist, ntl, nest=nest,inclusive=inclusive )

       ! merge new list with existing one
       nlost = ilist + ntl + 1 - list_size
       if ( nlost > 0 ) then
          print*,code//"> listpix is too short, it will be truncated at ",list_size
          print*,"                         pixels lost : ", nlost
          print*, list_size
          ntl = ntl - nlost
       endif
       do i = 0, ntl - 1
          ilist = ilist + 1
          listpix(ilist) = templist(i)
       enddo
       deallocate(templist)

       if (n_remain == 3) exit

       ! prune vertex list
       n_remain = n_remain - 1

    enddo
    deallocate(vvlist)

    !-------------------------
    ! make final pixel list
    !-------------------------
    nlist = ilist + 1

    ! sort final list
    call isort(nlist, listpix)

    ! remove redondant pixels
    ! (we keep 0th element of the list)
    k = 1
    do i = 1, nlist-1
       if (listpix(i) > listpix(i-1)) then
          listpix(k) = listpix(i)
          k = k + 1
       endif
    enddo
    nlist = k

    ! pad end of list with -1
    listpix(nlist:list_size-1) = -1

    return
  end subroutine query_polygon
  !=======================================================================
  subroutine query_triangle ( nside, v1, v2, v3, listpix, nlist, nest, inclusive)
    !=======================================================================
    !
    !    query_triangle ( nside, v1, v2, v3, listpix, nlist, nest, inclusive)
    !    --------------
    !     nside       = resolution parameter (a power of 2)
    !     v1, v2, v3  = 3D vector location of the 3 triangle vertices
    !     list_pix    = list of pixel lying in the triangle
    !     nlist       = number of pixels in the list
    !     nest  (OPT), :0 by default, the output list is in RING scheme
    !                  if set to 1, the output list is in NESTED scheme
    !     inclusive (OPT) , :0 by default, only the pixels whose center 
    !                       lie in the triangle are listed on output
    !                  if set to 1, all pixels overlapping the triangle are output
    !
    !     NB : the dimension of the listpix array is fixed in the calling 
    !     routine and should be large enough for the specific configuration
    !
    !
    ! v1.0, EH, Caltech, Nov-Dec-2001
    ! v1.1,  Aug-Sep-2002 : added nest and inclusive
    !=======================================================================
    integer(kind=I4B), intent(in) :: nside
    real(kind=DP),     intent(in),  dimension(1:)  :: v1, v2, v3
    integer(kind=I4B), intent(out), dimension(0:)  :: listpix
    integer(kind=I4B), intent(out)                 :: nlist
    integer(kind=I4B), intent(in), optional        :: nest
    integer(kind=I4B), intent(in), optional        :: inclusive

    integer(kind=I4B) :: npix, ilist !, ip1, ip2, ip3
    integer(kind=I4B) :: iz, irmin, irmax
    integer(kind=I4B) :: n12, n123a, n123b, ndom
    integer(kind=I4B), dimension(:),   allocatable  :: listir
    logical(kind=LGT) :: test1, test2, test3
    logical(kind=LGT) :: test1a, test1b, test2a, test2b, test3a, test3b
    real(kind=DP) :: dth1, dth2, determ, sdet
    real(kind=DP) :: zmax, zmin, z1max, z1min, z2max, z2min, z3max, z3min
    real(kind=DP) :: z, zz
    real(kind=DP) :: tgth, st
    real(kind=DP) :: offset, sin_off
    real(kind=DP), dimension(1:3,1:3) :: vv, vo
    real(kind=DP), dimension(1:3) :: sprod, sto, phi0i, tgthi
    real(kind=DP), dimension(1:3) :: dc
    real(kind=DP), dimension(1:2,1:3) :: dom
    real(kind=DP), dimension(1:4) :: dom12, dom123a, dom123b
    real(kind=DP), dimension(1:6) :: alldom
    real(kind=DP) :: a_i, b_i, phi0, dphiring
    integer(kind=I4B) :: idom, nir, ip
    integer(kind=I4B) :: status
    integer(kind=I4B) :: j
    character(len=*), parameter :: code = "QUERY_TRIANGLE"
    integer(kind=I4B) :: list_size, nlost
    logical(LGT)      :: do_inclusive

    !=======================================================================

    list_size = size(listpix)

    npix = nside2npix(nside)
    if (npix < 0) then 
       print*,"Invalid Nside = ",nside
       stop
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif

    !   ! sort pixels by number
    !   ip1 = MIN(ipix1, ipix2, ipix3)
    !   ip3 = MAX(ipix1, ipix2, ipix3)
    !   ip2 = ipix1 + ipix2 + ipix3 - ip1 - ip3

    !   !     ---------- check inputs ----------------
    !   if (ip1 < 0 .or. ip3 > npix-1) then
    !      write(unit=*,fmt="(a)") " > Non valid choice for pixel number :"
    !      write(unit=*,fmt="(a)") " > ",ipix1,ipix2,ipix3
    !      write(unit=*,fmt="(a,i2,i10)") " > valid range : ",0,npix-1
    !      nlist = 0
    !      listpix(0) = -1
    !   endif

    !   call pix2vec_ring( nside, ip1, vv(1:3,1))
    !   call pix2vec_ring( nside, ip2, vv(1:3,2))
    !   call pix2vec_ring( nside, ip3, vv(1:3,3))

    vv(1:3,1) = v1(1:3) / sqrt(dot_product(v1,v1))
    vv(1:3,2) = v2(1:3) / sqrt(dot_product(v2,v2))
    vv(1:3,3) = v3(1:3) / sqrt(dot_product(v3,v3))

    !     --------- allocate memory -------------
    ALLOCATE( listir(0: 4*nside-1), STAT = status)
    if (status /= 0) then
       write(unit=*,fmt="(a)") code//" > can not allocate memory for listir :"
       STOP " > program abort "
    endif

    dth1 = 1.0_dp / (3.0_dp*REAL(nside,kind=dp)**2)
    dth2 = 2.0_dp / (3.0_dp*REAL(nside,kind=dp))


    ! determ = (vect1 X vect2) . vect3
    ! determines the left(<0)/right(>0) handedness of the triangle
    determ = vv(1,1)*vv(2,2)*vv(3,3) + vv(1,2)*vv(2,3)*vv(3,1) + vv(1,3)*vv(2,1)*vv(3,2) &
         & - vv(3,1)*vv(2,2)*vv(1,3) - vv(3,2)*vv(2,3)*vv(1,1) - vv(3,3)*vv(2,1)*vv(1,2)

    if (abs(determ) < 1.e-20_dp) then 
       print*,' ************************************************************'
       print*,' The triangle is degenerated (2 of the vertices are antipodal)'
       print*,' The query can not be performed '
       print*,' ************************************************************'
       stop
    endif

    !   print*,determ
    sdet = SIGN(1.0_dp,determ) ! = +1 or -1, the sign of determ

    ! scalar product of vertices vectors
    sprod(1) = dot_product(vv(1:3,2),vv(1:3,3))
    sprod(2) = dot_product(vv(1:3,3),vv(1:3,1))
    sprod(3) = dot_product(vv(1:3,1),vv(1:3,2))

    ! vector orthogonal to the great circle containing the vertex doublet
    call vect_prod(vv(1:3,2), vv(1:3,3), vo(1:3,1))
    call vect_prod(vv(1:3,3), vv(1:3,1), vo(1:3,2))
    call vect_prod(vv(1:3,1), vv(1:3,2), vo(1:3,3))

    ! normalize the orthogonal vector
    vo(1:3,1) = vo(1:3,1) /  SQRT(SUM(vo(1:3,1)**2))
    vo(1:3,2) = vo(1:3,2) /  SQRT(SUM(vo(1:3,2)**2))
    vo(1:3,3) = vo(1:3,3) /  SQRT(SUM(vo(1:3,3)**2))

    ! test presence of poles in the triangle
    zmax = -1.0_dp
    zmin =  1.0_dp
    test1 = (vo(3,1) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 2-3
    test2 = (vo(3,2) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 1-2
    test3 = (vo(3,3) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 1-3
    if (test1 .and. test2 .and. test3) then 
       zmax = 1.0_dp ! north pole in the triangle
    endif
    if ((.not.test1) .and. (.not.test2) .and. (.not.test3)) then
       zmin = -1.0_dp ! south pole in the triangle
    endif

    ! look for northernest and southernest points in the triangle
    !   node(1,2) = vector of norm=1, in the plane defined by (1,2) and with z=0
    test1a = ((vv(3,3) - sprod(1) * vv(3,2)) >= 0.0_dp) ! segment 2-3 : -vector(3) . node(2,3)
    test1b = ((vv(3,2) - sprod(1) * vv(3,3)) >= 0.0_dp) !                vector(2) . node(2,3)
    test2a = ((vv(3,3) - sprod(2) * vv(3,1)) >= 0.0_dp) ! segment 1-3 : -vector(3) . node(1,3)
    test2b = ((vv(3,1) - sprod(2) * vv(3,3)) >= 0.0_dp) !                vector(1) . node(1,3)
    test3a = ((vv(3,2) - sprod(3) * vv(3,1)) >= 0.0_dp) ! segment 1-2 : -vector(2) . node(1,2)
    test3b = ((vv(3,1) - sprod(3) * vv(3,2)) >= 0.0_dp) !                vector(1) . node(1,2)

    ! sin of theta for orthogonal vector
    sto(1:3) = SQRT( (1.0_dp-vo(3,1:3))*(1.0_dp+vo(3,1:3)) )

    ! for each segment (=side of the triangle) the extrema are either 
    ! - the 2 vertices 
    ! - one of the vertices and a point within the segment

    ! segment 2-3
    z1max = vv(3,2)
    z1min = vv(3,3)
    if ( test1a .EQV. test1b ) then
       zz = sto(1)
       if ((vv(3,2)+vv(3,3)) >= 0.0_dp) then
          z1max =  zz
       else
          z1min = -zz
       endif
    endif

    ! segment 1-3
!     z2max = vv(3,1)
!     z2min = vv(3,3)
    z2max = vv(3,3)
    z2min = vv(3,1)
    if ( test2a .EQV. test2b ) then
       zz = sto(2)
       if ((vv(3,1)+vv(3,3)) >= 0.0_dp) then
          z2max =  zz
       else
          z2min = -zz
       endif
    endif

    ! segment 1-2
    z3max = vv(3,1)
    z3min = vv(3,2)
    if ( test3a .EQV. test3b ) then
       zz = sto(3)
       if ((vv(3,1)+vv(3,2)) >= 0.0_dp) then
          z3max =  zz
       else
          z3min = -zz
       endif
    endif

    zmax = MAX(z1max, z2max, z3max, zmax)
    zmin = MIN(z1min, z2min, z3min, zmin)

    ! if we are inclusive, move the upper point up, and the lower point down, by a half pixel size
    offset = 0.0_dp
    sin_off = 0.0_dp
    if (do_inclusive) then
       offset = PI / (4.0_dp*nside) ! half pixel size
       sin_off = sin(offset)
       zmax = min( 1.0_dp, cos( acos(zmax) - offset) )
       zmin = max(-1.0_dp, cos( acos(zmin) + offset) )
    endif

    !   print*,"zmin, zmax ",zmin,zmax

    ! northernest and sourthernest ring number
    irmin = ring_num(nside, zmax)
!!!  irmin = MAX(1, irmin - 1) ! start from a higher point, to be safe
    irmax = ring_num(nside, zmin)
!!!  irmax = MIN(4*nside-1, irmax + 1) ! go down to a lower point

    ilist = -1
    !    print*,"irmin, irmax ",irmin,irmax

    ! -------- loop on the rings -------------------------

    tgthi(1:3) = -1.0e30_dp * vo(3,1:3)
    phi0i(1:3) =  0.0_dp
    do j=1,3 
       if (sto(j) > 1.0e-10_dp) then
          tgthi(j) = -vo(3,j) / sto(j)  ! -cotan(theta_orth)
          phi0i(j) = ATAN2(vo(2,j),vo(1,j))
       endif
    enddo
    !   print*,tgthi,phi0i

    ! the triangle boundaries are geodesics : intersection of the sphere with plans going thru (0,0,0)
    ! if we are inclusive, the boundaries are the intersecion of the sphere with plans pushed outward
    ! by sin(offset)
    do iz = irmin, irmax
       if (iz <= nside-1) then      ! north polar cap
          z = 1.0_dp  - REAL(iz,kind=dp)**2 * dth1
       else if (iz <= 3*nside) then    ! tropical band + equat.
          z = REAL(2*nside-iz,kind=dp) * dth2
       else 
          z = - 1.0_dp + REAL(4*nside-iz,kind=dp)**2 * dth1
       endif
       ! computes the 3 intervals described by the 3 great circles
       st = SQRT((1.0_dp - z)*(1.0_dp + z))
       tgth = z / st ! cotan(theta_ring)
!        dc(1:3)  = tgthi(1:3) * tgth - sdet * sin_off / (sto(1:3) * st)
       dc(1:3)  = tgthi(1:3) * tgth - sdet * sin_off / ((sto(1:3)+1.e-30_dp) * st) ! sto is slightly offset to avoid division by 0

       do j=1,3
          if (dc(j)*sdet <= -1.0_dp) then  ! the whole iso-latitude ring is on the right side of the great circle
             dom(1:2, j) = (/ 0.0_dp, twopi /) 
          else if (dc(j)*sdet >= 1.0_dp) then ! all on the wrong side
             dom(1:2, j) = (/ -1.000001_dp, -1.0_dp /) * j 
          else ! some is good, some is bad
             dom(1:2, j) = MODULO( phi0i(j) + (ACOS(dc(j)) * sdet) * (/-1.0_dp, 1.0_dp /), twopi)
          endif
       enddo

       ! identify the intersections (0,1,2 or 3) of the 3 intervals
       call intrs_intrv( dom(1:2,1), dom(1:2,2), dom12, n12)
       if (n12 == 0) goto 20
       if (n12 == 1) then
          call intrs_intrv( dom(1:2,3), dom12, dom123a, n123a)
          if (n123a == 0) goto 20
          alldom(1:2*n123a) = dom123a(1:2*n123a)
          ndom = n123a ! 1 or 2
       endif
       if (n12 == 2) then
          call intrs_intrv( dom(1:2,3), dom12(1:2), dom123a, n123a)
          call intrs_intrv( dom(1:2,3), dom12(3:4), dom123b, n123b)
          ndom = n123a + n123b ! 0, 1, 2 or 3
          if (ndom == 0) goto 20
          if (n123a /= 0) alldom(1:2*n123a) = dom123a(1:2*n123a)
          if (n123b /= 0) alldom(2*n123a+1:2*ndom)  = dom123b(1:2*n123b)
          if (ndom > 3) then
             print*,code//"> too many intervals found"
          endif
       endif
       do idom=0,ndom-1
          a_i = alldom(2*idom+1)
          b_i = alldom(2*idom+2)
          phi0 = (a_i + b_i) * 0.5_dp
          dphiring = (b_i - a_i) * 0.5_dp
          if (dphiring < 0.0_dp) then
             phi0 = phi0 + pi
             dphiring = dphiring + pi
          endif

          !        ------- finds pixels in the triangle on that ring ---------
          call in_ring(nside, iz, phi0, dphiring, listir, nir, nest=nest)

          !        ----------- merge pixel lists -----------
          nlost = ilist + nir + 1 - list_size
          if ( nlost > 0 ) then
             print*,code//"> listpix is too short, it will be truncated at ",nir
             print*,"                         pixels lost : ", nlost
             print*, list_size
             nir = nir - nlost
          endif
          do ip = 0, nir-1
             ilist = ilist + 1
             listpix(ilist) = listir(ip)
          enddo
       enddo
20     continue
    enddo !-----------------------------------------

    !     ------ total number of pixel in the disc --------
    nlist = ilist + 1

    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir) 

    return
  end subroutine query_triangle

  !=======================================================================
  subroutine getdisc_ring ( nside, vector0, radius, listpix, nlist)
    integer(kind=I4B), intent(in)                 :: nside
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in)                 :: radius
    integer(kind=I4B), intent(out), dimension(0:) :: listpix
    integer(kind=I4B), intent(out)                :: nlist
  !=======================================================================
    
    print*,"-------------------------------------------------------------"
    print*,"WARNING : the routine getdisc_ring is now obsolete"
    print*,"  Use "
    print*," call query_disc(nside, vector0, radius_radian, listpix, nlist [, nest]) "
    print*,"  instead (same module)"
    print*," It lets you choose the indexing scheme (RING or NESTED) "
    print*," used for the outuput"
    print*,"-------------------------------------------------------------"

    call query_disc(nside, vector0, radius, listpix, nlist)

    return
  end subroutine getdisc_ring

  !=======================================================================
  subroutine query_disc ( nside, vector0, radius, listpix, nlist, nest, inclusive)
    !=======================================================================
    !
    !      query_disc (Nside, Vector0, Radius, Listpix, Nlist[, Nest, Inclusive])
    !      ----------
    !      routine for pixel query in the RING or NESTED scheme
    !      all pixels within an angular distance Radius of the center 
    !
    !     Nside    = resolution parameter (a power of 2)
    !     Vector0  = central point vector position (x,y,z in double precision)
    !     Radius   = angular radius in RADIAN (in double precision)
    !     Listpix  = list of pixel closer to the center (angular distance) than Radius
    !     Nlist    = number of pixels in the list
    !     nest  (OPT), :0 by default, the output list is in RING scheme
    !                  if set to 1, the output list is in NESTED scheme
    !     inclusive (OPT) , :0 by default, only the pixels whose center 
    !                       lie in the triangle are listed on output
    !                  if set to 1, all pixels overlapping the triangle are output
    !
    !      * all pixel numbers are in {0, 12*Nside*Nside - 1}
    !     NB : the dimension of the listpix array is fixed in the calling 
    !     routine and should be large enough for the specific configuration
    !
    !      lower level subroutines called by getdisc_ring : 
    !       (you don't need to know them)
    !      ring_num (nside, ir)
    !      --------
    !      in_ring(nside, iz, phi0, dphi, listir, nir, nest=nest)
    !      -------
    !
    ! v1.0, EH, TAC, ??
    ! v1.1, EH, Caltech, Dec-2001
    !=======================================================================
    integer(kind=I4B), intent(in)                 :: nside
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in)                 :: radius
    integer(kind=I4B), intent(out), dimension(0:) :: listpix
    integer(kind=I4B), intent(out)                :: nlist
    integer(kind=I4B), intent(in), optional       :: nest
    integer(kind=I4B), intent(in), optional       :: inclusive

    INTEGER(KIND=I4B) :: irmin, irmax, ilist, iz, ip, nir, npix
    REAL(KIND=DP) :: norm_vect0
    REAL(KIND=DP) :: x0, y0, z0, radius_eff
    REAL(KIND=DP) :: a, b, c, cosang
    REAL(KIND=DP) :: dth1, dth2
    REAL(KIND=DP) :: phi0, cosphi0, cosdphi, dphi
    REAL(KIND=DP) :: rlat0, rlat1, rlat2, zmin, zmax, z
    INTEGER(KIND=I4B), DIMENSION(:),   ALLOCATABLE  :: listir
    INTEGER(KIND=I4B) :: status
    character(len=*), parameter :: code = "QUERY_DISC"
    integer(kind=I4B) :: list_size, nlost
    logical(LGT) :: do_inclusive

    !=======================================================================

    list_size = size(listpix)
    !     ---------- check inputs ----------------
    npix = 12 * nside * nside

    if (radius < 0.0_dp .or. radius > PI) then
       write(unit=*,fmt="(a)") code//"> the angular radius is in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       STOP "> program abort "
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif
          
    !     --------- allocate memory -------------
    ALLOCATE( listir(0: 4*nside-1), STAT = status)
    if (status /= 0) then
       write(unit=*,fmt="(a)") code//"> can not allocate memory for listir :"
       STOP "> program abort "
    endif

    dth1 = 1.0_dp / (3.0_dp*real(nside,kind=dp)**2)
    dth2 = 2.0_dp / (3.0_dp*real(nside,kind=dp))

    radius_eff = radius
    if (do_inclusive) then
       ! increase radius by half pixel size
       radius_eff = radius + PI / (4.0_dp*nside)
    endif
    cosang = COS(radius_eff)

    !     ---------- circle center -------------
    norm_vect0 =  SQRT(DOT_PRODUCT(vector0,vector0))
    x0 = vector0(1) / norm_vect0
    y0 = vector0(2) / norm_vect0
    z0 = vector0(3) / norm_vect0

    phi0=0.0_dp
    if ((x0/=0.0_dp).or.(y0/=0.0_dp)) phi0 = ATAN2 (y0, x0)  ! in ]-Pi, Pi]
    cosphi0 = COS(phi0)
    a = x0*x0 + y0*y0

    !     --- coordinate z of highest and lowest points in the disc ---
    rlat0  = ASIN(z0)    ! latitude in RAD of the center
    rlat1  = rlat0 + radius_eff
    rlat2  = rlat0 - radius_eff
    if (rlat1 >=  halfpi) then
       zmax =  1.0_dp
    else
       zmax = SIN(rlat1)
    endif
    irmin = ring_num(nside, zmax)
    irmin = MAX(1, irmin - 1) ! start from a higher point, to be safe

    if (rlat2 <= -halfpi) then
       zmin = -1.0_dp
    else
       zmin = SIN(rlat2)
    endif
    irmax = ring_num(nside, zmin)
    irmax = MIN(4*nside-1, irmax + 1) ! go down to a lower point

    ilist = -1

    !     ------------- loop on ring number ---------------------
    do iz = irmin, irmax

       if (iz <= nside-1) then      ! north polar cap
          z = 1.0_dp  - real(iz,kind=dp)**2 * dth1
       else if (iz <= 3*nside) then    ! tropical band + equat.
          z = real(2*nside-iz,kind=dp) * dth2
       else 
          z = - 1.0_dp + real(4*nside-iz,kind=dp)**2 * dth1
       endif

       !        --------- phi range in the disc for each z ---------
       b = cosang - z*z0
       c = 1.0_dp - z*z
       if ((x0==0.0_dp).and.(y0==0.0_dp)) then
          cosdphi=-1.0_dp
          dphi=PI
          goto 500
       endif
       cosdphi = b / SQRT(a*c)
       if (ABS(cosdphi) <= 1.0_dp) then
          dphi = ACOS (cosdphi) ! in [0,Pi]
       else
          if (cosphi0 < cosdphi) goto 1000 ! out of the disc
          dphi = PI ! all the pixels at this elevation are in the disc
       endif
500    continue

       !        ------- finds pixels in the disc ---------
       call in_ring(nside, iz, phi0, dphi, listir, nir, nest)

       !        ----------- merge pixel lists -----------
       nlost = ilist + nir + 1 - list_size
       if ( nlost > 0 ) then
          print*,code//"> listpix is too short, it will be truncated at ",nir
          print*,"                         pixels lost : ", nlost
          nir = nir - nlost
       endif
       do ip = 0, nir-1
          ilist = ilist + 1
          listpix(ilist) = listir(ip)
       enddo

1000   continue
    enddo

    !     ------ total number of pixel in the disc --------
    nlist = ilist + 1


    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir) 

    return
  end subroutine query_disc
  !=======================================================================
  function ring_num (nside, z) result(ring_num_result)
    !=======================================================================
    !     returns the ring number in {1, 4*nside-1}
    !     from the z coordinate
    !=======================================================================
    INTEGER(KIND=I4B) :: ring_num_result
    REAL(KIND=DP), INTENT(IN) :: z
    INTEGER(KIND=I4B), INTENT(IN) :: nside

    INTEGER(KIND=I4B) :: iring
    !=======================================================================

    !     ----- equatorial regime ---------
    iring = NINT( nside*(2.0_dp-1.500_dp*z))

    !     ----- north cap ------
    if (z > twothird) then
       iring = NINT( nside* SQRT(3.0_dp*(1.0_dp-z)))
       if (iring == 0) iring = 1
    endif

    !     ----- south cap -----
    if (z < -twothird   ) then
       iring = NINT( nside* SQRT(3.0_dp*(1.0_dp+z)))
       if (iring == 0) iring = 1
       iring = 4*nside - iring
    endif

    ring_num_result = iring

    return
  end function ring_num
  !=======================================================================
  subroutine in_ring (nside, iz, phi0, dphi, listir, nir, nest)
    !=======================================================================
    !     returns the list of pixels in RING or NESTED scheme (listir) 
    !     and their number (nir)
    !     with latitude in [phi0-dphi, phi0+dphi] on the ring ir 
    !     (in {1,4*nside-1})
    !     the pixel id-numbers are in {0,12*nside^2-1}
    !     the indexing is RING, unless NEST is set to 1
    !=======================================================================
    integer(kind=i4b), intent(in)                 :: nside, iz
    integer(kind=i4b), intent(out)                :: nir
    real(kind=dp),     intent(in)                 :: phi0, dphi
    integer(kind=i4b), intent(out), dimension(0:) :: listir
    integer(kind=i4b), intent(in), optional       :: nest

!     logical(kind=lgt) :: conservative = .true.
    logical(kind=lgt) :: conservative = .false.
    logical(kind=lgt) :: take_all, to_top, do_ring

    integer(kind=i4b) :: ip_low, ip_hi, i, in, inext
    integer(kind=i4b) :: npix, nr, nir1, nir2, ir, ipix1, ipix2, kshift, ncap
    real(kind=dp)     :: phi_low, phi_hi, shift
    !=======================================================================

    take_all = .false.
    to_top   = .false.
    do_ring  = .true.
    if (present(nest)) then
       do_ring = (nest == 0)
    endif
    npix = 12 * nside * nside
    ncap  = 2*nside*(nside-1) ! number of pixels in the north polar cap
    listir = -1
    nir = 0

    phi_low = MODULO(phi0 - dphi, twopi)
    phi_hi  = MODULO(phi0 + dphi, twopi)
    if (ABS(dphi-PI) < 1.0e-6_dp) take_all = .true.

    !     ------------ identifies ring number --------------
    if (iz >= nside .and. iz <= 3*nside) then ! equatorial region
       ir = iz - nside + 1  ! in {1, 2*nside + 1}
       ipix1 = ncap + 4*nside*(ir-1) !  lowest pixel number in the ring
       ipix2 = ipix1 + 4*nside - 1   ! highest pixel number in the ring
       kshift = MODULO(ir,2)
       nr = nside*4
    else
       if (iz < nside) then       !    north pole
          ir = iz
          ipix1 = 2*ir*(ir-1)        !  lowest pixel number in the ring
          ipix2 = ipix1 + 4*ir - 1   ! highest pixel number in the ring
       else                          !    south pole
          ir = 4*nside - iz
          ipix1 = npix - 2*ir*(ir+1) !  lowest pixel number in the ring
          ipix2 = ipix1 + 4*ir - 1   ! highest pixel number in the ring
       endif
       nr = ir*4
       kshift = 1
    endif

    !     ----------- constructs the pixel list --------------
    if (take_all) then
       nir    = ipix2 - ipix1 + 1
       if (do_ring) then
          listir(0:nir-1) = (/ (i, i=ipix1,ipix2) /)
       else
          call ring2nest(nside, ipix1, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
       return
    endif

    shift = kshift * 0.5_dp
    if (conservative) then
       ! conservative : include every intersected pixels, 
       ! even if pixel CENTER is not in the range [phi_low, phi_hi]
       ip_low = nint (nr * phi_low / TWOPI - shift)
       ip_hi  = nint (nr * phi_hi  / TWOPI - shift)
       ip_low = modulo (ip_low, nr) ! in {0,nr-1}
       ip_hi  = modulo (ip_hi , nr) ! in {0,nr-1}
    else
       ! strict : include only pixels whose CENTER is in [phi_low, phi_hi]
       ip_low = ceiling (nr * phi_low / TWOPI - shift) 
       ip_hi  = floor   (nr * phi_hi  / TWOPI - shift) 
       if ((ip_low - ip_hi == 1) .and. (dphi*nr < PI)) then
          ! the interval is so small (and away from pixel center)
          ! that no pixel is included in it
          nir = 0
          return
       endif
       ip_low = min(ip_low, nr-1)
       ip_hi  = max(ip_hi , 0   )
    endif
    !
    if (ip_low > ip_hi) to_top = .true.
    ip_low = ip_low + ipix1
    ip_hi  = ip_hi  + ipix1

    if (to_top) then
       nir1 = ipix2 - ip_low + 1
       nir2 = ip_hi - ipix1  + 1
       nir  = nir1 + nir2
       if (do_ring) then
          listir(0:nir1-1)   = (/ (i, i=ip_low, ipix2) /)
          listir(nir1:nir-1) = (/ (i, i=ipix1, ip_hi) /)
       else
          call ring2nest(nside, ip_low, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
    else
       nir = ip_hi - ip_low + 1
       if (do_ring) then
          listir(0:nir-1) = (/ (i, i=ip_low, ip_hi) /)
       else
          call ring2nest(nside, ip_low, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
    endif

    return
  end subroutine in_ring
  !=======================================================================
  subroutine intrs_intrv( d1, d2, di, ni)
    !=======================================================================
    ! computes the intersection di 
    ! of 2 intervals d1 (= [a1,b1]) and d2 (= [a2,b2])
    ! on the periodic domain ( = [A,B], where A and B are arbitrary)
    ! ni is the resulting number of intervals (0,1, or 2)
    !
    ! if a1<b1 then d1 = {x | a1 <= x <= b1}
    ! if a1>b1 then d1 = {x | a1 <= x <= B  U  A <= x <= b1}
    !=======================================================================
    real(kind=DP), dimension(1:), INTENT(IN)  :: d1, d2
    real(kind=DP), dimension(1:), INTENT(OUT) :: di
    integer(kind=I4B), INTENT(OUT) :: ni

    real(kind=DP), dimension(1:4) :: dk
    integer(kind=I4B) :: ik
    logical(kind=LGT) :: tr12, tr21, tr34, tr43, tr13, tr31, tr24, tr42, tr14, tr32
    !=======================================================================

    tr12 = (d1(1) < d1(2))
    tr21 = .NOT. tr12
    tr34 = (d2(1) < d2(2))
    tr43 = .NOT. tr34
    tr13 = (d1(1) < d2(1))
    tr31 = .NOT. tr13
    tr24 = (d1(2) < d2(2))
    tr42 = .NOT. tr24
    tr14 = (d1(1) < d2(2))
    tr32 = (d2(1) < d1(2))

    ik = 0
    dk(1:4) = -1.0e9_dp


    if ((tr34.AND.tr31.AND.tr14) .OR. (tr43.AND.(tr31.OR.tr14))) then
       ik = ik + 1
       dk(ik) = d1(1)  ! a1
    endif
    if ((tr12.AND.tr13.AND.tr32) .OR. (tr21.AND.(tr13.OR.tr32))) then
       ik = ik + 1
       dk(ik) = d2(1)  ! a2
    endif
    if ((tr34.AND.tr32.AND.tr24) .OR. (tr43.AND.(tr32.OR.tr24))) then
       ik = ik + 1
       dk(ik) = d1(2)  ! b1
    endif
    if ((tr12.AND.tr14.AND.tr42) .OR. (tr21.AND.(tr14.OR.tr42))) then
       ik = ik + 1
       dk(ik) =  d2(2)  ! b2
    endif

    di(1:4) = 0.0_dp
    select case (ik)
    case (0)
       ni = 0
    case (2)
       ni = 1
       di(1:2) = (/ dk(1), dk(2) /) ! [a1,b1] or [a1,b2] or [a2,b1] or [a2,b2]
    case (4)
       ni = 2
       di(1:4) = (/ dk(1), dk(4), dk(2), dk(3) /) ! [a1,b2] U [a2,b1]
    case default
       print*,"error in intrs_intrv", ik
       print*,dk
       print*,d1,d2
    end select

    return
  end subroutine intrs_intrv

  !=======================================================================
  subroutine pix2ang_ring(nside, ipix, theta, phi)
    !=======================================================================
    !     renders theta and phi coordinates of the nominal pixel center
    !     for the pixel number ipix (RING scheme)  
    !     given the map resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: ipix, nside
    REAL(KIND=DP), INTENT(OUT) ::  theta, phi

    INTEGER(KIND=I4B) ::  nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
    REAL(KIND=DP) ::  fodd, hip, fihip
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    npix = 12*nside**2       ! total number of points
    if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

    ipix1 = ipix + 1 ! in {1, npix}
    nl2 = 2*nside
    ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1

    if (ipix1 <= ncap) then ! North Polar cap -------------

       hip   = ipix1*0.5_dp
       fihip = AINT ( hip , kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipix1 - 2*iring*(iring - 1)

       theta = ACOS( 1.0_dp - iring**2 / (3.0_dp*nside**2) )
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

    elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

       ip    = ipix1 - ncap - 1
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (1.5_dp*nside) )
       phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix1 + 1
       hip   = ip*0.5_dp
       fihip = AINT ( hip , kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

       theta = ACOS( -1.0_dp + iring**2 / (3.0_dp*nside**2) )
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

    endif

    return
  end subroutine pix2ang_ring ! pix2ang_ring
!   !=======================================================================
!   subroutine pix2vec0_ring(nside, ipix, vector)
!     !=======================================================================
!     !     renders vector (x,y,z) coordinates of the nominal pixel center
!     !     for the pixel number ipix (RING scheme)  
!     !     given the map resolution parameter nside
!     !=======================================================================
!     INTEGER(KIND=I4B), INTENT(IN) :: ipix, nside
!     REAL(KIND=DP), INTENT(OUT),dimension(1:) :: vector

!     INTEGER(KIND=I4B) :: nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
!     REAL(KIND=DP) ::  fact1, fact2, fodd, hip, fihip, z, sth, phi
!     !-----------------------------------------------------------------------
!     if (nside<1 .or. nside>ns_max) stop "nside out of range"
!     npix = 12*nside**2       ! total number of points
!     if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

!     ipix1 = ipix + 1 ! in {1, npix}
!     nl2 = 2*nside
!     nl4 = 4*nside
!     ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1
!     fact1 = 1.5_dp*nside
!     fact2 = 3.0_dp*nside**2

!     if (ipix1 <= ncap) then ! North Polar cap -------------

!        hip   = ipix1/2.0_dp
!        fihip = AINT ( hip ,kind=DP)
!        iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
!        iphi  = ipix1 - 2*iring*(iring - 1)

!        z =  1.0_dp - iring**2 / fact2 
!        phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

!     elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

!        ip    = ipix1 - ncap - 1
!        iring = INT( ip / nl4 ) + nside ! counted from North pole
!        iphi  = MODULO(ip,nl4) + 1

!        fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
!        z = (nl2 - iring) / fact1 
!        phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

!     else ! South Polar cap -----------------------------------

!        ip    = npix - ipix1 + 1
!        hip   = ip/2.0_dp
!        fihip = AINT ( hip ,kind=DP)
!        iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
!        iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

!        z = -1.0_dp + iring**2 / fact2 
!        phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

!     endif

!     sth = SQRT((1.0_dp-z)*(1.0_dp+z))
!     vector(1) = sth * COS(phi)
!     vector(2) = sth * SIN(phi)
!     vector(3) = z

!     return
!   end subroutine pix2vec0_ring
  !=======================================================================
  subroutine pix2vec_ring(nside, ipix, vector, vertex)
    !=======================================================================
    !     renders vector (x,y,z) coordinates of the nominal pixel center
    !     for the pixel number ipix (RING scheme)  
    !     given the map resolution parameter nside
    !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    !     in the order N,W,S,E
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN)                             :: ipix, nside
    REAL(KIND=DP),     INTENT(OUT),dimension(1:)              :: vector
    REAL(KIND=DP),     INTENT(OUT),dimension(1:,1:), optional :: vertex

    INTEGER(KIND=I4B) :: nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
    REAL(KIND=DP) ::  fact1, fact2, fodd, hip, fihip, z, sth, phi

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    npix = 12*nside**2       ! total number of points
    if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

    ipix1 = ipix + 1 ! in {1, npix}
    nl2 = 2*nside
    nl4 = 4*nside
    ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1
    fact1 = 1.5_dp*nside
    fact2 = 3.0_dp*nside**2

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          stop " pix2vec_ring : vertex array has wrong size "
       endif
    endif

    phi_nv = 0.0_dp
    phi_sv = 0.0_dp
    if (ipix1 <= ncap) then ! North Polar cap -------------

       hip   = ipix1/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipix1 - 2*iring*(iring - 1)

       z =  1.0_dp - iring**2 / fact2 
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = 1.0_dp - (iring-1)**2 / fact2 
          z_sv = 1.0_dp - (iring+1)**2 / fact2 
          iphi_mod = MODULO(iphi-1, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi-1) / iring      ! in {0,1,2,3}
          if (iring > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
          phi_sv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
       endif
 

    elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

       ip    = ipix1 - ncap - 1
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
       z = (nl2 - iring) / fact1 
       phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*nside)   ! half pixel width
          phi_nv = phi
          phi_sv = phi
          z_nv = (nl2 - iring +1) / fact1 
          z_sv = (nl2 - iring -1) / fact1 
          if (iring == nside) then ! northern transition
             z_nv = 1.0_dp - (nside-1)**2 / fact2 
             iphi_mod = MODULO(iphi-1, nside) ! in {0,1,... nside-1}
             iphi_rat = (iphi-1) / nside      ! in {0,1,2,3}
             if (nside > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(nside-1,kind=dp))
          elseif (iring == 3*nside) then ! southern transition
             z_sv = -1.0_dp + (nside-1)**2 / fact2 
             iphi_mod = MODULO(iphi-1, nside) ! in {0,1,... iring-1}
             iphi_rat = (iphi-1) / nside      ! in {0,1,2,3}
             if (nside > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(nside-1,kind=dp))
          endif
       endif

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix1 + 1
       hip   = ip/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

       z = -1.0_dp + iring**2 / fact2 
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = -1.0_dp + (iring+1)**2 / fact2 
          z_sv = -1.0_dp + (iring-1)**2 / fact2 
          iphi_mod = MODULO(iphi-1, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi-1) / iring      ! in {0,1,2,3}
          phi_nv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
          if (iring > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
       endif

    endif

    ! pixel center
    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z

    if (do_vertex) then
       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north vertex
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       vertex(1,1) = sth_nv * COS(phi_nv)
       vertex(2,1) = sth_nv * SIN(phi_nv)
       vertex(3,1) = z_nv

       ! south vertex
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       vertex(1,3) = sth_sv * COS(phi_sv)
       vertex(2,3) = sth_sv * SIN(phi_sv)
       vertex(3,3) = z_sv
    endif


    return
  end subroutine pix2vec_ring
  !=======================================================================
  subroutine ang2pix_ring(nside, theta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map 
    !     resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN) ::  theta, phi

    INTEGER(KIND=I4B) ::  nl4, jp, jm
    REAL(KIND=DP) ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(KIND=I4B) ::  ir, ip, kshift

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2PIX_RING: theta : ",theta," is out of range [0, Pi]"
       stop
    endif

    z = COS(theta)
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)


    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(.5_dp+tt)
       temp2 = nside*.75_dp*z
       jp = int(temp1-temp2) ! index of  ascending edge line 
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 1 - modulo(ir,2) ! kshift=1 if ir even, 0 otherwise

       nl4 = 4*nside
       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1) + nl4*(ir-1) + ip 

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT(tp          * tmp ) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir )     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 12*nside**2 - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix_ring

  !=======================================================================
  subroutine vec2pix_ring(nside, vector, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinate vector (=x,y,z), given the map 
    !     resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN), dimension(1:) :: vector

    INTEGER(KIND=I4B) :: nl2, nl4, ncap, npix, jp, jm, ipix1
    REAL(KIND=DP) ::  z, za, tt, tp, tmp, dnorm, phi
    INTEGER(KIND=I4B) :: ir, ip, kshift

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"

    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z = vector(3) / dnorm
    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[
    tt = phi / halfpi   ! in [0,4)

    nl2 = 2*nside
    nl4 = 4*nside
    ncap  = nl2*(nside-1) ! number of pixels in the north polar cap
    npix  = 12*nside**2

    if ( za <= twothird ) then ! Equatorial region ------------------

       jp = INT(nside*(0.5_dp + tt - z*0.75_dp)) ! index of  ascending edge line 
       jm = INT(nside*(0.5_dp + tt + z*0.75_dp)) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 0
       if (MODULO(ir,2) == 0) kshift = 1 ! kshift=1 if ir even, 0 otherwise

       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1 ! in {1,4n}
       if (ip > nl4) ip = ip - nl4

       ipix1 = ncap + nl4*(ir-1) + ip 

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT( nside * tp          * tmp ) ! increasing edge line index
       jm = INT( nside * (1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir ) + 1 ! in {1,4*ir}
       if (ip > 4*ir) ip = ip - 4*ir

       ipix1 = 2*ir*(ir-1) + ip
       if (z <= 0.0_dp) then
          ipix1 = npix - 2*ir*(ir+1) + ip
       endif

    endif

    ipix = ipix1 - 1 ! in {0, npix-1}

    return
  end subroutine vec2pix_ring
  !=======================================================================
  subroutine pix2ang_nest(nside, ipix, theta, phi)
    !=======================================================================
    !     renders theta and phi coordinates of the nominal pixel center
    !     for the pixel number ipix (NESTED scheme)  
    !     given the map resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT) :: theta, phi

    INTEGER(KIND=I4B) :: npix, npface, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4
    REAL(KIND=DP) :: z, fn, fact1, fact2

    INTEGER(KIND=I4B) :: ix, iy, face_num
!     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    npix = 12 * nside**2
    if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    fn = real(nside,kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)
    nl4   = 4*nside

    !     finds the face, and the number in the face
    npface = nside**2

    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    z  = (2*nside-jr)*fact2
    kshift = MODULO(jr - nside, 2)
    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*nr*fact1
       kshift = 0
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_dp + nr*nr*fact1
       kshift = 0
    endif
    theta = ACOS(z)

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

    return
  end subroutine pix2ang_nest
!   !=======================================================================
!   subroutine pix2vec_nest(nside, ipix, vector)
!     !=======================================================================
!     !     renders vector (x,y,z) coordinates of the nominal pixel center
!     !     for the pixel number ipix (NESTED scheme)  
!     !     given the map resolution parameter nside
!     !=======================================================================
!     INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
!     REAL(KIND=DP), INTENT(OUT), dimension(1:) :: vector

!     INTEGER(KIND=I4B) :: npix, npface, &
!          &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
!          &     jrt, jr, nr, jpt, jp, kshift, nl4
!     REAL(KIND=DP) :: z, fn, fact1, fact2, sth, phi

!     INTEGER(KIND=I4B) ::  ix, iy, face_num
! !     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

!     ! coordinate of the lowest corner of each face
!     INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
!     INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
!     !-----------------------------------------------------------------------
!     if (nside<1 .or. nside>ns_max) stop "nside out of range"
!     npix = 12 * nside**2
!     if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

!     !     initiates the array for the pixel number -> (x,y) mapping
!     if (pix2x(1023) <= 0) call mk_pix2xy()

!     fn = real(nside,kind=dp)
!     fact1 = 1.0_dp/(3.0_dp*fn*fn)
!     fact2 = 2.0_dp/(3.0_dp*fn)
!     nl4   = 4*nside

!     !     finds the face, and the number in the face
!     npface = nside**2

!     face_num = ipix/npface  ! face number in {0,11}
!     ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

!     !     finds the x,y on the face (starting from the lowest corner)
!     !     from the pixel number
!     ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
!     ip_trunc =   ipf/1024        ! truncation of the last 10 bits
!     ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
!     ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

!     ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
!     iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

!     !     transforms this in (horizontal, vertical) coordinates
!     jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
!     jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

!     !     computes the z coordinate on the sphere
!     jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

!     nr = nside                  ! equatorial region (the most frequent)
!     z  = (2*nside-jr)*fact2
!     kshift = MODULO(jr - nside, 2)
!     if (jr < nside) then     ! north pole region
!        nr = jr
!        z = 1.0_dp - nr*nr*fact1
!        kshift = 0
!     else if (jr > 3*nside) then ! south pole region
!        nr = nl4 - jr
!        z = - 1.0_dp + nr*nr*fact1
!        kshift = 0
!     endif
!     !ccc      theta = ACOS(z)

!     !     computes the phi coordinate on the sphere, in [0,2Pi]
!     jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
!     if (jp > nl4) jp = jp - nl4
!     if (jp < 1)   jp = jp + nl4

!     phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

!     sth = SQRT((1.0_dp-z)*(1.0_dp+z))
!     vector(1) = sth * COS(phi)
!     vector(2) = sth * SIN(phi)
!     vector(3) = z

!     return
!   end subroutine pix2vec_nest
  !=======================================================================
  subroutine pix2vec_nest(nside, ipix, vector, vertex)
    !=======================================================================
    !     renders vector (x,y,z) coordinates of the nominal pixel center
    !     for the pixel number ipix (NESTED scheme)  
    !     given the map resolution parameter nside
    !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    !     in the order N,W,S,E
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT), dimension(1:) :: vector
    REAL(KIND=DP),     INTENT(OUT),dimension(1:,1:), optional :: vertex

    INTEGER(KIND=I4B) :: npix, npface, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4
    REAL(KIND=DP) :: z, fn, fact1, fact2, sth, phi

    INTEGER(KIND=I4B) ::  ix, iy, face_num
!     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    npix = 12 * nside**2
    if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    fn = real(nside,kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)
    nl4   = 4*nside

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          stop " pix2vec_ring : vertex array has wrong size "
       endif
    endif

    !     finds the face, and the number in the face
    npface = nside**2

    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    z  = (2*nside-jr)*fact2
    kshift = MODULO(jr - nside, 2)
    if (do_vertex) then
       z_nv = (2*nside-jr+1)*fact2
       z_sv = (2*nside-jr-1)*fact2
       if (jr == nside) then ! northern transition
          z_nv =  1.0_dp - (nside-1)**2 * fact1
       elseif (jr == 3*nside) then  ! southern transition
          z_sv = -1.0_dp + (nside-1)**2 * fact1
       endif
    endif
    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = 1.0_dp - (nr-1)**2*fact1
          z_sv = 1.0_dp - (nr+1)**2*fact1
       endif
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_dp + nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = - 1.0_dp + (nr+1)**2*fact1
          z_sv = - 1.0_dp + (nr-1)**2*fact1
       endif
    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z

    if (do_vertex) then
       phi_nv = phi
       phi_sv = phi

       phi_up = 0.0_dp
       iphi_mod = MODULO(jp-1, nr) ! in {0,1,... nr-1}
       iphi_rat = (jp-1) / nr      ! in {0,1,2,3}
       if (nr > 1) phi_up = HALFPI * (iphi_rat +  iphi_mod   /real(nr-1,kind=dp))
       phi_dn             = HALFPI * (iphi_rat + (iphi_mod+1)/real(nr+1,kind=dp))
       if (jr < nside) then            ! North polar cap
          phi_nv = phi_up
          phi_sv = phi_dn
       else if (jr > 3*nside) then     ! South polar cap
          phi_nv = phi_dn
          phi_sv = phi_up          
       else if (jr == nside) then      ! North transition
          phi_nv = phi_up
       else if (jr == 3*nside) then    ! South transition
          phi_sv = phi_up
       endif

       hdelta_phi = PI / (4.0_dp*nr)

       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north vertex
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       vertex(1,1) = sth_nv * COS(phi_nv)
       vertex(2,1) = sth_nv * SIN(phi_nv)
       vertex(3,1) = z_nv

       ! south vertex
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       vertex(1,3) = sth_sv * COS(phi_sv)
       vertex(2,3) = sth_sv * SIN(phi_sv)
       vertex(3,3) = z_sv
    endif

    return
  end subroutine pix2vec_nest
  !=======================================================================
  subroutine ang2pix_nest(nside, theta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (NESTED scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map 
    !     resolution parametr nside
    !
    !     the computation is made to the highest resolution available (nside=8192)
    !     and then degraded to that required (by integer division)
    !     this doesn't cost more, and it makes sure 
    !     that the treatement of round-off will be consistent 
    !     for every resolution
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN) ::  theta, phi

    REAL(KIND=DP) ::  z, za, tt, tp, tmp
    INTEGER(KIND=I4B) :: jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2PIX_NEST: theta : ",theta," is out of range [0,Pi]"
       stop
    endif
    if (x2pix(128) <= 0) call mk_xy2pix()

    z  = COS(theta)
    za = ABS(z)
    tt = MODULO(phi, twopi) / halfpi  ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(ns_max*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(ns_max*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       ifp = jp / ns_max  ! in {0,4}
       ifm = jm / ns_max
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = MODULO(ifp,4)
       else                            ! (half-)faces 8 to 11
          face_num = MODULO(ifm,4) + 8
       endif

       ix = MODULO(jm, ns_max)
       iy = ns_max - MODULO(jp, ns_max) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( ns_max * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( ns_max * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(ns_max-1, jp) ! for points too close to the boundary
       jm = MIN(ns_max-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = ns_max - jm - 1
          iy = ns_max - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

       !         print*,z,face_num,ix,iy
    endif

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    ipf = ipf / ( ns_max/nside ) **2  ! in {0, nside**2 - 1}

    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
  end subroutine ang2pix_nest
  !=======================================================================
  subroutine vec2pix_nest(nside, vector, ipix)
    !=======================================================================
    !     renders the pixel number ipix (NESTED scheme) for a pixel which contains
    !     a point on a sphere at coordinate vector (=x,y,z), given the map 
    !     resolution parameter nside
    !
    !     the computation is made to the highest resolution available (nside=8192)
    !     and then degraded to that required (by integer division)
    !     this doesn't cost more, and it makes sure 
    !     that the treatement of round-off will be consistent 
    !     for every resolution
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN), dimension(1:) ::  vector

    REAL(KIND=DP) ::  z, za, tt, tp, tmp, dnorm, phi
    INTEGER(KIND=I4B) ::  jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    if (x2pix(128) <= 0) call mk_xy2pix()

    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z = vector(3) / dnorm
    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0)    phi = phi + twopi ! phi in [0,2pi[
    tt = phi / halfpi ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(ns_max*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(ns_max*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       ifp = jp / ns_max  ! in {0,4}
       ifm = jm / ns_max
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = MODULO(ifp,4)
       else                            ! (half-)faces 8 to 11
          face_num = MODULO(ifm,4) + 8
       endif

       ix = MODULO(jm, ns_max)
       iy = ns_max - MODULO(jp, ns_max) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( ns_max * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( ns_max * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(ns_max-1, jp) ! for points too close to the boundary
       jm = MIN(ns_max-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = ns_max - jm - 1
          iy = ns_max - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

       !         print*,z,face_num,ix,iy
    endif

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    ipf = ipf / ( ns_max/nside ) **2  ! in {0, nside**2 - 1}

    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
  end subroutine vec2pix_nest
  !=======================================================================
  subroutine convert_nest2ring(nside, map)
    !=======================================================================
    !     makes the conversion NEST to RING
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    REAL(KIND=SP), DIMENSION(0:), INTENT(INOUT) ::  map

    INTEGER(KIND=I4B) :: npix
    REAL(KIND=SP), DIMENSION(:), ALLOCATABLE :: map_tmp
    INTEGER(KIND=I4B) :: ipn, ipr
    !=======================================================================

    npix = 12*nside*nside
    ALLOCATE(map_tmp(0:npix-1))

    do ipn = 0, npix-1 
       call nest2ring(nside, ipn, ipr)
       map_tmp(ipr) = map(ipn)
    enddo

    do ipr = 0, npix - 1
       map(ipr) = map_tmp(ipr)
    enddo

    DEALLOCATE(map_tmp)

    return
  end subroutine convert_nest2ring
  !=======================================================================
  subroutine convert_ring2nest(nside, map)
    !=======================================================================
    !     makes the conversion NEST to RING
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    REAL(KIND=SP), DIMENSION(0:), INTENT(INOUT) ::  map

    INTEGER(KIND=I4B) :: npix
    REAL(KIND=SP), DIMENSION(:), ALLOCATABLE :: map_tmp
    INTEGER(KIND=I4B) :: ipn, ipr
    !=======================================================================

    npix = 12*nside*nside
    ALLOCATE(map_tmp(0:npix-1))

    do ipr = 0, npix-1 
       call ring2nest(nside, ipr, ipn)
       map_tmp(ipn) = map(ipr)
    enddo

    do ipn = 0, npix - 1
       map(ipn) = map_tmp(ipn)
    enddo

    DEALLOCATE(map_tmp)

    return
  end subroutine convert_ring2nest

  !====================================================================
  ! The following two routines convert both real
  ! and integer arrays between the
  ! NESTED and RING schemes. 
  !
  ! Author: Benjamin D. Wandelt October 1997
  !

  !====================================================================
  subroutine convert_inplace_real(subcall,map)
    !====================================================================
    !     Converts a single precision real map from RING to NESTED and vice versa 
    !     in place, ie. without allocating a temporary map (only a logical array). 
    !     This routine is more general, but slower than convert_nest2ring.
    !
    !     This is a wrapper for the functions "ring2nest" and 
    !     "nest2ring". Their names are supplied in the "subcall" 
    !     argument.
    !     
    !
    ! Benjamin D. Wandelt October 1997
    ! Added to pix_tools for version 1.00 in March 1999.
    !====================================================================

    external  subcall ! required by some compilers (DEC, Portland, ...)
    real(kind=sp), pointer, dimension(:) :: map

    integer(kind=i4b), parameter :: ns_max=8192
    integer(kind=i4b) :: npix,nside
    logical(kind=lgt), dimension(:), allocatable::check
    integer(kind=i4b) :: ilast,i1,i2
    real(kind=sp) :: pixbuf1,pixbuf2
    !------------------------------------------------------------------
    npix=size(map)
    nside=sqrt(real(npix/12))
    if(nside>ns_max) stop "convert_inplace_real: map too big"

    print*, "Convert: Converting map pixelisation scheme"
    allocate(check(0:npix-1))
    check=.false.

    ilast=0                   !start at first pixel
    do
       pixbuf2=map(ilast)      !initialise 
       i1=ilast
       call subcall(nside,i1,i2)
       do 
          if (check(i2)) exit
          pixbuf1=map(i2)
          map(i2)=pixbuf2
          check(i2)=.true.
          pixbuf2=pixbuf1
          i1=i2
          call subcall(nside,i1,i2)
       enddo
       do 
          if (.not. (check(ilast).and.ilast<npix-1)) exit ! npix-1 or npix
          ilast=ilast+1
       enddo
       if(ilast==npix-1) exit ! npix-1 or npix
    enddo
    deallocate(check)
    return
  end subroutine convert_inplace_real

  !====================================================================
  subroutine convert_inplace_int(subcall,map)
    !====================================================================
    !     Converts a 4 byte integer map from RING to NESTED and vice versa 
    !     in place, ie. without allocating a temporary map. This routine is more general,
    !     but slower than convert_nest2ring.

    !     This is a wrapper for the toolbox functions "ring2nest" and 
    !     "nest2ring". Their names are supplied in the "subcall" 
    !     argument.
    !
    ! Benjamin D. Wandelt October 1997
    ! Added to pix_tools for version 1.00 in March 1999.
    !====================================================================

    external  subcall ! required by some compiler (DEC, Portland, ...)
    integer(kind=i4b), dimension(:), pointer:: map

    integer(kind=i4b), parameter :: ns_max=8192
    integer(kind=i4b) :: npix,nside
    logical(kind=lgt), dimension(:), allocatable::check 
    integer(kind=i4b) :: ilast,i1,i2
    real(kind=sp) :: pixbuf1,pixbuf2
    !------------------------------------------------------------------
    npix=size(map)
    nside=sqrt(real(npix/12))
    if(nside>ns_max) stop "convert_inplace_int: map too big"

    print*, "Convert: Converting map pixelisation scheme"
    allocate(check(0:npix-1))
    check=.false.

    ilast=0                   !start at first pixel
    do
       pixbuf2=map(ilast)      !initialise 
       i1=ilast
       call subcall(nside,i1,i2)
       do 
          if (check(i2)) exit
          pixbuf1=map(i2)
          map(i2)=pixbuf2
          check(i2)=.true.
          pixbuf2=pixbuf1
          i1=i2
          call subcall(nside,i1,i2)
       enddo
       do 
          if (.not. (check(ilast).and.ilast<npix-1)) exit ! npix-1 or npix
          ilast=ilast+1
       enddo
       if(ilast==npix-1) exit ! npix-1 or npix
    enddo
    deallocate(check)
    return
  end subroutine convert_inplace_int

  !=======================================================================
  subroutine nest2ring(nside, ipnest, ipring)
    !=======================================================================
    !     performs conversion from NESTED to RING pixel number
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) ::  nside, ipnest
    INTEGER(KIND=I4B), INTENT(OUT) :: ipring
    INTEGER(KIND=I4B) ::  npix, npface, face_num, ncap, n_before, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    npix = 12 * nside**2
    if (ipnest<0 .or. ipnest>npix-1) stop "ipnest out of range"

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    ncap  = 2*nside*(nside-1) ! number of points in the North Polar cap
    nl4   = 4*nside

    !     finds the face, and the number in the face
    npface = nside**2
    !ccccc      ip = ipnest - 1         ! in {0,npix-1}

    face_num = ipnest/npface  ! face number in {0,11}
    ipf = MODULO(ipnest,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    n_before = ncap + nl4 * (jr - nside)
    kshift = MODULO(jr - nside, 2)
    if (jr < nside) then     ! north pole region
       nr = jr
       n_before = 2 * nr * (nr - 1)
       kshift = 0
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       n_before = npix - 2 * (nr + 1) * nr
       kshift = 0
    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}

    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    ipring = n_before + jp - 1 ! in {0, npix-1}
    return
  end subroutine nest2ring
  !=======================================================================
  subroutine ring2nest(nside, ipring, ipnest)
    !=======================================================================
    !     performs conversion from RING to NESTED pixel number
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipring
    INTEGER(KIND=I4B), INTENT(OUT) :: ipnest

    REAL(KIND=DP) :: fihip, hip
    INTEGER(KIND=I4B) :: npix, nl2, nl4, ncap, ip, iphi, ipt, ipring1, &
         &     kshift, face_num, nr, &
         &     irn, ire, irm, irs, irt, ifm , ifp, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    npix = 12*nside**2      ! total number of points
    if (ipring <0 .or. ipring>npix-1) stop "ipring out of range"
    if (x2pix(128) <= 0) call mk_xy2pix()

    nl2 = 2*nside
    nl4 = 4*nside
    ncap = nl2*(nside-1) ! points in each polar cap, =0 for nside =1
    ipring1 = ipring + 1

    !     finds the ring number, the position of the ring and the face number
    if (ipring1 <= ncap) then ! north polar cap

       hip   = ipring1/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       irn   = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipring1 - 2*irn*(irn - 1)

       kshift = 0
       nr = irn                  ! 1/4 of the number of points on the current ring
       face_num = (iphi-1) / irn ! in {0,3}

    elseif (ipring1 <= nl2*(5*nside+1)) then ! equatorial region

       ip    = ipring1 - ncap - 1
       irn   = INT( ip / nl4 ) + nside               ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       kshift  = MODULO(irn+nside,2)  ! 1 if irn+nside is odd, 0 otherwise
       nr = nside
       ire =  irn - nside + 1 ! in {1, 2*nside +1}
       irm =  nl2 + 2 - ire
       ifm = (iphi - ire/2 + nside -1) / nside ! face boundary
       ifp = (iphi - irm/2 + nside -1) / nside
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp + 1 == ifm) then ! (half-)faces 0 to 3
          face_num = ifp
       else if (ifp - 1 == ifm) then ! (half-)faces 8 to 11
          face_num = ifp + 7
       endif

    else ! south polar cap

       ip    = npix - ipring1 + 1
       hip   = ip/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       irs   = INT( SQRT( hip - SQRT(fihip) ) ) + 1  ! counted from South pole
       iphi  = 4*irs + 1 - (ip - 2*irs*(irs-1))

       kshift = 0
       nr = irs
       irn   = nl4 - irs
       face_num = (iphi-1) / irs + 8 ! in {8,11}

    endif

    !     finds the (x,y) on the face
    irt =   irn  - jrll(face_num+1)*nside + 1       ! in {-nside+1,0}
    ipt = 2*iphi - jpll(face_num+1)*nr - kshift - 1 ! in {-nside+1,nside-1}
    if (ipt >= nl2) ipt = ipt - 8*nside ! for the face #4

    ix =  (ipt - irt ) / 2
    iy = -(ipt + irt ) / 2

    ix_low = MODULO(ix,128)
    ix_hi  = ix/128
    iy_low = MODULO(iy,128)
    iy_hi  = iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))        ! in {0, nside**2 - 1}


    ipnest = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
  end subroutine ring2nest
  !=======================================================================
  subroutine xy2pix_nest(nside, ix, iy, face_num, ipix)
    !=======================================================================
    !     gives the pixel number ipix (NESTED) 
    !     corresponding to ix, iy and face_num
    !
    !     Benjamin D. Wandelt 13/10/97
    !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) ::  nside, ix, iy, face_num    
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    INTEGER(KIND=I4B) ::  ix_low, ix_hi, iy_low, iy_hi, ipf

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    if (ix<0 .or. ix>(nside-1))     stop "ix out of range"
    if (iy<0 .or. iy>(nside-1))     stop "iy out of range"
    if (x2pix(128) <= 0) call mk_xy2pix()

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}
    return
  end subroutine xy2pix_nest
  !=======================================================================
  subroutine pix2xy_nest(nside, ipf, ix, iy)
    !=======================================================================
    !     gives the x, y coords in a face from pixel number within the face (NESTED) 
    !
    !     Benjamin D. Wandelt 13/10/97
    !     
    !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipf
    INTEGER(KIND=I4B), INTENT(OUT) :: ix, iy    

    INTEGER(KIND=I4B) ::  ip_low, ip_trunc, ip_med, ip_hi

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    if (ipf <0 .or. ipf>nside*nside-1) &
         &     stop "ipix out of range"
    if (pix2x(1023) <= 0) call mk_pix2xy()

    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    return
  end subroutine pix2xy_nest
  !=======================================================================
  subroutine mk_pix2xy()
    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels
    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================
    INTEGER(KIND=I4B) ::  kpix, jpix, ix, iy, ip, id

    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do 
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31
    enddo

    return
  end subroutine mk_pix2xy
  !=======================================================================
  subroutine mk_xy2pix()
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    INTEGER(KIND=I4B):: k,ip,i,j,id
    !=======================================================================

    do i = 1,128           !for converting x,y into
       j  = i-1            !pixel numbers
       k  = 0
       ip = 1

       do
          if (j==0) then
             x2pix(i) = k
             y2pix(i) = 2*k
             exit
          else
             id = MODULO(J,2)
             j  = j/2
             k  = ip*id+k
             ip = ip*4
          endif
       enddo

    enddo

    RETURN
  END subroutine mk_xy2pix
  !=======================================================================

  ! The following is a routine which finds the 7 or 8 neighbours of 
  ! any pixel in the nested scheme of the HEALPIX pixelisation.
  !====================================================================
  subroutine neighbours_nest(nside,ipix,n,nneigh)
    !====================================================================
    !   Returns list n(8) of neighbours of pixel ipix (in NESTED scheme)
    !   the neighbours are ordered in the following way:
    !   First pixel is the one to the south (the one west of the south
    ! direction is taken
    ! for the pixels which don't have a southern neighbour). From
    ! then on the neighbours are ordered in the clockwise direction
    ! about the pixel with number ipix.
    !     
    !   nneigh is the number of neighbours (mostly 8, 8 pixels have 7 neighbours)
    !
    !   Benjamin D. Wandelt October 1997
    !   Added to pix_tools in March 1999
    !====================================================================
    use bit_manipulation
    integer(kind=i4b), intent(in)::nside, ipix
    integer(kind=i4b), intent(out), dimension(1:):: n
    integer(kind=i4b), intent(out):: nneigh

    integer(kind=i4b) :: npix,ipf,ipo,ix,ixm,ixp,iy,iym,iyp,ixo,iyo
    integer(kind=i4b) :: face_num,other_face
    integer(kind=i4b) :: ia,ib,ibp,ibm,ib2,icase,nsidesq
    integer(kind=i4b) :: local_magic1,local_magic2

!     integer(kind=i4b), intrinsic :: IAND

    !--------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    nsidesq=nside*nside
    npix = 12*nsidesq       ! total number of points
    if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

    !     initiates array for (x,y)-> pixel number -> (x,y) mapping
    if (x2pix(128) <= 0) call mk_xy2pix()

    local_magic1=(nsidesq-1)/3
    local_magic2=2*local_magic1
    face_num=ipix/nsidesq

    ipf=modulo(ipix,nsidesq)   !Pixel number in face

    call pix2xy_nest(nside,ipf,ix,iy)
    ixm=ix-1
    ixp=ix+1
    iym=iy-1
    iyp=iy+1

    nneigh=8                  !Except in special cases below

    !     Exclude corners
    if(ipf==local_magic2)     then !WestCorner
       icase=5             
       goto 100
    endif
    if(ipf==(nsidesq-1)) then !NorthCorner
       icase=6             
       goto 100
    endif
    if(ipf==0)           then !SouthCorner
       icase=7                
       goto 100
    endif
    if(ipf==local_magic1)     then !EastCorner
       icase=8      
       goto 100
    endif

    !     Detect edges
    if(IAND(ipf,local_magic1)==local_magic1) then !NorthEast
       icase=1 
       goto 100
    endif
    if(IAND(ipf,local_magic1)==0)      then !SouthWest
       icase=2    
       goto 100
    endif
    if(IAND(ipf,local_magic2)==local_magic2) then !NorthWest
       icase=3                 
       goto 100
    endif
    if(IAND(ipf,local_magic2)==0)      then !SouthEast
       icase=4                   
       goto 100
    endif

    !     Inside a face
    call xy2pix_nest(nside, ixm, iym, face_num, n(1))
    call xy2pix_nest(nside, ixm, iy , face_num, n(2))   
    call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
    call xy2pix_nest(nside, ix , iyp, face_num, n(4))
    call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
    call xy2pix_nest(nside, ixp, iy , face_num, n(6))
    call xy2pix_nest(nside, ixp, iym, face_num, n(7))
    call xy2pix_nest(nside, ix , iym, face_num, n(8))
    return

100 continue

    ia= face_num/4            !in {0,2}
    ib= modulo(face_num,4)       !in {0,3}
    ibp=modulo(ib+1,4)
    ibm=modulo(ib+4-1,4)
    ib2=modulo(ib+2,4)

    if(ia==0) then          !North Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))         
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1 , iyo, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=4+ib
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))   
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=4+ib
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=0+ibm
          n(3)=other_face*nsidesq+local_magic1
          n(4)=n(3)+2
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          n(8)=ipix-2
          other_face=0+ibm
          n(4)=other_face*nsidesq+nsidesq-1
          n(3)=n(4)-2
          other_face=0+ib2
          n(5)=other_face*nsidesq+nsidesq-1
          other_face=0+ibp
          n(6)=other_face*nsidesq+nsidesq-1
          n(7)=n(6)-1
       case(7)              !South corner
          other_face=8+ib
          n(1)=other_face*nsidesq+nsidesq-1
          other_face=4+ib
          n(2)=other_face*nsidesq+local_magic1
          n(3)=n(2)+2
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=4+ibp
          n(8)=other_face*nsidesq+local_magic2
          n(7)=n(8)+1
       case(8)              !East corner
          nneigh=7
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(5)=n(6)+1
          other_face=4+ibp
          n(7)=other_face*nsidesq+nsidesq-1
          n(1)=n(7)-1
       end select ! north

    elseif(ia==1) then      !Equatorial region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ib
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))         
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ib
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))   
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          other_face=8+ibm
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=4+ibm
          n(3)=other_face*nsidesq+local_magic1
          other_face=0+ibm
          n(4)=other_face*nsidesq
          n(5)=n(4)+1
          n(6)=ipix+1
          n(7)=ipix-1
          n(8)=ipix-2         
       case(6)              !North corner
          nneigh=7
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=0+ibm
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq+local_magic2
          n(6)=n(5)-2
          n(7)=ipix-2
       case(7)              !South corner
          nneigh=7
          other_face=8+ibm
          n(1)=other_face*nsidesq+local_magic1
          n(2)=n(1)+2
          n(3)=ipix+2
          n(4)=ipix+3
          n(5)=ipix+1
          other_face=8+ib
          n(7)=other_face*nsidesq+local_magic2
          n(6)=n(7)+1
       case(8)              !East corner
          other_face=8+ib
          n(8)=other_face*nsidesq+nsidesq-1
          n(1)=n(8)-1
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ib
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
          other_face=4+ibp
          n(7)=other_face*nsidesq+local_magic2
       end select ! equator
    else                    !South Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))         
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)        !W-E flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=4+ib
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))   
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=8+ibm
          n(2)=other_face*nsidesq+local_magic1
          n(1)=n(2)-1
          other_face=4+ib
          n(3)=other_face*nsidesq
          n(4)=n(3)+1
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=4+ib
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq
          other_face=4+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(7)=n(6)-2
          n(8)=ipix-2
       case(7)              !South corner
          other_face=8+ib2
          n(1)=other_face*nsidesq
          other_face=8+ibm
          n(2)=other_face*nsidesq
          n(3)=n(2)+1
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=8+ibp
          n(8)=other_face*nsidesq
          n(7)=n(8)+2
       case(8)              !East corner
          nneigh=7
          other_face=8+ibp
          n(7)=other_face*nsidesq+local_magic2
          n(1)=n(7)-2
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=4+ibp
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
       end select ! south
    endif

    return
  end subroutine neighbours_nest
  !====================================================================
  subroutine next_in_line_nest(nside, ipix, inext)
    !====================================================================
    !   given nside and a NESTED pixel number ipix, returns in inext
    !  the pixel that lies on the East side (and the same latitude) as ipix
    !
    !   Hacked by EH from BDW's neighbours_nest, 2001-12-18
    !====================================================================
    use bit_manipulation
    integer(kind=i4b), intent(in)::nside, ipix
    integer(kind=i4b), intent(out):: inext

    integer(kind=i4b) :: npix,ipf,ipo,ix,ixp,iy,iym,ixo,iyo
    integer(kind=i4b) :: face_num,other_face
    integer(kind=i4b) :: ia,ib,ibp,ibm,ib2,icase,nsidesq
    integer(kind=i4b) :: local_magic1,local_magic2

    !--------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    nsidesq=nside*nside
    npix = 12*nsidesq       ! total number of points
    if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"

    !     initiates array for (x,y)-> pixel number -> (x,y) mapping
    if (x2pix(128) <= 0) call mk_xy2pix()

    local_magic1=(nsidesq-1)/3
    local_magic2=2*local_magic1
    face_num=ipix/nsidesq

    ipf=modulo(ipix,nsidesq)   !Pixel number in face

    call pix2xy_nest(nside,ipf,ix,iy)
    ixp=ix+1
    iym=iy-1

    !     Exclude corners
    if(ipf==local_magic2)     then !WestCorner
       inext = ipix - 1         
       return
    endif
    if(ipf==(nsidesq-1)) then !NorthCorner
       icase=6             
       goto 100
    endif
    if(ipf==0)           then !SouthCorner
       icase=7                
       goto 100
    endif
    if(ipf==local_magic1)     then !EastCorner
       icase=8      
       goto 100
    endif

    !     Detect edges
    if(IAND(ipf,local_magic1)==local_magic1) then !NorthEast
       icase=1 
       goto 100
    endif
    if(IAND(ipf,local_magic2)==0)      then !SouthEast
       icase=4                   
       goto 100
    endif

    !     Inside a face
    call xy2pix_nest(nside, ixp, iym, face_num, inext)
    return

100 continue

    ia= face_num/4            !in {0,2}
    ib= modulo(face_num,4)       !in {0,3}
    ibp=modulo(ib+1,4)
    ibm=modulo(ib+4-1,4)
    ib2=modulo(ib+2,4)

    if(ia==0) then          !North Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ibp
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          inext = other_face*nsidesq+ipo         ! (6)
       case(4)              !SouthEast edge
          other_face=4+ibp
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, inext)
       case(6)              !North corner
          other_face=0+ibp
          inext=other_face*nsidesq+nsidesq-1
       case(7)              !South corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2+1
       case(8)              !East corner
          other_face=0+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! north

    elseif(ia==1) then      !Equatorial region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ib
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, inext)
       case(4)              !SouthEast edge
          other_face=8+ib
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, inext)
       case(6)              !North corner
          other_face=0+ib
          inext=other_face*nsidesq+local_magic2-2
       case(7)              !South corner
          other_face=8+ib
          inext=other_face*nsidesq+local_magic2+1
       case(8)              !East corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! equator
    else                    !South Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=4+ibp
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, inext)
       case(4)              !SouthEast edge
          other_face=8+ibp
          ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
          inext = other_face*nsidesq+ipo   ! (8)
       case(6)              !North corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2 -2
       case(7)              !South corner
          other_face=8+ibp
          inext=other_face*nsidesq
       case(8)              !East corner
          other_face=8+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! south
    endif

    return
  end subroutine next_in_line_nest
  !=======================================================================
  subroutine ang2vec(theta, phi, vector)
    !=======================================================================
    !     renders the vector (x,y,z) corresponding to angles 
    !     theta (co-latitude measured from North pole, in [0,Pi] radians) 
    !     and phi (longitude measured eastward, in radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    REAL(KIND=DP), INTENT(IN) :: theta, phi
    REAL(KIND=DP), INTENT(OUT), dimension(1:) :: vector

    REAL(KIND=DP) :: sintheta
    !=======================================================================
    
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2VEC: theta : ",theta," is out of range [0, Pi]"
       stop
    endif
    sintheta = SIN(theta)

    vector(1) = sintheta * COS(phi)
    vector(2) = sintheta * SIN(phi)
    vector(3) = COS(theta)

    return
  end subroutine ang2vec
  !=======================================================================
  subroutine vec2ang(vector, theta, phi)
    !=======================================================================
    !     renders the angles theta, phi corresponding to vector (x,y,z)
    !     theta (co-latitude measured from North pole, in [0,Pi] radians) 
    !     and phi (longitude measured eastward, in [0,2Pi[ radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    REAL(KIND=DP), INTENT(IN), dimension(1:) :: vector
    REAL(KIND=DP), INTENT(OUT) :: theta, phi

    REAL(KIND=DP) :: dnorm, z
    !=======================================================================
    
    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)

    z = vector(3) / dnorm
    theta = ACOS(z)

    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

    return
  end subroutine vec2ang
  !=======================================================================
  function npix2nside(npix) result(nside_result)
    !=======================================================================
    ! given npix, returns nside such that npix = 12*nside^2
    !  nside should be a power of 2 smaller than 8192
    !  if not, -1 is returned
    ! EH, Feb-2000
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: npix
!     INTEGER(KIND=I4B) :: npix2nside
    INTEGER(KIND=I4B) :: nside_result
    INTEGER(KIND=I4B) :: nside, ilog
    REAL(KIND=DP) :: fnside, fnpix, flog
    CHARACTER(LEN=*), PARAMETER :: code = "npix2nside"
    INTEGER(KIND=I4B), PARAMETER :: npix_max = 12*ns_max*ns_max
    !=======================================================================

    fnside = sqrt(npix/12.0_dp)
    nside = NINT(fnside)

    if (npix < 12) then
       print*, code//": Npix is too small :",npix
       print*, "                       < 12 "
       nside = -1
       goto 1       
    endif

    if (npix > npix_max) then
       print*, code//": Npix is too large :",npix
       print*, "                       > ",npix_max
       nside = -1
       goto 1
    endif

    fnpix = 12.0_dp*nside*nside
    if (abs(fnpix-npix) > 1.0e-2) then
       print*, code//": wrong Npix ",npix
       print*, "    not 12*nside*nside "
       nside = -1
       goto 1
    endif

    flog = log(real(nside,kind=dp))/log(2.0_dp)
    ilog = NINT(flog)
    if (abs(ilog-flog) > 1.0e-6_dp) then
       print*, code//": wrong Nside ",nside
       print*, "    not a power of 2 "
       nside = -1
       goto 1
    endif

1   continue
!     npix2nside = nside
    nside_result = nside
    return

  end function npix2nside
  !=======================================================================
  function nside2npix(nside) result(npix_result)
    !=======================================================================
    ! given nside, returns npix such that npix = 12*nside^2
    !  nside should be a power of 2 smaller than 8192
    !  if not, -1 is returned
    ! EH, Feb-2000
    !=======================================================================
!     INTEGER(KIND=I4B) :: nside2npix
    INTEGER(KIND=I4B) :: npix_result
    INTEGER(KIND=I4B), INTENT(IN) :: nside

    INTEGER(KIND=I4B), dimension(1:14) :: listnside
    INTEGER(KIND=I4B) :: npix
    CHARACTER(LEN=*), PARAMETER :: code = "nside2npix"
    !=======================================================================

    listnside = (/1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192/)
    npix = 12*nside*nside

    if (MINVAL(abs(listnside-nside)) > 0) then
       print*,code//": unvalid nside ",nside
!        write(unit=*,fmt="(a,4(i2),3(i3),3(i4),4(i5))") " Nside should be among ",listnside
       print "(a,4(i2),3(i3),3(i4),4(i5))", " Nside should be among ",listnside
       npix = -1
    endif

!     nside2npix = npix
    npix_result = npix

    return
  end function nside2npix
  !=======================================================================
  subroutine surface_triangle(vec1, vec2, vec3, surface)
    !=======================================================================
    ! returns the surface in steradians 
    !  of the spherical triangle with vertices vec1, vec2, vec3
    !
    ! algorithm : finds triangle sides and uses l'Huilier formula to compute
    ! "spherical excess" = surface area of triangle on a sphere of radius one
    ! see, eg Bronshtein, Semendyayev Eq 2.86
    !=======================================================================
    real(kind=dp), dimension(1:), intent(in) :: vec1, vec2, vec3
    real(kind=dp), intent(out) :: surface

    !   real(kind=dp), dimension(1:3) :: v1, v2, v3
    real(kind=dp), dimension(1:3) :: side
    !  real(kind=dp) :: hp
    real(kind=dp) :: x0, x1, x2, x3
    !=======================================================================

    ! half perimeter  
!   hp = 0.5_dp * (side1 + side2 + side3)

!   ! l'Huilier formula
!   x0 = tan( hp          * 0.5_dp)
!   x1 = tan((hp - side1) * 0.5_dp)
!   x2 = tan((hp - side2) * 0.5_dp)
!   x3 = tan((hp - side3) * 0.5_dp)

    ! find triangle sides
    call angdist(vec2, vec3, side(1))
    call angdist(vec3, vec1, side(2))
    call angdist(vec1, vec2, side(3))
    ! divide by 4
    side(1:3) = side(1:3) * 0.25_dp

    ! l'Huilier formula
    x0 = tan( side(1) + side(2) + side(3) )
    x1 = tan(-side(1) + side(2) + side(3) )
    x2 = tan( side(1) - side(2) + side(3) )
    x3 = tan( side(1) + side(2) - side(3) )
    surface = 4.0_dp * atan( sqrt(x0 * x1 * x2 * x3) )

    return
  end subroutine surface_triangle
  !=======================================================================
  subroutine angdist(v1, v2, dist)
    !=======================================================================
    ! call angdist(v1, v2, dist)
    ! computes the angular distance dist (in rad) between 2 vectors v1 and v2
    ! in general dist = acos ( v1 . v2 )
    ! except if the 2 vectors are almost aligned.
    !=======================================================================
    real(kind=DP), intent(IN), dimension(1:) :: v1, v2
    real(kind=DP), intent(OUT) :: dist

    real(kind=DP), dimension(1:3) :: r1, r2, vdiff
    real(kind=DP) :: diff, sprod
    !=======================================================================

    ! normalize both vectors
    r1(1:3) = v1(1:3) / sqrt(dot_product(v1,v1))
    r2(1:3) = v2(1:3) / sqrt(dot_product(v2,v2))

    sprod = DOT_PRODUCT(r1, r2)

    if (sprod > 0.999_dp) then
       ! almost colinear vectors
       vdiff(1:3) = r1(1:3) - r2(1:3)
       diff = sqrt(dot_product(vdiff,vdiff)) ! norm of difference
       dist = 2.0_dp * asin(diff * 0.5_dp)

    else if (sprod < -0.999_dp) then
       ! almost anti-colinear vectors
       vdiff(1:3) = r1(1:3) + r2(1:3)
       diff = sqrt(dot_product(vdiff,vdiff)) ! norm of sum
       dist = PI - 2.0_dp * asin(diff * 0.5_dp)

    else
       ! other cases
       dist = acos( sprod )
    endif


    return
  end subroutine angdist
!=======================================================================
  subroutine vect_prod(v1, v2, v3) !
    !=======================================================================
    !     returns in v3 the vectorial product of the 2 3-dimensional 
    !     vectors v1 and v2
    !=======================================================================
    real(kind=DP), dimension(1:), INTENT(IN)  :: v1, v2
    real(kind=DP), dimension(1:), INTENT(OUT) :: v3
    !=======================================================================
    
    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
    
    return
  end subroutine vect_prod

!=======================================================================
!=======================================================================
!=======================================================================

  subroutine chpix_nest (nsides_old,nsides_new,ipix_old,ipix_new)
  
  !---------------------------------------------------
  ! Added routine by Marcelo Lares.
  ! Given a pixel number in a pixelization squeme,
  ! returns the corresponding pixel number in another
  ! squeme (different nsides) both in the NESTED squeme.
  ! This code is a copy of routines Int2Binary and Binary2Int
  ! from odule GRID_tools.
  ! required: nsides_new < nsides_old
  !----------------------------------------------------

  implicit none
  integer,intent(in)   :: nsides_old,nsides_new,ipix_old
  integer,intent(out)  :: ipix_new

  integer              :: N, shift, i, k, j, s, ind
  integer,parameter    :: Nb=30
  integer,dimension(Nb):: B
  integer,dimension(1) :: L

  ! MAIN

  if(nsides_old<nsides_new) write(*,*) 'nsides_old must be > nsides_new'

  ! (1) Calculate jump in bits between the two squemes
 
  shift = nsides_old / nsides_new

  s = shift
  
  do j=1,Nb
    k   = mod(s,2)
    ind = Nb-j+1
    B(ind) = k
    s = (s-k)/2
  enddo
  L = maxloc(B)
  shift = (Nb - L(1))*2

  ! (2) Binary representation of the pixel (old)

  call Int2Binary(ipix_old,B)
 
  ! (3) Binary representation of the new pixel
  
  N = Nb-shift

  i = 0
  do ind=1,N
    k = N - ind
    i = i + B(ind) * 2**k
  enddo
  
  ipix_new = i

  end subroutine chpix_nest

!=======================================================================
  subroutine introspixtion (nsides_old,nsides_new,npixs)
  
     !---------------------------------------------------
     ! Added routine by Marcelo Lares.
     ! Given two pixel squemes, calculates how many pixels
     ! do a macro-pixel contain.
     ! required: nsides_new > nsides_old
     !----------------------------------------------------

     implicit none
     
     integer,intent(in)   :: nsides_old,nsides_new
     integer,intent(out)  :: npixs
     
     integer,parameter    :: N=30
     integer,dimension(N) :: B
     integer,dimension(1) :: L
     integer              :: shift, i, j, k, s, ind

     ! MAIN
     
     shift = nsides_new / nsides_old

     s = shift
     B = 0 
    
     do j=1,N
       k   = mod(s,2)
       ind = N-j+1
       B(ind) = k
       s = (s-k)/2
     enddo
     L = maxloc(B)
     npixs = N - L(1)
     npixs = 2**(2*npixs)

  end subroutine introspixtion
  !=======================================================================
  subroutine spreadpix_nest (pix_old,newpixs,ipixs)
  
  !---------------------------------------------------
  ! Added routine by Marcelo Lares.
  ! Given a pixel number in a pixelization squeme,
  ! returns the corresponding set of pixel numbers
  ! in another more refinated squeme, both in the NESTED squeme.
  ! This code is a copy of routines Int2Binary and Binary2Int
  ! from module GRID_tools.
  !----------------------------------------------------

    implicit none
    integer,intent(in)                      :: pix_old
    integer,intent(in)                      :: newpixs
    integer,dimension(newpixs),intent(out)  :: ipixs

    integer              :: i, k, j, s
    integer              :: shift, ind
    integer,parameter    :: N=30
    integer,dimension(N) :: B_new, B_old, B
    integer,dimension(1) :: L
    integer,allocatable  :: tail(:)

    ! MAIN

    ! (1) Binary representation of the integer number pix_old

    call Int2Binary(pix_old,B_old)

    ! (2) Shift the numbers

    call Int2Binary(newpixs,B)

    if(count(B==1)>1) write(*,*) &
    'Error@spreadpix_nest: newpixs is not a power of two'

    L = maxloc(B)
    shift = N - L(1)

    B_new = 0
    do i=1,N
       j = max(i-shift,1)
       B_new(j) = B_old(i)
    enddo

    ! (3) Binary representation of new pixels

    allocate(tail(shift))

    do i=1,newpixs
       call Int2Binary (i-1,tail)
       B_new(N-shift+1:) = tail
       call Binary2Int (B_new,ipixs(i))
    enddo

    deallocate(tail)

  end subroutine spreadpix_nest
  !======================================================================

  !-----------------------------------------------------------
   !this routine writes the binary representation
   !of an integer number <i> in a vector <B>.
   !Warning:the size of the vector must be big enough!
   ! NOTE:  this is a copy of the original routine from >>>>>> GRID_tools <<<<<<
   !-----------------------------------------------------------

  subroutine Int2Binary (i,B)

    !DIM&DEC
    implicit none
    integer,intent(in)               :: i
    integer,dimension(:),intent(out) :: B
 
    integer :: s, k, ind, j
    integer :: N(1)

    !MAIN

    N = size(B)

    if(i>=2**N(1)) &
    write(*,*) 'invalid number in Int2Binary',i

    s = i

    do j=1,N(1)

      k = mod(s,2)
      ind = N(1)-j+1
      B(ind) = k
      s = (s - k)/2

    enddo

  end subroutine Int2Binary

 !======================================================================

!this routine writes the binary representation
!of an integer number <i> in a vector <B>.
!Warning:the size of the vector must be big enough!

subroutine Binary2Int (B,i)

  !DIM&DEC
  implicit none
  integer,dimension(:),intent(in)  :: B
  integer,intent(out)              :: i

  integer :: ind, k
  integer :: N(1)

  !MAIN

  N = size(B)
  i = 0

  do ind=1,N(1)

    k = N(1) - ind
    i = i + B(ind) * 2**k

  enddo

end subroutine Binary2Int

!======================================================================
  
end module pix_tools
