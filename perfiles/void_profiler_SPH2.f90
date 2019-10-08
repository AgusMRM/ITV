module datas
  use healpix_types, only : DP, SP
  implicit none
  real(kind=DP),dimension(:),allocatable :: x,y,z
  integer,parameter::pr=selected_real_kind(10,20)
  real(pr),dimension(:),allocatable :: b2!,by,bz
  real(pr),dimension(:),allocatable :: rho!,by,bz
  real(pr),dimension(:),allocatable :: Hsml!,by,bz
  real(SP),dimension(:),allocatable :: xe,ye,ze,rv
  real       :: dbox
  integer(4) :: parbox
  integer :: ncen
end module datas

program cicjack_omp
  use datas
  use OMP_lib
  use healpix_types, only : DP, SP
  use pix_tools, only : pix2vec_nest,nside2npix
  implicit none
  type trivect
          real(DP),dimension(3) :: v
  end type trivect
  !!type pixshell
  !!        integer                               :: npix
  !!        type(trivect),dimension(:)allocatable :: pix
  !!end type pixshell
  !!type(pixshell),dimension(:),allocatable :: grid_pixels
  type(trivect),dimension(:),allocatable :: grid_pixels
  real(DP),dimension(3) :: v
  real(DP),parameter              :: rmin=1000.
  integer,parameter           :: nthread   = 32 
  integer,parameter           :: nside     = 60
  integer,parameter           :: nbindist  = 100
  integer,parameter           :: nsideheal = 2
  real,parameter              :: rmax_vu=10. , rmin_vu=0.05

  integer                     :: kk1,kk2,ipx,nwrites
  integer ,dimension(nthread) :: seed
  integer ,dimension(:),allocatable   :: nes
  real(pr),dimension(:,:),allocatable :: dens
  real(pr),dimension(:,:),allocatable :: dens_h
  real(pr),dimension(:,:),allocatable :: bfld
  real(pr),dimension(:),allocatable :: dens_pix
  real(pr),dimension(:),allocatable :: dens_h_pix
  real(pr),dimension(:),allocatable :: bfld_pix
  real(DP)       :: rdist,vol,rdist_
  integer    :: ibin,ibindist,isum
  integer    :: i,j,it,l,iix,iiy,iiz,ix,iy,iz,indx,indy,indz,ld
  real(DP)       :: remax, xc,yc,zc,uc,vc,wc,facbin       
  real(DP)       :: lrmax_vu,lrmin_vu,dbindist_vu       
  real(DP)       :: dbindist,distlog,dist,dc,rad
  real(pr)   :: distmedbin

  character(len=100) :: fort16, boxfilep 
  character(len=2) :: ixs,iys,izs,ibs,jbs,kbs 
  integer,dimension(nside,nside,nside):: lirst
  integer,dimension(:),allocatable :: ll
  real(DP) :: volbox,pi
  real(DP) :: xp,yp,zp
  real(pr):: bb
  real(pr):: rr
  logical :: ilogic
  real(DP) :: iHsml_min,iHsml_max
  real(pr) :: fac
  integer :: iHsml
  integer :: npix


  !call getreal(1, remax)
  !call getreal(2, dbindist)  ! ver q remax y dbindist sean multiplos 
  call read_voids() !xe,ye,ze,rv,ncen

  print*,'leyendo sim....'

  call read_sim() !x,y,z,bx,by,bz
  print*,'parbox',parbox
  print*,'dbox  ',dbox

  allocate(ll(size(x)))
  print*,'grideando....'
  facbin   = real(nside)/dbox
  remax    = dbox/real(nside)

  call grid(x,y,z,lirst,ll)

  allocate(dens  (ncen,nbindist))
  allocate(dens_h(ncen,nbindist))
  allocate(bfld  (ncen,nbindist))


  !! lo usare a esto si quiero diferentes pixeles por shell
  !!npix = 12*( 2 )**2  !12*nsideheal**2
  !!do i=1,nbindist
  !!    !quizas se puede aumentar el npix con el radio
  !!    !tengo mis dudas sobre la performance
  !!    grid_pixels%npix = npix
  !!    allocate(grid_pixels(i)%pix(npix))

  !!enddo

  npix = 12*( nsideheal )**2 
  allocate(grid_pixels(npix))

  do i=0,npix-1
      call pix2vec_nest(nsideheal,i,v)
      grid_pixels(i+1)%v=v
  enddo
  
  lrmax_vu=alog10(rmax_vu)
  lrmin_vu =alog10(rmin_vu)
  dbindist_vu = (lrmax_vu - lrmin_vu) / nbindist

  volbox   = dbox*dbox*dbox
  pi       = 4.*atan(1.)
  write(*,*)'nbinsit', nbindist
  !call OMP_set_dynamic(.true.)
  call OMP_set_num_threads(nthread)

  dens   = -999.999
  dens_h = -999.999
  bfld   = -999.999


!!!!!!mmdens=sum(rho)/size(rho)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(l,ix,iy,iz,iix,iiy,iiz, &
!$OMP indx,indy,indz,it,dc,distlog,dist,rad, &
!$OMP isum,iHsml,iHsml_min,iHsml_max,fac, &
!$OMP ibin,xp,yp,zp,xc,yc,zc,ilogic,bb,rr,distmedbin,ipx, &
!$OMP dens_pix,bfld_pix,dens_h_pix) &
!$OMP SHARED(xe,ye,ze,rv,x,y,z,nes,facbin,b2,hsml, &
!$OMP rho,bfld,dens_h,dens,grid_pixels,pi, &
!$OMP ll,lirst,remax,dbindist,ncen,lrmin_vu,lrmax_vu,dbindist_vu,npix)

  allocate(   dens_pix(npix))
  allocate( dens_h_pix(npix))
  allocate(   bfld_pix(npix))

!$OMP BARRIER

!--------------------------calculo part x esf------------------------
!$OMP DO SCHEDULE(DYNAMIC)
  do l=1, ncen  ! tomo centro random
            nes(l) = 0
            call Show_progress(l,ncen,60,':','ran sphere   ')
            radio: do ibin=1,nbindist
               !ibindist = floor((distlog-lrmin_vu)/dbindist_vu+ 0.5)+1
               !distlog = alog10(dist/rv(l))
  
               distlog = (ibin - 0.5)*dbindist_vu + lrmin_vu
               rad    = (10.0**distlog)*rv(l)
  
               dens_pix   = 0.0
               dens_h_pix = 0.0
               bfld_pix   = 0.0
  
               pixel: do ipx=1,npix
                  xc =xe(l) + grid_pixels(ipx)%v(1)*rad 
                  yc =ye(l) + grid_pixels(ipx)%v(2)*rad 
                  zc =ze(l) + grid_pixels(ipx)%v(3)*rad 
                  call periodoBX(xc)
                  call periodoBX(yc)
                  call periodoBX(zc)
  
                  ix = int(xc*facbin) + 1 
                  iy = int(yc*facbin) + 1 
                  iz = int(zc*facbin) + 1
  
                  pri: do iix = ix - 1,ix + 1
                    call periodoIND(iix,indx)
                    seg: do iiy = iy - 1,iy + 1
                      call periodoIND(iiy,indy)
                      ter: do iiz = iz - 1,iz + 1
                        call periodoIND(iiz,indz)
                        it = lirst(indx,indy,indz)
                        if(it /= 0) then
                          vecino: do
                            it = ll(it)
                            xp =x(it) ; yp=y(it) ; zp=z(it)
                            bb = b2(it)
                            rr = rho(it)
                            call periodoCM(xc,xp)
                            call periodoCM(yc,yp)
                            call periodoCM(zc,zp)
                            dc = (xc-xp)**2+(yc-yp)**2+(zc-zp)**2
                            dist = sqrt(dc)
  
                            if(dist < Hsml(it))then
                               call KernelGadget(dist, hsml(it),fac) 
                               bfld_pix  (ipx) = bfld_pix  (ipx)+       bb*fac
                               dens_pix  (ipx) = dens_pix  (ipx)+       rr*fac
                               dens_h_pix(ipx) = dens_h_pix(ipx)+ hsml(it)*fac
                            endif
                            if(it == lirst(indx,indy,indz))exit vecino
                          end do vecino
                        endif
  
                      end do ter
                    end do seg
                  end do pri 
               enddo pixel
  
               bfld  (l,ibin) = sum(  bfld_pix)/npix*4*pi
               dens  (l,ibin) = sum(  dens_pix)/npix*4*pi
               dens_h(l,ibin) = sum(dens_h_pix)/npix*4*pi
            enddo radio
  
  end do  ! termina centro voids
  !$OMP END DO
  !$OMP BARRIER
  deallocate(   dens_pix)
  deallocate( dens_h_pix)
  deallocate(   bfld_pix)

  !$OMP END PARALLEL 
  
  print*,'Escribiendo...'
  open(44,file='./profile_GAS.dat',status='unknown',access='stream')

  nwrites=0
  write(44)int(ncen,4)
  nwrites=nwrites+1
  write(44)int(nbindist,4)
  nwrites=nwrites+1
  write(44)real(10.**lrmin_vu,4)
  nwrites=nwrites+1
  write(44)real(10.**lrmax_vu,4)
  nwrites=nwrites+1
  do i=1,ncen
    write(44)real(rv(i),4)
      do j=1,nbindist
        write(44)real(dens(i,j),4)
  nwrites=nwrites+1
      enddo
      do j=1,nbindist
        write(44)real(dens_h(i,j),4)
  nwrites=nwrites+1
      enddo
      do j=1,nbindist
        write(44)real(bfld(i,j),4)
  nwrites=nwrites+1
      enddo
  enddo
  print*,nwrites,' <----- nwrites'
  deallocate(xe,ye,ze,rv)
  deallocate(nes)

contains
subroutine Show_progress(idx,nbt,howmany,deco,mess)
           ! uso: call show_progress(loop index, Nloops, rayitas)
           ! optional: caracter for decoration
           ! optional: message (max 10 characters)

           implicit none
           integer,intent(in) :: idx
           integer,intent(in) :: nbt

           integer,intent(in),optional          :: howmany
           character(len=1),intent(in),optional :: deco
           character(len=*),intent(in),optional :: mess
           
           integer            :: nprints, ntot, mark, i
           character(len=1)   :: A
           character(len=10)  :: screen
           character(len=50)  :: frmt 

           ! MAIN

           nprints = 10
           if(present(howmany)) nprints=max(howmany,10)
           ntot = nbt - mod(nbt,nprints)
           mark = max(ntot/nprints,1)

           if(nbt<nprints) then
              write(*,*) idx, nbt
              return
           endif

           if(idx==nbt) write(*,'(a/)',advance='yes') ' '
           if(idx>ntot) return

           A = ':'
           if(present(deco)) A = deco

           screen = 'progress '
           if(present(mess)) screen=mess(:10)

           if(idx==1) then
              if(nbt<99999) then
                write(frmt,'(''(t'',i2,'',f6.0,x,a1)'')') nprints+10
              else
                write(frmt,'(''(t'',i2,'',e8.2,x,a1)'')') nprints+10
              endif
              write(*,fmt=frmt,advance='yes') real(nbt),'|'
              write(*,'(a10,$)') screen
           endif
           if(mod(idx,mark)==0) write(*,'(a1,$)') A
           if(idx==ntot) write(*,'(a1,$)') A

end subroutine Show_progress  !}}}*/
    subroutine read_sim()
            implicit none
            real(4):: xtmp,ytmp,ztmp
            integer :: jj,i
            integer :: ppp
            character(len=100) :: path='./'

! Dark Gas 
              open(11,file=trim(path)//'magneticum_gas.dat',access='stream',status='old')
              parbox=204101470
       	      print*,'Gas mode'
!
            dbox=352000
!!            jj=0
!!            ppp=0
!!            do 
!!                ppp=ppp+1
!!                read(11,iostat=jj)xtmp 
!!                if(jj /=0)exit
!!            enddo
!!            print*,ppp
!!            ppp=ppp/3
!!            rewind(11)
!!            print*,ppp,parbox
   
            allocate(x(parbox),y(parbox),z(parbox))
            do i=1,parbox
                call Show_progress(i,parbox,60,':','lectura   ')
                read(11)xtmp 
                read(11)ytmp 
                read(11)ztmp 
                x(i)=xtmp 
                y(i)=ytmp 
                z(i)=ztmp 
            enddo
            close(11)

            allocate(b2(parbox))
            allocate(rho(parbox))
            allocate(hsml(parbox))
	    b2=0.0_pr
! Dark Gas 
            open(11,file=trim(path)//'magneticum_B.dat',access='stream',status='old')
            do i=1,parbox
                call Show_progress(i,parbox,60,':','lectura   ')
                read(11)xtmp 
!                read(11)ytmp 
!                read(11)ztmp 
                !b2(i)=(real(xtmp,pr)*real(xtmp,pr))!+real(ytmp,pr)*real(ytmp,pr)+real(ztmp,pr)*real(ztmp,pr))
                b2(i)=(real(xtmp,pr))
            enddo
            close(11)
            open(11,file=trim(path)//'magneticum_RHO.dat',access='stream',status='old')
            do i=1,parbox
                call Show_progress(i,parbox,60,':','lectura   ')
                read(11)xtmp 
                rho(i)=(real(xtmp,pr))
            enddo
            close(11)
            open(11,file=trim(path)//'magneticum_HSML.dat',access='stream',status='old')
            do i=1,parbox
                call Show_progress(i,parbox,60,':','lectura   ')
                read(11)xtmp 
                hsml(i)=(real(xtmp,pr))
            enddo
            close(11)
	    print*,"Valores B Min/Max",minval(b2),maxval(b2)
	    print*,"Valores Rho Min/Max",minval(rho),maxval(rho)
	    print*,"Valores hsml Min/Max",minval(hsml),maxval(hsml)
    end subroutine read_sim
    subroutine read_voids()
            implicit none
            real(4):: xtmp,ytmp,ztmp
            integer :: i
            character(len=100)::filename

            !call getstr (1, filename)
            filename='./CatalogoMag_144_GAS.txt'
            open(11,file=filename,form='formatted',status='old')
            call count_lines(11,ncen)
            allocate(xe(ncen),ye(ncen),ze(ncen),rv(ncen))
            allocate(nes(ncen))

            write(*,*)'leyendo voids...'
            read(11,*)(rv(i),xe(i),ye(i),ze(i),i=1,ncen)
            close(11)
    end subroutine read_voids
    subroutine count_lines(nunit,nlines)
      implicit none
      integer, intent(in)   :: nunit
      integer, intent(out)  :: nlines
      integer               :: check
      nlines = 0
      do
          read(nunit,*,iostat=check)
          if(check /= 0)exit
          nlines = nlines + 1
      end do
      rewind(nunit)
    end subroutine count_lines
    subroutine grid(x,y,z,lirst,ll)
     implicit none 
     integer                                            :: i,ind,ix,iy,iz
     integer,dimension(nside,nside,nside),intent(inout) :: lirst
     integer,dimension(:),intent(out)                   :: ll
     real(kind=DP),dimension(:),intent(in)              :: x,y,z
     
     lirst = 0
     ll    = 0
     print*,facbin,'--------------------------'
     if(facbin < 0.0) stop 'definir facbin antes de llamar grid'

     write(*,*) 'grid: start'
     do i = 1,size(ll)
       ix = int(x(i)*facbin)+1
       iy = int(y(i)*facbin)+1
       iz = int(z(i)*facbin)+1
        
       lirst(ix,iy,iz) = i
     end do

     do i = 1,size(ll)
       ix = int(x(i)*facbin)+1
       iy = int(y(i)*facbin)+1
       iz = int(z(i)*facbin)+1     
       ll(lirst(ix,iy,iz)) = i
       lirst(ix,iy,iz)   = i
     end do
     write(*,*) 'grid: end'

    end subroutine grid
    subroutine periodoCM(uc,u)
       implicit none
       real(kind=DP),intent(in)    :: uc
       real(kind=DP),intent(inout) :: u
       if((uc-u)> dbox/2.)u=u+dbox
       if((uc-u)<-dbox/2.)u=u-dbox
    end subroutine periodoCM
    subroutine periodoBX(u)
       implicit none
       real(kind=DP),intent(inout) :: u
       if(u>=dbox)u=u-real(dbox,4)
       if(u<0._4)u=u+real(dbox,4)
    end subroutine periodoBX
    subroutine periodoIND(inte,into)
       implicit none
       integer, intent(in)         :: inte
       integer, intent(inout)      :: into
       into = inte
       if(inte > nside)into = 1
       if(inte < 1)into = nside
    end subroutine periodoIND
    subroutine getstr(i,a)
      implicit none
      integer                        :: i
      character(len=*),intent(inout) :: a
      call getarg(i,a)
    end subroutine getstr
    subroutine getint(i,ii)
      implicit none
      integer           ::  i,ii
      character(len=20) ::  a
      call getarg(i,a)
      read(a,*) ii
    end subroutine getint
    subroutine getreal(i,r)
      implicit none
      integer           ::  i
      real              ::  r
      character(len=20) ::  a
      call getarg(i,a)
      read(a,*) r
    end subroutine getreal
    subroutine getreal8(i,r)
      implicit none
      integer           ::  i
      real(8)           ::  r
      character(len=20) ::  a
      call getarg(i,a)
      read(a,*) r
    end subroutine getreal8
    subroutine KernelGadget(dist,hsml,fac)
       implicit none
       real(pr),intent(in)    :: dist,hsml
       real(pr),intent(inout) :: fac
       real(pr)           ::  u,hinv,hinv3
       real(pr)           ::  COEF_1,COEF_2,COEF_5,NORM
       COEF_1 =  2.54647908947
       COEF_2 = 15.278874536822
       COEF_5 =  5.092958178941
!!!!!       NORM   = 1.0/3.1415 !!! si usas hinv3=1, recordar que q=2u (paper de price)
       NORM   = 8.0/3.1415

       if(dist>hsml) then 
          fac=0.0
       else
          hinv=1./hsml
          hinv3=(1./hsml)**3.
          u=dist*hinv
!!!!	 hinv3=1.0 ! NO damos en Kernel como Volumen
          if(u<0.5) then 
             fac=hinv3*(1.0 + 6.0 * ( u - 1.0)*u*u)*NORM 
!             fac=hinv3*(COEF_1 + COEF_2*(u-1.)*u*u) 
          else
             fac=hinv3 * 2.0* (1.-u) * (1.-u) * (1.-u) *NORM
          endif
       endif
     end subroutine KernelGadget
end program cicjack_omp

