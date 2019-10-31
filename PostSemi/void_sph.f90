
module datas
  use healpix_types, only : DP, SP
  implicit none
  integer,parameter::pr=selected_real_kind(10,20)
  real(pr),dimension(:),allocatable :: b2!,by,bz
!  real(pr),dimension(:),allocatable :: rho!,by,bz
!  real(pr),dimension(:),allocatable :: Hsml!,by,bz
  real       :: dbox
  integer(4) :: parbox
  integer :: ncen
  real(SP) :: xbox,ybox,zbox,xe,ye,ze,rv
end module datas
program cicjack_omp
  use modulos
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
 ! real(DP),parameter              :: rmin=1000.
  integer,parameter           :: nthread   =48, VOID=1373 
  integer,parameter           :: nside     = 100
  integer,parameter           :: nbindist  = 500
  integer,parameter           :: nsideheal = 2
  real,parameter              :: rmax_vu=25. , rmin_vu=0.2
  real,parameter              :: tlim=10**5  , dlim=100
  integer                     :: kk1,kk2,ipx,nwrites
  integer ,dimension(nthread) :: seed
  integer ,dimension(:),allocatable   :: nes
  real(pr),dimension(:,:),allocatable :: dens_h
  real(pr),dimension(:,:),allocatable :: bfld
  real(pr),dimension(:),allocatable :: dens_pix,dens1_pix,dens2_pix
  real(pr),dimension(:),allocatable :: dens3_pix,dens4_pix
  real(pr),dimension(:,:),allocatable ::dens2 ,dens_jack
  real(pr),dimension(:,:),allocatable ::dens2_1, dens2_2, dens2_3, dens2_4 
  real(pr),dimension(:,:),allocatable ::dens1_jack, dens2_jack, dens3_jack, dens4_jack 
  real(pr),dimension(:),allocatable :: dens_h_pix
  real(pr),dimension(:),allocatable :: bfld_pix
  real(pr),dimension(:),allocatable :: x,y,z
  real(DP)       :: rdist,vol,rdist_
  integer    :: ibin,ibindist,isum,o,ijack
  integer    :: it,l,iix,iiy,iiz,ix,iy,iz,indx,indy,indz,ld
  real(DP)       :: remax, xc,yc,zc,uc,vc,wc,facbin       
  real(DP)       :: lrmax_vu,lrmin_vu,dbindist_vu       
  real(DP)       :: dbindist,distlog,dist,dc,rad
  real(pr)   :: distmedbin

  character(len=100) :: fort16, boxfilep 
  character(len=2) :: ixs,iys,izs,ibs,jbs,kbs 
  integer,dimension(nside,nside,nside):: lirst
  integer,dimension(:),allocatable :: ll
  real(DP) :: volbox
  real(DP) :: xp,yp,zp
  real(pr):: bb
  real(pr):: rr
  logical,dimension(:),allocatable :: ilogic
  real(DP) :: iHsml_min,iHsml_max
  real(pr) :: fac
  integer :: iHsml
  integer :: npix
  real ::f1,f2,f3,f4,ftot
  real :: rhomin, rhomax,rad0,promedio,prom,densidad,dens3,error
  real :: pi, temp,densidad1,densidad2,densidad3,densidad4
  real :: promedio1,promedio2,promedio3,promedio4
  real ::dmean,error1,error2,error3,error4
  real,dimension(nbindist) :: jmed,jsig,sig
  real,dimension(nbindist) :: jmed1,jsig1,sig1
  real,dimension(nbindist) :: jmed2,jsig2,sig2
  real,dimension(nbindist) :: jmed3,jsig3,sig3
  real,dimension(nbindist) :: jmed4,jsig4,sig4
  pi       = 4.*atan(1.)
  dmean=(3*(100**2)*(0.045))/(8*pi*(4.3e-9))*1e-10
 ! dmean=1
     xbox=403.896
     ybox=459.8882
     zbox=440.9021 
      !   xbox=411.217
      !   ybox=162.1655
      !   zbox=453.0553 
        open(13,file='/mnt/is2/fstasys/ITV/base09/voids/voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo
        read(13,*) rv,xe,ye,ze
        xe=xe-xbox+250
        ye=ye-ybox+250
        ze=ze-zbox+250
  ncen = 1
close(13)

  print*,'leyendo sim....'

  call reader() !rho
!------------------------------------------------------------  
!------------------------------------------------------------  
!------------------------------------------------------------  
!------------------------------------------------------------  
f1=0;f2=0;f3=0;f4=0
write(*,*) 'CALCULO FRACCION DE CADA REGION'
do i=1,nall(0)
         
         dc=sqrt((pos(1,i)-xe)**2+(pos(2,i)-ye)**2+(pos(3,i)-ze)**2)
         if (dc<rmax_vu) then 
                call temperature(u(i),ne(i),temp)
                if (temp<tlim .and. (rho(i)/dmean)<dlim) then
                f1=f1+1
                elseif (temp<tlim .and. (rho(i)/dmean)>=dlim) then
                f2=f2+1
                elseif (temp>=tlim .and. (rho(i)/dmean)<dlim) then
                f3=f3+1
                elseif (temp>=tlim .and. (rho(i)/dmean)>=dlim) then
                f4=f4+1
                endif
         endif
enddo
print*, f1,f2,f3,f4
!------------------------------------------------------------  
!------------------------------------------------------------  
!------------------------------------------------------------  
!------------------------------------------------------------  
  dbox = 500     !al parecer, el dbox es el tamaño del box en Kpc
  parbox = nall(0) 

  allocate(x(nall(0)),y(nall(0)),z(nall(0)))
  allocate(ll(size(x)))
  print*,'grideando....'
  facbin   = real(nside)/dbox
  remax    = dbox/real(nside)
! grid me va a generar la linked list
 print*, int(dbox/rmax_vu)
! if (nside > int(dbox/rmax_vu)) stop 'nside muy chico' 
do i=1,nall(0)
        x(i)=pos(1,i)
        y(i)=pos(2,i)
        z(i)=pos(3,i)
enddo
  call grid(x,y,z,lirst,ll)
  allocate(dens2  (ncen,nbindist))
  allocate(dens2_1(ncen,nbindist),dens2_2(ncen,nbindist))
  allocate(dens2_3(ncen,nbindist),dens2_4(ncen,nbindist))
  npix = 12*( nsideheal )**2                 ! ESTO NO LO ENTIENDO
  allocate(grid_pixels(npix))
  allocate(ilogic(npix))
  allocate(dens_jack(nbindist,npix))
  allocate(dens1_jack(nbindist,npix))
  allocate(dens2_jack(nbindist,npix))
  allocate(dens3_jack(nbindist,npix))
  allocate(dens4_jack(nbindist,npix))
  
  do i=0,npix-1
      call pix2vec_nest(nsideheal,i,v)
      grid_pixels(i+1)%v=v
  enddo
!------------------------------------------------------------------------------
  
  lrmax_vu=alog10(rmax_vu)
  lrmin_vu =alog10(rmin_vu)
  dbindist_vu = (lrmax_vu - lrmin_vu) / nbindist

  volbox   = dbox*dbox*dbox
  write(*,*)'nbinsit', nbindist
! call OMP_set_num_threads(nthread)   !si no llamo esta subrutina la compu
                                         !elige aumeaticamente el numero de cores
  dens2   = -999.999
  dens2_1   = -999.999
  dens2_2   = -999.999
  dens2_3   = -999.999
  dens2_4   = -999.999
  do ibin=1,nbindist
               distlog = (ibin - 0.5)*dbindist_vu + lrmin_vu
               rad    = (10.0**distlog)!*rv
  enddo
  allocate(dens_pix(npix))
  allocate(dens1_pix(npix),dens2_pix(npix),dens3_pix(npix),dens4_pix(npix))
l=1
              dens_jack  = 0.0
              dens1_jack  = 0.0
              dens2_jack  = 0.0
              dens3_jack  = 0.0
              dens4_jack  = 0.0
 !do l=1, ncen  ! en mi caso, voy de uno a uno 
  !-------------------------------------------------------------------------------------------
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP PRIVATE(distlog,ibin,rad,rv,bfld_pix, &
      !$OMP xc,yc,zc,ix,iy,iz,iix,iiy,iiz,it,bb,rr,dc,dist,fac,xp,yp,zp,v,ipx, &
      !$OMP dens_pix,dens1_pix,dens2_pix,dens3_pix,dens4_pix, & 
      !$OMP ilogic,ijack,temp) & 
      !$OMP SHARED (lrmin_vu,xe,ye,ze,grid_pixels,lirst,ll,x,y,z,b2,rho,hsml,dens2, &
      !$OMP npix,facbin,pos,l,pi,dbindist_vu,mass,dens_jack,ne,u, &
      !$OMP dens2_1,dens2_2,dens2_3,dens2_4,dmean,dens1_jack,dens2_jack,dens3_jack,dens4_jack )
      !$OMP DO SCHEDULE(DYNAMIC)           !      L O   C A M B I E   D E    L U G A R 
            radio: do ibin=1, nbindist    
              dens_pix   = 0.0
              dens1_pix   = 0.0
              dens2_pix   = 0.0
              dens3_pix   = 0.0
              dens4_pix   = 0.0
               distlog = (ibin - 0.5)*dbindist_vu + lrmin_vu
               rad    = (10.0**distlog)         !*rv
               pixel: do ipx=1,npix
                  xc = xe + grid_pixels(ipx)%v(1)*rad 
                  yc = ye + grid_pixels(ipx)%v(2)*rad 
                  zc = ze + grid_pixels(ipx)%v(3)*rad 
  
                  ix = int(xc*facbin) + 1
                  iy = int(yc*facbin) + 1 
                  iz = int(zc*facbin) + 1
               
                  pri: do iix = ix - 1,ix + 1
                    seg: do iiy = iy - 1,iy + 1
                      ter: do iiz = iz - 1,iz + 1
                        it = lirst(iix,iiy,iiz)
               
                        if(it /= 0) then
                          vecino: do
                            it = ll(it)
                            call temperature(u(it),ne(it),temp)
                            xp =pos(1,it) ; yp=pos(2,it) ; zp=pos(3,it)
                        !    bb = b2(it)
                            rr = rho(it)
                            dc = (xc-xp)**2+(yc-yp)**2+(zc-zp)**2
                            dist = sqrt(dc)
                            if(dist < Hsml(it))then
                               call KernelGadget(dist, hsml(it),fac)

                                ! ASIGNO LAS PARTICULAS A CADA REGION Y CALCULO
                                ! LAS DENSIDAD RESPECTIVAS
                                if (temp<tlim .and. (rho(it)/dmean) <dlim) then
                                dens1_pix(ipx) = dens1_pix(ipx) +(mass(it))*fac
                                elseif (temp<tlim .and. (rho(it)/dmean) >=dlim) then
                                dens2_pix(ipx) = dens2_pix(ipx) +(mass(it))*fac
                                elseif (temp>=tlim .and. (rho(it)/dmean) <dlim) then
                                dens3_pix(ipx) = dens3_pix(ipx) +(mass(it))*fac
                                elseif (temp>=tlim .and. (rho(it)/dmean) >=dlim) then
                                dens4_pix(ipx) = dens4_pix(ipx) +(mass(it))*fac
                                endif


                            endif
                            if(it == lirst(iix,iiy,iiz)) exit vecino
                          end do vecino
                        endif
  
                      end do ter
                    end do seg
                  end do pri 
               enddo pixel
               dens2_1(l,ibin) = sum(dens1_pix)/npix   
               dens2_2(l,ibin) = sum(dens2_pix)/npix   
               dens2_3(l,ibin) = sum(dens3_pix)/npix   
               dens2_4(l,ibin) = sum(dens4_pix)/npix   
                do ijack=1,npix
                        ilogic=.TRUE.
                        ilogic(ijack)=.FALSE.
                        dens1_jack(ibin,ijack) = sum(dens1_pix,mask=ilogic)/npix
                        dens2_jack(ibin,ijack) = sum(dens2_pix,mask=ilogic)/npix
                        dens3_jack(ibin,ijack) = sum(dens3_pix,mask=ilogic)/npix
                        dens4_jack(ibin,ijack) = sum(dens4_pix,mask=ilogic)/npix
                enddo
            enddo radio
      !$OMP END DO
      !$OMP BARRIER
      !$OMP END PARALLEL 
                !write(*,*) dens_jack(10,:)
!  end do  ! termina centro voids
write(*,*) 'Escribiendo...'
 open(10,file='profile_gas.dat',status='unknown',action='write')
  deallocate(dens_pix)
  deallocate(dens1_pix,dens2_pix,dens3_pix,dens4_pix)
jmed1=0.0; jsig1=0.0; sig1=0.0 
jmed2=0.0; jsig2=0.0; sig2=0.0 
jmed3=0.0; jsig3=0.0; sig3=0.0 
jmed4=0.0; jsig4=0.0; sig4=0.0 
  
ftot=f1+f2+f3+f4
 do i=1,nbindist
                distlog = (i - 0.5)*dbindist_vu + lrmin_vu
                rad    = (10.0**distlog)         !*rv
               
                do ijack=1,npix
                jmed1(i)=jmed1(i)+dens1_jack(i,ijack)
                jmed2(i)=jmed2(i)+dens2_jack(i,ijack)
                jmed3(i)=jmed3(i)+dens3_jack(i,ijack)
                jmed4(i)=jmed4(i)+dens4_jack(i,ijack)
                enddo
                jmed1(i) = jmed1(i)/npix
                jmed2(i) = jmed2(i)/npix
                jmed3(i) = jmed3(i)/npix
                jmed4(i) = jmed4(i)/npix
                
                do ijack=1,npix
                jsig1(i)=jsig1(i)+(dens1_jack(i,ijack)-jmed1(i))**2
                jsig2(i)=jsig2(i)+(dens2_jack(i,ijack)-jmed2(i))**2
                jsig3(i)=jsig3(i)+(dens3_jack(i,ijack)-jmed3(i))**2
                jsig4(i)=jsig4(i)+(dens4_jack(i,ijack)-jmed4(i))**2
                enddo
                sig1(i)=sqrt(((npix-1)/real(npix))*jsig1(i))
                sig2(i)=sqrt(((npix-1)/real(npix))*jsig2(i))
                sig3(i)=sqrt(((npix-1)/real(npix))*jsig3(i))
                sig4(i)=sqrt(((npix-1)/real(npix))*jsig4(i))
     !          write(10,*) rad, dens2_1(1,i),dens2_2(1,i), &
     !                  dens2_3(1,i),dens2_4(1,i)
     !          write(10,*) rad, dens2_1(1,i)/(f1/real(ftot)),dens2_2(1,i)/(f2/real(ftot)), &
     !                  dens2_3(1,i)/(f3/real(ftot)),dens2_4(1,i)/(f4/real(ftot))
 enddo

j=0
distlog = (0+ 0.5)*dbindist_vu + lrmin_vu
rad    = (10.0**distlog)         !*rv
rad0 = rmin_vu  
densidad1 = 0 
densidad2 = 0 
densidad3 = 0 
densidad4 = 0 
do i=1,nbindist
       j=j+1
       distlog = (i+ 0.5)*dbindist_vu + lrmin_vu
       rad    = (10.0**distlog)         !*rv
       densidad1 = densidad1 + dens2_1(1,i)*rad
       densidad2 = densidad2 + dens2_2(1,i)*rad
       densidad3 = densidad3 + dens2_3(1,i)*rad
       densidad4 = densidad4 + dens2_4(1,i)*rad
       o = o + 1

       error1 = error1 + jsig1(i)**2
       error2 = error2 + jsig2(i)**2
       error3 = error3 + jsig3(i)**2
       error4 = error4 + jsig4(i)**2
       if (mod(i,20)==0) then
               !distlog = (i - 0.5)*dbindist_vu + lrmin_vu
               !rad    = (10.0**distlog)         !*rv
               promedio1 = 1./(rad-rad0)*densidad1*log(10.)*(log10(rad)-log10(rad0))/real(o)
               promedio2 = 1./(rad-rad0)*densidad2*log(10.)*(log10(rad)-log10(rad0))/real(o)
               promedio3 = 1./(rad-rad0)*densidad3*log(10.)*(log10(rad)-log10(rad0))/real(o)
               promedio4 = 1./(rad-rad0)*densidad4*log(10.)*(log10(rad)-log10(rad0))/real(o)
               error1 = sqrt(error1)/real(o)
               error2 = sqrt(error2)/real(o)
               error3 = sqrt(error3)/real(o)
               error4 = sqrt(error4)/real(o)
               rad0 = rad 
               densidad1 = 0   
               densidad2 = 0   
               densidad3 = 0   
               densidad4 = 0   
               o = 0             
       write(10,*) rad, promedio1,error1,promedio2,error2,promedio3,error3,promedio4,error4
               error1 = 0
               error2 = 0
               error3 = 0
               error4 = 0
       endif
 enddo
 close(10)
  deallocate(   dens2,dens_jack)
  deallocate(dens1_jack,dens2_jack,dens3_jack,dens4_jack)
  deallocate(dens2_1,dens2_2,dens2_3,dens2_4)
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
       real*4,intent(in)    :: hsml
       real*8,intent(in)    :: dist
       real(pr),intent(inout) :: fac
       real(pr)           ::  u,hinv,hinv3
       real(pr)           ::  NORM

       NORM   = 8.0/3.1415

          hinv=1./(hsml)      !hsml es el la dist a la ultima particula, en la
          hinv3=(1./(hsml))**3.
          u=dist*(hinv)
       if(dist>hsml) then 
          fac=0.0
       else
                                !formula entra q=d/(2*h), osea el diametro.
          !hinv3=1.0 ! NO damos en Kernel como Volumen
          if(u<=0.5) then 
             fac=hinv3*(1.0 + 6.0 * ( u - 1.0)*u*u)*NORM 
          else
             fac=hinv3 * 2.0* (1.-u) * (1.-u) * (1.-u) *NORM
          endif
       endif
     end subroutine KernelGadget
end program cicjack_omp
subroutine Kernel(r,h,w)
        implicit none
        real,parameter :: pi=acos(-1.)
        real*8 :: r,w,q
        real :: h
        q=r/(2.*h)
         if (q >= 0 .and. q <= .5) then
                  w = (8./pi)*(1.-6.*q**2+6.*q**3)
         elseif (q >.5 .and. q<=1. ) then
                  w = (8./pi)*(2.*(1.-q)**3)
          else 
                  w = 0
         endif
endsubroutine
