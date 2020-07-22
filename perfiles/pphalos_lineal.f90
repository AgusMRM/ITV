program first0
        implicit none
        integer,parameter :: bines=30, VOID= 443,halos=172548 -19  !106409 
        real,parameter :: rmax=3.5,rmin=.1, pi=acos(-1.),corte=100
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2,rad,r0
        integer :: bin,bin2  
        integer :: dmp,gsp,i 
        real :: x,y,z,rv,ri,vol2,masa,a,xc,yc,zc,xbox,ybox,zbox
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm
        integer,dimension(bines):: num_gs,num_dm,num_st
        real,dimension(3,halos) :: pos
        integer,dimension(halos) :: pnum
  !-----------------------------------------------------------------------------    
 !open(24,file='/mnt/is2/fstasys/ITV/base09/rockstar/halos_0.5.ascii',status='unknown',action='read')
 open(24,file='/home/arodriguez/Doctorado/snap50/halos_ 443.ascii',status='unknown',action='read')
 do i=1,19
        read(24,*)
 enddo
 do i=1,halos
        read(24,*) a,pnum(i),a,a,a,a,a,a, pos(1,i), pos(2,i), pos(3,i)
 enddo
 close(24) 

 open(13,file='/home/arodriguez/Doctorado/itv/voids_boxes.dat',status='unknown')

!   open(13,file='/mnt/is2/fstasys/ITV/base09/voids/voids_new.dat')
   do i=1,VOID-1
   read(13,*)
   enddo
   
   read(13,*) rv,x,y,z,xbox,ybox,zbox
   xc=x-xbox+250
   yc=y-ybox+250
   zc=z-zbox+250
   close(13)
   write(*,*) x,y,z
  !-----------------------------------------------------------------------------
        abin = ((rmax*rv)-(rmin*rv))/real(bines)        
        num_dm  = 0
  !---------------------------------------------- DM ------------------------
        do i=1,halos
                d=sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)
                if (d < rmax*rv .and. d>rmin*rv .and. pnum(i)>= corte) then
                bin = int((d-rmin)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                endif
        enddo
   !---------------------------------------------------------------------     
        open(17,file='halosprofiles.dat',status='unknown')
        dmp = 0
        r0=rmin*rv
        ri=0
        do i=1,bines
                dmp = dmp + num_dm(i)
                rad=i*abin+(rmin*rv)
                vol  = (4./3.)*pi*(rad**3-(rmin*rv)**3)
                vol2 = (4./3.)*pi*(rad**3-r0**3)
                write(17,*) ((i-.5)*abin+rmin*rv)/rv, dmp, vol,num_dm(i),vol2
                r0 = rad
        enddo
        close(17)
endprogram
