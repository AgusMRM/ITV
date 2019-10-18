program first0
        implicit none
        integer,parameter :: bines=30, VOID=1198,halos= 105390  !106409 
        real,parameter :: rmax=35,rmin=.1, pi=acos(-1.),corte=20
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2,rad,r0
        integer :: bin,bin2  
        integer :: dmp,gsp,i 
        real :: x,y,z,rv,ri,vol2,masa,a
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm
        integer,dimension(bines):: num_gs,num_dm,num_st
        real,dimension(3,halos) :: pos
        integer,dimension(halos) :: pnum
  !-----------------------------------------------------------------------------    
 !open(24,file='/mnt/is2/fstasys/ITV/base09/rockstar/halos_0.5.ascii',status='unknown',action='read')
 open(24,file='/home/arodriguez/halos.dat',status='unknown',action='read')
 do i=1,20
        read(24,*)
 enddo
 do i=1,halos
        read(24,*) a,pnum(i),a,a,a,a,a,a, pos(1,i), pos(2,i), pos(3,i)
 enddo
 close(24) 
 !------------------------------------------------------------------------------
print*, 'leo voids'
        open(13,file='voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo

        read(13,*) rv,x,y,z
        close(13)
print*, rv,x,y,z
  !-----------------------------------------------------------------------------
        abin = (rmax-rmin)/real(bines)        
        num_dm  = 0
  !---------------------------------------------- DM ------------------------
        do i=1,halos
                d=sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2)
                if (d < rmax .and. d>rmin .and. pnum(i)>= corte) then
                bin = int((d-rmin)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                endif
        enddo
   !---------------------------------------------------------------------     
        open(17,file='halosprofiles.dat',status='unknown')
        dmp = 0
        r0=rmin
        ri=0
        do i=1,bines
                dmp = dmp + num_dm(i)
                rad = i*abin+rmin
              !  vol2 = (4./3.)*pi*((10**(rad))**3-((10**(r0))**3))  ! ESTE ES
         !       PARA EL DIFERENCIAL
                vol  = (4./3.)*pi*((i*abin+rmin)**3-(rmin)**3)
                !write(17,*) ((i-.5)*abin+rmin), dmp, vol,num_dm(i),vol2
                write(17,*) ((i-0)*abin+rmin), dmp, vol,num_dm(i),vol2
                r0 = rad
        enddo
        close(17)
endprogram
