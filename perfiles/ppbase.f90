! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first0
        use modulos
        implicit none
        integer,parameter :: bines=20, VOID=1277
        real,parameter :: rmax=35,rmin=.5, pi=acos(-1.)
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2,rad,r0
        integer :: bin,bin2
        integer :: dmp,gsp 
        real :: x,y,z,rv,ri,vol2,vx,vy,vz
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm
        integer,dimension(bines):: num_gs,num_dm,num_st
  !-----------------------------------------------------------------------------     
        call reader()
print*, 'leo voids'
        open(13,file='voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo

        read(13,*) rv,x,y,z,vx,vy,vz
        close(13)
print*, rv,x,y,z
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
        ri=log10(rmin)

        abin = (rm-ri)/real(bines)        
        sf = 0
        num_dm  = 0
  !---------------------------------------------- DM ------------------------
        do i=nall(0)+1,nall(0)+nall(1)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm .and. d>ri) then
                bin = int((d-ri)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                endif
        enddo
   !---------------------------------------------------------------------     
        open(10,file='particlesprofileBASE.dat',status='unknown')
        dmp = 0
        r0=ri
        do i=1,bines
                dmp = dmp + num_dm(i)
                rad = i*abin+ri
                vol2 = (4./3.)*pi*((10**(rad))**3-((10**(r0))**3))
                vol  = (4./3.)*pi*((10**(i*abin+ri))**3-(10**ri)**3)
                write(10,*) 10**((i-.5)*abin+ri), dmp, vol,num_dm(i),vol2
                r0 = rad
        enddo
        close(10)
endprogram
subroutine titaangle(x,y,z,tita)
implicit none
real :: x,y,z,tita,pi
pi=atan(-1.)
if (z>0) then
        tita=atan(sqrt(x**2+y**2)/z)
elseif (z==0) then
        tita=pi/2.
else    
        tita=pi+atan(sqrt(x**2+y**2)/x)
endif        
endsubroutine        
subroutine  fiangle(x,y,fi)
implicit none
real :: x,y,z,fi,pi
pi=atan(-1.) 
if (x>0 .and. y>0) then
        fi=atan(y/x)
elseif (x>0 .and. y<0) then        
        fi=2*pi+atan(y/x)
elseif (x==0) then
        fi=pi/2.*sign(1.,y)
else
        fi=pi+atan(y/x)
endif
endsubroutine        
