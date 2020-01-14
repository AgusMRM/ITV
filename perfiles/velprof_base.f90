! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first0
!        use modulos
        implicit none
        integer,parameter :: bines=20, halos=850366-16,VOID=1198
        real,parameter :: rmax=100,rmin=5, pi=acos(-1.)
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2,rad,r0
        integer :: bin,bin2,i
        integer :: dmp,gsp 
        real :: x,y,z,rv,ri,vol2,vx,vy,vz,a
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm,vel_r
        integer,dimension(bines):: num_gs,num_dm,num_st
        real,dimension(3,halos) ::pos,vel
        real :: tita,fi,vr,dx,dy,dz
        integer :: np
  !-----------------------------------------------------------------------------     
!        call reader()
open(13,file='/mnt/is2/fstasys/ITV/base09/rockstar/out_0.list',action='read')
do i=1,16
read(13,*)
enddo
pos=0; vel=0
do i=1,halos
        read(13,*) a,a,a,a,a,a,a,np,x,y,z,vx,vy,vz
                if (np>19) then        
        
                pos(1,i)=x
                pos(2,i)=y
                pos(3,i)=z
                vel(1,i)=vx
                vel(2,i)=vy
                vel(3,i)=vz
        endif 
enddo
close(13)
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
        vel_r  = 0
  !---------------------------------------------- DM ------------------------
        do i=1,halos
                dx=pos(1,i)-x
                dy=pos(2,i)-y
                dz=pos(3,i)-z
                if (abs(dx)>250) dx=500-abs(dx)
                if (abs(dy)>250) dy=abs(dy)-500
                if (abs(dz)>250) dz=500-abs(dz)
                d=log10(sqrt(dx**2+dy**2+dz**2))
                if ((d) < rm .and. d>ri) then
                bin = int((d-ri)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                       ! call titaangle(pos(1,i),pos(2,i),pos(3,i),tita)
                       ! call fiangle(pos(1,i),pos(2,i),fi)
                        call titaangle(dx,dy,dz,tita)
                        call fiangle(dx,dy,fi)
                        vr = (vel(1,i)-vx)*sin(tita)*cos(fi) + &
                             (vel(2,i)-vy)*sin(tita)*sin(fi) + (vel(3,i)-vz)*cos(tita)
                        vel_r(bin) = vel_r(bin) + vr
                endif

        enddo
   !---------------------------------------------------------------------     
        open(10,file='halosprofileBASE.dat',status='unknown')
        dmp = 0
        r0=ri
        do i=1,bines
                dmp = dmp + num_dm(i)
                rad = i*abin+ri
                vol2 = (4./3.)*pi*((10**(rad))**3-((10**(r0))**3))
                vol  = (4./3.)*pi*((10**(i*abin+ri))**3-(10**ri)**3)
                write(10,*) 10**((i-.5)*abin+ri), dmp, &
                        vol,num_dm(i),vol2,vel_r(i)/num_dm(i)
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
