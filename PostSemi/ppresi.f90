! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=20, VOID=1198
        real,parameter :: rmax=30, pi=acos(-1.),rmin=.5
        real :: abin,d,vol,vol2,rm,rand,minmass,maxmass,minsfr,maxsfr,abin2,ri,rad,r0
        integer :: bin,bin2
        real ::   dgs,ddm,ddm2,dst,dtt,x,y,z,xbox,ybox,zbox,rv,rho,vr,tita,fi,vx0,vy0,vz0 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm,mass_tot,mass_gs,mass_dm2,mass_st
        integer,dimension(bines):: num_gs,num_dm,num_st
        real,dimension(bines) :: vel_dm, vel_gs
        !-----------------------------------------------------------------------------     
        call reader()

        xbox=411.217
        ybox=162.1655
        zbox=453.0553 
     !   xbox=403.896
     !   ybox=459.8882
     !   zbox=440.9021 

        open(13,file='/mnt/is2/fstasys/ITV/base09/voids/voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo

        read(13,*) rv,x,y,z,vx0,vy0,vz0
        x=x-xbox+250
        y=y-ybox+250
        z=z-zbox+250
        close(13)
        write(*,*) x,y,z
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
        ri=log10(rmin)
        abin = (rm-ri)/real(bines)        
        sf = 0
        num_gs  = 0
        num_dm  = 0
        num_st = 0
        wgth = 0
        minsfr=400
        maxsfr=0
        mass_tot=0
        mass_gs=0
        mass_dm=0
        mass_dm2=0
        mass_st=0
        vel_dm=0
        vel_gs=0
 !----------------------------------------------- GAS --------------------------       
        do i=1,nall(0)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm .and. d>ri) then
                        bin = int((d-ri)/abin) + 1
                        sf(bin) = sf(bin) + sfr(i)
                        if (bin==1 .and. sfr(i)>0) print*, sfr(i)
                        num_gs(bin)  = num_gs(bin) + 1
                      !  wgth(bin) = wgth(bin) + 1e10*mass(i)
                      !  mass_tot(bin) = mass_tot(bin) + mass(i)
                        mass_gs(bin) = mass_gs(bin) + mass(i)
                        call titaangle(pos(1,i),pos(2,i),pos(3,i),tita)
                        call fiangle(pos(1,i),pos(2,i),tita)
                        vr = (vel(1,i)-vx0)*sin(tita)*cos(fi) + &
                             (vel(2,i)-vy0)*sin(tita)*sin(fi) + (vel(3,i)-vz0)*cos(tita)
                        vel_gs(bin) = vel_gs(bin) + vr
                endif

        enddo
  !---------------------------------------------- DM ------------------------
        do i=nall(0)+1,nall(0)+nall(1)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm .and. d>ri) then
                bin = int((d-ri)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                mass_dm(bin) = mass_dm(bin) + mass(i)
                mass_dm2(bin) = mass_dm2(bin) + mass(i)
                mass_tot(bin) = mass_tot(bin) + mass(i)
                        call titaangle(pos(1,i),pos(2,i),pos(3,i),tita)
                        call fiangle(pos(1,i),pos(2,i),tita)
                        vr = (vel(1,i)-vx0)*sin(tita)*cos(fi) + &
                             (vel(2,i)-vy0)*sin(tita)*sin(fi) + (vel(3,i)-vz0)*cos(tita)
                        vel_dm(bin) = vel_dm(bin) + vr
                endif
        enddo
  !---------------------------------------------- TIDALES ------------------   
        do i=nall(0)+nall(1)+1,nall(0)+nall(1)+nall(2)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm) then
                bin = int((d-ri)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                mass_dm2(bin) = mass_dm2(bin) + mass(i)
                mass_tot(bin) = mass_tot(bin) + mass(i)
                endif
        enddo   
   !--------------------------------------------- ESTRELLAS -------------
        do i=nall(0)+nall(1)+nall(2)+1,nall(0)+nall(1)+nall(2)+nall(4)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm .and. d>ri) then
                bin = int((d-ri)/abin) + 1
                num_st(bin)  = num_st(bin) + 1
                mass_tot(bin) = mass_tot(bin) + mass(i)
                mass_st(bin) = mass_st(bin) + mass(i)
                endif
        enddo
   !---------------------------------------------------------------------     
        open(10,file='particlesprofile.dat',status='unknown')
        ddm = 0
        ddm2 = 0
        dgs = 0
        dst = 0
        dtt = 0
                r0 = ri
        do i=1,bines
                rho = mass_dm(i) + mass_gs(i) + mass_st(i)
                ddm = ddm + mass_dm(i)
                dgs = dgs + mass_gs(i)
                dst = dst + mass_st(i)
                dtt = dtt + mass_tot(i)
                rad = (i*abin+ri)
                vol2 = (4./3.)*pi*((10**rad)**3-(10**r0)**3)
                vol  = (4./3.)*pi*((10**rad)**3-(10**ri)**3)
       !        write(10,*) 10**(i*abin), num_gs(i)/vol, num_dm(i)/vol,num_st(i)/vol, wgth(i)/vol, sf(i)/vol, mass_dm(i)
             !   write(10,'(11(e25.10,1x))') 10**((i-.5)*abin+ri), dgs,ddm,ddm2,dst,dtt,vol, & 
             !                  mass_dm(i), mass_gs(i), mass_st(i),vol2,sf(i)
              !write(10,*) 10**((i-.5)*abin+ri), mass_gs(i), mass_dm(i),mass_st(i), sf(i), vol2,vol
              write(10,*) 10**((i-.5)*abin+ri), dgs, ddm ,mass_st(i), sf(i), vol2,vol,rho/vol2, &
!                      vlx_gs(i)/num_gs(i),vly_gs(i)/num_gs(i),vlz_gs(i)/num_gs(i),vlx_dm(i)/num_dm(i), &
!                      vly_dm(i)/num_dm(i),vlz_dm(i)/num_dm(i)
                      vel_gs(i)/num_dm(i),vel_dm(i)/num_dm(i)
                  r0=rad       
        enddo
        close(10)
        !print*, sf
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
