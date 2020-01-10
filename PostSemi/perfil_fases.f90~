! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=20, VOID=1373
        real,parameter :: rmax=30, pi=acos(-1.),rmin=.5, dth=116.24, Tc=1e5
        real :: abin,d,vol,vol2,rm,rand,minmass,maxmass,minsfr,maxsfr,abin2,ri,rad,r0
        integer :: bin,bin2
        real ::   dgs,ddm,ddm2,dst,dtt,x,y,z,xbox,ybox,zbox,rv,rho,vr,tita,fi,vx0,vy0,vz0 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm,mass_tot,mass_gs,mass_dm2,mass_st
        integer,dimension(bines):: num_gs1,num_gs2,num_gs3,num_gs4
        real,dimension(bines) :: vel_dm, vel_gs
        real :: xh,yhe,mu0,mp,kcgs,vv,te
        !-----------------------------------------------------------------------------     
        call reader()

      !  xbox=411.217
      !  ybox=162.1655
      !  zbox=453.0553 
        xbox=403.896
        ybox=459.8882
        zbox=440.9021 

        open(13,file='/mnt/is2/fstasys/ITV/base09/voids/voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo

        read(13,*) rv,x,y,z
        x=x-xbox+250
        y=y-ybox+250
        z=z-zbox+250
        close(13)
        write(*,*) x,y,z
        
print*, rmax
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
print*, rm
        ri=log10(rmin)
        abin = (rm-ri)/real(bines)    
     print*, 10**(bines*abin) 
  stop   
        sf = 0
        num_gs1  = 0
        num_gs2  = 0
        num_gs3  = 0
        num_gs4  = 0
 !----------------------------------------------- GAS --------------------------       
        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1e10
        do i=1,nall(0)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm .and. d>ri) then
                        mu0=(1.0-yHe)/(1+yHe+ne(i))
                        te=(5./3.-1.)*u(i)*vv*mu0*mp/kcgs
                        bin = int((d-ri)/abin) + 1
                        if (dens(i)<dth .and. te < Tc) then
                                num_gs1(bin)  = num_gs1(bin) + 1
                        elseif (dens(i)<dth .and. te >= Tc) then
                                num_gs2(bin)  = num_gs2(bin) + 1
                        elseif (dens(i)>=dth .and. te < Tc) then
                                num_gs3(bin)  = num_gs3(bin) + 1
                        elseif (dens(i)<dth .and. te < Tc) then
                                num_gs4(bin)  = num_gs4(bin) + 1
                        endif        
                endif

        enddo
   !---------------------------------------------------------------------     
        open(10,file='particlesprofile_fases.dat',status='unknown')
        dgs = 0
        r0 = rmin
        do i=1,bines
                dgs = dgs + mass_gs(i)
                rad = 10**(i*abin)+rmin
                print*, rad
                vol2 = (4./3.)*pi*((10**rad)**3-(10**r0)**3)
                vol  = (4./3.)*pi*((10**rad)**3-(10**ri)**3)
              write(10,*) 10**((i-.5)*abin+ri),num_gs1(i)/vol2, num_gs2(i)/vol2,num_gs3(i)/vol2,num_gs4(i)/vol2
                  r0=rad       
        enddo
        close(10)
        !print*, sf
endprogram
