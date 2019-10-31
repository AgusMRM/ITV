! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=30, VOID=1198
        real,parameter :: rmax=25, pi=acos(-1.),rmin=.5,dlim=10**2,tlim=10**5
        real :: abin,d,vol,vol2,rm,rand,minmass,maxmass,minsfr,maxsfr,abin2,ri,rad,r0
        integer :: bin,bin2
        real :: dgs,ddm,ddm2,dst,dtt,x,y,z,xbox,ybox,zbox,rv,dmean,temp 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm,mass_tot,mass_gs,mass_dm2,mass_st
        integer,dimension(bines):: num1_gs,num2_gs,num3_gs,num4_gs
  !-----------------------------------------------------------------------------     
  dmean=(3*(100**2)*(0.045))/(8*pi*(4.3e-9))*1e-10
        call reader()
print*, u(200), u(12121212)
        xbox=411.217
        ybox=162.1655
        zbox=453.0553 
   !     xbox=403.896
   !     ybox=459.8882
   !     zbox=440.9021 

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
        !x=270
        !y=270
        !z=270
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
        ri=log10(rmin)
        abin = (rm-ri)/real(bines)        
        num1_gs  = 0
        num2_gs  = 0
        num3_gs  = 0
        num4_gs  = 0
 !----------------------------------------------- GAS --------------------------       
        do i=1,nall(0)
                call temperature(u(i),ne(i),temp)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                !if (d < rm .and. d>ri) then
                        bin = int((d-ri)/abin) + 1
                        if (d < rm .and. d>ri .and. temp < tlim .and. (rho(i)/dmean)<dlim) then
                        num1_gs(bin)  = num1_gs(bin) + 1
                        elseif (d < rm .and. d>ri .and. temp<tlim .and. (rho(i)/dmean)>=dlim) then
                        num2_gs(bin)  = num2_gs(bin) + 1
                        elseif (d < rm .and. d>ri .and. temp>=tlim .and. (rho(i)/dmean)<dlim) then
                        num3_gs(bin)  = num3_gs(bin) + 1
                        elseif (d < rm .and. d>ri .and. temp>=tlim .and. (rho(i)/dmean)>=dlim) then
                        num4_gs(bin)  = num4_gs(bin) + 1
                        endif
                !endif
        enddo
        print*, dmean
        print*, sum(num1_gs)
        print*, sum(num2_gs)
        print*, sum(num3_gs)
        print*, sum(num4_gs)
   !---------------------------------------------------------------------     
        open(10,file='particlesprofile.dat',status='unknown')
        dgs = 0
        r0 = ri
        do i=1,bines
                rad = (i*abin+ri)
                vol = (4./3.)*pi*((10**rad)**3-(10**r0)**3)
              write(10,*) 10**((i-.5)*abin+ri),num1_gs(i), num2_gs(i), &
                     num3_gs(i), num4_gs(i), vol 
                  r0=rad       
        enddo
        close(10)
        !print*, sf
endprogram
