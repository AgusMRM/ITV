! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=100, VOID=1373
        real,parameter :: rmax=25, pi=acos(-1.),rmin=.2
        real :: abin,d,vol,vol2,rm,rand,minmass,maxmass,minsfr,maxsfr,abin2,ri,rad,r0
        integer :: bin,bin2
        real :: dgs,ddm,ddm2,dst,dtt,x,y,z,xbox,ybox,zbox,rv 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm,mass_tot,mass_gs,mass_dm2,mass_st
        integer,dimension(bines):: num_gs,num_dm,num_st
  !-----------------------------------------------------------------------------     
        call reader()

      !  xbox=411.217
      !  ybox=162.1655
      !  zbox=453.0553 
        xbox=403.896
        ybox=459.8882
        zbox=440.9021 

        open(13,file='voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo

        read(13,*) rv,x,y,z
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
        num_dm  = 0
        mass_dm=0
        mass_dm2=0
  !---------------------------------------------- DM ------------------------
        do i=nall(0)+1,nall(0)+nall(1)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm .and. d>ri) then
                bin = int((d-ri)/abin) + 1
                mass_dm(bin) = mass_dm(bin) + mass(i)
                num_dm(bin)  = num_dm(bin) + 1
                
                rad = (i*abin+ri)
                vol2 = (4./3.)*pi*((10**rad)**3-(10**r0)**3)


                endif
        enddo
   !---------------------------------------------------------------------     
        open(10,file='particlesprofile.dat',status='unknown')
        ddm = 0
        ddm2 = 0
                r0 = ri
        do i=1,bines
                ddm = ddm + mass_dm(i)
                ddm2 = ddm2 + mass_dm2(i)
                rad = (i*abin+ri)
                vol2 = (4./3.)*pi*((10**rad)**3-(10**r0)**3)
              write(10,*) 10**((i-.5)*abin+ri), num_dm(i), vol2, sqrt(real(num_dm(i)))
                  r0=rad       
        enddo
        close(10)
        !print*, sf
endprogram
