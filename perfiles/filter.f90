! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        real,parameter    :: masadm=0.0932880491
        integer,parameter :: bines=10, VOID=1373, n=20
        real,parameter :: rmax=15, pi=acos(-1.),rmin=.3
        real :: abin,d,vol,vol2,rm,minmass,maxmass,minsfr,maxsfr,abin2,ri,rad,r0
        integer :: bin,bin2,ijac,ijack,nparticulas
        real :: dgs,ddm,ddm2,dst,dtt,x,y,z,xbox,ybox,zbox,rv,rand 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm,mass_tot,mass_gs,mass_dm2,mass_st
        integer,allocatable :: jack(:)
        integer :: num_dm(bines,n)
        real,allocatable :: part(:,:)
        real,dimension(bines) :: jmed,jsig,sig
        logical,dimension(n):: ilogic
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

open(10,file='insidevoid_gas.dat',status='unknown')
open(11,file='insidevoid_dm.dat',status='unknown')
        mass_dm=0
        call reader()
       do i=1,nall(0)
        d=(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
        if (d < rmax) then
                write(10,*)  pos(1,i),pos(2,i),pos(3,i),vel(1,i),vel(2,i),vel(3,i),RHO(i),hsml(i),u(i)
        endif
       enddo
       
       do i=1+nall(0),nall(0)+nall(1)
        d=(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
        if (d < rmax) then
                write(11,*)  pos(1,i),pos(2,i),pos(3,i),vel(1,i),vel(2,i),vel(3,i)
        endif
        enddo 
       close(11)
       close(10)
 endprogram  
