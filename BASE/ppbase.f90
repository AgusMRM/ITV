! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=19, VOID=1198
        real,parameter :: rmax=40, pi=acos(-1.)
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2
        integer :: bin,bin2
        integer :: dmp,gsp 
        real :: x,y,z,rv
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

        read(13,*) rv,x,y,z
        close(13)
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
        

        abin = rm/real(bines)        
        sf = 0
        num_dm  = 0
  !---------------------------------------------- DM ------------------------
        do i=nall(0)+1,nall(0)+nall(1)
                d=log10(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm) then
                bin = int((d)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                endif
        enddo
   !---------------------------------------------------------------------     
        open(10,file='particlesprofileBASE.dat',status='unknown')
        dmp = 0
        do i=1,bines
                dmp = dmp + num_dm(i)
                vol = (4./3.)*pi*((10**(i*abin))**3-(10**((i-1)*abin))**3)
!                vol = (4./3.)*pi*((10**(i*abin))**3)
                write(10,*) 10**(i*abin) - (10**abin)/2, dmp, vol
        enddo
        close(10)
endprogram
