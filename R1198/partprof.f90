! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=19
        real,parameter :: rmax=200, pi=acos(-1.)
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2
        integer :: bin,bin2
        integer :: dmp,gsp 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm
        integer,dimension(bines):: num_gs,num_dm,num_st
  !-----------------------------------------------------------------------------     
        call reader()
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
        
        abin = rm/real(bines)        
        sf = 0
        num_gs  = 0
        num_dm  = 0
        wgth = 0
        minsfr=400
        maxsfr=0
 !----------------------------------------------- GAS --------------------------       
        do i=1,nall(0)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                        bin = int(d/abin) + 1
                        sf(bin) = sf(bin) + sfr(i)
                        if (bin==1 .and. sfr(i)>0) print*, sfr(i)
                        num_gs(bin)  = num_gs(bin) + 1
                        wgth(bin) = wgth(bin) + 1e10*mass(i)
                endif

        enddo
  !---------------------------------------------- DM ------------------------
        do i=nall(0)+1,nall(0)+nall(1)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                bin = int((d)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                mass_dm(bin) = mass_dm(bin) + mass(i)
                endif
        enddo
  !---------------------------------------------- TIDALES ------------------   
        do i=nall(0)+nall(1)+1,nall(0)+nall(1)+nall(2)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                bin = int((d)/abin) + 1
                mass_dm(bin) = mass_dm(bin) + mass(i)
                endif
        enddo   
   !--------------------------------------------- ESTRELLAS -------------
        do i=nall(0)+nall(1)+nall(2)+1,nall(0)+nall(1)+nall(2)+nall(4)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                bin = int((d)/abin) + 1
                num_st(bin)  = num_st(bin) + 1
                endif
        enddo
   !---------------------------------------------------------------------     
        open(10,file='particlesprofile.dat',status='unknown')
        dmp = 0
        gsp = 0
        do i=1,bines
                dmp = dmp + num_dm(i)
                gsp = gsp + num_gs(i)
                vol = (4./3.)*pi*((10**(i*abin))**3-(10**((i-1)*abin))**3)
        !        write(10,*) 10**(i*abin), num_gs(i)/vol, num_dm(i)/vol,num_st(i)/vol, wgth(i)/vol, sf(i)/vol, mass_dm(i)
!                write(10,*) 10**(i*abin) - (10**abin)/2, dmp, gsp, vol,mass_dm(i)
                write(10,*) 10**(i*abin) - (10**abin)/2, vol, num_dm(i), num_gs(i)
        enddo
        close(10)
        !print*, sf
endprogram
