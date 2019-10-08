! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=19
        real,parameter :: rmax=40, pi=acos(-1.)
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2
        integer :: bin,bin2
        integer :: dgs,ddm,ddm2,dst,dtt 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm,mass_tot,mass_gs,mass_dm2,mass_st
        integer,dimension(bines):: num_gs,num_dm,num_st
  !-----------------------------------------------------------------------------     
        call reader()
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
        
        abin = rm/real(bines)        
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
 !----------------------------------------------- GAS --------------------------       
        do i=1,nall(0)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                        bin = int(d/abin) + 1
                        sf(bin) = sf(bin) + sfr(i)
                        if (bin==1 .and. sfr(i)>0) print*, sfr(i)
                        num_gs(bin)  = num_gs(bin) + 1
                        wgth(bin) = wgth(bin) + 1e10*mass(i)
                        mass_tot(bin) = mass_tot(bin) + mass(i)
                        mass_gs(bin) = mass_gs(bin) + mass(i)
                endif

        enddo
  !---------------------------------------------- DM ------------------------
        do i=nall(0)+1,nall(0)+nall(1)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                bin = int((d)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                mass_dm(bin) = mass_dm(bin) + mass(i)
                mass_dm2(bin) = mass_dm2(bin) + mass(i)
                mass_tot(bin) = mass_tot(bin) + mass(i)
                endif
        enddo
  !---------------------------------------------- TIDALES ------------------   
        do i=nall(0)+nall(1)+1,nall(0)+nall(1)+nall(2)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                bin = int((d)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                mass_dm2(bin) = mass_dm2(bin) + mass(i)
                mass_tot(bin) = mass_tot(bin) + mass(i)
                endif
        enddo   
   !--------------------------------------------- ESTRELLAS -------------
        do i=nall(0)+nall(1)+nall(2)+1,nall(0)+nall(1)+nall(2)+nall(4)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if ((d) < rm) then
                bin = int((d)/abin) + 1
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
        do i=1,bines
                ddm = ddm + mass_dm(i)
                ddm2 = ddm2 + mass_dm2(i)
                dgs = dgs + mass_gs(i)
                dst = dst + mass_st(i)
                dtt = dtt + mass_tot(i)
              !  vol = (4./3.)*pi*((10**(i*abin))**3-(10**((i-1)*abin))**3)
                vol = (4./3.)*pi*((10**(i*abin))**3)
       !        write(10,*) 10**(i*abin), num_gs(i)/vol, num_dm(i)/vol,num_st(i)/vol, wgth(i)/vol, sf(i)/vol, mass_dm(i)
                write(10,*) 10**(i*abin) - (10**abin)/2, dgs,ddm,ddm2,dst,dtt,vol
        enddo
        close(10)
        !print*, sf
endprogram
