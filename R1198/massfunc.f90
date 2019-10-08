! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 
program first
        use modulos
        implicit none
        integer,parameter :: bines=19
        real,parameter :: rmax=40, pi=acos(-1.)
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2
        integer :: bin,bin2
        integer :: dmp,gsp 
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth
        integer,dimension(bines):: num_gs,num_dm,num_st
  !-----------------------------------------------------------------------------     
        call reader()
  !-----------------------------------------------------------------------------
        rm=log10(rmax)
        lim1=log10(6.)
        lim2=log10(12.)
        lim3=log10(18.)
        abin = rm/real(bines)        
        !abin2= rmax/real(bines)
 !----------------------------------------------- GAS --------------------------       
        do i=1,nall(0)
                d=log10(sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2))
                if (d < limit1) then
                        bin = int(d/abin) + 1
                        mf1(bin) = sf(bin) + sfr(i)
                        if (bin==1 .and. sfr(i)>0) print*, sfr(i)
                endif
        enddo
  !---------------------------------------------- DM ------------------------
endprogram
