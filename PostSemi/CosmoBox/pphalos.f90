! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modV.f90 ppresi.f90 
program first
        implicit none
        integer,parameter :: bines=100, halos=87780-19,pcorte=1231,VOID=1198
       ! integer,parameter :: bines=100, halos=71417-19,pcorte=2*615,VOID=1373
        real,parameter :: rmax=30, pi=acos(-1.),rmin=0.001
        real :: abin,d,vol,vol2,rm,rand,minmass,maxmass,minsfr,maxsfr,abin2,ri,rad,r0
        integer :: bin,bin2,i
        real :: dgs,ddm,ddm2,dst,dtt,x,y,z,xbox,ybox,zbox,rv ,a
        real,dimension(bines):: sf,hsml0
        real*8,dimension(bines):: wgth,mass_dm,mass_tot,mass_gs,mass_dm2,mass_st
        integer,dimension(bines):: num_gs,num_dm,num_st
        real, allocatable :: pos(:,:)
        integer, allocatable :: nump(:)
  !-----------------------------------------------------------------------------     

        allocate(pos(3,halos),nump(halos))      

       xbox=411.217
       ybox=162.1655
       zbox=453.0553 
       !  xbox=403.896
       !  ybox=459.8882
       !  zbox=440.9021 

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
        open(14,file='/mnt/is2/dpaz/ITV/R1198/halos/halos_50.ascii')
      !  open(14,file='/mnt/is2/dpaz/ITV/S1373/halos/halos_50.ascii')
        do i=1,19
                read(14,*) 
        enddo

        do i=1,halos
                read(14,*) a, nump(i), a,a,a,a,a,a,pos(1,i), pos(2,i), pos(3,i) 

        enddo
        close(14)
  !-----------------------------------------------------------------------------
        rm=(rmax)
        ri=(rmin)
        abin = (rm-ri)/real(bines)        
        num_dm  = 0
  !---------------------------------------------- DM ------------------------
        do i=1,halos
                d=(sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2))
                if ((d) < rm .and. d>ri .and. nump(i) > pcorte) then
                bin = int((d-ri)/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
               ! mass_dm(bin) = mass_dm(bin) + mass(i)
               ! mass_dm2(bin) = mass_dm2(bin) + mass(i)
               ! mass_tot(bin) = mass_tot(bin) + mass(i)
                endif
        enddo
   !---------------------------------------------------------------------     
        open(10,file='halosprofile.dat',status='unknown')
        ddm = 0
        r0 = rmax
        do i=1,bines
                ddm = ddm + num_dm(i)
                rad = (i*abin+ri)
                vol2 = (4./3.)*pi*((rad)**3-(rmin)**3)
       !        write(10,*) 10**(i*abin), num_gs(i)/vol, num_dm(i)/vol,num_st(i)/vol, wgth(i)/vol, sf(i)/vol, mass_dm(i)
             !   write(10,'(11(e25.10,1x))') 10**((i-.5)*abin+ri), dgs,ddm,ddm2,dst,dtt,vol, & 
             !                  mass_dm(i), mass_gs(i), mass_st(i),vol2,sf(i)
              !write(10,*) 10**((i-.5)*abin+ri), mass_gs(i), mass_dm(i),mass_st(i), sf(i), vol2,vol
              write(10,*) (i-0.5)*abin + rmin, ddm , vol2
                  !r0=rad       
        enddo
        close(10)
        !print*, sf
endprogram
