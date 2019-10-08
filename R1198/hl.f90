program asd
        implicit none
        integer,parameter :: halos=64439, bines=18,VOID=1373
        real, parameter  :: lim1=6, lim2=12, lim3=18

        integer :: i,bin1,bin2,bin
        real :: abin1,abin2,abin3,d,xbox,ybox,zbox,l1,l2,l3
        real :: x,y,z,rvir,np,npgs,npdm,npst,nptd,mgs,mdm,mst,mtd,lgs,ldm
        real ::ssfr,e,h,hsml,lx,ly,lz,spin,maxdm,x0,y0,z0,rv
        real,dimension(bines) ::mass_dm1,mass_dm2,mass_dm3
        
        xbox=403.8960 
        ybox=459.8882
        zbox= 440.9021 
        open(13,file='voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo
        read(13,*) rv,x,y,z
        x=x-xbox+250
        y=y-ybox+250
        z=z-zbox+250
        
        abin1=log10(34e3)/real(bines) ! 2.8 es el halo de maxima masa
        maxdm=0
        open(10,file='halosprop_2rvir.dat',status='unknown')
        do i=1,halos
                read(10,*) x0,y0,z0,rvir,np,npgs,npdm,npst,nptd,mgs,mdm,mst,lgs,ldm,&
                        ssfr,e,h,hsml,lx,ly,lz,spin
                d=log10(sqrt((x0-x)**2+(y0-y)**2+(z0-z)**2))
                if (d<=log10(lim1)) then
                        bin=int(log10(mdm)/abin1)+1
                        mass_dm1(bin)=mass_dm1(bin) + 1
                elseif (d >log10(lim1) .and. d<=log10(lim2)) then
                        bin=int(log10(mdm)/abin1)+1
                        mass_dm2(bin)=mass_dm2(bin) + 1
                elseif (d>log10(lim2) .and. d<log10(lim3)) then 
                        bin=int(log10(mdm)/abin1)+1
                        mass_dm3(bin)=mass_dm3(bin) + 1
                endif
                 if (mdm >= maxdm) maxdm=mdm             
        enddo
print*, maxdm
        close(10)
        open(11,file='haloes_mf.dat',status='unknown')
        do i=1,bines
                 write(11,*) 10**(i*abin1), mass_dm1(i),mass_dm2(i),mass_dm3(i)
        enddo
        close(11)
        
endprogram
