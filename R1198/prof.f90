program asd
        implicit none
        integer,parameter :: halos=64439, bines=12,VOID=1198
        real, parameter  :: lim=40
        integer :: i,bin
        real :: abin,d,xbox,ybox,zbox,l1,l2,l3
        real :: x,y,z,rvir,mgs,mdm,mst,mtd,lgs,ldm
        integer ::np,npgs,npdm,npst,nptd
        real ::ssfr,e,h,hsml,lx,ly,lz,spin,maxdm,x0,y0,z0,rv
        integer,dimension(bines) ::part_gs,part_dm,part_st
        real,dimension(bines) ::mass_gs,mass_dm,mass_st
        real,dimension(bines) ::sfrs,e_gs,h_gs,hsml_gs,l_hl,spin_hl
        
        xbox=411.2170 
        ybox=162.1655
        zbox=453.0553 
        open(13,file='voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo
        read(13,*) rv,x,y,z
        x=x-xbox+250
        y=y-ybox+250
        z=z-zbox+250
        
        abin=log10(lim)/bines
        open(10,file='halosprop_2rvir.dat',status='unknown')
        
        mass_gs = 0
        mass_dm =0
        mass_st =0
        part_gs = 0
        part_dm =0
        part_st =0
        sfrs   =0
        e_gs   =0
        h_gs    =0
        hsml_gs =0
        l_hl   =0
        spin_hl =0
        do i=1,halos
                read(10,*) x0,y0,z0,rvir,np,npgs,npdm,npst,nptd,mgs,mdm,mst,lgs,ldm,&
                        ssfr,e,h,hsml,lx,ly,lz,spin
                d=log10(sqrt((x0-x)**2+(y0-y)**2+(z0-z)**2))
                bin = int(d/abin) + 1
                mass_gs(bin) = mass_gs(bin) + mgs
                mass_dm(bin) = mass_dm(bin) + mdm
                mass_st(bin) = mass_st(bin) + mst
                part_gs(bin) = part_gs(bin) + npgs
                part_dm(bin) = part_dm(bin) + npdm
                part_st(bin) = part_st(bin) + npst
                sfrs(bin)    = sfrs(bin)    + ssfr
!                e_gs(bin)    = e_gs(bin)    + e

!                h_gs(bin)    = h_gs(bin)    + h
!                hsml_gs(bin) = hsml_gs(bin) + hsml
!                l_hl(bin)    = l_hl(bin)    + sqrt(lx**2+ly**2+lz**2)
!                spin_hl(bin) = spin_hl(bin) + spin
        enddo
        close(10)

        open(10,file='perfiles.dat',status='unknown')
        do i=1,bines
                write(10,*) 10**((i-.5)*abin),mass_gs(i),mass_dm(i),mass_st(i),part_gs(i),part_dm(i),part_st(i),sfrs(i) !,e_gs(i),h_gs(i),hsml_gs(i),l_hl(i),spin_hl(i)
        enddo
        close(10)
        
endprogram
