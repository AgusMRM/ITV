program asd
        implicit none
        integer,parameter :: halos=52126,bines=18,VOID=1373
        real, parameter  :: lim1=6, lim2=12, lim3=18

        integer :: i,bin1,bin2,bin
        real :: abin1,abin2,abin3,d,xbox,ybox,zbox,l1,l2,l3
        real :: x,y,z,rvir,np,npgs,npdm,npst,nptd,mgs,mdm,mst,mtd,lgs,ldm
        real ::ssfr,e,h,hsml,lx,ly,lz,spin,maxdm,x0,y0,z0,rv
        real,dimension(bines) ::mass_dm1,mass_dm2,mass_dm3
        
        xbox=403.8960
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
        
        maxdm=0
        open(10,file='halosprop_2rvir.dat',status='unknown')
        open(11,file='halosprop_int.dat',status='unknown')
        open(12,file='halosprop_wll.dat',status='unknown')
        open(13,file='halosprop_ext.dat',status='unknown')
        open(14,file='halosconESTRELLAS.dat',status='unknown')
        open(15,file='halossinESTRELLAS.dat',status='unknown')
        open(16,file='halosGRANDES.dat',status='unknown')
        open(17,file='halosCHICOS.dat',status='unknown')
        do i=1,halos
                read(10,*) x0,y0,z0,rvir,np,npgs,npdm,npst,nptd,mgs,mdm,mst,lgs,ldm,&
                        ssfr,e,h,hsml,lx,ly,lz,spin
                d=(sqrt((x0-x)**2+(y0-y)**2+(z0-z)**2))
                if (d<=(lim1)) then
                write(11,*) d, mgs,mdm,mst, lgs, ldm, ssfr, spin
 
                elseif (d >(lim1) .and. d<=(lim2)) then
                write(12,*) d, mgs,mdm,mst, lgs, ldm, ssfr, spin
 
                elseif (d>(lim2) .and. d<(lim3)) then 
                write(13,*) d, mgs,mdm,mst, lgs, ldm, ssfr, spin
 
                endif
                if (mst>0) then
                        write(14,*) d,npgs, npdm, npst, mgs,mdm, mst, ssfr
                else
        
                        write(15,*) d,npgs, npdm, npst, mgs,mdm, mst, ssfr
                endif
                
                if (npdm>50) then
                        write(16,*) d,mgs, mdm, mst, ssfr
                else
                        write(17,*) d,mgs,mdm,mst,ssfr
                endif
                 
         enddo
        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)
        close(17)
        
endprogram
