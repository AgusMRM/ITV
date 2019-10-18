program hp
        use modulos
        implicit none
        integer,parameter :: halo=4332
        integer :: ngas,ndm,ntid,nst
        integer :: k,m,l,p,particles,part_gs,part_dm,part_st,part_td
        real :: nn, rmax,rmx,abin,masa_dm,masa_gs,masa_st,masa_t,masa_td,rho,rvir
        integer:: bxg,byg,bzg,cand,bx,by,bz,bxdm,bydm,bzdm,partcls
        real :: rx,ry,rz,x,y,z,velx,vely,velz,G,L_dm,L_gs,L_st,L_td,espin_dm,ssfr,e,h

        call reader()

!---------------------------------------------------------------------------------------------------        
        open(7,file='halos_R1198.dat',status='unknown')  !saltear 19 filas
        do i=1,19 + halo -1 
                read(7,*)
        enddo
        read(7,*) a,a,a,a,rvir,a,a,a,x,y,z
        close(7)
        print*, rvir, x,y,z
        open(12,file='halossinstar.dat',status='unknown')

        do i=1,nall(0)
                if ((sqrt(pos(1,i)-x)**2 + (pos(2,i)-y)**2 + (pos(3,i)-z)**2)<2*rvir*1e-3) then
                write(12,*) pos(1,i), pos(2,i), pos(3,i), &
                                                    vel(1,i), vel(2,i), vel(3,i), &   
                                                              u(i), dens(i), ne(i), nh(i),sfr(i)                       
                endif
        enddo 
        close(12)
endprogram 
