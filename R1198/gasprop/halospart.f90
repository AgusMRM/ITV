program hp
        use modulos
        implicit none
        integer,parameter :: halo=1480, box=500, cell=50
        integer :: ngas,ndm,ntid,nst
        integer :: k,m,l,p,particles,part_gs,part_dm,part_st,part_td
        real,allocatable :: pos_gs(:,:), vel_gs(:,:)
        real,allocatable :: u_gs(:), mass_gs(:)
        integer,allocatable :: id_gs(:)
        real,allocatable :: pos_dm(:,:), vel_dm(:,:)
        real,allocatable :: mass_dm(:)
        integer,allocatable :: id_dm(:)
        real,allocatable :: pos_td(:,:), vel_td(:,:)
        real,allocatable :: mass_td(:)
        integer,allocatable :: id_td(:)
        real,allocatable :: pos_st(:,:), vel_st(:,:)
        real,allocatable :: mass_st(:)
        integer,allocatable :: id_st(:)
        real :: nn, rmax,rmx,abin,masa_dm,masa_gs,masa_st,masa_t,masa_td,rho,rvir
        integer,dimension(cell,cell,cell) :: tot_gs, head_gs
        integer,dimension(cell,cell,cell) :: tot_dm, head_dm
        integer,dimension(cell,cell,cell) :: tot_td, head_td
        integer,dimension(cell,cell,cell) :: tot_st, head_st
        integer,allocatable :: link_gs(:)
        integer,allocatable :: link_dm(:)
        integer,allocatable :: link_td(:)
        integer,allocatable :: link_st(:)
        integer:: bxg,byg,bzg,cand,bx,by,bz,bxdm,bydm,bzdm,partcls
        real :: rx,ry,rz,x,y,z,velx,vely,velz,G,L_dm,L_gs,L_st,L_td,espin_dm,ssfr,e,h

        call reader()

        ngas=nall(0)
        ndm=nall(1)
        ntid=nall(2)
        nst=nall(4)
        allocate(pos_gs(3,ngas),vel_gs(3,ngas))
        allocate(u_gs(ngas),mass_gs(ngas),id_gs(ngas),link_gs(ngas)) 
        allocate(pos_dm(3,ndm),vel_dm(3,ndm))
        allocate(mass_dm(ndm),id_dm(ndm),link_dm(ndm)) 
        allocate(pos_td(3,ntid),vel_td(3,ntid))
        allocate(mass_td(ntid),id_td(ntid),link_td(ntid)) 
        allocate(pos_st(3,nst),vel_st(3,nst))
        allocate(mass_st(nst),id_st(nst),link_st(nst)) 

        do i=1,nall(0)
                pos_gs(1,i) = pos(1,i)
                pos_gs(2,i) = pos(2,i)
                pos_gs(3,i) = pos(3,i)
                vel_gs(1,i) = vel(1,i)
                vel_gs(2,i) = vel(2,i)
                vel_gs(3,i) = vel(3,i)
                mass_gs(i)  = mass(i)
                u_gs(i)     = u(i)

        enddo
        j=0
        do i=1+nall(0),nall(0)+nall(1)
                j=j+1
                pos_dm(1,j) = pos(1,i)
                pos_dm(2,j) = pos(2,i)
                pos_dm(3,j) = pos(3,i)
                vel_dm(1,j) = vel(1,i)
                vel_dm(2,j) = vel(2,i)
                vel_dm(3,j) = vel(3,i)
                mass_dm(j)  = mass(i)
        enddo
        j=0
        do i=1+nall(0)+nall(1),nall(0)+nall(1)+nall(2)
                j=j+1
                pos_td(1,j) = pos(1,i)
                pos_td(2,j) = pos(2,i)
                pos_td(3,j) = pos(3,i)
                vel_td(1,j) = vel(1,i)
                vel_td(2,j) = vel(2,i)
                vel_td(3,j) = vel(3,i)
                mass_td(j)  = mass(i)
        enddo
        j=0
        do i=1+nall(0)+nall(1)+nall(2)+nall(3),nall(0)+nall(1)+nall(2)+nall(3)+nall(4)
                j=j+1
                pos_st(1,j) = pos(1,i)
                pos_st(2,j) = pos(2,i)
                pos_st(3,j) = pos(3,i)
                vel_st(1,j) = vel(1,i)
                vel_st(2,j) = vel(2,i)
                vel_st(3,j) = vel(3,i)
                mass_st(j)  = mass(i)
        enddo
!---------------------------------------------------------------------------------------------------        
        open(7,file='halos_R1198.dat',status='unknown')  !saltear 19 filas
        do i=1,19 + halo -1 
                read(7,*)
        enddo
        read(7,*) a,a,a,a,rvir,a,a,a,x,y,z
        close(7)
abin    = real(box)/real(cell)
tot_gs  = 0
head_gs = 0
link_gs = 0
tot_dm  = 0       
head_dm = 0
link_dm = 0
tot_td  = 0
head_td = 0
link_td = 0
tot_st  = 0
head_st = 0
link_st = 0
write(*,*) 'LiNKEANDO DM'
call linkedlist(ndm,abin,cell,pos_dm,head_dm,tot_dm,link_dm)
write(*,*) 'LiNKEANDO GAS'
call linkedlist(ngas,abin,cell,pos_gs,head_gs,tot_gs,link_gs)
write(*,*) 'LiNKEANDO TIDALES'
call linkedlist(ntid,abin,cell,pos_td,head_td,tot_td,link_td)
write(*,*) 'LiNKEANDO ESTRELLAs'
call linkedlist(nst,abin,cell,pos_st,head_st,tot_st,link_st)

        bx = int(x/abin) + 1
        by = int(y/abin) + 1
        bz = int(z/abin) + 1

        open(12,file='halosconstar.dat',status='unknown')
        call cuentasgas(x,y,z,bx,by,bz,ngas,cell,tot_gs,head_gs,pos_gs,vel_gs,&
                link_gs,rvir,dens,sfr,u,ne,nh)
        close(12)
endprogram 
subroutine cuentasgas(x,y,z,bx,by,bz,n,cell,tot,head,pos,vel,link,rvir,dens,sfr,u,ne,nh)
                                    
        !use modulos
        implicit none
        integer:: cell,n
        integer,dimension(cell,cell,cell) :: tot, head
        real,dimension(3,n):: pos,vel 
        real,dimension(n) :: mass,sfr,ne,nh,hsml,dens,u
        integer,dimension(n) :: link       
        integer ::j, o,k,m,l,bx,by,bz,p
        real::x,y,z,rx,ry,rz,rvir,masa,Lx,Ly,Lz,Ltot,vx,vy,vz,velx,vely,velz
        real::ssfr,rho,maxhsml
                do k = bx-1,bx+1
                 do m = by-1,by+1
                  do l = bz-1,bz+1
                        if (tot(k,m,l) == 0) cycle
                        p = head(k,m,l)
                        do j = 1,tot(k,m,l)
                                rx = pos(1,p) - x 
                                ry = pos(2,p) - y 
                                rz = pos(3,p) - z
                                vx = vel(1,p) - velx
                                vy = vel(2,p) - vely
                                vz = vel(3,p) - velz
                                if (sqrt(rx**2+ry**2+rz**2)<2*rvir*1e-3) then
                                        write(12,*) pos(1,p), pos(2,p), pos(3,p), &
                                                    vel(1,p), vel(2,p), vel(3,p), &   
                                                    u(p), dens(p), ne(p), nh(p),sfr(p)                       
                                endif
                               p = link(p)       
                        enddo
                  enddo      
                 enddo
                enddo
                rho=rho/real(o)
endsubroutine
