! el radio virial mas grande es de ~1.5 Mpc, entonces yo voy a buscar
! particulas que esten a 2 radios virial es del centro del halos, esto
! es 3 Mpc. El entero mas cercano que es multiplo de 500 es 4, entonces
! hice 500/4 = 125 celdas (cada una de 4 Mpc de lado. 
!---------------------------------------------------------------------
! hay que compilar el modulo y subroutina de el archivo modulos.f90
! ej:  gfortran -o pp modulos.f90 program.f90 

program first
        use modulos
        implicit none
        integer,parameter::halos=71417-19, box=500, cell=125
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
        real,dimension(3,halos) :: pos_hl, vel_hl
        real,dimension(halos) :: rvir,spin,lx_hl,ly_hl,lz_hl
        integer,dimension(halos) :: id_hl,np_hl
        real :: nn, rmax, rmx,abin,masa_dm,masa_gs,masa_st,masa_t,masa_td
        integer,dimension(cell,cell,cell) :: tot_gs, head_gs
        integer,dimension(cell,cell,cell) :: tot_dm, head_dm
        integer,dimension(cell,cell,cell) :: tot_td, head_td
        integer,dimension(cell,cell,cell) :: tot_st, head_st
        integer,allocatable :: link_gs(:)
        integer,allocatable :: link_dm(:)
        integer,allocatable :: link_td(:)
        integer,allocatable :: link_st(:)
        integer:: bxg,byg,bzg,cand,bx,by,bz,bxdm,bydm,bzdm,partcls
        real :: rx,ry,rz,x,y,z,velx,vely,velz,G,L_dm,L_gs,L_st,L_td,espin_dm,ssfr,e,h,hsml_gs,mas
  !-----------------------------------------------------------------------------     
        call reader()

        STOP

        !mas=0
        !do i=1,nall(0)
        !        mas=mas+mass(i)
        !enddo
        !print*, mass(2), mas/real(nall(0))

        !stop
  !-----------------------------------------------------------------------------
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
        write(*,*) 'PROGRAMA PRINCIPAL'
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
        rmax = 0

        open(7,file='halos_S1373.dat',status='unknown')  !saltear 19 filas
        do i=1,19
                read(7,*)
        enddo
        print*, 'empiezo a leer halos'
        do i = 1, halos
                read(7,*) id_hl(i), np_hl(i), nn, nn, rvir(i), nn, nn, nn,pos_hl(1,i), pos_hl(2,i), pos_hl(3,i), &
                                vel_hl(1,i), vel_hl(2,i), vel_hl(3,i),nn,nn,nn,nn,spin(i)
                if (rvir(i) >= rmax) rmax = rvir(i)
        enddo 
        close(7)
write(*,*) 'MAXIMO RADIO VIRIAL:',rmax*1e-3,'Mpc'
        abin    = real(box)/real(cell)
      ! stop 
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

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
write(*,*) 'COMENZANDO CALCULOS'
!open(10,file='propiedades_halos.dat',status='unknown')
open(10,file='halosprop_2rvir.dat',status='unknown')
        do i=1,halos                            
                x = pos_hl(1,i)
                y = pos_hl(2,i)
                z = pos_hl(3,i)
                velx = vel_hl(1,i)
                vely = vel_hl(2,i)
                velz = vel_hl(3,i)

                bx = int(x/abin) + 1
                by = int(y/abin) + 1
                bz = int(z/abin) + 1
             !------------------------------------   
                part_td = 0 
                call chequeo(x,y,z,bx,by,bz,part_td,ntid,cell,tot_td, &
                               head_td,pos_td,vel_td,link_td,rvir(i))
                if (part_td>0) cycle
             !-----------------------------------
                part_gs = 0 
                call cuentasgas(x,y,z,velx,vely,velz,bx,by,bz,part_gs,ngas,cell,tot_gs, &
                                head_gs,pos_gs,vel_gs,link_gs,rvir(i),mass_gs,sfr,ne,nh,hsml,masa_gs,L_gs,ssfr,e,h,hsml_gs) 
                part_dm = 0 
                call cuentas(x,y,z,velx,vely,velz,bx,by,bz,part_dm,ndm,cell,tot_dm, &
                               head_dm,pos_dm,vel_dm,link_dm,rvir(i),mass_dm,masa_dm,L_dm) 
                part_st = 0 
                call cuentas(x,y,z,velx,vely,velz,bx,by,bz,part_st,nst,cell,tot_st, &
                               head_st,pos_st,vel_st,link_st,rvir(i),mass_st,masa_st,L_st)
                if (part_dm>25) then
                        !espin_dm = L_dm/(masa_dm*rvir(i)**2)*sqrt(G*masa_dm/(rvir(i)**3))
                        write(10,*) x,y,z,rvir(i),np_hl(i),part_gs,part_dm,part_st,part_td,masa_gs, &
                                        masa_dm,masa_st,L_gs,L_dm,ssfr,e,h,hsml_gs,lx_hl(i),ly_hl(i),lz_hl(i),spin(i)
                endif
              !  endif
        enddo
       close(10) 
endprogram first

subroutine cuentas(x,y,z,velx,vely,velz,bx,by,bz,o,n,cell,tot,head,pos,vel,link,rvir,mass,masa,Ltot)
        !use modulos
        implicit none
        integer:: cell,n
        integer,dimension(cell,cell,cell) :: tot, head
        real,dimension(3,n):: pos,vel 
        real,dimension(n) :: mass
        integer,dimension(n) :: link       
        integer ::j, o,k,m,l,bx,by,bz,p
        real ::x,y,z,rx,ry,rz,rvir,masa,Lx,Ly,Lz,Ltot,vx,vy,vz,velx,vely,velz
                masa = 0
                o    = 0
                Lx   = 0
                Ly   = 0
                Lz   = 0
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
                                        Lx = Lx + ry*vz - rz*vy
                                        Ly = Ly - rx*vz + rz*vx
                                        Lz = Lz + rx*vy - ry*vx                 
                                        masa = masa + mass(p)    
                                        o = o + 1
                                endif
                               p = link(p)       
                        enddo
                  enddo      
                 enddo
                enddo
                Ltot=sqrt(Lx**2+Ly**2+Lz**2)
endsubroutine
subroutine cuentasgas(x,y,z,velx,vely,velz,bx,by,bz,o,n,cell,tot,head,pos,vel,link,rvir,mass,sfr,ne,nh,hsml,masa, &
                                      Ltot,ssfr , e,h,hsml_gs)
        !use modulos
        implicit none
        integer:: cell,n
        integer,dimension(cell,cell,cell) :: tot, head
        real,dimension(3,n):: pos,vel 
        real,dimension(n) :: mass,sfr,ne,nh,hsml
        integer,dimension(n) :: link       
        integer ::j, o,k,m,l,bx,by,bz,p
        real::x,y,z,rx,ry,rz,rvir,masa,Lx,Ly,Lz,Ltot,vx,vy,vz,velx,vely,velz,ssfr,e,h,hsml_gs
                masa = 0
                o    = 0
                Lx   = 0
                Ly   = 0
                Lz   = 0
                ssfr = 0
                e    = 0
                h    = 0
                hsml_gs = 0
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
                                        Lx = Lx + ry*vz - rz*vy
                                        Ly = Ly - rx*vz + rz*vx
                                        Lz = Lz + rx*vy - ry*vx                 
                                        masa = masa + mass(p) 
                                        ssfr = ssfr + sfr(p)  
                                        e    = e + ne(p)
                                        h    = h + nh(p)
                                        hsml_gs=hsml_gs + hsml(p) 
                                        o = o + 1
                                endif
                               p = link(p)       
                        enddo
                  enddo      
                 enddo
                enddo
                Ltot=sqrt(Lx**2+Ly**2+Lz**2)
                e = e/real(o)
                h = h/real(masa)
                hsml_gs=hsml_gs/real(o)         ! TENGO EL HSML PROMEDIO DEL HALO               
endsubroutine
subroutine chequeo(x,y,z,bx,by,bz,o,n,cell,tot,head,pos,vel,link,rvir)
        !use modulos
        implicit none
        integer:: cell,n
        integer,dimension(cell,cell,cell) :: tot, head
        real,dimension(3,n):: pos,vel 
        real,dimension(n) :: mass
        integer,dimension(n) :: link       
        integer ::j, o,k,m,l,bx,by,bz,p
        real ::x,y,z,rx,ry,rz,rvir,masa,Lx,Ly,Lz,Ltot,vx,vy,vz,velx,vely,velz
                masa = 0
                o    = 0
                Lx   = 0
                Ly   = 0
                Lz   = 0
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
                                if (sqrt(rx**2+ry**2+rz**2)<3*rvir*1e-3) then
                                        o = o + 1
                                endif
                               p = link(p)       
                        enddo
                  enddo      
                 enddo
                enddo
endsubroutine
