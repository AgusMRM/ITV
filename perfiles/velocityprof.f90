program first0
        implicit none
        integer,parameter :: bines=30, VOID=443,halos=172548 -19  !106409 
        real,parameter :: rmax=35,rmin=3, pi=acos(-1.),corte=100!1232
        real :: abin,d,vol,rm, rand,minmass,maxmass,minsfr,maxsfr,abin2,rad,r0
        integer :: bin,bin2  
        integer :: dmp,gsp,i 
        real :: x,y,z,rv,ri,vol2,masa,a,xc,yc,zc,xbox,ybox,zbox
        real :: rx,ry,rz,tita,fi,vr,velx,vely,velz,vx,vy,vz
        real,dimension(bines):: sf,hsml0
        real,dimension(bines):: wgth,mass_dm
        integer,dimension(bines):: num_gs,num_dm,num_st
        real,dimension(bines):: vel_dm
        real,dimension(3,halos) :: pos,vel
        integer,dimension(halos) :: pnum
  !-----------------------------------------------------------------------------    
 !open(24,file='/mnt/is2/fstasys/ITV/base09/rockstar/halos_0.5.ascii',status='unknown',action='read')
 open(24,file='/home/arodriguez/Doctorado/snap50/halos_ 443.ascii',status='unknown',action='read')
 do i=1,19
        read(24,*)
 enddo
 do i=1,halos
        read(24,*) a,pnum(i),a,a,a,a,a,a, pos(1,i), pos(2,i), pos(3,i),&
                vel(1,i), vel(2,i), vel(3,i)
 enddo
 close(24) 

 open(13,file='/home/arodriguez/Doctorado/itv/voids_boxes.dat',status='unknown')
   do i=1,VOID-1
   read(13,*)
   enddo
   
   read(13,*) rv,x,y,z,xbox,ybox,zbox
   xc=x-xbox+250
   yc=y-ybox+250
   zc=z-zbox+250
   close(13)
   write(*,*) x,y,z
 open(14,file='/mnt/is2/fstasys/ITV/base09/voids/voids_new.dat')
   do i=1,VOID-1
   read(14,*)
   enddo
   read(14,*) a,a,a,a,velx,vely,velz
 close(14)
  !-----------------------------------------------------------------------------
        abin = (log10(rmax)-log10(rmin))/real(bines)        
        num_dm  = 0
        vel_dm  = 0
  !---------------------------------------------- DM ------------------------
        do i=1,halos
                d=sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)
                if (d < rmax .and. d>rmin .and. pnum(i)>= corte) then
               
                rx = (pos(1,i)-xc)
                ry = (pos(2,i)-yc)
                rz = (pos(3,i)-zc)
                vx = vel(1,i) - velx
                vy = vel(2,i) - vely
                vz = vel(3,i) - velz
                call tita_ang(rx,ry,rz,tita)
                call fi_ang(rx,ry,rz,fi)                           
                vr = vx*sin(tita)*cos(fi) + &
                    vy*sin(tita)*sin(fi) + vz*cos(tita) 
                bin = int((log10(d)-log10(rmin))/abin) + 1
                num_dm(bin)  = num_dm(bin) + 1
                vel_dm(bin)  = vel_dm(bin) + vr
                endif
        enddo
   !---------------------------------------------------------------------     
        open(17,file='velprof.dat',status='unknown')
        dmp = 0
        r0=rmin
        ri=0
        print*, vel_dm
        do i=1,bines
                dmp = dmp + num_dm(i)
                rad=10**((i)*abin+log10(rmin))
              !  vol2 = (4./3.)*pi*((10**(rad))**3-((10**(r0))**3))  ! ESTE ES
         !       PARA EL DIFERENCIAL
                vol  = (4./3.)*pi*(rad**3-rmin**3)
                vol2 = (4./3.)*pi*(rad**3-r0**3)
                !write(17,*) ((i-.5)*abin+rmin), dmp, vol,num_dm(i),vol2
                write(17,*) rad, dmp,&
                        vol,num_dm(i),vol2,vel_dm(i)/float(num_dm(i))
                r0 = rad
        enddo
        close(17)
endprogram
subroutine tita_ang(x,y,z,ang)
        implicit none
        real :: x,y,z,ang,pii
        pii = acos(-1.)
        if (z > 0 ) then
                ang = atan(sqrt(x**2+y**2)/z)
        elseif (z==0) then
                ang = pii/2.
        else        
                ang = pii + atan(sqrt(x**2+y**2)/z)
        endif
endsubroutine
subroutine fi_ang(x,y,z,ang)
        implicit none
        real :: x,y,z,ang,pii
        pii = acos(-1.)
        if (x>0 .and. y>0) then
                ang = atan(y/x)
        elseif (x>0 .and. y<0) then
                ang = atan(y/x)
        elseif (x==0) then
                ang = (pii/2.)*sign(1.,y)
        else
                ang = pii + atan(y/x)
        endif        
endsubroutine 
