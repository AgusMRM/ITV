program PERFILGAS
        use modulos
!        use OMP_lib
        implicit none
        integer,parameter:: bines=20, VOID=1373
        real,parameter :: rmax=20, rmin=0, pi=acos(-1.)
        real, allocatable:: d(:),te(:)
        integer,allocatable :: bin(:),seccion(:)
        integer,dimension(bines) :: p_difu, p_whim, p_cond, p_hot
        real :: xh,yhe,mu,mp,kcgs,vv,abin,rad,distancia
        real :: xbox,ybox,zbox,xc,yc,zc,vol,x,y,z,rv
        abin=log10(rmax-rmin)/real(bines)
        call reader()
        allocate(d(nall(0)),bin(nall(0)),seccion(nall(0)),te(nall(0)))
        
        d=0; bin=0; seccion=0; te=0
        p_difu=0; p_whim=0; p_cond=0; p_hot=0

!  xbox=411.217
!  ybox=162.1655
!  zbox=453.0553 
  
  xbox=403.896
  ybox=459.8882
  zbox=440.9021 
   
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
       ! xc=408.205481 - xbox + 250 
       ! yc=457.777839 - ybox + 250   
       ! zc=441.538681 - zbox + 250
        
        xh=0.76
        yhe=(1.0-xh)/(4.0*xh)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1e10

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE (i,mu)  &
        !$OMP SHARED (d,nall,pos,xc,yc,zc,xh,yhe,ne,mp,kcgs,vv,bin,te,u,abin,rho,seccion)
        !$OMP DO SCHEDULE (DYNAMIC)
        do i=1,nall(0)
                d(i) = sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)
                if ( d(i)<rmax ) then
                mu = (1.0-yhe)/(1+yhe+ne(i))
                te(i)=(5./3.-1.)*u(i)*vv*mu*mp/kcgs
                bin(i) = int(log10(d(i)-rmin)/abin) + 1 
                        if (te(i) <= 10**(5) .and. rho(i) < 116) then
                               seccion(i) = 1 
                                p_difu(bin(i)) = p_difu(bin(i)) + 1
                        elseif (te(i) <= 10**(5) .and. rho(i) >= 116) then
                               seccion(i) = 2 
                                p_cond(bin(i)) = p_cond(bin(i)) + 1
                        elseif (te(i) > 10**(5) .and. rho(i) < 116) then
                               seccion(i) = 3 
                                p_whim(bin(i)) = p_whim(bin(i)) + 1
                        elseif (te(i) > 10**(5) .and. rho(i) >=116) then
                               seccion(i) = 4 
                                p_hot(bin(i)) = p_hot(bin(i)) + 1
                        endif 
                      
                endif 
        enddo
        !$OMP END DO
        !$OMP BARRIER
        !$OMP END PARALLEL
       ! do i=1,nall(0)
       !         if (seccion(i) == 1 ) then
       !                 p_difu(bin(i)) = p_difu(bin(i)) + 1
       !         elseif (seccion(i) == 2 ) then
       !                 p_cond(bin(i)) = p_cond(bin(i)) + 1
       !         elseif (seccion(i) == 3 ) then
       !                 p_whim(bin(i)) = p_whim(bin(i)) + 1
       !         elseif (seccion(i) == 4) then
       !                 p_hot(bin(i)) = p_hot(bin(i)) + 1
       !         endif        
       ! enddo 

        open(10,file='perfiles_FasesGAS.dat',status='unknown')
        do i =1,bines
                vol = (4./3.)*pi*(10**(i*abin)**3)
                write(10,*) 10**(i*abin), p_difu(i)/vol, p_cond(i)/vol, p_whim(i)/vol, p_hot(i)/vol 
        enddo
        close(10)


endprogram
