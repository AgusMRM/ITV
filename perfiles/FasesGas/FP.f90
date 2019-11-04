program PERFILGAS
        use modulos
        implicit none
        integer,parameter:: bines=20, VOID=1198
        real,parameter :: rmax=25, rmin=0, pi=acos(-1.)
        real, allocatable:: d(:),te(:)
        integer :: bin
        integer,dimension(bines) :: p_difu, p_whim, p_cond, p_hot
        real :: xh,yhe,mu,mp,kcgs,vv,abin,rad,distancia
        real :: xbox,ybox,zbox,xc,yc,zc,vol
       
        abin=log10(rmax-rmin)/real(bines)
        call reader()
        allocate(d(nall(0)),te(nall(0)))
        
        d=0; te=0
        p_difu=0; p_whim=0; p_cond=0; p_hot=0

        if (VOID==1198) then        
        xbox=411.2170 
        ybox=162.1655
        zbox=453.0553
        elseif (VOID==1373) then
        xbox=403.8960 
        ybox=459.8882
        zbox=440.9021 
        xc=408.205481 - xbox + 250 
        yc=457.777839 - ybox + 250   
        zc=441.538681 - zbox + 250
        endif
       stop 
        xh=0.76
        yhe=(1.0-xh)/(4.0*xh)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1e10

        do i=1,nall(0)
                d(i) = sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)
                mu = (1.0-yhe)/(1+yhe+ne(i))
                te(i)=(5./3.-1.)*u(i)*vv*mu*mp/kcgs
                        if (te(i) <= 10**(4.5) .and. rho(i) < 10**2 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int(log10(d(i)-rmin)/abin) + 1 
                                p_difu(bin) = p_difu(bin) + 1
                        elseif (te(i) <= 10**(4.5) .and. rho(i) >= 10**2 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int(log10(d(i)-rmin)/abin) + 1 
                                p_cond(bin) = p_cond(bin) + 1
                        elseif (te(i) > 10**(4.5) .and. rho(i) < 10**2 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int(log10(d(i)-rmin)/abin) + 1 
                                p_whim(bin) = p_whim(bin) + 1
                        elseif (te(i) > 10**(4.5) .and. rho(i) >= 10**2 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int(log10((d(i)-rmin))/abin) + 1 
                                p_hot(bin) = p_hot(bin) + 1
                        endif 
        enddo
        open(10,file='perfiles_FasesGAS.dat',status='unknown')
        do i =1,bines
                vol = (4./3.)*pi*(10**(i*abin)**3)
                write(10,*) 10**(i*abin), p_difu(i)/vol, p_cond(i)/vol, p_whim(i)/vol, p_hot(i)/vol 
        enddo
        close(10)


endprogram