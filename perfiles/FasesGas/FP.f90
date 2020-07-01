program PERFILGAS
        use modulos
        implicit none
        integer,parameter:: bines=20, VOID=832
        real,parameter :: rmax=35, rmin=2, pi=acos(-1.)
        real, allocatable:: d(:),te(:)
        integer :: bin,tot
        integer*8,dimension(bines) :: p_difu, p_whim, p_cond, p_hot,p_st
        real :: xh,yhe,mu,mp,kcgs,vv,abin,rad,distancia
        real :: xbox,ybox,zbox,xc,yc,zc,vol,r0
        real :: x,y,z ,rv,dist      
        abin=(log10(rmax)-log10(rmin))/real(bines)
        print*, abin
        call reader()
        allocate(d(nall(0)),te(nall(0)))
        
        d=0; te=0
        p_difu=0; p_whim=0; p_cond=0; p_hot=0; p_st=0

        open(13,file='/home/arodriguez/Doctorado/itv/voids_boxes.dat',status='unknown')

        !   open(13,file='/mnt/is2/fstasys/ITV/base09/voids/voids_new.dat')
           do i=1,VOID-1
           read(13,*)
           enddo
   
           read(13,*) rv,x,y,z,xbox,ybox,zbox
           xc=x-xbox+250
           yc=y-ybox+250
           zc=z-zbox+250
           close(13)
           write(*,*) x,y,z

        xh=0.76
        yhe=(1.0-xh)/(4.0*xh)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1e10

        do i=1,nall(0)
                d(i) = sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)
                mu = (1.0-yhe)/(1+yhe+ne(i))
                te(i)=(5./3.-1.)*u(i)*vv*mu*mp/kcgs
                        if (te(i) <= 10**(5) .and. rho(i) < 116.24 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int((log10(d(i))-log10(rmin))/abin) + 1 
                                p_difu(bin) = p_difu(bin) + 1
                        elseif (te(i) <= 10**(5) .and. rho(i) >= 116.24 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int((log10(d(i))-log10(rmin))/abin) + 1 
                                p_cond(bin) = p_cond(bin) + 1
                        elseif (te(i) > 10**(5) .and. rho(i) < 116.24 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int((log10(d(i))-log10(rmin))/abin) + 1 
                                p_whim(bin) = p_whim(bin) + 1
                        elseif (te(i) > 10**(5) .and. rho(i) >= 116.24 .and. d(i) > rmin .and. d(i)<= rmax) then
                                bin = int((log10(d(i))-log10(rmin))/abin) + 1 
                                p_hot(bin) = p_hot(bin) + 1
                        endif 
        enddo
        do i=1+nall(0)+nall(1)+nall(2),nall(0)+nall(1)+nall(2)+nall(4)
                dist = sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)
                if (dist > rmin .and. dist<= rmax) then
                bin = int((log10(dist)-log10(rmin))/abin) + 1 
                p_st(bin) = p_st(bin) + 1
                endif
        enddo
        open(10,file='perfiles_FasesGAS.dat',status='unknown')
        r0=rmin
        do i =1,bines
                tot = p_difu(i) + p_cond(i) + p_whim(i) + p_hot(i)        
                rad=10**((i)*abin+log10(rmin))
              !  vol = (4./3.)*pi*(10**(i*abin)**3)
                vol = (4./3.)*pi*(rad**3-r0**3)
                write(10,*) rad, p_difu(i)/vol, p_cond(i)/vol, p_whim(i)/vol, p_hot(i)/vol, &
                        sqrt(real(p_difu(i)))/vol, sqrt(real(p_cond(i)))/vol, sqrt(real(p_whim(i)))/vol, &
                        sqrt(real(p_hot(i)))/vol, p_difu(i)/real(tot), &
                        p_cond(i)/real(tot), p_whim(i)/real(tot), p_hot(i)/real(tot)
                r0=rad 
        enddo
        close(10)
        open(11,file='perfiles_FasesBARIONES.dat',status='unknown')
        r0=rmin
        do i =1,bines
                tot = p_difu(i) + p_cond(i) + p_whim(i) + p_hot(i) + p_st(i)     
                rad=10**((i)*abin+log10(rmin))
                vol = (4./3.)*pi*(rad**3-r0**3)
                write(11,*) rad, p_difu(i)/vol, p_cond(i)/vol, p_whim(i)/vol, p_hot(i)/vol, p_st(i)/vol, &
                        sqrt(real(p_difu(i)))/vol, sqrt(real(p_cond(i)))/vol, sqrt(real(p_whim(i)))/vol, &
                        sqrt(real(p_hot(i)))/vol,sqrt(real(p_st(i)))/vol, p_difu(i)/real(tot), &
                        p_cond(i)/real(tot), p_whim(i)/real(tot), p_hot(i)/real(tot), p_st(i)/real(tot)
                r0=rad 
        enddo
        close(11)


endprogram
