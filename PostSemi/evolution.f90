program recorrido
        use modulos 
        implicit none
        real,parameter:: rmax=9.57, rmin=8.08 
        character(len=200)::snumber
        real:: xbox,ybox,zbox,xc,yc,zc,xh,yhe,mp,mu,kcgs,vv,te,d
        integer::p,dmn,snapshot
        integer,allocatable:: idp(:)
        real,allocatable:: hist(:,:) 
        write(snumber,'(I3)') 50
        call reader(snumber)
        print*, redshift
        xbox=403.8960 
        ybox=459.8882
        zbox=440.9021
        xc=408.205481-xbox+250
        yc=457.777839-ybox+250
        zc=441.538681-zbox+250 
        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1e10
        
        p=0
        call contador(p,rmin,rmax,xc,yc,zc,mu,yhe,te,vv,mp,kcgs)
        dmn=p
        allocate(idp(dmn))
        allocate(hist(10,dmn))
        idp= 0 
        p=0
        do i=1,nall(0)
                d=sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)       
                if (d<rmax .and. d>rmin) then 
                 mu=(1.0-yHe)/(1+yHe+ne(i))
                 te=(5./3.-1.)*u(i)*vv*mu*mp/kcgs
                 if (te>10**(5.5) .and. dens(i)<10**(1.)) then
                        p = p +1
                        idp(p) = id(i)
                 endif         
                endif 
        enddo

        do i=45,50
                call Id_Finder(dmn,idp,snapshot,hist)         
        enddo 
        deallocate(hist,idp)
endprogram
