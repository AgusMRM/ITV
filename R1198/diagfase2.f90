program df
        use modulos
        implicit none
        integer,parameter :: VOID=1373
        real :: xh,yhe,mp,kcgs,vv,mu,te,d
        real,allocatable :: dist(:), dist2(:),dist3(:)
        integer,allocatable :: o2(:) 
        open(10,file='diagfases_ord.dat',status='unknown',action='write')
        call reader() 
        allocate(dist(nall(0)),dist2(nall(0)),dist3(nall(0)))
        do i=1,nall(0)
                dist(i)=sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2)
        enddo
        dist2=dist
        dist3=dist
        
        call sort2(nall(0),dist,dens)
        call sort2(nall(0),dist2,u) 
        call sort2(nall(0),dist3,ne) 
        write(*,*) 'VECTORES YA ORDENADOS'
        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1E10
        do i=1,nall(0)
                
                if (dist2(i) > 15) cycle
                mu=(1.0-yHe)/(1+yHe+ne(i))
                te=(5./3.-1.)*u(i)*vv*mu*mp/kcgs

                write(10,*) dens(i), te, dist(i)
                
        enddo


        close(10) 
        close(11) 
        close(12) 
endprogram 

