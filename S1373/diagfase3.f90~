! EN ESTE PROGRAMA QUIERO VER DE DONDE VIENEN LAS PARTICULAS DE GAS QUE
! ESTAN EN LA ZONA CALIENTE PARA VER SI SON DEL FEEDBACK EN LA PARED
program df
        use modulos
        implicit none
        integer,parameter :: VOID=1373
        real :: xh,yhe,mp,kcgs,vv,mu,te,d,x,y,z,rv,xbox,ybox,zbox
        real,allocatable :: dist(:)
        integer,allocatable :: o2(:)
        open(10,file='diagfases_in_snp25.dat',status='unknown',action='write')
        open(11,file='diagfases_wall1_snp25.dat',status='unknown',action='write')
        open(12,file='diagfases_wall2_snp25.dat',status='unknown',action='write')
        open(14,file='diagfases_wall3_snp25.dat',status='unknown',action='write')
        open(15,file='diagfases_wall4_snp25.dat',status='unknown',action='write')
        open(16,file='diagfases_out_snp25.dat',status='unknown',action='write')
        call reader() 

        xbox=403.896 
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
        
        allocate(dist(nall(0)))
        
        do i=1,nall(0)
                dist(i)=sqrt((pos(1,i)-x)**2+(pos(2,i)-y)**2+(pos(3,i)-z)**2)
        enddo

        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1E10
        do i=1,nall(0)
                mu=(1.0-yHe)/(1+yHe+ne(i))
                te=(5./3.-1.)*u(i)*vv*mu*mp/kcgs
                if ( dist(i) < 6 ) then
                        write(10,*)  te, dens(i), dist(i), id(i)
        
                elseif (dist(i) > 6 .and. dist(i) <= 8)  then          
                        write(11,*)  te, dens(i), dist(i), id(i)
                elseif (dist(i) > 8 .and. dist(i) <= 9)  then          
                        write(12,*)  te, dens(i), dist(i), id(i)
                elseif (dist(i) > 9 .and. dist(i) <= 10)  then          
                        write(14,*)  te, dens(i), dist(i), id(i)
                elseif (dist(i) > 10 .and. dist(i) <= 12)  then          
                        write(15,*)  te, dens(i), dist(i), id(i)
                elseif (dist(i) > 12 .and. dist(i) <= 18 ) then
                        write(16,*)  te, dens(i), dist(i), id(i)
                endif
        enddo
        close(10) 
        close(11) 
        close(12) 
        close(14) 
        close(15) 
        close(16) 

        

endprogram 

