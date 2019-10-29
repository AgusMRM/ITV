program recorrido
        use modulos
        use OMP_lib  
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
   ! LA SUBRUTINA CONTADOR ME DICE CUANTAS PARTICULAS HAY EN EL RANGO QUE ME
   ! INTERESA, DE ESTA MANERA OBTENGO LAS DIMENSION DE LOS VECTORES CON LOS QUE
   ! VOY A TRABAJAR--------     
        write(*,*) 'PARTICULAS IDENTIFICADAS',dmn
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
   ! EL VECTOR IDP TIENE LOS ID DE LAS PARTICULAS QUE ME INTERESA RASTREAR,
   ! AHORA LO QUE TENGO QUE HACER ES ORDENARLO EN FORMA CRECIENTE   
        call orden(dmn,idp) 
        write(*,*) 'VECTOR ID ORDENADO: READY'
        snapshot=49
        deallocate(pos,vel,id,idch,idgn,u,dens,mass,ne)
        
   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP PRIVATE (i,snapshot,pos,vel,id,idch,idgn,u,dens,mass,ne) &
   !$OMP SHARED (dmn,idp,hist)   
   !$OMP DO SCHEDULE (DYNAMIC)
        
        do i=40,50
                snapshot=i
                call Id_Finder(dmn,idp,snapshot,hist)         
        enddo

   !$OMP END DO
   !$OMP BARRIER
   !$OMP END PARALLEL
        open(20,file='evolucion.dat',status='unknown')
        do i=1,dmn
                write(20,*) hist(:,i)
        enddo
        close(20)
        deallocate(hist,idp)

endprogram
