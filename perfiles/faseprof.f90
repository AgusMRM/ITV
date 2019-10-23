program
        use modulos
        use openmp
        implicit none
        integer,parameter:: bines=20
        real,parameter :: rmax=25, rmin=0
        integer :: j
        real, allocatable(:) :: d
        integer,allocatable(:) :: bin,seccion
        real :: xh,yhe,mu,mp,kcgs,vv

        abin=log10(rmax-rmin)/real(bines)
        call reader()
        allocate(d(nall(0)),bin(nall(0)),seccion(nall(0)))

        xbox=403.8960 
        ybox=459.8882
        zbox=440.9021 
        xc=408.205481 - xbox + 250 
        yc=457.777839 - ybox + 250   
        zc=441.538681 - zbox + 250
        
        xh=0.76
        yhe=(1.0-xh)/(4.0*xh)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1e10

        !$OMP PARALLEL DEFAULT(NONE)
        !$OMP PRIVATE (i   ) &
        !$OMP SHARED d,nall,pos,xc,yc,zc,xh,yhe,ne,mp,kcgs,vv,bin,abin,dens,seccion)
        !$OMP DO SCHEDULE (DYNAMIC)
        do i=1,nall(0)
                d(i) = sqrt((pos(1,i)-xc)**2+(pos(2,i)-yc)**2+(pos(3,i)-zc)**2)
               
                if (d(i)<rmax) then
                mu = (1.0-yhe)/(1+yhe+ne(i))
                temp(i)=(5./3.-1.)*u(i)*vv*mu*mp/kcgs 
                bin(i) = int(log10(d(i)-rmin)/abin) + 1 
                        
                        if (temp(i) < 10**5 & dens(i) < 10**2) then
                               seccion(i) = 1 
                        elseif (temp(i) < 10**5 & dens(i) > 10**2) then
                               seccion(i) = 2 
                        elseif (temp(i) > 10**5 & dens(i) < 10**2) then
                               seccion(i) = 3 
                        else 
                               seccion(i) = 4 
                        endif 
                      
                endif 
        enddo
        !$OMP END DO
        !$OMP BARRIER
        !$OMP END PARALLEL

        do i=1,nall(0)

        enddo 




endprogram
