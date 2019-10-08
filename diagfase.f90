rogram
        use modulos
        implicit none
        integer :: i
        real :: xh,yhe,mp,kcgs,vv,mu,te
        open(10,file='diagfases.dat',status='unknown')
        
        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1E10
        
        do i=1,nall(0)
                mu=(1.0-yHe)/(1+yHe+ne(i))
                te=(5./3.-1.)*u*vv*mu*mp/kcgs
                write(10,*) te, dens(i)
        enddo


        close(10) 
endprogram 
