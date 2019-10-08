program df
        use modulos
        implicit none
        integer,parameter :: VOID=1373
        real :: xh,yhe,mp,kcgs,vv,mu,te,d

        open(10,file='df_raster.dat',status='unknown',action='write')
        call reader() 
       

        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1E10
        do i=1,nall(0),200000
                
                mu=(1.0-yHe)/(1+yHe+ne(i))
                te=(5./3.-1.)*u(i)*vv*mu*mp/kcgs
                if (te > 1e6 .and. dens(i) < 1e2) write(10,*) id(i) 
        enddo


        close(10) 
endprogram 
