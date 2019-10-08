program tmep
        use modulos
        implicit none
        real:: xh,yhe,mp,kcgs,vv,mu,temperatura,rand
        
        call reader()
        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1E10
       open(10,file='rhoTemp.dat',status='unknown')
        do i=1,nall(0)
                call random_number(rand)
                if (rand<.005 ) then        
                mu=(1.0-yHe)/(1+yHe+ne(i))
                temperatura=(5./3.-1)*u(i)*vv*mu*mp/kcgs
                write(10,*) temperatura, dens(i)
                endif
        enddo
        
       close(10) 
endprogram
