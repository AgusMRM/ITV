program df
        use modulos
        implicit none
        integer,parameter :: VOID=1373
        real :: xh,yhe,mp,kcgs,vv,mu,te,d

        !open(10,file='diagfases_int.dat',status='unknown',action='write')
        !open(11,file='diagfases_wll.dat',status='unknown',action='write')
        !open(12,file='diagfases_ext.dat',status='unknown',action='write')
        !open(13,file='diagfases_sf.dat',status='unknown',action='write')
        open(13,file='diagfases_VOID_fueraRE.dat',status='unknown',action='write')
        call reader() 
       

        xH=0.76
        yHe=(1.0-xH)/(4.0*xH)
        mp=1.6726E-24
        kcgs=1.3807E-16
        vv=1E10
        do i=1,nall(0)

                d=sqrt((pos(1,i)-250)**2+(pos(2,i)-250)**2+(pos(3,i)-250)**2)
                mu=(1.0-yHe)/(1+yHe+ne(i))
                te=(5./3.-1.)*u(i)*vv*mu*mp/kcgs
                !if(d<=6)                write(10,*) te, dens(i), sfr(i)
                !if(d>6 .and. d<=12)     write(11,*) te, dens(i), sfr(i)
                !if(d>12 .and. d<=18)    write(12,*) te, dens(i), sfr(i)
                if (d>20 .and. d<25)        write(13,*) te, dens(i), sfr(i)
        enddo


        close(10) 
        close(11) 
        close(12) 
endprogram 
