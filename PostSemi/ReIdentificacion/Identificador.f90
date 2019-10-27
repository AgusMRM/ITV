program IdEnTiFiCaToR
        !use Identikit
        implicit none
        integer, parameter :: halos=71417-19, cell=16, bines=10000
        real,parameter :: DELTA=-0.9,rmax=12, rmin=3, halos_tot=850350
        real,parameter :: pi=acos(-1.)
        real ::a
        integer ::i,j
        integer,allocatable:: id(:),num_p(:)
        real,allocatable::    rvir(:)
        real,allocatable::    pos(:,:)
        real,dimension(bines) :: D_r
        integer,dimension(bines) :: cont
        real::x0,y0,z0,vol,vol0,abin,mp,d,dmean,radio
        real::dx,dy,dz,xv,yv,zv,r_void
        real::px,py,pz
        integer::part,bin,dimen,num
        dimen=0
        open(10,file='/mnt/is2/dpaz/ITV/S1373/halos/halos_50.ascii',status='unknown')
        do i=1,19
        read(10,*)
        enddo
        do i=1,halos
                read(10,*)a, num
                if  (num>1529) dimen=dimen+1
        enddo   
        close(10)
        
        allocate(id(dimen),num_p(dimen),rvir(dimen),pos(3,dimen))
        
        open(10,file='/mnt/is2/dpaz/ITV/S1373/halos/halos_50.ascii',status='unknown')
        j=0
        do i=1,19
        read(10,*)
        enddo
        do i=1,halos
                
                read(10,*) a,num,a,a,a,a,a,a,px,py,pz        
                if  (num>1529) then
                j=j+1
                num_p(j)=num
                pos(1,j)=px
                pos(2,j)=py
                pos(3,j)=pz 
                endif
        enddo   
        close(10) 
       ! call linkedlist(halos,abin,cell,pos,head,tot,link) 
        x0=254
        y0=247
        z0=250
       ! x0=252
       ! y0=245
       ! z0=250
        abin=(rmax-rmin)/real(bines)
        dmean=halos_tot/500**3
        r_void=0
        do j=1,100000
                cont=0
                do i=1,dimen
                d= sqrt((pos(1,i)-x0)**2+(pos(2,i)-y0)**2+(pos(3,i)-z0)**2)
                        if (d<rmax .and. d>rmin) then 
                        bin=int((d-rmin)/abin) + 1
                        cont(bin) = cont(bin) + 1  
                        endif    
                enddo 
                
                vol=0
                part=0
                D_r=0 
                do i=1,bines
                        radio=i*abin
                        vol=(4./3.)*pi*((i*abin)**3-rmin**3)
                        part= part + cont(i)
                        D_r(i) = (part/vol-dmean)/dmean
                 !       print*, D_r(i)
                        if (D_r(i)>=-0.9) cycle
                enddo
           !    print*, radio
                if (radio>r_void) then
                        write(*,*) 'CAMBIANDO VOID',j
                        r_void = radio
                        xv = x0
                        yv = y0
                        zv = z0
                endif
                call random_number(dx)        
                call random_number(dy)        
                call random_number(dz)
                x0 = x0 + 0.1*(0.5-dx)        
                y0 = y0 + 0.1*(0.5-dy)        
                z0 = z0 + 0.1*(0.5-dz)
                if (abs(x0-250)<10) x0 = 250
                if (abs(y0-250)<10) y0 = 250
                if (abs(z0-250)<10) z0 = 250
        enddo               
        write(*,*) radio,xv,yv,zv 
       stop 
endprogram IdEnTiFiCaToR
