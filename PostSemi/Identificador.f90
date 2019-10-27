program IdEnTiFiCaToR
        implicit none
        integer, parameter :: halos=60000,bines=100
        character(len=200) :: filename
        real,parameter :: pi=acos(-1.)
        real ::a
        integer ::i,j
        integer,dimension(halos)::id,num_p
        real,dimension(halos)::rvir
        real,dimension(3,halos)::pos
        open(10,file='halos_S1373.dat',status='unknown')
        do i=1,19
        read(10,*)
        enddo
        do i=1,halos
                read(10,*)id(i), num_p(i),a, a, rvir(i), a,a,a,pos(1,i),pos(2,i),pos(3,i) 
        enddo   

        x0=250
        y0=250
        z0=250
        r=0
        
        do i=1,1000
                
        enddo      
        

        close(10) 
        
endprogram
