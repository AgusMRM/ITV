program df
        use modulos
        implicit none
        real :: rr,tt,dd
        call reader()
        open(10,file='diagfases_random.dat',status='unknown')   
        open(11,file='diagfases.dat',status='unknown')
        do i=1,nall(0)
                read(11,*) tt,dd
                call random_number(rr)
                if (rr<.1) write(10,*) tt,dd
        enddo
        close(10)
        close(11)
        

endprogram
