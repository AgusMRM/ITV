program MUSICREADER
        implicit none
        integer :: i
        real :: np1,np2,np3,iseed
        character(len=200), parameter :: path='/mnt/is2/fstasys/ITV/ICs/Exclusive/'
        character(len=200) :: filename

        filename=trim(path)//'wnoise_0008.bin'
        open(11,file=filename,form='unformatted')
        read(11,*) 
        close(11)
endprogram
