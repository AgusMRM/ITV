program nose
implicit none
real,dimension(2,5) :: x
integer:: i
real:: med
logical,dimension(5) :: ilogic
print*, 'ingreso 5 valores x'
do i=1,5
        read(*,*) x(1,i)
enddo
do i=1,5
        read(*,*) x(2,i)
enddo
ilogic=.TRUE.
ilogic(1)=.FALSE.
print*, sum(x(2,:),mask=ilogic)

endprogram
