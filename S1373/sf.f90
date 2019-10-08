program SF
       use modulos
       use OMP_lib 
       implicit none
       integer,parameter :: VOID = 1373, nthread = 6
       integer :: SNAPSHOT,c1,c2,k     
       real :: d,xbox,ybox,zbox,x,y,z,rv
       character(len=200) :: snumber
        xbox=403.8960
        ybox=459.8882
        zbox=440.9021 
        open(13,file='voids_new.dat')
        do i=1,VOID-1
        read(13,*)
        enddo
        read(13,*) rv,x,y,z
        x=x-xbox+250
        y=y-ybox+250
        z=z-zbox+250
        
       open(10,file='sfh_in.dat.out',status='unknown')
       open(11,file='sfh_wall.dat.out',status='unknown')
       
      call OMP_set_num_threads(nthread)
       do j=30,50
       c1=0
       c2=0
       SNAPSHOT=j
       write(*,*) snapshot
       write(snumber,'(I3)') SNAPSHOT
       call reader(snumber)
        
!$OMP PARALLEL DEFAULT(NONE) &       
!$OMP PRIVATE( k, d) &
!$OMP SHARED(x,y,z,pos,c2,nall)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = nall(0)+nall(1)+nall(2)+1, nall(0)+nall(1)+nall(2)+nall(4)
               d=sqrt((pos(k,1)-x)**2 + (pos(k,2)-y)**2 + (pos(k,3)-z)**2)
         !       if ( d<5 ) c1=c1+1
                if ( d>5 .and. d <10) c2=c2+1
        enddo
!$OMP END DO      
!$OMP END PARALLEL   
        print*, c1, c2
        !write(10,*) c1, redshift
        write(11,*) c2, redshift
        deallocate(pos)
   
    enddo     
   close(10) 
   close(11) 
endprogram 
