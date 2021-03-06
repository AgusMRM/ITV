!SUBRUTINAS UTILIZADAS PARA LA RE-IDENTIFICACION DE LOS VOIDS
!SOBRE LOS HALOS
subroutine linkedlist(n,cmin,abin,cell,pos,head,tot,link)
        implicit none
        integer :: i,n,bx,by,bz,cell
        integer,dimension(cell,cell,cell) :: tot, head
        real,dimension(3,n) :: pos
        integer,dimension(n) :: link
        real :: abin,cmin
        
        do i = 1,n
               bx = int((pos(1,i)-cmin)/abin) + 1
               by = int((pos(2,i)-cmin)/abin) + 1
               bz = int((pos(3,i)-cmin)/abin) + 1
               if (bx==0)     bx=cell
               if (bx==cell+1)bx=1
               if (by==0)     by=cell
               if (by==cell+1)by=1
               if (bz==0)     bz=cell
               if (bz==cell+1)bz=1

               tot(bx,by,bz) = tot(bx,by,bz) + 1
               head(bx,by,bz)   = i
       enddo

       do i = 1,n
               bx = int((pos(1,i)-cmin)/abin) + 1
               by = int((pos(2,i)-cmin)/abin) + 1
               bz = int((pos(3,i)-cmin)/abin) + 1

                if (bx==0)     bx=cell
                if (bx==cell+1)bx=1
                if (by==0)     by=cell
                if (by==cell+1)by=1
                if (bz==0)     bz=cell
                if (bz==cell+1)bz=1
               
                link(head(bx,by,bz)) = i
                head(bx,by,bz)       = i
        enddo 
endsubroutine 
