program grid_mesh
        use modulos
        use OMP_lib
        implicit none
        integer, parameter :: cell=20, vec=5,vec2=50 ,nd=2000,bines=20, halos=71417-19
        real, parameter :: box=500, rmax=20, rmin=0.0! elegir el lado del box que voy a armar
       ! real, parameter :: distmin=0.001, distmax=4
        real, parameter :: distmin=0.5, distmax=5
        integer, dimension(cell,cell) :: tot, head
        !real*8,dimension(12232826) :: distance
        real, allocatable :: distance(:)
        real, allocatable :: dist(:)
        integer,allocatable,dimension(:,:,:) :: tot_gs, head_gs
        integer,allocatable,dimension(:,:,:) :: tot_dm, head_dm
        integer,allocatable,dimension(:,:,:) :: tot_st, head_st
        integer,allocatable,dimension(:,:,:) :: tot_hl, head_hl
        integer,allocatable :: link_dm(:)
        integer,allocatable :: link_st(:)
        integer,allocatable :: link_gs(:)
        integer,allocatable :: link_hl(:)
        real :: abin,xbox,ybox,zbox,xc,yc,zc,d
        integer :: ndm, ngs, nst,k,nst2
        real,allocatable :: pos_hl2(:,:)
        real,allocatable :: pos_dm(:,:), pos_gs(:,:), pos_st(:,:)
        real,allocatable :: pos_hl(:,:)
        integer :: grid,ibin
        integer,dimension(bines) :: partbin_st, partbin_gs, partbin_dm
        real ::abin0,radio
        partbin_st=0
        partbin_gs=0
        partbin_dm=0
 xbox=403.8960 
 ybox=459.8882
 zbox=440.9021
 xc=408.205481-xbox+250
 yc=457.777839-ybox+250
 zc=441.538681-zbox+250
 
        abin = box/real(cell)
        allocate(tot_gs(cell,cell,cell), head_gs(cell,cell,cell))
        allocate(tot_dm(cell,cell,cell), head_dm(cell,cell,cell))
        allocate(tot_st(cell,cell,cell), head_st(cell,cell,cell))
        allocate(tot_hl(cell,cell,cell), head_hl(cell,cell,cell))
        allocate(pos_hl(3,halos))
        
        open(21,file='/mnt/is2/dpaz/ITV/S1373/halos/halos_50.ascii',status='unknown')
        do i=1,19
        read(21,*)
        enddo
        do i=1,halos
                read(21,*) a,a,a,a,a,a,a,a,a,pos_hl(1,i),pos_hl(2,i),pos_hl(3,i)
        enddo
        close(21)
        print*, 'ya lei'

        do i=1,halos
                pos_hl(1,i)=pos_hl(1,i)
                pos_hl(2,i)=pos_hl(2,i)
                pos_hl(3,i)=pos_hl(3,i)
        enddo
        
        print*, 'ya allocatie las cosas para la linked list'
        call reader()
        allocate(pos_gs(3,nall(0)), pos_dm(3,nall(1)), pos_st(3,nall(4)))
        do i=1,nall(0)
        pos_gs(1,i)=pos(1,i)
        pos_gs(2,i)=pos(2,i)
        pos_gs(3,i)=pos(3,i)
        enddo   
        k=0
        do i=1+nall(0),nall(0)+nall(1)
        k=k+1
        pos_dm(1,k)=pos(1,i)
        pos_dm(2,k)=pos(2,i)
        pos_dm(3,k)=pos(3,i)
        enddo  
        k=0 
        do i=1+nall(0)+nall(1)+nall(2)+nall(3),nall(0)+nall(1)+nall(2)+nall(3)+nall(4)
        k=k+1
        pos_st(1,k)=pos(1,i)
        pos_st(2,k)=pos(2,i)
        pos_st(3,k)=pos(3,i)
        enddo   
        deallocate(pos,vel)



        allocate(link_gs(nall(0)),link_dm(nall(1)),link_st(nall(4)))
        allocate(link_hl(halos))
        ngs=nall(0)
        ndm=nall(1)
        nst=nall(4)

        tot_gs  = 0
        head_gs = 0
        link_gs = 0
        tot_dm  = 0       
        head_dm = 0
        link_dm = 0
        tot_st  = 0
        head_st = 0
        link_st = 0
        
        k=0
        do i=1,halos
                d=sqrt((pos_hl(1,i)-xc)**2+(pos_hl(2,i)-yc)**2+(pos_hl(3,i)-zc)**2)
                if (d<rmax .and. d>rmin) then
                        k=k+1
                endif
        enddo
        print*, 'tracers:', k
        allocate(pos_hl2(3,k))
        nst2=k
        k=0
        
        do i=1,halos
                d=sqrt((pos_hl(1,i)-xc)**2+(pos_hl(2,i)-yc)**2+(pos_hl(3,i)-zc)**2)
                if (d<rmax .and. d>rmin) then
                        k=k+1
                        pos_hl2(1,k)=pos_hl(1,i)
                        pos_hl2(2,k)=pos_hl(2,i)
                        pos_hl2(3,k)=pos_hl(3,i)
                endif
        enddo

        write(*,*) 'LiNKEANDO DARK MATER'
     !   call linkedlist(ndm,abin,cell,pos_dm,head_dm,tot_dm,link_dm)
        call linkedlist(halos,abin,cell,pos_hl,head_hl,tot_hl,link_hl)
       
        abin0 = (log10(distmax)-log10(distmin))/real(bines)
        call correlacion(halos,pos_hl,vec,nst2,pos_hl2,abin,cell,box,tot_hl,head_hl,link_hl,distmin,distmax, &
              bines,abin0,partbin_dm)!,distance)
     
        print*, 'DM READY'
     
        open(22,file='correlacion_st.dat',status='unknown')
        print*, distmin, distmax
        do i=1,bines
                radio = 10**(i*abin0 + log10(distmin))
                write(22,*) radio, partbin_st(i),partbin_gs(i),partbin_dm(i)

        enddo
        close(22)
endprogram grid_mesh
subroutine correlacion(ngs,pos2,vec,n,pos1,abin,cell,box,tot,head,link,distmin,distmax,bines,abin0,partbin)!,d)
        implicit none
        integer :: i,j,k,m,cand,bxx,byy,bzz,u,p,l
        integer :: bx,by,bz,n,cell,vec,nd,ngs
        real :: x0,y0,z0,x,y,z,abin,box
        real*8 :: r,dx,dy,dz
        real,dimension(3,n):: pos1
        real,dimension(3,ngs):: pos2
        integer,dimension(ngs) :: link
        integer,dimension(cell,cell,cell):: head, tot
        real,dimension(vec) :: dist
        integer :: vec0,o,q
        integer :: ibin,bines
        integer,dimension(bines) :: partbin
        real :: distmin,distmax,lmin,lmax,abin0
        cand = 0
        l = 0 
        lmin=log10(distmin)
        lmax=log10(distmax)
        partbin=0
  !$OMP PARALLEL DEFAULT (NONE) &
  !$OMP PRIVATE (i,bx,by,bz,x,y,z,j,k,m,u,p,bxx,byy,bzz, &
  !$OMP cand,dist,l,x0,y0,z0,dx,dy,dz,r,vec0,o,q,ibin) &
  !$OMP SHARED (pos1,abin,tot,head,link,n,cell,vec,box,nd,pos2,lmin , &
  !$OMP lmax,bines,partbin,abin0)
  !$OMP DO SCHEDULE (DYNAMIC)
        do i=1,n
                l=0
                x = pos1(1,i)
                y = pos1(2,i)
                z = pos1(3,i)
                bx = int(x/abin) + 1 
                by = int(y/abin) + 1 
                bz = int(z/abin) + 1
                
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1
                if (bz == 0)      bz=cell
                if (bz == cell+1) bz=1
                do j= bx-1,bx+1
                 if (j == 0) then 
                        bxx=cell
                 elseif (j == cell+1) then
                        bxx=1
                 else
                        bxx=j
                 endif      
                 do k= by-1,by+1
                    if (k == 0) then
                        byy=cell
                    elseif (k == cell+1) then
                        byy=1
                    else
                        byy=k
                    endif  
                 do m= bz-1,bz+1
                    if (m == 0) then
                        bzz=cell
                    elseif (m == cell+1) then
                        bzz=1
                    else
                        bzz=m
                    endif
           
                     p  = head(bxx,byy,bzz)
                    do u=1,tot(bxx,byy,bzz)
                        l  = l + 1
                        x0 = pos2(1,p)
                        y0 = pos2(2,p)
                        z0 = pos2(3,p)
                        dx=abs(x-x0)
                        dy=abs(y-y0)
                        dz=abs(z-z0)
                       ! if (dx > box/2.) dx = box - dx
                       ! if (dy > box/2.) dy = box - dy
                       ! if (dz > box/2.) dz = box - dz

                        r  = sqrt(dx**2+dy**2+dz**2)
                       ! r=log10(r)
                       if (r>distmin .and. r<distmax) then
                        ibin=int((log10(r)-log10(distmin))/abin0)+1  
      !                  if (ibin > 20) print*, r
      !                  if (ibin > 20) go to 4

                        partbin(ibin)=partbin(ibin)+1                 
                       endif 
                        p = link(p)
                    enddo
                 enddo    
                  enddo    
                enddo
                  l = 0
        enddo
   !$OMP END DO
   !$OMP BARRIER
   !$OMP END PARALLEL
endsubroutine        
subroutine ordenador(r,vec,dist)
    implicit none
    ! integer,parameter :: vec=3,n=20 
    integer ::i,o,q,vec,activador
    real,dimension(vec) ::dist
    real :: r,dist_vieja, dist_vieja2
   ! dist=99999
    activador=1
    if (r<dist(vec)) then
            do o=1,vec-1      
                if (r>dist(vec-o) .and. r<= dist(vec+1-o)) then   
                   activador=0     
                    dist_vieja=dist(vec-o)
                    dist(vec-o)=r
                    do q=vec-o,vec
                    dist_vieja2=dist(q)
                    dist(q)=dist_vieja
                    dist_vieja=dist_vieja2
                    enddo
                    goto 4
                endif
            enddo
      4      if (activador==1) then
                    dist_vieja=dist(1)
                dist(1)=r    
                    do q=2,vec
                    dist_vieja2=dist(q)
                    dist(q)=dist_vieja
                    dist_vieja=dist_vieja2
                    enddo
            endif    
    endif   
endsubroutine
SUBROUTINE ORDEN(NELEM,ARREG)
! -----------------------------------------------------
!ORDENACION POR BURBUJA ("buble sort") de un arreglo
! unidimensional, de menor a mayor.
!
! NELEM = Número de elementos del arreglo
! ARREG = Arreglo unidimensional a ordenar
! -----------------------------------------------------
IMPLICIT NONE
INTEGER :: NELEM
REAL :: ARREG(*)
!-----------------------------------------------------
INTEGER:: I,J
REAL:: AUX
!-----------------------------------------------------
IF (NELEM.LT.2) RETURN
DO I=1,NELEM-1
DO J=1,NELEM-I
IF (ARREG(J).GT.ARREG(J+1)) THEN
AUX = ARREG(J)
ARREG(J) = ARREG(J+1)
ARREG(J+1) = AUX
ENDIF
ENDDO
ENDDO
RETURN
END
SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
