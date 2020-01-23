program grid_mesh
        use modulos
        use OMP_lib
        implicit none
        integer, parameter :: cell=200, vec=5,vec2=50 ,nd=2000,bines=20
        real, parameter :: box=200, rmax=30, rmin=25! elegir el lado del box que voy a armar
        real, parameter :: distmin=0.001, distmax=1
        integer, dimension(cell,cell) :: tot, head
        !real*8,dimension(12232826) :: distance
        real, allocatable :: distance(:)
        real, allocatable :: dist(:)
        integer,allocatable,dimension(:,:,:) :: tot_gs, head_gs
        integer,allocatable,dimension(:,:,:) :: tot_dm, head_dm
        integer,allocatable,dimension(:,:,:) :: tot_st, head_st
        integer,allocatable :: link_gs(:)
        integer,allocatable :: link_dm(:)
        integer,allocatable :: link_st(:)
        real :: abin,xbox,ybox,zbox,xc,yc,zc,d
        integer :: ndm, ngs, nst,k,nst2
        real,allocatable :: pos_st2(:,:)
        real,allocatable :: pos_st(:,:), pos_gs(:,:), pos_dm(:,:)
        integer :: grid,ibin
        integer,dimension(bines) :: partbin_st, partbin_gs, partbin_dm
        real ::abin0,radio
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  ! NECESITO CREAR UN BOX MAS CHICO PARA NO TENER CELDAS AL PEDO AL LOS COSTADOS
  ! (DONDE ESTAN LAS TIDALES QUE NO ME INTERESAN) ENTONCES LO QUE HAGO ES CORRER
  ! LAS PARTICULAS DE LA RESIMULACION Y EL CENTRO DEL BOX 150 MPC PARA DE ESTA 
  ! MANERA TRABAJAR CON UN BOX DE 200 DE LADO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
 
        xbox=411.2170 
        ybox=162.1655 
        zbox=453.0553 
        xc=413.621475 - xbox + 250 -150  !estoy restando 150 para tener un box
        yc=162.604601 - ybox + 250 -150  !mas chico de la resi, de 150**3 (maso) 
        zc=448.953638 - zbox + 250 -150
 !xbox=403.8960 
 !ybox=459.8882
 !zbox=440.9021
 !xc=408.205481-xbox+250-150
 !yc=457.777839-ybox+250-150
 !zc=441.538681-zbox+250-150

        abin = box/real(cell)
        allocate(tot_gs(cell,cell,cell), head_gs(cell,cell,cell))
        allocate(tot_dm(cell,cell,cell), head_dm(cell,cell,cell))
        allocate(tot_st(cell,cell,cell), head_st(cell,cell,cell))
        call reader()
        allocate(pos_gs(3,nall(0)), pos_dm(3,nall(1)), pos_st(3,nall(4)))
        do i=1,nall(0)
        pos_gs(1,i)=pos(1,i)-150
        pos_gs(2,i)=pos(2,i)-150
        pos_gs(3,i)=pos(3,i)-150
        enddo   
        k=0
        do i=1+nall(0),nall(0)+nall(1)
        k=k+1
        pos_dm(1,k)=pos(1,i)-150
        pos_dm(2,k)=pos(2,i)-150
        pos_dm(3,k)=pos(3,i)-150
        enddo  
        k=0 
        do i=1+nall(0)+nall(1)+nall(2)+nall(3),nall(0)+nall(1)+nall(2)+nall(3)+nall(4)
        k=k+1
        !print*, pos(1,i)
        pos_st(1,k)=pos(1,i)-150
        pos_st(2,k)=pos(2,i)-150
        pos_st(3,k)=pos(3,i)-150
        enddo   
        deallocate(pos,vel)



        allocate(link_gs(nall(0)),link_dm(nall(1)),link_st(nall(4)))
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
        do i=1,nst
                d=sqrt((pos_st(1,i)-xc)**2+(pos_st(2,i)-yc)**2+(pos_st(3,i)-zc)**2)
                if (d<rmax .and. d>rmin) then
                        k=k+1
                endif
        enddo
        allocate(pos_st2(3,k))
        print*, 'stars particles in interests area:', k
        !dpares = k*
        nst2=k
        k=0
        do i=1,nst
                d=sqrt((pos_st(1,i)-xc)**2+(pos_st(2,i)-yc)**2+(pos_st(3,i)-zc)**2)
                if (d<rmax .and. d>rmin) then
                        k=k+1
                        pos_st2(1,k)=pos_st(1,i)
                        pos_st2(2,k)=pos_st(2,i)
                        pos_st2(3,k)=pos_st(3,i)
                endif
        enddo
        write(*,*) 'LiNKEANDO GAS'
        call linkedlist(ngs,abin,cell,pos_gs,head_gs,tot_gs,link_gs)
        write(*,*) 'LiNKEANDO DARMATE'
        call linkedlist(ndm,abin,cell,pos_dm,head_dm,tot_dm,link_dm)
        write(*,*) 'LiNKEANDO ESTRELLAs'
        call linkedlist(nst,abin,cell,pos_st,head_st,tot_st,link_st)
     !************************************************************   
     !*********** VECTOR DE DISTANCIAS ***************************
    ! allocate(distance(nst))
       write(*,*) 'VECINAS STARS'
        abin0 = (log10(distmax)-log10(distmin))/real(bines)
        print*, 'ESTRELLAS READY'
       call correlacion(nd,nst,pos_st,vec,nst2,pos_st2,abin,cell,box,tot_st,head_st,link_st,distmin,distmax, &
               bines,abin0,partbin_st)!,distance)
        print*, 'GAS READY'
       call correlacion(nd,ngs,pos_gs,vec2,nst2,pos_st2,abin,cell,box,tot_gs,head_gs,link_gs,distmin,distmax, &
               bines,abin0,partbin_gs)!,distance)
        print*, 'DM READY'
       call correlacion(nd,ndm,pos_dm,vec2,nst2,pos_st2,abin,cell,box,tot_dm,head_dm,link_dm,distmin,distmax, &
               bines,abin0,partbin_dm)!,distance)
     open(12,file='correlacion_st.dat',status='unknown')
        print*, abin0
     print*, distmin, distmax
     do i=1,bines
                radio = 10**(i*abin0 + log10(distmin))
                write(12,*) radio, partbin_st(i),partbin_gs(i),partbin_dm(i)

     enddo
     close(12)
      ! abin=box/real(cell)
       write(*,*) 'VECINAS GAS'
     !open(12,file='correlacion_gs.dat',status='unknown')
     !  call correlacion(nd,ngs,pos_gs,vec2,nst2,pos_st2,abin,cell,box,tot_gs,head_gs,link_gs)!,distance)
     !close(12)
     !  write(*,*) 'VECINAS DM'
     !open(12,file='correlacion_dm.dat',status='unknown')
     !  call correlacion(nd,ndm,pos_dm,vec2,nst2,pos_st2,abin,cell,box,tot_dm,head_dm,link_dm)!,distance)
     !close(12)
       
    !enddo
    ! deallocate(distance)
endprogram grid_mesh
subroutine correlacion(nd,ngs,pos2,vec,n,pos1,abin,cell,box,tot,head,link,distmin,distmax,bines,abin0,partbin)!,d)
        implicit none
        integer :: i,j,k,m,cand,bxx,byy,bzz,u,p,l
        integer :: bx,by,bz,n,cell,vec,nd,ngs
        real :: x0,y0,z0,x,y,z,r,abin,box,dx,dy,dz
        real,dimension(3,n):: pos1
        real,dimension(3,ngs):: pos2
        integer,dimension(ngs) :: link
      !  real,dimension(n) :: d
        integer,dimension(cell,cell,cell):: head, tot
        real,dimension(vec) :: dist
        integer :: vec0,o,q
        integer :: ibin,bines
        integer,dimension(bines) :: partbin
        real :: distmin,distmax,lmin,lmax,abin0
        cand = 0
        !d    = 0
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
        !write(*,*) ii
               ! dist=999999
               ! vec0=vec
                l=0
               ! cand=0
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
                        r=log10(r)
                       if (r>log10(distmin) .and. r<log10(distmax)) then
                        ibin=int((r-log10(distmin))/abin0)+1   
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
