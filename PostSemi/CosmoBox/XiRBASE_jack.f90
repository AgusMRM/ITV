program grid_mesh
        use modulos
        use OMP_lib
        implicit none
        integer, parameter :: cell=10, vec=5,vec2=50, realjack=10
      !  integer,parameter  :: nd=2000,bines=20,halos=10682-19 !12317-19!
      !  integer,parameter  :: nd=2000,bines=20,halos=10682-19 !12317-19!
        integer,parameter  :: nd=2000,bines=20,halos=850366-16 !12317-19!
        real, parameter :: box=500,pi=acos(-1.), rmax=155, rmin=0! elegir el lado del box que voy a armar
        real*8, parameter :: distmin=1, distmax=30
        real,parameter :: boxt=125 
        integer*8,parameter :: estrellas=10439184 , darkmatter=134217728 , &
                            gas=134217728 - 10439184

        integer, dimension(cell,cell) :: tot, head
        real, allocatable :: distance(:)
        real, allocatable :: dist(:)
        integer,allocatable,dimension(:,:,:) :: tot_gs, head_gs
        integer,allocatable,dimension(:,:,:) :: tot_dm, head_dm
        integer,allocatable,dimension(:,:,:) :: tot_st, head_st
        integer,allocatable,dimension(:,:,:) :: tot_hl, head_hl
        integer,allocatable :: link_gs(:)
        integer,allocatable :: link_dm(:)
        integer,allocatable :: link_st(:)
        integer,allocatable :: link_hl(:)
        real :: abin,xbox,ybox,zbox,xc,yc,zc,d,rand,vs
        real :: dpares_st, dpares_gs, dpares_dm, dpares_hl
        integer :: ndm, ngs, nst,k,trazers_st,trazers_hl,ijack,trazers
        real,allocatable :: pos_st2(:,:),pos_hl2(:,:)
        real,allocatable :: pos_st(:,:), pos_gs(:,:), pos_dm(:,:),pos_hl(:,:)
        integer :: grid,ibin
        integer*8,dimension(bines,realjack) :: partbin_st, partbin_gs, partbin_dm
        integer*8,dimension(bines,realjack) :: partbin_hl
        real*8 ::abin0,radio, radio0
        integer*8,allocatable :: binpart(:,:)
        integer, allocatable,dimension(:)  :: jack_dm, jack_gs, jack_st, jack_hl
        real,dimension(bines) :: pares_teoricos_st, pares_teoricos_dm,pares_teoricos_gs,pares_teoricos_hl
        real,dimension(bines) :: u_st,u_dm,u_gs,u_hl, var_st,var_dm,var_gs, var_hl
        real,dimension(bines,realjack) :: xi_st,xi_dm,xi_gs, xi_hl
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  ! NECESITO CREAR UN BOX MAS CHICO PARA NO TENER CELDAS AL PEDO AL LOS COSTADOS
  ! (DONDE ESTAN LAS TIDALES QUE NO ME INTERESAN) ENTONCES LO QUE HAGO ES CORRER
  ! LAS PARTICULAS DE LA RESIMULACION Y EL CENTRO DEL BOX 150 MPC PARA DE ESTA 
  ! MANERA TRABAJAR CON UN BOX DE 200 DE LADO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
 
  !      xbox=411.2170 
  !      ybox=162.1655 
  !      zbox=453.0553 
  !      xc=413.621475 - xbox + 250  !estoy restando 150 para tener un box
  !      yc=162.604601 - ybox + 250  !mas chico de la resi, de 150**3 (maso) 
  !      zc=448.953638 - zbox + 250
 xbox=403.8960 
 ybox=459.8882
 zbox=440.9021
 xc=408.205481-xbox+250
 yc=457.777839-ybox+250
 zc=441.538681-zbox+250

xc=0
yc=0
zc=0

        abin = box/real(cell)
        allocate(tot_gs(cell,cell,cell), head_gs(cell,cell,cell))
        allocate(tot_dm(cell,cell,cell), head_dm(cell,cell,cell))
        allocate(tot_st(cell,cell,cell), head_st(cell,cell,cell))
       allocate(tot_hl(cell,cell,cell), head_hl(cell,cell,cell))
        call reader()
       allocate(pos_gs(3,nall(0)), pos_dm(3,nall(1)), &
               pos_st(3,nall(4)),pos_hl(3,halos))
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
        !print*, pos(1,i)
        pos_st(1,k)=pos(1,i)
        pos_st(2,k)=pos(2,i)
        pos_st(3,k)=pos(3,i)
        enddo   
        deallocate(pos,vel)



        allocate(link_gs(nall(0)),link_dm(nall(1)),link_st(nall(4)),link_hl(halos))
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
        k=10000
        allocate(pos_st2(3,k))
        allocate(binpart(k,bines))
        print*, 'stars particles in interests area:', k
        !dpares = k*
        trazers_st=k
        k=0
 !       do i=1,nst
 !          !     call random_number(rand)
 !             !  if (rand < 0.05 .and. k<10000) then
 !              d=sqrt((pos_st(1,i)-xc)**2+(pos_st(2,i)-yc)**2+(pos_st(3,i)-zc)**2)
 !              if (d<rmax .and. d>rmin) then
 !                       k=k+1
 !                       pos_st2(1,k)=pos_st(1,i)
 !                       pos_st2(2,k)=pos_st(2,i)
 !                       pos_st2(3,k)=pos_st(3,i)
 !               endif
 !       enddo
! AHORA SELECCIONO LOS HALOS
open(123,file='/mnt/is2/fstasys/ITV/base09/rockstar/out_0.list',status='unknown')
!open(123,file='/home/arodriguez/halosprueba/halos_0.0.ascii',status='unknown')
        do i=1,16
        read(123,*)
        enddo
        do i=1,halos
                read(123,*) a,a,a,a,a,a,a,a,pos_hl(1,i),pos_hl(2,i),pos_hl(3,i)
        enddo
  close(123)
        k=0
        do i=1,halos
                d=sqrt((pos_hl(1,i)-xc)**2+(pos_hl(2,i)-yc)**2+(pos_hl(3,i)-zc)**2)
                if (d<rmax .and. d>rmin) then
                        k=k+1
                endif
        enddo
        k=10000
        allocate(pos_hl2(3,k))
       ! allocate(binpart(k,bines))
        print*, 'halos particles in interests area:', k
        !dpares = k*
        trazers_hl=k
        k=0
        do i=1,halos
                call random_number(rand)
                if (rand < 0.05 .and. k<10000) then
                d=sqrt((pos_hl(1,i)-xc)**2+(pos_hl(2,i)-yc)**2+(pos_hl(3,i)-zc)**2)
              !  if (d<rmax .and. d>rmin) then
                        k=k+1
                        pos_hl2(1,k)=pos_hl(1,i)
                        pos_hl2(2,k)=pos_hl(2,i)
                        pos_hl2(3,k)=pos_hl(3,i)
                endif
        enddo        
        write(*,*) 'LiNKEANDO GAS'
   !     call linkedlist(ngs,abin,cell,pos_gs,head_gs,tot_gs,link_gs)
        write(*,*) 'LiNKEANDO DARMATE'
   !     call linkedlist(ndm,abin,cell,pos_dm,head_dm,tot_dm,link_dm)
        write(*,*) 'LiNKEANDO ESTRELLAs'
   !     call linkedlist(nst,abin,cell,pos_st,head_st,tot_st,link_st)
        call linkedlist(halos,abin,cell,pos_hl,head_hl,tot_hl,link_hl)
     !************************************************************   
     !*********** VECTOR DE DISTANCIAS ***************************
    ! allocate(distance(nst))
       write(*,*) 'VECINAS STARS'
        
       
       abin0 = (log10(distmax)-log10(distmin))/real(bines)
   ! ARMO LOS VECTORES JACK, ESTOS TIENE ENTEROS DE ENTRE 1 Y EL NUMERO DE
   ! REALIZACIONES QUE VAYA A HACER  
  allocate(jack_gs(ngs), jack_dm(ndm),jack_st(nst), jack_hl(halos))   
   do i=1,ngs
        call random_number(rand)
        jack_gs(i)=int(rand*realjack) + 1
   enddo 
   do i=1,ndm
        call random_number(rand)
        jack_dm(i)=int(rand*realjack) + 1
   enddo 
   do i=1,nst
        call random_number(rand)
        jack_st(i)=int(rand*realjack) + 1
   enddo 
   do i=1,halos
        call random_number(rand)
        jack_hl(i)=int(rand*realjack) + 1
   enddo
 !##########################################################################33   
 !##########################################################################33   
 !##########################################################################33   
  partbin_gs=0; partbin_dm=0; partbin_st=0; partbin_hl=0 
 ! binpart_dm=0; binpart_gs=0; binpart_st=0; binpart_hl=0
    trazers = trazers_hl
    
  do ijack=1,realjack 
   write(*,*) 'realizacion:', ijack
    !  call correlacion(nst,pos_st,vec,trazers,pos_hl2,abin,cell,box,tot_st,head_st,link_st,distmin,distmax, &
    !         bines,abin0,partbin_st,realjack,jack_st,ijack)!,distance)
    !  call correlacion(ngs,pos_gs,vec2,trazers,pos_hl2,abin,cell,box,tot_gs,head_gs,link_gs,distmin,distmax, &
    !          bines,abin0,partbin_gs,realjack,jack_gs,ijack)!,distance)

   !   call correlacion(ndm,pos_dm,vec2,trazers,pos_hl2,abin,cell,box,tot_dm,head_dm,link_dm,distmin,distmax, &
   !          bines,abin0,partbin_dm,realjack,jack_dm,ijack)!,distance)
      call correlacion(halos,pos_hl,vec2,trazers,pos_hl2,abin,cell,box,tot_hl,head_hl,link_hl,distmin,distmax, &
          bines,abin0,partbin_hl,realjack,jack_hl,ijack)!,distance)
   enddo  


 !##########################################################################33   
 !##########################################################################33   
 !##########################################################################33   
    dpares_st=trazers*estrellas/(boxt**3) 
    dpares_gs=trazers*darkmatter/(boxt**3) 
    dpares_dm=trazers*gas/(boxt**3) 
    dpares_hl=trazers*real(halos)/(real(box)**3)
   write(*,*) 'trazers:', trazers 
   write(*,*) 'dpares:', dpares_hl
  radio0 = distmin 
    do i=1,bines
                radio = 10**(i*abin0 + log10(distmin))
              !  vs=(4./3.)*pi*(radio**3-distmin**3)
                vs=(4./3.)*pi*(radio**3-radio0**3)
                pares_teoricos_st(i)=dpares_st*vs
                pares_teoricos_dm(i)=dpares_dm*vs
                pares_teoricos_gs(i)=dpares_gs*vs
                pares_teoricos_hl(i)=dpares_hl*vs
                radio0 = radio
    enddo 
   
   do i=1,bines
        write(*,*) partbin_hl(i,1)
   enddo
   do i=1,bines
        write(*,*) 'teoricos', pares_teoricos_hl(i)
   enddo
   
    do j=1,realjack
        do i=1,bines
                xi_st(i,j)=partbin_st(i,j)/pares_teoricos_st(i) - 1
                xi_gs(i,j)=partbin_gs(i,j)/pares_teoricos_gs(i) - 1
                xi_dm(i,j)=partbin_dm(i,j)/pares_teoricos_dm(i) - 1
                xi_hl(i,j)=real(partbin_hl(i,j))/pares_teoricos_hl(i) - 1
        enddo 
    enddo    
! CALCULO LA MEDIA JACKNIFE DE CADA BIN        
u_st = 0
u_dm = 0
u_gs = 0
u_hl = 0
   do i=1,bines
        do j=1,realjack
                u_st(i) = u_st(i) + xi_st(i,j)
                u_dm(i) = u_dm(i) + xi_dm(i,j)
                u_gs(i) = u_gs(i) + xi_gs(i,j)
                u_hl(i) = u_hl(i) + xi_hl(i,j)
        enddo
        u_st(i)=u_st(i)/realjack
        u_dm(i)=u_dm(i)/realjack
        u_gs(i)=u_gs(i)/realjack
        u_hl(i)=u_hl(i)/real(realjack)
   enddo
! CALCULO LA VARIANZA
var_st = 0
var_dm = 0
var_gs = 0
var_hl = 0
   do i=1,bines
        do j=1,realjack
                var_st(i) = (xi_st(i,j) - (u_st(i)))**2
                var_dm(i) = (xi_dm(i,j) - (u_dm(i)))**2
                var_gs(i) = (xi_gs(i,j) - (u_gs(i)))**2
                var_hl(i) = (xi_hl(i,j) - (u_hl(i)))**2
              !  var_st(i) = xi_st(i,j)**2 - (u_st(i)**2)*realjack
              !  var_dm(i) = xi_dm(i,j)**2 - (u_dm(i)**2)*realjack
              !  var_gs(i) = xi_gs(i,j)**2 - (u_gs(i)**2)*realjack
        enddo
        
               var_st(i) = ((real(realjack)-1.)/real(realjack))*var_st(i) 
               var_dm(i) = ((real(realjack)-1.)/real(realjack))*var_dm(i)
               var_gs(i) = ((real(realjack)-1.)/real(realjack))*var_gs(i)
               var_hl(i) = ((real(realjack)-1.)/real(realjack))*var_hl(i)
   enddo


     open(12,file='correlacion.dat',status='unknown')
     do i=1,bines
                radio = 10**(i*abin0 + log10(distmin))
!               write(12,*) radio, partbin_st(i,1),partbin_gs(i,1),partbin_dm(i,1), &
!                               partbin_hl(i,2)
        write(12,*) radio, u_hl(i), u_dm(i), u_gs(i), var_hl(i), var_dm(i), var_gs(i)
                
     enddo
     close(12)
endprogram grid_mesh
subroutine correlacion(npart,pos2,vec,n,pos1,abin,cell,box,tot,head,link,distmin,distmax,bines,abin0,partbin,&
                realjack,jack,ijack)!,d)
        implicit none
        integer :: i,j,k,m,cand,bxx,byy,bzz,u,p,l,realjack,ijack
        integer :: bx,by,bz,n,cell,vec,nd,npart
        real*8 :: x0,y0,z0,x,y,z,r
        real :: abin,box,dx,dy,dz
        real,dimension(3,n):: pos1
        real,dimension(3,npart):: pos2
        integer,dimension(npart) :: link,jack
        integer,dimension(cell,cell,cell):: head, tot
        real,dimension(vec) :: dist
        integer*8,dimension(n,bines) ::binpart
        integer :: vec0,o,q
        integer :: ibin,bines
        integer*8 ::pares
        integer*8,dimension(bines,realjack) :: partbin
        real*8 :: distmin,distmax,lmin,lmax,abin0
        cand = 0
        l = 0 
        lmin=log10(distmin)
        lmax=log10(distmax)
        pares=0
        binpart=0
  !$OMP PARALLEL DEFAULT (NONE) &
  !$OMP PRIVATE (i,bx,by,bz,x,y,z,j,k,m,u,p,bxx,byy,bzz, &
  !$OMP cand,dist,l,x0,y0,z0,dx,dy,dz,r,vec0,o,q,ibin) &
  !$OMP SHARED (pos1,abin,tot,head,link,n,cell,vec,box,nd,pos2,lmin , &
  !$OMP lmax,bines,partbin,abin0,distmin,distmax,binpart,ijack,jack)
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
                        if (dx > box/2.) dx = box - dx
                        if (dy > box/2.) dy = box - dy
                        if (dz > box/2.) dz = box - dz

                        r  = sqrt(dx**2+dy**2+dz**2)
                        r=log10(r)
                       if (r>log10(distmin) .and. r<log10(distmax) .and. &
                               jack(p) /= ijack ) then
                        ibin=int((r-log10(distmin))/abin0)+1  
                        binpart(i,ibin)=binpart(i,ibin)+1 
!                        partbin(ibin)=partbin(ibin)+1                 
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
   do i=1,n
        do j=1,bines
                pares = binpart(i,j)
                partbin(j,ijack)=partbin(j,ijack) + pares
        enddo
   enddo
    
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
