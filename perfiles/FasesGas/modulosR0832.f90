module modulos
   implicit none
   integer, parameter :: SNAPSHOT = 50     ! number of dump
   integer, parameter :: FILES = 1         ! number of files per snapshot
!   character(len=200), parameter :: path='/mnt/is2/dpaz/ITV/S1373/out'
   character(len=200), parameter :: path='/mnt/is2/fstasys/ITV/R0832/out'

   character(len=200) :: filename, snumber, fnumber

   integer(4),dimension(0:5) :: npart, nall
   real(8),dimension(0:5)    :: massarr
   real*8    a
   real*8    redshift
   integer*4 fn,i,j,nstart,flag_sfr,flag_feedback
   integer*4 N,Ntot
   integer(4),dimension(0:5) :: Newnpart

   integer*4 flag_cool, num_files
   real*8    boxsize
   real*8 newbox
   real*8    omega_m,omega_l
   integer*4 unused(26)

   real*4,allocatable    :: pos(:,:),vel(:,:)
   integer*4,allocatable :: id(:),idch(:),idgn(:)
   real*4, allocatable,dimension(:)    :: mass, u, rho, ne, &
                                        nh,hsml,sfr,abvc,temp 
   real*8 n0,mol,R,vx,vy,vz
   real :: prom1,prom2,prom3
       
   character(len=4) :: blckname 
   integer(kind=4) :: hwm
   real(kind=4)  :: hwm2 
   real :: xmin,ymin,zmin
   real :: xmax,ymax,zmax
   integer(4) :: i1,i2

   logical :: ilogic
endmodule
subroutine reader()
        use modulos
   write(snumber,'(I3)') SNAPSHOT
   write(fnumber,'(I3)') 1
   do i=1,3 
      if (snumber(i:i) .eq. ' ') then 
          snumber(i:i)='0'
      end if
   end do

   !filename= path(1:len_trim(path)) // '/snapshot_' // snumber(1:len_trim(snumber)) !!// '.' // fnumber(verify(fnumber,' '):3)
   filename= trim(path)//'/snapshot_'//trim(snumber) !!// '.' // fnumber(verify(fnumber,' '):3)

   print *,'opening...  '//filename

   ! now, read in the header

   open (1, file=filename, access='stream',action='read')


   read(1)hm1
   read(1)blckname,hwm
   read(1)hm2
   if(hm1/=hm2) stop 'apa la papa'
   write(*,*)'leyendo ',blckname,hwm

 
   read(1)hm1
   read (1) npart, massarr, a, redshift, flag_sfr,flag_feedback, nall,flag_cool, num_files, boxsize,omega_m,omega_l,unused
   read(1)hm2
   if(hm1/=hm2) stop 'apa la papa'
   write(*,*) a
   write(*,*) 'redshift:', redshift
   write(*,*) nall


   do i=0,5
      print *,'Type=',i,'    Particles=', nall(i)
   end do
   ! Now we just read in the coordinates, and only keep those of type=1
   ! The array PartPos(1:3,1:Ntot) will contain all these particles
   Ntot= sum(nall)
   write(*,*)"---------> ", Ntot

   allocate(pos(3,nall(0)+nall(1)+nall(2)+nall(4)))
   allocate(vel(3,Ntot))
   !allocate(id(Ntot))
   !allocate(mass(Ntot))
   !allocate(u(nall(0))) !+nall(4)))
   !allocate(dens(nall(0))) !+nall(4)))
   N=sum(npart)

   print*, (3*sum(nall)*4), '<--------- 3*sum(nall)*4 '
  print*, 'Ntot', ntot 
   !hwm contiene el numero de bytes del bloque de datos mas 8 bytes
  
  
   read(1) hm1 
   read (1)blckname,hwm
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   print*, blckname, hwm
   read(1) hm1 
   read (1)pos
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   write(*,*)'leyendo ',blckname,hwm-8
   
   
   
   read(1) hm1
   read (1)blckname,hwm
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   write(*,*)'leyendo ',blckname,hwm-8
   read(1) hm1
   read (1) vel
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 

   read(1) hm1
   read (1)blckname,hwm
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   write(*,*)'leyendo ',blckname,hwm-8
   read(1) hm1
   allocate(id(ntot))
   read (1)id
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
 
   read(1) hm1 
   read (1)blckname,hwm
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   write(*,*)'leyendo ',blckname,hwm-8
  
   read(1) hm1 
   allocate(idch(ntot))
   read (1)idch
   read(1) hm2 
   read(1) hm1 
   read (1)blckname,hwm
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 
   
   read(1) hm1 
   write(*,*)'leyendo ',blckname,hwm-8
   allocate(idgn(ntot))
   read (1)idgn
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 

   read(1) hm1 
   read (1)blckname,hwm
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 
   write(*,*)'leyendo ',blckname,hwm-8
   allocate(mass(ntot))
   read(1) hm1 
   read (1)mass
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 

   read(1) hm1 
   read (1)blckname,hwm
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 
   write(*,*)'leyendo ',blckname,hwm-8
   allocate(u(nall(0)))
   !write(*,*) 'allocate con', int(real(hwm-8)/4.)
   read(1) hm1 
   read (1)u
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 
   
   read(1) hm1 
   read (1)blckname,hwm
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 
   write(*,*)'leyendo ',blckname,hwm-8
   allocate(rho(nall(0)))
   read(1) hm1 
   read (1)rho
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 
   
   read(1) hm1 
   read (1)blckname,hwm
   read(1) hm2 
   if (hm1/=hm2) stop 'apalapapa' 
   allocate(ne(nall(0)))
   read(1) hm1 
   read(1)ne
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 


   write(*,*)'leyendo ',blckname,hwm-8
   read(1) hm1
   read (1)blckname,hwm
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   allocate(nh(nall(0)))
   read(1) hm1
   read(1)nh
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 

   write(*,*)'leyendo ',blckname,hwm-8
   read(1) hm1
   read (1)blckname,hwm
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   allocate(hsml(nall(0)))
   read(1) hm1
   read(1) hsml
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   
   write(*,*)'leyendo ',blckname,hwm-8
   read(1) hm1
   read (1)blckname,hwm
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
   allocate(sfr(nall(0)))
   read(1) hm1
   read(1) sfr
   read(1) hm2
   if (hm1/=hm2) stop 'apalapapa' 
close(1)  
!   write(*,*)'leyendo ',blckname,hwm-8
!   read(1) hm1
!   read (1)blckname,hwm
!   read(1) hm2
!   if (hm1/=hm2) stop 'apalapapa' 
!   allocate(abvc(nall(0)))
!   read(1) hm1
!   read(1) abvc
!   read(1) hm2
!   if (hm1/=hm2) stop 'apalapapa' 
endsubroutine
subroutine linkedlist(n,abin,cell,pos,head,tot,link)
        implicit none
        integer :: i,n,bx,by,bz,cell
        integer,dimension(cell,cell,cell) :: tot, head
        real,dimension(3,n) :: pos
        integer,dimension(n) :: link
        real :: abin
        
        do i = 1,n
               bx = int(pos(1,i)/abin) + 1
               by = int(pos(2,i)/abin) + 1
               bz = int(pos(3,i)/abin) + 1
               if (bx==0)     bx=cell
               if (bx==cell+1)bx=1
               if (by==0)     by=cell
               if (by==cell+1)by=1
               if (bz==0)     bz=cell
               if (bz==cell+1)bz=1

               if (bx > cell) write(*,*) pos(1,i),i
               tot(bx,by,bz) = tot(bx,by,bz) + 1
               head(bx,by,bz)   = i
       enddo

       do i = 1,n
               bx = int(pos(1,i)/abin) + 1
               by = int(pos(2,i)/abin) + 1
               bz = int(pos(3,i)/abin) + 1

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


