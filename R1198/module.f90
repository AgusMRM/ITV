module modulos
   implicit none
   integer, parameter :: FILES = 1         ! number of files per snapshot
   character(len=200), parameter :: path='/mnt/is2/dpaz/ITV/R1198/out'
   !character(len=200), parameter :: path='/mnt/is2/fstasys/ITV/S1050/out'

   character(len=200) :: filename, fnumber

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
   real*4, allocatable,dimension(:)    :: mass, u, dens, ne, &
                                        nh,hsml,sfr,abvc,temp 
   real*8 n0,mol,R,vx,vy,vz
   real :: prom1,prom2,prom3
       
   character(len=4) :: blckname 
   integer(kind=4) :: hwm
   real(kind=4)  :: hwm2 
   real :: xmin,ymin,zmin
   real :: xmax,ymax,zmax
   integer(4) :: i1,i2,npgs

   logical :: ilogic
endmodule

subroutine reader(snumber)
        use modulos
       character(len=200) :: snumber

   write(fnumber,'(I3)') 1
   do i=1,3 
      if (snumber(i:i) .eq. ' ') then 
          snumber(i:i)='0'
      end if
   end do

   !filename= path(1:len_trim(path)) // '/snapshot_' // snumber(1:len_trim(snumber)) !!// '.' // fnumber(verify(fnumber,' '):3)
   filename= trim(path)//'/snapshot_'//trim(snumber) !!// '.' // fnumber(verify(fnumber,' '):3)

   ! now, read in the header

   open (1, file=filename, form='unformatted')
   !open (1, file=filename, access='stream')

   read (1)blckname,hwm

   read (1) npart, massarr, a, redshift, flag_sfr,flag_feedback, nall,flag_cool, num_files, boxsize,omega_m,omega_l,unused

 !  do i=0,5
 !     print *,'Type=',i,'    Particles=', nall(i)
 !  end do
   ! Now we just read in the coordinates, and only keep those of type=1
   ! The array PartPos(1:3,1:Ntot) will contain all these particles
   Ntot= sum(nall)

   allocate(pos(3,Ntot))
   allocate(vel(3,Ntot))
   N=sum(npart)

   
   !hwm contiene el numero de bytes del bloque de datos mas 8 bytes
   read (1)blckname,hwm
   read (1)pos
 !  write(*,*) blckname,hwm
 !  write(*,*)'leyendo ',blckname,hwm-8
   read (1)blckname,hwm
   read (1) vel
   read (1)blckname,hwm
 !  write(*,*) blckname,hwm
 !  write(*,*)'leyendo ',blckname,hwm-8
   allocate(id(int((hwm-8)/4)))
   read (1)id
        
 
   read (1)blckname,hwm
 !  write(*,*)'leyendo ',blckname,hwm-8
   allocate(idch(int((hwm-8)/4)))
   read (1)idch
   read (1)blckname,hwm
 !  write(*,*)'leyendo ',blckname,hwm-8
   allocate(idgn(int((hwm-8)/4)))
   read (1)idgn

   read (1)blckname,hwm
 !  write(*,*)'leyendo ',blckname,hwm-8
 !       print*, 'allocateo con ',int(real((hwm-8)/4.)-1)
   allocate(mass(int((hwm-8)/4)))
   read (1)mass
   read (1)blckname,hwm
 !  write(*,*)'leyendo ',blckname,hwm-8
   allocate(u(int((hwm-8)/4)))
   !write(*,*) 'allocate con', int(real(hwm-8)/4.)
   read (1)u
   read (1)blckname,hwm
 !  write(*,*)'leyendo ',blckname,hwm-8
   allocate(dens(int(hwm-8)/4))
   read (1)dens
 
   read (1)blckname,hwm
   allocate(ne(int(hwm-8)/4))
   read(1)ne
   close(1)
   npgs = nall(0)
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


