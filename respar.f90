!........................................................................
!...                       Subroutine respar                          ...
!........................................................................

subroutine respar(x,y,z,nx,ny,nz,nfun,name1,name2,fun1,fun2)
implicit none
integer (kind=4)  :: nfun
integer (kind=4)  :: nx,ny,nz
integer (kind=4)  :: ix,iy,iz
real    (kind=8)  :: x(nx),y(ny),z(nz)
real    (kind=8)  :: fun1(nx,ny,nz)
real    (kind=8)  :: fun2(nx,ny,nz)
real    (kind=8)  :: aux
character (len=*) :: name1,name2

open(10,file=name1//'-x.dat')
do ix=1,nx
  aux=fun1(ix,ny/2+1,nz/2+1)
  if(abs(aux).lt.1.d-40) aux=0.0d0
  write(10,7020) x(ix),aux
end do
close(10)
open(10,file=name1//'-y.dat')
do iy=1,ny
  aux=fun1(nx/2+1,iy,nz/2+1)
  if(abs(aux).lt.1.d-40) aux=0.0d0
  write(10,7020) y(iy),aux
end do
close(10)
open(10,file=name1//'-z.dat')
do iz=1,nz
  aux=fun1(nx/2+1,ny/2+1,iz)
  if(abs(aux).lt.1.d-40) aux=0.0d0
  write(10,7020) z(iz),aux
end do
close(10)
if(nfun.eq.2) then
   open(11,file=name2//'-x.dat')
   do ix=1,nx
     aux=fun2(ix,ny/2+1,nz/2+1)
     if(abs(aux).lt.1.d-40) aux=0.0d0
     write(11,7020) x(ix),aux
   end do
   close(11)
   open(11,file=name2//'-y.dat')
   do iy=1,ny
     aux=fun2(nx/2+1,iy,nz/2+1)
     if(abs(aux).lt.1.d-40) aux=0.0d0
     write(11,7020) y(iy),aux
   end do
   close(11)
   open(11,file=name2//'-z.dat')
   do iz=1,nz
     aux=fun2(nx/2+1,ny/2+1,iz)
     if(abs(aux).lt.1.d-40) aux=0.0d0
     write(11,7020) z(iz),aux
   end do
   close(11)
end if
return
7020 format(1x,1P,2E14.5)
end
