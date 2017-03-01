subroutine potenimpini
use interpol !, only:DelInter,potpi,potdel,npot,vpi,delta
use grid
use classicimp, only: rimp, lselection
implicit none
integer (kind=4) :: i
real    (kind=8) :: r,rmax,V_gs
real    (kind=8) :: r1,r2


npot =100*max(nx,ny,nz)+1
rmax = 4.d0*max(xmax,ymax,zmax)

DelInter = rmax/dfloat(npot-1)

!.....................................
! Interpolation for V_gs
!.....................................

allocate(potion(npot))

do i=1,npot
  r = dfloat(i-1)*DelInter
  potion(i) = V_gs(r)
enddo

call updatepoten(rimp)

end subroutine potenimpini

