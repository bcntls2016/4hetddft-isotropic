subroutine potenimpini
use interpol !, only:DelInter,potpi,potdel,npot,vpi,delta
use grid
use classicimp!, only: rimp, lselection, N_imp
implicit none
integer (kind=4) :: i,k
real    (kind=8) :: r,rmax
real    (kind=8) :: r1,r2,Select_pot


npot =100*max(nx,ny,nz)+1
rmax = 4.d0*max(xmax,ymax,zmax)

DelInter = rmax/dfloat(npot-1)

!.....................................
! Interpolation for V_gs
!.....................................

allocate(potion(N_imp,npot))
do k=1,N_imp
  do i=1,npot
    r = dfloat(i-1)*DelInter
	potion(k,i) = Select_pot(selec_gs_k(k),r,r_cutoff_gs_k(k),umax_gs_k(k)) 
  enddo
enddo

call updatepoten()

end subroutine potenimpini

