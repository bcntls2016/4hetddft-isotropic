subroutine potenimp(rimp)
!
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
!
use classicimp , only: uimp!,Also2!,pairpot
use grid
use interpol !, only: potdel,potpi,DelInter
implicit none
integer (kind=4) :: ix,iy,iz
integer (kind=4) :: ir
real    (kind=8) :: zt,yt,r
real    (kind=8) :: xx,yy,zz
real    (kind=8) , intent(in) :: rimp(3)
real    (kind=8)              :: rmod

call updatepoten(rimp)

end subroutine potenimp

!---------------------------------------------------------------------------!

subroutine forceimp(rimp,F)
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
use classicimp , only: uimp
use deriva
use grid
use rho
implicit none
integer (kind=4) :: ix,iy,iz
real    (kind=8) :: zt,yt,zt2,yt2,r,d
real    (kind=8) :: aux1,aux2,aux3,drV_HeXor,dzV_XTiO,dV
real    (kind=8) , intent(in) :: rimp(3)
real    (kind=8) , intent(out):: F(3)

!.................................!
!... First, the term due to He ...!
!.................................!


F(1) = -sum(dxden*uimp)*dxyz
F(2) = -sum(dyden*uimp)*dxyz
F(3) = -sum(dzden*uimp)*dxyz

!...........................................!
!... Second, the term due to the surface ...!
!...........................................!

!F(3) = F(3) - dzV_XTiO(rimp(3)) 


end subroutine forceimp


!double precision function V_ion(x)
!use impur, only :r_cutoff,selec,umax
!implicit none
!Real (Kind=8)  :: r_cutoff=2.0d0, umax=7476.405d0 
!Character  (Len=80) :: selec='Rb_plus_Fausto'
!real (kind=8) :: x, Select_Pot
!V_ion = Select_Pot(selec,x,r_cutoff,umax)
!   V_ion = min(15000.d0,V_ion)

!end function



subroutine updatepoten(rimp)
use classicimp , only: uimp!,Also2!,pairpot
use grid
use interpol
implicit none
real    (kind=8) , intent(in) :: rimp(3)
real    (kind=8)              :: dist(3)
real    (kind=8)              :: r,yt,zt,rmod
integer (kind=8)              :: ix,iy,iz,ir


 do iz=1,nz
  zt = (z(iz)-rimp(3))**2
  do iy=1,ny
   yt = (y(iy)-rimp(2))**2 + zt
   do ix=1,nx
    r = dsqrt((x(ix)-rimp(1))**2 + yt)
    ir = int(r/DelInter)+1
    rmod = mod(r,DelInter)/DelInter
      uimp(ix,iy,iz) =  potion(ir)*(1.d0-rmod) +  potion(ir+1)*rmod
    enddo
   enddo
  enddo

end subroutine updatepoten
