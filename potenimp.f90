subroutine potenimp()
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
real    (kind=8) :: rmod

call updatepoten()

end subroutine potenimp

!---------------------------------------------------------------------------!

subroutine forceimp()
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
use classicimp!, only: uimp_k, rimp, F, F_ij, N_imp
use deriva
use grid
use rho
implicit none
integer (kind=4) :: ix,iy,iz,k,m
real    (kind=8) :: zt,yt,zt2,yt2,r,d,Select_pot
real    (kind=8) :: aux1,aux2,aux3,drV_HeXor,dzV_XTiO,dV

!.................................!
!... First, the term due to He ...!
!.................................!

F_ij = 0
do k=1,N_imp
  do m=k+1,N_imp
  aux1 = dsqrt((rimp(k,1)-rimp(m,1))**2 + (rimp(k,2)-rimp(m,2))**2 + (rimp(k,3)-rimp(m,3))**2)
  F_ij(k,m,:) = -Select_pot(drselec_gs_k_k(k,m),aux1,drr_cutoff_gs_k_k(k,m),drumax_gs_k_k(k,m)) * (rimp(k,:)-rimp(m,:))/aux1
  F_ij(m,k,:) = -F_ij(k,m,:)
  enddo
enddo

do k=1,N_imp
  F(k,1) = -sum(dxden(:,:,:)*uimp_k(k,:,:,:))*dxyz
  F(k,2) = -sum(dyden(:,:,:)*uimp_k(k,:,:,:))*dxyz
  F(k,3) = -sum(dzden(:,:,:)*uimp_k(k,:,:,:))*dxyz
  do m=1,N_imp
    F(k,1) = F(k,1) + F_ij(k,m,1)
    F(k,2) = F(k,2) + F_ij(k,m,2)
    F(k,3) = F(k,3) + F_ij(k,m,3)
  enddo
enddo

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



subroutine updatepoten()
use classicimp , only: uimp_k, uimp, rimp, N_imp!,Also2!,pairpot
use grid
use interpol
implicit none
real    (kind=8)              :: dist(3)
real    (kind=8)              :: r,yt,zt,rmod
integer (kind=4)              :: ix,iy,iz,ir,k,m


do k=1,N_imp
  do iz=1,nz
    zt = (z(iz)-rimp(k,3))**2
    do iy=1,ny
      yt = (y(iy)-rimp(k,2))**2 + zt
      do ix=1,nx
        r = dsqrt((x(ix)-rimp(k,1))**2 + yt)
        ir = int(r/DelInter)+1
        rmod = mod(r,DelInter)/DelInter
        uimp_k(k,ix,iy,iz) =  potion(k,ir)*(1.d0-rmod) +  potion(k,ir+1)*rmod
	  enddo
    enddo
  enddo
enddo

uimp = 0d0
do k=1,N_imp
	uimp(:,:,:) = uimp(:,:,:) + uimp_k(k,:,:,:)
enddo

end subroutine updatepoten
