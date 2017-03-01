      SUBROUTINE STEPPC(deltat,errHe,errimp)
!
!      Predictor-Modifier-Corrector method
!       (Ralston & Wilf Vol I, pag. 99)
!
use Para_DerivnD ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
! use impur  !
use classicimp
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc)

implicit none

real (kind=8) :: c112=112.d0/121.d0
real (kind=8) :: c9  =9.d0/121.d0
real (kind=8) :: c1o3=1.d0/3.d0
real (kind=8) :: c4o3=4.d0/3.d0
real (kind=8) :: c5o3=5.d0/3.d0

integer (kind=4) :: ix,iy,iz,iaux
real    (kind=8) :: deltat
real    (kind=8) :: errHe, errimp
real    (kind=8) :: auxr(3)
complex (kind=8) :: auxc(6)
complex (kind=8) :: aux1c,aux2c,aux3c,aux4c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)

! include 'interface_derivnD.include'  ! Per fer servir les derivades generiques

!call pdergc(2,npd,nn,hx,psi,sto1c,3,1,mmx,iw,icon)
!call pdergc(2,npd,nn,hy,psi,sto2c,3,2,mmx,iw,icon)
!call pdergc(2,npd,nn,hz,psi,sto3c,3,3,mmx,iw,icon)

!  Call deriv3Dc(2,nn,hx,1,psi,sto1c,Icon)
!  Call deriv3Dc(2,nn,hy,2,psi,sto2c,Icon)
!  Call deriv3Dc(2,nn,hz,3,psi,sto3c,Icon)

  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

!
!   We compute H·Psi
!
!
!      Predictor:
!
Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 - pot4*psi) - ci*uimp*psi
Sto1c = psiold(:,:,:,ioldp(3)) + c4o3*deltat*(2.d0*Sto4c-hpsiold(:,:,:,ioldh(1))           &
        + 2.d0*hpsiold(:,:,:,ioldh(2)))
!
!      Modificador:
!
psiold(:,:,:,ioldp(3)) = psi
psi = Sto1c - c112*pc
pc  = Sto1c
hpsiold(:,:,:,ioldh(2)) = Sto4c
den = Abs(psi)**2
!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
iaux=ioldh(2)
ioldh(2)=ioldh(1)
ioldh(1)=iaux

!................!
!... position ...!
!................!
! Predictor
 auxr(:) = rimpold(:,ioldr(3)) + c4o3*deltat*(2.d0*vimp(:) - vimpold(:,ioldv(1)) + 2.d0*vimpold(:,ioldv(2)))
! Modificador
 rimpold(:,ioldr(3)) = auxr(:) - c112*pcr(:)
 pcr = auxr

!..................!
!... velocities ...!
!..................!
! Predictor
 auxr(:) = vimpold(:,ioldv(3)) + c4o3*deltat*(2.d0*aimp(:) - aimpold(:,iolda(1)) + 2.d0*aimpold(:,iolda(2)))
! Modificador
 vimpold(:,ioldv(3)) = auxr(:) - c112*pcv(:)
 pcv = auxr



 aimpold(:,iolda(2)) = aimp(:)
  ! Reubicacion indices
 iaux=iolda(2)  ; iolda(2)=iolda(1)   ; iolda(1)=iaux


!........................................
! call poten()
    auxr =  rimpold(:,ioldr(3))
! if(newpoten)then
!     call superpoten(auxr,auxc,aimp,Hinvar)
! else
    call potenimp(auxr)
    call poten()
    call forceimp(auxr,aimp)
    aimp(:) = aimp(:)/mAg_u
! endif
!........................................

!call pdergc(2,npd,nn,hx,sto4c,sto1c,3,1,mmx,iw,icon)
!call pdergc(2,npd,nn,hy,sto4c,sto2c,3,2,mmx,iw,icon)
!call pdergc(2,npd,nn,hz,sto4c,sto3c,3,3,mmx,iw,icon)

!  Call deriv3Dc(2,nn,hx,1,psi,sto1c,Icon)
!  Call deriv3Dc(2,nn,hy,2,psi,sto2c,Icon)
!  Call deriv3Dc(2,nn,hz,3,psi,sto3c,Icon)

  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4  - pot4*psi) -ci*uimp*psi
Sto5c = 0.125d0*( 9.d0*psiold(:,:,:,ioldp(3)) - psiold(:,:,:,ioldp(2))   &
      +3.d0*deltat*(Sto4c + 2.d0*hpsiold(:,:,:,ioldh(1)) - hpsiold(:,:,:,ioldh(2))  ))
pc = pc - Sto5c
!
!     Valor final:
!
psi = Sto5c + c9*pc
errHe = c9*Sum(Abs(pc))
den = Abs(psi)**2

errHe=errHe/nxyz

!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
      iaux=ioldp(3)
      ioldp(3)=ioldp(2)
      ioldp(2)=ioldp(1)
      ioldp(1)=iaux

!.................!
!... positions ...!
!.................!
! Corrector:
auxr(:) = 0.125d0*( 9.d0*rimp(:) - rimpold(:,ioldr(2))     &
                 +3.d0*deltat*(vimpold(:,ioldv(3)) + 2.d0*vimp(:) - vimpold(:,ioldv(1)) ))
! vpold3 is actually the v_temporal just computed, so it is the 'newest'.
! The combination of v and vold is different because THE INDEXS HAVE NOT BEEN REALLOCATED YET.
pcr = pcr -auxr
rimpold(:,ioldr(3)) = rimp(:)
! Valor final:
  rimp = auxr + c9*pcr
errimp = sum(Abs(c9*pcr))*0.3333333333d0
! Reubicacion
iaux=ioldr(3) ; ioldr(3)=ioldr(2) ; ioldr(2)=ioldr(1) ; ioldr(1)=iaux

!..................!
!... velocities ...!
!..................!
! Corrector:
auxr(:) = 0.125d0*( 9.d0*vimp(:) - vimpold(:,ioldv(2))     &
                 +3.d0*deltat*(aimp(:) + 2.d0*aimpold(:,iolda(1)) - aimpold(:,iolda(2)) ))
pcv = pcv -auxr
vimpold(:,ioldv(3)) = vimp(:)
! Valor final:
  vimp = auxr + c9*pcv
! errvimp = sum(Abs(c9*pcv))*0.3333333333d0


! Reubicacion
iaux=ioldv(3) ; ioldv(3)=ioldv(2) ; ioldv(2)=ioldv(1) ; ioldv(1)=iaux


return
end
