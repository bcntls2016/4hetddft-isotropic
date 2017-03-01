      SUBROUTINE STEPRK(deltat)
!
!      Runge-Kutta-Gill method
!       (Ralston & Wilf Vol I, pag. 117)
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
use rkpc   ! (Storage for Steprk & Steppc rutines)

implicit none

real (kind=8) :: arun(4),brun(4),crun(4)

integer (kind=4) :: ix,iy,iz,jrun
real    (kind=8) :: deltat
real    (kind=8) :: aux2r(3)
complex (kind=8) :: auxc(6)
complex (kind=8) :: aux1c,aux2c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)

! include 'interface_derivnD.include'  ! Per fer servir les derivades generiques

arun(1)=0.5d0
arun(2)=1.0d0-1.d0/dsqrt(2.d0)
arun(3)=1.0d0+1.d0/dsqrt(2.d0)
arun(4)=1.d0/6.d0

brun(1)=2.0d0
brun(2)=1.0d0
brun(3)=1.0d0
brun(4)=2.0d0

 crun(1)=0.5d0
 crun(2)=1.0d0-1.d0/dsqrt(2.d0)
 crun(3)=1.0d0+1.d0/dsqrt(2.d0)
 crun(4)=0.5d0


!........................................
!.. Laplacian of H^{iorder-1}ï¿½Psi (0) ...
!........................................
!
!   icon   =  0 ! Take the derivative.
!   icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!   icon   =  8 ! Take the derivative. Use periodic conditions.
!

do jrun=1,4

!  call pdergc(2,npd,nn,hx,psi,sto1c,3,1,mmx,iw,icon)
!  call pdergc(2,npd,nn,hy,psi,sto2c,3,2,mmx,iw,icon)
!  call pdergc(2,npd,nn,hz,psi,sto3c,3,3,mmx,iw,icon)

!  Call deriv3Dc(2,nn,hx,1,psi,sto1c,Icon)
!  Call deriv3Dc(2,nn,hy,2,psi,sto2c,Icon)
!  Call deriv3Dc(2,nn,hz,3,psi,sto3c,Icon)

  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

!
!   We compute H·Psi
!
Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 - pot4*psi) - ci*uimp*psi
Sto1c = arun(jrun)*(Sto4c - brun(jrun)*q)
q = q + 3.*Sto1c - crun(jrun)*Sto4c
if(jrun.eq.1)then
  hpsiold(:,:,:,2) = hpsiold(:,:,:,1)
  hpsiold(:,:,:,1) = Sto4c
  psiold(:,:,:,3) = psiold(:,:,:,2)
  psiold(:,:,:,2) = psiold(:,:,:,1)
  psiold(:,:,:,1) = psi
endif
psi = psi + deltat*Sto1c
den = Abs(psi)**2

!...............................................................test
! print*,'inside steprk, jrun=',jrun
 do iz=1,nz ; do iy=1,ny ; do ix=1,nx 
  if(.not.(den(ix,iy,iz).gt.0))then ; print*,ix,iy,iz,den(ix,iy,iz) 
  call respar(x,y,z,nx,ny,nz,1,'uimp','den',uimp,den); STOP ; endif
 end do ; enddo ; enddo
!...............................................................test

!
! Impurity evolution if it is necessary
!
!.................!
!... positions ...!
!.................!
  aux2r(:) = arun(jrun)*(vimp(:) - brun(jrun)*qr(:))
     qr(:) = qr(:) + 3.*aux2r(:) - crun(jrun)*vimp(:)
  if(jrun.eq.1)then
!    vpold(:,:,2) = vpold(:,:,1) 
!    vpold(:,:,1) =    vp(:,:)            
   rimpold(:,3) = rimpold(:,2) 
   rimpold(:,2) = rimpold(:,1) 
   rimpold(:,1) = rimp(:) 
  endif

  rimp(:) = rimp(:) + deltat*aux2r(:)

!..................!
!... velocities ...!
!..................!
  aux2r(:) = arun(jrun)*(aimp(:) - brun(jrun)*qv(:))
  qv(:) = qv(:) + 3.*aux2r(:) - crun(jrun)*aimp(:)
  if(jrun.eq.1)then
   aimpold(:,2) = aimpold(:,1) 
   aimpold(:,1) =    aimp(:)            
   vimpold(:,3) = vimpold(:,2) 
   vimpold(:,2) = vimpold(:,1) 
   vimpold(:,1) =    vimp(:) 
  endif
  vimp(:) = vimp(:) + deltat*aux2r(:)

! ...........................................................

  if(jrun.le.3)then
! if(newpoten)then
!     call superpoten(rimp,invar,aimp,Hinvar)
! else
    call potenimp(rimp)
    call poten()
    call forceimp(rimp,aimp)
    aimp(:) = aimp(:)/mAg_u
! endif
  endif
enddo

do ix=1,3
  ioldp(ix)=ix
  ioldr(ix)=ix
  ioldv(ix)=ix
enddo

do ix=1,2
  ioldh(ix)=ix
  iolda(ix)=ix
enddo

return
end
