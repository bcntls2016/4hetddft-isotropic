!............................................................
!...                      Subroutine term_alfa            ...
!............................................................

subroutine term_alfa()

use alphasterm ! (
               !  denalf,falfs,kalfs,intxalf,intyalf,intzalf)
use Para_DerivnD ! (npd,dxden,dyden,dzden,icon)
use deriva     ! (npd,dxden,dyden,dzden,icon)
use grid       ! (nx,ny,nz,hx,hy,hz,nxyz,dxyz)
use lenard4    ! (core4)
use he4        ! (h2o2m4,den0s)
use rho        ! (den)
use util1      ! (nn,mmx,iw)
use work1      ! (wk1,wk2,wk3,workx,worky,workz,sto1....sto7)

implicit none
integer  (kind=4) :: ix,iy,iz

!include 'interface_derivnD.include'  ! Per fer servir les derivades generiques

!.............................................................
!... Calculate  f(r), integral common to U_alfs,and E_alfs ...
!.............................................................


forall(ix=1:nx,iy=1:ny,iz=1:nz)
  falfs(ix,iy,iz) = (1.d0-denalf(ix,iy,iz)/den0s)   ! Obtain f(r) function.
end forall

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  sto1(ix,iy,iz) = falfs(ix,iy,iz)*dxden(ix,iy,iz)  ! Obtain 'f(r)*grad(rho)'    (x-component)
  sto2(ix,iy,iz) = falfs(ix,iy,iz)*dyden(ix,iy,iz)  !   "                        (y-component)
  sto3(ix,iy,iz) = falfs(ix,iy,iz)*dzden(ix,iy,iz)  !   "                        (z-component)
end forall

! call fftfw(sto1,wk1) ! Obtain FFT[f(r)*grad(rho)] (x-component)
! call fftfw(sto2,wk2) !   "                        (y-component)
! call fftfw(sto3,wk3) !   "                        (z-component)
call fftfw_123() !   Obtain FFT[f(r)*grad(rho)] (all)

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
  wk1(ix,iy,iz) = kalfs(ix,iy,iz)*wk1(ix,iy,iz)  ! Calculate FFT [ F(r,r') ] * FFT[ f(r)*grad(rho) ] (x-component)
  wk2(ix,iy,iz) = kalfs(ix,iy,iz)*wk2(ix,iy,iz)  !   "                                               (y-component)
  wk3(ix,iy,iz) = kalfs(ix,iy,iz)*wk3(ix,iy,iz)  !   "                                               (z-component)
end forall

! call fftbk(wk1,intxalf)  ! Obtain FFT^-1(exp(-(pi*l*p)**2 * FFT(....) )  x-component
! call fftbk(wk2,intyalf)  !  "                                            y-component
! call fftbk(wk3,intzalf)  !  "                                            z-component
call fftbk_xyz() ! Obtain FFT^-1(exp(-(pi*l*p)**2 * FFT(....) ) (all)

if(core4.eq.'OTC') then 


!.....................................
!.. Contribution to the mean field ...
!.....................................

!
!... Term with the double integral.
!

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  sto1(ix,iy,iz) = (dxden(ix,iy,iz)*intxalf(ix,iy,iz)) +    &
                   (dyden(ix,iy,iz)*intyalf(ix,iy,iz)) +    &
                   (dzden(ix,iy,iz)*intzalf(ix,iy,iz))
end forall

! call fftfw(sto1,wk1)         ! Obtain FFT[f(r)*grad(rho)]
call fftfw_1()         ! Obtain FFT[f(r)*grad(rho)]

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
  wk1(ix,iy,iz) = kalfs(ix,iy,iz)*wk1(ix,iy,iz)
end forall

! call fftbk(wk1,ualphas)      ! Obtain first term of the U_alphas mean field.
!                              ! actually is ualphas*den0s
call fftbk_ua()      ! Obtain first term of the U_alphas mean field.
                             ! actually is ualphas*den0s

!....................................................
!... Term with the gradient  and a single integral
!....................................................

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  ualphas(ix,iy,iz) = ualphas(ix,iy,iz)/den0s
  sto1(ix,iy,iz)    = falfs(ix,iy,iz)*intxalf(ix,iy,iz)
  sto2(ix,iy,iz)    = falfs(ix,iy,iz)*intyalf(ix,iy,iz)
  sto3(ix,iy,iz)    = falfs(ix,iy,iz)*intzalf(ix,iy,iz)
end forall

!call pderg(1,npd,nn,hx,sto1,sto4,3,1,mmx,iw,icon)   ! derivative respect X
!call pderg(1,npd,nn,hy,sto2,sto5,3,2,mmx,iw,icon)   ! derivative respect Y
!call pderg(1,npd,nn,hz,sto3,sto6,3,3,mmx,iw,icon)   ! derivative respect Z

!Call deriv3D_p(1,nn,hx,1,Sto1,Sto4,Icon)
!Call deriv3D_p(1,nn,hy,2,Sto2,Sto5,Icon)
!Call deriv3D_p(1,nn,hz,3,Sto3,Sto6,Icon)

!Call deriv3D(1,nn,hx,1,Sto1,Sto4,Icon)
!Call deriv3D(1,nn,hy,2,Sto2,Sto5,Icon)
!Call deriv3D(1,nn,hz,3,Sto3,Sto6,Icon)

Call derivnD(1,nn,hx,1,Sto1,Sto4,Icon)
Call derivnD(1,nn,hy,2,Sto2,Sto5,Icon)
Call derivnD(1,nn,hz,3,Sto3,Sto6,Icon)

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  ualphas(ix,iy,iz) = (ualphas(ix,iy,iz)  + sto4(ix,iy,iz)+sto5(ix,iy,iz)+sto6(ix,iy,iz)) &
                      *h2o2m4*alphas
end forall

return
end if
end
