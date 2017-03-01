!............................................................
!...                      Subroutine energy               ...
!............................................................

!
! In order to have actualized arrays. It is very convient
! that before to call this subroutine calls the POTEN
! subroutine (for delj4)  and derden (for derivadores of 
! the helium4 density)subroutines


subroutine energy()

use alphasterm
use Para_DerivnD
use deriva
use energies
use field
! use impur
use classicimp
use lenard4
use grid
use gridk
use he4
use rho
use util1
use work1

implicit none

real (kind=8)    ::  aux1,aux2,aux3,aux4,aux5 ! Auxiliar variables
complex (kind=8) :: invars(6)
integer (kind=4) :: ix,iy,iz

! include 'interface_derivnD.include'  ! Per fer servir les derivades generiques

!................................ Use derivatives for (grad(den))**2
!
!  icon   =  0 ! Take the derivative.
!  icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!  icon   =  8 ! Take the derivative. Use periodic conditions.

!call pdergc(1,npd,nn,hx,psi,sto1c,3,1,mmx,iw,icon)   ! derivative respect X
!call pdergc(1,npd,nn,hy,psi,sto2c,3,2,mmx,iw,icon)   ! derivative respect Y
!call pdergc(1,npd,nn,hz,psi,sto3c,3,3,mmx,iw,icon)   ! derivative respect Z

!  Call deriv3Dc(1,nn,hx,1,psi,sto1c,Icon)
!  Call deriv3Dc(1,nn,hy,2,psi,sto2c,Icon)
!  Call deriv3Dc(1,nn,hz,3,psi,sto3c,Icon)

  Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(1,nn,hz,3,psi,sto3c,Icon)

aux1   = 0.5d0*cp4
aux2   = (1.0d0/3.0d0)*cpp4

!..................................................................
!... Calculate the density of energy (kinetic, and correlation) ...
!..................................................................

ekin4 = 0.0d0
elj4  = 0.0d0
ecor4 = 0.0d0

do iz=1,nz
   do iy=1,ny
      do ix=1,nx
          aux3  = den(ix,iy,iz)
          aux4  = dencg(ix,iy,iz)
          aux5  = aux3*aux4**2*(aux1+aux2*aux4)
          ekin4 = ekin4 + abs(sto1c(ix,iy,iz))**2+abs(sto2c(ix,iy,iz))**2+abs(sto3c(ix,iy,iz))**2
          elj4  = elj4  + delj4(ix,iy,iz)*aux3   
          ecor4 = ecor4 + aux5
      end do
   end do
end do

! if(limp) then
!    call pdergc(1,npd,nn,hx,psix,sto1c,3,1,mmx,iw,icon)   ! derivative respect X
!    call pdergc(1,npd,nn,hy,psix,sto2c,3,2,mmx,iw,icon)   ! derivative respect Y
!    call pdergc(1,npd,nn,hz,psix,sto3c,3,3,mmx,iw,icon)   ! derivative respect Z
!    eimpu = 0.0d0             
!    ekinx = 0.0d0
!    do iz=1,nz
!       do iy=1,ny
!          do ix=1,nx
!           ekinx = ekinx + abs(sto1c(ix,iy,iz))**2+abs(sto2c(ix,iy,iz))**2+abs(sto3c(ix,iy,iz))**2
!           eimpu = eimpu + potx4(ix,iy,iz)*den(ix,iy,iz)
!          end do
!       end do
!    end do
! end if

!......................................................
!... Calculate the density of energy (alpha_s term) ...
!......................................................

select case(core4)
   case('OTE','OTC')
     ealphas = sum(falfs*( dxden*intxalf + dyden*intyalf + dzden*intzalf) )
     ealphas = -h2o2m4*0.5d0*alphas*ealphas*dxyz
   case default
     continue
end select



esolid=0.0d0
if(lsolid)esolid = C*sum(den*(1.d0+dtanh(beta*(den-den_m))))*dxyz

ekin4 = h2o2m4*ekin4*dxyz           ! TOTAL Kinetic energy for 4He
elj4  = 0.5d0 *elj4 *dxyz           ! TOTAL Lennard-Jones energy
ecor4 =        ecor4*dxyz           ! TOTAL Correlation energy for 4He
etot4 = ekin4+elj4+ecor4 + esolid   ! TOTAL ENERGY without impurity

select case(core4)
   case('OTE','OTC')
     etot4 = etot4+ealphas          ! TOTAL ENERGY including Alpha_s term
   case default
     continue
end select

etot4    =  etot4 + sum(uext*den)*dxyz

! Classic vecotrial particle energy:
 ekinx = 0.5d0*mAg_u*sum(vimp*vimp)
 eHeX = sum(uimp*den)*dxyz

 eimpu = ekinx + eHeX 

etot   =  etot4 + eimpu

! if(limp) then                       ! ..... Term due to impurities
!   ekinx = h2o2mx*ekinx*dxyz         ! TOTAL Kinetic energy for the impurity
!   eimpu =        eimpu*dxyz         ! TOTAL energy due to the impurity
!   etot4 = etot4+eimpu               ! TOTAL ENERGY including impurity
!   etot  = etot4+ekinx
! end if


return

end
