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
integer (kind=4) :: ix,iy,iz,k,m
real (kind=8) :: r_ij(3), Select_pot

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
 ekinx = 0d0
 do k=1, N_imp
   ekinx = ekinx + 0.5d0*m_imp(k)*sum(vimp(k,:)*vimp(k,:))
 enddo
 eHeX = sum(uimp*den)*dxyz

 eimpu_impu = 0
 do k=1,N_imp
   do m=k+1,N_imp
     r_ij(:) = rimp(k,:)-rimp(m,:)
	 aux1 = dsqrt(sum(r_ij(:)**2))
	 eimpu_impu = eimpu_impu + Select_pot(selec_gs_k_k(k,m),aux1,r_cutoff_gs_k_k(k,m),umax_gs_k_k(k,m))
   enddo
 enddo

 eimpu = ekinx + eHeX + eimpu_impu
 etot   =  etot4 + eimpu

return

end
