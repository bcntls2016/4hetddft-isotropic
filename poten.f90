!---------------------------------------------------------------------
!---                             Subroutine Poten                  ---
!---------------------------------------------------------------------

subroutine poten()

use alphasterm ! (ualphas)
use grid     ! (nxyz)
use he4      ! (cp4,cpp4)
! use impur    ! (vq,)
use classicimp
use lenard4  ! (wk2,pelj4,fvlj4,delj4,core4,lalphas)
use field    ! (pot4,limp)
use rho      ! (dencg,den,wcgk,)
use work1    ! (sto1,sto2,sto3,sto4)

implicit none

real    (kind=8) :: a0,a1
integer (kind=4) :: ix,iy,iz

!...............................................
!... Calculate of Fourier Transforms for 4He ...
!...............................................

! call fftfw(den,fden)
 call fftfw_den()

!..............................
!.. Coarse-graining density ...
!.. Alfa_s density          ...
!....................................

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
   wk1(ix,iy,iz) = fden(ix,iy,iz)*wcgk(ix,iy,iz)
end forall
! call fftbk(wk1,dencg)  ! get Coarse graining density
call fftbk_cg()  ! get Coarse graining density

call derden() ! Calculate derivatives of the density

!..............................
! Lennard-Jones contribution  .
!..............................

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
   wk1(ix,iy,iz) = fden(ix,iy,iz)*fvlj4(ix,iy,iz)
end forall 
! call fftbk(wk1,delj4) ! Get delj4 -> (   int{ rho_4*V_4 dr'}  )
call fftbk_lj() ! Get delj4 -> (   int{ rho_4*V_4 dr'}  )

!........................
! Correlation terms.  ...
!........................

a0 = cp4 /2.d0               ! Auxiliar variable useful for saving operations
a1 = cpp4/3.d0               ! Auxiliar variable useful for saving operations

!   The rest of the correlation contribution is calculated as
!   a convolution product.

forall(ix=1:nx,iy=1:ny,iz=1:nz)
   sto1(ix,iy,iz) = den(ix,iy,iz)*dencg(ix,iy,iz)*(cp4+cpp4*dencg(ix,iy,iz))
end forall 

! call fftfw(sto1,wk1)
call fftfw_1()

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
   wk1(ix,iy,iz) = wk1(ix,iy,iz)*wcgk(ix,iy,iz)
end forall 

! call fftbk(wk1,sto1)
call fftbk_1()


!..........................
!... Solid penalty term ...
!..........................
If(lsolid)Then
  penalty = dtanh(beta*(den-den_m))
  penalty = C*(1.d0 + penalty + beta*den*(1.d0 - penalty**2) )
Endif

!.......................
!.. Final 'Potential' ...
!........................


forall(ix=1:nx,iy=1:ny,iz=1:nz)
   pot4(ix,iy,iz) = delj4(ix,iy,iz) +                                 &   ! Lennard-Jones
                    dencg(ix,iy,iz)**2*(a0+a1*dencg(ix,iy,iz)) +      &   ! Correlation
                    sto1(ix,iy,iz)                                        ! Correlation
end forall

if(core4.eq.'OTC') then
    forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
       wk1(ix,iy,iz)   = fden(ix,iy,iz)*kalfs(ix,iy,iz)
    end forall
!     call fftbk(wk1,denalf) ! Get Alfa_s density
    call fftbk_as() ! Get Alfa_s density
    call term_alfa()    ! Calculates alpha_s contribution to the field
    forall(ix=1:nx,iy=1:ny,iz=1:nz)
       pot4(ix,iy,iz) = pot4(ix,iy,iz)  +                                 &   ! Mean field...
                        ualphas(ix,iy,iz)                                     ! Alfa_s term
    end forall
end if


If(lsolid)Then
  forall(ix=1:nx,iy=1:ny,iz=1:nz)
    pot4(ix,iy,iz) = pot4(ix,iy,iz) + penalty(ix,iy,iz)
  end forall 
Endif

return

end
