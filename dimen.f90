!..............................................................
!...                     Subroutine dimen                   ...
!..............................................................
!
! This subroutine allocates almost all the arrays....
!
!
subroutine dimen()

use deriva
use field
use grid
use gridk
use he4
use classicimp
! use impur
use lenard4
use rho
use util1
use work1
use rkpc

implicit none

!.........................................
!.. Arrays for real and momentum grids ...
!.........................................
!
allocate (x(nx)) ; allocate (px(nx)) ;  
allocate (y(ny)) ; allocate (py(ny)) ;
allocate (z(nz)) ; allocate (pz(nz)) ;
allocate ( pmod(nx/2+1,ny,nz))    ! Array with the modules of p.

!............................................
!.. Arrays for Lennard-Jones calculations ...
!............................................
!
allocate (fvlj4(nx/2+1,ny,nz))    ! Array with the FFT of the Kernel of LJ potential  
allocate (delj4(nx    ,ny,nz))    ! Array with rhe energy-density of LJ

!..................................................
!.. Arrays for partial derivatives of densities ...
!..................................................
!
allocate(dxden(nx,ny,nz))       ! Array for partial derivatives.
allocate(dyden(nx,ny,nz))       !        "
allocate(dzden(nx,ny,nz))       !        "

!..........................
!.. Arrays for Helium 4 ...
!..........................
!
allocate (     pot4(nx,ny,nz))    ! Array with al the 'potential'
allocate (     hpsi(nx,ny,nz))    ! Array with H PSi
allocate (  den(nx    ,ny,nz))    ! Array with the density
allocate (  psi(nx    ,ny,nz))    ! Array with the field
allocate ( psiold(nx   ,ny,nz,3))  ! Array with the field 1-step old iteration
allocate ( hpsiold(nx  ,ny,nz,2))  ! Array with the field 1-step old iteration
allocate (dencg(nx    ,ny,nz))    ! Array with the coarse-graining density.
allocate ( fden(nx/2+1,ny,nz))    ! Array with the density in p-space
!allocate ( fpsi(nx/2+1,ny,nz))    ! Array with the FFT-field
allocate ( wcgk(nx/2+1,ny,nz))    ! Array with the kernel of the coarse-graining

!..............................
!.. Arrays for the impurity ...
!..............................
!
   allocate( uext(nx,ny,nz))      ! Fourier transform of the impurity external potential
   If(Lsolid)allocate(penalty(nx,ny,nz))      ! Potential of penalty term, for the solid functional
   allocate( uimp(nx,ny,nz))      ! Fourier transform of the impurity external potential
!  allocate(pairpot(nx,ny,nz,6))      ! Fourier transform of the impurity external potential
!    allocate( uextimp(nx,ny,nz))      ! Fourier transform of the impurity external potential

! if(limp) then                     ! Things for impurities.
!    allocate(hpsix(nx,ny,nz))      ! Array with H PSi_x
!    allocate(upotx(nx,ny,nz))      ! Mean field for the impurity
!    allocate(potx4(nx,ny,nz))      ! 
!    allocate( psix(nx,ny,nz))      ! Wave function for the impurity
!    allocate(psixold(nx,ny,nz,3))  ! Wave function for the impurity
!    allocate(hpsixold(nx,ny,nz,2))  ! Wave function for the impurity
!    allocate( denx(nx,ny,nz))      ! density function for the impurity
! !   allocate(fpsix(nx/2+1,ny,nz))  ! Fourier Transform of the Wave function for the impurity
!    allocate(fdenx(nx/2+1,ny,nz))  ! Fourier Transform of the Wave function for the impurity
!    allocate(vq(nx/2+1,ny,nz))     ! 
!    allocate(Epsix(nx,ny,nz))      ! Envelop function
!    allocate(rmod0(nx,ny,nz))      ! For the  envelop function
! end if


!....................................................
!.. Arrays for temporal storage and working areas ...
!....................................................
!
allocate(sto1(nx,ny,nz))       ! Real*8 Array for temporal calculations...
allocate(sto2(nx,ny,nz))       !    "
allocate(sto3(nx,ny,nz))       !    "
allocate(sto4(nx,ny,nz))       !    "
allocate(sto5(nx,ny,nz))       !    "
allocate(sto6(nx,ny,nz))       !    "
allocate(wk1(nx/2+1,ny,nz))    ! Complex array for FFT temporal calculations
allocate(wk2(nx/2+1,ny,nz))    !    "
allocate(wk3(nx/2+1,ny,nz))    !    "
allocate(sto1c(nx,ny,nz))      ! Complex*8 Array for temporal calculations...
allocate(sto2c(nx,ny,nz))      !    "
allocate(sto3c(nx,ny,nz))      !    "
allocate(sto4c(nx,ny,nz))      !    "
allocate(sto5c(nx,ny,nz))      !    "
allocate(sto6c(nx,ny,nz))      !    "
allocate(sto7c(nx,ny,nz))      !    "
allocate(sto8c(nx,ny,nz))      !    "
!...................................................................
!.. Array for Runge-Kutta-Gill & Predictor-Corrector-Modificator ...
!...................................................................
!
allocate(q(nx,ny,nz))
allocate(pc(nx,ny,nz))         
q  = 0.d0
pc = 0.d0
! if(limp) then                     ! Things for impurities.
! allocate(qx(nx,ny,nz))
! allocate(pcx(nx,ny,nz))
! qx  = 0.d0
! pcx = 0.d0
! end if


allocate(timec(nx,ny,nz))
return

end
