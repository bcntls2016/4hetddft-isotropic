!------------------------------------------------------------------
!---                    Subroutine readenc                      ---
!------------------------------------------------------------------

subroutine readenc(npart,densat,fileden,fileimp,mode,rimpur,r_clust)

! INPUT QUANTITIES:
!
! From list:
!
! npart   ----> Number of particles
! densat  ----> Density of saturation (useful for fermi distributions)
! fileden ----> Name of the file with the density
! mode    ----> Select the density:
!                       0 = continue a calculation
!                       1 = Build from scratch a new density
!                       2 = Build from scratch a new density with impurity
! rimpur  ----> Radius of the impurity
!
!
! OUTPUT QUANTITIES:
!
!  psi    ----> He Wave function)
!  den    ----> Density (module rho)
!  r_clust ---> Size of the cluster
!  psix ------> Wave function for the impurity
!
!
! NOTE:
! This subroutine check the consistency of namelist input data when an 
! external density is used as a input. The quantities to check consistency
! came from modules
!-------------------------------------------------------------------------------

use rho
use field
use grid
! use impur
use classicimp
use util1

implicit none
integer   (kind=4), intent(in)    :: npart
real      (kind=8)                :: nph     ! Num. part in hole
real      (kind=8), intent(in)    :: densat
character (len=60), intent(in)    :: fileden,fileimp
integer   (kind=4), intent(in)    :: mode
integer	(kind=4)					:: ninvar
complex (kind=8), 	allocatable	::	invar(:)
real      (kind=8), intent(out)   :: r_clust
real      (kind=8), intent(in)    :: rimpur

real      (kind=8) :: xmaxp,ymaxp,zmaxp,xcp,ycp,zcp,hxp,hyp,hzp
real      (kind=8) :: aux1, aux2, aux3
real      (kind=8) :: aux1b,aux2b,aux3b
! real      (kind=8) :: ximp,yimp,zimp

integer   (kind=4) :: nxp,nyp,nzp
integer   (kind=4) :: ix,iy,iz,isalto
!logical            :: limp,Ldensity=.true.
logical            :: limp
real      (kind=8) :: rr

real      (kind=8) :: lambda,norm

! real      (kind=8)    :: Oneosq3,Oneosq2! = dsqrt(1.d0/3.d0)
! real      (kind=8)    :: tre,fiv,norm,rr
! real      (kind=8) , parameter :: p1=0.573585,p2=0.584807
! parameter(Oneosq3 = dsqrt(1.d0/3.d0))
! parameter(Oneosq2 = dsqrt(1.d0/2.d0))
! parameter(tre=0.297044d0)
! parameter(fiv=0.495074d0)
!...........................................
!With 'mode' select the kind of density...
!...........................................
!

r_clust=0.d0

select case(mode)
!-------------------------------------------------------------------
  case(0)   ! Continue a calculation (read a previous density or wave function (vortex))
!-------------------------------------------------------------------

     open(unit=1,file=fileden,status='old')
     Go to 20
10   Continue
       Ldensity=.false.
20   Continue     
     Rewind(1)
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
     If(Ldensity)Then
       read(1,*,Err=10) den
       Write(6,'("From Readenc: We have read a density")')
     Else
       read(1,*)psi
       Write(6,'("From Readenc: We have read a complex w.f.")')
       den=Abs(psi)**2
     Endif      
     close(1)

!---------------------------------------------------------------
 case(2)   ! Continue a calculation (read a previous wave functon)
!---------------------------------------------------------------
   open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp
     read(1,*) rimp
     read(1,*) vimp
     read(1,*) psi
     den=Conjg(psi)*psi
   close(1)

!---------------------------------------------------------------
 case(3)   ! Continue a calculation (read a previous wave functon)
!---------------------------------------------------------------
   open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp
     read(1,*) rimp
     read(1,*) vimp
     read(1,*) ninvar
     allocate(invar(ninvar))
     read(1,*) invar
     read(1,*) psi
     den=Conjg(psi)*psi
   close(1)
!
!-------------------------------------------------
   case default ! Non programmed situation
!-------------------------------------------------
!
      stop 'ERROR READEN: This mode is still not programmed'
end select


return
end
