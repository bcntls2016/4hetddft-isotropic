!.....................................................................
!...                      Subroutine initcg                        ...
!.....................................................................

subroutine initcg(h,wcgk)

! This subrotuine computes the Analitical Fourier Transform of
! the kernel used for coarse graining density.
!
! INPUT:
!   List:
!        h -------> value for the kernel
!   Modules:
!        nx,ny,nz -> Number of points          (in grid)
!        px,py,pz -> Grid in k-space           (in gridk)
!        pmod     -> Value of the module of p  (in grid)
!        twopi  ---> 2*pi                      (in util1)
! OUTPUT:
!        wcgk -----> FFT(discrete) of the kernel of the Coarse-Graining

use grid
use gridk
use util1
 
implicit none

integer    (kind=4)  :: ipx,ipy,ipz
real       (kind=8)  :: aux
real       (kind=8)  :: p,arg,h
real       (kind=8)  :: wcgk(nx/2+1,ny,nz)

do ipz=1,nz
  do ipy=1,ny
    do ipx=1,nx/2+1
       p    = pmod(ipx,ipy,ipz)
       if(p.gt.0.0d0) then
         arg = twopi*p*h
         aux = sin(arg)-arg*cos(arg)
         if(abs(aux).lt.1.d-10) then
             wcgk(ipx,ipy,ipz) = 1.0d0
         else
             wcgk(ipx,ipy,ipz) = 3.0d0*aux/(arg**3)
         end if
       else
         wcgk(ipx,ipy,ipz) = 1.0d0
       end if
    end do
  end do
end do
return
end
