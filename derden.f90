!............................................................
!...                      Subroutine derden               ...
!............................................................

subroutine derden()

use Para_DerivnD     ! (npd,dxden,dyden,dzden,icon)
use deriva     ! (npd,dxden,dyden,dzden,icon)
use grid       ! (nx,ny,nz,hx,hy,hz,nxyz,dxyz)
use rho        ! (den)
use util1      ! (nn,mmx,iw)


implicit none

!include 'interface_derivnD.include'  ! Per fer servir les derivades generiques

!call pderg(1,npd,nn,hx,den,dxden,3,1,mmx,iw,icon)   ! derivative respect X
!call pderg(1,npd,nn,hy,den,dyden,3,2,mmx,iw,icon)   ! derivative respect Y
!call pderg(1,npd,nn,hz,den,dzden,3,3,mmx,iw,icon)   ! derivative respect Z

!  Call deriv3D(1,nn,hx,1,den,dxden,Icon)
!  Call deriv3D(1,nn,hy,2,den,dyden,Icon)
!  Call deriv3D(1,nn,hz,3,den,dzden,Icon)

  Call derivnD(1,nn,hx,1,den,dxden,Icon)
  Call derivnD(1,nn,hy,2,den,dyden,Icon)
  Call derivnD(1,nn,hz,3,den,dzden,Icon)

return
end
