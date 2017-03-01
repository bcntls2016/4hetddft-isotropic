!--------------------------------------------------------------------
!---                    Subroutine r_cm                         ---
!--------------------------------------------------------------------
!
! Gives the z-component of the center-of-mass for the rho4 density
!
! INPUT:
!
! Npart -> Number of particles
! Rho  --> Density or Psi**2
!
! Through Module grid: 
!    arrays of x(nx), y(ny) and z(nz)
!    dimensions nx,ny and nz
!    dxyz -> volume element
!
! OUTPUT:
!
! (xcm,ycm,zcm) ---> coordinates of the center of mass
!
!
!
subroutine r_cm(rho,npart,xcm,ycm,zcm)

use grid !(x,y,z,dxyz,nx,ny,nz)

implicit none
integer (kind=4)                      :: ix,iy,iz,npart
real    (kind=8)                      :: xcm,ycm,zcm
real    (kind=8), dimension(nx,ny,nz) :: rho

real    (kind=8) :: xx,yy,zz,aux

xcm=0.0d0
ycm=0.0d0
zcm=0.0d0

do iz=1,nz
  zz = z(iz)
  do iy=1,ny
    yy = y(iy)
    do ix=1,nx
      xx  = x(ix)
      aux = rho(ix,iy,iz)
      xcm = xcm+aux*xx
      ycm = ycm+aux*yy
      zcm = zcm+aux*zz
    end do
  end do
end do

xcm = xcm*dxyz/npart
ycm = ycm*dxyz/npart
zcm = zcm*dxyz/npart

return

end
