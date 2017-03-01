!--------------------------------------------------------------------
!---                  Subroutine FForma                           ---
!--------------------------------------------------------------------

! Gives the FFT of lennard-Jones potential.
! Thre are two cores possible: 
!    OP -> Orsay-Paris  core.
!    OT -> Orsay-Trento core.


subroutine fforma(core,b,eps,sigma,hc,nx,ny,nz,pmod,fvlj)
implicit none

integer   (kind=4) :: nx,ny,nz
real      (kind=8) :: eps,b,sigma,hc,p
real      (kind=8) :: pmod(nx/2+1,ny,nz)
real      (kind=8) :: fvlj(nx/2+1,ny,nz)
character  (len=2) :: core

real      (kind=8)  :: pi,twopi
real      (kind=8)  :: aux1,aux2,aux3,sigma6,v0
real      (kind=8)  :: t1,t2,t3
integer   (kind=4)  :: ipx,ipy,ipz

!........................
!... Useful constants ...
!........................

pi    = 4.0d0*datan(1.0d0)
twopi = 2.0d0*pi

!...............................................................................
!... calculo de la transformada del potencial en los puntos de la grid de las p
!...............................................................................

sigma6 = sigma**6
aux1   = 8.d0*eps
aux3   = sigma6/hc**6

do ipz=1,nz
  do ipy=1,ny
    do ipx=1,nx/2+1
       p = pmod(ipx,ipy,ipz)
       if(abs(p).lt.1.d-10) then
          fvlj(ipx,ipy,ipz) = bforce(core,eps,sigma,hc)
       else
          aux2 = twopi*p
          t1   =  sinx11(aux2,hc)
          t2   =   sinx5(aux2,hc)
          if(core.eq.'OP') then    !........................ Core Orsay-Paris
             v0                = aux1*aux3*(aux3-1.0d0)/hc**4
             t3                = intcore(aux2,hc)
             fvlj(ipx,ipy,ipz) = (aux1*sigma6*(sigma6*t1-t2)+v0*t3)/p
          else                     !........................ Core Orsay-Trento
             fvlj(ipx,ipy,ipz) = (aux1*sigma6*(sigma6*t1-t2))/p
          end if
       end if
    end do
  end do
end do
return

!
!.. Subrutinas internas del calculo de la transformada del Lennard-Jones
!   (sinxn, cosxn, sinint, cosint, factorial)

contains

!--------------------------------------------------------------------
!                  Function sinx11
!--------------------------------------------------------------------

function sinx11(a,x) result(res)
 
! Esta funcion calcula la integral \int{ \frac{\sin(ax)}{x^11} }
! entre x e infinito. La expresion fue calculada con Mathematica
 
implicit none

real    (kind=8)             :: a,x
real    (kind=8)             :: b(5),c(5),d(5)
real    (kind=8)             :: res,res1,res2
real    (kind=8)             :: y,y2
integer (kind=4)             :: ifail


y = a*x
y2= y*y

b(1) = 40320.d0 ; c(1) = 362880.d0  ; d(1) =   1.d0   
b(2) =  -720.d0 ; c(2) =  -5040.d0  ; d(2) =   y2
b(3) =    24.d0 ; c(3) =    120.d0  ; d(3) =   y2**2
b(4) =    -2.d0 ; c(4) =     -6.d0  ; d(4) =   y2**3
b(5) =     1.d0 ; c(5) =      1.d0  ; d(5) =   y2**4

res1 = (a*dot_product(b,d)*cos(y) + dot_product(c,d)*sin(y)/x) &
       /(3628800*x**9)
res2 = -(a**10*sinint(a,x))/3628800.d0

res = res1+res2
return
end function sinx11

!--------------------------------------------------------------------
!                  Function sinx5
!--------------------------------------------------------------------

function sinx5(a,x) result(res)
 
! Esta funcion calcula la integral \int{ \frac{\sin(ax)}{x^11} }
! entre x e infinito. La expresion fue calculada con Mathematica
 
implicit none

real    (kind=8), intent(in) :: a,x
real    (kind=8)             :: res1,res2,res
real    (kind=8)             :: y,y2
real    (kind=8)             :: b(2),c(2),d(2)
integer (kind=4) :: ifail


y  = a*x
y2 = y*y

b(1) = -2.d0   ; c(1) = -6.d0 ; d(1) = 1.d0
b(2) =  1.d0   ; c(2) =  1.d0 ; d(2) = y2


res1 = -(a*dot_product(b,d)*cos(y) + dot_product(c,d)*sin(y)/x) &
       /(24*x**3)
res2 =  (a**4*sinint(a,x))/24.d0


res  =  res1+res2

return
end function sinx5

!--------------------------------------------------------------------
!                  Function sinint
!--------------------------------------------------------------------

function sinint(a,x) result(res)

! Funcion seno integral definido como Schaum (viejo) 1.285 p230

implicit none

real    (kind=8)             :: res
real    (kind=8)             :: a,x
real    (kind=8)             :: y,y2,yact
real    (kind=8)             :: term,signo
real    (kind=8), parameter  :: aux0 = -1
real    (kind=8)             :: halfpi,pi
integer (kind=4), parameter  :: maxn   = 500! 'Solo 51 terminos'
integer (kind=4)             :: n,aux1
integer (kind=4)             :: ifail 

interface
  double precision function s13adf(x,ifail)
         real (kind=8) :: x
         integer (kind=4) :: ifail
  end function s13adf
end interface

pi     = 4.0d0*datan(1.0d0)
halfpi = pi*0.5d0

y     = (a*x)


!y2    = y*y
!yact  = y
!res   = 0.0d0
 
!signo = -1.0d0
!do n=1,maxn
! signo  = signo*aux0
!  aux1   = (2*n-1)
!  term   = signo*yact/(aux1*factorial(aux1))
!  res    = res+term
!  if(abs(term/res).lt.1.d-16) exit
!  yact   = yact*y2
!end do


ifail=0
res = s13adf(y,ifail)

if(a.gt.0.0d0) then
  res = halfpi-res
else
  res = halfpi+res
end if
return
end function sinint

!--------------------------------------------------------------------
!                  Function factorial
!--------------------------------------------------------------------

! Ojo no tiene limitaciones. para 100! el orden es 10**157

recursive function factorial(n) result(res)
implicit none
integer, intent(in) :: n
real    (kind=8)    :: res
if(n==1) then
  res = 1.0d0
else
  res = n*factorial(n-1)
end if
return
end function factorial


!--------------------------------------------------------------------
!                  Function intcore
!--------------------------------------------------------------------

function intcore(a,x) result(res)
 
! Esta funcion calcula la integral la transformada de fourier del 
! core Orsay-Paris del potencial de Lennard-Jones
 
implicit none

real (kind=8) :: a,x
real (kind=8) :: b(3),c(3),d(3),e(3)
real (kind=8) :: y,y2
real (kind=8) :: res

y  = a*x
y2 = y*y

b(1) = -120.d0 ; c(1) = 120.d0  ;  d(1) = x/a**5   ; e(1) = 1.d0/a**6
b(2) =   20.d0 ; c(2) = -60.d0  ;  d(2) = d(1)*y2  ; e(2) = e(1)*y2
b(3) =   -1.d0 ; c(3) =   5.d0  ;  d(3) = d(2)*y2  ; e(3) = e(2)*y2


res  = (dot_product(b,d)*cos(y) + dot_product(c,e)*sin(y) )

return

end function intcore

!-----------------------------------------------------------------
!--               Functio bforce                               ---
!-----------------------------------------------------------------
!
!.. Esta funcion integral el Potencial de Lennard-Jones con Core
!   Orsay Paris y con el core de Orsay-Trento
!   en coordenadas esfericas a todo el espacio.
!   (sirve para comprobar el valor de 'b')

function bforce(core,eps,sigma,hc) result(res)
implicit none
character (len=2) :: core
real (kind=8) :: eps4,xx,v0,s6,s12,bb1,bb2,bb3,pi,pi4,res
real (kind=8) :: eps,hc,sigma

pi      = 4.0d0*datan(1.0d0)
pi4     = 4.0d0*pi
eps4    = 4.0d0*eps
xx      = (sigma/hc)**6
v0      = xx*(xx-1.0d0)/hc**4
s6      = sigma**6
s12     = s6*s6
if(core.eq.'OP') then
  bb1     =  pi4*eps4*(v0*hc**7)/7.d0
  bb2     =  pi4*eps4*s12/(9.d0*hc**9)
  bb3     = -pi4*eps4*s6 /(3.d0*hc**3)
  res     = bb1+bb2+bb3
else
  bb2     =  pi4*eps4*s12/(9.d0*hc**9)
  bb3     = -pi4*eps4*s6 /(3.d0*hc**3)
  res     = bb2+bb3
end if
return
end function bforce

end 
