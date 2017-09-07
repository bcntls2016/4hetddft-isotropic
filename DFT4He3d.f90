program DFT4HeImpd

!
!                         ******************
!                         *** DFT4HeImpd ***
!                         ******************
!
!

!
! This program computes time dependent Helium_4 drops with/without impurities using
! Density Functional Theory.
! The interation can be Orsay-Paris or Orsay-Trento
! The code is fully 3-dimensional
!
!-----------------------------------------------------------------------------
! Version 0  (alpha)   Barcelona April-19,2004   R. Mayol & M. Pi
! Version 99 (Epsilon) Barcelona April-19,2006   R. Mayol & M. Pi
!-----------------------------------------------------------------------------
!07oct2004/11:30AM  
!06oct2004/11:30AM
!23sep2004/11:30AM
!Feb2005
! (LE FALTA: Ficheros binarios)
!-----------------------------------------------------------------------------
use seleccio_de_potencial
use Impressio
use alphasterm
use Para_DerivnD
use deriva
use energies
use rho
use field
use fftmodule
use grid 
use gridk
use classicimp
use rkpc
use impur, only:  r_cutoff, selec, umax
use lenard4
use he4
use util1
use work1

implicit none

Type (info_printout) pr


real (kind=4) :: t0,t1,t2,t3,t4,t5,t6   ! Variables used to printout run time execution.

logical              :: lfilepv         ! T-> print save file when change Paflo parameter.
logical              :: lpaflv=.false.  ! T-> allows change of Paflov coeffient
logical              :: lrkpc=.true.    ! T-> allows to use diferent evolution procedures
logical              :: lrk=.false.     ! T-> allows to only Runge-Kutta method
integer    (kind=4)  :: ndmax=2         ! maxima derivada a calcular
integer    (kind=4)  :: naux            ! Auxiliar variable
integer    (kind=4)  :: nstepp=1        ! Number of 'Paflov parameter'.
integer    (kind=4)  :: ix ,iy ,iz      ! Used as indexs of arrays
integer    (kind=4)  :: ipx,ipy,ipz     ! Used as indexs of arrays
integer    (kind=4)  :: iter,niter      ! Control the number of iterations
integer    (kind=4)  :: pener=50        ! Computes energy each 'pener' iterations
integer    (kind=4)  :: pdenpar=50      ! Writes partial densities each 'pdenpar' iterations
integer    (kind=4)  :: pchem=50        ! Writes partial chemical potential each 'pchem' iter.
integer    (kind=4)  :: ppot=1          ! Computes the potential each 'ppot' iterations
integer    (kind=4)  :: tstgrid
integer    (kind=4)  :: norder          ! Order of Taylor expansion for time evolution
integer    (kind=4)  :: nv_iter         ! Number of iterationsfor the first time evolution
integer    (kind=4)  :: iteraux		! Iter
real       (kind=8)  :: Ktops=7.63822291d0
real       (kind=8)  :: pstoK=1.d0/7.63822291d0
real       (kind=8)  :: norma           ! Use for normalization
real       (kind=8)  :: denxmax         ! Maximum value for density of the impurity
real       (kind=8)  :: p,px2,py2,pz2   ! Temporary variables for momentum values
real       (kind=8)  :: epsx,epsxerr    ! Value of autovalue and associated error
real       (kind=8)  :: errmu4          ! Relative change betwen iteration for chemical potential
real       (kind=8)  :: deltat=0.d0     ! Time step
real       (kind=8)  :: deltatps=0.d0   ! Time step in picoseconds
real       (kind=8)  :: xlamda=0.0d0    ! Monopol costraint
real       (kind=8)  :: xlamdax=0.0d0    ! He Dipol costraint x direction
real       (kind=8)  :: xlamday=0.0d0    ! He Dipol costraint y direction
real       (kind=8)  :: xlamdaz=0.0d0    ! He Dipol costraint z direction
real       (kind=8)  :: xlamdax_x=0.0d0    ! Impurity Dipol costraint x direction
real       (kind=8)  :: xlamdax_y=0.0d0    ! Impurity Dipol costraint y direction
real       (kind=8)  :: xlamdax_z=0.0d0    ! Impurity Dipol costraint z direction
real       (kind=8)  :: eold            ! Auxiliar variables
real       (kind=8)  :: aux,aux1,aux2,aux3             ! Auxiliar variables
real       (kind=8)  :: temps           ! Absolute time
real       (kind=8)  :: time0 = -1      !Starting time
real       (kind=8)  :: aux4,aux5,aux6                 ! Auxiliar variables
integer    (kind=4),allocatable  :: nitera(:) ! In wich iteration change to the
real       (kind=8),allocatable  :: paflv(:)  !    corresponging Paflov coeffiecient
real       (kind=8)  :: cnorm           !    Normalization constant
real       (kind=8)  :: pmod1        !    Work variable for controlling p-values
real       (kind=8)  :: pmod2        !    Work variable for controlling p-values
complex    (kind=8)  :: ci=(0.d0,1.d0), caux !  Work complex variables
real       (kind=8)  :: xcm4,ycm4,zcm4,xcmx,ycmx,zcmx ! Center of mass Drop and Impurity
real       (kind=8)  :: xcm,ycm,zcm ! Center of mass Drop 
real       (kind=8)  :: xlx,xly,xlz ! Angular momentum 
real       (kind=8)  :: distx,disty,distz  ! Distance between center of masses
real       (kind=8)  :: errHe, errimp,errvimp     ! Error evolution (only form Predictor-Corrector-Modificator)
real       (kind=8)  :: Zsurf = -25.d0, Morse_HeTiO2_1D, Morse_HeTiO2_3D, auxn4 ! Position of the surface

integer    (kind=4)  :: n4=300         ! Number of helium_4 atoms
integer    (kind=4)  :: mode=0          ! Way to start the program (see readen subroutine)
integer    (kind=4)  :: instate=0          ! Way to start the program (see readen subroutine)
integer    (kind=4)  :: iter0=0         ! Starting point for iteration procedure
integer    (kind=4)  :: ncurr,pcurr,k,m
integer    (kind=4)  :: icurr = 0

character  (len=40)  :: title         = 'Helium4 - 3dim.  '
character  (len=60)  :: fileout       = 'DFT4He3d.res'
character  (len=60)  :: filedenin     = 'he4.input.density'
character  (len=60)  :: filedenout    = 'he4.output.density'
character  (len=60)  :: fileimpin     = 'X.input.wf'
character  (len=60)  :: fileimpout    = 'X.output.wf'
character  (len=60)  :: filemonopol   = 'monopol.out'
character  (len=60)  :: filerimp      = 'rimp.out'
character  (len=60)  :: filevimp      = 'vimp.out'
character  (len=60)  :: fileaimp      = 'aimp.out'
character  (len=60)  :: filequadrupol = 'quadrupol.out'
character  (len=60)  :: namefile,namefile1

character  (len=23)  :: curvfile
character  (len=4)   :: chariter


logical              :: lsurf=.false.         ! include TiO2 surface or not
logical              :: lsurf3D=.false.         ! include TiO2 surface or not

real       (kind=8)  :: Lambdah=2.d0,txmean=1000.d0,txsurf=2.d0      &
                                    ,tymean=1000.d0,tysurf=2.d0      &
                                    ,tzmean=1000.d0,tzsurf=2.d0

! added by manu (see the removed stuff file modules.f90)
integer (kind=4)              :: nthread   ! Number of threads

!Variables mesures du temps
integer                      ::  ir,t_1, t_2
real(kind=4)                 :: temps_ela, t_cpu, t_cpu_0, t_cpu_1


!....................Variables read in a NAMELIST statement ..............................

namelist /input/title,fftwplan,nthread,nsfiles,                         &
                fileout,filemonopol, filerimp,filevimp,fileaimp,        &
                filedenin,filedenout,                                   &
                fileimpin,fileimpout,                                   &
                n4,N_imp,mode,                                          &
                nx,ny,nz,xmax,ymax,zmax,                                &
                xc,yc,zc,afermi,                                        &
                eps4,sigma4,core4,l,                                    &
                cp4,cpp4,den4c,alphas,h2o2m4,                           &
                denmin,psimin,npd,ndmax,icon,                           &
                niter,printpot,pchem,irespar,                           &
                pdenpar,pener,ppot,lsolid,                              &
                deltat,lrkpc,lrk,xlamda,lselection,                     &
                xlamdax,xlamday,xlamdaz,deltatps,                       &
                xlamdax_x,xlamdax_y,xlamdax_z,iter0,time0,Zsurf,        &
                pcurr,icurr,lsurf,lsurf3D,Lambdah,                      &
                txmean,txsurf,tymean,tysurf,tzmean,tzsurf

namelist /imp/rimp,vimp,m_imp_u,selec_gs_k,r_cutoff_gs_k,umax_gs_k,		&
				selec_gs_k_k,r_cutoff_gs_k_k,umax_gs_k_k,				&
				drselec_gs_k_k,drr_cutoff_gs_k_k,drumax_gs_k_k,			&
				filerimp_k, filevimp_k, fileaimp_k

!................................ Start main Program ..............................
call timer(t0)
!.............................................
!... Inicializate some numerical constants ...
!.............................................

pi     = 4.0d0*datan(1.0d0) ! Initialization of pi
twopi  = 2.0d0*pi
fourpi = 4.0d0*pi
piq    = pi*pi


!...............................
!... Read  master-input file ...
!...............................
read(5,nml=input,end=999)

open(10,file="DFT4He3d.namelist.read")
write(10,nml=input)


Allocate(pr%psi(nx,ny,nz))
Allocate(pr%selec_gs_k(N_imp))
Allocate(pr%rimp(N_imp,3))
Allocate(pr%vimp(N_imp,3))
Allocate(pr%selec_gs_k_k(N_imp,N_imp))
Allocate(pr%drselec_gs_k_k(N_imp,N_imp))
Allocate(pr%r_cutoff_gs_k(N_imp))
Allocate(pr%r_cutoff_gs_k_k(N_imp,N_imp))
Allocate(pr%drr_cutoff_gs_k_k(N_imp,N_imp))
Allocate(pr%umax_gs_k(N_imp))
Allocate(pr%umax_gs_k_k(N_imp,N_imp))
Allocate(pr%drumax_gs_k_k(N_imp,N_imp))


nn(1)  = nx ; nn(2)  = ny ; nn(3)  = nz;                ! Initialize things for PDERG
mmx(1) = nx ; mmx(2) = ny ; mmx(3) = nx ; mmx(4) = ny   ! (NO NOT MOVE THAT!!!!!!!)

!.............................................................
!.. Check if the size of the grid is among the valid values ..
!.............................................................

if(tstgrid(nx).ne.0 ) stop 'SEVERE ERROR: NX is not correct'
if(tstgrid(ny).ne.0 ) stop 'SEVERE ERROR: NY is not correct'
if(tstgrid(nz).ne.0 ) stop 'SEVERE ERROR: NZ is not correct'

!...................................................
!.. Controls Paflov parameters (read and storage ...
!...................................................

close(5)
close(10)

!.........................................
!.. Some consistency check on Delta t  ...
!.........................................
If(deltat.eq.0.d0)then
 if(deltatps.eq.0.d0)then    ! deltat ==  0, deltatps ==  0
  print*,'You must specify either Deltat (in kelvin) or Deltatps (in picosecond)'
  STOP
 else                        ! deltat ==  0, deltatps =/= 0
  deltat = deltatps/7.63822291d0
 endif
Else
 if(deltatps.eq.0.d0)then    ! deltat =/= 0, deltatps ==  0
  deltatps = deltat*7.63822291d0
 else                        ! deltat =/= 0, deltatps =/= 0
  if(deltatps.ne.(deltat*7.63822291d0))then
   print *,'Inconsistent deltat - deltatps'
   STOP
  endif
 endif
Endif
write(*,*)'Time step is ',deltat,' kelvins or ',deltatps,' picoseconds.'


!................................................
!.. Some consistency check on input variables ...
!................................................

nthread=abs(nthread)

Call Init_deriv_p(npd,ndmax,nthread)

hx    = 2.0d0*abs(xmax)/(nx)  ! Step in x-grid
hy    = 2.0d0*abs(ymax)/(ny)  ! Step in y-grid
hz    = 2.0d0*abs(zmax)/(nz)  ! Step in z-grid

dxyz  = hx*hy*hz              ! Element of volum in real space
nxyz  = nx*ny*nz              ! Total number of points

hpx   = 1.0d0/(nx*hx)         ! Step in px-grid
hpy   = 1.0d0/(ny*hy)         ! Step in py-grid
hpz   = 1.0d0/(nz*hz)         ! Step in pz-grid

pmaxx = 1.0d0/(2.0d0*hx)      ! Maximum 'frequency' in X-grid
pmaxy = 1.0d0/(2.0d0*hy)      ! Maximum 'frequency' in Y-grid
pmaxz = 1.0d0/(2.0d0*hz)      ! Maximum 'frequency' in Z-grid


!...............................
!.. Dimensionate main ARRAYS ...
!...............................

call dimen()
open(1,file="imp.input")
read(1,nml=imp)
close(1)

write(6,nml=imp)

Write(6,'("Used potentials....:")')
do k=1,N_imp
  pr%selec_gs_k(k)=selec_gs_k(k)
  pr%r_cutoff_gs_k(k)=r_cutoff_gs_k(k)
  pr%umax_gs_k(k)=umax_gs_k(k)
  Write(6,'(A/,"selc_gs,r_cutoff_gs,umax_gs..:"1p,2E15.6)')selec_gs_k(k),r_cutoff_gs_k(k),umax_gs_k(k)
  m_imp(k) = m_imp_u(k)*mp_u
  do m=k+1,N_imp
    pr%selec_gs_k_k(k,m)=selec_gs_k_k(k,m)
	pr%r_cutoff_gs_k_k(k,m)=r_cutoff_gs_k_k(k,m)
	pr%umax_gs_k_k(k,m)=umax_gs_k_k(k,m)
    pr%drselec_gs_k_k(k,m)=drselec_gs_k_k(k,m)
	pr%drr_cutoff_gs_k_k(k,m)=drr_cutoff_gs_k_k(k,m)
	pr%drumax_gs_k_k(k,m)=drumax_gs_k_k(k,m)
    Write(6,'(A/,"selec_gs_k_k,r_cutoff_gs_k_k,umax_gs_k_k..:"1p,2E15.6)')selec_gs_k_k(k,m),r_cutoff_gs_k_k(k,m),umax_gs_k_k(k,m)
    Write(6,'(A/,"drselec_gs_k_k,drr_cutoff_gs_k_k,drumax_gs_k_k..:"1p,2E15.6)')drselec_gs_k_k(k,m),drr_cutoff_gs_k_k(k,m),drumax_gs_k_k(k,m)
  enddo
enddo
!.....................................
!.. Initial value of rimp and vimp ...
!.....................................
vimp = vimp*7.63822291d0 ! Because it is given in A/ps, we have to transform to A*K


!................................
!... Build grid in real space ...
!................................

do ix=1,nx  !.................... Grid X
 x(ix) = -xmax+hx*(ix-1)
end do
do iy=1,ny  !.................... Grid Y
 y(iy) = -ymax+hy*(iy-1)
end do
do iz=1,nz  !.................... Grid  Z
 z(iz) = -zmax+hz*(iz-1)
end do

!....................................
!... Build grid in momentum space ...
!....................................

!.... Build p-grid. In order to use FFTW the grid must
!     start from frequency zero to the maximum and then continue from
!     (-maximum) to zero (negative).

!............................................ grid Px
do ipx=1,nx/2+1
   px(ipx) =        hpx*(ipx-1)
end do
do ipx=nx/2+2,nx
   px(ipx) = -pmaxx+hpx*(ipx-(nx/2+1))
end do

!............................................ grid Py
do ipy=1,ny/2+1
   py(ipy) =        hpy*(ipy-1)
end do
do ipy=ny/2+2,ny
   py(ipy) = -pmaxy+hpy*(ipy-(ny/2+1))
end do

!............................................ grid Pz
do ipz=1,nz/2+1
   pz(ipz) =        hpz*(ipz-1)
end do
do ipz=nz/2+2,nz
   pz(ipz) = -pmaxz+hpz*(ipz-(nz/2+1))
end do

!............................................ Compule modulus of p
do ipz=1,nz
  pz2=pz(ipz)**2
  do ipy=1,ny
    py2=py(ipy)**2
    do ipx=1,nx/2+1
      px2               = px(ipx)**2
      pmod(ipx,ipy,ipz) = sqrt(px2+py2+pz2)
    end do
  end do
end do

pmod1=maxval(pmod)
pmod2=sqrt(pmaxx**2+pmaxy**2+pmaxz**2)


write(6,*) '    Initialize Linear Interpolation for V_Pi and V_Sig'
call flush(8)




call potenimpini() ! interpolation + first call to updatepoten

!................................
!... read density or build-it ...
!................................

call readenc(n4,densat4,filedenin,fileimpin,mode)
!...............................................................test
 do iz=1,nz ; do iy=1,ny ; do ix=1,nx 
  if(.not.(den(ix,iy,iz).gt.0))print*,ix,iy,iz,den(ix,iy,iz)
 end do ; enddo ; enddo
!...............................................................test


call respar(x,y,z,nx,ny,nz,1,'hedenini','hedenini',den,den)

!....................................
!.. Print-out initial quantities  ...
!....................................

!open(6,file=fileout)

write(6,6010) title

select case(mode)
   case(0) !................................... Start a dynamic form static calcultaion
         write(6,6011) filedenin,filedenout
   case(2) !................................... Continue a dynamic calculation from a prevous one
      write(6,6013) filedenout,fileimpout
   case(3) !................................... Continue a dynamic calculation from impurity with excited internal state
      write(6,6013) filedenout,fileimpout
   case(4) !................................... Continue a dynamic calculation from impurity with excited internal state
      write(6,6013) filedenout,fileimpout
   case default !.............................. Start still not programed.
      write(6,*) ' '
      write(6,*) ' The variable mode has no acceptable value.'
      write(6,*) ' '
      stop
end select

!...............................................................
!................................ Write the input parameters ...
!...............................................................

write(6,6018) nthread,niter
if(mode.ne.0) then
  write(6,6020) n4
else
  write(6,6025) n4
end if
write(6,6030) nx,ny,nz,hx, hy, hz, x(1) ,y(1) ,z(1) ,x(nx),y(ny),z(nz)
write(6,6035)          hpx,hpz,hpz,px(1),py(1),pz(1),pmaxx,pmaxy,pmaxz,&
                       pmod1,pmod2
write(6,6037) cp4,cpp4,den4c,alphas,l,den0s,h2o2m4


!...............
!.. Impurity ...
!...............

!...................................................................
!... Compute the FFT of Lennard-Jones                            ...
!... Prepara \alpha_s term in case of Orsay-Trento Interaction.  ...
!...................................................................

If(Lsolid)core4='OT '

select case(core4)
   case('OP ')
     h4=h4op
     write(6,*) '    Use Orsay-Paris-Barcelona Interaction.'
     write(6,6040) core4,h4,eps4,sigma4,b4
   case('OT ')
     h4=h4ot
     write(6,*) '    Use Orsay-Trento Interaction. (ONLY CORE)'
     write(6,*) '    Do not calculate Alpha_s in the field neither energy.'
     If(Lsolid)write(6,'("    We will use the solid functional (see: PRB72, 214522(2005)) ")')
     write(6,6040) core4,h4,eps4,sigma4,b4
   case('OTC')
     h4=h4ot
     write(6,*) '    Use Orsay-Trento Interaction.'
     write(6,*) '    Full Orsay-Trento calculation. (Field and Energy)'
     write(6,6040) core4,h4,eps4,sigma4,b4
     allocate( denalf(nx  ,ny,nz))                                                            
     allocate(  falfs(nx  ,ny,nz))
     allocate(kalfs(nx/2+1,ny,nz))
     allocate(intxalf(nx  ,ny,nz))
     allocate(intyalf(nx  ,ny,nz))
     allocate(intzalf(nx  ,ny,nz))
     allocate(ualphas(nx  ,ny,nz))
   case default
     print *,' ***************** WARNING ************************'
     print *,' '
     print *,' I do not know how to work with the core: ',core4
     print *,' '
     print *,' **************************************************'
     print *,' '
     STOP ' ... The program stops due to a severe error.'
end select



!...............................
!... Prepare plans for FFTWs ...
!...............................
write(6,*) '    Initialize Plans for FFT.'
!call fftini(nx,ny,nz)
call fftini(nx,ny,nz,nthread)

!...........................................................
!... Form  factor for Lennard-Jones and for the impurity ...
!...........................................................

write(6,*) '    Compute the FFT of the kernel of Lennard-Jones integrals.'

call fforma(core4,b4,eps4,sigma4,h4,nx,ny,nz,pmod,fvlj4)


!........................................
!.. Initialize coarse-graining kernel ...
!........................................


If(Lsolid)Then

  aux=h4*1.065d0
  write(6,'("    Initialize Coarse-graining kernel, for Solid DF, h_cg=",1p, E15.6)')aux
  call initcg(aux,wcgk)

Else

  write(6,'("    Initialize Coarse-graining kernel, h_cg=",1p, E15.6)')h4
  call initcg(h4,wcgk)

Endif

if(core4.eq.'OTC') then
   forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
      kalfs(ix,iy,iz) = exp(-(pi*l*pmod(ix,iy,iz))**2)
   end forall
end if

!
!    We compute the perturbate wave function
!
! Initial velocity : Altough xlamdax_i is in K units, xlamda
! is introduced in Angstrom/picosecond for the sake of lazyness.
auxn4 = sum(den)*dxyz

if(xlamda.ne.0.d0)then
write(*,*)'xlambda not equal zero'
xlamdaz = xlamda*7.63822291d0/(2.d0*h2o2m4)
endif

If(mode.eq.0.And.Ldensity)then
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         aux=x(ix)*xlamdax       &
            +y(iy)*xlamday       &
            +z(iz)*xlamdaz
         psi(ix,iy,iz) = sqrt(den(ix,iy,iz)) &
                       * cmplx(cos(aux),sin(aux))
       end do
     end do
   end do
Endif
!.................................
!.. First call to total energy ...
!.................................

den=Conjg(psi)*psi
call potenimp()
call poten()              ! First Potential  (for Lagrange Equation)
call forceimp()
aimp(:,1) = F(:,1)/m_imp(:)
aimp(:,2) = F(:,2)/m_imp(:)
aimp(:,3) = F(:,3)/m_imp(:)
call energy()             ! Calculate energies


! TEST: Treu el valor de uimp
call respar(x,y,z,nx,ny,nz,1,'uimp','den',uimp,den)
write(6,'("Number of He4 atoms",1P,E15.6)')auxn4
write(6,6050) etot4,etot4/n4,ekin4,elj4,ealphas,esolid,ecor4
  write(6,6060) eimpu_impu,eimpu,ekinx,eHeX,0.d0,etot
eold = etot4
call flush(6)

!-------------------------------------------------------------------------------
!---                            Iterative procedure                           --
!-------------------------------------------------------------------------------

! TIME CONSTANT !
! This time it's a cylinder.
do iz=1,nz
 do iy=1,ny
  do ix=1,nx
   timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((abs(x(ix))-txmean)/txsurf)            &
                                 +1.d0+tanh((abs(y(iy))-tymean)/tysurf)            &
                                 +1.d0+tanh((abs(z(iz))-tzmean)/tzsurf))           &
                                 ,1.d0)
  enddo
 enddo
enddo


!plot it
! call respar(x,y,z,nx,ny,nz,1,'timec','den',timec,den)
open(unit=32,file='timec.dat')
do iz=1,nz
 do ix=1,nx
  write(32,*)x(ix),z(iz),real(timec(ix,ny/2+1,iz))
 enddo
 write(32,*)''
enddo
close(32)

do k=1,N_imp
  open(unit=(100+k), file=filerimp_k(k))
  open(unit=(200+k), file=filevimp_k(k))
  open(unit=(300+k), file=fileaimp_k(k))
enddo

lfilepv  = .false.

call timer(t5)

open(unit=42,file='uext.dat')
 do iz=1,nz
     write(42,*)z(iz),uext(3,5,iz)
 enddo
close(42)

qr = 0.d0
qv = 0.d0
rimpold = 0.d0
vimpold = 0.d0
aimpold = 0.d0

do k=1,N_imp
  write((100+k),'("# Tiempo(ps), x(AA), y(AA), z(AA)")')
  write((100+k),'(1x,1p,E15.6,3E18.10)')time0, rimp(k,1), rimp(k,2), rimp(k,3)
  write((200+k),'("# Tiempo(ps), Vx(AA/ps), Vy(AA/ps), Vz(AA/ps)")')
  write((200+k),'(1x,1p,E15.6,3E18.10)')time0, vimp(k,1)*pstoK, vimp(k,2)*pstoK, vimp(k,3)*pstoK
  write((300+k),'("# Tiempo(ps), Ax(AA/ps**2), Ay(AA/ps**2), Az(AA/ps**2)")')
  write((300+k),'(1x,1p,E15.6,3E18.10)')time0, aimp(k,1)*pstoK*pstoK, aimp(k,2)*pstoK*pstoK, aimp(k,3)*pstoK*pstoK
enddo

pr%it0   = iter0
pr%time0 = time0
pr%dtps  = deltatps
pr%nx    = nx
pr%ny    = ny
pr%nz    = nz
pr%hx    = hx
pr%hy    = hy
pr%hz    = hz
pr%xmax  = xmax
pr%ymax  = ymax
pr%zmax  = zmax
pr%N_imp = N_imp

!mesure du temps
! Temps CPU de calcul initial.
call cpu_time(t_cpu_0)
! Temps CPU de restitution initial.
call system_clock(count=t_1, count_rate=ir)

iter0=iter0+1
do iter=iter0,niter  ! <--------- Iterative procedure starts here.

Iteraux = iter - iter0 + 1

pr%it = iter

    if((iter-iter0+1).le.3.Or.lrk)then
      call steprk(deltat)

    else
      call steppc(deltat,errHe,errimp,errvimp)

      write(6,'(" Error( He, imp) (From Steppc)...",1p,3E15.6)')errHe,errimp,errvimp
    endif

    call potenimp()
    call poten()
    call forceimp()
    aimp(:,1) = F(:,1)/m_imp(:)
    aimp(:,2) = F(:,2)/m_imp(:)
    aimp(:,3) = F(:,3)/m_imp(:)

    aux1 = time0 + Iteraux*deltatps
    temps = aux1

	do k=1,N_imp
  	  write((100+k),'(1x,1p,E15.6,3E18.10)')aux1, rimp(k,1), rimp(k,2), rimp(k,3)
  	  write((200+k),'(1x,1p,E15.6,3E18.10)')aux1, vimp(k,1)*pstoK, vimp(k,2)*pstoK, vimp(k,3)*pstoK
  	  write((300+k),'(1x,1p,E15.6,3E18.10)')aux1, aimp(k,1)*pstoK*pstoK, aimp(k,2)*pstoK*pstoK, aimp(k,3)*pstoK*pstoK
	enddo

   if(mod(iter,pener).eq.0) then          ! Compute New energy and max of density

      Write(6,'(" Iteration....:",I10)')iter
      auxn4 =sum(den)*dxyz
      Write(6,'(" Number of particles....:",1p,E18.10)')auxn4

      call energy()

      write(6,7010) etot4,(etot4-eold),etot4/n4,ekin4,elj4,ealphas,esolid,ecor4
      write(6,7015) eimpu_impu,eimpu,ekinx,eHeX,0.d0,etot

      eold = etot4

        call r_cm(den,n4,xcm4,ycm4,zcm4)    ! Center of mass of 4He Drop

          xcm = xcm4; ycm=ycm4; zcm=zcm4

          Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
          Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
          Call derivnD(1,nn,hz,3,psi,sto3c,Icon)
!          
! Z Component of angular momentum 
!          
          caux = (0.d0, 0.d0)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                ((y(iy)-ycm)*sto1c(ix,iy,iz) - (x(ix)-xcm)*sto2c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          xlz = caux*dxyz
!          
! Y Component of angular momentum 
!          
          caux = (0.d0, 0.d0)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                ((x(ix)-xcm)*sto3c(ix,iy,iz) - (z(iz)-zcm)*sto1c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          xly = caux*dxyz
!          
! X Component of angular momentum 
!          
          caux = (0.d0, 0.d0)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                ((z(iz)-zcm)*sto2c(ix,iy,iz) - (y(iy)-ycm)*sto3c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          xlx = caux*dxyz
          
        Write(6,'("<Lx,Ly,Lz>.......:",1p,3E20.11)')xlx,xly,xlz
        
          aux1 = 0.d0
          aux2 = 0.d0
          aux3 = 0.d0
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                aux1 = aux1 + den(ix,iy,iz)*x(ix)**2
                aux2 = aux2 + den(ix,iy,iz)*y(iy)**2
                aux3 = aux3 + den(ix,iy,iz)*z(iz)**2
              EndDo
            EndDo
          EndDo
          aux1 = aux1*dxyz
          aux2 = aux2*dxyz
          aux3 = aux3*dxyz
          
        pr%r2(1)   = aux1
        pr%r2(2)   = aux2
        pr%r2(3)   = aux3
        pr%ang(1)  = xlx
        pr%ang(2)  = xly
        pr%ang(3)  = xlz
        pr%cm(1)   = xcm4
        pr%cm(2)   = ycm4
        pr%cm(3)   = zcm4
        pr%ekin    = ekin4
        pr%elj     = elj4
        pr%ealphas = ealphas
        pr%esolid  = esolid
        pr%ecor    = ecor4
        pr%auxn4   = auxn4
        pr%ekinx   = ekinx
        pr%evx     = eimpu
        pr%etot    = etot
        pr%time    = temps
		pr%rimp(:,:)  = rimp(:,:)
		pr%vimp(:,:)    = vimp(:,:)
        write(6,7100) xcm4,ycm4,zcm4

   end if

!..............................................................................

   if(mod(iteraux,pcurr).eq.0) then        ! Save wavefunction for current

     ncurr = iteraux/pcurr + icurr
      select case (ncurr)
       case(1:9)
!      write(chariter,8011)ncurr
         write(chariter,'("000",I1)')ncurr
       case(10:99)
!      write(chariter,8012)ncurr
         write(chariter,'("00",I2)')ncurr
       case(100:999)
!      write(chariter,8013)ncurr
         write(chariter,'("0",I3)')ncurr
       case(1000:9999)
         write(chariter,'(I4)')ncurr
      end select
      namefile='density.'//chariter//'.dat'
      pr%namefile = namefile
      pr%psi(:,:,:) = psi(:,:,:)
       call printoutc(pr)
   endif

!..............................................................................

!..............................................................................

   call timer(t6)                         ! Compute use time
   t5=t6
end do

!Fin mesure temps
! Temps CPU de restitution final
call system_clock(count=t_2, count_rate=ir)
temps_ela=real(t_2 - t_1,kind=4)/real(ir,kind=4)
! Temps CPU de calcul final
call cpu_time(t_cpu_1)
t_cpu = t_cpu_1 - t_cpu_0
! affichage temps 
print*, "elapsed partie iterative : ", temps_ela
print*, "cpu_time partie iterative : ", t_cpu
      pr%namefile = filedenout
      pr%psi(:,:,:) = psi(:,:,:)
       call printoutc(pr)
call timer(t4)
print *,' Total  ',t4-t0

stop
999 stop 'DFT3He3d. Error in input master file. Too short'

!...............
!... Formats ...
!...............

3100 format(3x,0P,f9.4,2x,1P,E13.5)

3156 format(10E13.5)

6010 format(//,&
T10,'   ######  ####### ####### #       #     #          #####          ',/,  &
T10,'   #     # #          #    #    #  #     #  ###### #     #  #####  ',/,  &
T10,'   #     # #          #    #    #  #     #  #            #  #    # ',/,  &
T10,'   #     # #####      #    #    #  #######  #####   #####   #    # ',/,  &
T10,'   #     # #          #    ####### #     #  #            #  #    # ',/,  &
T10,'   #     # #          #         #  #     #  #      #     #  #    # ',/,  &
T10,'   ######  #          #         #  #     #  ######  #####   #####  ',//, &
T6,'Title of the run: ',A)

6011 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input  densitity file: ',A,/,&
               T6,'Output densitity file: ',A)

6111 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input file with helium densitity    : ',A,/,&
               T6,'Input file with impurity wave func. : ',A,/,&
               T6,'Output file with helium densitity   : ',A,/,&
               T6,'Output file with impurity wave func.: ',A,/,' ')

6012 format(//,T6,'Start a new calculation:',//,&
               T6,'Output densitity file: ',A)
6013 format(//,T6,'Start a new calculation with an impurity:',//,&
               T6,'Output file for Helium density: ',A,/,        &
               T6,'Output file for the impurity wave function: ',A)
6018 format(//,T6,'Number of threads:    ',I6,/,&
               T6,'Number of iterations: ',i6)
6020 format(//,T6,'Number of particles:    ',0P,I10,/,&
               T6,'Radius of the cluster : ',F10.3,' A')
6025 format(//,T6,'Number of particles:    ',0P,I10,/,' ')
6030 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| REAL GRID         |     X-grid       Y-grid       Z-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Number of points  |',0P,T32,I4,T45,I4,T58,I4,T66,' |',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+')
6035 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| MOMEMTUM GRID     |    Px-grid      Py-grid      Pz-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+' ,//,&
               T6,'Maximum modulus of p ',0P,F11.6,3x,F11.6,/,' ')
6037 format(//,T6,'Parameters for the funcional:',//,          &
               T8,'cp4 ...... ',1P,E13.5 ,' K \AA**6'      ,/, &
               T8,'cpp4 ..... ',1P,E13.5 ,' K \AA**9'      ,/, &
               T8,'den4c .... ',0P,F13.3 ,' \AA**{-3}'     ,/, &
               T8,'Alphas ... ',0P,F13.3 ,' K ^-1 \AA**3'  ,/, &
               T8,'L ........ ',0P,F13.2 ,' \AA'           ,/, &
               T8,'den0s .... ',0P,F13.2 ,' \AA**-3'       ,/, &
               T8,'h2o2m4 ... ',0P,F14.11,' hbar**2 / (2 m_4)' )
6138 format(' ',T8,'h2o2mx ... ',0P,F14.11,' hbar**2 / (2 m_x)' )
6038 format(//,T6,'Change of Paflov parameter allowed: ',//, &
     T18,'From     to      iter     Factor',/, &
     T19,'------  ------  ------  -----------')

6039 format(1x,0P,T17,i6,T25,I6,T33,i6,t42,f11.7)

6150 format(//,T6,'Pavlov parameter fixed for all the run to: ',F8.4)

6040 format( /,T6,'Lennard-Jones parameters:',//,&
               T10,'Core    ',A3,/,&
               T10,'h     ',F11.7,' A',/,&
               T10,'eps   ',F11.7,' K',/,&
               T10,'sigma ',F11.7,' A',/,&
               T10,'b     ',F11.3,' K A**3 '//,' ')
6050 format(//,T5,'FIRST ENERGY BALANCE: ',                    //     &
              ,T5,'TOTAL   energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Energy per particle (He) .....: ',F18.6,' K',/,    &
               T5,'Kinetic energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Lennard-Jones energy (He) ....: ',F18.6,' K',/,    &
               T5,'Alpha_s term  energy (He) ....: ',F18.6,' K',/,    &
               T5,'Solid energy (He)  ...........: ',F18.6,' K',/,    &
               T5,'Correlation energy   (He) ....: ',F18.6,' K')
6060 format(1x,T5,'Impurity-impurity energy .....: ',F18.6,' K',/,    &
			   T5,'Impurity energy (X) ..........: ',F18.6,' K',/,    &
               T5,'Kinetic energy (X) ...........: ',F18.6,' K',/,    &
               T5,'Interaction energy (X-He) ....: ',F18.6,' K',/,    &
               T5,'Spin-Orbit energy (X) ........: ',F18.6,' K',/,    &
               T5,'TOTAL ENERGY (He+X) ..........: ',F18.6,' K',/)
6065 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' A',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' A',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' A',/)
7000 format(//,1x,T2,'Iter     Mu(K)      Err(Mu)    Ttime  / Lap Time',/,&
 '--------------------------------------------------')
7001 format(//,1x,T2, &
     'Iter     Mu(K)      Err(Mu)    Autovalue(K)   err(K)   ETtime  / Lap Time',&
     /,74('-'))

7010 format(//,T5,'ITERATIVE PROCEDURE ',                                    //  &
              ,T5,'Total Energy (He).......... ',0P,F18.6,' K +- ',1P,e12.4,' K',&
             /,T5,'Energy per particle (He)... ',0P,F18.6,' K',/,                &
             /,T5,'Kinetic Energy (He)........ ',0P,F18.6,' K',                  &
             /,T5,'Lennard-Jones Energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Alpha_s term  energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Solid energy (He) ......... ',0P,F18.6,' K',                  &
             /,T5,'Correlation Energy  (He)... ',0P,F18.6,' K')
7015 format(   T5,'Impurity-impurity energy... ',0P,F18.6,' K',/,&
			   T5,'Impurity energy (X->He) ... ',0P,F18.6,' K',/,&
               T5,'Kinetic energy (X) ........ ',0P,F18.6,' K',/,&
               T5,'Interaction energy (X-He) . ',0P,F18.6,' K',/,    &
               T5,'Spin-Orbit energy (X) ..... ',0P,F18.6,' K',/,    &
               T5,'TOTAL energy (He+X) ....... ',0P,F18.6,' K',/,' ')

7016 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' K',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' K',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' K',/,' ')

7017 format(   T5,'Chemical Potential ........ ',0P,F18.6,' K +- ',1P,e12.4,'K')
7018 format(   T5,'Autovalue (impurity) ...... ',0P,F18.6,' K +- ',1P,e12.4,'K')

! 7100 format(/1x,T5,'Center of Mass of the Helium ...(', &
!                      0P,F10.6,',',F10.6,',',F10.6,') A')
7100 format(/1x,T5,'Center of Mass of the Helium ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A')
7110 format(/1x,T5,'Center of Mass of the Helium ........(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Center of Mass of the Impurity ......(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Distances between centers of mass ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,' ')

7020 format(1x,1P,2E14.5)
7030 format(0P,I5,T7,F13.7,T21,1P,E9.2,T32,0P,F8.0,'/',F8.2)
7035 format(0P,I5,T7,F13.7,T21,1P,E9.2,T32,0P,F13.7,T46,1P,E9.2,T57,0P,F8.0,'/',F8.2)


8010 format('partial.density' ,i1)
8015 format('partial.densityx',i1)
8020 format('partial.density' ,i2)
8025 format('partial.densityx',i2)
8030 format('partial.density' ,i3)
8035 format('partial.densityx',i3)

8011 format('00',i1)
8012 format('0',i2)
8013 format(i3)

5010 format('density.',SS,i1,'.out')
5020 format('density.',SS,i2,'.out')
5030 format('density.',SS,i3,'.out')
5040 format('density.',SS,i4,'.out')
5050 format('density.',SS,i5,'.out')
5060 format('density.',SS,i6,'.out')

5015 format('densityx.',SS,i1,'.out')
5025 format('densityx.',SS,i2,'.out')
5035 format('densityx.',SS,i3,'.out')
5045 format('densityx.',SS,i4,'.out')
5055 format('densityx.',SS,i5,'.out')
5065 format('densityx.',SS,i6,'.out')
!         1         2         3         4         5         6         7         8
!|2345678901234567890123456789012345678901234567890123456789012345678901234567890

end program 
