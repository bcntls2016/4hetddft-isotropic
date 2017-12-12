Module Impressio
Character (len=60) :: outdir = 'he-wfs'
Type, Public  :: info_printout
  Character (len=60) :: namefile
  Character (Len=80), Allocatable :: selec_gs_k(:), selec_gs_k_k(:,:), drselec_gs_k_k(:,:)
  Integer   (Kind=4) :: it, it0, nx, ny, nz, N_imp
  Real      (Kind=8) :: r2(3), cm(3), ang(3), dtps, Ekin, Elj, Ecor, auxn4
  Real      (Kind=8) :: Ealphas, Esolid, Ekinx, EVx, Etot
  Real      (Kind=8) :: Time0, Time, hx, hy, hz, xmax, ymax, zmax
  Real      (Kind=8), Allocatable :: r_cutoff_gs_k(:), r_cutoff_gs_k_k(:,:), drr_cutoff_gs_k_k(:,:)
  Real      (Kind=8), Allocatable :: umax_gs_k(:),umax_gs_k_k(:,:),drumax_gs_k_k(:,:)
  Real		(Kind=8), Allocatable :: rimp(:,:), vimp(:,:)
  Complex   (Kind=8), Allocatable :: psi(:,:,:)
End type info_printout

End Module Impressio        
!------------------------------------------------------------------
!---                    MODULES                                 ---
!------------------------------------------------------------------
module alphasterm
real    (kind=8) , allocatable  :: denalf(:,:,:)  ! Gaussian weighted density
real    (kind=8) , allocatable  :: falfs(:,:,:)   ! f(r) function for alpha_s term
real    (kind=8) , allocatable  :: intxalf(:,:,:) ! Intermediate integral in alpha_s terms
real    (kind=8) , allocatable  :: intyalf(:,:,:) !    "
real    (kind=8) , allocatable  :: intzalf(:,:,:) !    " 
real    (kind=8) , allocatable  :: kalfs(:,:,:)   ! FFT-kernel for the gaussian.
real    (kind=8) , allocatable  :: ualphas(:,:,:)  ! Piece of field due to the alpha_s term.
end module alphasterm
!------------------------------------------------------------------
module deriva
integer (kind=4)              :: npd=13             ! Number of points for derivatives
real    (kind=8), allocatable :: dxden(:,:,:)      ! Partial derivative in X for den in Real-space
real    (kind=8), allocatable :: dyden(:,:,:)      ! Partial derivative in X for den in Real-space
real    (kind=8), allocatable :: dzden(:,:,:)      ! Partial derivative in X for den in Real-space
integer (kind=4)              :: icon=13            ! Bounday conditions for derivatives
end module deriva
!------------------------------------------------------------------
module energies
real    (kind=8)              :: etot             ! Total energy (system)
real    (kind=8)              :: etot4            ! Total energy (helium)
real    (kind=8)              :: ekin4            ! Kinetic energy (helium)
real    (kind=8)              :: ekinx            ! Kinetic energy (impurity)
real    (kind=8)              :: eimpu_impu       ! Impurity-impurity interaction energy (X-X)
real    (kind=8)              :: elj4             ! Lennard-Jones energy
real    (kind=8)              :: ealphas          ! Alpha_s term
real    (kind=8)              :: ecor4            ! Correlation contribution
real    (kind=8)              :: eimpu            ! Impurity energy
real    (kind=8)              :: eHeX             ! Interaction He-X energy
real    (kind=8)              :: esolid             ! Solid term
end module energies
!------------------------------------------------------------------
module field
real    (kind=8), allocatable :: hpsi(:,:,:)     ! H Psi
real    (kind=8), allocatable :: hpsix(:,:,:)     ! H Psi
real    (kind=8), allocatable :: pot4(:,:,:)     ! 'Potential' of lagrange equation
real    (kind=8), allocatable :: uext(:,:,:)     ! 'Potential due to impurities
real    (kind=8), allocatable :: penalty(:,:,:)     ! 'Potential due to impurities
integer (kind=4)              :: iron=1          ! Number of Smoothings
integer (kind=4)              :: ironx=1          ! Number of Smoothings
end module field
!------------------------------------------------------------------
module fftmodule

character (len=15)            :: fftwplan="FFTW_PATIENT"
! real    (kind=8), allocatable :: fin(:,:,:)  ! Work Array for FFT 
! complex (kind=8), allocatable :: fout(:,:,:) ! Work Array for FFT 
! integer (kind=8)              :: pfftfw      ! Pointer for FFT forward
! integer (kind=8)              :: pfftbk      ! Pointer for FFT bakward
integer (kind=8)              :: pfftfw_den  ! Pointer for FFT forward
integer (kind=8)              :: pfftfw_1    ! Pointer for FFT forward
integer (kind=8)              :: pfftfw_2    ! Pointer for FFT forward
integer (kind=8)              :: pfftfw_3    ! Pointer for FFT forward
integer (kind=8)              :: pfftbk_cg   ! Pointer for FFT bakward
integer (kind=8)              :: pfftbk_lj   ! Pointer for FFT bakward
integer (kind=8)              :: pfftbk_as   ! Pointer for FFT bakward
integer (kind=8)              :: pfftbk_ua   ! Pointer for FFT bakward
integer (kind=8)              :: pfftbk_1    ! Pointer for FFT bakward
integer (kind=8)              :: pfftbk_1x   ! Pointer for FFT bakward
integer (kind=8)              :: pfftbk_2y   ! Pointer for FFT bakward
integer (kind=8)              :: pfftbk_3z   ! Pointer for FFT bakward
! Removed by Manu/Francois
!integer (kind=4)              :: nthread=1   ! Number of threads
integer (kind=4)              :: npx,npy,npz ! Number of of points (axis)
real    (kind=8)              :: renor       ! Inverse of (nx*ny*nz)
end module fftmodule
!------------------------------------------------------------------
module grid
real    (kind=8), allocatable :: x(:)            ! Values in X
real    (kind=8), allocatable :: y(:)            !   Id.     Y
real    (kind=8), allocatable :: z(:)            !   Id.     Z
real    (kind=8)              :: xmax=64.0d0     ! Last  point of the grid
real    (kind=8)              :: ymax=64.0d0     ! Last  point of the grid
real    (kind=8)              :: zmax=64.0d0     ! Last  point of the grid
real    (kind=8)              :: xc=0.0d0        ! Center of the cluster. X
real    (kind=8)              :: yc=0.0d0        ! Center of the cluster. Y
real    (kind=8)              :: zc=0.0d0        ! Center of the cluster. Z
real    (kind=8)              :: hx,hy,hz        ! Steps in X, Y and Z
integer (kind=4)              :: nx=64           ! Number of points in X
integer (kind=4)              :: ny=64           ! Number of points in Y
integer (kind=4)              :: nz=64           ! Number of points in Z
integer (kind=4)              :: nxyz            ! nxyz=nx*ny*nz
real    (kind=8)              :: dxyz            ! dxyz=hx*hy*hz
complex (kind=8), allocatable :: timec(:,:,:)   ! time constant, =i for real-time evol, = 1 for minimization
end module grid
!------------------------------------------------------------------
module gridk
real    (kind=8), allocatable :: px(:)          ! Values in PX
real    (kind=8), allocatable :: py(:)          !   Id.     PY
real    (kind=8), allocatable :: pz(:)          !   Id.     PZ
real    (kind=8), allocatable :: pmod(:,:,:)    ! Module of p
real    (kind=8)              :: pmaxx,pmaxy,pmaxz
real    (kind=8)              :: hpx,hpy,hpz
end module gridk
!------------------------------------------------------------------
module he3
real    (kind=8)              :: densat3=0.0163 ! Density of saturation for 3He
!
!.......... Density functional parameter  (Orsay-Trento)
!           (From Barranco et al. PRB 56(1997)8997-9003
!
!..Helium-3  (Refs. 1 and 2)
real      (kind=8) ::    cp3 =1.588790d6   ! K \AA**{3+\gamma_3}
real      (kind=8) ::   cpp3 =    -3.5d4   ! K \AA**{3+\gamma_3}
real      (kind=8) ::   gam3 =  2.1251d0
real      (kind=8) ::  den3c =  0.0406d0   ! \AA**{-3}
real      (kind=8) :: h2o2m3 =  8.041775d0 ! \hbar**2 / (2 m_3)
real      (kind=8) ::   beta =5.55555555555556d-02 !  = 1/18 See Ref.2 p-5026
real      (kind=8) ::  betap =3.33333333333333d-01 !  = 1/3  See Ref.2 p-5026
!real     (kind=8) ::  betap =     0.0d0
end module he3
!------------------------------------------------------------------
module he4

real    (kind=8)              :: densat4=0.02184 ! Density of saturation for 4He
!
!.......... Density functional parameter  (Orsay-Trento)
!           (From Barranco et al. PRB 56(1997)8997-9003
!
!..Helium-4  (Refs. 1 and 2)
real      (kind=8)               ::    cp4=-2.41186d4       ! K \AA**6
real      (kind=8)               ::   cpp4= 1.85850d6       ! K \AA**9
real      (kind=8)               ::  den4c= 0.062d0         ! \AA**{-3}
real      (kind=8)               :: alphas= 54.31d0         ! K ^-1 \AA**3
real      (kind=8)               ::      l= 1.0d0           ! \AA
real      (kind=8)               ::  den0s= 0.04d0          ! \AA**-3
real      (kind=8)               :: h2o2m4= 6.05969638298d0 ! \hbar**2 / (2 m_4)
logical                          :: lsolid=.true. ! Sumem els termes de la funcional solida
real      (kind=8)               :: C = 3.1577504d4
real      (kind=8)               :: beta= 40.d0
real      (kind=8)               :: den_m = 0.37d0
!
end module he4
!------------------------------------------------------------------
module impur
! real      (kind=8)               :: ximp=0.0d0   ! Position of the impurity  X-ccordinate.
! real      (kind=8)               :: yimp=0.0d0   !  "        "  "   "        Y-Coordinate
! real      (kind=8)               :: zimp=0.0d0   !  "        "  "   "        Z-Coordinate
! integer   (kind=4)               :: rimp(3)      ! Coordinates of the impurity.
logical                          :: limp=.false. ! T-> Impurity, F->Pure
real      (kind=8)               :: umax=400.d0  ! Maximum value for the potential
real      (kind=8)               :: r_cutoff=2.0d0 ! Cut_off del potencial
Character (Len=80)               :: selec=''
Integer   (kind=8)               :: nq=1000      ! Number of q-values for Patil
real      (kind=8)               :: tol=1.d-7    ! Tolerance for Romberg
real      (kind=8)               :: rmin=3.d0    ! Core For Patil
complex   (kind=8),allocatable   :: vq(:,:,:)    ! Fourier of Patil Potential
complex   (kind=8),allocatable   :: psix(:,:,:)  ! Wave function for the impurity
complex   (kind=8),allocatable   :: Psixold(:,:,:,:)   ! Old values of Psix
complex   (kind=8),allocatable   :: hpsixold(:,:,:,:)  ! Old values of H�Psix
real      (kind=8),allocatable   :: denx(:,:,:)  ! Density for the impurity
real      (kind=8),allocatable   :: upotx(:,:,:) ! Mean field for the impurity
real      (kind=8),allocatable   :: potx4(:,:,:) !
real      (kind=8),allocatable   :: epsix(:,:,:) ! Enveloping impurity wave function
real      (kind=8),allocatable   :: rmod0(:,:,:) ! Fo the enveloping function
real      (kind=8)               :: x0=0.0d0     ! Position of enveloping  X-ccordinate.
real      (kind=8)               :: y0=0.0d0     !  "        "  "   "      Y-Coordinate
real      (kind=8)               :: z0=0.0d0     !  "        "  "   "      Z-Coordinate
real      (kind=8)               :: Renv=5.d0    !  Parameter for the envelop function
real      (kind=8)               :: Aenv=1.d-2   !  "          "   "      "
real      (kind=8)               :: h2o2mx       ! Term hbar**2/(2*m_x)
real      (kind=8)               :: rinfi=500        ! Infinite for Patil
real      (kind=8)               :: gwf=5.69     ! Parameter for starting gaussian
complex   (kind=8), allocatable  :: fdenx(:,:,:) ! FFT of denx


character (len=3)                :: elem        ! Chemical Symbol for Alkali
character (len=3)                :: leepot='NO '! YES/NO Read VXQ external
character (len=40)               :: vxpot       ! Name of External file with VXQ

end module impur
!------------------------------------------------------------------
module lenard3
real      (kind=8), allocatable  ::   fvlj3(:,:,:)
real      (kind=8)               ::   h3   ! \AA

real      (kind=8)               ::     h3op =2.356415d0   ! \AA
real      (kind=8)               ::     h3ot =2.356415d0   ! \AA
real      (kind=8)               ::   eps3   =   10.22d0   ! K
real      (kind=8)               :: sigma3   =   2.556d0   ! \AA
real      (kind=8)               ::     b3   =-684.676d0   ! K \AA**3
character (len=2)                ::  core3   ='OT'

real      (kind=8), allocatable  :: delj3(:,:,:) ! Density of energy-lennard-Jones
end module lenard3
!------------------------------------------------------------------
module lenard4
real      (kind=8), allocatable  ::   fvlj4(:,:,:)
real      (kind=8)               ::   h4  

real      (kind=8)               ::     h4op =2.359665d0   ! \AA
real      (kind=8)               ::     h4ot =2.190323d0   ! \AA
real      (kind=8)               ::     eps4 =   10.22d0   ! K
real      (kind=8)               ::   sigma4 =   2.556d0   ! \AA
real      (kind=8)               ::       b4 =-718.99d0    ! K \AA**3
character (len=3)                :: core4 ='OT '
real      (kind=8), allocatable  :: delj4(:,:,:) ! Density of energy-lennard-Jones

end module lenard4

!------------------------------------------------------------------
module rho
real    (kind=8), allocatable :: den(:,:,:)      ! Density in Real-space
complex (kind=8), allocatable :: Psi(:,:,:)      ! density = Mod(Psi)**2
complex (kind=8), allocatable :: Psiold(:,:,:,:) ! Old values of Psi
complex (kind=8), allocatable :: hpsiold(:,:,:,:) ! Old values of H�Psi
complex (kind=8), allocatable :: fden(:,:,:)     ! Density in K-space
real    (kind=8), allocatable :: dencg(:,:,:)    ! Coarse-Graining density
real    (kind=8), allocatable :: wcgk(:,:,:)     ! Kernel of coarse graining
                                                 ! in fourier space
real    (kind=8)              :: denmin=1.d-60   ! Minimum value for densities.
real    (kind=8)              :: psimin=1.d-30   ! Minimum value for densities.
Logical                       :: Ldensity=.true.

end module rho

!------------------------------------------------------------------
module util1
character (len=1)  :: cchar="#"  ! Used in routine 'Titols'
real      (kind=8) :: pi         ! Pi=3.141592.... value
real      (kind=8) :: twopi      ! twopi  = 2*pi
real      (kind=8) :: fourpi     ! fourpi = 4*pi
real      (kind=8) :: piq        ! piq    = pi*pi
real      (kind=8) :: afermi     ! Parameter for the initial fermi-distribution
real      (kind=8) :: rfermi     ! Parameter for the initial fermi-distribution
real      (kind=8) :: vdt(2)     ! Speeds for imaginary step-time method.
integer   (kind=4) :: nn(3)      ! Auxiliar array for pderg
integer   (kind=4) :: mmx(4)     ! Auxiliar array for pderg
integer   (kind=4) :: iw(11)     ! Auxiliar array for pderg
integer   (kind=4) :: nsfiles=1 ! Number-of-save-files
integer   (kind=4) :: nsfaux=0   ! Actual generation of backup file
integer   (kind=4) :: nsfaux2=0   ! Actual generation of backup file
integer   (kind=4) :: irespar=0  ! Does not write partial auxiliar plot files...
logical            :: printpot=.false. ! Prints (no) impurity potential

end module util1
!------------------------------------------------------------------
module work1
complex (kind=8), allocatable :: wk1(:,:,:)
complex (kind=8), allocatable :: wk2(:,:,:)
complex (kind=8), allocatable :: wk3(:,:,:)
real    (kind=8), allocatable :: sto1(:,:,:)
real    (kind=8), allocatable :: sto2(:,:,:)
real    (kind=8), allocatable :: sto3(:,:,:)
real    (kind=8), allocatable :: sto4(:,:,:)
real    (kind=8), allocatable :: sto5(:,:,:)
real    (kind=8), allocatable :: sto6(:,:,:)
complex (kind=8), allocatable :: sto1c(:,:,:)
complex (kind=8), allocatable :: sto2c(:,:,:)
complex (kind=8), allocatable :: sto3c(:,:,:)
complex (kind=8), allocatable :: sto4c(:,:,:)
complex (kind=8), allocatable :: sto5c(:,:,:)
complex (kind=8), allocatable :: sto6c(:,:,:)
complex (kind=8), allocatable :: sto7c(:,:,:)
complex (kind=8), allocatable :: sto8c(:,:,:)
end module work1

!------------------------------------------------------------------
module rkpc
integer (kind=4) :: ioldp(3),ioldh(2),ioldpx(3),ioldhx(2)      ! Auxiliar arrays for steppc
complex (kind=8), allocatable ::  q(:,:,:),  qx(:,:,:)         ! Auxiliar arrays for steprk
complex (kind=8), allocatable :: pc(:,:,:), pcx(:,:,:)         ! Auxiliar arrays for steppc
end module rkpc
!------------------------------------------------------------------
module classicimp
real    (kind=8), allocatable :: stor(:,:)
real    (kind=8), allocatable :: rimp(:,:)
real    (kind=8), allocatable :: vimp(:,:)
real    (kind=8), allocatable :: aimp(:,:)
real    (kind=8), allocatable :: F(:,:)
real    (kind=8), allocatable :: F_ij(:,:,:)
real    (kind=8), allocatable :: uimp(:,:,:)  !  Potencial from imp for He
real    (kind=8), allocatable :: uimp_k(:,:,:,:)  !  Potencial from imp_k for He
logical                       :: lselection = .true.
real    (kind=8), allocatable              :: qr(:,:)
real    (kind=8), allocatable              :: qv(:,:)
real    (kind=8), allocatable              :: rimpold(:,:,:)
real    (kind=8), allocatable              :: vimpold(:,:,:)
real    (kind=8), allocatable              :: aimpold(:,:,:)
integer (kind=4)              :: ioldr(3)
integer (kind=4)              :: ioldv(3)
integer (kind=4)              :: iolda(2)
integer (kind=4)              :: N_imp=1 ! Number of impurities. Default is one impurity
real    (kind=8), allocatable              :: pcr(:,:)
real    (kind=8), allocatable              :: pcv(:,:)
real    (kind=8) , parameter   :: mp_u = 0.02061484d0  ! proton mass in Angs**-2 * Kelvin**-1, go figure!
real    (kind=8), allocatable :: m_imp_u(:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
real    (kind=8), allocatable :: m_imp(:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
character (len=80), allocatable :: selec_gs_k(:)
character (len=80), allocatable :: selec_gs_k_k(:,:)
character (len=80), allocatable :: drselec_gs_k_k(:,:)
character (len=80), allocatable :: filerimp_k(:)
character (len=80), allocatable :: filevimp_k(:)
character (len=80), allocatable :: fileaimp_k(:)
real    (kind=8), allocatable :: r_cutoff_gs_k(:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
real    (kind=8), allocatable :: r_cutoff_gs_k_k(:,:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
real    (kind=8), allocatable :: drr_cutoff_gs_k_k(:,:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
real    (kind=8), allocatable :: umax_gs_k(:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
real    (kind=8), allocatable :: umax_gs_k_k(:,:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
real    (kind=8), allocatable :: drumax_gs_k_k(:,:)  ! argon mass in Angs**-2 * Kelvin**-1, go figure!
!


! logical                        :: newpoten = .true.
end module classicimp
!------------------------------------------------------------------
module Aziz
!.......................................!
!... Values of Aziz parameter for Ag ...!
!.......................................!
real      (kind=8) :: A    = 78741.3877
real      (kind=8) :: alpha= 0.65299
real      (kind=8) :: beta = 0.36367
real      (kind=8) :: D    = 12.093
real      (kind=8) :: C6   = 137042.07d0
real      (kind=8) :: C8   = 189386.217d0
real      (kind=8) :: C10  = 111358462.d0
real      (kind=8) :: C12  = 1.17670377d10
! Groundstate
real      (kind=8) :: Ags    = 78741.3877
real      (kind=8) :: alphags= 0.65299
real      (kind=8) :: betags = 0.36367
real      (kind=8) :: Dgs    = 12.093
real      (kind=8) :: C6gs   = 137042.07d0
real      (kind=8) :: C8gs   = 189386.217d0
real      (kind=8) :: C10gs  = 111358462.d0
real      (kind=8) :: C12gs  = 1.17670377d10
! Pi
real      (kind=8) :: Api    = 0.7095092d6
real      (kind=8) :: alphapi= 1.9928
real      (kind=8) :: betapi = 0.45500
real      (kind=8) :: Dpi    = 4.0562
real      (kind=8) :: C6pi   = 0.21661134d6
real      (kind=8) :: C8pi   = 123454.7952
real      (kind=8) :: C10pi  = 0.307155024d6
real      (kind=8) :: C12pi  = 0.740464032d-10
! Sigma
real      (kind=8) :: Asi    = 40005.2
real      (kind=8) :: alphasi= 0.41529
real      (kind=8) :: betasi = 0.14888
real      (kind=8) :: Dsi    = 25.730
real      (kind=8) :: C6si   = 1.97443d-9
real      (kind=8) :: C8si   = 3.37393d-6
real      (kind=8) :: C10si  = 0.00788363
real      (kind=8) :: C12si  = 7.52624d12
end module Aziz
!------------------------------------------------------------------
module interpol
real      (kind=8) :: DelInter
real      (kind=8) , allocatable :: potion(:,:) 
integer   (kind=4) :: npot
!real      (kind=8) , allocatable :: Vion(:,:,:) 
!real      (kind=8) :: lastr(3)
end module interpol
!------------------------------------------------------------------
