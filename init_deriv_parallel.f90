Module Para_DerivnD
Interface DerivnD
  Subroutine Deriv1D_p(id,n,h,f,df,icon)
    Implicit none
    Integer (kind=4), Intent (in) :: id, n, icon
    Real    (Kind=8), Intent (in) :: h
    Real    (Kind=8), Intent (in) :: f(n)
    Real    (Kind=8), Intent(out) :: df(n)
  End Subroutine Deriv1D_p
  Subroutine Deriv2D_p(id,n,h,iv,f,df,icon)
    Implicit none
    Integer (kind=4), Intent(in)  :: id, n(2), icon, iv
    Real    (Kind=8), Intent(in)  :: h
    Real    (Kind=8), Intent(in)  :: f(n(1),n(2))
    Real    (Kind=8), Intent(out) :: df(n(1),n(2))
  End Subroutine Deriv2D_p
  Subroutine Deriv3D_p(id,n,h,iv,f,df,icon)
    Implicit none
    Integer (kind=4), Intent(in)  :: id, n(3), icon, iv
    Real    (Kind=8), Intent(in)  :: h
    Real    (Kind=8), Intent(in)  :: f(n(1),n(2),n(3))
    Real    (Kind=8), Intent(out) :: df(n(1),n(2),n(3))
  End Subroutine Deriv3D_p
  Subroutine Deriv1D_pc(id,n,h,f,df,icon)
    Implicit none
    Integer (kind=4), Intent(in)  :: id, n, icon
    Real    (Kind=8), Intent(in)  :: h
    Complex (Kind=8), Intent(in)  :: f(n)
    Complex (Kind=8), Intent(out) :: df(n)
  End Subroutine Deriv1D_pc
  Subroutine Deriv2D_pc(id,n,h,iv,f,df,icon)
    Implicit none
    Integer (kind=4), Intent(in)  :: id, n(2), icon, iv
    Real    (Kind=8), Intent(in)  :: h
    Complex (Kind=8), Intent(in)  :: f(n(1),n(2))
    Complex (Kind=8), Intent(out) :: df(n(1),n(2))
  End Subroutine Deriv2D_pc
  Subroutine Deriv3D_pc(id,n,h,iv,f,df,icon)
    Implicit none
    Integer (kind=4), Intent(in)  :: id, n(3), icon, iv
    Real    (Kind=8), Intent(in)  :: h
    Complex (Kind=8), Intent(in)  :: f(n(1),n(2),n(3))
    Complex (Kind=8), Intent(out) :: df(n(1),n(2),n(3))
  End Subroutine Deriv3D_pc
End Interface
End Module Para_DerivnD
Module Deriv_p
  Real    (Kind=8), Allocatable :: Cdc(:,:,:)
  Integer (Kind=4)              :: npd, ndmax
End Module Deriv_p
SubRoutine Init_deriv_p(k,kmax,Number_of_Threads)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)    :: k,kmax
Integer (Kind=4), Intent(in)    :: Number_of_Threads
Integer (Kind=4), Save   :: Icon_save=-1
Logical (Kind=4)   :: OMP_Dynamic_Enable=.false.
Real (Kind=8), Allocatable      :: f(:), Caux(:)
Real (Kind=8) ::tmp_sum

If(k.Lt.kmax)Then
  Write(0,'(" From Init_deriv_p: k no pot ser mes petit que kmax",2i5)')k,kmax
  Stop '001'
Endif
npd=k
ndmax=kmax
kkp1=k*k+1
If(Allocated(Cdc))Deallocate(Cdc)
Allocate(Cdc(kmax,k,k))
Allocate(Caux(kkp1))
Allocate(f(kkp1))

!$ Call Omp_Set_Dynamic(OMP_Dynamic_Enable)
!$ Call Omp_Set_Num_Threads(Number_of_Threads)

Do id=1,kmax
  call newder(id,K,kkp1,h,f,Caux,1)
  ForAll(i=1:k,j=1:k)
      cdc(id,i,j)=Caux(j+(i-1)*k)
  EndForAll
Enddo
Deallocate(f)
Deallocate(Caux)
Return
End
!
!   El parametre Icon es el que controla les condicions de contorn de totes 
!   les subroutines, per veure com s'ha de definir aneu a veure la capsalera
!   de la Subroutine Newder i/o Pderg
!
Subroutine Deriv1D_p(id,n,h,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,n,icon     ! Ordre de la derivada, nombre de punts, condicions de contorn, dummy argument
Real    (Kind=8), Intent(in)  :: h             ! Interval en la definicio de la funcio
Real    (Kind=8), Intent(in)  :: f(n)          ! Funcio (entrada)
Real    (Kind=8), Intent(out) :: df(n)         ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
NMKO2=N-KO2
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv1D_p: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(n.Lt.npd)Then
  If(n.Ge.3)Write(*,'("From Deriv1D_p: n ha de ser mes gran que npd",2i5)')n,npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon

  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  !$Omp Parallel
  !$OMP Do
     Do i=1,Ko2
       df(i)=Sum(Cd(i,:)*f(1:k))
      End Do
  !$OMP End Do Nowait

  !$OMP Do
    Do i=ko2p1,NmKo2
      df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
    EndDo
  !$OMP End Do Nowait

  !$OMP Do
    Do i=Nmko2+1,N
      df(i)=Sum(Cd(i-N+k,:)*f(N-k+1:N))
    EndDo
  !$OMP End Do
  !$Omp End Parallel
Else
  is=1
  If(Icon.Eq.18)is=-1

  !$Omp Parallel
  !$OMP Do
     Do i=1,Ko2
       df(i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:N))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i))
      End Do
  !$OMP End Do Nowait

  !$OMP Do
    Do i=ko2p1,NmKo2
      df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
    EndDo
  !$OMP End Do Nowait

  !$OMP Do
    Do i=Nmko2+1,N
      df(i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:N))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2))
    EndDo
  !$OMP End Do
  !$Omp End Parallel
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 2D
!
Subroutine Deriv2D_p(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer (Kind=4), Intent(in)  :: iv,nn(2) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real    (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Real    (Kind=8), Intent(in)  ::  f(nn(1),nn(2)) ! Funcio (entrada)
Real    (Kind=8), Intent(out) :: df(nn(1),nn(2)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv2D_p: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv2D_p: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon        

  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
!$Omp Parallel 
!$OMP Do 
    Do i2=1, n2
      Do i=1,Ko2
        df(i,i2)=Sum(Cd(i,:)*f(1:k,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=ko2p1,NmKo2
        df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=Nmko2+1,n1
        df(i,i2)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(i,:)*f(i1,1:k))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv2D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Else

  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2

!$Omp Parallel 
!$OMP Do 
    Do i2=1, n2
      Do i=1,Ko2
        df(i,i2)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=ko2p1,NmKo2
        df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=Nmko2+1,n1
        df(i,i2)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i1=1, n1
        df(i1,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv2D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 3D
!
Subroutine Deriv3D_p(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer (Kind=4), Intent(in)  :: iv,nn(3) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real    (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Real    (Kind=8), Intent(in)  ::  f(nn(1),nn(2),nn(3)) ! Funcio (entrada)
Real    (Kind=8), Intent(out) :: df(nn(1),nn(2),nn(3)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)
integer ik

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2); n3=nn(3)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv3D_p: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv3D_p: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon
  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then

  If(iv.Eq.1)Then
    NMKO2=n1-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=1,Ko2
          df(i,i2,i3)=Sum(Cd(i,:)*f(1:k,i2,i3))
        EndDo
    EndDo
  EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=ko2p1,NmKo2
          df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=Nmko2+1,n1
          df(i,i2,i3)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i=1,Ko2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(i,:)*f(i1,1:k,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i=ko2p1,NmKo2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i=Nmko2+1,n2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(i,:)*f(i1,i2,1:k))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i2=1, n2
        Do i1=1, n1
	  tmp_sum=0.d0
	  !DIR$ VECTOR ALWAYS
	  !DIR$ DISTRIBUTE POINT 
!         do ik=0,2*ko2
!                tmp_sum=tmp_sum+Cd(ko2p1,ik+1)*f(i1,i2,i-ko2+ik)
!         enddo   
!         df(i1,i2,i)=tmp_sum
          df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n3
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(i-n3+k,:)*f(i1,i2,n3-k+1:n3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv3D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Else

  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=1,Ko2
          df(i,i2,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2,i3))
        EndDo
    EndDo
  EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=ko2p1,NmKo2
          df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=Nmko2+1,n1
          df(i,i2,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i=1,Ko2
        Do i1=1, n1
          df(i1,i,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait
!$OMP Do 

    Do i3=1, n3
      Do i=ko2p1,NmKo2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i=Nmko2+1,n2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,i2,Nmko2+i:n3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,i2,1:ko2+i))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n3
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i2,i-ko2:n3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,i2,1:i-Nmko2))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv3D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
Subroutine Deriv1D_pc(id,n,h,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,n,icon     ! Ordre de la derivada, nombre de punts, condicions de contorn, dummy argument
Real    (Kind=8), Intent(in)  :: h             ! Interval en la definicio de la funcio
Complex    (Kind=8), Intent(in)  :: f(n)       ! Funcio (entrada)
Complex    (Kind=8), Intent(out) :: df(n)      ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
NMKO2=N-KO2
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv1D_p: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(n.Lt.npd)Then
  If(n.Ge.3)Write(*,'("From Deriv1D_p: n ha de ser mes gran que npd",2i5)')n,npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon

  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  !$Omp Parallel
  !$OMP Do
     Do i=1,Ko2
       df(i)=Sum(Cd(i,:)*f(1:k))
      End Do
  !$OMP End Do Nowait

  !$OMP Do
    Do i=ko2p1,NmKo2
      df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
    EndDo
  !$OMP End Do Nowait

  !$OMP Do
    Do i=Nmko2+1,N
      df(i)=Sum(Cd(i-N+k,:)*f(N-k+1:N))
    EndDo
  !$OMP End Do
  !$Omp End Parallel
Else
  is=1
  If(Icon.Eq.18)is=-1

  !$Omp Parallel
  !$OMP Do
     Do i=1,Ko2
       df(i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:N))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i))
      End Do
  !$OMP End Do Nowait

  !$OMP Do
    Do i=ko2p1,NmKo2
      df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
    EndDo
  !$OMP End Do Nowait

  !$OMP Do
    Do i=Nmko2+1,N
      df(i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:N))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2))
    EndDo
  !$OMP End Do
  !$Omp End Parallel
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 2D
!
Subroutine Deriv2D_pc(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer (Kind=4), Intent(in)  :: iv,nn(2) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real    (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Complex    (Kind=8), Intent(in)  ::  f(nn(1),nn(2)) ! Funcio (entrada)
Complex    (Kind=8), Intent(out) :: df(nn(1),nn(2)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv2D_p: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv2D_p: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon

  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
!$Omp Parallel 
!$OMP Do 
    Do i2=1, n2
      Do i=1,Ko2
        df(i,i2)=Sum(Cd(i,:)*f(1:k,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=ko2p1,NmKo2
        df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=Nmko2+1,n1
        df(i,i2)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(i,:)*f(i1,1:k))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv2D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Else

  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2

!$Omp Parallel 
!$OMP Do 
    Do i2=1, n2
      Do i=1,Ko2
        df(i,i2)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=ko2p1,NmKo2
        df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i2=1, n2
      Do i=Nmko2+1,n1
        df(i,i2)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i1=1, n1
        df(i1,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n2
      Do i1=1, n1
        df(i1,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2))
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv2D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 3D
!
Subroutine Deriv3D_pc(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer (Kind=4), Intent(in)  :: iv,nn(3) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real    (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Complex    (Kind=8), Intent(in)  ::  f(nn(1),nn(2),nn(3)) ! Funcio (entrada)
Complex    (Kind=8), Intent(out) :: df(nn(1),nn(2),nn(3)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2); n3=nn(3)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv3D_p: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv3D_p: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon

  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then

  If(iv.Eq.1)Then
    NMKO2=n1-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=1,Ko2
          df(i,i2,i3)=Sum(Cd(i,:)*f(1:k,i2,i3))
        EndDo
    EndDo
  EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=ko2p1,NmKo2
          df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=Nmko2+1,n1
          df(i,i2,i3)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i=1,Ko2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(i,:)*f(i1,1:k,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i=ko2p1,NmKo2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i=Nmko2+1,n2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(i,:)*f(i1,i2,1:k))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n3
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(i-n3+k,:)*f(i1,i2,n3-k+1:n3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv3D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Else

  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=1,Ko2
          df(i,i2,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2,i3))
        EndDo
    EndDo
  EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=ko2p1,NmKo2
          df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i2=1, n2
        Do i=Nmko2+1,n1
          df(i,i2,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2

!$Omp Parallel 
!$OMP Do 
    Do i3=1, n3
      Do i=1,Ko2
        Do i1=1, n1
          df(i1,i,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait
!$OMP Do 

    Do i3=1, n3
      Do i=ko2p1,NmKo2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i3=1, n3
      Do i=Nmko2+1,n2
        Do i1=1, n1
          df(i1,i,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2,i3))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2

!$Omp Parallel 
!$OMP Do 
    Do i=1,Ko2
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,i2,Nmko2+i:n3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,i2,1:ko2+i))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=ko2p1,NmKo2
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
        EndDo
      EndDo
    EndDo
!$OMP End Do Nowait

!$OMP Do 
    Do i=Nmko2+1,n3
      Do i2=1, n2
        Do i1=1, n1
          df(i1,i2,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i2,i-ko2:n3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,i2,1:i-Nmko2))
        EndDo
      EndDo
    EndDo
!$OMP End Do
!$OMP End Parallel

  Else
    Write(*,'("From Deriv3D_p: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
Subroutine Deriv1D(id,n,h,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,n,icon  ! Ordre de la derivada, nombre de punts, condicions de contorn, dummy argument
Real    (Kind=8), Intent(in)  :: h             ! Interval en la definicio de la funcio
Real    (Kind=8), Intent(in)  :: f(n)          ! Funcio (entrada)
Real    (Kind=8), Intent(out) :: df(n)         ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
NMKO2=N-KO2
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv1D: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(n.Lt.npd)Then
  If(n.Ge.3)Write(*,'("From Deriv1D: n ha de ser mes gran que npd",2i5)')n,npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon

  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
!OMP Parallel Shared(Df,Cd,f)
!OMP Whorkshare
  ForAll (i=1:Ko2)
     df(i)=Sum(Cd(i,:)*f(1:k))
  End ForAll
!OMP End Whorkshare Nowait

!OMP Whorkshare
  ForAll(i=ko2p1:Nmko2)
    df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
  End ForAll
!OMP End Whorkshare Nowait

!OMP Whorkshare
  ForAll(i=Nmko2+1:n)
    df(i)=Sum(Cd(i-N+k,:)*f(N-k+1:N))
  End ForAll
!OMP End Whorkshare
!OMP End Parallel
Else
  is=1
  If(Icon.Eq.18)is=-1
!OMP Parallel Shared(Df,Cd,f,is)
!OMP Whorkshare
  ForAll (i=1:Ko2)
     df(i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:N))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i))
  End ForAll
!OMP End Whorkshare Nowait

!OMP Whorkshare
  ForAll(i=ko2p1:Nmko2)
    df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
  End ForAll
!OMP End Whorkshare Nowait

!OMP Whorkshare
  ForAll(i=Nmko2+1:n)
    df(i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:N))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2))
  End ForAll
!OMP End Whorkshare
!OMP End Parallel
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 2D
!
Subroutine Deriv2D(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer (Kind=4), Intent(in)  :: iv,nn(2) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real    (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Real    (Kind=8), Intent(in)  ::  f(nn(1),nn(2)) ! Funcio (entrada)
Real    (Kind=8), Intent(out) :: df(nn(1),nn(2)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv2D: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv2D: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon
  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll
  Call Redef(npd,Caux,Icon)
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)
Else
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll
Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i2=1:n2,i=1:Ko2)
      df(i,i2)=Sum(Cd(i,:)*f(1:k,i2))
    End ForAll
    ForAll(i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
    End ForAll
    ForAll(i2=1:n2,i=Nmko2+1:n1)
      df(i,i2)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i1=1:n1,i=1:Ko2)
      df(i1,i)=Sum(Cd(i,:)*f(i1,1:k))
    End ForAll
    ForAll(i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
    End ForAll
    ForAll(i1=1:n1,i=NmKo2+1:n2)
      df(i1,i)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2))
    End ForAll
  Else
    Write(*,'("From Deriv2D: el parametre iv",i5,"no es valid")')iv
  Endif
Else
  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i2=1:n2,i=1:Ko2)
      df(i,i2)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2))
    End ForAll
    ForAll(i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
    End ForAll
    ForAll(i2=1:n2,i=Nmko2+1:n1)
      df(i,i2)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i1=1:n1,i=1:Ko2)
      df(i1,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i))
    End ForAll
    ForAll(i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
    End ForAll
    ForAll(i1=1:n1,i=NmKo2+1:n2)
      df(i1,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2))
    End ForAll
  Else
    Write(*,'("From Deriv2D: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 3D
!
Subroutine Deriv3D(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer (Kind=4), Intent(in)  :: iv,nn(3) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real    (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Real    (Kind=8), Intent(in)  ::  f(nn(1),nn(2),nn(3)) ! Funcio (entrada)
Real    (Kind=8), Intent(out) :: df(nn(1),nn(2),nn(3)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2); n3=nn(3)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv3D: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv3D: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon
  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll
  Call Redef(npd,Caux,Icon)
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)
Else
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll
Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i3=1:n3,i2=1:n2,i=1:Ko2)
      df(i,i2,i3)=Sum(Cd(i,:)*f(1:k,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=Nmko2+1:n1)
      df(i,i2,i3)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2,i3))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i3=1:n3,i1=1:n1,i=1:Ko2)
      df(i1,i,i3)=Sum(Cd(i,:)*f(i1,1:k,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=Nmko2+1:n2)
      df(i1,i,i3)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2,i3))
    End ForAll
  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2
    ForAll(i2=1:n2,i1=1:n1,i=1:Ko2)
      df(i1,i2,i)=Sum(Cd(i,:)*f(i1,i2,1:k))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=Nmko2+1:n3)
      df(i1,i2,i)=Sum(Cd(i-n3+k,:)*f(i1,i2,n3-k+1:n3))
    End ForAll
  Else
    Write(*,'("From Deriv3D: el parametre iv",i5,"no es valid")')iv
  Endif
Else
  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i3=1:n3,i2=1:n2,i=1:Ko2)
      df(i,i2,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=Nmko2+1:n1)
      df(i,i2,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2,i3))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i3=1:n3,i1=1:n1,i=1:Ko2)
      df(i1,i,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=Nmko2+1:n2)
      df(i1,i,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2,i3))
    End ForAll
  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2
    ForAll(i2=1:n2,i1=1:n1,i=1:Ko2)
      df(i1,i2,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,i2,Nmko2+i:n3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,i2,1:ko2+i))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=Nmko2+1:n3)
      df(i1,i2,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i2,i-ko2:n3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,i2,1:i-Nmko2))
    End ForAll
  Else
    Write(*,'("From Deriv3D: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
!
! Derivades de funcions comlexes
!
Subroutine Deriv1Dc(id,n,h,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer (Kind=4), Intent(in)  :: id,n,icon  ! Ordre de la derivada, nombre de punts, condicions de contorn, dummy argument
Real       (Kind=8), Intent(in)  :: h          ! Interval en la definicio de la funcio
Complex    (Kind=8), Intent(in)  :: f(n)       ! Funcio (entrada)
Complex    (Kind=8), Intent(out) :: df(n)      ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
NMKO2=N-KO2
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv1D: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon

  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll

  Call Redef(npd,Caux,Icon)

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)

Else

  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll

Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  ForAll (i=1:Ko2)
     df(i)=Sum(Cd(i,:)*f(1:k))
  End ForAll

  ForAll(i=ko2p1:Nmko2)
    df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
  End ForAll

  ForAll(i=Nmko2+1:n)
    df(i)=Sum(Cd(i-N+k,:)*f(N-k+1:N))
  End ForAll
Else
  is=1
  If(Icon.Eq.18)is=-1
  ForAll (i=1:Ko2)
     df(i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:N))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i))
  End ForAll

  ForAll(i=ko2p1:Nmko2)
    df(i)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2))
  End ForAll

  ForAll(i=Nmko2+1:n)
    df(i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:N))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2))
  End ForAll
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 2D de funcions comlexes
!
Subroutine Deriv2Dc(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer    (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer    (Kind=4), Intent(in)  :: iv,nn(2) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real       (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Complex    (Kind=8), Intent(in)  ::  f(nn(1),nn(2)) ! Funcio (entrada)
Complex    (Kind=8), Intent(out) :: df(nn(1),nn(2)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv2D: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv2Dc: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon
  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll
  Call Redef(npd,Caux,Icon)
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)
Else
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll
Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i2=1:n2,i=1:Ko2)
      df(i,i2)=Sum(Cd(i,:)*f(1:k,i2))
    End ForAll
    ForAll(i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
    End ForAll
    ForAll(i2=1:n2,i=Nmko2+1:n1)
      df(i,i2)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i1=1:n1,i=1:Ko2)
      df(i1,i)=Sum(Cd(i,:)*f(i1,1:k))
    End ForAll
    ForAll(i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
    End ForAll
    ForAll(i1=1:n1,i=NmKo2+1:n2)
      df(i1,i)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2))
    End ForAll
  Else
    Write(*,'("From Deriv2D: el parametre iv",i5,"no es valid")')iv
  Endif
Else
  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i2=1:n2,i=1:Ko2)
      df(i,i2)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2))
    End ForAll
    ForAll(i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2))
    End ForAll
    ForAll(i2=1:n2,i=Nmko2+1:n1)
      df(i,i2)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i1=1:n1,i=1:Ko2)
      df(i1,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i))
    End ForAll
    ForAll(i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2))
    End ForAll
    ForAll(i1=1:n1,i=NmKo2+1:n2)
      df(i1,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2))
    End ForAll
  Else
    Write(*,'("From Deriv2D: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
!
!   Derivades 3D de funcions comlexes
!
Subroutine Deriv3Dc(id,nn,h,iv,f,df,icon)
Use Deriv_P
Implicit Real*8(A-H,O-Z)
Integer    (Kind=4), Intent(in)  :: id,icon  ! Ordre de la derivada i condicions de contorn
Integer    (Kind=4), Intent(in)  :: iv,nn(3) ! Variable de la que calculem la derivada, nombre de punts en cada direccio
Real       (Kind=8), Intent(in)  :: h        ! Interval en la definicio de la funcio
Complex    (Kind=8), Intent(in)  ::  f(nn(1),nn(2),nn(3)) ! Funcio (entrada)
complex    (Kind=8), Intent(out) :: df(nn(1),nn(2),nn(3)) ! Derivada (sortida)
Real    (Kind=8), Allocatable :: Cd(:,:),Caux(:)

K=npd; Kk=K*K; Kkp1=KK+1
KO2=K/2; KO2P1=KO2+1; K1=K-1; Ko2p2=ko2+2
n1=nn(1); n2=nn(2); n3=nn(3)
Allocate(Cd(k,k))

If(id.Gt.ndmax)Then
  Write(*,'("From Deriv3D: Has de aumentar kmax",2i5)')id,kmax
  return
Endif

If(nn(iv).Lt.npd)Then
  If(nn(iv).Ge.3)Write(*,'("From Deriv3Dc: n ha de ser mes gran que npd",2i5)')nn(iv),npd
  df=0.d0
  return
Endif

If(Icon.Ne.0.And.Icon.Ne.8.And.Icon.Ne.18.And.Icon_save.Ne.Icon)Then
  Icon_save=Icon
  Allocate(Caux(kkp1))
  ForAll(i=1:k,j=1:k)
    Caux(j+(i-1)*npd)=Cdc(id,i,j)
  End ForAll
  Call Redef(npd,Caux,Icon)
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Caux(j+(i-1)*npd)*Sto
  End ForAll
  DeAllocate(Caux)
Else
  Sto=1.0d0/h**id
  ForAll(i=1:k,j=1:k)
    Cd(i,j)=Cdc(id,i,j)*Sto
  End ForAll
Endif
If(Icon.Ne.8.And.Icon.Ne.18)Then
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i3=1:n3,i2=1:n2,i=1:Ko2)
      df(i,i2,i3)=Sum(Cd(i,:)*f(1:k,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=Nmko2+1:n1)
      df(i,i2,i3)=Sum(Cd(i-n1+k,:)*f(n1-k+1:n1,i2,i3))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i3=1:n3,i1=1:n1,i=1:Ko2)
      df(i1,i,i3)=Sum(Cd(i,:)*f(i1,1:k,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=Nmko2+1:n2)
      df(i1,i,i3)=Sum(Cd(i-n2+k,:)*f(i1,n2-k+1:n2,i3))
    End ForAll
  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2
    ForAll(i2=1:n2,i1=1:n1,i=1:Ko2)
      df(i1,i2,i)=Sum(Cd(i,:)*f(i1,i2,1:k))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=Nmko2+1:n3)
      df(i1,i2,i)=Sum(Cd(i-n3+k,:)*f(i1,i2,n3-k+1:n3))
    End ForAll
  Else
    Write(*,'("From Deriv3D: el parametre iv",i5,"no es valid")')iv
  Endif
Else
  is=1
  If(Icon.Eq.18)is=-1
  If(iv.Eq.1)Then
    NMKO2=n1-KO2
    ForAll(i3=1:n3,i2=1:n2,i=1:Ko2)
      df(i,i2,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(Nmko2+i:n1,i2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(1:ko2+i,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=ko2p1:NmKo2)
      df(i,i2,i3)=Sum(Cd(ko2p1,:)*f(i-ko2:i+ko2,i2,i3))
    End ForAll
    ForAll(i3=1:n3,i2=1:n2,i=Nmko2+1:n1)
      df(i,i2,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i-ko2:n1,i2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(1:i-Nmko2,i2,i3))
    End ForAll
  ElseIf(iv.Eq.2)Then
    NMKO2=n2-KO2
    ForAll(i3=1:n3,i1=1:n1,i=1:Ko2)
      df(i1,i,i3)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,Nmko2+i:n2,i3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,1:ko2+i,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i,i3)=Sum(Cd(ko2p1,:)*f(i1,i-ko2:i+ko2,i3))
    End ForAll
    ForAll(i3=1:n3,i1=1:n1,i=Nmko2+1:n2)
      df(i1,i,i3)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i-ko2:n2,i3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,1:i-Nmko2,i3))
    End ForAll
  ElseIf(iv.Eq.3)Then
    NMKO2=n3-KO2
    ForAll(i2=1:n2,i1=1:n1,i=1:Ko2)
      df(i1,i2,i)=is*Sum(Cd(ko2p1,1:ko2p1-i)*f(i1,i2,Nmko2+i:n3))+sum(Cd(ko2p1,ko2p2-i:k)*f(i1,i2,1:ko2+i))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=ko2p1:NmKo2)
      df(i1,i2,i)=Sum(Cd(ko2p1,:)*f(i1,i2,i-ko2:i+ko2))
    End ForAll
    ForAll(i2=1:n2,i1=1:n1,i=Nmko2+1:n3)
      df(i1,i2,i)=Sum(Cd(ko2p1,1:K-i+Nmko2)*f(i1,i2,i-ko2:n3))+is*sum(Cd(ko2p1,K-i+Nmko2+1:k)*f(i1,i2,1:i-Nmko2))
    End ForAll
  Else
    Write(*,'("From Deriv3D: el parametre iv",i5,"no es valid")')iv
  Endif
Endif
DeAllocate(Cd)
Return
End
Subroutine dhms(id,ih,im,s,ic,icd)
Implicit None
Integer (Kind=4), Intent(in)     ::  ic,icd
Integer (Kind=4), Intent(out)    ::  id,ih,im
Real    (Kind=8), Intent(out) ::  s
s=Dfloat(ic)/Dfloat(icd)
id=s/86400.
s=s-id*86400.
ih=s/3600.
s=s-ih*3600.
im=s/60.
s=s-im*60.
Return
End
