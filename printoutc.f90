!------------------------------------------------------------------
!---                    Subroutine printout                     ---
!------------------------------------------------------------------

!    iprint = 1 ---> Last printout before exit.
!    iprint = 2 ---> Print due a change of Paflov parameter.
!    iprint = 3 ---> Security copy.
!
!  If we use negative values the output files will be written
!  in binary format.


!subroutine printoutc(time,iprint,namefile,namefile1,             &
!                    psi,nx,ny,nz,hx,hy,hz,       &
!                    xmax,ymax,zmax,rimp,vimp,psix,    &
!                    paflv,iter)
subroutine printoutc(pr)
!subroutine printoutc(namefile,nx,ny,nz,psi,rimp,vimp,pr)

Use Impressio

implicit none

integer (kind=4) :: k,m,N_imp

Type(info_printout) pr
  N_imp=pr%N_imp
  open(10,file=pr%namefile,form='FORMATTED',BUFFERED='yes')
    write(10,'("#  ")')
    write(10,'("#  Density after ",I15," iterations")') pr%it
    write(10,'("#  Actual deltatps....:",1p,E20.10)') pr%dtps
    Write(10,'("#  Total evolution time(ps)...:",1p,E20.12)')pr%time
	do k=1,N_imp
      Write(10,'("#  Impurity-helium potential.........:  ",A)')pr%selec_gs_k(k)
      Write(10,'("#  r_cutoff_gs:",1p,E16.8)')pr%r_cutoff_gs_k(k)
      Write(10,'("#  umax_gs....:",1p,E16.8)')pr%umax_gs_k(k)
	  do m=k+1,N_imp
        Write(10,'("#  Impurity-impurity potential.........:  ",A)')pr%selec_gs_k_k(k,m)
        Write(10,'("#  r_cutoff_gs:",1p,E16.8)')pr%r_cutoff_gs_k_k(k,m)
        Write(10,'("#  umax_gs....:",1p,E16.8)')pr%umax_gs_k_k(k,m)
        Write(10,'("#  Derv. Impurity-impurity potential.........:  ",A)')pr%drselec_gs_k_k(k,m)
        Write(10,'("#  r_cutoff_gs:",1p,E16.8)')pr%drr_cutoff_gs_k_k(k,m)
        Write(10,'("#  umax_gs....:",1p,E16.8)')pr%drumax_gs_k_k(k,m)
	  enddo
	enddo
    Write(10,'("#  xcm,ycm,zcm:",1p,3E16.8)')pr%cm
    Write(10,'("#  NParticles.:",1p,E16.8)')pr%auxn4
    Write(10,'("#  Ekin.......:",1p,E16.8)')pr%ekin
    Write(10,'("#  Elj........:",1p,E16.8)')pr%elj
    Write(10,'("#  Ealphas....:",1p,E16.8)')pr%ealphas
    Write(10,'("#  Ecorr......:",1p,E16.8)')pr%ecor
    Write(10,'("#  Esolid.....:",1p,E16.8)')pr%esolid
    Write(10,'("#  Ekinx......:",1p,E16.8)')pr%ekinx
    Write(10,'("#  Eimpu......:",1p,E16.8)')pr%evx
    Write(10,'("#  Etot.......:",1p,E16.8)')pr%etot
    Write(10,'("#  <Lx,y,z>...:",1p,3E16.8)')pr%ang
    Write(10,'("#  <x2,y2,z2>.:",1p,3E16.8)')pr%r2
    write(10,'("#  ")')
    write(10,*) pr%xmax,pr%ymax,pr%zmax,pr%hx,pr%hy,pr%hz,pr%nx,pr%ny,pr%nz
    write(10,*) pr%N_imp
    write(10,*) pr%rimp
    write(10,*) pr%vimp
    write(10,*) pr%psi
  close(10)


return

end
