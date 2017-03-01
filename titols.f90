!----------------------------------------------------------------------
!--                 Subroutine TITOLS                               ---
!----------------------------------------------------------------------
!
!   Esta rutina posiciona el puntero de lectura de un fichero
!   hasta la primera linea cuya primera columna NO comience por el
!   caracter  'cchar' y devuelve en la variable isalto el numero de
!   lineas que se ha saltado.
!
!   La variable ulog indica el numero de unidad logica a utilizar.
!
subroutine titols(ulog,cchar,isalto)
implicit none

character (LEN=1) :: pchar
character (LEN=1) :: pcolumn
character (LEN=1) :: cchar

integer :: nl
integer :: i
integer :: ulog
integer :: isalto

nl    = 0
pcolumn=cchar
do while(pcolumn.eq.cchar)
  read(ulog,5000,end=9000) pchar
  pcolumn=pchar
  nl=nl+1
end do
rewind(ulog)
isalto=nl-1
do i=1,isalto
  read(ulog,*)
end do

return

9000 print *,' Ey Mister this file is too short ...'  
     stop 'STOP due to severe errors in routine TITOLS'
5000 format(A1)
end
