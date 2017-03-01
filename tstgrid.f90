integer function tstgrid(n)

implicit none

integer (kind=4), intent(in)  :: n
integer (kind=4)              :: i

integer (kind=4), parameter, dimension(53) ::                         &
        nvalid=(/1,2,3,4,5,6,8,9,10,12,15,16,18,20,24,25,27,30,32,36, &
                40,45,48,50,54,60,64,72,75,80,81,90,96,100,108,120,   &
               125,128,135,144,150,160,162,180,192,200,216,225,240,   &
               243,250,256,320/)

do i=53,1,-1
  if(n.eq.nvalid(i)) then
       tstgrid = 0
       return
  end if
end do
tstgrid = 1
return
end function tstgrid
