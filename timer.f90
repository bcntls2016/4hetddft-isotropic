!--------------------------------------------------------------------
!--                      Subroutine TIMER                         ---
!--------------------------------------------------------------------
!
! This subroutine gives the time in seconds (REAL*4) measured from
! the begining of the run
! Should be adapted for each operating system.
!
! Version for Intel compiler. for Other UNIX* see manual
!
subroutine timer(secs)
use ifport
!real :: etime
real (kind=4) :: secs
real (kind=4) :: ta(2)
secs=etime(ta)
return
end
