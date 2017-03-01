129c129,130
< character  (len=3)   :: chariter
---
> !character  (len=3)   :: chariter
> character  (len=4)   :: chariter


135c136,138
< real       (kind=8)  :: Lambdah,tzmean,tzsurf,rt
---
> real       (kind=8)  :: Lambdah=2.d0,txmean=1000.d0,txsurf=2.d0      &
>                                     ,tymean=1000.d0,tysurf=2.d0      &
>                                     ,tzmean=1000.d0,tzsurf=2.d0


164c167,168
<                 pcurr,icurr,lsurf,lsurf3D,Lambdah,tzmean,tzsurf,        &
---
>                 pcurr,icurr,lsurf,lsurf3D,Lambdah,                      &
>                 txmean,txsurf,tymean,tysurf,tzmean,tzsurf,              &


415a420,421
>    case(4) !................................... Continue a dynamic calculation from impurity with excited internal state
>       write(6,6013) filedenout,fileimpout


691,694c697,701
< !    rt = dsqrt(x(ix)*x(ix)+y(iy)*y(iy)+z(iz)*z(iz))
<    rt = dsqrt(x(ix)*x(ix)+y(iy)*y(iy))
< !    timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((rt-tmean)/tsurf)),1.d0)
<    timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((abs(z(iz))-tzmean)/tzsurf)),1.d0)
---
> !   timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((abs(z(iz))-tzmean)/tzsurf)),1.d0)
>    timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((abs(x(ix))-txmean)/txsurf)            &
>                                  +1.d0+tanh((abs(y(iy))-tymean)/tysurf)            &
>                                  +1.d0+tanh((abs(z(iz))-tzmean)/tzsurf))           &
>                                                                         ,1.d0)

966c973,974
<       write(chariter,8011)ncurr
---
> !      write(chariter,8011)ncurr
>          write(chariter,'("000",I1)')ncurr


968c976,977
<       write(chariter,8012)ncurr
---
> !      write(chariter,8012)ncurr
>          write(chariter,'("00",I2)')ncurr


970c979,982
<       write(chariter,8013)ncurr
---
> !      write(chariter,8013)ncurr
>          write(chariter,'("0",I3)')ncurr
>        case(1000:9999)
>          write(chariter,'(I4)')ncurr
