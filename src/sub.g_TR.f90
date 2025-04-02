! Y. Yagi
subroutine getTR(vr0,xx,yy,mn,nn,m0,n0,tr,rtime)
  implicit none
  integer,intent(in) :: mn,nn,m0,n0
  real,intent(in)    :: vr0,xx,yy,rtime
  real,intent(out)   :: tr(mn,nn)
  real :: tr_r(mn,nn)
  integer :: m,n,ios
  real :: wkm,wkn,dist,nw
  do m=1,mn
     do n=1,nn
        if(m > m0 ) then
          wkm = - xx
        else
          wkm = + xx
        endif
        if(n > n0 ) then
          wkn = - yy
        else
          wkn = + yy
        endif
        if(m == m0 ) wkm = 0.
        if(n == n0 ) wkn = 0.
        dist = sqrt((xx*(m-m0)+wkm)**2+(yy*(n-n0)+wkn)**2) !-wk
        tr(m,n)= dist/vr0 + 0.1
        tr(m0,n0) = 0.
     end do
  end do
  !=============================================
  do m=1,mn
     do n=1,nn
       nw = nint(tr(m,n)/rtime -0.49999)
       tr(m,n) = nw * rtime
     end do
  end do
  !=============================================
  open(11, file = 'TR_fault_segments.txt', status = 'old', iostat = ios)
  if (ios == 0) then
    do n = nn, 1, -1
      read(11,*, iostat=ios)  (tr_r(m, n), m = 1, mn)
      if (ios < 0) goto 99
    end do
    tr = tr_r
  endif
99 close(11)
  return
END subroutine getTR
!----
subroutine getTR0(vr0,xx,yy,mn,nn,m0,n0,tr)
  implicit none
  integer,intent(in) :: mn,nn,m0,n0
  real,intent(in)    :: vr0,xx,yy
  real,intent(out)   :: tr(mn,nn)
  integer :: m,n
  real :: dist
  do m=1,mn
     do n=1,nn
        dist = sqrt((xx*(m-m0))**2+(yy*(n-n0))**2)
        tr(m,n)= dist/vr0
     end do
  end do
  return
END subroutine getTR0
!----
