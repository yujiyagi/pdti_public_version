!
! Waveform processing subroutines (Y. Yagi)
!
!==========================================================
subroutine taper_all2(x,nterm,m1)
  !-------------------------------------
  !       Taper Wavefrome Program
  !       cut Cosain taper
  !       x ; waveform
  !       nterm ; number of waveform
  !       widtap ; reate of taper in waveform
  !------------------------------------
  implicit none
  integer,intent(in)  :: nterm, m1
  real,intent(inout)  :: x(nterm)
  ! 
  real :: pi,ca,cb
  integer :: m2,i
  pi=3.141593
  m2=m1
  do i=1,nterm
     if(i.le.1+m1) then
       ca=0.5*(1.0+cos(pi*real(m1+1-i)/real(m1)))
       x(i)=x(i)*ca
     endif
     if(i.ge.nterm-m2) then
       cb=0.5*(1.0+cos(pi*real(i-nterm+m2)/real(m2)))
       x(i)=x(i)*cb
     endif
  enddo
  return
end subroutine taper_all2
!==========================================================
subroutine taper_tail(x,nterm,m1)
  !-------------------------------------
  !       Taper Wavefrome Program
  !       cut Cosain taper
  !       x ; waveform
  !       nterm ; number of waveform
  !       widtap ; reate of taper in waveform
  !------------------------------------
  implicit none
  integer,intent(in)  :: nterm, m1
  real,intent(inout)  :: x(nterm)
  ! 
  real :: pi,ca,cb
  integer :: m2,i
  pi=3.141593
  m2=m1
  do i=1,nterm
     if(i.ge.nterm-m2) then
       cb=0.5*(1.0+cos(pi*real(i-nterm+m2)/real(m2)))
       x(i)=x(i)*cb
     endif
  enddo
  return
end subroutine taper_tail
!==========================================================
subroutine taper_all(x,nterm,widtap)
  !-------------------------------------
  !       Taper Wavefrome Program
  !       cut Cosain taper
  !       x ; waveform
  !       nterm ; number of waveform
  !       widtap ; reate of taper in waveform
  !------------------------------------
  implicit none
  integer,intent(in)  :: nterm
  real,intent(in)     :: widtap
  real,intent(inout)  :: x(nterm)
  !
  real :: pi,ca,cb
  integer  ::  m1,m2,i
  pi=3.141593
  m1=int(widtap*real(nterm))
  m2=int(widtap*real(nterm))
  do i=1,nterm
     if(i.le.1+m1) then
       ca=0.5*(1.0+cos(pi*real(m1+1-i)/real(m1)))
       x(i)=x(i)*ca
     endif
     if(i.ge.nterm-m2) then
       cb=0.5*(1.0+cos(pi*real(i-nterm+m2)/real(m2)))
       x(i)=x(i)*cb
     endif
  enddo
  return
end subroutine taper_all
!==========================================================
subroutine diff_wave(x,dt0,in0)
!
!   Differential Calculus
!
  implicit none
  integer,intent(in) :: in0
  real,intent(in)    :: dt0
  real,intent(inout) :: x(in0)
  !---
  integer :: i
  real,allocatable:: dx(:)
  allocate(dx(in0))
  dx = 0.
!----
  do i=1,in0-1
    dx(i) = (x(i+1)-x(i))/dt0
  enddo
!----
  x = dx
  deallocate(dx)
end subroutine diff_wave
!==========================================================
subroutine int_wave(dt,nn,ddy,ivd)
!
!  Integral Calculus
!  (ivd=1 : single integrate, ivd=2 : twice integrate )
!
  implicit none
  integer,intent(in) :: nn,ivd
  real,intent(in)    :: dt
  real,intent(inout) :: ddy(nn)
  !--
  real,allocatable:: dy(:),y(:)
  integer :: m
  allocate(dy(nn),y(nn))
  dy(1)=0.0
  y(1)=0.0
  do m=2,nn
     dy(m)=dy(m-1)+(ddy(m-1)+ddy(m))*dt/2.
     y(m)=y(m-1)+dy(m-1)*dt+(ddy(m-1)/3.+ddy(m)/6.)*dt**2
  end do
  if(ivd.eq.1) ddy=dy
  if(ivd.eq.2) ddy=y
  deallocate(dy,y)
  return
END subroutine int_wave
!==========================================================
subroutine rmean(x,nterm)
  !-------------------------------------
  !       remove mean Wavefrome Program
  !       x ; waveform
  !       nterm ; number of waveform
  !------------------------------------
  implicit none
  integer,intent(in) :: nterm
  real,intent(inout) :: x(nterm)
  !
  real :: pi, dm
  integer :: i
  pi=3.141593
  dm= 0.
  do i=1,nterm
     dm = dm + x(i)
  end do
  dm = dm / nterm
  x = x - dm
  return
end subroutine rmean
!==========================================================
subroutine offset(x,nd,m)
  !-----------------
  !  x  : waveform ( input/output)
  !  nd : maximum number of point (input)
  !  m  : we calculate offset using x(0:m)  (input)
  !-----------------
  implicit none
  integer,intent(in) :: nd,m
  real,intent(inout) :: x(nd)
  !--
  real :: xm
  integer :: i
  !--
  xm = 0.
  do i=1,m
     xm = xm + x(i)
  end do
  xm = xm/real(m)
  x = x -xm
  return
end subroutine offset
!==========================================================
subroutine resample_shift(w1,imax1,dt1,w2,imax2,dt2,shift)
 implicit none
 integer, intent(in) :: imax1, imax2
 real,    intent(in) :: dt1,   dt2
 real,    intent(in) :: shift
 real,    intent(in) :: w1(imax1)
 real,    intent(out):: w2(imax2)
 !-----
 real :: trend, t1, t2
 integer :: i, i1
 w2 = 0.
 do i=1,imax2
   t2 = dt2*(i-1) - shift
   i1 = int(t2/dt1) + 1
   t1 = dt1*(i1-1)
   if(i1 <= 0  .or.  i1 >= imax1) cycle
   if(t2 < 0.) then
     trend = (w1(1) - 0. )/dt1
     w2(i) =  w1(1) + trend*(t2)
   else
     trend = (w1(i1+1)-w1(i1))/dt1
     w2(i) = w1(i1) + trend*(t2-t1)
   endif
 end do
 return
end subroutine resample_shift
!==========================================================
subroutine readsac(fname,x,nm,nd,beg,pick,dt,nerr)
  !  Read SAC file (ASCII)
  !----------------
  !  fname : sac file name  (input)
  !  nm : maximum number of point (input)
  !----------------
  !  x :  waveform (output)
  !  nd : maximum number of point (output)
  !  beg : start point (output)
  !  pick : p-pick point (output)
  !  dt : sampling interval (output)
  !  nerr  :  mode (output)
  !----------------
  implicit none
  character,intent(in) :: fname*40
  integer,intent(in)   :: nm
  real,intent(out)     :: x(nm),beg,pick,dt
  integer,intent(out)  :: nd,nerr
  character            :: stcd*8,comp*8
  real                 :: wk
  integer              :: i
  integer              :: iyear,imd,ihour,imin,isec,idsec,iw
  x = 0.
  open(11,file=fname)
  read(11,*)dt
  read(11,*)beg,wk,wk,pick
  do i=1,12
     read(11,*)
  end do
  read(11,*)iyear,imd,ihour,imin,isec
  read(11,*)idsec,iw,iw,iw,nd
  do i=1,6
     read(11,*)
  end do
  read(11,*) stcd
  do i=1,5
     read(11,*)
  end do
  read(11,*)iw,iw,comp
  read(11,*)
  read(11,*)(x(i),i=1,nd)
  close(11)
  nerr = 0
  return
end subroutine
!==========================================================
subroutine stime_m(x,n,dt,t1)
  !   Source time function
  implicit none
  real,intent(in)     :: dt ,t1
  integer,intent(in)  :: n
  real,intent(out)   :: x(n)
  !---
  real :: t2
  integer :: i, it1,it2
  !---
  t2 = t1 * 2.
  do i=1,n
     x(i) = 0.
  end do
  it1=t1/dt+1.5
  it2=T2/dt+1.5
  do i=2,it1
     x(i)=2.*(i-1)*dt/(t1*t2)
  end do
  do i=it1+1,it2
     x(i)=2.*(t2-(i-1)*dt)/(t2-t1)/t2
  end do
  return
END subroutine stime_m
!==========================================================
