!  The sourutines were developed based on sourutines written by 
!  Prof. Kikuchi and Prof. Kanamori. Prof. Kanamori has kindly agreed to 
!  release to the public a modified version. Modified subroutines for ease 
!  of use and changed to fortran90 format (Yuji Yagi)
!
!  The core of sourutines were written by Prof. Kikuchi and Prof. Kanamori.
!  The original programs are available at "Note on Teleseismic Body-Wave 
!     Inversion Program (https://wwweic.eri.u-tokyo.ac.jp/ETAL/KIKUCHI/).
!==========================================================
subroutine lp_filter(x,dt,dtj,nm)
  !
  !     Anti-aliasing filter
  !     input x (raw data)
  !     output x (filtering data)
  !
  implicit none
  real,parameter :: pai= 3.1415926
  integer,intent(in) :: nm
  real,intent(in)    :: dt,dtj
  real,intent(inout) :: x(nm)
  complex,allocatable:: z(:)
  complex:: zc
  !----
  integer :: i
  integer :: np,nn
  real    :: amp
  real    :: df,dw,f,f3
  !----
  ! set parameter
  np = 2
  amp = 10.
  if(dtj.gt.2.0) then
     np = 10
     amp = 25.
  endif
  !
  f3  = 10.**(alog10(0.5/dtj)-1./np*alog10(amp-1.))
  nn=log(real(nm))/log(2.)
  nn=2**nn
  if(nn.lt.nm) nn=nn*2
  allocate(z(nn))
  z = 0.
  do i=1,nm
     z(i) = x(i)
  enddo
  call cfft(z,nn,-1)
  df = 1./(dt*nn)
  dw = df*pai*2
  do i = 2,nn/2
     f=df*(i-1)
     zc = 1/(1+(f/f3)**np)
     z(i) = z(i) *zc /nn
     z(nn-i+2)=conjg(z(i))
  enddo
  close(99)
  z(1) = 0.
  z(nn/2+1)=0.
  call cfft(z,nn,1)
  do i=1,nm
     x(i) = z(i)
  enddo
  deallocate(z)
  return
end subroutine lp_filter
!==========================================================
subroutine bp_filter(x,f1,f2,dt,nm)
  !
  !     input x (raw data)
  !     output x (filtering data)
  !
  implicit none
  real,parameter :: pai=3.1415926, amp=20.
  integer,parameter :: np=6
  integer,intent(in) :: nm
  real,intent(in)    :: f1,f2,dt
  real,intent(inout) :: x(nm)
  complex,allocatable ::z(:)
  complex:: zc1,zc2
  integer :: i,nn
  real    :: df,dw,f
  !
  nn=log(real(nm))/log(2.)
  nn=2**nn
  if(nn.lt.nm) nn=nn*2
  allocate(z(nn))
  z = 0.
  do i=1,nm
     z(i) = x(i)
  enddo
  call cfft(z,nn,-1)
  df = 1./(dt*nn)
  dw = df*pai*2
  do i = 2,nn/2
     f=df*(i-1)
     zc1 =    1./(1.+(f/f2)**np)
     zc2=  1.-1./(1.+(f/f1)**np)
     z(i) = z(i) *zc1*zc2 /nn
     z(nn-i+2)=conjg(z(i))
  enddo
  close(99)
  z(1) = 0.
  z(nn/2+1)=0.
  call cfft(z,nn,1)
  do i=1,nm
     x(i) = z(i)
  enddo
  deallocate(z)
  return
end subroutine bp_filter
!==========================================================
subroutine deconv(y,dt,zpole,zero,sc,nd,npole,nzero)
  !
  !     remove seismograph respons
  !     input y : seismic wavefrom
  !     output y : real velocity
  !
  implicit none
  integer,intent(in)    :: nd,npole,nzero
  real,intent(in)       :: sc,dt
  complex,intent(in)    :: zpole(npole),zero(nzero)
  real,intent(inout)    :: y(nd)
  real,parameter ::  pi=3.1415926
  complex,allocatable:: z(:)
  complex(16) :: zpole_d(npole),zero_d(nzero),zjw,zc1,zt1,zt2,zt,zc
  complex :: zr
  integer :: n,nn,i
  real    :: f,df,dw
  !
  do  n=1,npole
     zpole_d(n) = zpole_d(n)
  enddo
  do n=1,nzero
     zero_d(n) = zero(n)
  enddo
  nn=log(real(nd))/log(2.)
  nn=2**nn
  if(nn.lt.nd) nn=nn*2
  allocate(z(nn))
  do  i=1,nn
     z(i)=0.
  end do
  do  i=1,nd
     z(i)=y(i)
  end do
  call cfft(z,nn,-1)
  df=1/(dt*nn)
  dw=df*pi*2.
  do  i=2,nn/2
     f=df*(i-1)
     zjw=cmplx(0.,dw*(i-1))
     !=============== filter ==================
     !     zc1=1.-1./(1.+(f/0.001)**12)
     zc1=1.
     !-------------------
     zt1=1.
     do  n=1,npole
        zt1 = zt1 * (zjw-zpole(n))
     end do
     zt2=1.
     do n=1,nzero
        zt2 = zt2 * (zjw-zero(n))
     end do
     zt = zt1/zt2
     zc=zc1/nn*zt*sc
     !   diff waveform
     !   zc=zc*cmplx(0,2*pi*f)
     !-----------------------------------------
     zr = zc
     z(i)=z(i)*zr
     z(nn-i+2)=conjg(z(i))
  end do
  z(1)=0.
  z(nn/2+1)=0.
  call cfft(z,nn,1)
  do  i=1,nd
     y(i)=z(i)
  end do
  deallocate(z)
  return
end subroutine deconv
!==========================================================
subroutine cfft(x,n,id)
  ! <  fft for complex variables >
  implicit none
  integer,intent(in)    :: n,id
  complex,intent(inout) :: x(n)
  real :: pi
  integer :: i,j,k,l,m,ni,n2
  complex :: z,zk,xi,xj
  pi=sign(3.141593,id*1.)
  n2=n/2
  l=n
1 ni=l/2
  z=cexp(cmplx(0.,pi/ni))
  zk=1
  do  k=1,ni
     do  m=1,n/l
        j=k+2*(m-1)*ni
        xj=x(j)
        x(j)=xj+x(j+ni)
        x(j+ni)=(xj-x(j+ni))*zk
     enddo
     zk=zk*z
  enddo
  l=ni
  if(l.gt.1) goto 1
  do i=1,n
     j=i-1
     m=n2
     l=1
5    l=l+m*(j-2*(j/2))
     m=m/2
     j=j/2
     if(j.ge.1) goto 5
     if(i.gt.l) goto 10
     xi=x(i)
     x(i)=x(l)
     x(l)=xi
10   continue
  enddo
end subroutine cfft
