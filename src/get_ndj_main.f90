Program get_ndj_main
  !============================================================
  !  The length of the waveform at each observation used for inversion 
  !   is determined from the Green's function at the edge of the model plane. (Y. Yagi)
  !============================================================
  implicit none
  integer,parameter :: imax=2501,jmax=201
  !  Inversion Parameter
  !  Station Parameter
  character :: stcd(jmax)*10,comp(jmax)*4,title*50
  real,allocatable ::  Wobs(:,:)
  real ::  cp(jmax),t1(jmax),dt(jmax),sigw(jmax),sigwd(jmax),tlength(jmax)
  integer ndj(jmax)
  integer,allocatable:: jtn(:,:)
  real,allocatable :: tr(:,:)
  !--------
  integer :: j,jn
  integer :: mn,nn,m0,n0,icmn,jtn0
  real    :: rtime,vr0,shift,para1,para2,para3,st_max
  real    :: xx,yy
  real    :: p_sec
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  !  Read Input Parameter
  read(5,'(a50)')title
  read(5,*)rtime,jtn0,vr0,shift,para1,para2,para3,st_max
  read(5,*)
  do j=1,jmax
     read(5,*,err=998,end =999)stcd(j),comp(j),cp(j),sigw(j),tlength(j),dt(j)
     cp(j)=cp(j)+shift
     jn = j
998  continue
  end do
999 continue
  !---
  open(10,file="fault.dat")
  read(10,*)xx,yy,mn,nn,m0,n0,icmn
  allocate(tr(mn,nn))
  call getTR(vr0,xx,yy,mn,nn,m0,n0,tr,rtime)   !Tr(m,n); Start time for each knot
  allocate(jtn(mn,nn))
  call getJTN(jtn,tr,mn,nn,jtn0,rtime,st_max)

  close(10)
  allocate(Wobs(imax,jn))
  p_sec =  10.
  call readOBS_f(jn,imax,stcd,comp,Wobs,t1,dt,ndj,sigw,sigwd,tlength,p_sec)
  call get_ndj(mn,nn,jtn,icmn,tr,jn,stcd,comp,dt,cp,ndj,rtime)
  !------------------------------------------------------------------
  open(60,file=".station.abic")
  write(60,'(i10,f10.2)') jn,p_sec
  do j=1,jn
    write(60,'(a10,1x,a4,1x,f7.3,1x,e15.4,1x,e15.4,1x,i10,1x,f10.2)')  &
                  stcd(j),comp(j),dt(j),sigw(j),sigwd(j),ndj(j),cp(j)
  enddo
  close(60)
  !------------------------------------------------------------------
  stop
END Program get_ndj_main
