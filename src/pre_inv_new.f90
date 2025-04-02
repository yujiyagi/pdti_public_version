
Program Pre_inv_new
  !============================================================
  !  This program obtains the matrices H and the data vector d. (Y. Yagi)
  !============================================================
  !  Parameter
  !     H(1:kmax,1:lmax): Kannel matrix 
  !     d(1:kmax): Data vector 
  !  Waveform and Fault Parameter
  !     stcd(j): Component code for J wave
  !     comp(j): Component code for J wave
  !     wobs(i,j): Observation data at J station
  !     wcal(i,j): Caluculation data at J station
  !     g(i,m,n,ic): Green function  for ic vector at a subfault (m,n,)
  !     cp(j): delay time for J station (Tobs - Tcal)
  !     t1(j): Start time ; DT(j): Sampling Frequ. ; ndj(j): Number of data
  !     tlength(j) : Length of wavefrom used inversion (sec)
  !     sigw(j): Waight for J wave ; STCD(j): Staction code for J wave
  !  judgment parameter
  !     nsurface : if input 1, fault zone extend to free surface
  !============================================================
  implicit none
  integer,parameter :: imax=2501, jmax=201
  !  Inversion Parameter
  real,allocatable :: H(:,:),d(:)
  !  Station Parameter
  character :: stcd(jmax)*10,comp(jmax)*4,title*50
  real,allocatable ::  Wobs(:,:)
  real :: cp(jmax),t1(jmax),dt(jmax),sigw(jmax),sigwd(jmax),tlength(jmax)
  integer ndj(jmax)
  !---
  integer :: j,k
  integer :: lmax,kmax,ndata,nmodel,jn
  integer :: mn,nn,m0,n0,icmn,nsurface
  integer :: jtn0
  real :: xx,yy,depth,rslip,ylength
  real :: p_sec
  real :: rtime,vr0,shift,para1,para2,para3,st_max
  !---
  real,allocatable :: tr(:,:)
  integer,allocatable :: jtn(:,:),l_id_m(:,:,:,:)
  !  working space
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
  !------------------------------------------------------------------
  open(10,file="fault.dat")
  read(10,*)xx,yy,mn,nn,m0,n0,icmn,depth,rslip,nsurface,ylength
  allocate(tr(mn,nn))
  call getTR(vr0,xx,yy,mn,nn,m0,n0,tr,rtime)   !Tr(m,n); Start time for each knot
  close(10)
  if(icmn.eq.1) rslip = 0.
  !------------------------------------------------------------------
  open(60,file=".station.abic")
  read(60,*) jn,p_sec
  allocate(Wobs(imax,jn))
  call readOBS_f(jn,imax,stcd,comp,Wobs,t1,dt,ndj,sigw,sigwd,tlength,p_sec)
  ndata = 0
  do j=1,jn
     read(60,*) stcd(j),comp(j),dt(j),sigw(j),sigwd(j),ndj(j)
     ndata = ndata + ndj(j)
  enddo
  close(60)
  !------------------------------------------------------------------
  allocate(jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn))
  call getJTN(jtn,tr,mn,nn,jtn0,rtime,st_max)
  call get_l_id(mn,nn,jtn,jtn0,icmn,l_id_m)
  nmodel = maxval(l_id_m)
  write(6,'(" N_model:",i9," , min_jtn:",i4,"/",i4," , st_max:",f5.1)')  nmodel,minval(jtn),jtn0,st_max
  kmax = ndata
  lmax = nmodel
  if (nmodel .gt. 150000) then
     write(6,*) ; write(6,*)" ****************************************** "
     write(6,*)"   Warning: Too many model parameters. (> 150,000)"
     write(6,*) ; write(6,*)" ****************************************** "
  end if
  !------------------------------------------------------------------
  allocate(H(kmax,lmax),d(kmax))
  H = 0.; d=0.
  call vector_d(imax,jn,Wobs,ndj,dt,d)  ! b(k): observation vector
  !---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  open(10,file="d_H.matrix",form='unformatted')
  write(10) ndata,nmodel
  do k = 1, ndata
    write(10)d(k)!,(H(k,l),l=1,nmodel)
  enddo
  close(10)
  !---
  deallocate(H,d)
  stop
END Program pre_inv_new
