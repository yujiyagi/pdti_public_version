! Y. Yagi
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
subroutine  get_ndj(mn,nn,jtn,icmn,tr,jn,stcd,comp,dt,cp,ndj,rtime)
!  parameter (p_sec = 8.0)
  implicit none
  integer,intent(in)   :: mn,nn,icmn,jn
  integer,intent(in)   :: jtn(mn,nn)
  real,intent(in)      :: tr(mn,nn)
  real,intent(in)      :: rtime
  character,intent(in) :: stcd(jn)*10,comp(jn)*4
  real,intent(in)      :: dt(jn),cp(jn)
  integer,intent(inout):: ndj(jn)
  !  working space
  real,allocatable     :: green(:),so(:),gb(:,:,:),Umnji(:,:)
  character            :: cha*1
  !----
  integer :: nf0(mn*nn),m_f(mn*nn),n_f(mn*nn)
  integer :: i,j,m,n,icm,iw,num,ndmax,ndg,nwk,i_mn
  integer :: ic1,ic2,nflag,ndj_0,ndj_get
  real :: dtmin,tmax,tw,shift
  real :: tg0,dtg
  real  :: work
  !--- Modified to determine lengths of used waveforms
  !--- based on the four sides of the model plane, not the four corners.
  !--- Ohara, Yagi & Yamashita 2021/03/11
  integer :: count
  num=mn*nn
  do m=1,mn;do n=1,nn
    i_mn = (n-1)*mn + m
    m_f(i_mn) = m
    n_f(i_mn) = n
  enddo;enddo
  ndmax = maxval(ndj)
  dtmin = minval(dt)
  do j=1,jn
    !================================
    tmax = ndj(j) * dt(j)
    ic1=index(stcd(j),' ') - 1 ; ic2 =index(comp(j),' ') - 1
    open(21,file='wave.grn/'//stcd(j)(1:ic1)//comp(j)(1:ic2)//cha(1), &
                                         form='unformatted')
    read(21,err=999);read(21,err=999);read(21,err=999)
    read(21)tg0,dtg,work,ndg
    iw  =  ndg *2
    allocate(green(iw),so(iw),gb(iw,mn,nn),Umnji(iw,num)) ! added by Ohara, Yagi & Yamashita 2021/03/11
    !================================
    call stime_m(so,iw,dtg,rtime)
    nflag = 0
    nf0   = 0
    ndj_0 = ndg*(dtg/dt(j))
    gb = 0.
    do n=1,nn;do m=1,mn
      read(21,err=999)(gb(i,m,n),i=1,ndg)
    enddo;enddo
    close(21)
    !$omp parallel do private(m,n,green,shift)
    do i_mn = 1, mn*nn
      m = m_f(i_mn) 
      n = n_f(i_mn) 
      if (jtn(m,n)>0) then
!      if( m==1 .or.  m==mn .or. n==1 .or. n==nn) then
        call conv_y(gb(1,m,n),So,green,iw)
        shift = Tr(m,n) + rtime*real(jtn(m,n)) + dt(j) + cp(j) + tg0  
        if(shift >= ndj(j)*dt(j))  nflag = 1
        call resample_shift(green,iw,dtg,Umnji(1,i_mn),iw,dt(j),shift)
        nf0(i_mn) = 1  
!      endif
      endif
    enddo
    nwk = ndj_get(umnji,iw,icmn,ndj_0,num,nf0) ! added by Ohara, Yagi & Yamashita 2021/03/11 - end
    if(nwk <= 3 ) nwk = ndj(j)
    if(nflag == 0 ) ndj(j) = min(ndj(j),nwk)
    deallocate(green,so,gb,Umnji)
  enddo
  goto 998
999 continue
  write(6,'(a)')" ++++ Error in Read Green Function ++++ "
  write(6,'(" Please cheack :",a20)')'wave.grn/'//stcd(j)(1:ic1)//comp(j)//cha(icm)
  write(6,'(a)')" ++++++++++++++++++++++++++++++++++++++ "
  stop
998 continue
  close(30)
  return
END subroutine get_ndj
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
function ndj_get(umnji,iw,icmn,nd,num,nf0)
  ! real:: Umnji(iw,4)
  implicit none
  integer :: ndj_get
  integer,intent(in)   :: iw,icmn,nd,num
  real,intent(in)      :: Umnji(iw,num) ! added by Ohara, Yagi & Yamashita 2021/03/11
  integer,intent(in)   :: nf0(num)
  integer:: ni         ! added by Ohara, Yagi & Yamashita 2021/03/11
  integer :: i
  real :: sum_pu,wk
  real :: pu(nd)
  real :: ndj_w(num)
  !---------------------------
  ndj_w = 0
  ndj_get = 0.
  !$omp parallel do private(i,pu,sum_pu,wk)
  do ni=1,num ! added by Ohara, Yagi & Yamashita 2021/03/11
    if(nf0(ni) == 0 ) cycle
    do i=1,nd
      pu(i) = Umnji(i,ni) ** 2
    enddo
    sum_pu = sum(pu)
    wk = 0.
    do i=1,nd
      wk = wk + pu(i)
      !  Note: If the tail of the Green's function is long (thick crust setting),
      !        you should consider using a value such as 0.7 instead of 0.90.
      !        0.95 -> 0.9  2024/02/26  Yuji Yagi
      if(wk >= sum_pu*0.95 ) then
        ndj_w(ni) = i
        exit
      endif
    enddo
  enddo
  ndj_get = maxval(ndj_w)
  return
end function ndj_get
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
