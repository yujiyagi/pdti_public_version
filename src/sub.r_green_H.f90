!Y. Yagi
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
subroutine  r_green_H(mn,nn,jtn,jtn0,l_id_m,icmn,tr,rtime, &
               jn,stcd,comp,dt,ndj,cp,  H,kmax)
  implicit none
  integer,intent(in)  :: mn,nn,jtn0,icmn,jn,kmax
  integer,intent(in)  :: jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn)
  real,intent(in)     :: rtime
  real,intent(in)     :: tr(mn,nn)
  character,intent(in):: stcd(jn)*10,comp(jn)*4
  real,intent(in)     :: dt(jn),cp(jn)  
  integer,intent(in)  :: ndj(jn)
  real,intent(out)    :: H(kmax,*)
  !  staion information
  !  working space
  real :: slip1(5)
  real,allocatable :: green(:),Umn(:),so(:),gb(:,:,:,:),so2(:)
  real    :: strike,dip,slip,rigid
  real    :: az(jn),del(jn)
  character :: cha*1
  integer :: m,n,j,icm,jt,i_mn
  integer :: l,k,i,ic1,ic2,iw,ndg
  real    :: shift,tg0,dtg,tmax
  real    :: work
  integer :: k_id
  integer :: nwk
  integer :: m_f(mn*nn*icmn),n_f(mn*nn*icmn),icm_f(mn*nn*icmn)
  !- Preparation for efficient parallelization
  do icm=1,icmn; do n=1,nn; do m=1,mn
    i_mn =(icm-1)*nn*mn + (n-1)*mn + m
    m_f(i_mn) = m
    n_f(i_mn) = n
    icm_f(i_mn) = icm
  enddo; enddo; enddo
  !-
  do j=1,jn
    tmax = ndj(j) * dt(j)
    ic1=index(stcd(j),' ') - 1 ; ic2 =index(comp(j),' ') - 1
    open(21,file='wave.grn/'//stcd(j)(1:ic1)//comp(j)(1:ic2)//cha(1), &
                                         form='unformatted')
    read(21,err=999);read(21,err=999);read(21,err=999)
    read(21)tg0,dtg,work,ndg
    close(21)
    iw  =  ndg *2
    allocate(green(iw),Umn(iw),so(iw),gb(iw,mn,nn,icmn),so2(iw))
    gb = 0.
    do icm=1,icmn
      open(20,file='wave.grn/'//stcd(j)(1:ic1)//comp(j)(1:ic2)//cha(icm), &
                                          form='unformatted')
      read(20,err=999)
      read(20,err=999)
      read(20,err=999)strike,dip,slip1(icm),rigid
      read(20)tg0,dtg,work,ndg,del(j),az(j)
      call stime_m(so,iw,dtg,rtime)
      do n=1,nn;do m=1,mn
        read(20,err=999)(gb(i,m,n,icm),i=1,ndg)
      enddo; enddo
      close(20)
    enddo
    nwk = max(nint(ndj(j)*0.025),nint(2.0/dt(j))) 
    !$omp parallel 
    !$omp do private(m,n,icm,green,jt,i,Umn)
    do i_mn = 1,mn*nn*icmn
      m = m_f(i_mn)
      n = n_f(i_mn)
      icm = icm_f(i_mn)
      call conv_y(gb(1:ndg,m,n,icm),So,green,iw)
      !======   Creating Kernel Matrix   ======
      do jt=1,jtn(m,n)
        call resample_shift(green,iw,dtg,Umn,iw,dt(j),Tr(m,n)+rtime*(jt-1)+cp(j)+tg0)
        call taper_tail(Umn,ndj(j),nwk)  !To stably obtain product of covariance inverse matrix and Kernel matrix (2024/05/30)
        do i=1,ndj(j)
          H(k_id(i,j,jn,ndj),l_id_m(m,n,jt,icm))=Umn(i)
        enddo
      enddo
    enddo
    !$omp end do
    !$omp end parallel
    deallocate(green,Umn,so,gb,so2)
  enddo       ! end j (station) loop
  if(icmn.eq.5)  slip = 90.
  if(icmn.eq.2)  slip=(slip1(1)+slip1(2))/2.
  if(icmn.eq.1)  slip=slip1(1)
  goto 998
999 continue
  write(6,'(a)')" ++++ Error in Read Green Function ++++ "
  write(6,'(" Please cheack :",a20)')'wave.grn/'//stcd(j)(1:ic1)//comp(j)//cha(icm)
  write(6,'(a)')" ++++++++++++++++++++++++++++++++++++++ "
  stop
998 continue
  return
END subroutine r_green_H
