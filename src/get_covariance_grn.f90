program get_covariance_grn
!========================================
! This progamme calculates the covariance matrix of the 
!  Green's function from the potency-rate density (slip-rate) function. (Y. Yagi)
!  (eq15-YF2011)
!========================================
  implicit none
  real,allocatable:: tr(:,:)
  real,allocatable:: dt(:)
  integer,allocatable:: ndj(:)
  character,allocatable:: stcd(:)*10,comp(:)*4
  !---------------------------------------
  character :: title
  integer   :: mn,nn,m0,n0,icmn,jtn0,nsurface
  integer   :: j, jn
  real      :: rtime,vr0
  real      :: xx,yy,depth,rslip,ylength
  real      :: p_sec
  real      :: wk1,wk2
  !---------------------------------------
  read(5,'(a50)')title
  read(5,*)rtime,jtn0,vr0
  !---------------------------------------
  open(10,file="fault.dat")
  read(10,*)xx,yy,mn,nn,m0,n0,icmn,depth,rslip,nsurface,ylength
  close(10)
  allocate(tr(mn,nn))
  call getTR(vr0,xx,yy,mn,nn,m0,n0,tr,rtime)
  !---------------------------------------
  open(60,file=".station.abic")
  read(60,*) jn,p_sec
  allocate(stcd(jn),comp(jn),dt(jn),ndj(jn))
  do j=1,jn
    read(60,*) stcd(j),comp(j),dt(j),wk1,wk2,ndj(j)
  enddo
  close(60)
  !---------------------------------------
  call covariance_grn(stcd,comp,ndj,dt,jn,tr,mn,nn,m0,n0,rtime,jtn0,icmn)
  deallocate(tr,dt,ndj)
  stop
end program get_covariance_grn
!========================================
!========================================
!========================================
subroutine covariance_grn(stcd,comp,ndj,dt,jn,tr,mn,nn,m0,n0,rtime,jtn0,icmn)
  !
  implicit none
  integer,intent(in)   :: jn,mn,nn,m0,n0,jtn0,icmn
  real,intent(in)      :: rtime
  real,intent(in)      :: tr(mn,nn)
  character,intent(in) :: stcd(jn)*10,comp(jn)*4
  real,intent(in)      :: dt(jn)
  integer,intent(in)   :: ndj(jn)
  !----------
  real,allocatable:: slip_rate(:,:,:,:),dly(:,:,:),   &
         wv(:),gb(:),eg(:,:),cgj(:,:,:),xmax(:,:,:,:),td_p(:,:,:,:)
  !---------
  integer :: imax,lkmax
  integer :: i,j,l,m,n,icm,l1,l2,lk
  integer :: nwk,ic1,ic2
  integer :: nds,ndt
  real    :: dtg,dts,dly0
  real    :: tw,xwk
  !section #### [ Get Slip_rate at earth space grid point ] ####
  !===========================================
  !######### Add Qeffect ####################
  open(20,file="Qeffect.dat")
  read(20,*) nds,dts
  imax = nint( ( maxval(Tr) + (jtn0+2)*rtime + 5.) * 2. / dts )
  allocate(wv(max(imax,nds)),gb(imax))
  do i=1,min(nds,imax)
    read(20,*) wv(i)
  enddo
  call diff_wave(wv,dts,imax)   ! We use velocity waveforms
  ! Waveform is shifted forward by 0.5*dts.
  close(20)
  !-----------
  allocate(slip_rate(imax,mn,nn,icmn))
  call get_slip_rate(slip_rate,imax,mn,nn,icmn,dts)
  !-----------
  do m=1,mn
    do n=1,nn; do icm = 1,icmn
      call conv_y(slip_rate(1,m,n,icm),wv,gb,imax)  !eq.9 in YF2011
      slip_rate(1:imax,m,n,icm) = gb(1:imax)
    enddo; enddo
  enddo
  deallocate(wv,gb)
  !===========================================
  open(61,file="Cg.matrix",form='unformatted')
  allocate(dly(mn,nn,jn),xmax(mn,nn,icmn,jn),td_p(mn,nn,icmn,jn))
  do j=1,jn
    ic1=index(stcd(j),' ') - 1 ; ic2 =index(comp(j),' ') - 1
    open(21,file='wave.grn/'//stcd(j)(1:ic1)//comp(j)(1:ic2)//'info', &
                                         form='unformatted')
    do n=1,nn;do m=1,mn
      read(21) dly(m,n,j)
    enddo;enddo
    do icm = 1,icmn
      do n=1,nn;do m=1,mn
        read(21) xmax(m,n,icm,j),td_p(m,n,icm,j)
      enddo;enddo
      if(icm == 1 ) then
        read(21) nwk
        read(21) (xwk,i=1,nwk)
        read(21) (xwk,i=1,nwk)
        read(21) (xwk,i=1,nwk)
        read(21) (xwk,i=1,nwk)
      endif
    enddo
    close(21)
    dly0 = dly(m0,n0,j)
    do n = 1,nn
      dly(1:mn,n,j) = dly(1:mn,n,j) - dly0
    enddo
  enddo
  !---
  allocate(cgj(maxval(ndj),maxval(ndj),jn))
  !$omp parallel do  private(dtg,ndt,lkmax,Eg,gb,lk,icm,m,n,l,tw)
  do j=1,jn
    dtg  =  dt(j)/2.
    ndt = nint(ndj(j)*(dt(j)/dtg)) + 10
    lkmax = icmn * mn * nn * ndt + 1
    allocate(Eg(lkmax,ndj(j)),gb(imax))
    !----------------------------- 
    !  The effect of the filter was excluded because putting the effect of the filter 
    !  in the covariance matrix corresponds to applying the filter to the waveform and 
    !  then deconvolving the applied filter. 
    !  Y. Yagi, Filtering Effect in Waveform Data Inversion for Seismic Source Process,
    !      Zisin, 66(4), p. 147-149, 2014, https://doi.org/10.4294/zisin.66.147
    !  Eq. 13 in YF2011 is skipped.
    !----------------------------- 
    lk = 0
    eg = 0.
    do icm =1,icmn     ! icm -> q  in YF2011
      do m=1,mn;  do n=1,nn   ! space(m,n)  -> k  in YF2011
        do l=1,ndt
          tw = dly(m,n,j) + (l-1) * dtg  !No error in the Green's function until P arrives. 
          tw = tw + 0.5*dts              !Corrects time deviations due to numerical differentiation (2024/06/13)
          tw = tw - 0.5*dt(j)            !The first time point contains the GF error and for stabilization (2024/06/13)
          if( tw > ndj(j)*dt(j) ) exit   ! exit l (time) loop
          call resample_shift(slip_rate(1,m,n,icm),imax,dts,gb,imax,dt(j),tw)
          lk = lk + 1
          eg(lk,1:ndj(j)) = gb(1:ndj(j)) * xmax(m,n,icm,j)   ! S_qkj*P_qkj(a) in YF2011
        enddo
      enddo; enddo
    enddo
    call ssyrk('u','t',ndj(j),lkmax,1.,Eg,lkmax,0.,Cgj(1,1,j),maxval(ndj))
    do l=1,ndj(j)
      Cgj(l,1:l-1,j) = Cgj(1:l-1,l,j) 
    enddo
    deallocate(Eg,gb)
  enddo       ! end j (station) loop
  !$omp end parallel do 
  do j=1,jn
    write(61)j,ndj(j)
    write(61)((cgj(l1,l2,j),l1=1,ndj(j)),l2=1,ndj(j))
  enddo       ! end j (station) loop
  close(61)
  return
end subroutine covariance_grn
!========================================
