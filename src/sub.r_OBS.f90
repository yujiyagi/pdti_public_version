! Y. Yagi
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
subroutine vector_d(imax,jn,Wobs,ndj,dt,b)
  implicit none
  integer,intent(in)  :: imax, jn
  integer,intent(in)  :: ndj(jn)
  real,intent(in)     :: Wobs(imax,jn),dt(jn)
  real,intent(out)    :: b(*)
  integer :: iwork,i,j,nwk
  iwork=0
  do j=1,jn
     ! Taper should be applied only when the low-pass filter is applied,
     ! so the below script was commented out (2019/04/09 yagi & smz)
     ! nwk = max(nint(ndj(j)*0.05),nint(1.5/dt(j)))
     ! call taper_all2(Wobs(1,j),ndj(j),nwk)
     ! If "dt" is set small enough to represent the characteristics of 
     !  the source time function, no filter is needed. (2023/03/26 Yagi)
     ! if (dt(j) > 1. ) call taper_all2(Wobs(1,j),ndj(j),nwk)
     nwk = max(nint(ndj(j)*0.025),nint(2.0/dt(j))) 
     call taper_tail(Wobs(1,j),ndj(j),nwk)  !To stably obtain product of covariance inverse matrix and obs. vector
     do i=1,ndj(j)
        iwork=iwork+1
        b(iwork)=Wobs(i,j)
     end do
  end do
  return
END subroutine vector_d
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
subroutine readOBS_f(jn,imax,stcd,comp,Wobs,t1,dt,ndj,sigw,sigwd,tlength,p_sec)
  implicit none
  integer,intent(in)      :: imax
  integer,intent(in)      :: jn
  character,intent(in)    :: stcd(jn)*10,comp(jn)*4
  real,intent(in)         :: dt(jn),tlength(jn)
  real,intent(in)         :: p_sec
  real,intent(out)        :: Wobs(imax,jn),t1(jn),sigw(jn),sigwd(jn)
  integer,intent(out)     :: ndj(jn)
  real,allocatable        :: wv(:),wvv(:)
  !---
  integer :: i,j,nstcd,nend,np0,nd0,nd1
  real :: sig,sig2,shift,dt0
  real :: wk1,wk2,wk3,wsum,wmax
  !---   
  Wobs = 0.
  !---
  do j=1,jn
     nstcd = index(stcd(j),' ') - 1
     open(10,file='wave.obs/'//stcd(j)(1:nstcd)//comp(j))
     read(10,*)t1(j),dt0,nd0 ,wk1,wk2,wk3,nend,sig
     nend = nend - nint(10./dt0)
     np0 = nint(t1(j)/dt0)-1
     allocate(wv(nd0),wvv(nd0))
     !-------
     read(10,*)(wv(i),i=1,nd0)
     !-------
     call diff_wave(wv,dt0,nd0)
     ! Waveform is shifted forward by 0.5*dt0.
     !
     call offset(wv,nd0,np0)
!     call bp_filter(wv,0.1,2.0,dt0,nd0)
     !--------
     !  The effect of the filter was excluded because putting the effect of the filter 
     !  in the covariance matrix corresponds to applying the filter to the waveform and 
     !  then deconvolving the applied filter. 
     !  Y. Yagi, Filtering Effect in Waveform Data Inversion for Seismic Source Process,
     !      Zisin, 66(4), p. 147-149, 2014, https://doi.org/10.4294/zisin.66.147
     !  Eq. 12 in YF2011 was skipped. (2024/03/28 Y.Yagi)
     !---------
     nd1 = nd0 
     if(nd0 > nend) nd1 = nend
     !-------
     shift = - p_sec
     shift = shift + dt0*0.5     !Corrects time deviations due to numerical differentiation (2024/06/13)
     ndj(j) = nint( (nd1*dt0 - p_sec + dt0*0.5 ) / dt(j) ) ! Correction (2025/01/07)
     call resample_shift(wv,nd0,dt0,wvv,nd0,dt(j),shift)
     !-------
     if( tlength(j) < ndj(j)*dt(j) ) ndj(j)=nint(tlength(j)/dt(j))
     !-------
     wsum=0.;do i=1,ndj(j); wsum=wsum+wvv(i)**2; enddo
     sig2=sqrt(wsum/ndj(j))
     !-------
     wmax=0.;do i=1,ndj(j); if(wmax<abs(wvv(i))) wmax=abs(wvv(i));enddo
     !-------
     Wobs(1:ndj(j),j) = wvv(1:ndj(j)) 
     !-------
     sigwd(j) = sig     ! Sigma of velocity wavefrom before P-wave arrival time
     sigw(j)  = wmax
     !---------------
     close(10)
     deallocate(wv,wvv)
  end do
  return
END subroutine readOBS_f
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
