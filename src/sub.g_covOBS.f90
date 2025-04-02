!
!  A program to calculate the covariance component of the observation error.  (Y. Yagi)
!  The effect of the filter was excluded because putting the effect of the filter 
!  in the covariance matrix corresponds to applying the filter to the waveform and 
!  then deconvolving the applied filter. 
!  Y. Yagi, Filtering Effect in Waveform Data Inversion for Seismic Source Process,
!      Zisin, 66(4), p. 147-149, 2014, https://doi.org/10.4294/zisin.66.147
!
subroutine get_covariance_obs_r1(wwm,dt,nd,stcd,comp)
  !  We assumed observation noise = back ground noise (no mechanical noise)
  !  implicit complex (z)
  implicit none
  character,intent(in) ::  stcd*10,comp*4
  integer,intent(in)   ::  nd
  real,intent(in)      ::  dt
  real,intent(out)     :: wwm(nd,nd)
  integer :: l
  !  real(4),allocatable:: gb(:),f(:,:),x(:),wv(:)
  !  complex:: zero(50),zpole(50)
  !------------------------------------------------------
  !ic1 = index(stcd," ")-1 ; ic2 = index(comp," ")-1
  !open(12,file="wave.obs/"//stcd(1:ic1)//comp(1:ic2)//".info")
  !read(12,*) dt0
  !read(12,*) nzero
  !do n=1,nzero
  !  read(12,*)a1,a2
  !   zero(n)=cmplx(a1,a2)
  !enddo
  !read(12,*) npole
  !do n=1,npole
  !  read(12,*)a1,a2
  !   zpole(n)=cmplx(a1,a2)
  !end do
  !close(12)
  !nwk = nint(nd*dt/dt0)*4+100
  !nn=log(real(nwk))/log(2.)
  !nn=2**nn
  !if(nn.lt.nwk) nn=nn*2
  !nd2 = nn
  !allocate(gb(nd2))
  !nwk = nint(dt/dt0)
  !do i= -9,9 
  !  n_input = nd2/2 + i
  !  if(mod(n_input,nwk) == 0) exit
  !enddo
  !!gb = 0.  ; gb(n_input:nd2) = 1.
  !! Velocity 
  !gb = 0.  ; gb(n_input) = 1.
  !! Velocity 
  !!call deconv(gb,dt0,zpole,zero,1.,nd2,npole,nzero-1)
  !!call taper_all(gb,nd2,0.2)
  !if(dt > 1.0 ) then 
  !   call lp_filter(gb,dt0,dt,nd2)
  !   call offset(gb,nd2,n_input-50)
  !   nd3 = nint(nd*dt/dt0) + 2
  !   !nd3 = nint(nd*dt/dt0) 
  !   allocate(f(nd3,nd),wv(nd2))
  !   nwk = nint(n_input * dt0/dt)
  !   f = 0.
  !   do l=1,nd3
  !     tw = (l-1)*dt0 - (n_input-2)*dt0
  !     !tw = (l-1)*dt0 - (n_input-1)*dt0
  !     call resample_shift(gb,nd2,dt0,wv,nd2,dt,tw)
  !     f(l,1:nd) = wv(1:nd)    !get f'
  !   enddo
  !   call multi_ata(nd3,nd,f,nd3,wwm,nd)
  !   sig_obs = maxval(wwm) * 1.e-7      ! To obtaine stable invers matrix
  !   forall(l=1:nd) wwm(l,l) = wwm(l,l) + sig_obs
  !else
!    neglecting filter effect when dt is enough small (2016/12/21)
  wwm = 0.0
  forall(l=1:nd) wwm(l,l) = 1.0
  !endif
  return
end subroutine get_covariance_obs_r1
