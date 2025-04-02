subroutine get_slip_rate(R3,imax,mn,nn,icmn,dt)
  ! The program computes a potency-rate density (slip-rate) function 
  ! to estimate the covariance matrix of the Green's function. (Y. Yagi)
  !
  implicit none
  integer,intent(in) :: imax,mn,nn,icmn
  real,intent(in)    :: dt
  real,intent(out)   :: R3(imax,mn,nn,icmn)
  real,allocatable :: PD_dist(:,:),PDT_comp(:,:,:),PDT_rate(:,:,:,:),Tr(:,:), &
                      slip(:,:)
  real,allocatable :: Ro(:),Ro2(:)
  integer  :: jt,m,n,icm
  integer  :: m0,n0,jtn
  integer  :: nwk
  real     :: xx,yy,rtime,ylength
  real     :: shift
  !----------------------------------------------
  read(40,*); read(40,*)
  read(40,*); read(40,*)
  read(40,*); read(40,*)xx,yy,nwk,nwk,m0,n0,rtime,jtn,nwk,ylength
  read(40,*); read(40,*)
  read(40,*) 
  allocate (PD_dist(mn,nn),PDT_comp(mn,nn,icmn),PDT_rate(mn,nn,jtn,icmn), &
            Tr(mn,nn), slip(mn,nn))
  do n=nn,1,-1 ; read(40,*)(PD_dist(m,n),m=1,mn) ; end do
  read(40,*) 
  do n=nn,1,-1 ; read(40,*)(slip(m,n),m=1,mn) ; end do
  read(40,*) 
  do n=nn,1,-1 ; read(40,*)(Tr(m,n),m=1,mn)   ; end do
  read(40,*) 
  do icm=1,icmn ; read(40,*) 
    do n=nn,1,-1 ; read(40,*)(PDT_comp(m,n,icm),m=1,mn) ; end do
  end do
  do jt=1,jtn ; do icm=1,icmn ;read(40,*) 
    do n=nn,1,-1 ; read(40,*)(PDT_rate(m,n,jt,icm),m=1,mn) ; end do
  end do ;  end do
  !---------------------------------------------
  ! Slip-rate Time Function
  !---------------------------------------------
  allocate (Ro(imax),Ro2(imax))
  Ro=0.; R3 = 0.
  call stime_m(Ro,imax,dt,rtime)
  do icm=1,icmn ; do m=1,mn ; do n=1,nn 
    do jt=1,jtn
      shift = rtime * (jt - 1) + tr(m, n)
      call resample_shift(Ro, imax, dt, Ro2, imax, dt, shift)
      R3(1:imax,m,n,icm) = R3(1:imax,m,n,icm) + PDT_rate(m,n,jt,icm)*Ro2(1:imax) 
    end do
  end do; end do; end do
  return 
end subroutine get_slip_rate
