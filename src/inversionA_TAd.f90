Program inversionA_TAd
  !============================================================
  ! Potency density tensor inversion  (Y. Yagi)
  ! To be run after the optimal values of the hyperparameters are estimated by ABIC_new2. 
  ! Ifalg is set to 1, the LBFGSB can be used.
  !  Parameter
  !     imax:  The maximun number of Waveform
  !     jmax:  The maximun number of Station
  !     i_LBFGSB : Set i_LBFGSB to 1 to use LBFGSB. (default = 0)
  !     a_low_b  : Minimum value of displacement(icm=1) (only used when i_LBFGSB=1) 
  !-----------
  !     H(1:kmax,1:lmax): Kannel matrix 
  !     d(1:kmax): Data vector 
  !     a(1:lmax): Solution vecoter 
  !     G1(1:lmax,1:lmax) : time and space constrain matrix     
  !     Cm(1:kmax,1:kmax) : Covariance matrix of data
  !     d_* : double precision
  !  Estation Error problem
  !     The estimation error in the output seems to be not meaningful 
  !     because the non-diagonal components of the model covariance matrix are ignored.
  !     The calculation was changed so that estimation errors are not calculated.
  !============================================================
!  use, non_intrinsic:: optimize_lib, only: nnls_lbfgsb
  implicit none
  integer,parameter  :: imax=2501,jmax=201
  integer,parameter  :: i_LBFGSB=0
  real,parameter     :: a_low_b=-0.2
  !  Inversion Parameter
  real,allocatable   :: H(:,:), d(:), a(:), siga(:)
  real,allocatable   :: G1(:,:)
  real, allocatable  :: grn_m(:,:,:),obs_m(:,:,:)
  double precision, allocatable :: d_HtH(:,:), d_Htd(:), d_a(:), d_siga(:)
  double precision, allocatable :: d_H(:,:), d_d(:), d_Cm(:,:)
  double precision   :: d_btb, d_det_A
  !  Station Parameter
  character stcd(jmax)*10,comp(jmax)*4,title*50,cmode*1
  real,allocatable ::  Wobs(:,:)
  real cp(jmax),t1(jmax),dt(jmax),sigw(jmax),sigwd(jmax),cp1(jmax),tlength(jmax)
  real az(jmax),del(jmax)
  integer ndj(jmax)
  !  ABIC Parameter
  real beta(3),ABIC
  real,allocatable :: tr(:,:)
  integer,allocatable :: jtn(:,:),l_id_m(:,:,:,:)
  !  working space
  real,allocatable :: zz1(:),zz2(:)
  !-------
  integer :: i,j,k,l,n,l1,l2
  integer :: kmax,lmax,jn,ndmax,nmodel,ndata
  integer :: mn,nn,m0,n0,icmn
  integer :: jtn0,itera,nfabic,nsurface,nsurface_o,nflag
  integer :: nwk
  real    :: xx,yy,vr0,ylength,depth
  real    :: strike,dip,slip,rslip,rtime,st_max,r_s_t,cr
  real    :: shift,para1,para2,para3
  real    :: alpha1,alpha2,dump,s,var
  real    :: p_sec
  real    :: wk1,wk2,f_o,f_g
  real    :: det_Cm
  real    :: f_nrom2
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  !  Read Input Parameter
  read(5,'(a50)')title
  read(5,*)rtime,jtn0,vr0,shift,para1,para2,para3,st_max, r_s_t, itera, cr
  read(5,*)(beta(i),i=1,3),nfabic,nsurface_o,nflag,cmode,alpha1,alpha2,dump,s
  do j=1,jmax
     read(5,*,err=998,end =999)stcd(j),comp(j),cp(j),sigw(j),tlength(j),dt(j)
     cp(j)=cp(j)+shift
     jn = j
998  continue
  end do
999 continue
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  open(10,file="fault.dat")
  read(10,*)xx,yy,mn,nn,m0,n0,icmn,depth,rslip,nsurface,ylength,strike,dip,slip
  close(10)
  if(icmn.eq.1) rslip = 0.
  !---
  allocate(Wobs(imax,jn))
  !p_sec =  get_pre_sec(tlength,jn)
  open(60,file=".station.abic")
  read(60,*) jn,p_sec
  call readOBS_f(jn,imax,stcd,comp,Wobs,t1,dt,ndj,sigw,sigwd,tlength,p_sec)
  do j=1,jn
    read(60,*) stcd(j),comp(j),dt(j),sigw(j),sigwd(j),ndj(j)
  enddo
  close(60)
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  !---
  allocate(tr(mn,nn))
  call getTR(vr0,xx,yy,mn,nn,m0,n0,tr,rtime)   !Tr(m,n); Start time for each knot
  allocate(jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn))
  call getJTN(jtn,tr,mn,nn,jtn0,rtime,st_max)
  call get_l_id(mn,nn,jtn,jtn0,icmn,l_id_m)
  call get_stinfo(jn,stcd,comp,az,del)
  !----
  open(10,file="d_H.matrix",form='unformatted')
  read(10) ndata,nmodel
  kmax = ndata
  lmax = nmodel
  allocate(H(kmax,lmax),d(kmax))
  do k = 1, ndata
    read(10)d(k) !,(H(k,l),l=1,nmodel)
  enddo
  close(10)
  call r_green_H(mn,nn,jtn,jtn0,l_id_m,icmn,tr,rtime,jn,stcd,comp,dt,ndj,cp,H,kmax)
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  ndmax = maxval(ndj(1:jn))
  allocate(grn_m(ndmax,ndmax,jn),obs_m(ndmax,ndmax,jn))
  allocate(d_H(kmax,lmax), d_d(kmax), d_Cm(kmax,kmax))
  open(15,file="Cg.matrix",form='unformatted')
  do j=1,jn
     call get_covariance_obs_r1(obs_m(1:ndmax,1:ndmax,j),dt(j),ndmax,stcd(j),comp(j))
     read(15)nwk,nwk
     read(15)((grn_m(l1,l2,j),l1=1,ndj(j)),l2=1,ndj(j))
  enddo
  close(15)
  call correction_obs_m(ndata,d,obs_m,ndmax,jn,ndj,sigw,sigwd,f_o)
  call correction_grn_m(ndata,d,grn_m,ndmax,jn,ndj,f_g)
  call get_dCm(alpha2, obs_m, grn_m, ndmax, ndj, jn, d_Cm, kmax, det_cm)   !eq16-YF2011
  d_H = dble(H)
  d_d = dble(d)
  allocate(d_HtH(lmax,lmax), d_Htd(lmax))
  call get_HtH_Htd_dtd_d(kmax, lmax, d_H, d_d, d_Cm, d_HtH, d_Htd, d_btb)
  deallocate(d_H, d_d, d_Cm)
  deallocate(grn_m,obs_m)
  allocate(G1(lmax,lmax))
  call get_G12_matrix(mn,nn,jtn,icmn,jtn0,vr0,xx,yy,tr,l_id_m,rtime,st_max,r_s_t,itera,cr,nsurface,nsurface_o,G1,lmax)
!  open(13,file="G_12.matrix",form='unformatted')
!  do l1 = 1, nmodel
!    read(13)(G1(l1,l2),l2=l1,nmodel)
!  enddo
!  close(13)
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
!  call symetric_m(lmax,G1,lmax)
  !$omp parallel do
  do l=1,lmax
    d_HtH(1:lmax,l) = d_HtH(1:lmax,l) + dble(G1(1:lmax,l)*beta(1))
  enddo
  !$omp end parallel do
  deallocate(G1)
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  allocate(a(lmax),siga(lmax))
  allocate(d_a(lmax) )
  if( i_LBFGSB == 1 ) then
!    d_a(:) = nnls_lbfgsb(d_HtH, d_Htd, 5, 1.0d+7, 1.0d-5,dble(a_low_b))
    siga = 0.
  else
    call get_sol_posv_d(lmax, d_HtH, lmax, d_Htd, d_a, d_det_A)   !eq29-YF2011
    allocate(d_siga(lmax))
!----------------
!    Since estimation errors that do not take into account the covariance component 
!    of the model variance shoud be not meaningful, the calculation is skipped.
!    call get_sigx_LU_d(lmax,d_HtH,lmax,d_siga)                
!    siga = d_siga 
    siga = 0.    
!----------------
    deallocate(d_siga)
  endif
  a = d_a
  siga = sqrt(siga * (s/ndata))             ! It's not a good approximation.
  !-----------------------------------------
  !  allocate(d_x(lmax) )
  !  call get_sol_svd_d(lmax, d_HtH, lmax, d_Htd,d_x,1.0d-07)
  !  x = d_x
  !  siga = 0.
  !  deallocate(d_x)
  !-----------------------------------------
!  if( dump  /= 0. ) call get_dump_sol(a,nmodel,dump)    !eq32-YF2011
  open(16,file="x.vector",form='unformatted')
  write(16)(a(n),n=1,nmodel)
  close(16)
  wk1 = f_nrom2(ndata,d)
  allocate(zz1(max(lmax,kmax)),zz2(max(lmax,kmax)))
!  call multi_ab(ndata,nmodel,H,kmax,a,zz1)
  call SGEMV('N', ndata, nmodel, 1.0, H, kmax, a, 1, 0.0, zz1, 1)
  d(1:ndata) = d(1:ndata) - zz1(1:ndata)    ! d(k) = d(k) - H(k,l)*a(l)
  wk2 = f_nrom2(ndata,d)
  !---
  var = wk2/wk1
  write(6,'("check: ",3(f11.4,1x),f9.5,1x,f13.1)') para1,para2,para3,var,s
  !---
  !     Out Put Final Solution
  ABIC  = 0.
  if(cmode.eq."F")then
    cp1 = 0.
    call writesol(mn,nn,jtn,jtn0,l_id_m,icmn,a,tr,var,abic,strike,dip,slip,    &
          m0,n0,xx,yy,rtime,beta,siga,depth,stcd,comp,sigw,dt,   &
          cp,jn,vr0,nsurface,rslip,ylength,alpha1,alpha2)
  endif
  deallocate(wobs,zz1)
  !---
  stop
END Program inversionA_TAd
