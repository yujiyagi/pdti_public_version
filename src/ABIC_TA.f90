Program ABIC_TA
  !============================================================
  !  Calculation ABIC with proper Covariance Matrix  ( Y. Yagi)
  !     H(1:kmax,1:lmax): Kannel matrix 
  !     d(1:kmax): Data vector 
  !     a(1:lmax): Solution vecoter 
  !     G1(1:lmax,1:lmax) : time and space constrain matrix     
  !     Cm(1:kmax,1:kmax) : Covariance matrix of data
  !     d_* : double precision
  !============================================================
  implicit none
  real,allocatable :: H(:,:),d(:),G1(:,:)
  real,allocatable :: obs_m(:,:,:),grn_m(:,:,:) 
  !  ABIC Parameter
  real beta(3)
  real,allocatable :: abic(:),abic_min(:),beta_min(:)
  !  Station Parameter
  character,allocatable:: stcd(:)*10,comp(:)*4
  real,allocatable:: dt(:),sigw(:),sigwd(:),cp(:)
  integer,allocatable:: ndj(:)
  character:: cmode*1
  real, parameter:: pi=3.14159265
  real,allocatable :: tr(:,:)
  integer,allocatable :: jtn(:,:),l_id_m(:,:,:,:)
  !--------------------
  real, allocatable :: Cm(:,:), HtH(:,:) , Htd(:), HtH_r(:,:)
  real, allocatable :: a(:),siga(:)
  real(8), allocatable :: d_HtH_r(:,:),d_Htd(:),d_a(:)
  real(8)              :: d_det_HtH_G
  real :: dtd
  integer :: info
  !--------------------
  integer :: i,j,k,l,l1,l2,ia2,ia2_0,ita2,ib1,itb1
  integer :: ia2_1,ia2_2,ib1_2
  integer :: ndata,nmodel,kmax,lmax,jn,ndmax,nfabic,itera,jtn0
  integer :: mn,nn,m0,n0,icmn,nsurface,nsurface_o
  integer :: nwk,nflag_a,nflag_w,nflag_b
  real    :: shift,para1,para2,para3,alpha1,alpha2,dump
  real    :: st_max,r_s_t,cr
  real    :: xx,yy,depth,rslip,ylength,rtime,vr0
  real    :: a2_0,d_a2,b1_0,d_b1,b1_m
  real    :: abice,abic_wk,beta1,abic_1,abic_2,beta_wk
  real    :: det_HtH_G, det_Cm, det_G, s 
  real    :: sig_s_g,sig_s_o,f_o,f_g
  real    :: get_hyper_para
  double precision :: d_beta1
  real    :: strike,dip,rake
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  !  Read Input Parameter
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  read(5,*)
  read(5,*)rtime,jtn0,vr0,shift,para1,para2,para3,st_max, r_s_t, itera, cr
  para2 = 0.0
  read(5,*)(beta(i),i=1,3),nfabic,nsurface_o,nwk,cmode,alpha1,alpha2,dump
  !---
  open(10,file="fault.dat")
  read(10,*)xx,yy,mn,nn,m0,n0,icmn,depth,rslip,nsurface,ylength,strike,dip,rake
  close(10)
  open(60,file=".station.abic")
  read(60,*) jn
  allocate(stcd(jn),comp(jn),dt(jn),sigw(jn),sigwd(jn),cp(jn),ndj(jn))
  do j=1,jn
    read(60,*) stcd(j),comp(j),dt(j),sigw(j),sigwd(j),ndj(j),cp(jn)
  enddo
  close(60)
  ndmax = maxval(ndj)
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  !  Read Matrix 
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  open(10,file="d_H.matrix",form='unformatted')
  read(10) ndata,nmodel
  kmax = ndata
  lmax = nmodel
  !---
  allocate(H(kmax,lmax),d(kmax)) 
  do k = 1, ndata
     read(10)d(k) !,(H(k,l),l=1,nmodel)
  enddo
  close(10)
  !--------
  allocate(tr(mn,nn),jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn))
  call getTR(vr0,xx,yy,mn,nn,m0,n0,tr,rtime)  
  call getJTN(jtn,tr,mn,nn,jtn0,rtime,st_max)
  call get_l_id(mn,nn,jtn,jtn0,icmn,l_id_m)
  call r_green_H(mn,nn,jtn,jtn0,l_id_m,icmn,tr,rtime,jn,stcd,comp,dt,ndj,cp,H,kmax)
  !--------
  allocate(G1(lmax,lmax))
  call get_G12_matrix(mn,nn,jtn,icmn,jtn0,vr0,xx,yy,tr,l_id_m,rtime,st_max,r_s_t,itera,cr,nsurface,nsurface_o,G1,lmax)
  allocate(grn_m(ndmax,ndmax,jn),obs_m(ndmax,ndmax,jn))
  open(15,file="Cg.matrix",form='unformatted')
  do j=1,jn
     call get_covariance_obs_r1(obs_m(1,1,j),dt(j),ndmax,stcd(j),comp(j))
     read(15)nwk,nwk
     read(15)((grn_m(l1,l2,j),l1=1,ndj(j)),l2=1,ndj(j))
  enddo
  close(15)
  call correction_obs_m(ndata,d,obs_m,ndmax,jn,ndj,sigw,sigwd,f_o)
  call correction_grn_m(ndata,d,grn_m,ndmax,jn,ndj,f_g)
  !-------
  allocate(HtH(lmax,lmax),Htd(lmax),a(lmax))
  !-------
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  !  Grid Search covariance matrix where  Min. ABIC
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  call set_alpha_para0(ita2,a2_0,d_a2)
  nflag_a = 0
  allocate(beta_min(100))
801 continue 
  !=+=+=+=+=+  TMP
  nflag_w = 0
  !=+=+=+=+=+  TMP
  allocate(abic_min(ita2))
  abic_min = 1.e15
  abice = 99999999.
  ia2_0 = 1
  !=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  do ia2 = ia2_0,ita2          ! Alpha2 Loop  (Uncertainty of Green's function)
     alpha2 = get_hyper_para(a2_0,d_a2,ia2)
     if(nflag_a >= 1 .and. ia2==2) then
        abic_min(ia2) = abic_wk
        beta_min(ia2) = beta_wk
        goto 802
     endif
     allocate(Cm(kmax,kmax))
     call get_Cm(alpha2, obs_m, grn_m, ndmax, ndj, jn, Cm, kmax, det_Cm)  !eq16-YF2011
     call get_HtH_Htd_dtd(kmax, lmax, H, d, Cm, HtH, Htd, dtd)
     deallocate(Cm)
     !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
     !  Grid Search hyper parameter of Constraint where Min. ABIC
     !-------------------
     call set_beta_para0(itb1,b1_0,d_b1)
     if(nflag_w == 1) b1_0 =  b1_m
     nflag_b = 0
701  continue
     !-------------------
     allocate(abic(itb1));  abic = 1.e15
     !-----
     allocate(HtH_r(lmax,lmax))
     do ib1 = 1,itb1                   !   beta1 Loop (Smoothing)
        beta1 = get_hyper_para(b1_0,d_b1,ib1)
        !-------------------
        if( nflag_b > 0 .and. ib1 == 2 ) then
           abic(ib1) = abic_min(ia2)
           if(abic(ib1) - abic(ib1-1) > 2. ) exit
           cycle
        endif
        det_G = lmax * log(beta1)   ! Single hyper-para of smoothing.
        d_beta1 = dble(beta1)
!----
        !$omp parallel do
        do l = 1, lmax
          call scopy(lmax, HtH(1,l) , 1, HtH_r(1,l) , 1)
          call saxpy(lmax, beta1, G1(1,l),1, HtH_r(1,l),1)
        enddo
        !$omp end parallel do
        call get_sol_posv   (lmax, HtH_r, lmax, Htd, a, det_HtH_G,info) !eq29-YF2011
        !$omp parallel do
        do l = 1, lmax
          call scopy(lmax, HtH(1,l) , 1, HtH_r(1,l) , 1)
          call saxpy(lmax, beta1, G1(1,l),1, HtH_r(1,l),1)
        enddo
        !$omp end parallel do
        call get_residual (lmax, HtH_r, Htd, dtd, a, s)          !eq28-YF2011
!----
        abic(ib1) = ndata*alog(s) - det_G + det_HtH_G + det_Cm          !eq31-YF2011
        if(abs(abic(ib1)) > abice)  abic(ib1) = abice
        sig_s_g = sqrt(f_g * (s/ndata))
        sig_s_o = sqrt(f_o*(s/ndata)*alpha2)
        if(info .ne. 0 ) then
          abic(ib1) = 9999999000.
          s = 9999999000.
          sig_s_g = 9999999000.
          sig_s_o = 9999999000.
        endif
        write(6,'("ABIC :",e12.5,3f4.1,1x,e12.5,3(1x,f7.1),1x,f20.2,1x,4f18.2,2e10.3)') &
             beta1,0.0,0.0,0.0,alpha2,para1,para2,para3,abic(ib1), &
             det_G,det_HtH_G,det_Cm,s,sig_s_g,sig_s_o
        if(ib1 > 1 ) then
           if(abic(ib1) - abic(ib1-1) > 2. ) exit
        endif
     enddo  ! ib1 Loop
     deallocate(HtH_r)
     !---
     call get_min_value_1(abic,itb1,abic_2,ib1_2)
     abic_min(ia2) = min(abic_min(ia2), abic_2)
     beta_min(ia2) = get_hyper_para(b1_0,d_b1,ib1_2)
     deallocate(abic)
     if(nflag_b == 0 ) then
        if(ib1_2==1) then
           b1_0 = b1_0 - 2.
           goto 701
        else
           b1_m = b1_0 + real(ib1_2) - 2.
           nflag_w = 1
        endif
     endif
     call set_beta_para1(itb1,b1_0,d_b1,ib1_2)
     nflag_b = nflag_b + 1
     if (nflag_b <= 3 ) goto 701
     !---
802  continue
     call get_min_value_1(abic_min,ita2,abic_1,ia2_1)
     write(6,'("alpha2, abic :",1e10.3,2f14.2)')  alpha2,abic_min(ia2),abic_1
     if(abic_min(ia2) - abic_1 > 2. ) exit
  enddo      ! Alpha2 Loop
  call get_min_value_1(abic_min,ita2,abic_2,ia2_2)
  !----
  deallocate(abic_min)
  abic_wk = abic_2
  beta_wk = beta_min(ia2_2)
  nflag_a = nflag_a + 1
  if (nflag_a <= 1 ) then 
    call set_alpha_para1(ita2,a2_0,d_a2,ia2_2)
    goto 801
  endif
!----------------------------
  beta1 = beta_min(ia2_2)
  alpha2 =  get_hyper_para(a2_0,d_a2,ia2_2)
  allocate(Cm(kmax,kmax))
  call get_Cm(alpha2, obs_m, grn_m, ndmax, ndj, jn, Cm, kmax, det_Cm)  !eq16-YF2011
  call get_HtH_Htd_dtd(kmax, lmax, H, d, Cm, HtH, Htd, dtd)
  deallocate(Cm,obs_m,grn_m,H,d)
  allocate(HtH_r(lmax,lmax))
  !$omp parallel do   
  do l = 1, lmax     
    call scopy(lmax, HtH(1,l) , 1, HtH_r(1,l) , 1)   
    call saxpy(lmax, beta1, G1(1,l),1, HtH_r(1,l),1) 
  enddo 
  !$omp end parallel do 
  deallocate(G1)
  !-----------
  !Reduce the effect of errors that occur when calculating with single precision, which can have a growing effect on each iteration.
  allocate(d_HtH_r(lmax,lmax),d_Htd(lmax),d_a(lmax))
  !$omp parallel do   
  do l = 1, lmax
    d_HtH_r(1:lmax,l) = dble(HtH_r(1:lmax,l))
    d_Htd(l) = dble(Htd(l))
  enddo
  !$omp end parallel do 
  call get_sol_posv_d   (lmax, d_HtH_r, lmax, d_Htd, d_a, d_det_HtH_G)  !eq29-YF2011
  a = real(d_a)
  !-----------
  call get_residual (lmax, HtH_r, Htd, dtd, a, s)           !eq28-YF2011
  write(6,*)beta1,alpha2,s
  allocate(siga(lmax))
  siga = 0.
  beta = 0.
  if(itera>0) call get_dump_sol(a,lmax,dump)  
  open(16,file="x.vector",form='unformatted') 
  write(16)(a(l),l=1,lmax) 
  close(16)
    call writesol(mn,nn,jtn,jtn0,l_id_m,icmn,a,tr,0.,0.,strike,dip,rake,    & 
          m0,n0,xx,yy,rtime,beta,siga,depth,stcd,comp,sigw,dt,   & 
          cp,jn,vr0,nsurface,rslip,ylength,alpha1,alpha2) 
  stop
END Program ABIC_TA
