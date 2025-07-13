! S. Yamashita & Y. Yagi
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_G12_matrix(mn,nn,jtn,icmn,jtn0,vr0,xx,yy,tr,l_id_m,rtime,st_max,r_s_t,itera,cr,nsurface,nsurface_o,G1,lmax)
  implicit none
  integer,intent(in) :: mn, nn, jtn0, icmn, itera, lmax, nsurface, nsurface_o
  integer,intent(in) :: jtn(mn,nn), l_id_m(mn,nn,jtn0,icmn)
  real,intent(in)    :: xx,yy,rtime,vr0,st_max, r_s_t, cr
  integer,intent(in) :: tr(mn,nn)
  real,intent(out)   :: G1(lmax,lmax)
  !--------------
  integer :: ns,l,m,n,ios
  integer :: icm,jt
  real,allocatable :: Wwm1(:,:), G2(:,:)
  real,allocatable :: weight(:,:,:,:), tr_r(:,:)
  integer,allocatable :: nfault(:,:)
  !-------------
  ns = 0
  if(nsurface_o == 1 .and. nsurface == 1 ) ns  = 1
  if(nsurface_o == 2 ) ns  = 1
  allocate(weight(mn,nn,jtn0,icmn))
  call get_tmp_weight(itera, cr, weight, mn, nn, jtn0, icmn, rtime)
  allocate(Wwm1(lmax,lmax),G2(lmax,lmax))
  !-------------
  Wwm1 = 0. 
  call tmp_get_lap_time(mn,nn,jtn,jtn0,l_id_m,icmn,1.0,0,Wwm1,lmax,rtime,weight)
  call ssyrk('u','t',lmax,lmax,1.,Wwm1,lmax,0.,G1,lmax)
  !-------------
  Wwm1 = 0. 
  allocate(tr_r(mn,nn),nfault(mn,nn))
  open(11, file = 'TR_fault_segments.txt', status = 'old', iostat = ios)
  if (ios == 0) then
    do n = nn, 1, -1
      read(11,*, iostat=ios)  (tr_r(m, n), m = 1, mn)
      if (ios < 0) goto 98
    end do
    read(11,*)
    do n = nn, 1, -1
      read(11,*, iostat=ios)  (nfault(m, n), m = 1, mn)
      if (ios < 0) goto 98
    end do
    write(6,*) "*** TR_fault_segmetns.txt was loaded properly. -> Multi-Fault Mode ***"
    call tmp_get_lap_space_fs(ns,xx,yy,mn,nn,jtn,jtn0,l_id_m,rtime,icmn,tr,r_s_t,1.0,0,Wwm1,lmax,weight,nfault)
    goto 99
  endif
98 call tmp_get_lap_space(ns,xx,yy,mn,nn,jtn,jtn0,l_id_m,rtime,icmn,tr,r_s_t,1.0,0,Wwm1,lmax,weight)
99 call ssyrk('u','t',lmax,lmax,1.,Wwm1,lmax,0.,G2,lmax)
  deallocate(tr_r,nfault)
  close(11)
  deallocate(Wwm1)
  !$omp parallel do
  do l=1,lmax
    call SAXPY( lmax, 1.0, G2(1,l), 1, G1(1,l), 1 )
  end do
  !$omp end parallel do
  deallocate(G2,weight)
  !----- Fixing bugs that do not affect results (2025/03/13) ---
  !$omp parallel do
  do l=2,lmax
    G1(l,1:l-1) = G1(1:l-1,l) 
  enddo
  !$omp end parallel do
  !-------------------------------------------------------------
!--------TEST DUMPING--------
!  do m=1,mn ;do n=1,nn ;do icm = 2, icmn; do jt = 1, jtn(m,n)
!    l = l_id_m(m,n,jt,icm)
!    G1(l,l) = G1(l,l)  + 1.e-5
!  enddo; enddo; enddo; enddo
!--------TEST DUMPING--------
!----------------------------
  return
end subroutine get_G12_matrix
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tmp_get_lap_time(mn,nn,jtn,jtn0,l_id_m,icmn,beta,ndata,a,kmax,rtime,w)
  !
  !     (d^2/dt^2)slip-rate(x,y,t) : Laplasian
  !
  implicit none
  real,parameter :: amagic = 1.0 
  integer,intent(in)  :: mn,nn,jtn0,icmn,ndata,kmax
  integer,intent(in)  :: jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn)
  real,intent(in)     :: beta,rtime
  real,intent(in)     :: w(mn,nn,jtn0,icmn)
  real,intent(out)    :: a(kmax,*)
  !  
  integer :: m,n,jt,icm,k,l
  integer :: k_id_c
  real    :: work1,sigc
  !
  integer :: i_mn
  integer :: m_f(mn*nn*icmn),n_f(mn*nn*icmn),icm_f(mn*nn*icmn) 
  !- Preparation for efficient parallelization  
  do icm=1,icmn; do n=1,nn; do m=1,mn          
    i_mn =(icm-1)*nn*mn + (n-1)*mn + m        
    m_f(i_mn) = m                            
    n_f(i_mn) = n                           
    icm_f(i_mn) = icm                      
  enddo; enddo; enddo 
  !-----
  sigc = 1.
  work1=beta * (amagic/rtime)**2
  !------
  !$omp parallel 
  !$omp do private(m,n,icm,jt)
  do i_mn = 1,mn*nn*icmn
    m = m_f(i_mn) 
    n = n_f(i_mn)    
    icm = icm_f(i_mn)
    do jt=1,jtn(m,n)
      if( jt > 1) then
        a(l_id_m(m,n,jt,icm),l_id_m(m,n,jt-1,icm)) = - work1 * w(m,n,jt,icm) 
      endif
      if( jt < jtn(m,n) ) then
        a(l_id_m(m,n,jt,icm),l_id_m(m,n,jt+1,icm)) = - work1 * w(m,n,jt,icm)
      endif
      a(l_id_m(m,n,jt,icm),l_id_m(m,n,jt,icm)) = 2.* work1 * w(m,n,jt,icm) 
    enddo
  enddo
  !$omp end do 
  !$omp end parallel 
  return
END subroutine tmp_get_lap_time
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tmp_get_lap_space(ns,xx,yy,mn,nn,jtn,jtn0,l_id_m,rtime,icmn,tr,r_s_t,beta,ndata,a,kmax,w)
  implicit none
  integer,intent(in)  :: ns,mn,nn,jtn0,icmn,ndata,kmax
  integer,intent(in)  :: jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn)
  real,intent(in)     :: xx,yy,rtime,r_s_t,beta
  real,intent(in)     :: tr(mn,nn)
  real, intent(in)    :: w(mn,nn,jtn0,icmn)
  real,intent(out)    :: a(kmax,*)
  integer :: m,n,jt,jt2,icm,k,l
  real    :: wt1,wt2,wk,tw
  integer :: k_id_c
  real    :: disp_amp
  !-------
  integer :: i_mn
  integer :: m_f(mn*nn*icmn),n_f(mn*nn*icmn),icm_f(mn*nn*icmn)
  !- Preparation for efficient parallelization  
  do icm=1,icmn; do n=1,nn; do m=1,mn
    i_mn =(icm-1)*nn*mn + (n-1)*mn + m
    m_f(i_mn) = m
    n_f(i_mn) = n
    icm_f(i_mn) = icm
  enddo; enddo; enddo
  !-----
  !
  !  Spatial Constraint: Laplasian
  write(6,'(" 1/x = r2(space)/r1(time) =   ",f9.2,"   km/sec")')  r_s_t
  wt1 = beta * (r_s_t/xx)**2
  wt2 = beta * (r_s_t/yy)**2
  !$omp parallel 
  !$omp do private(m,n,icm,wk,jt,k,l,tw,jt2)
  do i_mn = 1,mn*nn*icmn
    m = m_f(i_mn)
    n = n_f(i_mn)
    icm = icm_f(i_mn)
    wk = 2.  !; if(ns == 1  .and.  n==nn) wk=1.
    !---
    if(ns == 1 ) then
      if(n==nn) then
        wk = wk - 1.
      else
        if(jtn(m,n+1) == 0) wk = wk - 1.
      endif
    endif
    !---
    do jt = 1,jtn(m,n)
      k = l_id_m(m,n,jt,icm)
      l = l_id_m(m,n,jt,icm)
      a(k,l) = -2 * wt1 - wk * wt2
      a(k,l) = a(k,l) * w(m,n,jt,icm) ! add for adaptive smoothing
!--------
      tw = tr(m,n) + rtime * jt
!--------
!        tw = tr(m,n) + rtime * jt + xx / r_s_t
!--------
      if(m > 1) then
        if(jtn(m-1,n) .ne. 0) then 
          do jt2 = 1,jtn(m-1,n)
            l = l_id_m(m-1,n,jt2,icm)
            a(k,l) =  wt1 *  disp_amp(jt2,rtime,tr(m-1,n),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
          enddo
        else
          if (n > 1) then
            if(jtn(m-1,n-1) .ne. 0) then
              do jt2 = 1,jtn(m-1,n-1)
                l = l_id_m(m-1,n-1,jt2,icm)
                a(k,l) =  wt1 *  disp_amp(jt2,rtime,tr(m-1,n-1),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
              enddo
            endif
          endif
        endif
      endif
      if(m < mn) then
        if(jtn(m+1,n) .ne. 0) then
          do jt2 = 1,jtn(m+1,n)
            l = l_id_m(m+1,n,jt2,icm)
            a(k,l) =  wt1 *  disp_amp(jt2,rtime,tr(m+1,n),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
          enddo
        else
          if (n > 1) then
            if(jtn(m+1,n-1) .ne. 0) then
              do jt2 = 1,jtn(m+1,n-1)
                l = l_id_m(m+1,n-1,jt2,icm)
                a(k,l) =  wt1 *  disp_amp(jt2,rtime,tr(m+1,n-1),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
              enddo
            endif
          endif
        endif
      endif
      if(n > 1) then
        do jt2 = 1,jtn(m,n-1)
          l = l_id_m(m,n-1,jt2,icm)
          a(k,l) =  wt2 *  disp_amp(jt2,rtime,tr(m,n-1),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
        enddo
      endif
      if(n < nn) then
        do jt2 = 1,jtn(m,n+1)
          l = l_id_m(m,n+1,jt2,icm)
          a(k,l) =  wt2 *  disp_amp(jt2,rtime,tr(m,n+1),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
        enddo
      endif
    enddo            !Loop jt   (time1)
  enddo        
  !$omp end do 
  !$omp end parallel 
  return
end subroutine tmp_get_lap_space
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tmp_get_lap_space_fs(ns,xx,yy,mn,nn,jtn,jtn0,l_id_m,rtime,icmn,tr,r_s_t,beta,ndata,a,kmax,w,nfault)
  implicit none
  integer,intent(in)  :: ns,mn,nn,jtn0,icmn,ndata,kmax
  integer,intent(in)  :: jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn),nfault(mn,nn)
  real,intent(in)     :: xx,yy,rtime,r_s_t,beta
  real,intent(in)     :: tr(mn,nn)
  real, intent(in)    :: w(mn,nn,jtn0,icmn)
  real,intent(out)    :: a(kmax,*)
  integer :: m,n,jt,jt2,icm,k,l
  real    :: wt1,wt2,wk,wk1,wk2,tw
  integer :: k_id_c
  real    :: disp_amp
  !-------
  integer :: i_mn
  integer :: m_f(mn*nn*icmn),n_f(mn*nn*icmn),icm_f(mn*nn*icmn)
  !- Preparation for efficient parallelization  
  do icm=1,icmn; do n=1,nn; do m=1,mn
    i_mn =(icm-1)*nn*mn + (n-1)*mn + m
    m_f(i_mn) = m
    n_f(i_mn) = n
    icm_f(i_mn) = icm
  enddo; enddo; enddo
  !-----
  !
  !  Spatial Constraint: Laplasian
  write(6,'(" 1/x = r2(space)/r1(time) =   ",f9.2,"   km/sec")')  r_s_t
  wt1 = beta * (r_s_t/xx)**2
  wt2 = beta * (r_s_t/yy)**2
  !$omp parallel 
  !$omp do private(m,n,icm,wk1,wk2,jt,k,l,tw,jt2)
  do i_mn = 1,mn*nn*icmn
    m = m_f(i_mn)
    n = n_f(i_mn)
    icm = icm_f(i_mn)
    wk1 = 2. 
    wk2 = 2. 
    if(ns == 1 ) then
      if(n==nn) then
        wk2 = wk2 - 1.
      else
        if(jtn(m,n+1) == 0) wk2 = wk2 - 1.
      endif
!      if(m > 1 .and. m< mn) then
!         if(jtn(m-1,n) == 0) wk1 = wk1 - 1.
!         if(jtn(m+1,n) == 0) wk1 = wk1 - 1.
!      endif
    endif
    do jt = 1,jtn(m,n)
      k = l_id_m(m,n,jt,icm)
      l = l_id_m(m,n,jt,icm)
      a(k,l) = -wk1 * wt1 - wk2 * wt2
      a(k,l) = a(k,l) * w(m,n,jt,icm) ! add for adaptive smoothing
!--------
      tw = tr(m,n) + rtime * jt
!--------
!      tw = tr(m,n) + rtime * jt + xx / r_s_t
!--------
      if(m > 1) then
        if( nfault(m-1,n) == nfault(m,n) ) then
          do jt2 = 1,jtn(m-1,n)
            l = l_id_m(m-1,n,jt2,icm)
            a(k,l) =  wt1 *  disp_amp(jt2,rtime,tr(m-1,n),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
          enddo
        endif
      endif
      if(m < mn) then
        if( nfault(m+1,n) == nfault(m,n) ) then
          do jt2 = 1,jtn(m+1,n)
            l = l_id_m(m+1,n,jt2,icm)
            a(k,l) =  wt1 *  disp_amp(jt2,rtime,tr(m+1,n),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
          enddo
        endif
      endif
      if(n > 1) then
        if( nfault(m,n-1) == nfault(m,n) ) then
          do jt2 = 1,jtn(m,n-1)
            l = l_id_m(m,n-1,jt2,icm)
            a(k,l) =  wt2 *  disp_amp(jt2,rtime,tr(m,n-1),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
          enddo
        endif
      endif
      if(n < nn) then
        if( nfault(m,n+1) == nfault(m,n) ) then
          do jt2 = 1,jtn(m,n+1)
            l = l_id_m(m,n+1,jt2,icm)
            a(k,l) =  wt2 *  disp_amp(jt2,rtime,tr(m,n+1),tw) * w(m,n,jt,icm) ! add for adaptive smoothing
          enddo
        endif
      endif
    enddo            !Loop jt   (time1)
  enddo
  !$omp end do 
  !$omp end parallel 
  return
end subroutine tmp_get_lap_space_fs
!++++++++++++++++++++++++++++++++++++++++++++++++
function disp_amp(jt,rtime,ts,t)
  implicit none
  real :: disp_amp
  integer,intent(in) :: jt
  real,intent(in) :: rtime,ts,t
  if(t <= ts + (jt-1)*rtime ) then
    disp_amp = 0.
    return
  endif
  if(t >= ts + (jt+1)*rtime ) then
    disp_amp = 0.
    return
  endif
  if(t <= ts + (jt)*rtime ) then
    disp_amp = (t-ts-(jt-1)*rtime)/rtime
    return
  endif
  if(t > ts + (jt)*rtime ) then
    disp_amp = 1.- (t-ts-(jt)*rtime)/rtime
    return
  endif
end function disp_amp
!++++++++++++++++++++++++++++++++++++++++++++++++
