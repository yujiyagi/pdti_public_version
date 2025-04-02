program GreenpointSource
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  Green's function calculation program (Y. Yagi)
  !  Developed to be able to calculate Green's function for each space-knot
  !  to adapt for multiple fault segment by introduting the knot.dat (2018.03 smz)
  !  Modified to calculate the geometrical spreading factor using the TauP ray 
  !  parameter of the ak135 model. (2024.07.01 Y. Yagi)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none
  integer,parameter :: nl0 = 50, ncal = 8 
  real,parameter    :: pi = 3.141593, rad = .0174533
  character stcd * 8, comp * 4, cha * 1 
  real, allocatable :: gf(:, :, :), y(:), depth(:), rigid(:), dly(:, :), ww(:)
  complex, allocatable :: zi(:), zQp(:), zqs(:)
  complex                 zz (50), zp0 (50)
  real                    slip1 (6)
  real(8)                 gfact (500), pp0 (500)
  real                    vp(nl0),  vs(nl0),  den(nl0),  dep(nl0) ! near-source
  real                    vp1(nl0), vs1(nl0), den1(nl0), dep1(nl0)! near-station
  real                    vp2(nl0), vs2(nl0), den2(nl0), dep2(nl0)! for PP-wave
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  real, allocatable :: sstk(:, :), sdip(:, :), srake(:, :), lat(:, :), lon(:, :), sdep(:, :), amp(:, :)
  real    :: stktmp, diptmp, raketmp, lattmp, lontmp, depthtmp
  real    :: b_strike(1:6), b_dip(1:6), b_rake(1:6)
  real    :: tQp,tqs,epilat,epilon,stk,dip,slip0,slip,vr_i,vr0
  real    :: slat, slon, az, del, delta, az2, fc, p, tp0
  !----------
  integer :: i,l,k,icm
  integer :: nl,nl1,nl2,nt,nt0,nd0,ll,nq
  integer :: ktag,ltag,ios,nline
  integer :: ic1,ic2,i1,iwork,ier,nwk,nwk1
  real    :: dt,df,dw,f,dtq,ds
  real    :: vps,vss,vhypo,vobs,vp_p
  real    :: wk,wk1,wk2,fwork,rigidm,gf_var,gf_var2,td_p,gmax
  !----------
  integer :: nlen,nk,l0,k0,icmn,nflag,nlen_t,nk_t,l0_t,k0_t,icmn_t,nflag_t
  real    :: dl,dk,h0,rslip,ylength,dl_t,dk_t,h0_t,rslip_t,ylength_t
  !----------
  integer :: lk
  integer,allocatable :: k_f(:),l_f(:)
  !---------- function ----
  integer :: l_source
  real    :: getDelta,getAZ,getTP,getG,getP
  character :: name
  !----------
  real    :: get_TauP
  real :: rayp(110,19), TPtable(110,19)
  call  load_taup_file(rayp, "rayP" )
  call  load_taup_file(TPtable, "tP" )
  !----------  2024.07.01 Yagi
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  open(11, file="Rotation_basis_DC.dat", status="old", err = 888 )
  read(11,*, err=888)(b_strike(i),i=1,6)
  read(11,*, err=888)(b_dip(i),i=1,6)
  read(11,*, err=888)(b_rake(i),i=1,6)
  goto 887
888 continue
  b_strike(1) = 0.
  b_strike(2) = 135.
  b_strike(3) = 180.
  b_strike(4) = 90.
  b_strike(5) = 90.
  b_strike(6) = 0.
  b_dip(1) = 90.
  b_dip(2) = 90.
  b_dip(3) = 90.
  b_dip(4) = 90.
  b_dip(5) = 45.
  b_dip(6) = 500.
  b_rake(1) = 0.
  b_rake(2) = 0.
  b_rake(3) = 90.
  b_rake(4) = 90.
  b_rake(5) = 90.
  b_rake(6) = 0.
887 close(11)
  write(6,*) " basis double couple component. "
  write(6,'("strike: ", 6f7.1)') ( b_strike(i),i=1,6)
  write(6,'("dip:    ", 6f7.1)') ( b_dip(i),i=1,6)
  write(6,'("rake:   ", 6f7.1)') ( b_rake(i),i=1,6)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  open (2, file = "structure.dat") 
  read (2, '(a40)') name 
  read (2, * ) tQp, tqs, nl, (vp (l), vs (l), den (l), dep (l), wk1, wk2, l = 1, nl)
  read (2, * ) nl1, (vp1 (l), vs1 (l), den1 (l), dep1 (l), l = 1,  nl1) 
  read (2, * ) nl2, (vp2 (l), vs2 (l), den2 (l), dep2 (l), l = 1,  nl2) 
  close (2) 
  open (3, file = "epicenter.dat") 
  read (3, * ) epilat, epilon
  close(3)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! < from unit 5 >
  !         nt should be powers of 2.
  read (5, * ) nt0, dt, stk, dip, slip0, h0, vr_i
  vr0 = vr_i
  read (5, * ) dl, dk, nlen, nk, l0, k0, icmn, rslip, nflag
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  open (20, file = "fault.dat_in")
  read(20, *, iostat = ios)dl_t, dk_t, nlen_t, nk_t, l0_t, k0_t, icmn_t, h0_t ,rslip_t, nflag_t, ylength_t
  close(20)
  open (10, file = "fault.dat")
  if (ios == 0) then
     ylength = ylength_t
     nk = nk_t
     write (10, '(f8.2,1x,f8.2,1x,4(i3,1x),i2,1x,f5.1,1x,f10.2,i3,1x,f7.2,3f9.3)')  &
          dl_t, dk_t, nlen_t, nk, l0_t, k0_t, icmn_t, h0_t ,rslip_t, nflag_t, ylength, stk, dip, slip0
  else
     ylength = dk
     write (10, '(f8.2,1x,f8.2,1x,4(i3,1x),i2,1x,f5.1,1x,f10.2,i3,1x,f7.2,3f9.3)')  &
          dl, dk, nlen, nk, l0, k0, icmn, h0 ,rslip, nflag, ylength, stk, dip, slip0
  endif
  close (10)
  allocate(amp(1:nlen, 1:nk))
  !When multiple faults are included, the value of amp varies depending on the fault geometory.
  amp = 1.
  amp(1:nlen,nk) = 1. - 0.5 * ((dk-ylength)**2 / dk **2) ! Relative area of the upper knot
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allocate(sstk(1:nlen, 1:nk), sdip(1:nlen, 1:nk), srake(1:nlen, 1:nk))
  allocate(sdep(1:nlen, 1:nk), dly(1:nlen, 1:nk), lat(1:nlen, 1:nk), lon(1:nlen, 1:nk))
  allocate(k_f(nlen*nk),l_f(nlen*nk))
  open(430,file='knot.dat_in',status='old')
  nline = 0
  do
     read(430, *, iostat = ios)
     if (ios == 0) then
        nline = nline + 1
     elseif (ios < 0) then
        exit
     end if
  end do
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (nline /= nk * nlen) then
     write(6, *)" Stop. Please check the number of space knots"
     stop
  endif
  rewind(430)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do k = 1, nk
     do l = 1, nlen
        read(430, *, iostat = ios)ktag, ltag, lattmp, lontmp, depthtmp, stktmp, diptmp, raketmp
        if (ios == 0) then
           lat(ltag, ktag)   = lattmp
           lon(ltag, ktag)   = lontmp
           sdep(ltag, ktag)  = depthtmp
           sstk(ltag, ktag)  = stktmp
           sdip(ltag, ktag)  = diptmp
           srake(ltag, ktag) = raketmp
        else
           write(6, *) " Stop. Failed to read knot.dat_in"
           stop
        end if
     enddo
  enddo
  close(430)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slip = slip0 + 45.
  rslip = 45.
  slip1 = 0.
  if (icmn == 5) then
     slip1(1:icmn) = b_rake(1:icmn)
  else if (icmn == 2) then
     slip1 (1) = slip - rslip
     slip1 (2) = slip + rslip
  else
     slip1 (1) = slip0
  end if
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nl = l_source(h0+100., dep, nl) ! # of layers in structure.dat
  dep(nl) = 0.
  !---   set dimension
  nt = int(nt0 * ncal)
  nd0 = nt  + int(30. / dt)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allocate (gf(1:nd0, 1:nlen, 1:nk), depth(nk + 1), ww(nd0),  &
       rigid(nk + 1), zi(nd0), zQp(nt), zqs(nt), stat=ier)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  open (96, file = 'rigid.dat') 
  rigidm = 0.
  nwk = 0
  do k = 1, nk
     !do l = 1, nlen
     depth(k) = sum(sdep(1:nlen, k))/real(nlen) ! assuming the average depth at k
     wk = 0.0
     if (vs(1) == 0.) wk = dep(1) + 0.0
     if(depth(k) < wk )  then
        nwk = max(nwk, k)
        cycle
     endif
     ll = l_source(depth(k), dep, nl)
     if(ll == 0 ) then
        rigid (k) = rigid (k + 1) 
     else
        rigid (k) = vs (ll) * vs (ll) * den (ll)
        vp_p = vp(ll)
     endif
     write(96,'(3f10.2,i5,f10.2)') rigid(k), sum(amp(1:nlen,k))/real(nlen), depth(k), 1, vp_p
     rigidm = rigid (k) + rigidm 
     !enddo
  enddo
  close(96)
  rigidm = rigidm / nk
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  open (97, file = 'rigid_amp.info') 
  do k=1,nk
     do l = 1,nlen
        ll =  l_source(sdep(l,k), dep, nl)
        write(97,'(2i6,3f12.5)') k, l, vs (ll) * vs (ll) * den (ll), amp(l,k), sdep(l,k)
     enddo
  enddo
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  df = 1. / (dt * nt)
  dw = df * 2. * pi
  ! < q-response >
  zQp = 0.
  zqs = 0.
  call qf (zQp, nt, tQp, df) 
  call qf (zqs, nt, tqs, df) 
  !-----  diff   
  ! To simplify the covariance matrix, velocity waveforms are used for inversion.
  do i = 2,nt/2
     f = df * (i - 1)
     zQp(i) = zQp(i) * cmplx(0,2. * 3.1415926 * f)
     zQp(nt - i + 2) = conjg(zQp(i))
     zqs(i) = zqs(i) * cmplx(0,2. * 3.1415926 * f)
     zqs(nt - i + 2) = conjg(zqs(i))
  enddo
  !-----  diff
  ! < source layer ll >
  ll = l_source(h0, dep, nl)
  vps = vp (ll) 
  vss = vs (ll) 
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call instg (zi, nt, dw, 0, zp0, zz, 0, 0, 1., 0) 
  ds = 1.e4
  call system("mkdir -p wave.grn")
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  open (44, file = "station.list") 
100 read (44, *, end = 999) stcd, comp, slat, slon, az, del 
  do icm = 1, icmn
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     do k = 1, nk
        do l = 1, nlen
           if(srake(l, k) .le. -180.) srake(l, k) =   360. + srake(l, k)
           if(srake(l, k) .ge.  180.) srake(l, k) = - 360. + srake(l, k)
        enddo
     enddo
     if (icmn .eq. 5) then
        sstk(1:nlen, 1:nk)  = b_strike(icm)
        sdip(1:nlen, 1:nk)  = b_dip(icm)
        srake(1:nlen, 1:nk) = b_rake(icm)
     end if
     !----- get time lag
     delta = getDelta(epilat, epilon, slat, slon) !epicentral distance of the station
     az2 = getAZ(epilat, epilon, slat, slon)      !azimuth from epicenter
!     tp0 = getTP(delta, h0)                 !travel time from epi. to station
     tp0 = get_TauP(delta, h0, TPtable,"tP")   ! 2024.07.01 Yagi
     do k = 1, nk
        do l = 1, nlen
           delta = getDelta(lat(l, k), lon(l, k), slat, slon) !distan. from eachknot to the station
           if(sdep(l, k) > 0. )then
!              dly(l, k) =  getTP(delta, sdep(l,k)) - tp0   !difference of the travel time between
              dly(l, k) =  get_TauP(delta, sdep(l,k), TPtable,"tP")  - tp0   ! 2024.07.01 Yagi
           else                                   ! epicenter-station and each-knot-station
!              dly(l, k) =  getTP(delta, 0.) - tp0
              dly(l, k) =   get_TauP(delta, 0., TPtable,"tP") - tp0   ! 2024.07.01 Yagi
           endif
        enddo
     enddo
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     fwork = -minval(dly) + 10.
     iwork = nint(fwork / dt)
     dly = dly + fwork
     gf = 0. ! green's function at each knot
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     do k=1,nk
       do l=1,nlen
         lk = (k-1)*nlen + l
         k_f(lk) = k
         l_f(lk) = l
!         write(6,*) lk,k,l
       enddo
     enddo
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     !$omp parallel 
     !$omp do private(k,l,delta,az2,fc,p,ll,ww,gfact,pp0)
     do lk = 1, nlen*nk
        k = k_f(lk)
        l = l_f(lk)
        delta = getDelta(lat(l,k), lon(l,k), slat, slon)
        az2 = getAZ(lat(l,k), lon(l,k), slat, slon) 
!---
!        call geom_m (gfact, pp0, dble(sdep(l,k)), vhypo, vobs)
!        fc = getG (delta, gfact) * ds
!        p = getp (delta, pp0)
        fc = get_TauP(delta, sdep(l,k), rayp,"G") * ds
        p  = get_TauP(delta, sdep(l,k), rayp,"rayP")
!---  2024.07.01 Y. Yagi
        call bodyw (ww, nt, dt, 1, 1, sstk(l, k), sdip(l, k), srake(l, k), sdep(l, k), &
             az2, p, fc, zQp, zQs, zi, nl, vp, vs, den, dep, &
             nl1, vp1, vs1, den1, dep1, nl2, vp2, vs2, den2, dep2)
        ll = l_source(sdep(l, k), dep, nl)
        ww = ww  * vs (ll) * vs (ll) * den (ll) * amp(l, k)
        call resample_shift(ww,nt,dt,gf(1, l, k),nt,dt, dly(l, k)+dt*0.5) !(+dt*0.5)P-wave pick position is just after increasing and for improved stability (2024/06/13)
     enddo  
     !$omp end do
     !$omp end parallel
     !------------------
     !
     !    - - out put green function - -
     !
     ic1 = index(stcd, " ") - 1
     ic2 = index(comp, " ") - 1
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! open (99, file = "wave.grn_syn/"//stcd(1:ic1)//comp(1:ic2)//cha(icm),form='unformatted')
     open (99, file = "wave.grn/"//stcd(1:ic1)//comp(1:ic2)//cha(icm),form='unformatted')
     write (99) stcd(1:ic1) , comp(1:ic2)
     write (99) dl,  dk,  nlen,  nk,  l0,  k0,  icmn,  h0
     write (99) stk, dip,  slip1(icm),  rigidm
     write (99) -fwork,  dt,  0.,  nt/ncal + iwork,  az,  del,  slat,  slon
     gf = gf * (dt * dl * dk / 1000.)
     do k = 1, nk
        do l = 1, nlen
           write (99) (gf(i,l,k), i = 1, nt/ncal+iwork)
        enddo
     enddo
     !- - - check - - -
     ! open(77, file = " CheckGF.txt", status='replace')
     ! do n = 1, nk
     !    write(77,'(">-Z",e8.2)') sdep(1,n)
     !    do i = 1, nt/ncal+iwork
     !       write(77,'(i6, f10.5)') i, gf(i,1,n)
     !    enddo
     !    !write(77,'(">")')
     ! end do
     ! close(77)
     !- - - - - - - - - 
     close (99)
     if(icm == 1) then
        open (98, file = "wave.grn/"//stcd(1:ic1)//comp(1:ic2)//'info', &
             form='unformatted')
        do k = 1, nk
           do l = 1, nlen
              write (98)  dly(l, k)
           enddo
        enddo
     endif
     do k = 1, nk
        do l = 1, nlen
           gf_var = 0.
           nwk = nt/ncal + iwork
           i1 = nint(dly(l, k)/dt)
           do i = i1, nwk
              gf_var = gf_var + gf(i, l, k) ** 2
           enddo
           gf_var2 = 0.
           do i = i1, nwk
              gf_var2 = gf_var2 + gf(i, l, k) ** 2
              if(gf_var2 >= gf_var * 0.90) then
                 td_p = i * dt - dly(l, k)
                 nwk1 = i
                 exit
              endif
           enddo
           gf_var = sqrt(gf_var2 / (nwk1 - i1 + 1) )
           gmax = 0.
           do i = 1, nwk
              if(abs(gf(i, l, k)) > gmax ) gmax = abs(gf(i, l, k))
           enddo
           write (98) gmax, td_p
           !----
        enddo
     enddo
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     if(icm == 1) then
        write(98) nt / ncal + iwork
        write (98) (gf(i, 1, nk), i = 1, nt/ncal+iwork)
        write (98) (gf(i, 1, 1), i = 1, nt/ncal+iwork)
        write (98) (gf(i, nlen, nk), i = 1, nt/ncal+iwork)
        write (98) (gf(i, nlen, 1), i = 1, nt/ncal+iwork)
     endif
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     if (icmn == 5) then
        if (icm == 5) then
           close(98)
        end if
     else
        if (icm == 2) then
           close(98)
        end if
     end if
699  format(5(e14.7,1x))
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  enddo
  goto 100
999 continue
  close (44)
  open(97, file = "Qeffect.dat")
  dtq = 0.1
  nq = nt0 * ncal * 4
  deallocate(zQp)
  allocate(zQp(nq), y(nq))
  y = 0.
  df = 1. / (dtq * nq) 
  call Qf (zQp, nq, tQp, df) 
  do i = 2, nq/2
     zQp(nq - i + 2) = conjg(zQp(i))
  enddo
  call cfft (zQp, nq, 1)
  do  i = 1, nq
     y(i) = zQp(i) * df
  enddo
  deallocate(zQp)
  write(97,'(i7,1x,f10.5)') nq/ncal, dtq
  write (97, '(f10.7)') (y(i), i = 1, nq/ncal)
  close(97)
998 stop
end program GreenpointSource
