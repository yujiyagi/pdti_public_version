! S. Yamashita and Y. Yagi

subroutine slip_rate2fort40 (imax, dt, mn, nn, icmn, jtn0, rtime, slip_rate, Tr)
  implicit none
  integer, intent(in) :: imax, mn, nn, icmn, jtn0
  real, intent(in)    :: dt, rtime
  real, intent(in)    :: slip_rate(imax,mn,nn,icmn), Tr(mn,nn)
  integer             :: jt, icm, m, n, i, j, ir
  real                :: dm, xmo(1:6), dum(1:3), disp_length
  real, allocatable   :: stf(:,:,:,:), PDT_comp(:,:,:), PDT_rate(:,:,:,:), slip(:,:)
  integer             :: m0, n0, ns
  real                :: elon, elat, strike, dip, rake_in, depth, vr, xx, yy, rslip
  real, parameter     :: ylen=0.
  !--- get parameters
  open(17, file="epicenter.dat", status="old")
  read(17, *) elat, elon
  close(17)
  open(16, file="i_greenf", status="old")
  read(16, *) dm, dm, strike, dip, rake_in, depth, vr
  read(16, *) xx, yy, dm, dm, m0, n0, dm, rslip, ns
  close(16)
  !--- resample slip_rate from dt to rtime
  ir=nint(imax*dt/rtime)
  allocate(stf(mn,nn,ir,icmn))
  stf=0.
  do jt=1,ir
    i=nint(rtime / dt)
    j=nint(((jt - 1) * rtime) / dt)
    do icm=1,icmn
      do n=1,nn; do m=1,mn
        stf(m,n,jt,icm)=slip_rate(i+j,m,n,icm)
      enddo; enddo
    enddo
  enddo
  !--- shift to front
  allocate(PDT_rate(mn,nn,icmn,jtn0))
  PDT_rate=0.
  do n=1,nn; do m=1,mn
    do jt=1,jtn0
      j=nint(Tr(m,n) / rtime)+jt
      do icm=1,icmn
        PDT_rate(m,n,icm,jt)=stf(m,n,j,icm)
      end do
    end do
  end do; end do
  !--- get total slip
  allocate(PDT_comp(mn,nn,icmn))
  do n=1,nn; do m=1,mn
    do icm=1,icmn
      PDT_comp(m,n,icm)=sum(PDT_rate(m,n,icm,1:jtn0))
    end do
  end do; end do
  allocate(slip(mn,nn))
  do n=1,nn; do m=1,mn
    if (icmn .eq. 5) then
      xmo=0.
      xmo(1:icmn)=PDT_comp(m,n,1:icmn)
      call d_cp(xmo, dm, dm, dm, slip(m, n), dm, dum)
    else if (icmn .eq. 2) then
      slip(m,n)=disp_length(PDT_comp(m,n,1), PDT_comp(m,n,2), rake_in+45., rslip)
    else if (icmn .eq. 1) then
      slip(m,n)=PDT_comp(m,n,1)
    end if
  end do; end do
  !--- generate fort.40
  call writesol(mn, nn, jtn0, icmn, tr, 0., 0., strike, dip, rake_in, rake_in, &
     0., m0, n0, xx, yy, rtime, 0., depth, vr, ns, rslip, ylen, 0., &
     0., 0., elat, elon, slip, PDT_comp, PDT_rate)
  !---
!  open(15, file="original_slip_rate.dat", status="replace")
!  do i=1,imax
!    do icm=1,icmn
!      write(15, *) "JT, ICM: ", i, icm
!      do n=nn,1,-1
!        write(15,'(30(f9.4,x))')  (slip_rate(i,m,n,icm), m=1,mn)
!      end do
!    end do
!  end do
!  close(15)
!  !---
!  open(14, file="resampled_slip_rate.dat", status="replace")
!  do jt=1,ir
!    do icm=1,icmn
!      write(14, *) "JT, ICM: ", jt, icm
!      do n=nn,1,-1
!        write(14,'(30(f9.4,x))')  (stf(m,n,jt,icm), m=1,mn)
!      end do
!    end do
!  end do
!  close(14)
!  !---
!  open(13, file="d_tq_slip_rate.dat", status="replace")
!  do jt=1,ir
!    write(13,'(i4, f7.2, 5f10.4)') jt, jt*rtime, ( sum(stf(1:mn,1:nn,jt,icm)), icm=1,icmn)
!  end do
!  close(13)
!  !---
  deallocate(stf, PDT_comp, PDT_rate, slip)
end subroutine

!------------------------------------------------------------

subroutine  writesol(mn, nn, jtn,icmn,tr,variance,abic,strike,dip,rake,rake_in, &
     rigid, m0, n0, xx, yy, rtime, beta1, depth, vr, nsurface, rslip, ylen,alpha2, &
     moment, Mw, lat, lon, slip, PDT_comp, PDT_rate)
  implicit none
  integer,intent(in) :: mn, nn, m0, n0, jtn, nsurface, icmn
  real,intent(in) :: moment, Mw, rigid, lat, lon, depth, vr, xx, yy, rtime, ylen
  real,intent(in) :: variance, abic, beta1, alpha2, strike, dip, rake, rake_in, rslip
  real,intent(in) :: tr(1:mn,1:nn), slip(1:mn, 1:nn), PDT_comp(1:mn, 1:nn, 1:icmn)
  real,intent(in) :: PDT_rate(1:mn, 1:nn, 1:icmn, 1:jtn)
  integer :: m, n, icm, jt
  real :: slipd(1:mn, 1:nn), disp_slip
  if (icmn .eq. 2) then
     do m = 1, mn
        do n = 1, nn
           slipd(m, n) = disp_slip(PDT_comp(m, n, 1), PDT_comp(m, n, 2), rake_in+45, rslip)
        end do
     end do
  else
     slipd(1:mn, 1:nn) = rake_in
  endif
  open(40, file = 'fort.40', status = 'replace', action = 'write')
  write(40, *)"Moment(Nm),Mw,Rigid(GPa),Lat.,Lon.,Depth(km),Vr(km/sec),nsurface"
  write(40,'(e11.4,1x,f5.2,1x,f7.2,1x,2f10.3,1x,f7.1,1x,e10.3,1x,i2)') &
       moment, Mw, rigid, lat, lon, depth, vr, nsurface
10 Format(" Seismic Moment : ",e11.4," [Nm],  Mw : ",f5.2)
  write(40,*)" Strike, Dip, Slip, slip_input, rslip  "
  write(40,'(6(f7.2,3x))')strike, dip, rake, rake_in, rslip
  write(40,*)"xx,  yy,  mn,  nn,  m0, n0, Raisetime, jtn icmn "
  write(40,'(2(f6.2,x),4(i3,x),f6.2,x,i3,x,i3,x,f7.2)') &
       xx, yy, mn, nn, m0, n0, rtime, jtn, icmn, ylen
  write(40,*)" Varicance, ABIC, beta(1),beta(2),beta(3),alpha1,alpha2 "
  write(40,'(f8.5,x,f13.1,x,3(f11.5,x),2(e12.4,x))') &
       variance, abic, beta1, 0., 0., 0., alpha2
  write(40,*)" Total Slip for each sub-fault "
  if(mn .eq. 1 .and. nn .eq. 1) then
     do n = nn, 1, -1
        write(40,'(e12.3,x)') (slip(m, n)/500., m = 1, mn)
     end do
  else
     do n = nn, 1, -1
        write(40,'(30(f9.4,x))') (slip(m, n), m = 1, mn)
     end do
  endif
  write(40, *)" Slip angel for each sub-fault "
  do n = nn, 1, -1
     write(40,'(30(f7.1,1x))')(slipd(m, n), m = 1, mn)
  end do
  !
  write(40, *)" Start Time for each sub-fault "
  do n = nn, 1, -1
     write(40,'(30(f6.2,1x))')(Tr(m, n), m = 1, mn)
  end do
  write(40,*)" Total Slip Vecoter for each sub-fault "
  do icm = 1, icmn
     write(40, *)" Vector ",icm
     do n = nn, 1, -1
        if(mn.eq.1 .and. nn.eq.1) then
           write(40,'(e12.3)') (PDT_comp(m, n, icm)/500., m = 1, mn)
        else
           write(40,'(30(f9.4,x))') (PDT_comp(m, n, icm), m = 1, mn)
        endif
     end do
  end do
  !
  do jt = 1,jtn
     do icm = 1, icmn
        write(40,*)"JT, ICM : ",jt,",  ",icm
        do n = nn, 1, -1
           if(mn .eq. 1 .and. nn .eq. 1) then
              write(40,'(e12.3)')  (PDT_rate(m, n, icm, jt),m = 1, mn)
           else
              write(40,'(30(f9.4,x))')  (PDT_rate(m, n, icm, jt), m = 1, mn)
           endif
        end do
     end do
  end do
  do jt = 1,jtn
     do icm = 1, icmn
        write(40,*)"Error JT, ICM : ",jt,",  ",icm
        do n = nn, 1, -1
           write(40,'(30(f9.4,x))')  (real(m), m = 1, mn)
        end do
     end do
  end do
  close(40)
  return
END subroutine writesol


