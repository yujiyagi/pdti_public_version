!  S. Yamashita & Y. Yagi
program get_init_model_para_TmpSm
  implicit none
  integer            :: jtn0, mn, nn, m0, n0, icmn
  real               :: rtime, vr0, xx, yy
  real, allocatable  :: tr(:,:), slip_rate(:,:,:,:), so(:)
  integer            :: imax, m, n, jt, nst, i
  real               :: wdt
  !---
  read(5, *)
  read(5, *) rtime, jtn0, vr0
  !---
  open(10, file="fault.dat", status="old")
  read(10,*) xx, yy, mn, nn, m0, n0, icmn
  close(10)
  !---
  allocate(tr(mn, nn))
  call getTR(vr0, xx, yy, mn, nn, m0, n0, tr, rtime)
  !---
  write(6,*) "Inital Slip Model : Uniform slip Model "
  imax = 5000
  wdt  = 0.1
  allocate(slip_rate(imax,mn,nn,icmn),so(imax))
  call stime_m(so,imax,wdt,rtime)
  !slip_rate = 0.
  slip_rate = 0.1e-8
  do m=1,mn ;  do n=1,nn; do jt = 1,jtn0
    nst = nint((Tr(m,n)+real(jt-1)*rtime)/wdt)
    do i=1,min(imax-nst,imax)
      slip_rate(i+nst,m,n,1) = slip_rate(i+nst,m,n,1)+ So(i)*exp(-(jt-1)*5./jtn0)
     ! slip_rate(i+nst,m,n,2) = slip_rate(i+nst,m,n,2)+ So(i)*1.e-8
    enddo
  enddo; enddo; enddo
  deallocate(so)
  ! write(6,*) "Vertical-striped Model"
  ! do m=1,mn
  !   if ( mod(m,2) .eq. 0 ) slip_rate(1:imax,m,1:nn,1:icmn)=0.
  ! end do
  call slip_rate2fort40 (imax,wdt,mn,nn,icmn,jtn0,rtime,slip_rate,tr) ! add for adaptive smoothing
end program
