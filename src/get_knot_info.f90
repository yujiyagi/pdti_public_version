program get_knot_info
  !-----------
  ! output knot.dat_in for calculation of green's function 
  !        fault.dat (The number of knots in the depth direction has been updated.)
  !-----------
  implicit none
  real, allocatable :: lat(:,:), lon(:,:), sdep(:,:)
  integer :: nt0, nlen, nk, l0, k0, icmn, nflag
  integer :: k, l
  real    :: stk, dip, rake, h0, vr_i, dl, dk, dt, rslip
  real    :: wx, wy, hypo_lat, hypo_lon
  real    :: tqp, tqs, nl, vp , vs , den , dep
  real    :: dwork, dpoint, ylength 
  real    :: rad = 0.0174533
  !------------
  read (5, * ) nt0, dt, stk, dip, rake, h0, vr_i
  read (5, * ) dl, dk, nlen, nk, l0, k0, icmn, rslip, nflag
  open (3, file = "epicenter.dat") 
  read (3, * ) hypo_lat,hypo_lon
  close(3)
  !------------
  allocate(lat(nlen,nk),lon(nlen,nk),sdep(nlen,nk))
  open (2, file ="structure.dat")
  read (2, * )
  read (2, * ) tqp, tqs, nl, vp , vs , den , dep
  close(2)
  dwork = 0.
  if (vs .eq. 0.) dwork = dwork + dep
  do k = 1, nk
    if(h0 - (k - k0) * dk * sin (dip * rad)  < dwork ) then
      nk = k -1
      exit
    endif
    do l = 1, nlen
      wx =  dl * (l -l0)
      wy =  dk * (k -k0)
      call xy2geo(hypo_lat,hypo_lon,wx,wy,stk,dip,lat(l,k),lon(l,k))
      sdep(l,k) = h0 - (k - k0) * dk * sin (dip * rad)
    enddo
  enddo
  !------------
  open(20,file='knot.dat_in',status='replace')
  do k = 1, nk
    do l = 1, nlen
       write(20,'(2i5,3f12.5,3f10.2)') k,l,lat(l,k),lon(l,k),sdep(l,k),stk,dip,rake
    enddo
  enddo
  !------------
  dpoint = h0 + real(k0 - nk ) * dk * sin (rad * dip) 
  ylength = (dpoint - dwork)/sin(rad*dip)
  ylength = min(ylength,dk)
  open (10, file = "fault.dat_in") 
  write (10, '(f8.2,1x,f8.2,1x,4(i3,1x),i2,1x,f5.1,1x,f10.2,i3,1x,f7.2)')  &
           dl, dk, nlen, nk, l0, k0, icmn, h0 ,rslip, nflag, ylength
  close (10)
  !------------
  stop
end program get_knot_info
