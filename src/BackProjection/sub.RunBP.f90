subroutine get_nj(nj)
  implicit none
  integer             :: ios
  integer, intent(out):: nj
  open(99, file='station.list', status='old', action='read')
  nj = 0
  do
     read(99, *, iostat = ios)
     if (ios == 0) then
        nj = nj + 1
     else
        exit
     end if
  end do
  close(99)
  return
end subroutine get_nj

subroutine read_station_info(station_code, comp, weight, nj)
  implicit none
  integer, intent(in)            :: nj
  character(len=*), intent(out)  :: station_code(nj),comp(nj)
  real, intent(out)              :: weight(nj)      !station weight (station)
  integer                        :: j, j1, j2, n_wei
  real, allocatable              :: slat(:),slon(:)
  real                           :: delta, getDelta, sweight
  allocate(slat(nj),slon(nj))
  open(99, file='station.list', status='old', action='read')
  do j = 1, nj
     read(99, *) station_code(j), comp(j), slat(j), slon(j)
  end do
  close(99)
  do j1 = 1, nj
     n_wei = 1
     do j2 = 1, nj
        if (j1 == j2 ) cycle
        delta = getDelta(slat(j1), slon(j1), slat(j2), slon(j2))
        if (delta .lt. 20.) n_wei = n_wei + 1
     end do
     weight(j1) = 1.0/n_wei
  end do
  sweight = sum(weight(1:nj))
  weight(1:nj) = weight(1:nj) / sweight 
  deallocate(slat,slon)
  return
end subroutine read_station_info

FUNCTION getDelta (alat1, alon1, alat2, alon2)
  PARAMETER (rad = 0.017453292)
  REAL getDelta, alat1, alon1, alat2, alon2
  t0 = alat1 * rad
  t1 = alat2 * rad
  Dphi = (alon2 - alon1) * rad
  t00 = t0 - 11.55 / 60 * rad * sin (2 * t0)
  t11 = t1 - 11.55 / 60 * rad * sin (2 * t1)
  C00 = COS (T00)
  C11 = COS (T11)
  S00 = SIN (T00)
  S11 = SIN (T11)
  C2 = COS (Dphi)
  d11 = c00 * c11 * c2 + s00 * s11
  getDelta = acos (d11) / rad
  RETURN
END FUNCTION getDelta

function to_radian(degree) result(rad)
  ! degrees to radians
  real,intent(in) :: degree
  real, parameter :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
  real :: rad
  rad = degree*deg_to_rad
end function to_radian
function haversine(deglat1,deglon1,deglat2,deglon2) result (dist)
  ! great circle distance -- adapted from Matlab 
  real,intent(in) :: deglat1,deglon1,deglat2,deglon2
  real :: a,c,dist,dlat,dlon,lat1,lat2
  real,parameter :: radius = 6372.8
  dlat = to_radian(deglat2-deglat1)
  dlon = to_radian(deglon2-deglon1)
  lat1 = to_radian(deglat1)
  lat2 = to_radian(deglat2)
  a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
  c = 2*asin(sqrt(a))
  !  dist = radius*c       !output in km
  dist = c*180/3.1415927 !output in degree
end function haversine

subroutine read_obs(station_code, comp, f1, f2, dt, Fj, alpha, nt, mbpm_flag, potflag)
  implicit none
  integer,intent(in):: nt, mbpm_flag, potflag
  character(len = *), intent(in):: station_code,comp
  real, intent(in) :: f1, f2, dt
  real, intent(out):: Fj(nt),alpha
  real, allocatable:: obs(:)
  integer:: obs_nt, nend
  real:: dto, amp, w 
  integer:: t    !  j: station, t: time
  open(99, file='wave.obs/'//trim(station_code)//trim(comp), status='old', action='read')
  read(99, *) w, dto, obs_nt, w, w, w, nend
  allocate(obs(obs_nt))
  do t = 1, obs_nt
     read(99, *) obs(t)
  enddo
  close(99)
  call diff_wave(obs, dto, obs_nt)
  call rmean(obs, obs_nt)
  call taper_all(obs, obs_nt, 0.025)
  call bp_filter(obs, f1, f2, dto,obs_nt)
  if(mbpm_flag == 0) then
     call resample_shift(obs, obs_nt, dto, Fj(1:nt), nt, dt, -10.)
  else
     call resample_shift(obs, obs_nt, dto, Fj(1:nt), nt, dt, -5.)
  end if
  deallocate(obs)
  nend = nend - nint(5./0.05)
  if (nt > nend) Fj(nend:nt) = 0.   ! Cut PP phase
  alpha = sqrt(sum(Fj(1:nt) **2))
  if (potflag == 1) then
     Fj(1:nt) = Fj(1:nt) !observed waveforms will later be normalized by ampGij (P-phase amp. in GF)
  else if (potflag == 0) then
     Fj(1:nt) = Fj(1:nt) / alpha
  else
     write(6, '("Error! Specify potflag either 0 (OriginalBP/HBP) or 1 (PotencyBP/HBP)")')
     stop
  endif
  return
end subroutine read_obs

subroutine diff_wave(x,dt0,in0)
  real             :: x(in0)
  real, allocatable:: dx(:)
  allocate(dx(in0))
  dx(1) = 0.
  do i=2,in0
     dx(i) = (x(i)-x(i-1))/dt0
  enddo
  do i=1,in0
     x(i) = dx(i)
  enddo
end subroutine diff_wave

subroutine rmean(disp,nterm)
  !       remove mean Wavefrome Program
  !       disp ; waveform
  !       nterm ; number of waveform in points
  dimension disp(nterm)
  pi=3.141593
  dm= 0.
  do i=1,nterm
     dm = dm + disp(i)
  end do
  dm = dm / nterm
  disp = disp - dm
  return
end subroutine rmean

subroutine taper_all(disp,nterm,widtap)
  !       Taper Wavefrome Program
  !       cut Cosain taper
  !       disp ; waveform
  !       nterm ; number of waveform
  !       widtap ; reate of taper in waveform
  dimension disp(nterm)
  pi=3.141593
  m1=ifix(widtap*float(nterm))
  m2=ifix(widtap*float(nterm))
  do i=1,nterm
     if(i.le.1+m1) then
        ca=0.5*(1.0+cos(pi*float(m1+1-i)/float(m1)))
        disp(i)=disp(i)*ca
     endif
     if(i.ge.nterm-m2) then
        cb=0.5*(1.0+cos(pi*float(i-nterm+m2)/float(m2)))
        disp(i)=disp(i)*cb
     endif
  enddo
  return
end subroutine taper_all

subroutine bp_filter(x,f1,f2,dt,nm)
  !     The original code was written by Professor Kikuchi.
  !      https://wwweic.eri.u-tokyo.ac.jp/ETAL/KIKUCHI/
  !     input x (raw data)
  !     output x (filtering data)
  !
  parameter(pai= 3.1415926,np = 6,amp = 20.)
  !parameter(pai= 3.1415926,np = 6,amp = 10.)
  real x(nm)
  complex(8),allocatable ::z(:)
  complex(8):: zc1,zc2
  !
  !f3  = 10.**(alog10(0.5/dtj)-1./np*alog10(amp-1.))
  !write(6,*) dtj,dt,f3
  !
  nn=log(real(nm))/log(2.)
  nn=2**nn
  if(nn.lt.nm) nn=nn*2
  allocate(z(nn))
  z = 0.
  do i=1,nm
     z(i) = x(i)
  enddo
  call cfft(z,nn,-1)
  df = 1./(dt*nn)
  dw = df*pai*2
  do i = 2,NN/2
     f=df*(i-1)
     zc1 =    1./(1.+(f/f2)**np)
     zc2=  1.-1./(1.+(F/f1)**np)
     z(i) = z(i) *zc1*zc2 /nn
     z(nn-i+2)=conjg(z(i))
  enddo
  !  close(99) ! commented out by Ryo Okuwaki, 2017-06-15
  !  Above unnecessary `close(99)` could return error such as
  !  `double free or corruption (out)`
  z(1) = 0.
  z(NN/2+1)=0.
  call cfft(z,nn,1)
  do i=1,nm
     x(i) = z(i)
  enddo
  deallocate(z)
  return
END subroutine bp_filter

SUBROUTINE CFFT(X,N,ID)
  !     The original code was written by Professor Kikuchi.
  !      https://wwweic.eri.u-tokyo.ac.jp/ETAL/KIKUCHI/
  ! <  FFT for complex variables >
  IMPLICIT COMPLEX(8) (X-Z)
  DIMENSION:: X(N)
  PI=SIGN(3.141593,ID*1.)
  N2=N/2
  L=N
1 NI=L/2
  Z=CEXP(CMPLX(0.,PI/NI))
  ZK=1
  DO  K=1,NI
     DO  M=1,N/L
        J=K+2*(M-1)*NI
        XJ=X(J)
        X(J)=XJ+X(J+NI)
        X(J+NI)=(XJ-X(J+NI))*ZK
     enddo
     ZK=ZK*Z
  enddo
  L=NI
  IF(L.GT.1) GOTO 1
  DO I=1,N
     J=I-1
     M=N2
     L=1
5    L=L+M*(J-2*(J/2))
     M=M/2
     J=J/2
     IF(J.GE.1) GOTO 5
     IF(I.GT.L) GOTO 10
     XI=X(I)
     X(I)=X(L)
     X(L)=XI
10   CONTINUE
  enddo
END SUBROUTINE CFFT

subroutine resample_shift(w1,imax1,dt1,w2,imax2,dt2,shift)
  implicit none
  integer, intent(in) :: imax1, imax2
  real,    intent(in) :: dt1,   dt2
  real,    intent(in) :: shift
  real,    intent(in) :: w1(imax1)
  real,    intent(out):: w2(imax2)
  !-----
  real :: trend, t1, t2
  integer :: i, i1
  w2 = 0.
  do i=1,imax2
     t2 = dt2*(i-1) - shift
     i1 = int(t2/dt1) + 1
     t1 = dt1*(i1-1)
     if(i1 <= 0  .or.  i1 >= imax1) cycle
     if(t2 < 0.) then
        trend = (w1(1) - 0. )/dt1
        w2(i) =  w1(1) + trend*(t2)
     else
        trend = (w1(i1+1)-w1(i1))/dt1
        w2(i) = w1(i1) + trend*(t2-t1)
     endif
     !   write(6,'(2f10.3,2f10.3,1i4,f10.3)') trend, w2(i),t1,t2,i1,(i-1)*dt2
  end do
  return
end subroutine resample_shift

subroutine read_green (station_code,comp,f1,f2,TempGij,Tij,polarity,nt,nx,ny,ampGij,bpflag,potflag,acvelflag)
  implicit none
  integer, intent(in) :: nt, nx, ny
  character(len = *), intent(in):: station_code,comp
  integer, intent(in):: acvelflag, bpflag, potflag
  real, intent(in) :: f1, f2
  real, intent(out):: TempGij(nt,nx,ny), ampGij(nx,ny)
  real, intent(out):: polarity
  integer, intent(out):: Tij(nx,ny)
  real :: dtg, w, minT, maxT, peakmax, peakmin, SearchWindow, tmpamp
  integer :: t, x, y
  integer :: green_nt, tmpt
  character :: cwk*9
  if (bpflag==0 .and. potflag==1) then
     open(99, file='wave.grn_withoutQ/'//trim(station_code)//trim(comp)//"1", &
          status='old',form='unformatted', action='read')
  else
     open(99, file='wave.grn/'//trim(station_code)//trim(comp)//"1", &
          status='old',form='unformatted', action='read')
  end if
  read(99)
  read(99)
  read(99) 
  read(99) w, dtg, w, green_nt
  do y=ny, 1, -1; do x=1,nx
     read(99) ( TempGij(t,x,y), t=1,green_nt)
  enddo; enddo
  close(99)
  do x=1,nx ; do y=1,ny
     t=1
     do
        if ( TempGij(t,x,y) .ne. 0 ) then
           Tij(x,y)=t ; exit
        else
           t=t+1
        endif
        if ( t>green_nt ) then 
           write(6,*) 'Travel time of  ', station_code, 'could not be read'
           !check green error
           write(6, *) maxval(tempGij(1:green_nt, x, y)), x, y, 'check'          
           write(88, *) x, y, station_code
           exit
        endif
     enddo
  enddo; enddo
  
  ampGij=0
  !Below is used for calculating normalizing factor for Potency BP g_{ij}^{p}
  !This routine fails to search local min. and max. if peak and trough of GF are located next to each other.
  !Written on 180918171517 by Ryo Okuwaki
  !SearchWindow = 0.1 / dtg !0.1 [s] / 0.05 [sampling sec] 
  !do x = 1, nx
     !do y = 1, ny
        !minT=Tij(x,y) - SearchWindow
        !maxT=Tij(x,y) + SearchWindow
        !peakmax=maxval(TempGij(minT:maxT,x,y))
        !peakmin=minval(TempGij(minT:maxT,x,y))
        !if (abs(peakmax) >= abs(peakmin)) then
           !ampGij(x,y) = peakmax
        !else
           !ampGij(x,y) = peakmin
        !end if
     !enddo
  !enddo

  !New routine for searching local max. or min. for normalizing factor g_{ij}^{p}
  !Written on 180918171517 by Ryo Okuwaki
  do x = 1, nx
     do y = 1, ny
        if (sign(1.0, TempGij(Tij(x, y), x, y)) > 0) then !polarity is up
           do t = 1, 11
              tmpt = Tij(x, y) - 2 + t
              if (TempGij(tmpt, x, y) .gt. TempGij(tmpt+1, x, y)) then
                 ampGij(x, y) = TempGij(tmpt, x, y)
                 exit
              end if
           end do
        else !polarity is down
           do t = 1, 11
              tmpt = Tij(x, y) - 2 + t
              if (TempGij(tmpt, x, y) .le. TempGij(tmpt+1, x, y)) then
                 ampGij(x, y) = TempGij(tmpt, x, y)
                 exit
              end if
           end do
        end if
     end do
  end do
  
  x = nint(nx/2.)
  y = nint(ny/2.)
  !NOTE!
  !polarity should be variable on each source knot? (2017-08-09)
  polarity = sign(1.0, sum ( TempGij (Tij(x,y):Tij(x,y)+10, x, y) ) )
  if (acvelflag == 0) then
     !$omp parallel do private(x, y)
     do x = 1, nx
        do y = 1, ny
           call diff_wave(TempGij(1:green_nt,x,y),dtg,green_nt)
           call taper_all(TempGij(1:green_nt,x,y),green_nt,0.05)
           call bp_filter(TempGij(1:green_nt,x,y),f1,f2,dtg,green_nt)
        end do
     end do
  else if (acvelflag == 1) then
     !$omp parallel do private(x, y)
     do x = 1, nx
        do y = 1, ny
           !call diff_wave(TempGij(1:green_nt,x,y),dtg,green_nt)
           call taper_all(TempGij(1:green_nt,x,y),green_nt,0.05)
           call bp_filter(TempGij(1:green_nt,x,y),f1,f2,dtg,green_nt)
        end do
     end do
  else
     write(6, *) 'Input either 0 or 1 for acvelflag'
     stop
  end if
  return
end subroutine read_green

subroutine shift_normalized_green(TempGij,Tij,alpha,Gij,nt,nx,ny,nw,mbpm_flag,potflag, ampGij)
  implicit none
  integer, intent(in) :: nt, nx, ny, nw, mbpm_flag, potflag
  real, intent(in):: TempGij(nt,nx,ny)
  integer, intent(in):: Tij(nx,ny)
  real, intent(in):: alpha
  real, intent(out):: Gij(nw,nx,ny)
  real, intent(inout):: ampGij(nx,ny)
  real :: dtg, w
  integer :: t, x, y
!  integer,parameter :: nshift = 0  
  integer,parameter :: nshift = 5 * 20 ! 
                  ! 1 sec shift for adjust obs shift (9 sec) 20 Hz
  integer :: green_nt
!++++++++++++++++++++ Sifting Gij(t) ++++++++++++++++++++++
  forall (t=1:nw, x=1:nx, y=1:ny)  &
       Gij(t,x,y) = TempGij(t+Tij(x,y)-nshift,x,y) 
!++++++++++++++++++++ Nomalizing Gij(t) ++++++++++++++++++++++
  if (mbpm_flag == 1) then
     if (potflag == 0) then
        do x=1,nx
           do y=1,ny
              ampGij = sqrt ( sum( Gij(1:nw,x,y)**2 ) )
              Gij(1:nw,x,y) = Gij(1:nw,x,y) / ampGij(x,y)
           end do
        end do
     else if (potflag == 1) then
        do x=1,nx
           do y=1,ny
              ampGij(x,y) = sum( Gij(1:nw,x,y)**2 )
              Gij(1:nw,x,y) = Gij(1:nw,x,y) / ampGij(x,y)
           end do
        end do
     else
        write(6, '("Error! Specify potflag either 0 (OriginalBP/HBP) or 1 (PotencyBP/HBP)")')
        stop
     end if
  else if (mbpm_flag == 0) then
     !ampGij for PotBP was already defined in read_green subroutine
     !ampGij is not used for conventional BP
     do x=1,nx
        do y=1,ny
           ampGij(x,y) = ampGij(x,y)
        end do
     end do
  endif
  
  return
end subroutine shift_normalized_green

subroutine combining_Uij(Fj,Gij,Uij,polarity,nt,ng,nx,ny,mbpm_flag,potflag,ampGij)
  use, intrinsic:: iso_fortran_env, only: ERROR_UNIT
  implicit none
  !-----
  integer :: ng,nt,nx,ny,mbpm_flag, potflag
  real, intent(in)  :: Fj (nt), Gij(ng,nx,ny), ampGij(nx,ny)
  real, intent(out) :: Uij(nt,nx,ny)
  real, intent(in) :: polarity
  !----
  real :: pol, maxval_Uij
  integer :: t,x,y,j,ti
  !-----
  
  if(mbpm_flag==1) then !HBP
     !$omp parallel do private(x,y,ti)
     do y=1,ny
        do x=1,nx
           do ti=1,nt-ng ! LOOP for calculating uij(t=ti)
              Uij(ti,x,y) = dot_product(Fj(ti:ti+ng-1),Gij(1:ng, x, y))
           enddo
        enddo
     end do
  else
     !BP
     if (potflag == 1) then !kBP
        !$omp parallel do private(x,y,t)
        do y = 1, ny
           do x = 1, nx
              do t = 1, nt
                 Uij(t,x,y) = Fj(t) / ampGij(x, y) ! Suggestion from A.Kasahara (Tue Dec  6 13:41:51 JST 2016), ampGij is amplitude of first impulse (direct P) of calculated GFs. We do not have to consider polarity of observed waveform since ampGij has polatity information, and true polarity is retained by multiplying polarities of Fj and BPSij.
                 !Uij(t, x, y) = Fj(t) / maxval(Fj(1:nt))
              end do
           end do
        end do
        ! Below is used for prohibiting overflow by stacking large number
        ! Ryo Okuwaki 2017-10-23 (below fix was suggested by Yagi-sensei)
        maxval_Uij = maxval(Uij(1:nt, 1:nx, 1:ny))
        Uij = Uij/maxval_Uij
     else if (potflag == 0) then
        !$omp parallel do private(x,y,t)
        do y = 1, ny
           do x = 1, nx
              do t = 1, nt
                 Uij(t,x,y) = Fj(t) * polarity
              end do
           end do
        end do
     else
        write(6, '("Error! Specify potflag either 0 (OriginalBP/HBP) or 1 (PotencyBP/HBP)")')
        stop
     end if
  end if
  return
end subroutine combining_Uij

!no londer used... 2018-09-24 Ryo Okuwaki
subroutine nth_root_stacking_without_power_n(Uij,Tij,weight,stack,ns,nt,nx,ny,np,x0,y0)
  implicit none
  !-----
  integer :: ns,nt,nx,ny
  real:: np
  real, intent(in):: Uij(nt,nx,ny), weight
  integer, intent(in):: Tij(nx,ny)
  real(kind=8), intent(inout):: stack(ns,nx,ny)
  real(8) :: dwk
  real :: wk, wk2
  integer :: t,tw,ix,iy,x0,y0
  Integer:: it_min, it_max, delta_it
  interface nth_root
     elemental function nth_root(x, n) result(ret)
       Real, intent(in):: x, n
       !Integer, intent(in):: n
       Real(kind=kind(x)):: ret
     end function nth_root
  end interface nth_root
  !-----
  do iy = 1, ny
     do ix = 1, nx
        delta_it = Tij(ix, iy) - Tij(x0, y0)
        it_min = max(1 + delta_it, 1)
        it_max = min(ns + delta_it, nt)
        ! todo: this is the correct code. Uncomment after test.
        ! stack(it_min:it_max, x, y) = stack(it_min:it_max, x, y) + weight*dble(nth_root(Uij(it_min:it_max, x, y), np))
        !        write(0, *) Uij(it_min:it_max, ix, iy)
        ! if(ix == nx/2 .and. iy == ny/2)then
        !    write(0, *) Uij(it_min:it_max, ix, iy)
        !    write(0, *) weight*Uij(it_min:it_max, ix, iy)
        !    write(0, *) nth_root(weight*Uij(it_min:it_max, ix, iy), np) ! todo: this is wrong`
        ! end if
        stack(it_min-delta_it:it_max-delta_it, ix, iy) = stack(it_min-delta_it:it_max-delta_it, ix, iy) + nth_root(weight*Uij(it_min:it_max, ix, iy), np) ! todo: this is wrong`
        !stack(it_min-delta_it:it_max-delta_it, ix, iy) = stack(it_min-delta_it:it_max-delta_it, ix, iy) + weight*nth_root(Uij(it_min:it_max, ix, iy), np) 
     end do
  end do
  return
end subroutine nth_root_stacking_without_power_n

elemental function nth_root(x, n) result(ret)
  Real, intent(in):: x, n
!  real, intent(in):: n
  Real(kind=kind(x)):: ret
  ret = sign(abs(x)**(1.0/n), x)
end function nth_root

subroutine NthRoot(nx, ny, Tij, x0, y0, ns, nt, weight, Uij, np, stack, dt)
  !Written by Ryo Okuwaki (rokuwaki@gmail.com) on 2017-10-23
  implicit none
  integer, intent(in)        :: nx, ny, ns, nt, x0, y0, Tij(nx, ny)
  real, intent(in)           :: weight, Uij(nt, nx, ny), np, dt
  real(kind=8), intent(inout):: stack(ns, nx, ny)
  integer                    :: x, y, it_min, it_max, s, e, delta_it, t, intt
  do x = 1, nx
     do y = 1, ny
        delta_it = Tij(x, y) - Tij(x0, y0)
        it_min = max(1 + delta_it, 1)
        it_max = min(ns + delta_it, nt)
        s = it_min - delta_it
        e = it_max - delta_it
        intt = 0
        do t = s, e
           !stack(t, x, y) = stack(t, x, y) + &
           !     weight * sign(1.0, Uij(it_min+intt, x, y)) * abs(Uij(it_min+intt, x, y)) ** (1.0/np)
           !Below function is conventional form used in Yagi+2012, Oku+2014, Oku+2016, OkuYagi2017
           !, but it does not equivalent to the formulations in paper (Yagi+2012)
           stack(t, x, y) = stack(t, x, y) + &
                sign(1.0, Uij(it_min+intt, x, y)) * abs(weight*Uij(it_min+intt, x, y)) ** (1.0/np)
           intt = intt + 1
        end do
     end do
  end do
end subroutine NthRoot

subroutine NthRootNTTfluct(nx, ny, Tij, x0, y0, ns, nt, weight, Uij, np, stack, NTTfluct, dt)
  !Written by Ryo Okuwaki (rokuwaki@gmail.com) on 2017-08-08
  implicit none
  integer, intent(in)        :: nx, ny, ns, nt, x0, y0, Tij(nx, ny), NTTfluct
  real, intent(in)           :: weight, Uij(nt, nx, ny), np, dt
  real(kind=8), intent(inout):: stack(ns, nx, ny, NTTfluct)
  integer                    :: x, y, it_min, it_max, s, e, delta_it, t, intt, TTfluct, f
  do x = 1, nx
     do y = 1, ny
        do f = 1, NTTfluct
           TTfluct = (f - int(NTTfluct/2) - 1)*nint(1.0/dt)
           !if NTTfluct = 5, then, fluctuations will be +-2s
           delta_it = Tij(x, y) - Tij(x0, y0) + TTfluct
           it_min = max(1 + delta_it, 1)
           it_max = min(ns + delta_it, nt)
           s = it_min - delta_it
           e = it_max - delta_it
           intt = 0
           do t = s, e
              !stack(t, x, y, f) = stack(t, x, y, f) + &
              !     weight * sign(1.0, Uij(it_min+intt, x, y)) * abs(Uij(it_min+intt, x, y)) ** (1.0/np)
              !Below function is conventional form used in Yagi+2012, Oku+2014, Oku+2016, OkuYagi2017
              !, does not equivalent to the formulations in paper (Yagi+2012)
              stack(t, x, y, f) = stack(t, x, y, f) + &
                   sign(1.0, Uij(it_min+intt, x, y)) * abs(weight*Uij(it_min+intt, x, y)) ** (1.0/np)
              
!              stack(t, x, y, f) = stack(t, x, y, f) + &
!                   sign(1.0, Uij(it_min+intt, x, y)) * abs(weight * Uij(it_min+intt, x, y))

              
              intt = intt + 1
           end do
        end do
     end do
  end do
end subroutine NthRootNTTfluct

subroutine NthPower(stack, ns, nx, ny, np, nj)
  implicit none
  integer,intent(in)        :: ns, nx, ny, nj
  real, intent(in)          :: np
  real(kind=8),intent(inout):: stack(ns, nx, ny)
  integer                   :: t, x, y
  !$omp parallel do private(t, x, y)
  do t = 1, ns
     do x = 1, nx
        do y = 1, ny
           stack(t, x, y) = sign(1.0, stack(t, x, y))*abs(stack(t, x, y))**np
        end do
     end do
  end do
  return
end subroutine NthPower

subroutine SortMedian(x, n, median)
  implicit none
  real(kind=8), intent(in) :: x(n)
  integer, intent(in)      :: n
  real(kind=8), intent(out):: median
  
  integer            :: i1, i2, k, m
  real               :: wr
  real, allocatable  :: s(:)
  allocate(s(1:n))
  
  do k = 1, n
     s(k) = x(k)
  end do
  
  do k = 1, n - 1
     do m = k + 1, n
        if (s(k) > s(m)) then
           wr = s(k)
           s(k) = s(m)
           s(m) = wr
        end if
     end do
  end do
  
  if (mod(n, 2) == 0) then
     i1 = nint(0.5 * n); i2 = nint((0.5 * n) + 1)
     median = 0.5 * (s(i1) + s(i2))
  else
     i1 = nint(0.5 * (n + 1))
     median = s(i1)
  end if
  
  return
end subroutine SortMedian

subroutine OutputSol(stack,ns,nx,ny)
  implicit none
  integer,intent(in)     :: ns, nx, ny
  real(kind=8),intent(in):: stack(ns,nx,ny)
  !----  Working Space
  integer            ::  t, x, y
  real               ::  xp, yp, dx, dy, x1, x2, y1, y2
  real               ::  slip_xy, wk
  !----
  open(88,file='OutputSol.dat',form='unformatted', status='replace', action='write')
  write(88) nx, ny, ns
  do x = 1, nx
     do y = 1, ny
        write(88)  x,y
        write(88) (stack(t,x,y), t = 1, ns)
     enddo
  enddo
  write(6, *) "Max stacked value =", maxval(stack)
  close(88)
  return
end subroutine OutputSol
