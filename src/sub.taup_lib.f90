!---
! The AK135 travel-time and ray-parameter tables calculated using TauP are used.
!
! This program was originally written by R. Yamaguchi. 
! Y. Yagi redesigned and rewrote it to improve accuracy and support parallelisation.
!  (2024.6.30)
!---
! Usage 
!- 
! real :: rayp(110,19)
! call load_taup_file(rayp, "rayP" )
! p  = get_tauP(delta, depth, rayp, "rayP") ! ray parameter at source point
! g  = get_taupP(delta, depth, rayp, "G")   ! geometrical spreading factor
!-
! real :: TPtable(110,19)
! call load_taup_file(TPtable, "tP" )
! tp = get_tauP(delta, depth, TPtable,"tP" ) ! P-wave arrival time 
!-
! real :: TPPtable(110,19)
! call load_taup_file(TPPtable, "tPP" )
! tp = get_tauP(delta, depth, TPPtable, "tPP") ! PP-wave arrival time 
!--- 
subroutine load_taup_file(para, cflag)
  implicit none
  integer,parameter   :: n_del=110, n_dep=19
  real,intent(out)    :: para(n_del,n_dep)
  character(*),intent(in)  :: cflag
  !-------
  integer:: i_del, i_dep, i
  real,parameter ::  rad = 0.017453292
  character :: homedir*20
  !-------
  CALL GETENV( 'HOME', homedir )
  i = index(homedir,' ') -1 
  if(cflag == 'rayP' ) open(20, file=homedir(1:i)//'/data/ak135_taup.rayp', status='old', action='read')
  if(cflag == 'tP'   ) open(20, file=homedir(1:i)//'/data/ak135_taup.ptime', status='old', action='read')
  if(cflag == 'tPP'  ) open(20, file=homedir(1:i)//'/data/ak135_taup.pptime', status='old', action='read')
  para = 0.
  do i = 1, n_del
    read(20,*, end=100) i_del, (para (i_del,i_dep), i_dep = 1, n_dep)
  enddo
100 close(20)
  if(cflag == 'rayP' )  para = para / rad    ! [input s/deg]  -> [rayp = r*sin(i)/v :r(km),vp(km/s)]
  if(maxval(para) == 0. ) write(6,*) "(ERROR) subroutine load_taup_file : cflag = " ,cflag 
  return
end subroutine load_taup_file
!----------
!----------
function get_tauP(delta, depth, para, cflag)
  implicit none
  real :: get_tauP
  real,intent(in)     :: delta, depth
  integer,parameter   :: n_del=110
  real,intent(in)    :: para(n_del,*)
  character(*),intent(in)  :: cflag
  real :: get_G_taup, get_para_taup
  real,parameter :: R0 = 6371.
  get_tauP = 99999.
  if(cflag == 'rayP' ) get_tauP = get_para_taup(delta, depth, para)/(R0-depth)
  if(cflag == 'G'    ) get_tauP = get_G_taup(delta, depth, para) 
  if(cflag == 'tP'   ) get_tauP = get_para_taup(delta, depth, para)
  if(cflag == 'tPP'  ) then 
      get_tauP = get_para_taup(delta, depth, para) 
      if(get_tauP .gt. 1150. ) get_tauP = 9999.
  endif
  if(get_tauP == 99999. ) write(6,*) "(ERROR) function get_tauP : cflag = " ,cflag ,delta,depth
  return
end function
!----------
!----------
function get_G_taup(delta, depth, rayp)
  implicit none
  real :: get_G_taup
  real,intent(in)     :: delta, depth
  integer,parameter   :: n_del=110
  real,intent(in)     :: rayp(n_del,*)
  !-------
  real :: dtih, vp, p_0, p_1, p_2, R1
  real :: get_para_taup, get_velp, get_p_angle
  real,parameter :: R0 = 6371., rad = 0.017453292, vp0 = 5.8
  real :: dd 
  real :: d_min, d_max, G_min, G_max, get_G_taup_interpolation
  !-------
  vp = get_velp (depth)
  R1 = R0 - depth
  !---
  dd = 1.0 
  !--- Mitigate abrupt drop in the geometrical spreading factor (2025/2/21 Yagi)
  d_min = 85.
  d_max = 100.
  if(delta >= d_min .and. delta < d_max ) then 
    p_1 = get_para_taup(d_min-dd,depth,rayp)
    p_2 = get_para_taup(d_min+dd,depth,rayp)
    dtih = abs(asin(p_1/R1*vp)/rad - asin(p_2/R1*vp)/rad) /(2.*dd)
    p_1 = get_para_taup(d_min,depth,rayp)/R1 * vp   ! sin(takeoff)
    p_2 = get_para_taup(d_min,depth,rayp)/R0 * vp0  ! sin(inci)
    G_min = sqrt( p_1 * dtih / (1.-p_2**2) /sin(d_min*rad) )  / r0 * 1.e5
    !---
    p_1 = get_para_taup(d_max-dd,depth,rayp)
    p_2 = get_para_taup(d_max+dd,depth,rayp)
    dtih = abs(asin(p_1/R1*vp)/rad - asin(p_2/R1*vp)/rad) /(2.*dd)
    p_1 = get_para_taup(d_max,depth,rayp)/R1 * vp   ! sin(takeoff)
    p_2 = get_para_taup(d_max,depth,rayp)/R0 * vp0  ! sin(inci)
    G_max = sqrt( p_1 * dtih / (1.-p_2**2) /sin(d_max*rad) )  / r0 * 1.e5
    !---
    get_G_taup_interpolation = G_min +  (delta-d_min) * (G_max-G_min)/(d_max-d_min) 
    !-------------------------------------
    p_1 = get_para_taup(delta-dd,depth,rayp)
    p_2 = get_para_taup(delta+dd,depth,rayp)
    dtih = abs(asin(p_1/R1*vp)/rad - asin(p_2/R1*vp)/rad) /(2.*dd)
    p_1 = get_para_taup(delta,depth,rayp)/R1 * vp   ! sin(takeoff)
    p_2 = get_para_taup(delta,depth,rayp)/R0 * vp0  ! sin(inci)
    get_G_taup = sqrt( p_1 * dtih / (1.-p_2**2) /sin(delta*rad) )  / r0 * 1.e5
    get_G_taup = max(get_G_taup_interpolation,get_G_taup)
  else
    p_1 = get_para_taup(delta-dd,depth,rayp) 
    p_2 = get_para_taup(delta+dd,depth,rayp) 
    dtih = abs(asin(p_1/R1*vp)/rad - asin(p_2/R1*vp)/rad) /(2.*dd)
    !---
    p_1 = get_para_taup(delta,depth,rayp)/R1 * vp   ! sin(takeoff)
    p_2 = get_para_taup(delta,depth,rayp)/R0 * vp0  ! sin(inci)
    get_G_taup = sqrt( p_1 * dtih / (1.-p_2**2) /sin(delta*rad) )  / r0 * 1.e5
  endif
  ! Velocity and density contrasts are not considered here because they are corrected in the Haskell matrix
  ! by A. Kasahara (2016).
  return
end function get_G_taup
!----------
!----------
function get_velp ( depth )
  implicit none
  real :: get_velp
  real,intent(in) :: depth
  integer :: n
  !---
  integer,parameter :: nl = 25
  real :: dl(nl),vp(nl)
  !---  AK 135 model
  data dl /0.,20.,20.,35.,35.,77.5,120.,120.,165.,210.,210.,260.,310.,360.,410.,410., &
          460.,510.,560.,610.,660.,660.,710.,760.,809.5 /
  data vp /5.8000,5.8000,6.5000,6.5000,8.0400,8.0450,8.0500,8.0505,8.1750,8.3007,8.3007,&
           8.4822,8.6650,8.8476,9.0302,9.3601,9.5280,9.6962,9.8640,10.0320,10.2000,&
          10.7909,10.9222,11.0553,11.1355 /
  !---
  do n = 2, nl
    if( dl(n)  > depth ) then
       get_velp = vp(n-1) + (depth-dl(n-1))*(vp(n)-vp(n-1))/(dl(n)-dl(n-1))
       exit
    endif
  enddo
  return
end function get_velp
!----------
!----------
!----------
function get_para_taup(delta, depth, para)
  implicit none
  real :: get_para_taup
  real,intent(in) :: delta, depth
  integer,parameter  :: n_del=110, n_dep=19
  real,intent(in) :: para(n_del, n_dep)
  !-------
  integer :: i,idelta,idepth
  real :: p1, p2
  real :: h0(n_dep)
  data h0 /  0., 20., 35., 78., 120., 165., 210., 260., 310., 360., 410., 460., 510., 560., 610., 660., 710., 760., 809. /
  !-------
  idelta = delta
  do i=2,n_dep
    if(depth < h0(i) ) then
      idepth = i-1
      exit
    endif
  enddo
  p1 = para(idelta,idepth)   +  (delta - idelta) * (para(idelta+1,idepth)   - para(idelta,idepth)  )
  p2 = para(idelta,idepth+1) +  (delta - idelta) * (para(idelta+1,idepth+1) - para(idelta,idepth+1))
  get_para_taup = p1 + (depth - h0(idepth)) * (p2 - p1) / (h0(idepth+1) - h0(idepth))
  return
end function get_para_taup
!----------
!----------
