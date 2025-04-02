! Y. Yagi
!---
function disp_slip(D1,D2,slip,rslip)
  implicit none
  real disp_slip
  real,parameter :: rad =  0.017453292
  real,intent(in)  :: D1,D2,slip,rslip
  real :: x,y
  real :: disp_str,disp_dip
  x = disp_str(D1,D2,slip,rslip)
  y = disp_dip(D1,D2,slip,rslip)
  if(x.gt.0. .and. y.gt.0.)then
     disp_slip = atan(y/x)/rad
     return
  endif
  if(x.lt.0. .and. y.gt.0.)then
     disp_slip = 180. - atan(abs(y/x))/rad
     return
  endif
  if(x.lt.0. .and. y.lt.0.)then
     disp_slip =  atan(abs(y/x))/rad -180.
     return
  endif
  if(x.gt.0. .and. y.lt.0.)then
     disp_slip =  - atan(abs(y/x))/rad
     return
  endif
  if(x.eq.0. .and. y.gt.0.)then
     disp_slip =   90.
     return
  endif
  if(x.eq.0. .and. y.lt.0.)then
     disp_slip =   -90.
     return
  endif
  if(x.gt.0. .and. y.eq.0.)then
     disp_slip =   0.
     return
  endif
  if(x.lt.0. .and. y.eq.0.)then
     disp_slip =   180.
     return
  endif
  if(x.eq.0. .and. y.eq.0.)then
     disp_slip =  0.
     return
  endif
end function disp_slip
!--
function disp_length(D1,D2,slip,rslip)
  implicit none
  real :: disp_length
  real,intent(in) :: D1,D2,slip,rslip
  real :: disp_str,disp_dip
  if(d1==0. .and.  d2==0.) then
    disp_length = 0.
    return
  endif
  if(d1/=0. .and.  d2== 0.) then
    disp_length = abs(d1)
    return
  endif
  if(d1==0. .and.  d2/= 0.) then
    disp_length = abs(d2)
    return
  endif
  if(rslip == 45.) then
    disp_length = sqrt(D1**2 + D2 **2)
    return
  endif  
  disp_length = sqrt(disp_str(D1,D2,slip,rslip)**2+disp_dip(D1,D2,slip,rslip)**2)
  return
end function disp_length
!---
function disp_str(D1,D2,slip,rslip)
  implicit none
  real :: disp_str
  real,parameter :: rad =  0.017453292
  real,intent(in)  ::  D1,D2,slip,rslip
  disp_str = D1*cos((slip-rslip)*rad) + D2*cos((slip+rslip)*rad)
  return
end function disp_str
!---
function disp_dip(D1,D2,slip,rslip)
  implicit none
  real :: disp_dip
  real,parameter :: rad =  0.017453292
  real,intent(in)  ::  D1,D2,slip,rslip
  disp_dip =  D1*sin((slip-rslip)*rad) + D2*sin((slip+rslip)*rad)
  return
end function disp_dip
