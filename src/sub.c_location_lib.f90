!Subroutines to calculate the location of the space knot point (Y. Yagi)
!------------------------------------------------                       
function getdeglat2 (alat, alon)   
  !-
  !-  partial(latitude)/partial(km)
  !-
  implicit none
  real :: getdeglat2
  real,parameter ::rad = 0.017453292, delta = 0.1
  real,intent(in) :: alat, alon
  real :: a,b,e,rc,dist
  !------ GRS80
  a =6378.137
  b = 6356.752
  e = sqrt(1. - (b/a)**2)
  rc = a*(1.-e**2)/sqrt( (1.-(e*sin( (alat+delta/2.)*rad))**2)**3 )
  dist = rc * delta * rad
  getdeglat2 = delta/dist 
  return 
end function getdeglat2
!------------------------------------------------                       
function getdeglon2 (alat, alon) 
  !-
  !-  partial(longitude)/partial(km)
  !-
  implicit none
  real :: getdeglon2
  real,parameter ::rad = 0.017453292, delta = 0.1
  real,intent(in) :: alat, alon
  real :: a,b,e,rc,dist
  !------ GRS80
  a =6378.137
  b =6356.752
  e = sqrt(1. - (b/a)**2)
  rc = (a/sqrt( 1.- (e* sin(alat*rad))**2)) * cos(alat*rad)
  dist = rc * delta * rad
  getdeglon2 = delta / dist
  return 
end function getdeglon2
!---------------------------------------------                          
subroutine xy2geo(olat,olon,xd,yd,str,dip,wlat,wlon)
!	IF You Put Strike, Dip, Epicenter and center of 
!	Subfaulte location, To output Center of Subfaulte
!	Longitude and Latitude.
  implicit none
  real,parameter :: rad = 0.017453292
  real,intent(in)  :: olat,olon,xd,yd,str,dip
  real,intent(out) :: wlat,wlon
  real :: str1,dip1,work1,work2
  real :: getdeglat2, getdeglon2
  real :: deglat,deglon
  deglat = getdeglat2 (olat, olon) 
  deglon = getdeglon2 (olat, olon) 
  dip1=dip*rad
  str1=(str-90.)*rad
  work1 =  xd*cos(str1) + yd*sin(str1)*cos(dip1)
  work2 = -xd*sin(str1) + yd*cos(str1)*cos(dip1)
  wlon = olon + deglon * work1
  wlat = olat + deglat * work2
  if(wlon > 180. ) wlon = wlon - 360.
  return
end subroutine xy2geo
