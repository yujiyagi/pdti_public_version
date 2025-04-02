!  The sourutines were developed based on sourutines written by 
!  Prof. Kikuchi and Prof. Kanamori. Prof. Kanamori has kindly agreed 
!  to release to the public a modified version. 
!  Conversion to Fortran90 and minor modifications avoid unexpected 
!    problems during parallelization and optimisation.  (Yuji Yagi)
!---
!  The original programs are available at "Note on Teleseismic Body-Wave 
!     Inversion Program (https://wwweic.eri.u-tokyo.ac.jp/ETAL/KIKUCHI/).
!---------------------------------                                      
function l_source(h0, dep, nl)
  implicit none
  integer :: l_source
  real,intent(in)    :: h0
  real,intent(inout) :: dep(nl)
  integer,intent(in) :: nl
  integer :: l
  real :: hl,dh
  !-------------
  hl = 0. 
  if(dep(nl) /= 0. ) dep (nl) = 0. 
  do l = 1, nl - 1 
     hl = hl + dep (l) 
     dh = h0 - hl 
     if (dh.lt.0.) exit  
  enddo
  l_source = l 
end function l_source
!---------------------------------                                      
function getdelta (alat1, alon1, alat2, alon2) 
  implicit none
  real :: getdelta
  real, intent(in) ::  alat1, alon1, alat2, alon2 
  real :: rad = 0.017453292
  real :: t0,t1,dphi,t00,t11,c00,c11,s00,s11,c2,d11
  !-------------
  t0 = alat1 * rad 
  t1 = alat2 * rad 
  dphi = (alon2 - alon1) * rad 
  t00 = t0 - 11.55 / 60 * rad * sin (2 * t0) 
  t11 = t1 - 11.55 / 60 * rad * sin (2 * t1) 
  c00 = cos (t00) 
  c11 = cos (t11) 
  s00 = sin (t00) 
  s11 = sin (t11) 
  c2 = cos (dphi) 
  d11 = c00 * c11 * c2 + s00 * s11 
  getdelta = acos (d11) / rad 
  return 
end function getdelta                         
!---------------------------------                                      
function getaz (alat1, alon1, alat2, alon2) 
  implicit none
  real :: getaz
  real,intent(in)::  alat1, alon1, alat2, alon2 
  real :: rad = 0.017453292
  real :: t0,t1,dphi,c0,c1,s0,s1,c2,s2,a1,a2
  t0 = alat1 * rad 
  t1 = alat2 * rad 
  dphi = (alon2 - alon1) * rad 
  c0 = cos (t0) 
  c1 = cos (t1) 
  s0 = sin (t0) 
  s1 = sin (t1) 
  c2 = cos (dphi) 
  s2 = sin (dphi) 
  a1 = c1 * s2 
  a2 = s1 * c0 - s0 * c1 * c2 
  getaz = atan2 (a1, a2) / rad 
  return 
end function getaz                            
!---------------------------------                                      
