!  The sourutines were developed based on sourutines written by 
!  Prof. Kikuchi and Prof. Kanamori. Prof. Kanamori has kindly agreed to 
!  release to the public a modified version. 
!  Changed to fortran 90 format, fixed a bug that originated from numerical errors, 
!  and changed to use lapack for eigenvalue estimation. (Yuji Yagi)
!----
!  The original programs are available at "Note on Teleseismic Body-Wave 
!     Inversion Program (https://wwweic.eri.u-tokyo.ac.jp/ETAL/KIKUCHI/).
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
subroutine d_cp(v,f1,d1,a1,sm,dsm,ax) 
  !     << to interpret a moment-tensor in terms of a double-couple. >>   
  implicit none
  real,intent(in)  :: v(6)
  real,intent(out) :: f1,d1,a1,sm,ax(3)
  !------
  real,parameter   :: pi=3.1415926, rad=0.017453292, root2=1.41421356
  integer :: ii,icon
  real    :: eg1(3),eg2(3),Ev(3,3),Mij(3,3) 
  real    :: dsm,cf,sf,c,s
  real    :: acos_y
  !     output                                                            
  !     moment      sm +- dsm                                             
  !     fault mechanism (f1,t1,a1)                                        
  !     tension axis   ax()                                               
  !     input                                                             
  !     vector v(6)                                                       
  !     mij(*) = moment tensor   (x:north, y:east, z:downward vertical)   
  !     eg(i)  = eigen values of moment-tensor (to be modified later)     
  !     ev(*,i)= eigen vector for i-th eigen-value                        
  call Mtrx(v,Mij) 
  ! call eig1(mij,3,3,2,eg1,eg2,ev,wk1,icon) 
  call eig1_y(Mij,3,eg1,Ev,icon)         ! changed to use lapack 
  if(icon.ne.0) print *,'  seig1.cond= ',icon 
  sm=(eg1(1)-eg1(3))/2                   ! M_0 for best-fitting double couple 
  dsm=sm-eg1(1)
  do ii=1,3 
     ax(ii)=Ev(ii,1) 
     eg1(ii)=(Ev(ii,3)-Ev(ii,1))/root2   !slip vector  
     eg2(ii)=(Ev(ii,1)+Ev(ii,3))/root2   !normal vector for a fault plane
  end do
  d1=acos_y(eg2(3))                      ! Dealing with numerical error 
  f1=atan2(eg2(2),eg2(1))+pi/2 
  cf=cos(f1) 
  sf=sin(f1) 
  c=eg1(1)*cf+eg1(2)*sf 
  if(d1.ne.0.) s=eg1(3)/sin(d1) 
  if(d1.eq.0.) s=eg1(1)*sf-eg1(2)*cf 
  f1=f1/rad 
  d1=d1/rad 
  a1=-atan2(s,c)/rad 
  if(f1.lt.0.)   f1=f1+360.
  if(d1.lt.90.) return 
  a1=-a1 
  d1=180.-d1 
  f1=f1+180. 
  if(f1.gt.360.) f1=f1-360.
end subroutine d_cp
!----------------------------------------------
subroutine Mtrx(v,M) 
  !     < basis tensor:v(n) to moment-tensor:M(i,j)                       
  implicit none
  real,intent(in)  :: v(6)
  real,intent(out) :: M(3,3) 
  M(1,1)=v(2)-v(5)+v(6) 
  M(1,2)=v(1) 
  M(2,1)=M(1,2) 
  M(1,3)=v(4) 
  M(3,1)=M(1,3) 
  M(2,2)=-v(2)+v(6) 
  M(2,3)=v(3) 
  M(3,2)=M(2,3) 
  M(3,3)=v(5)+v(6) 
end subroutine Mtrx
!----------------------------------------------
subroutine conj(f1,d1,a1,f2,d2,r2) 
  !******************************************                             
  !*      conjugate nodal plane            **                             
  !******************************************                             
  implicit none
  real,parameter :: pi=3.1415926, rad=0.017453292
  real,intent(in)  :: f1,d1,a1
  real,intent(out) :: f2,d2,r2
  !---
  real :: rf1,rd1,ra1,rd2,rf2,ra2 
  real :: acos_y
  !---
  if(f1.eq.0..and.d1.eq.0..and.a1.eq.0.) stop 
  rf1=rad*f1 
  rd1=rad*d1 
  ra1=rad*a1 
  if(a1.ge.0.) then 
     rd2=acos_y( sin(rd1)*sin(ra1) )              ! Dealing with numerical error
     rf2=rf1+atan2(cos(ra1),sin(ra1)*cos(rd1))+pi 
     ra2=pi-acos_y( cos(ra1)*sin(rd1)/sin(rd2) )  ! Dealing with numerical error
  else 
     rd2=pi-acos_y( sin(rd1)*sin(ra1) )           ! Dealing with numerical error
     rf2=rf1+atan2(cos(ra1),sin(ra1)*cos(rd1)) 
     ra2=-pi+acos_y( cos(ra1)*sin(rd1)/sin(rd2) ) ! Dealing with numerical error
  end if
  if(rf2.gt.2*pi) rf2=rf2-2*pi 
  f2 = rf2/rad 
  d2 = rd2/rad 
  r2 = ra2/rad 
  if(f2.lt.0.) f2 = 360. + f2 
  return 
end subroutine conj
!----------------------------------------------
subroutine Mxy_Mrf(Mxy,Mrf) 
  implicit none
  real,intent(in)  :: Mxy(3,3)
  real,intent(out) :: Mrf(3,3) 
  Mrf(1,1) =  Mxy(3,3) 
  Mrf(1,2) =  Mxy(3,1) 
  Mrf(1,3) = -Mxy(3,2) 
  Mrf(2,1) =  Mxy(1,3) 
  Mrf(2,2) =  Mxy(1,1) 
  Mrf(2,3) = -Mxy(1,2) 
  Mrf(3,1) = -Mxy(2,3) 
  Mrf(3,2) = -Mxy(2,1) 
  Mrf(3,3) =  Mxy(2,2) 
  return 
end subroutine Mxy_Mrf
!----------------------------------------------
!  Add subroutine to change to use Lapack (Y. Yagi)
subroutine eig1_y(Mij,nm,eg1,Ev,info)
  implicit none
  integer,intent(in) :: nm
  real,intent(in)    :: Mij(nm,nm)
  real,intent(out)   :: eg1(nm), Ev(nm,nm)
  integer,intent(out):: info
  !-----
  real:: wkv(nm*3-1)
  integer :: i,j
  real:: wk
  !-----
  Ev = Mij
  call SSYEV ('V','L', nm, Ev , nm, eg1, wkv, nm*3-1, info)
  do i=1,nm-1 
     do j=i+1,nm
        if(eg1(i) < eg1(j) ) then
           wk = eg1(i)
           eg1(i)=eg1(j)
           eg1(j)=wk
           wkv(1:nm)=Ev(1:nm,i)
           Ev(1:nm,i)=Ev(1:nm,j)
           Ev(1:nm,j)=wkv(1:nm)
        endif
     enddo
  enddo
  return
end subroutine eig1_y
!-----------------------------------------------
!  Add subroutines to address numerical error issues (Y. Yagi)
function acos_y(x)
  implicit none
  real :: acos_y
  real,intent(in) :: x
  if( x > 1. ) then
     acos_y = acos(1.)
     return
  endif
  if( x < -1. ) then
     acos_y = acos(-1.)
     return
  endif
  acos_y = acos(x)
  return
end function acos_y
