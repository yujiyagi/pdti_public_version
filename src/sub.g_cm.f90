!  Subroutines to create the data covariance matrix (Y. Yagi)
!
subroutine get_Cm(alpha2,obs_m,grn_m,ndmax,ndj,jn,Cm,kmax,det_cm)
  !----------------------------------------------------
  implicit none
  real,intent(in)   :: alpha2
  integer,intent(in):: ndmax,jn,kmax
  integer,intent(in):: ndj(jn)
  real,intent(in):: obs_m(ndmax,ndmax,jn)
  real,intent(in):: grn_m(ndmax,ndmax,jn)
  real,intent(out) :: Cm(kmax,kmax) ! Invers of Data Covariance 
  real,intent(out) :: det_cm  ! Determinant of Data Covariance Matrix
  ! working space
  real,allocatable :: Wm(:,:)
  real:: wk
  integer :: j,i1,i2,k1,k2
  integer :: k_id
  !----------------------------------------------------
  Cm = 0.
  det_Cm = 0.
  do j=1,jn
    allocate(Wm(1:ndj(j),1:ndj(j)))
    do i1=1,ndj(j)
      do i2 = 1,ndj(j)
        Wm(i1,i2) = obs_m(i1,i2,j)*alpha2  +  grn_m(i1,i2,j)  ! Unnecessary "dble" removed.(2025/03/13)
      enddo
    enddo
    call get_inv_LU(ndj(j),Wm,ndj(j),wk)
    det_cm = wk + det_cm
    do i1=1,ndj(j)     ;   k1=k_id(i1,j,jn,ndj)
      do i2 = 1,ndj(j) ;   k2=k_id(i2,j,jn,ndj)
        Cm(k1,k2) = Wm(i1,i2)
      enddo
    enddo
    deallocate(Wm)
  enddo    !  J (station) loop
  return
end subroutine get_cm
!------------------------------
!------------------------------
!------------------------------
subroutine get_dCm(alpha2,obs_m,grn_m,ndmax,ndj,jn,dcm,kmax,det_cm)
  !----------------------------------------------------
  implicit none
  real,intent(in)   :: alpha2
  integer,intent(in):: ndmax,jn,kmax
  integer,intent(in):: ndj(jn)
  real,intent(in):: obs_m(ndmax,ndmax,jn)
  real,intent(in):: grn_m(ndmax,ndmax,jn)
  real(8),intent(out) :: dcm(kmax,kmax) ! Invers of Data Covariance 
  real,intent(out) :: det_cm  ! Determinant of Data Covariance Matrix
  ! working space
  real(8),allocatable :: dwwm(:,:)
  real(8):: dwk
  integer :: j,i1,i2,k1,k2
  integer :: k_id
  !----------------------------------------------------
  dCm = 0.
  det_Cm = 0.
  do j=1,jn
    allocate(dwwm(1:ndj(j),1:ndj(j)))
    do i1=1,ndj(j)
      do i2 = 1,ndj(j)
        dwwm(i1,i2) = dble(obs_m(i1,i2,j)*alpha2  +  grn_m(i1,i2,j) )
      enddo
    enddo
    call get_inv_LU_d(ndj(j),dwwm,ndj(j),dwk)
    det_cm = real(dwk) + det_cm
     do i1=1,ndj(j)     ;   k1=k_id(i1,j,jn,ndj)
       do i2 = 1,ndj(j) ;   k2=k_id(i2,j,jn,ndj)
         dCm(k1,k2) = dwwm(i1,i2)
       enddo
     enddo
     deallocate(dwwm)
  enddo    !  J (station) loop
  return
end subroutine get_dcm
!------------------------------
!------------------------------
!------------------------------
subroutine correction_grn_m(ndata,b,grn_m,ndmax,jn,ndj,wk)
  implicit none
  integer,intent(in) :: ndata,ndmax,jn
  integer,intent(in) :: ndj(jn)
  real,intent(in)    :: b(ndata)
  real,intent(inout) :: grn_m(ndmax,ndmax,jn)
  !---
  integer :: i,j
  real    :: traceb,traceg,wk
  traceb = 0.
  do i=1,ndata
    traceb = traceb + b(i)**2
  enddo
  traceg = 0.
  do j=1,jn
    do i=1,ndj(j)
      traceg = traceg + grn_m(i,i,j)
    enddo
  enddo
  wk =  traceb / traceg
  grn_m = grn_m * wk
  return
end subroutine correction_grn_m
!------------------------------
!------------------------------
!------------------------------
subroutine correction_obs_m(ndata,b,obs_m,ndmax,jn,ndj,sigw,sigwd,wk)
  implicit none
  integer,intent(in) :: ndata,ndmax,jn
  integer,intent(in) :: ndj(jn)
  real,intent(in)    :: b(ndata),sigw(jn),sigwd(jn)
  real,intent(inout) :: obs_m(ndmax,ndmax,jn)
  integer :: j,i1,i2
  real    :: wk,wk1
  do j=1,jn
    !  sigwd : variance before P-wave
    wk1 = sigwd(j) ** 2
    do i1=1,ndj(j)    
      do i2 = 1,ndj(j) 
        obs_m(i1,i2,j) = obs_m(i1,i2,j)*wk1
      enddo
    enddo
  enddo    !  J (station) loop
  call correction_grn_m(ndata,b,obs_m,ndmax,jn,ndj,wk)
  return
end subroutine correction_obs_m

