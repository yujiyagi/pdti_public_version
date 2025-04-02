! Y. Yagi
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
function l_id(m,n,jt,icm,mn,nn,jtn)
  implicit none
  integer :: l_id
  integer,intent(in) :: m,n,jt,icm,mn,nn,jtn
  l_id=(mn*nn*jtn)*(icm-1)+mn*nn*(jt-1)+mn*(n-1)+m
  RETURN
END function l_id
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
function k_id(ik,jk,jn,ndj)
  implicit none
  integer :: k_id
  integer,intent(in) :: ik,jk,jn
  integer,intent(in) :: ndj(jn)
  integer :: kwork1,j
  kwork1=0
  do j=1,jk-1
     kwork1=kwork1+ndj(j)
  end do
  k_id=kwork1+ik
  return
END function k_id
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
function k_id_c(ik,ndata)
  implicit none
  integer :: k_id_c
  integer,intent(in) :: ik,ndata
  k_id_c=ndata+ik
  return
END function k_id_c
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
subroutine get_l_id(mn,nn,jtn,jtn0,icmn,l_id_m)
  implicit none
  integer,intent(in) :: mn,nn,jtn0,icmn
  integer,intent(in) :: jtn(mn,nn)
  integer,intent(out):: l_id_m(mn,nn,jtn0,icmn)
  integer :: l0,icm,m,n,jt
  l0 = 0
  l_id_m = l0
  do icm = 1,icmn
    do m = 1,mn;do n= 1,nn
      do jt = 1,jtn(m,n)
        l0 = l0 + 1
        l_id_m(m,n,jt,icm) = l0
      enddo
    enddo;enddo
  enddo
  !---
  return
end subroutine get_l_id
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
subroutine getJTN(jtn,tr,mn,nn,jtn_max,rtime,st_max)
  implicit none
  integer,intent(in)  :: mn,nn,jtn_max
  real,intent(in)     :: tr(mn,nn)
  real,intent(in)     :: rtime,st_max
  integer,intent(out) :: jtn(mn,nn)
  integer :: nflag(mn,nn)
  integer :: m,n,jt,ios
  do m=1,mn
    do n=1,nn
      do jt=1,jtn_max
        jtn(m,n) = jt
        if(st_max <= tr(m,n)) jtn(m,n) = 0.
        if(st_max <= jt*rtime+tr(m,n)) exit
      enddo
    enddo
  enddo
!-----
!  If a configuration file (zeroone.txt) is available, non-rectangular faults mode.
!-----
 open(4, file = 'zeroone.txt', status = 'old', iostat = ios)
 if (ios == 0) then
    do n = nn, 1, -1
       read(4,*, iostat=ios)  (nflag(m, n), m = 1, mn)
       if (ios < 0) goto 99
    end do
    write(6,*) "*** zeroone.txt was loaded properly. -> Non-Rectangular Model Plane ***"
    do m=1,mn ; do n=1,nn
      if(nflag(m, n) == 0 ) jtn(m, n) = 0
    enddo; enddo
  endif
99 close(4)
  return
end subroutine getJTN
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

