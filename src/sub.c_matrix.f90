! Subroutines for matrix operations  (Y. Yagi)
! It should not be used because of its slow speed.
!-------------------------
function get_sig(km,lm,a,kmax,b,x)
  implicit none
  real get_sig
  integer,intent(in) :: km,lm,kmax
  real,intent(in) :: a(kmax,*),b(*),x(*)
  real(8)         ::zz,wk1
  integer :: k,l
  !---
  wk1 = 0.d0
  do k=1,km
    zz = 0.d0
    do l=1,lm
       zz = zz + dprod(a(k,l),x(l))
    enddo
    wk1 = wk1 + (dble(b(k)) - zz)**2
  enddo
  get_sig = wk1
  return
end function get_sig
!-------------------------
!-------------------------
subroutine multi_atb(km,lm,a,kmax,b,atb)
  implicit none
  integer,intent(in) :: km,lm,kmax
  real,intent(in)  :: a(kmax,*),b(*)
  real,intent(out) :: atb(*)
  real(8):: wk1
  integer :: k,l
  !$omp parallel do shared(a,b,atb) private(k,wk1)
  do l =  1, lm
    wk1 = 0.0d0
    do k = 1,km
       wk1 = wk1 +  dprod(a(k,l),b(k))
    end do
    atb(l) =  wk1
  enddo
  return
END subroutine multi_atb
!-------------------------
!-------------------------
subroutine multi_ab(km,lm,a,kmax,b,ab)
  !  ab(k) = a(k,l) * b(l)
  implicit none
  integer,intent(in) :: km,lm,kmax
  real,intent(in)  :: a(kmax,lm),b(lm)
  real,intent(out) :: ab(km)
  real(8):: wk1
  integer :: k,l
  !$omp parallel do shared(a,b,ab) private(l,wk1)
  do k = 1,km
    wk1 = 0.0d0
    do l = 1,lm
      wk1 = wk1 +  dprod(a(k,l),b(l))
    end do
    ab(k) = wk1
  enddo
  return
END subroutine multi_ab
!-------------------------
subroutine symetric_m(lm,a,lmax)
  implicit none
  integer,intent(in)    :: lm, lmax
  real,intent(inout) :: a(lmax,*)
  integer :: l1,l2
  do l1 = 2,lm
     do l2 = 1,l1-1
        a(l1,l2) = a(l2,l1)
     end do
  end do
  return
end subroutine symetric_m
!-------------------------
!-------------------------
subroutine symetric_m_d(lm,a,lmax)
  implicit none
  integer,intent(in)    :: lm, lmax
  real(8),intent(inout) :: a(lmax,*)
  integer :: l1,l2
  do l1 = 2,lm
     do l2 = 1,l1-1
        a(l1,l2) = a(l2,l1)
     end do
  end do
  return
end subroutine symetric_m_d
!-------------------------
!-------------------------
subroutine multi_ata(km,lm,a,kmax,ata,lmax)
  !     ata(l1,l2) =  a(k,l1)' * a(k,l1)
  implicit none
  integer,intent(in)    :: km,lm,kmax,lmax
  real,intent(in)       :: a(kmax,*)
  real,intent(out)      :: ata(lmax,*)
  real(8):: wk1
  integer:: l1,l2,k
  do l1 = 1, lm 
    !$omp parallel do shared(a,ata,l1) private(k,wk1)
    do l2 = l1, lm
      wk1=0.d0
      do k = 1,km
        wk1 = wk1 +  dprod(a(k,l1),a(k,l2))
      enddo
      ata(l1,l2) = wk1
    end do
  end do
  !======================================
  do l1 = 2,lm
     do l2 = 1,l1-1
        ata(l1,l2) = ata(l2,l1)
     end do
  end do
  return
END subroutine multi_ata
!-------------------------
!-------------------------
subroutine multi_ata_u(km,lm,a,kmax,ata,lmax)
  !     ata(l1,l2) =  a(k,l1)' * a(k,l1)
  implicit none
  integer,intent(in)  :: km,lm,kmax,lmax
  real,intent(in)    ::  a(kmax,*)
  real,intent(out)   ::  ata(lmax,*)
  real(8):: wk1
  integer :: l1,l2,k
  do l1 = 1, lm
    !$omp parallel do shared(a,ata,l1) private(k,wk1)
    do l2 = l1, lm
      wk1=0.d0
      do k = 1,km
        wk1 = wk1 +  dprod(a(k,l1),a(k,l2))
      enddo
      ata(l1,l2) = wk1
    end do
  end do
  return
END subroutine multi_ata_u
!-------------------------
!-------------------------
subroutine multi_m(km,lm,nm,r1,im1,r2,im2,r3,im3)
  !  r3(k,n) = r1(k,l) * r2(l,n)
  !--
  implicit none
  integer,intent(in) :: km,lm,nm,im1,im2,im3
  real,intent(in)    :: r1(im1,*),r2(im2,*)
  real,intent(out)   :: r3(im3,*)
  real(8):: wk1
  integer :: k,n,l
  !$omp parallel do shared(r1,r2,r3) private(n,l,wk1)
  do k =  1, km
    do n =  1, nm
      wk1 = 0.0d0
      do l = 1,lm
        wk1 = wk1 +  dprod(r1(k,l),r2(l,n))
      end do
      r3(k,n) = wk1
    enddo
  enddo
  return
END subroutine multi_m
!-------------------------
!-------------------------
function f_in_product(km,x,y)
  implicit none
  real:: f_in_product
  integer,intent(in) :: km
  real,intent(in)    :: x(km),y(km)
  real(8):: wk1
  integer :: k
  wk1 = 0.0d0
  do k = 1,km
    wk1 = wk1 + dprod(x(k),y(k))
  enddo
  f_in_product = wk1
  return
end function f_in_product
!-------------------------
!-------------------------
function f_nrom2(km,x)
  implicit none
  real:: f_nrom2
  integer,intent(in) :: km
  real,intent(in)    :: x(km)
  real(8):: wk1
  integer :: k
  wk1 = 0.0d0
  do k = 1,km
    wk1 = wk1 + dprod(x(k),x(k))
  enddo
  f_nrom2 = wk1
  return
end function f_nrom2
!-------------------------
!-------------------------
function amaxval_trace(ndata,c)
  implicit none
  real :: amaxval_trace
  integer,intent(in) :: ndata
  real,intent(in)    :: c(ndata,*)
  integer :: i
  amaxval_trace  = -1.e20
  do i=1,ndata
   if(c(i,i).gt.amaxval_trace) amaxval_trace = c(i,i)
  enddo
  return
end function amaxval_trace
!-------------------------
!-------------------------
real function get_trace(ata,lmax)
  implicit none
  integer, intent(in):: lmax
  real,    intent(in):: ata(lmax,lmax)
  integer:: l 
  get_trace = 0.
  do l = 1, lmax
    get_trace = get_trace + ata(l,l)
  enddo
  return
end function get_trace
