! Y. Yagi
!===================================================
subroutine get_sol_posv_d(n,A,nmax,b,x,det_m)
  implicit none
  integer,intent(in)               :: n, nmax
  double precision, intent(in)     :: A(nmax,n),b(n)
  double precision, intent(out)    :: x(n)
  double precision, intent(out)    :: det_m
  !----  Working Space
  double precision, allocatable    :: Wm(:,:)
  integer, parameter               :: nrhs = 1
  integer                          :: l
  integer                          :: INFO
  !---------------------
  external :: DCOPY, DPOSV
  !---------------------
  allocate(Wm(n,nrhs))
  call DCOPY( n, b, 1, Wm(1,1), 1 )
  call DPOSV ('U', n, NRHS, A, nmax, Wm, n, INFO)
  det_m = 0.d0
  do l=1,n
    det_m = det_m + dlog(dabs(A(l,l)))
  enddo
  det_m  = det_m *2.d0
  x(1:n) = Wm(1:n,1)
  deallocate(Wm)
  return
end subroutine get_sol_posv_d
!===================================================
subroutine get_sol_posv(n,A,nmax,b,x,det_m,info)
  implicit none
  integer,intent(in)   :: n, nmax
  real, intent(in)     :: A(nmax,n),b(n)
  real, intent(out)    :: x(n)
  real, intent(out)    :: det_m
  !----  Working Space
  real, allocatable    :: Wm(:,:)
  double precision     :: dwk
  integer, parameter   :: nrhs = 1
  integer              :: l
  integer              :: INFO
  !---------------------
  external :: SCOPY, SPOSV
  !---------------------
  allocate(Wm(n,nrhs))
  call SCOPY( n, b, 1, Wm(1,1), 1 )
  call SPOSV ('U', n, NRHS, A, nmax, Wm, n, INFO)
  dwk = 0.d0 
  do l=1,n
    dwk = dwk + dble(alog(abs(A(l,l))))
  enddo
  det_m  = real(dwk) * 2.
  call SCOPY( n, Wm(1,1), 1, x, 1 )
  deallocate(Wm)
  return
end subroutine get_sol_posv
!===================================================
subroutine get_sol_LU_d(n,a,nmax,b,x,det_m)
  implicit none
  !--------------------
  integer, intent(in)              :: n, nmax
  double precision, intent(in)     :: a(nmax,n),b(n)
  double precision, intent(out)    :: x(n)
  double precision, intent(out)    :: det_m
  !----  Working Space
  double precision, allocatable    :: wm(:,:), Atmp(:,:)
  integer, allocatable             :: ipiv(:)
  integer, parameter               :: nrhs = 1
  integer                          :: l
  integer                          :: INFO
  !--------------------
  external :: DGETRF, DGETRS
  !--------------------
  allocate(wm(n,nrhs),ipiv(n),Atmp(nmax,n))
  wm(1:n,1) = b(1:n)
  Atmp = A
  call DGETRF( n, n, Atmp, nmax, ipiv, INFO )
  if(info /= 0) then
    det_m = 0.d0
    return
  endif
  det_m = 0.d0
  do l=1,n
     det_m = det_m + dlog(dabs(Atmp(l,l)))
  end do
  call DGETRS( "N", n, nrhs, Atmp, nmax, ipiv, wm, n, INFO )
  if(info /= 0) then
    det_m = 0.d0
    return
  endif
  x(1:n) = wm(1:n,1)
  deallocate(wm,ipiv,Atmp)
  return
end subroutine get_sol_LU_d
!===================================================
subroutine get_sigx_LU_d(n,a,nmax,sigx)
  implicit none
  !--------------------
  integer, intent(in)              :: n, nmax
  double precision, intent(in)     :: a(nmax,n)
  double precision, intent(out)    :: sigx(n)
  !----  Working Space
  double precision, allocatable    :: Atmp(:,:)
  double precision, allocatable    :: work(:)
  integer, allocatable             :: ipiv(:)
  integer                          :: l, lwork
  integer                          :: INFO
  integer                          :: NB, ILAENV
  !--------------------
  external :: ILAENV, DGETRF, DGETRI
  !--------------------
  !lwork = 4*nmax+n*(nmax+n)
  NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
  IF( NB.LE.1 ) NB = MAX( 1, N )
  lwork = n*NB
  allocate(ipiv(n),Atmp(nmax,n),work(lwork))
  Atmp = A
  call DGETRF( n, n, Atmp, nmax, ipiv, info )
  call DGETRI( n, Atmp, nmax, ipiv, work, lwork, info )
  !-------------
  do l = 1, n
     sigx(l) = Atmp(l,l)
  enddo
  !test
  deallocate(ipiv,Atmp,work)
  return
end subroutine get_sigx_LU_d
!===================================================
function get_det_m(lm,a,lmax)
  ! get determinant of matrix A(l,l)
  implicit none
  real get_det_m
  integer,intent(in) :: lm,lmax
  real,intent(in)    :: A(lmax,lmax)
  !
  integer :: info, l
  real,allocatable:: atmp(:,:)
  integer,allocatable:: ipiv(:)
  external :: SGETRF
  allocate(atmp(lmax,lmax))
  atmp = a
  allocate(ipiv(lm))
  call SGETRF( lm, lm, atmp, lmax, ipiv, info )
  get_det_m = 0.
  do l=1,lm
     get_det_m = get_det_m + alog(abs(atmp(l,l)))
  end do
  deallocate(ipiv,atmp)
  return
end function get_det_m
!===================================================
!===================================================
subroutine get_inv_LU(n,A,nmax,det_a)
  !get inverse matrix of A
  !  a (input: matrix /output: inverse matrix)
  implicit none
  integer,intent(in) :: n,nmax
  real,intent(inout) :: A(nmax,n)
  real,intent(out)   :: det_a
  real,allocatable:: ww(:)
  integer,allocatable:: ipiv(:)
  integer :: ILAENV, nb, info, l
  double precision :: dwk
  external :: ILAENV, DGETRF, DGETRI
  NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
  allocate(ww(n*nb),ipiv(n))
  call SGETRF( n, n, A, nmax, ipiv, INFO )
  dwk = 0.d0
  do l=1,n
     dwk = dwk + dble(alog(abs(a(l,l))))
  end do
  det_a = real(dwk)
  call SGETRI( n, A, nmax, ipiv, ww, n*nb, INFO )
  deallocate(ww,ipiv)
  return
end subroutine get_inv_LU
!===================================================
!===================================================
subroutine get_inv_LU_d(n,a,nmax,det_a)
  !get inverse matrix of A
  !  a (input: matrix /output: inverse matrix)
  implicit none
  integer,intent(in) :: n,nmax
  double precision,intent(inout) :: A(nmax,n)
  double precision,intent(out)   :: det_a
  real(8),allocatable:: ww(:)
  integer,allocatable:: ipiv(:)
  integer :: ILAENV, nb, info, l
  external :: ILAENV, DGETRF, DGETRI
  NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
  allocate(ww(n*nb),ipiv(n))
  call DGETRF( n, n, a, nmax, ipiv, INFO )
  det_a = 0.d0
  do l=1,n
     det_a = det_a + dlog(dabs(a(l,l)))
  end do
  call DGETRI( n, a, nmax, ipiv, ww, n*nb, INFO )
  deallocate(ww,ipiv)
  return
end subroutine get_inv_LU_d
!===================================================
!===================================================
subroutine get_HtH_Htd_dtd_d(k, l, H, d, Cm, HtH, Htd, dtd)
  implicit none
  integer, intent(in):: k, l
  double precision, intent(in)  :: H(K,L), d(K), Cm(K,K)
  double precision, intent(out) :: HtH(L,L), Htd(L), dtd
  ! working space
  double precision, allocatable :: HtCm(:,:), Cmd(:)
  double precision :: one, zero
  double precision :: ddot    ! blas funtion
  external :: DGEMM, DGEMV, ddot
  !-------------
  allocate(HtCm(L,K), Cmd(K))
  one  = 1.0d0
  zero = 0.0d0
  !------------
  call DGEMM('T', 'N', l, k, k, one,    H, k, Cm, k, zero, HtCm, l )
  call DGEMM('N', 'N', l, l, k, one, HtCm, l,  H, k, zero,  HtH, l )
  call DGEMV('N', l, k, one, HtCm, l, d, 1, zero, Htd, 1 )
  call DGEMV('N', k, k, one,   Cm, k, d, 1, zero, Cmd, 1 )
  dtd = ddot(k, d, 1, Cmd, 1)
  return
end subroutine get_HtH_Htd_dtd_d
!===================================================
subroutine get_HtH_Htd_dtd(k, l, H, d, Cm, HtH, Htd, dtd)
  implicit none
  integer, intent(in):: k, l
  real, intent(in)  :: H(K,L), d(K), Cm(K,K)
  real, intent(out) :: HtH(L,L), Htd(L), dtd
  ! working space
  real, allocatable :: HtCm(:,:), Cmd(:)
  real :: one, zero
  real :: sdot    ! blas funtion
  external :: SGEMM, SGEMV, sdot
  !-------------
  allocate(HtCm(L,K), Cmd(K))
  one  = 1.0
  zero = 0.0
  !------------
  call SGEMM('T', 'N', l, k, k, one,    H, k, Cm, k, zero, HtCm, l )
  call SGEMM('N', 'N', l, l, k, one, HtCm, l,  H, k, zero,  HtH, l )
  call SGEMV('N', l, k, one, HtCm, l, d, 1, zero, Htd, 1 )
  call SGEMV('N', k, k, one,   Cm, k, d, 1, zero, Cmd, 1 )
  dtd = sdot(k, d, 1, Cmd, 1)
  return
end subroutine get_HtH_Htd_dtd
!==================================================================
subroutine  get_residual_d (L, AtA,  Atb, btb, x, s)
  implicit none
  integer, intent(in)            :: L
  double precision, intent (in)  :: AtA(L,L), Atb(L), btb, x(L)
  double precision, intent (out) :: s
  ! working space
  double precision, parameter    :: one=1.0d0, zero = 0.0d0
  double precision, allocatable  :: AtAx(:)
  double precision, allocatable  :: wv(:)
  double precision               :: ddot  ! BLAS function
  external :: DGEMV, ddot
  allocate(AtAx(L),wv(L))
  !-----------
  call DGEMV("N", L, L, one, AtA, L, x, 1, zero,  AtAx, 1)
  wv =  AtAx - 2*Atb
  s = ddot ( L, x, 1, wv, 1 ) + btb
  return
end subroutine get_residual_d
!==================================================================
subroutine  get_residual (L, AtA,  Atb, btb, x, s)
  implicit none
  integer, intent(in) :: L
  real, intent (in)   :: AtA(L,L), Atb(L), btb, x(L)
  real, intent (out)  :: s
  ! working space
  real, parameter     :: one=1.0, zero = 0.0
  real, allocatable   :: AtAx(:)
  real, allocatable   :: wv(:)
  real                :: sdot  ! BLAS function
  external :: SGEMV, sdot
  allocate(AtAx(L),wv(L))
  !-----------
  call SGEMV("N", L, L, one, AtA, L, x, 1, zero,  AtAx, 1)
  wv =  AtAx - 2 * Atb
  s = sdot ( L, x, 1, wv, 1 ) + btb
  return
end subroutine get_residual
!==================================================================
subroutine conv_y(x,y,z,nm)
!----------------------------
! Convolution (*)
!       z(output) = x(input) * y(input) 
!----------------------------
  implicit none
  integer,intent(in) :: nm
  real,intent(in)    :: x(nm),y(nm)
  real,intent(out)   :: z(nm)
  integer :: i, i1, i2
  external :: saxpy
  z = 0.
  do i=1,nm
    i1 = i
    i2 = nm -i + 1
!    call saxpy(i2, x(i1), y(1:i2),1, z(i1:nm),1 )
    call saxpy(i2, x(i1), y(1),1, z(i1),1 )
  enddo
  return
end subroutine conv_y
