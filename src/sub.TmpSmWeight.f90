! S. Yamashita and Y. Yagi
subroutine get_tmp_weight(itera, cr, w, mn, nn, jtn, icmn, rtime)
  implicit none
  integer, intent(in) :: itera, mn, nn, jtn, icmn
  real, intent(in)    :: cr, rtime
  real, intent(out)   :: w(mn,nn,jtn,icmn)
  character(50)       :: itera_str
  integer             :: i, m, n, jt, icm, jn, j
  real                :: cPmax
  real, allocatable   :: PDT_rate(:,:,:,:),Tr(:,:)
  real, allocatable   :: d_tq(:,:)
  real, allocatable   :: dummy(:,:),dummy2(:,:)
  w = 1.
  ! if (itera > 0) then
    !--- read fort.40
    open(50,file="input_fort.40")
    do i = 1, 9
      read(50,*)
    end do
    allocate (PDT_rate(1:mn,1:nn,1:jtn,1:icmn),Tr(1:mn,1:nn),dummy(1:mn,1:nn))
    do n=nn,1,-1 ; read(50,*)(dummy(m,n),m=1,mn) ; end do
    read(50,*)
    do n=nn,1,-1 ; read(50,*)(dummy(m,n),m=1,mn) ; end do
    read(50,*)
    do n=nn,1,-1 ; read(50,*)(Tr(m,n),m=1,mn) ; end do
    read(50,*)
    do icm=1,icmn ; read(50,*)
      do n=nn,1,-1 ; read(50,*)(dummy(m,n),m=1,mn) ; end do
    end do
    do jt=1,jtn ; do icm=1,icmn ;read(50,*)
      do n=nn,1,-1 ; read(50,*)(PDT_rate(m,n,jt,icm),m=1,mn) ; end do
    end do ;  end do
    deallocate(dummy)
    !--- get temporal total moment tensor
    jn = nint(maxval(Tr)/rtime) + jtn
    allocate (d_tq(jn,icmn),dummy(jn,icmn),dummy2(jn,icmn))
    d_tq =  0.
    do m=1, mn
      do n=1,nn
        do jt=1,jtn
          do icm=1,icmn
            j = nint( Tr(m,n) / rtime) + jt
            d_tq(j,icm) = d_tq(j,icm) + PDT_rate(m,n,jt,icm)
          enddo
        enddo
      enddo
    enddo
    dummy = d_tq ! dummy is now a pure potency-rate function for each icmn
    !--- get each temporal weight & allocate
    d_tq = abs(d_tq)
    if ( cr .eq. 0. .and. minval(d_tq) .eq. 0. ) call gimmick_avoid_zero(jn, icmn, d_tq)
    d_tq = d_tq / maxval(d_tq)
    do j=1,jn
      do icm = 1, icmn
        d_tq(j,icm) = max(d_tq(j,icm), cr)
      enddo
    enddo
    do m=1, mn
      do n=1,nn
        do jt=1,jtn
          do icm=1,icmn
            j = nint( Tr(m,n) / rtime ) + jt
            w(m,n,jt,icm) = 1./d_tq(j,icm)
          enddo
        enddo
      enddo
    enddo
    w = w / maxval(w)
    !--- output for check
    write(itera_str,*) itera
    do j=1,jn
      do icm=1,icmn
        d_tq(j,icm) = 1./d_tq(j,icm)
      end do
    end do
    d_tq = d_tq / maxval(d_tq) ! d_tq is now an actual given weight
    dummy2 = abs(dummy)
    cPmax = cr*maxval(dummy2) ! cPmax is Eq. (8)
    do j=1,jn
      do icm=1,icmn
        dummy2(j,icm) = max( dummy2(j,icm), cPmax ) ! dummy2 is now an adjusted absolute potency-rate function in Eq. (8)
      end do
    end do
    open(11, file = 'tmp_total_meca_'//trim(adjustl(itera_str))//'.dat', status = 'replace')
    do j = 1, jn
      write(11, '(i4, f7.2, 15f10.4)') j, j*rtime, (d_tq(j,icm), icm=1,icmn), (dummy2(j,icm), icm=1,icmn), (dummy(j,icm), icm=1,icmn)
    end do
    close(11)
    open(13, file='tmp_weight_'//trim(adjustl(itera_str))//'.dat', status='replace')
    do m = 1, mn
      do n = 1, nn
        do jt = 1, jtn
          do icm = 1, icmn
            write(13,'(f10.4,5i6)') w(m,n,jt,icm), m, n, jt, icm
            if (w(m,n,jt,icm) > 1.) write(6, *) " !!WARNING!! Smoothness weight to be more than 1"
            if (w(m,n,jt,icm) < cr-0.0001) write(6, '(" !!WARNING!! Smoothness weight to be smaller than",f5.2)') cr
          end do
        end do
      end do
    end do
    close(13)
    deallocate(PDT_rate,Tr,d_tq,dummy,dummy2)
  ! end if
  return
end subroutine get_tmp_weight

subroutine gimmick_avoid_zero(jn, icmn, d_tq)
  implicit none
  integer, intent(in) :: jn, icmn
  real, intent(inout) :: d_tq(jn, icmn)
  integer             :: j, icm
  real                :: sec_min
  !--- seeking the second minimum
  sec_min = 1.0e+33
  do j=1,jn
    do icm=1,icmn
      if ( d_tq(j,icm) .gt. 0. ) sec_min = min( sec_min, d_tq(j,icm) )
    end do
  end do
  !--- replace zero with the second minimum
  do j=1,jn
    do icm=1,icmn
      if ( d_tq(j,icm) .eq. 0. ) d_tq(j,icm) = sec_min
    end do
  end do
  write(6, '("TmpSm weight replaces 0. with ", f9.4)'), sec_min
  return
end subroutine
