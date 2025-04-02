! Y. Yagi
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
subroutine  writesol(mn,nn,jtn,jtn0,l_id_m,icmn,x,tr,var,abic,strike,dip,slip, &
               m0,n0,xx,yy,rtime,beta,sigx,depth,stcd,comp,sigw,  &
               dt,cp,jn,vr,nsurface,rslip,ylength,alpha1,alpha2)
  !
  ! The format of the output file was set referring to the format of FFI program
  !  developed by Prof. S. Yoshida (ERI, University of Tokyo).
  !
  implicit none
  integer,intent(in) :: mn,nn,jtn0,icmn,m0,n0,jn,nsurface
  integer,intent(in) :: jtn(mn,nn),l_id_m(mn,nn,jtn0,icmn)
  real,intent(in)    :: x(*),tr(mn,nn),var,abic,strike,dip,slip,xx,yy,rtime,beta(3), &
                        sigx(*),depth,cp(jn),sigw(jn),dt(jn),          &
                        vr,rslip,ylength,alpha1,alpha2
  character,intent(in) :: stcd(jn)*10,comp(jn)*4
  real:: xmo(6),TAXS(3)
  real:: moment,Mw
  real,allocatable:: PDT_comp(:,:,:),PD_dist(:,:),PDT_rate(:,:,:,:),E_PDT(:,:,:,:),  &
                      rake(:,:)
  real,allocatable :: srigid(:,:),samp(:,:)
  !----
  integer :: m,n,icm,jt,l,i,j
  real :: disp_length, disp_slip
  real :: dis1,dis2,wk,slip0,alat,alon
  real :: strike_s,dip_s,slip_s,rigid,AMO_s,DMO_s
  !----
  allocate ( PDT_comp(mn,nn,icmn),PD_dist(mn,nn),PDT_rate(mn,nn,icmn,jtn0),  &
                         E_PDT(mn,nn,icmn,jtn0),rake(mn,nn))
  allocate (srigid(mn,nn),samp(mn,nn))
  !----
  PDT_comp = 0.
  PDT_rate = 0.
  E_PDT = 0.
  do m=1,mn;do n=1,nn;
     do icm=1,icmn;do jt=1,jtn(m,n)
        l=l_id_m(m,n,jt,icm)
        PDT_comp(m,n,icm)=PDT_comp(m,n,icm)+x(l)
        PDT_rate(m,n,icm,jt)=x(l)
        E_PDT(m,n,icm,jt)=sigx(l)
     end do;end do
     if(icmn.eq.5)then
        xmo = 0.
        xmo(1:icmn) = PDT_comp(m,n,1:icmn)
        CALL D_CP(xmo,strike_s,dip_s,slip_s,AMO_s,DMO_s,TAXS)
        PD_dist(m,n) = AMO_s
        rake(m,n) = slip_s
     endif
     if(icmn.eq.2)then
        PD_dist(m,n)= disp_length(PDT_comp(m,n,1),PDT_comp(m,n,2),slip,rslip)
        rake(m,n) = disp_slip(PDT_comp(m,n,1),PDT_comp(m,n,2),slip,rslip)
     endif
     if(icmn.eq.1)then
        PD_dist(m,n)= PDT_comp(m,n,1)
        rake(m,n) = slip
     end if
  end do;end do
  !------------
  open(38,file='rigid_amp.info')
  do i=1,mn*nn
    read(38,*) n, m, srigid(m,n),samp(m,n)
  enddo
  close(38)
  rigid = sum(srigid) / real(mn*nn)
  !------------
  dis1 = 0.
  dis2 = 0.
  if(icmn.eq.5)then
    xmo = 0.
    do n=1,nn
       do m=1,mn
          wk = xx * yy * samp(m,n) * srigid(m,n) * 1.e15
          xmo(1:icmn) = xmo(1:icmn) + PDT_comp(m,n,1:icmn)* wk
       enddo
    enddo
    CALL D_CP(xmo,strike_s,dip_s,slip_s,AMO_s,DMO_s,TAXS)
    moment = AMO_s
    slip0  = 0.
  else
    do n=1,nn
       do m=1,mn
          wk = xx * yy * samp(m,n)
          !if(nsurface.eq.1 .and. n .eq. nn ) wk = xx * yy * 0.5 + xx * ylength
          dis1=dis1+PDT_comp(m,n,1)*srigid(m,n) * wk
          if(icmn.eq.2)then
             dis2=dis2+PDT_comp(m,n,2)*srigid(m,n) * wk
          end if
       end do
    end do
    slip0 = disp_slip(dis1,dis2,slip,rslip)
    moment=disp_length(dis1,dis2,slip,rslip)*1.e15
  endif
  Mw=alog10(Moment) / 1.5 - 6.06667
  !----
  alat= 0. ;alon = 0.
  open(38,file="epicenter.dat")
  read(38,*,err=997)alat,alon
997 close(38)

  !----
  write(40,*)"Moment(Nm),Mw,Rigid(GPa),Lat.,Lon.,Depth(km),Vr(km/sec),nsurface"
  write(40,'(e11.4,1x,f5.2,1x,f7.2,1x,2f10.3,1x,f7.1,1x,e10.3,1x,i2)') Moment,Mw,Rigid,alat,alon,depth ,vr, nsurface
  if(var>0.) write(6,10) Moment,Mw
10 Format(" Seismic Moment : ",e11.4," [Nm],  Mw : ",f5.2)
  write(40,*)" Strike, Dip, Slip, slip_input, rslip  "
  write(40,'(6(f7.2,3x))')strike,dip,slip0,slip,rslip
  write(40,*)"xx,  yy,  mn,  nn,  m0, n0, Raisetime, jtn icmn "
  write(40,'(2(f6.2,x),4(i3,x),f6.2,x,i3,x,i3,x,f7.2)') xx,yy,mn,nn,m0,n0,rtime,jtn0,icmn,ylength
  write(40,*)" Variance, ABIC, beta(1),beta(2),beta(3),alpha1,alpha2 "
  write(40,'(f8.5,x,f13.1,x,3(e13.5,x),2(e12.4,x))') var,abic,beta(1),beta(2),beta(3),alpha1,alpha2
  if(var > 0.)   write(6,*)" "
  if(var > 0.)   write(6,*)" Total Slip for each sub-fault "
  if(var > 0.)  write(6,*)" "
  write(40,*)" Total Slip for each sub-fault "
  if(mn.eq.1 .and. nn.eq.1) then
     do n=nn,1,-1
        write(40,'(e12.3,x)') (PD_dist(m,n)/500.,m=1,mn)
     end do
  else
     do n=nn,1,-1
   if(var > 0.)  write(6,'(30(f5.2,x))')  (PD_dist(m,n),m=1,mn)
        write(40,'(30(f9.4,x))') (PD_dist(m,n),m=1,mn)
     end do
  endif
  write(40,*)" Slip angel for each sub-fault "
  do n=nn,1,-1
     !         write(6,'(30(f5.1,x))')
     !     &        (Tr(m,n),m=1,mn)
     write(40,'(30(f7.1,1x))')(rake(m,n),m=1,mn)
  end do
  !
  write(40,*)" Start Time for each sub-fault "
  do n=nn,1,-1
     !         write(6,'(30(f5.1,x))')
     !     &        (Tr(m,n),m=1,mn)
     write(40,'(30(f6.2,1x))')(Tr(m,n),m=1,mn)
  end do
  write(40,*)" Total Slip Vecoter for each sub-fault "
  do icm=1,icmn
     !         write(6,*)" Vecoter ",icm
     write(40,*)" Vector ",icm
     do n=nn,1,-1
        if(mn.eq.1 .and. nn.eq.1) then
           write(40,'(e12.3)') (PDT_comp(m,n,icm)/500.,m=1,mn)
        else
           write(40,'(30(f9.4,x))') (PDT_comp(m,n,icm),m=1,mn)
        endif
     end do
  end do
  !
  do jt=1,jtn0
     do icm=1,icmn
        write(40,*)"JT, ICM : ",jt,",  ",icm
        do n=nn,1,-1
           if(mn.eq.1 .and. nn.eq.1) then
              write(40,'(e12.3)')  (PDT_rate(m,n,icm,jt),m=1,mn)
           else
              write(40,'(30(f9.4,x))')  (PDT_rate(m,n,icm,jt),m=1,mn)
           endif
        end do
     end do
  end do
  do jt=1,jtn0
     do icm=1,icmn
        write(40,*)"Error JT, ICM : ",jt,",  ",icm
        do n=nn,1,-1
           write(40,'(30(f9.5,x))')                                 &
                (E_PDT(m,n,icm,jt),m=1,mn)
        end do
     end do
  end do
  write(40,'(a)')" Station information "
  write(40,*)"stcd, comp, cp,  sigw,  dt, f1, f2, f3"

  do j = 1,jn
     write(40,'(a10,1x,a3,2f6.2,f10.4,f7.1,f7.3)')    &
          stcd(j),comp(j),cp(j),sigw(j),dt(j)
  end do
  close(40)
  deallocate (PDT_comp,PD_dist,PDT_rate,E_PDT,srigid,samp)
  return
END subroutine writesol
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
