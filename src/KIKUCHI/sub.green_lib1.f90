!  The sourutines were developed based on sourutines written by 
!  Prof. Kikuchi and Prof. Kanamori. Prof. Kanamori has kindly agreed 
!  to release to the public a modified version. 
!  Conversion to Fortran90 and minor modifications avoid unexpected 
!     problems during parallelization and optimisation. (Yuji Yagi) 
!
!  The original programs are available at "Note on Teleseismic Body-Wave 
!     Inversion Program(https://wwweic.eri.u-tokyo.ac.jp/ETAL/KIKUCHI/).
!---------
subroutine bodyw(x,n,dt,ib,ic,f1,d1,a1,h,az,p,fc0,zQp,zQs,zi,   &
     nl,vp,vs,den,dep,nl1,vp1,vs1,den1,dep1,nl2,vp2,vs2,den2,dep2)
  !===================================================*
  !  Calculate synthetic waveforms                    *
  !    Type of body wave  IB= 1/2/3/4: P/SV/SH/PP     *
  !    Component          IC= 1/2/3: UD/NS/EW (IB=1/4)*
  !                       IC= 1/2  : UD/HR    (IB=2)  *
  !                       IC= any  : SH       (IB=3)  *
  !============================ Ver.900715 ===========*
  ! Modification:                                     *
  ! 1) An isotropic component of M.T. is added        *
  !         for d1 > 360(degree)        -900531      *
  ! 2) Source layer # is determined in this subroutine*
  !      independently from a reference point -900715 *
  !===================================================*
  ! Modefication:
  !    Changed to improve usability(Y. Yagi)
  !===================================================*
  implicit none
  real :: pi = 3.141593
  !-------------------
  integer,intent(in)  :: n, ib, ic, nl, nl1, nl2
  real,intent(in)     :: dt, f1, d1, a1, h, az, p, fc0
  complex,intent(in)  :: zQp(n),zQs(n),zi(n)
  real,intent(in)     :: vp(nl), vs(nl), den(nl), dep(nl)
  real,intent(in)     :: vp1(nl1), vs1(nl1), den1(nl1), dep1(nl1)
  real,intent(in)     :: vp2(nl2), vs2(nl2), den2(nl2), dep2(nl2)
  real,intent(out)    :: x(n)
  !-------------------
  real :: df,dw,tl,hl,dh,ys,ts,yp,tp,dely,w,fc
  real :: pd,pu,svd,svu,shd,shu
  real :: tr0,tr1,tr2
  integer :: l, ll,i
  complex,allocatable :: z(:),zr0(:),zrpu(:),zrpd(:),zrsu(:),zrsd(:),zrpp(:),zdm(:)
  !-------------------
  allocate(z(n),zr0(n),zrpu(n),zrpd(n),zrsu(n),zrsd(n),zrpp(n),zdm(n))
  df = 1 /(dt * n)
  dw = df * 2 * pi
  tl = dt * n
  ! < source layer # >
  hl = 0.
  do l = 1, nl - 1
     hl = hl + dep(l)
     dh = h - hl
     if(dh.lt.0.) exit
  enddo
  ll = l
  hl = hl - dep(ll)
  ! < radiation pattern >
  call radp(az-f1,d1,a1,p,pd,pu,svd,svu,shd,shu,vp(ll),vs(ll))
  ! < structure effects: near-source & near-reciever >
  dh = h - hl
  if(ib.eq.3) then
     call reflsh(zrsu,zrsd,dw,n,p,vs,den,dep,nl,ll,tr1)
     call cnvrsh(zr0,dw,n,p,vs1,den1,dep1,nl1,tr2)
  else
     call refl(zrpu,zrpd,zrsu,zrsd,dw,n,p,vp,vs,den,dep,nl,ll,ib,tr1)
     if(ic.eq.1) call cnvr(zdm,zr0,dw,n,p,vp1,vs1,den1,dep1,nl1,ib,tr2)
     if(ic.ne.1) call cnvr(zr0,zdm,dw,n,p,vp1,vs1,den1,dep1,nl1,ib,tr2)
     !   PP-reflector
     call refl(zrpp,zdm,zdm,zdm,dw,n,p,vp2,vs2,den2,dep2,nl2,nl2,ib,tr0)
  endif
  ! < delay for depth phase >
  ys = sqrt(1 / vs(ll) **2 - p**2)
  ts = ys * dh
  if(ib.eq.1.or.ib.eq.4) then
     yp = sqrt(1 / vp(ll) **2 - p**2)
     tp = yp * dh
     if(ib.eq.1) dely = tr1 + tr2 - tp
     ! an additional delay of 10 sec is put for PP wave:
     if(ib.eq.4) dely = tr1 + tr2 - tp - 10.
  else
     dely = tr1 + tr2 - ts
  endif
  !-------
  do i = 1, n / 2
     w = dw *(i - 1)
     ! P or PP wave
     if(ib.eq.1.or.ib.eq.4) then
        fc = fc0 /(4 * pi * den(nl) * vp(nl) **3)
        z(i) = zrpd(i) * pd * cexp(cmplx(0., + w * tp) )
        !   exclude pP & sP phases outside the time window
        if(ll.eq.nl.and.tp * 2.0.ge.tl) goto 11
        z(i) = z(i) + zrpu(i) * pu * cexp(cmplx(0., - w * tp) ) &
             - zrsd(i) * svd * cexp(cmplx(0., + w * ts) ) - zrsu (i)   &
             * svu * cexp(cmplx(0., - w * ts) )
11      z(i) = z(i) * zQp(i)
        if(ib.eq.1) goto 50
        !  PP-reflector & additional q & hilbert-transform
        z(i) = z(i) * zrpp(i) * zQp(i) * cmplx(0., 1.)
        ! SV-wave
     elseif(ib.eq.2) then
        fc = fc0 /(4 * pi * den(nl) * vs(nl) **3) 
        z(i) = zrsd(i) * svd * cexp(cmplx(0., + w * ts) ) 
        !   exclude pS & sS phases outside the time window
        if(ll.eq.nl.and.ts * 2.0.ge.tl) goto 21
        z(i) = z(i) - zrpu(i) * pu * cexp(cmplx(0., - w * tp) ) &
             - zrpd(i) * pd * cexp(cmplx(0., + w * tp) ) + zrsu (i)    &
             * svu * cexp(cmplx(0., - w * ts) )
21      z(i) = z(i) * zQs(i)
        ! SH-wave
     elseif(ib.eq.3) then
        fc = fc0 /(4 * pi * den(nl) * vs(nl) **2 * vs(ll) ) 
        z(i) = zrsd(i) * shd * cexp(cmplx(0., + w * ts) ) 
        !   exclude ss phase outside the time window
        if(ll.eq.nl.and.ts * 2.0.ge.tl) goto 31
        z(i) = z(i) + zrsu(i) * shu * cexp(cmplx(0., - w * ts) ) 
31      z(i) = z(i) * zQs(i)
     endif
50   z(i) = fc * z(i) * zr0(i) * zi(i) * cexp(cmplx (0., w * dely) )
     if(i.eq.1) cycle
     z(n + 2 - i) = conjg(z(i) ) 
  end do
  z(n/2+1) = 0
  call cfft(z,n,1)
  x(1:n) = z(1:n) * df
  return
end subroutine bodyw
!                                                                       
!                                                                       
!                                                                       
subroutine cnvr(zr1, zr2, dw, n, p, vp, vs, den, dep, nl, ib, tr) 
  !  < Displacement-to-displacement conversion near the free surface >    
  !      P  = Ray parameter: sin(ih)/v                                    
  !  --------------------------------------                               
  !   For IB=1/4(incident P)                                              
  !      ZR1 = u(0)/uP(nl)(P to Horizontal component)                    
  !      ZR2 =-w(0)/uP(nl)(P to Vertical component(up))                  
  !   For IB=2(incident SV)                                             
  !      ZR1 = u(0)/uS(nl)(SV to Horizontal component)                   
  !      ZR2 =-w(0)/uS(nl)(SV to Vertical component(up))                 
  !  --------------------------------------                               
  !      TR = P-wave travel time from nl-th layer to top                  
  !  --------------------------------------    
  implicit none
  integer,intent(in) :: n, nl, ib
  real,intent(in)    :: dw, p
  real,intent(in)    :: vp(nl), vs(nl), den(nl),dep(nl) 
  real,intent(out)   :: tr
  complex,intent(out):: zr1(n), zr2(n)
  !-------------------
  integer  :: i,j,m
  real     :: p2,gm,vp2,y1,ra,vs2,y2,rb,w,pm,cpm,rmc2,qm,cqm
  complex  :: zspm,zsqm,zdet,zj1,zj2,zj3,zj4
  complex  :: ze(4,4), za(4,4), zaa(4,4), za1(4,4), zj(4,4)                                             
  !-------------------
  ze = 0 ; za= 0 ; zaa = 0 ; za1 = 0 ; zj  = 0
  p2 = p**2 
  tr = 0 
  do m = 1, nl - 1 
     tr = tr + dep(m) * sqrt(1 / vp(m) **2 - p2) 
  end do
  ! e-1 matrix for nl layers                                              
  gm = 2 *(p * vs(nl) ) **2 
  vp2 = 1 / vp(nl) **2 
  y1 = sqrt(vp2 - p2) 
  ra = y1 / p 
  vs2 = 1 / vs(nl) **2 
  y2 = sqrt(vs2 - p2) 
  rb = y2 / p 
  ze(1, 1) = - 2 *(vs(nl) / vp(nl) ) **2 
  ze(1, 3) = 1 /(den(nl) * vp(nl) **2) 
  ze(2, 2) = (gm - 1) /(vp(nl) **2 * ra * p2) 
  ze(2, 4) = 1 /(den(nl) * vp(nl) **2 * ra) 
  ze(3, 1) = (gm - 1) /(gm * rb) 
  ze(3, 3) = - p2 /(den(nl) * gm * rb) 
  ze(4, 2) = 1 
  ze(4, 4) = p2 /(den(nl) * gm) 
  !                                                                       
  do i = 1, n / 2 
     w = (i - 1) * dw 
     zaa = 0
     do  j = 1, 4 
        zaa(j, j) = 1
     end do
     do m = 1, nl - 1 
        vp2 = 1 / vp(m) **2 
        y1 = sqrt(vp2 - p2) 
        pm = y1 * w * dep(m) 
        cpm = cos(pm) 
        zspm = cmplx(0., sin(pm) ) 
        ra = y1 / p 
        rmc2 = den(m) / p**2 
        if(vs(m) .ne.0.) then 
           vs2 = 1 / vs(m) **2 
           y2 = sqrt(vs2 - p2) 
           gm = 2 *(p * vs(m) ) **2 
           qm = y2 * w * dep(m) 
           cqm = cos(qm) 
           zsqm = cmplx(0., sin(qm) ) 
           rb = y2 / p 
           ! a-matrix                                                              
           za(1, 1) = gm * cpm -(gm - 1) * cqm 
           za(1, 2) = (gm - 1) / ra * zspm + gm * rb * zsqm 
           za(1, 3) = -(1 / rmc2) *(cpm - cqm) 
           za(1, 4) = (1 / rmc2) *(zspm / ra + rb * zsqm) 
           za(2, 1) = -(gm * ra * zspm +(gm - 1) / rb * zsqm) 
           za(2, 2) = -(gm - 1) * cpm + gm * cqm 
           za(2, 3) = (1 / rmc2) *(ra * zspm + zsqm / rb) 
           za(2, 4) = za(1, 3) 
           za(3, 1) = rmc2 * gm *(gm - 1) *(cpm - cqm) 
           za(3, 2) = rmc2 *((gm - 1) **2 / ra * zspm + gm**2 * rb * zsqm)                                               
           za(3, 3) = za(2, 2) 
           za(3, 4) = za(1, 2) 
           za(4, 1) = rmc2 *(gm**2 * ra * zspm +(gm - 1) **2 /  rb * zsqm)                                               
           za(4, 2) = za(3, 1) 
           za(4, 3) = za(2, 1) 
           za(4, 4) = za(1, 1) 
           ! water layer                                                           
        else 
           za(1, 1) = 1 
           za(1, 2) = - zspm / ra 
           za(1, 3) = -(1 / rmc2) * cpm 
           za(1, 4) = (1 / rmc2) * zspm / ra 
           za(2, 1) = 0 
           za(2, 2) = cpm 
           za(2, 3) = ra / rmc2 * zspm 
           za(2, 4) = za(1, 3) 
           za(3, 1) = 0 
           za(3, 2) = rmc2 / ra * zspm 
           za(3, 3) = za(2, 2) 
           za(3, 4) = za(1, 2) 
           za(4, 1) = 0 
           za(4, 2) = 0 
           za(4, 3) = 0 
           za(4, 4) = 0 
        endif
        !                                                                       
        za1 = zaa 
        call prod(za, za1, zaa, 4) 
     end do
!                                                                       
     ! j-matrix for nl layers                                                
     call prod(ze, zaa, zj, 4) 
     !                                                                       
     zj1 = zj(4, 2) - zj(3, 2) 
     zj2 = zj(3, 1) - zj(4, 1) 
     zj3 = zj(2, 2) - zj(1, 2) 
     zj4 = zj(1, 1) - zj(2, 1) 
     zdet = zj4 * zj1 - zj3 * zj2 
     ! conversion at the free surface                                        
     if(ib.eq.1.or.ib.eq.4) then 
        !   for p-wave                                                          
        zr1(i) = - zj1 * 2 /(zdet * vp(nl) * p) 
        zr2(i) = + zj2 * 2 /(zdet * vp(nl) * p) 
     else 
        !   for s-wave                                                          
        zr1(i) = zj3 /(zdet * vs(nl) * p) 
        zr2(i) = zj4 /(zdet * vs(nl) * p) 
     endif
     if(i.eq.1) cycle
     zr1(n + 2 - i) = conjg(zr1(i) ) 
     zr2(n + 2 - i) = conjg(zr2(i) ) 
  end do
  zr1(n / 2 + 1) = 0 
  zr2(n / 2 + 1) = 0
  return
end subroutine cnvr
!                                                                       
!                                                                       
!                                                                       
subroutine cnvrsh(zrc, dw, n, p, vs, den, dep, nl, tr) 
  !  < Displacement-to-displacement conversion near the free surface >    
  !      P  = Ray parameter: sin(ih)/v                                    
  !  --------------------------------------                               
  !      ZRC = v(0)/vSH(nl)(SH to Horizontal component)                  
  !  --------------------------------------                               
  !      TR = S-wave travel time from nl-th layer to top                  
  !  -------------------------------------- 
  implicit none
  integer,intent(in)   :: n, nl
  real,intent(in)      :: dw, p
  real,intent(in)      :: vs(nl), den(nl), dep(nl) 
  real,intent(out)     :: tr
  complex,intent(out)  :: zrc(n)
  !-------------------
  integer :: i,j,m
  real    :: p2,w,vs2,y2,rgb,qm,cqm,zsqm
  complex :: zdet
  complex :: za(2, 2), zaa(2, 2), za1(2, 2) 
  !-------------------
  za= 0 ; zaa = 0 ; za1 = 0
  p2 = p**2 
  tr = 0 
  do m = 1, nl - 1 
     if(vs(m) .eq.0.) cycle
     tr = tr + dep(m) * sqrt(1 / vs(m) **2 - p2) 
  end do
  do i = 1, n / 2 
     w = (i - 1) * dw 
     zaa = 0
     do j = 1, 2 
        zaa(j, j) = 1
     end do
     do m = 1, nl - 1 
        vs2 = 1 / vs(m) **2 
        y2 = sqrt(vs2 - p2) 
        rgb = den(m) * vs(m) **2 * y2 / p 
        qm = y2 * w * dep(m) 
        cqm = cos(qm) 
        zsqm = cmplx(0., sin(qm) ) 
        ! a-matrix                                                              
        za(1, 1) = cqm 
        za(1, 2) = 1 / rgb * zsqm 
        za(2, 1) = rgb * zsqm 
        za(2, 2) = cqm 
        !                                                                       
        za1 = zaa 
        call prod(za, za1, zaa, 2) 
     end do
     !                                                                       
     vs2 = 1 / vs(nl) **2 
     y2 = sqrt(vs2 - p2) 
     rgb = den(nl) * vs(nl) **2 * y2 / p 
     zdet = zaa(1, 1) + zaa(2, 1) / rgb 
     zrc(i) = 2 / zdet 
     if(i.eq.1) cycle
     zrc(n + 2 - i) = conjg(zrc(i) ) 
  end do
  zrc(n / 2 + 1) = 0 
end subroutine cnvrsh
!                                                                       
!                                                                       
!                                                                       
subroutine refl(zrpu, zrpd, zrsu, zrsd, dw, n, p, vp, vs, den, dep, nl, l, ib, tr)                                               
!  < Haskell(1953,BSSA;1962;JGR)'s matrix >                             
!      IB = 1/2 for P/SV                                                
!      P  = Ray parameter: sin(ih)/v                                    
!  --------------------------------------                               
!   For IB=1/4(incident P)                                              
!     ZRPU = delta(l)' /delta(nl)''                                     
!     ZRPD = delta(l)''/delta(nl)''                                     
!     ZRSU =2omega(l)' /delta(nl)''                                     
!     ZRSD =2omega(l)''/delta(nl)''                                     
!   For IB=2(incident SV)                                               
!     ZRPU = delta(l)' /2omega(nl)''                                    
!     ZRPD = delta(l)''/2omega(nl)''                                    
!     ZRSU = omega(l)' / omega(nl)''                                    
!     ZRSD = omega(l)''/ omega(nl)''                                    
!  --------------------------------------                               
!   delta,omega are coefficients of <<potentials>>                      
!  --------------------------------------                               
!     TR = P-wave travel time from l to nl layer                        
!  -------------------------------------- 
  implicit none
  integer,intent(in)     :: n, nl, l, ib
  real,intent(in)        :: dw, p
  real,intent(out)       :: tr
  real,intent(in)        :: vp(nl), vs(nl), den(nl), dep(nl) 
  complex,intent(out)    :: zrpu(n), zrpd(n), zrsu(n), zrsd(n)
  !-------------------
  integer  :: i,j,m
  real     :: p2,gm,vp2,y1,ra,vs2,y2,rb,w,pm,cpm,rmc2,qm,cqm
  complex  :: zspm,zsqm,zj1,zj2,zj3,zj4,zdet
  complex  :: zel(4,4), ze(4,4), za(4,4), zaa(4,4), za1(4,4), zj(4,4), zjl(4,4) 
  !-------------------
  zel = 0 ; ze = 0 ; za= 0 ; zaa = 0 ; za1 = 0 ; zj = 0 ; zjl = 0
  p2 = p**2 
  tr = 0 
  do m = l, nl - 1 
     tr = tr + dep(m) * sqrt(1 / vp(m) **2 - p2)
  end do
  ! e-1 matrix for l layers                                               
  gm = 2 *(p * vs(l) ) **2 
  vp2 = 1 / vp(l) **2 
  y1 = sqrt(vp2 - p2) 
  ra = y1 / p 
  vs2 = 1 / vs(l) **2 
  y2 = sqrt(vs2 - p2) 
  rb = y2 / p 
  zel(1, 1) = - 2 *(vs(l) / vp(l) ) **2 
  zel(1, 3) = 1 /(den(l) * vp(l) **2) 
  zel(2, 2) = (gm - 1) /(vp(l) **2 * ra * p2) 
  zel(2, 4) = 1 /(den(l) * vp(l) **2 * ra) 
  zel(3, 1) = (gm - 1) /(gm * rb) 
  zel(3, 3) = - p2 /(den(l) * gm * rb) 
  zel(4, 2) = 1 
  zel(4, 4) = p2 /(den(l) * gm) 
  ! e-1 matrix for nl layers                                              
  gm = 2 *(p * vs(nl) ) **2 
  vp2 = 1 / vp(nl) **2 
  y1 = sqrt(vp2 - p2) 
  ra = y1 / p 
  vs2 = 1 / vs(nl) **2 
  y2 = sqrt(vs2 - p2) 
  rb = y2 / p 
  ze(1, 1) = - 2 *(vs(nl) / vp(nl) ) **2 
  ze(1, 3) = 1 /(den(nl) * vp(nl) **2) 
  ze(2, 2) = (gm - 1) /(vp(nl) **2 * ra * p2) 
  ze(2, 4) = 1 /(den(nl) * vp(nl) **2 * ra) 
  ze(3, 1) = (gm - 1) /(gm * rb) 
  ze(3, 3) = - p2 /(den(nl) * gm * rb) 
  ze(4, 2) = 1 
  ze(4, 4) = p2 /(den(nl) * gm) 
  !                                                                       
  do i = 1, n / 2 
     w = (i - 1) * dw 
     zaa = 0
     do j = 1, 4 
        zaa(j, j) = 1
     end do
     if(l.eq.1) call prod(zel, zaa, zjl, 4) 
     do m = 1, nl - 1 
        vp2 = 1 / vp(m) **2 
        y1 = sqrt(vp2 - p2) 
        pm = y1 * w * dep(m) 
        cpm = cos(pm) 
        zspm = cmplx(0., sin(pm) ) 
        ra = y1 / p 
        rmc2 = den(m) / p2 
        if(vs(m) .ne.0.) then 
           vs2 = 1 / vs(m) **2 
           y2 = sqrt(vs2 - p2) 
           gm = 2 *(p * vs(m) ) **2 
           qm = y2 * w * dep(m) 
           cqm = cos(qm) 
           zsqm = cmplx(0., sin(qm) ) 
           rb = y2 / p 
           ! a-matrix                                                              
           za(1, 1) = gm * cpm -(gm - 1) * cqm 
           za(1, 2) = (gm - 1) / ra * zspm + gm * rb * zsqm 
           za(1, 3) = -(1 / rmc2) *(cpm - cqm) 
           za(1, 4) = (1 / rmc2) *(zspm / ra + rb * zsqm) 
           za(2, 1) = -(gm * ra * zspm +(gm - 1) / rb * zsqm) 
           za(2, 2) = -(gm - 1) * cpm + gm * cqm 
           za(2, 3) = (1 / rmc2) *(ra * zspm + zsqm / rb) 
           za(2, 4) = za(1, 3) 
           za(3, 1) = rmc2 * gm *(gm - 1) *(cpm - cqm) 
           za(3, 2) = rmc2 *((gm - 1) **2 / ra * zspm + gm**2 *  &
                rb * zsqm)                                               
           za(3, 3) = za(2, 2) 
           za(3, 4) = za(1, 2) 
           za(4, 1) = rmc2 *(gm**2 * ra * zspm +(gm - 1) **2 /   &
                rb * zsqm)                                               
           za(4, 2) = za(3, 1) 
           za(4, 3) = za(2, 1) 
           za(4, 4) = za(1, 1) 
           ! water layer                                                           
        else 
           za(1, 1) = 1 
           za(1, 2) = - zspm / ra 
           za(1, 3) = -(1 / rmc2) * cpm 
           za(1, 4) = (1 / rmc2) * zspm / ra 
           za(2, 1) = 0 
           za(2, 2) = cpm 
           za(2, 3) = ra / rmc2 * zspm 
           za(2, 4) = za(1, 3) 
           za(3, 1) = 0 
           za(3, 2) = rmc2 / ra * zspm 
           za(3, 3) = za(2, 2) 
           za(3, 4) = za(1, 2) 
           za(4, 1) = 0 
           za(4, 2) = 0 
           za(4, 3) = 0 
           za(4, 4) = 0 
        endif
        !                                                                       
        za1 = zaa 
        call prod(za, za1, zaa, 4) 
        ! j-matrix for l layers                                                 
        if(m.eq.l - 1) call prod(zel, zaa, zjl, 4) 
     end do
     !                                                                       
     ! j-matrix for nl layers                                                
     call prod(ze, zaa, zj, 4) 
     !                                                                       
     zj1 = zj(4, 2) - zj(3, 2) 
     zj2 = zj(3, 1) - zj(4, 1) 
     zj3 = zj(2, 2) - zj(1, 2) 
     zj4 = zj(1, 1) - zj(2, 1) 
     zdet = zj4 * zj1 - zj3 * zj2 
     ! propagator coefficients                                               
     if(ib.eq.1.or.ib.eq.4) then 
        zrpu(i) = (zj1 *(zjl(1, 1) + zjl(2, 1) ) + zj2 *        &
          (zjl (1, 2) + zjl(2, 2) ) ) / zdet                         
        zrpd(i) = (zj1 *(zjl(1, 1) - zjl(2, 1) ) + zj2 *        &
          (zjl (1, 2) - zjl(2, 2) ) ) / zdet                         
        zrsu(i) = 2 *(zj1 *(zjl(3, 1) + zjl(4, 1) ) + zj2 *    &
          (zjl (3, 2) + zjl(4, 2) ) ) / zdet                         
        zrsd(i) = 2 *(zj1 *(zjl(4, 1) - zjl(3, 1) ) + zj2 *    &
          (zjl (4, 2) - zjl(3, 2) ) ) / zdet                         
     else 
        zrpu(i) = (zj3 *(zjl(1, 1) + zjl(2, 1) ) + zj4 *        &
          (zjl (1, 2) + zjl(2, 2) ) ) / zdet / 2                     
        zrpd(i) = (zj3 *(zjl(1, 1) - zjl(2, 1) ) + zj4 *        &
          (zjl (1, 2) - zjl(2, 2) ) ) / zdet / 2                     
        zrsu(i) = (zj3 *(zjl(3, 1) + zjl(4, 1) ) + zj4 *        &
          (zjl (3, 2) + zjl(4, 2) ) ) / zdet                         
        zrsd(i) = (zj3 *(zjl(4, 1) - zjl(3, 1) ) + zj4 *        &
          (zjl (4, 2) - zjl(3, 2) ) ) / zdet                         
     endif
     if(i.eq.1) cycle
     zrpu(n + 2 - i) = conjg(zrpu(i) ) 
     zrpd(n + 2 - i) = conjg(zrpd(i) ) 
     zrsu(n + 2 - i) = conjg(zrsu(i) ) 
     zrsd(n + 2 - i) = conjg(zrsd(i) ) 
  end do
  zrpu(n / 2 + 1) = 0 
  zrpd(n / 2 + 1) = 0 
  zrsu(n / 2 + 1) = 0 
  zrsd(n / 2 + 1) = 0 
end subroutine refl
!
!
!
subroutine reflsh(zrsu, zrsd, dw, n, p, vs, den, dep, nl, l, tr) 
  !  < Haskell(1953,BSSA;1960,JGR) matrix for SH wave >                   
  !      P  = Ray parameter: sin(ih)/v                                    
  !  --------------------------------------                               
  !     ZRSU = v(l)' / v(nl)''                                            
  !     ZRSD = v(l)''/ v(nl)''                                            
  !  v is coefficient of <<displacement>>                                 
  !  --------------------------------------                               
  !     TR = S-wave travel time from l to nl layer                        
  !  --------------------------------------    
  !  --------------------------------------                               
  implicit none
  integer,intent(in)  :: n, nl, l
  real,intent(in)     :: dw, p
  real,intent(in)     :: vs(nl), den(nl), dep(nl) 
  real,intent(out)    :: tr
  complex,intent(out) :: zrsu(n), zrsd(n)
  !-------------------
  integer :: i,j,m
  real    :: p2,w,vs2,y2,rgb,qm,cqm,rgbl
  complex :: zsqm,zdet
  complex :: za(2,2), zaa(2,2), za1(2,2), zal(2,2)                                                        
  !-------------------
  za= 0 ; zaa = 0 ; za1 = 0 ; zal = 0
  p2 = p**2 
  tr = 0 
  do m = l, nl - 1 
     if(vs(m) .eq.0.) cycle
     tr = tr + dep(m) * sqrt(1 / vs(m) **2 - p2) 
  end do
  do i = 1, n / 2 
     w = (i - 1) * dw 
     zaa = 0
     do j = 1, 2 
        zaa(j, j) = 1
     end do
     if(l.ne.1) goto 3 
     zal(1, 1) = 1. 
     zal(2, 2) = 1. 
     zal(1, 2) = 0. 
     zal(2, 1) = 0. 
3    do m = 1, nl - 1 
        if(vs(m) .ne.0.) then 
           vs2 = 1 / vs(m) **2 
           y2 = sqrt(vs2 - p2) 
           rgb = den(m) * vs(m) **2 * y2 / p 
           qm = y2 * w * dep(m) 
           cqm = cos(qm) 
           zsqm = cmplx(0., sin(qm) ) 
           ! a-matrix                                                              
           za(1, 1) = cqm 
           za(1, 2) = 1 / rgb * zsqm 
           za(2, 1) = rgb * zsqm 
           za(2, 2) = cqm 
           ! water layer                                                           
        else 
           za(1, 1) = 1 
           za(1, 2) = 0 
           za(2, 1) = 0 
           za(2, 2) = 0 
        endif
        !                                                                       
        za1 = zaa 
        call prod(za, za1, zaa, 2) 
        ! j-matrix for l layers                                                 
        if(m.ne.l - 1) cycle
        zal = zaa  
     end do
     !                                                                       
     vs2 = 1 / vs(nl) **2 
     y2 = sqrt(vs2 - p2) 
     rgb = den(nl) * vs(nl) **2 * y2 / p 
     zdet = zaa(1, 1) + zaa(2, 1) / rgb 
     vs2 = 1 / vs(l) **2 
     y2 = sqrt(vs2 - p2) 
     rgbl = den(l) * vs(l) **2 * y2 / p 
     ! propagator coefficients                                               
     zrsu(i) = (zal(1, 1) - zal(2, 1) / rgbl) / zdet 
     zrsd(i) = (zal(1, 1) + zal(2, 1) / rgbl) / zdet 
     if(i.eq.1) cycle
     zrsu(n + 2 - i) = conjg(zrsu(i) ) 
     zrsd(n + 2 - i) = conjg(zrsd(i) ) 
  end do
  zrsu(n / 2 + 1) = 0 
  zrsd(n / 2 + 1) = 0 
end subroutine reflsh
!                                                                       
!                                                                       
!                                                                       
subroutine radp(az, d0, a0, p, pd, pu, svd, svu, shd, shu, vp, vs)
  !  < Radiation pattern >                                                
  !  -----------------------------------                                  
  !        PD,PU for P,p                                                  
  !        SD,SU for S,s(SV&SH)                                           
  !     Polarity                                                          
  !        SH : Clockwise                                                 
  !        SV : SH x P                                                    
  !  -----------------------------------                                  
  !  May 29,1990                                                          
  !    adding radiation pattern for isotropic moment-tensor               
  !        for D0 > 360                                                   
  !  ----------------------------------- 
  implicit none
  real,intent(in)  :: az, d0, a0, p, vp, vs
  real,intent(out) :: pd, pu, svd, svu, shd, shu
  real :: rad = .0174533
  real :: dr0,ar0,thet,sih,cih,sih2,cih2
  real :: c1,c2,c3,sv1,sv2,sv3,sh1,sh2
  real :: the2,dr02,sth,cth,sth2,cth2,sd0,cd0,sd02,cd02,ca0,sa0
  real :: a1,a2,a3,a4,a5
  if(d0.gt.360.) then
     pd = 1
     pu = 1
     svd = 0
     svu = 0
     shd = 0
     shu = 0
     return
  endif
  dr0 = d0 * rad
  ar0 = a0 * rad
  thet = az * rad
  sih = p * vp
  if(sih.ge.1.) cih = 0.
  if(sih.lt.1.) cih = sqrt(1 - sih**2)
  sih2 = p * vs
  cih2 = sqrt(1 - sih2**2)
  !     p2=p*p                                                            
  !   *** c parameters                                                    
  c1 = sih**2
  c2 = - 2. * sih * cih
  c3 = 2 - 3 * sih**2
  !    *** sv parameters                                                  
  sv1 = sih2 * cih2
  sv2 = 2 * sih2**2 - 1
  sv3 = - 3 * sih2 * cih2
  !    *** sh parameters                                                  
  sh1 = sih2
  sh2 = cih2
  !    ***** a parameters                                                 
  the2 = thet * 2.
  dr02 = dr0 * 2.
  sth = sin(thet)
  cth = cos(thet)
  sth2 = sin(the2)
  cth2 = cos(the2)
  sd0 = sin(dr0)
  cd0 = cos(dr0)
  sd02 = sin(dr02)
  cd02 = cos(dr02)
  ca0 = cos(ar0)
  sa0 = sin(ar0)
  a1 = sth2 * ca0 * sd0 + .5 * cth2 * sa0 * sd02
  a2 = cth * ca0 * cd0 - sth * sa0 * cd02
  a3 = .5 * sa0 * sd02
  a4 = cth2 * ca0 * sd0 - .5 * sth2 * sa0 * sd02
  a5 = sth * ca0 * cd0 + cth * sa0 * cd02
  !                                                                       
  pd = a1 * c1 + a2 * c2 + a3 * c3
  pu = a1 * c1 - a2 * c2 + a3 * c3
  svd = a1 * sv1 + a2 * sv2 + a3 * sv3
  svu = - a1 * sv1 + a2 * sv2 - a3 * sv3
  shd = a4 * sh1 + a5 * sh2
  shu = a4 * sh1 - a5 * sh2
  return
end subroutine radp
!
!
!
subroutine instg(z, n, dw, id, zp, zz, izp, izz, a0, ip)
  !===============================================                        
  !   << Response of GDSN instrument >>          *                        
  !     ID = 0/1/2 for delta/step/ramp function  *  
  !===============================================                        
  implicit none
  integer,intent(in)  :: n, id, izp, izz, ip
  real,intent(in)     :: dw, a0
  complex,intent(in)  :: zp(50), zz(50)
  complex,intent(out) :: z(n)
  complex :: zw0,zzz
  integer :: m2,m22,i,j
  !----
  m2 = n / 2
  m22 = n + 2
  do i = 2, m2
     zw0 = cmplx(0., dw *(i - 1) )
     zzz = 1.
     do j = 1, izp
        zzz = zzz /(zw0 - zp(j) )
     end do
     do j = 1, izz
        zzz = zzz *(zw0 - zz(j) )
     end do
     z(i) = a0 * zzz * zw0**(ip - id)
     z(m22 - i) = conjg(z(i) )
  end do
  return
end subroutine instg
!                                                                       
!                                                                       
!                                                                       
subroutine prod(za, zb, zc, n)
  ! product of matrix                                                     
  implicit none
  integer,intent(in)   :: n
  complex,intent(in)   :: za(n,n), zb(n,n)
  complex,intent(out)  :: zc(n,n)
  integer :: j,j1,j2
  !-------------------
  do j1 = 1, n
     do j2 = 1, n
        zc(j1, j2) = 0
        do j = 1, n
           zc(j1, j2) = zc(j1, j2) + za(j1, j) * zb(j, j2)
        end do
     end do
  end do
  return
end subroutine prod
!                                                                       
!                                                                       
!                                                                       
subroutine Qf(z, n, tq, df)
  !  < Q-filter >                                                         
  !   modification is made for imaginary part
  !    following Prof. Kanamori's suggestion  -01/06/29
  implicit none
  integer,intent(in)  :: n
  real,intent(in)     :: tq, df
  complex,intent(out) :: z(n)
  !-------------------
  integer:: i
  real :: pi = 3.141593
  real :: s2,f
  complex  :: z1
  !-------------------
  ! fn = df * n / 2. 
  s2 = 3.*(2.*pi/tq)
  do i = 2, n / 2
     f = df *(i - 1)
     ! z1 = cmplx(0., 2.*f*tq) * clog(cmplx(0., f/fn) ) 
     z1 = cmplx(0., 2.*f*tq) * clog(cmplx(0., f/s2) )
     z(i) = cexp(z1)
     z(n + 2 - i) = conjg(z(i) )
  enddo
  z(1) = 1.
  z(n / 2 + 1) = 0.
  return
end subroutine Qf
