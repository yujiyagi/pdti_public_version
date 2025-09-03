program RunBP
  implicit none
  real                          :: f1, f2, np, dt, wr
  real, allocatable             :: weight(:), Fj(:), alpha(:), TempGij(:,:,:)
  real, allocatable             :: polarity(:), ampGij(:,:), Uij(:,:,:), Gij(:,:,:)
  real(kind=8), allocatable     :: stack(:,:,:)
  integer                       :: bpflag, potflag, nx, ny, x0, y0, wi, nj, j, t, x, y, acvelflag
  integer, allocatable          :: Tij(:,:)
  integer, parameter            :: nt=400*20, ng=40*20, ns=300*20!! nt > ns, assuming dt = 0.05 (20Hz)
  character(len=20), allocatable:: station_code(:), comp(:)  
  read(5, *) f1, f2, bpflag, potflag, np, acvelflag
  write(6, '("F1:", f5.2, "   F2:", f5.2, "   bpflag:", i2, "   potflag:", i2, "   &
       np:", f5.2, "   acvelflag:", i2)') f1, f2, bpflag, potflag, np, acvelflag
  open(10, file='i_greenf', status='old', action='read')
  read(10, *) wi, dt
  if (dt /= 0.05) then
     write(6, *) "Input dt = 0.05 in i_greenf. Stopped."
     stop
  end if
  read(10, *) wr, wr, nx, ny, x0, y0
  close(10)
  call get_nj(nj)
  allocate(station_code(1:nj), comp(1:nj), weight(1:nj))
  call read_station_info(station_code, comp,  weight, nj)
  allocate(alpha(1:nj), polarity(1:nj), stack(1:ns, 1:nx, 1:ny))
  alpha = 0
  polarity = 0
  stack = 0
  do j = 1, nj
     allocate(Fj(1:nt), TempGij(1:nt, 1:nx, 1:ny), Tij(1:nx, 1:ny), ampGij(1:nx, 1:ny), Gij(1:ng, 1:nx, 1:ny))
     Fj = 0
     TempGij = 0
     Tij = 0
     ampGij = 0
     Gij = 0
     call read_obs(station_code(j), comp(j), f1, f2, dt, Fj, alpha(j), nt, bpflag, potflag)
     call read_green(station_code(j),comp(j),f1,f2, TempGij, Tij, polarity(j), &
          nt, nx, ny, ampGij, bpflag, potflag, acvelflag)
     call shift_normalized_green(TempGij, Tij, alpha(j), Gij, nt, nx, ny, ng, bpflag, potflag, ampGij)
     deallocate(TempGij)
     allocate(Uij(1:nt, 1:nx, 1:ny))
     Uij = 0
     call combining_Uij(Fj, Gij, Uij, polarity(j), nt, ng, nx, ny, bpflag, potflag, ampGij)
     call NthRoot(nx, ny, Tij, x0, y0, ns, nt, weight(j), Uij, np, stack, dt)
     deallocate(Fj, Tij, ampGij, Gij, Uij)
     write(6, '(a5, i3, "/", i3)') station_code(j), j, nj
  end do
  call NthPower(stack, ns, nx, ny, np, nj)
  call OutputSol(stack, ns, nx, ny)
end program RunBP
