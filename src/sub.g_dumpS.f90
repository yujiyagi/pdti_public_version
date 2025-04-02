! Y. Yagi
subroutine get_dump_sol(x,nmodel,dump)
  implicit none
  integer,intent(in) :: nmodel
  real,intent(in)    :: dump
  real,intent(inout) :: x(nmodel)
  real,allocatable:: zz1(:)
  real,parameter :: threshold = 0.02
  integer :: n,i
  real    :: xdif,znorm
  !----
  allocate(zz1(nmodel))
  open(15,file="x.vector",form='unformatted')
  read(15)(zz1(n),n=1,nmodel)
  close(15)
  !---------
  xdif = 0. ;  znorm = 0.
  do n = 1,nmodel
     xdif = xdif + (x(n)-zz1(n))**2
     znorm =znorm +  zz1(n)**2
  enddo
  write(6,'("Xdif : ",f10.4)')xdif/znorm
  !---------
  if(xdif/znorm <= threshold) then
     open(17,file=".nstep.inf")
     write(17,*)"STOP"
     do i=1,nmodel
        write(17,'(2f10.3)')x(i),zz1(i)
     enddo
     close(17)
     !---------
  ! if(xdif/znorm > 0.005 ) then ! comment out and correct as below (2019.06.14 smz)
  elseif(xdif/znorm > threshold ) then
     x(1:nmodel) = zz1(1:nmodel) + (x(1:nmodel) - zz1(1:nmodel)) * dump
  endif
  deallocate(zz1)
  return
end subroutine get_dump_sol
