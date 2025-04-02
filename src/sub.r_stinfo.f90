! Y. Yagi
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
subroutine  get_stinfo(jn,stcd,comp,az,del)
  !  staion information
  implicit none
  integer,intent(in)  :: jn
  character,intent(in):: stcd(jn)*10,comp(jn)*4
  real(4),intent(out)::  az(jn),del(jn)
  !  working space
  character:: cha*1
  integer :: j,ndg,ic1,ic2
  real :: tg0,dtg,work

  !--  
  do j=1,jn
    ic1=index(stcd(j),' ') - 1 ; ic2 =index(comp(j),' ') - 1
    open(21,file='wave.grn/'//stcd(j)(1:ic1)//comp(j)(1:ic2)//cha(1), &
                                         form='unformatted')
    read(21);read(21);read(21)
    read(21)tg0,dtg,work,ndg,az(j),del(j)
    close(21)
  enddo
  return
END subroutine get_stinfo
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
