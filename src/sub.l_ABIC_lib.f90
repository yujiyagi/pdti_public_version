! Y. Yagi
!==============================================
function get_hyper_para(x0,dx,i)
  implicit none
  real get_hyper_para
  real,intent(in)    :: x0, dx
  integer,intent(in) :: i
  real :: wk
  wk = x0 + dx * (i-1)
  get_hyper_para = 10. ** wk
  return
end function get_hyper_para
!==============================================
subroutine get_min_value_1(abic,it1,abic_m,i1_m)
  implicit none
  integer,intent(in)  :: it1
  real,intent(in)     :: abic(it1)
  real,intent(out)    :: abic_m
  integer,intent(out) :: i1_m
  integer :: i1
  abic_m = 1.e15
  do i1 = 1,it1
    if(abic_m > abic(i1) ) then
      i1_m = i1
      abic_m =  abic(i1)
    endif
  enddo
  return
end subroutine get_min_value_1
!==============================================
subroutine set_alpha_para0(ita2,a2_0,d_a2)
  implicit none
  integer,intent(out) :: ita2
  real,intent(out)    :: a2_0,d_a2
  ita2 = 9
  a2_0 =  -3. ; d_a2 = 1.
  return
end subroutine set_alpha_para0
!==============================================
subroutine set_beta_para0(itb1,b1_0,d_b1)
  implicit none
  integer,intent(out) :: itb1
  real,intent(out)    :: b1_0,d_b1
  itb1 = 15 
  b1_0 =  -2. ; d_b1 = 1.
  return
end subroutine set_beta_para0
!==============================================
!==============================================
subroutine set_beta_para1(itb1,b1_0,d_b1,ib1)
  implicit none
  integer,intent(out)  :: itb1
  integer,intent(in)   :: ib1
  real,intent(inout)   :: b1_0,d_b1
  itb1 = 3 
  b1_0 = b1_0 + d_b1*(ib1-1)
  d_b1 = d_b1 / 2.
  b1_0 = b1_0 - d_b1
  return
end subroutine set_beta_para1
!==============================================
subroutine set_alpha_para1(ita2,a2_0,d_a2,ia2)
  implicit none
  integer,intent(out) :: ita2
  integer,intent(in)  :: ia2
  real,intent(inout)  :: a2_0,d_a2
  ita2 = 3
  a2_0 = a2_0 + d_a2*(ia2-1)
  d_a2 = d_a2 / 2.
  a2_0 = a2_0 - d_a2
  return
end subroutine set_alpha_para1
