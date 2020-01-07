program test
  implicit none
  integer,parameter::N=100
  double precision y(N)
  double precision,parameter::x_begin=0.0,x_end=1.0,y_init=1.0
  integer i

  call AdamsMoulton(y,N,y_init,x_begin,x_end,f)
contains
  double precision function f(x,y)
    double precision,intent(in)::x,y
    f=x
  end function

end program test
