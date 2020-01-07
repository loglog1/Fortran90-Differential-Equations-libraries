subroutine RungeKutta(y,N,init,x_begin,x_end,f)
  implicit none

  interface
    real function f(x_,y_)
      real, intent(in) :: x_,y_
    end function
  end interface

  real,intent(in)::init,x_begin,x_end
  real y(N)
  real k1,k2,k3,k4,x,Div
  integer i,N

  if(x_begin>x_end) then
    write(*,*) "[RungeKutta] error: x_begin is more than x_end. you must select x_begin < x_end."
    return
  end if

  Div=(x_end-x_begin)/N
  x=x_begin

  y(1)=init
  do i=2,N
    x=x+Div
    k1=f(x , y(i-1))
    k2=f(x+Div/2. , y(i-1)+Div*k1/2.)
    k3=f(x+Div/2. , y(i-1)+Div*k2/2.)
    k4=f(x+Div , y(i-1)+Div*k3)
    y(i)=y(i-1)+(k1+2.*k2+2.*k3+k4)*Div/6
  end do
end subroutine
