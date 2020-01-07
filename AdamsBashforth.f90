subroutine AdamsBashforth(res, N, y_init, x_begin, x_end, f)
  implicit none
  integer,intent(in)::N
  double precision,intent(out)::res(N)
  double precision,intent(in)::x_begin,x_end,y_init
  double precision fx(N),Div,k1,k2,k3,k4
  double precision,parameter::a1=55./24.,a2=-59./24.,a3=37./24.,a4=-9./24.
  integer i

  interface
    double precision function f(x_,y_)
      double precision x_,y_
    end function
  end interface

  if(x_begin>x_end) then
    write(*,*) "[AdamsBashforth] error: x_begin is more than x_end. you must select x_begin < x_end."
    return
  end if
  Div= (x_end-x_begin)/N

  res(1)=y_init
  do i=1,3
    k1=f(Div*i,res(i));fx(i)=k1
    k2=f(Div*i+Div/2.,res(i)+k1*Div/2.)
    k3=f(Div*i+Div/2.,res(i)+k2*Div/2.)
    k4=f(Div*i+Div,res(i)+k3*Div)
    res(i+1)=res(i)+Div/6.*(k1+2.*k2+2.*k3+k4)
  end do

  do i=4,N-1
    fx(i)=f(Div*i,res(i))
    res(i+1)=res(i)+Div*( a1*fx(i) + a2*fx(i-1) + a3*fx(i-2) + a4*fx(i-3) )
  end do
end subroutine AdamsBashforth
