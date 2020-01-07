subroutine AdamsMoulton(res, N, y_init, x_begin, x_end, f)
  implicit none
  integer,intent(in)::N
  double precision res(N),y_init,x_begin,x_end,fx(N)
  double precision Div,k1,k2,k3,k4
  double precision,parameter::a1=9./24.,a2=19./24.,a3=-5./24.,a4=1./24.
  integer i

  interface
    double precision function f(x_,y_)
      double precision,intent(in)::x_,y_
    end function
  end interface

  if(x_begin>x_end) then
    write(*,*) "[AdamsMoulton] error: x_begin is more than x_end. you must select x_begin < x_end."
    return
  end if

  Div= (x_end-x_begin)/N
  res(1)=y_init

  do i=1,2
    k1=f(Div*i,res(i)); fx(i)=k1
    k2=f(Div*i+Div/2.,res(i)+k1*Div/2.)
    k3=f(Div*i+Div/2.,res(i)+k2*Div/2.)
    k4=f(Div*i+Div,res(i)+k3*Div)
    res(i+1)=res(i)+Div/6.*(k1+2.*k2+2.*k3+k4)
  end do

  do i=3,N-1
        k1=f(Div*i,res(i)); fx(i)=k1
        k2=f(Div*i+Div/2.,res(i)+k1*Div/2.)
        k3=f(Div*i+Div/2.,res(i)+k2*Div/2.)
        k4=f(Div*i+Div,res(i)+k3*Div)
        res(i+1)=res(i)+Div/6.*(k1+2.*k2+2.*k3+k4)
        fx(i+1)=f(Div*(i+1),res(i+1))
        res(i+1)=res(i)+Div*( a1*fx(i+1) + a2*fx(i) + a3*fx(i-1) + a4*fx(i-2) )
  end do

end subroutine AdamsMoulton
