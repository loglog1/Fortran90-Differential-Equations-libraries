subroutine PredictorCorrector(res, N, y_init, x_begin, x_end, f)
  implicit none
  integer,intent(in)::N
  double precision,intent(in)::y_init,x_begin,x_end
  double precision,intent(out)::res(N)
  double precision fx(N),Div,k1,k2,k3,k4,yDummy
  double precision,parameter::AB_a1=55./24.,AB_a2=-59./24.,AB_a3=37./24.,AB_a4=-9./24.
  double precision,parameter::AM_a1=9./24.,AM_a2=19./24.,AM_a3=-5./24.,AM_a4=1./24.
  integer i

  interface
    double precision function f(x_,y_)
      double precision x_,y_
    end function
  end interface

  res(1)=y_init
  Div= (x_end-x_begin)/N
  if(x_begin>x_end) then
    write(*,*)"[PredictorCorrector] error: x_begin is more than x_end. you must select x_begin < x_end."
    return
  end if

    do i=1,3
      k1=f(Div*i,res(i)); fx(i)=k1;
      k2=f(Div*i+Div/2.,res(i)+k1*Div/2.)
      k3=f(Div*i+Div/2.,res(i)+k2*Div/2.)
      k4=f(Div*i+Div,res(i)+k3*Div)
      res(i+1)=res(i)+Div/6.*(k1+2.*k2+2.*k3+k4)
    end do

    yDummy=0.
    do i=4,N-1
      fx(i)=f(Div*i,res(i));
      yDummy=res(i)+Div*( AB_a1*fx(i) + AB_a2*fx(i-1) + AB_a3*fx(i-2) + AB_a4*fx(i-3) );
      do while(abs(yDummy-res(i+1))>1e-5)
        fx(i+1)=f(Div*(i+1),yDummy)
        res(i+1)=res(i)+Div*( AM_a1*fx(i+1) + AM_a2*fx(i) + AM_a3*fx(i-1) + AM_a4*fx(i-2) )
        yDummy=res(i+1)
      end do
    end do

end subroutine
