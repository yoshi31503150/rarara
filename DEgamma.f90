   module DATA01
     implicit none
     double precision ka,theta,pi
     complex(8),parameter::i=(0.d0,1.d0)
     end module DATA01

	 module DATA02
	 implicit none
	 complex(8),parameter::u=(0.5d0, 0.d0)
	 end module DATA02

   program main
	 use DATA01
	 implicit none
	 double precision the1,A
	 complex(8) z,DE
	 integer n,m
	 external DE

	 pi=3.14159265358979323846264338327950288d0
	 ka=10*pi

	 OPEN(10,file='DEgamma.TXT')

	  do n=-1799,1799,1
      the1=dfloat(n)*0.1d0
      theta=the1*(pi/180.d0)
      z=-2.d0*i*ka*(1.d0-cos(theta))

	  write(10,*) the1,z,DE(z)

	  end do
	  close(10)
	  stop
	  end program main

   function DE(z)
	 use DATA01
	 implicit none
	 integer j,num,m
	 real(8) h,nu,xjh,dxjh,jj
	 complex(8) DE,z,total

	 m=1
	 nu=1.d0
	 num=226
	 h=4.d-2
	 total=(0.d0,0.d0)
	 do j=-num/2,num/2,1
	  jj=dfloat(j)
	  xjh=dexp(pi*dsinh(jj*h)/2.d0)
	  dxjh=(pi*dcosh(jj*h)/2)*dexp(pi*dsinh(jj*h)/2.d0)
	  total=total+((((xjh)**(-nu/2.d0))*dexp(-xjh)*dxjh)/((xjh+z)**m))
	 end do

	 DE=total*h

	 return
	 end
