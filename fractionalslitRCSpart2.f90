  module DATA01
  implicit none
  double precision theta,pi,f,epsilon0,mu0,omega,the0,theta0,ka,k,a,ab
  complex(8) A_1,A_2,B_1,B_2,C_1_nu0,C_2_nu0,C_1_nu1,C_2_nu1
  complex(8) ab01,ab11,ab2,ab3,ab4
  complex(8),parameter::i=(0.d0,1.d0)
  end module DATA01

	 module DATA02
	 implicit none
	 real(8) nu
	 end module DATA02

	 module DATA03
	 implicit none
	 integer n,m
	 real(8) c(100),d(100),xi(100),wi(100),jj,nn,cmax
	 complex(8) x0(100),x(100)
	 end module DATA03

	 module DATA04
	 implicit none
	 double precision x_shisu_nu0,x_shisu_nu1
	 integer z_shisu
   end module DATA04

   program main
	 use DATA01
	 use DATA02
	 use DATA03
	 use DATA04
	 implicit none
	 double precision the1,db,length,Sigma,alpha
	 complex(8) Phi_1,Phi_2
   complex(8) Uplus,Uminus,Vplus,Vminus
   complex(8) chi_1_nu0,chi_2_nu0,chi_1_nu1,chi_2_nu1
   complex(8) P_1_nu0,P_2_nu0,Q_1_nu1,Q_2_nu1
   complex(8) Mplus,Mminus,Kplus,Kminus
   complex(8) ETA_1_nu0,ETA_2_nu0,ETA_1_nu1,ETA_2_nu1
   complex(8) ZAI_nu0,ZAI_nu1
   complex(8) z,DE_nu0,DE_nu1

   integer the

   external Phi_1,Phi_2
   external Uplus,Uminus,Vplus,Vminus
   external chi_1_nu0,chi_2_nu0,chi_1_nu1,chi_2_nu1
   external P_1_nu0,P_2_nu0,Q_1_nu1,Q_2_nu1
   external Mplus,Mminus,Kplus,Kminus
   external ETA_1_nu0,ETA_2_nu0,ETA_1_nu1,ETA_2_nu1
   external ZAI_nu0,ZAI_nu1
   external DE_nu0,DE_nu1

	 f=1.d9
	 epsilon0=8.854187817d-12
	 mu0=1.2566370614d-6
	 pi=3.14159265358979323846264338327950288d0
	 omega=pi*f*2.d0
	 k=omega*dsqrt(epsilon0*mu0)
	 the0=60.00000000000000000000d0
	 theta0=the0*(pi/180.00000000000000000000d0)
	 ab=k*dcos(theta0)   !kcosθ0

	  write(*,*) 'nu=?(0<=nu<=1)'
	  read(*,*) nu
	  write(*,*) 'ka=?'
	  read(*,*) length
	  write(*,*)'m=?'
	  read(*,*) z_shisu

    ka=pi*length
	  a=ka/k
	  n=z_shisu+1
	  x_shisu_nu0=((nu-1.d0)*0.50d0)
	  x_shisu_nu1=nu/2.d0

  open(10,file='SlitBistaticRCSpart2(nu=0~1).TXT')

   do the=-1799,1799,1

   the1=dfloat(the)*0.1d0
   theta=the1*(pi/180.00000000000000000000d0)
   alpha=-k*dcos(theta)
   ab01=dcmplx(-ka*dcos(theta),0.00000000000000000000d0)
   !ab11=dcmplx(-k*dcos(theta),0.00000000000000000000d0)

 if(the.gt.-1800 .and. the.lt.(-1800+int(10.d0*the0)))then
		 Sigma=(cdabs(Phi_2(alpha)*k*dsin(theta)))**2
		 db=10.d0*dlog10(sigma)
	  else if(the.eq.(-1800+int(10.d0*the0)))then
		 Sigma=0.d0
		else if(the.gt.(-1800+int(10.d0*the0)) .and. the.lt.0)then
		 Sigma=(cdabs(Phi_2(alpha)*k*dsin(theta)))**2
		 db=10.d0*dlog10(sigma)
		else if(the.eq.0)then
		 Sigma=0.d0
		else if(the.gt.0 .and. the.lt.(1800-int(10.d0*the0)))then
		 Sigma=(cdabs(Phi_1(alpha)*k*dsin(theta)))**2
		 db=10.d0*dlog10(sigma)
		else if(the.eq.(1800-int(10.d0*the0)))then
		 Sigma=0.d0
		else
		 Sigma=(cdabs(Phi_1(alpha)*k*dsin(theta)))**2
		 db=10.d0*dlog10(sigma)
	  end if

   write(10,*) db

   end do
   close(10)

   open(18,file='Slitkakunin.dat',status='replace')

   write(18,*) Phi_1(alpha)
   write(18,*) Phi_2(alpha)
   write(18,*) Uplus(alpha)
   write(18,*) Uminus(alpha)
   write(18,*) Vplus(alpha)
   write(18,*) Vminus(alpha)
   write(18,*) C_1_nu0
   write(18,*) C_2_nu0
   write(18,*) C_1_nu1
   write(18,*) C_2_nu1
   write(18,*) chi_1_nu0(alpha)
   write(18,*) chi_2_nu0(alpha)
   write(18,*) chi_1_nu1(alpha)
   write(18,*) chi_2_nu1(alpha)
   write(18,*) A_1
   write(18,*) A_2
   write(18,*) B_1
   write(18,*) B_2
   write(18,*) P_1_nu0(alpha)
   write(18,*) P_2_nu0(alpha)
   write(18,*) Q_1_nu1(alpha)
   write(18,*) Q_2_nu1(alpha)
   write(18,*) Mplus(alpha)
   write(18,*) Mminus(alpha)
   write(18,*) Kplus(alpha)
   write(18,*) Kminus(alpha)
   write(18,*) ETA_1_nu0(alpha)
   write(18,*) ETA_2_nu0(alpha)
   write(18,*) ETA_1_nu1(alpha)
   write(18,*) ETA_2_nu1(alpha)
   write(18,*) ZAI_nu0(alpha)
   write(18,*) ZAI_nu1(alpha)
   write(18,*) DE_nu0(z)
   write(18,*) DE_nu1(z)

   close(18)

  write(*,*) alpha,x_shisu_nu0,x_shisu_nu1

   stop
   end program main

!----------Phi_1(alpha)------------------------
function Phi_1(alpha)
	   use DATA01
	   use DATA02
	   implicit none
     double precision alpha

	   complex(8) Phi_1,Uplus,Uminus,Vplus,Vminus,MPgamma,Utotal,Vtotal
	   external Uplus,Uminus,Vplus,Vminus

	   Utotal=Uplus(alpha)*cdexp(i*a*alpha)+Uminus(alpha)*cdexp(-i*a*alpha)
	   Vtotal=Vplus(alpha)*cdexp(i*a*alpha)+Vminus(alpha)*cdexp(-i*a*alpha)
	   MPgamma=cdexp(-i*pi*0.50d0)*dsqrt(k**2-alpha**2)

	    !式(5.6)参照
	   if(nu.eq.0)then
	    Phi_1=(Utotal*dsin(theta))/(-MPgamma)
	   else if(nu>0 .and. nu<1)then
	    Phi_1=((Utotal*dsin(theta))/(-MPgamma))+Vtotal
	   else
	    Phi_1=Vtotal
	    end if

	   return
	   end

!----------Phi_2(alpha)------------------------
      function Phi_2(alpha)
	   use DATA01
	   use DATA02
	   implicit none
     double precision alpha


	   complex(8) Phi_2,Uplus,Uminus,Vplus,Vminus,MPgamma,Utotal,Vtotal
	   external Uplus,Uminus,Vplus,Vminus

     Utotal=Uplus(alpha)*cdexp(i*a*alpha)+Uminus(alpha)*cdexp(-i*a*alpha)
	   Vtotal=Vplus(alpha)*cdexp(i*a*alpha)+Vminus(alpha)*cdexp(-i*a*alpha)
	   MPgamma=cdexp(-i*pi*0.50d0)*dsqrt(k**2-alpha**2)

	   !式(5.6)参照
	   if(nu.eq.0)then
	    Phi_2=(Utotal*dsin(theta))/(-MPgamma)
	   else if(nu>0 .and. nu<1)then
	    Phi_2=((Utotal*dsin(theta))/(-MPgamma))-Vtotal
	   else
	    Phi_2=-Vtotal
	    end if

	   return
	   end

!----------Uplus(alpha)------------------------
     function Uplus(alpha)
	   use DATA01
	   implicit none
     real(8) alpha

	   complex(8) Uplus,Mplus,ETA_2_nu0,ZAI_nu0
	   external Mplus,ETA_2_nu0,ZAI_nu0

	  call A1teisu
	  call A2teisu
	  call C1_nu0_teisu

	   !式(4.26)参照
	   Uplus=Mplus(alpha)*((-A_1/(Mplus(ab)*(alpha-ab)))+A_2*ETA_2_nu0(alpha)+C_1_nu0*ZAI_nu0(alpha))

	   return
	   end

!----------Uminus(alpha)------------------------
    function Uminus(alpha)
	   use DATA01
	   implicit none
     real(8) alpha

	   complex(8) Uminus,Mplus,Mminus,ETA_1_nu0,ZAI_nu0
	   external Mplus,Mminus,ETA_1_nu0,ZAI_nu0

	  call A1teisu
	  call A2teisu
	  call C2_nu0_teisu

	   !式(4.26)参照
	   Uminus=Mminus(alpha)*((A_2/(Mminus(ab)*(alpha-ab)))+A_1*ETA_1_nu0(-alpha)+C_2_nu0*ZAI_nu0(-alpha))

	   return
	   end

!----------Vplus(alpha)------------------------
    function Vplus(alpha)
	   use DATA01
	   use DATA02
	   implicit none
     real(8) alpha

	   complex(8) Vplus,Kplus,ETA_2_nu1,ZAI_nu1,ETA_1_nu1
	   external Kplus,ETA_2_nu1,ZAI_nu1,ETA_1_nu1

	  call B1teisu
	  call B2teisu
	  call C1_nu1_teisu

	   !式(4.50)参照
	   Vplus=Kplus(alpha)*(-B_1/(Kplus(ab)*(alpha-ab))+B_2*ETA_2_nu1(alpha)+C_1_nu1*ZAI_nu1(alpha))

	   return
	   end

!----------Vminus(alpha)------------------------
    function Vminus(alpha)
	   use DATA01
	   use DATA02
	   implicit none
     real(8) alpha

	   complex(8) Vminus,Kplus,Kminus,ETA_1_nu1,ZAI_nu1,ETA_2_nu1
	   external Kplus,Kminus,ETA_1_nu1,ZAI_nu1,ETA_2_nu1

	  call B1teisu
	  call B2teisu
	  call C2_nu1_teisu

	   !式(4.50)参照
	   Vminus=Kminus(alpha)*((B_2/(Kminus(ab)*(alpha-ab)))+B_1*ETA_1_nu1(-alpha)+C_2_nu1*ZAI_nu1(-alpha))

	   return
	   end


!----------C1_nu0_teisu(nu=0)------------------------
    subroutine C1_nu0_teisu
	   use DATA01
	   implicit none
	   complex(8) Mplus,chi_1_nu0,chi_2_nu0,ZAI_nu0
	   external Mplus,chi_1_nu0,chi_2_nu0,ZAI_nu0

	   !式(4.27)参照
	   C_1_nu0=(Mplus(k)*(chi_2_nu0(k)+Mplus(k)*ZAI_nu0(k)*chi_1_nu0(k)))/(1.d0-(Mplus(k)*ZAI_nu0(k))**2)

	   return
	   end

!----------C2_nu0_teisu(nu=0)------------------------
    subroutine C2_nu0_teisu
	   use DATA01
	   implicit none
	   complex(8) Mplus,chi_1_nu0,chi_2_nu0,ZAI_nu0
	   external Mplus,chi_1_nu0,chi_2_nu0,ZAI_nu0

	   !式(4.27)参照
	   C_2_nu0=(Mplus(k)*(chi_1_nu0(k)+Mplus(k)*ZAI_nu0(k)*chi_2_nu0(k)))/(1.d0-(Mplus(k)*ZAI_nu0(k))**2)

	   return
	   end

!----------C1_nu1_teisu(nu=1)------------------------
    subroutine C1_nu1_teisu
	   use DATA01
	   implicit none
	   complex(8) Kplus,chi_1_nu1,chi_2_nu1,ZAI_nu1
	   external Kplus,chi_1_nu1,chi_2_nu1,ZAI_nu1

!式(4.51)参照
	   C_1_nu1=(Kplus(k)*(chi_2_nu1(k)+Kplus(k)*ZAI_nu1(k)*chi_1_nu1(k)))/(1.d0-((Kplus(k)*ZAI_nu1(k))**2))

	   return
	   end

!----------C2_nu1_teisu(nu=1)------------------------
    subroutine C2_nu1_teisu
	   use DATA01
	   implicit none
	   complex(8) Kplus,chi_1_nu1,chi_2_nu1,ZAI_nu1
	   external Kplus,chi_1_nu1,chi_2_nu1,ZAI_nu1

	  !式(4.51)参照
	   C_2_nu1=(Kplus(k)*(chi_1_nu1(k)+Kplus(k)*ZAI_nu1(k)*chi_2_nu1(k)))/(1.d0-((Kplus(k)*ZAI_nu1(k))**2))

	   return
	   end

!----------chi_1_nu0(alpha)------------------------
    function chi_1_nu0(alpha)
	   use DATA01
	   implicit none
     real(8) alpha

	   complex(8) chi_1_nu0,P_1_nu0,ETA_2_nu0
	   external P_1_nu0,ETA_2_nu0

	   call A1teisu
	   call A2teisu

	   !式(4.23)参照
	   chi_1_nu0=A_1*P_1_nu0(alpha)+A_2*ETA_2_nu0(alpha)

	   return
	   end

!----------chi_2_nu0(alpha)------------------------
    function chi_2_nu0(alpha)
	   use DATA01
	   implicit none
     real(8) alpha

	   complex(8) chi_2_nu0,P_2_nu0,ETA_1_nu0
	   external P_2_nu0,ETA_1_nu0

	   call A1teisu
	   call A2teisu

	   !式(4.23)参照
	   chi_2_nu0=A_2*P_2_nu0(alpha)+A_1*ETA_1_nu0(alpha)

	   return
	   end

!----------chi_1_nu1(alpha)------------------------
    function chi_1_nu1(alpha)
	   use DATA01
	   implicit none
     real(8) alpha
	   complex(8) chi_1_nu1,Q_1_nu1,ETA_2_nu1
	   external Q_1_nu1,ETA_2_nu1

	   call B1teisu
	   call B2teisu

	   !式(4.47)参照
	   chi_1_nu1=B_1*Q_1_nu1(alpha)+B_2*ETA_2_nu1(alpha)

	   return
	   end

!----------chi_2_nu1(alpha)------------------------
    function chi_2_nu1(alpha)
	   use DATA01
	   implicit none
     real(8) alpha
	   complex(8) chi_2_nu1,Q_2_nu1,ETA_1_nu1
	   external Q_2_nu1,ETA_1_nu1

	   call B1teisu
	   call B2teisu

!式(4.47)参照
	   chi_2_nu1=B_2*Q_2_nu1(alpha)+B_1*ETA_1_nu1(alpha)


	   return
	   end

!----------A1teisu------------------------
    subroutine A1teisu
	  use DATA01
	  use DATA02
    implicit none

!式(2.102)参照
	   A_1=(cdexp(-i*ka*dcos(theta0))*(-k*dsin(theta0))&
     *(1.d0+cdexp(i*pi*nu)))/(2.d0*dsqrt(2.d0*pi))

     !write(*,*) A_1

	   return
	   end

!----------A2teisu------------------------
    subroutine A2teisu
	   use DATA01
	   use DATA02
     implicit none

	   !式(2.102)参照
	   A_2=(cdexp(i*ka*dcos(theta0))*(-k*dsin(theta0))&
     *(1.d0+cdexp(i*pi*nu)))/(2.d0*dsqrt(2.d0*pi))

	   return
	   end

!----------B1teisu------------------------
    subroutine B1teisu
	   use DATA01
	   use DATA02
     implicit none

!式(2.103)参照
	   B_1=(cdexp(-i*ka*dcos(theta0))*(1.d0-cdexp(i*pi*nu)))&
     /(i*2.d0*dsqrt(2.d0*pi))

	   return
	   end

!----------B2teisu------------------------
    subroutine B2teisu
	   use DATA01
	   use DATA02
     implicit none

!式(2.103)参照
	   B_2=(cdexp(i*ka*dcos(theta0))*(1.d0-cdexp(i*pi*nu)))&
     /(i*2.d0*dsqrt(2.d0*pi))

	   return
	   end

!----------P_1_nu0(alpha)------------------------
	  function P_1_nu0(alpha)
	   use DATA01
	   implicit none
     real(8) alpha
	   complex(8) P_1_nu0,Mplus
	   external Mplus

	   !式(4.2)参照
	   P_1_nu0=((1.d0/Mplus(alpha))-(1.d0/Mplus(ab)))/(alpha-ab)

	   return
	   end

!----------P_2_nu0(alpha)------------------------
	  function P_2_nu0(alpha)
	   use DATA01
	   implicit none
     real(8) alpha
	   complex(8) P_2_nu0,Mplus,Mminus
	   external Mplus,Mminus

	   !式(4.3)参照
	   P_2_nu0=((1.d0/Mplus(alpha))-(1.d0/Mminus(ab)))/(alpha+ab)

	   return
	   end

!----------Q_1_nu1(alpha)------------------------
	  function Q_1_nu1(alpha)
	   use DATA01
	   implicit none
     real(8) alpha
	   complex(8) Q_1_nu1,Kplus
	   external Kplus

	   !式(4.29)参照
	   Q_1_nu1=((1.d0/Kplus(alpha))-(1.d0/Kplus(ab)))/(alpha-ab)

	   return
	   end

!----------Q_2_nu1(alpha)------------------------
	  function Q_2_nu1(alpha)
	   use DATA01
	   implicit none
     real(8) alpha
	   complex(8) Q_2_nu1,Kplus,Kminus
	   external Kplus,Kminus

	   !式(4.30)参照
	   Q_2_nu1=((1.d0/Kplus(alpha))-(1.d0/Kminus(ab)))/(alpha+ab)

	   return
	   end

!----------Mplus------------------------
    function Mplus(alpha)
	  use DATA01
	  use DATA02
    implicit none
    real(8) alpha
	  complex(8) Mplus

	  !式(3.2)参照
	  Mplus=cdsqrt(1.d0+cdexp(i*pi*nu))/(2.d0*cdexp((i*pi*(nu+1.d0))*0.25d0)&
    *((k+alpha)**((nu-1.d0)*0.50d0)))

    !write(*,*) Mplus

	  return
	  end

!----------Mminus------------------------
    function Mminus(alpha)
	  use DATA01
	  use DATA02
    implicit none
    real(8) alpha
	  complex(8) Mminus

	  !式(3.2)参照
	  Mminus=cdsqrt(1.d0+cdexp(i*pi*nu))/(2.d0*cdexp((i*pi*(nu+1.d0))*0.25d0)&
    *((k-alpha)**((nu-1.d0)*0.50d0)))

	  return
	  end

!----------Kplus------------------------
   function Kplus(alpha)
	  use DATA01
	  use DATA02
    implicit none
    real(8) alpha
	  complex(8) Kplus

	  !式(3.40)参照
	  Kplus=i*cdsqrt(1.d0-cdexp(i*pi*nu))/(2*cdexp(i*pi*nu/4.d0)*(dsqrt(k+alpha)**nu))

	  return
	  end

!----------Kminus------------------------
   function Kminus(alpha)
	  use DATA01
	  use DATA02
    implicit none
    real(8) alpha
	  complex(8) Kminus

	  !式(3.40)参照
	  Kminus=i*cdsqrt(1.d0-cdexp(i*pi*nu))/(2*cdexp(i*pi*nu/4.d0)*(dsqrt(k-alpha)**nu))

	  return
	  end

!---------- ETA_1_nu0------------------------
   function ETA_1_nu0(alpha)
	  use DATA01
	  implicit none
    real(8) alpha
	  complex(8) ZAI_nu0,ETA_1_nu0
	  external ZAI_nu0
	  !abはkcosθ0
	  !式(4.8)参照
	  ETA_1_nu0=(ZAI_nu0(alpha)-ZAI_nu0(-ab))/(alpha+ab)

	  return
	  end

!---------- ETA_2_nu0------------------------
    function ETA_2_nu0(alpha)
	  use DATA01
	  implicit none
    real(8) alpha
	  complex(8) ZAI_nu0,ETA_2_nu0
	  external ZAI_nu0
	  !式(4.8)参照
	  ETA_2_nu0=(ZAI_nu0(alpha)-ZAI_nu0(ab))/(alpha-ab)

	  return
	  end

!---------- ETA_1_nu1------------------------
   function ETA_1_nu1(alpha)
	  use DATA01
	  implicit none
    real(8) alpha
	  complex(8) ZAI_nu1,ETA_1_nu1
	  external ZAI_nu1
	  !式(4.51)参照
	  ETA_1_nu1=(ZAI_nu1(alpha)-ZAI_nu1(-ab))/(alpha+ab)

	  return
	  end

!---------- ETA_2_nu1------------------------
    function ETA_2_nu1(alpha)
	  use DATA01
	  implicit none
    real(8) alpha
	  complex(8) ZAI_nu1,ETA_2_nu1
	  external ZAI_nu1
	  !式(4.51)参照
	  ETA_2_nu1=(ZAI_nu1(alpha)-ZAI_nu1(ab))/(alpha-ab)

	  return
	  end

!---------- ZAI_nu0 ------------------------
   function ZAI_nu0(alpha)
	  use DATA01
	  use DATA02
	  implicit none
    real(8) alpha
	  complex(8) ZAI_nu0,DE_nu0,z
	  external DE_nu0

	  z=-2.d0*i*a*(alpha+k)

	  !式(4.15)参照
	  ZAI_nu0=((1.d0+cdexp(i*pi*nu)**0.5)*cdexp(i*ka*2.d0)&
    *cdexp(i*pi*(nu-1.d0))*DE_nu0(z))/(cdexp(i*pi*nu)*pi*(dsqrt(a*2.d0)**(nu-1.d0)))

	  return
	  end

!---------- ZAI_nu1 ------------------------
   function ZAI_nu1(alpha)
	  use DATA01
	  use DATA02
	  implicit none
    real(8) alpha
	  complex(8) ZAI_nu1,DE_nu1,z
	  external DE_nu1

	  z=-2.d0*i*a*(k+alpha)

	  !式(4.39)参照
	  ZAI_nu1=((1.d0-cdexp(i*pi*nu)**0.5)*cdexp(i*pi*nu)*cdexp(i*ka*2.d0)*DE_nu1(z))/(pi*(dsqrt(a*2.d0)**nu))

	  return
	  end

!ここからはDE公式のプログラムを用いている
!---------- DE_nu0(z)------------------------
	 function DE_nu0(z)
	 use DATA01
	 use DATA04
	 implicit none
	 integer j,num
	 real(8) h,xjh,dxjh,jj
	 complex(8) DE_nu0,z,total

	 num=226
	 h=4.d-2
	 total=(0.d0,0.d0)

	 !式(4.15)参照
	 do j=-num/2,num/2,1
	  jj=dfloat(j)
	  xjh=dexp(pi*dsinh(jj*h)/2.d0)
	  dxjh=(pi*dcosh(jj*h)/2.d0)*dexp(pi*dsinh(jj*h)/2.d0)
	  total=total+((((xjh)**x_shisu_nu0)*dexp(-xjh)*dxjh)/((xjh+z)**z_shisu))
	 end do

	 DE_nu0=total*h

	 return
	 end

!---------- DE_nu1(z)------------------------
	 function DE_nu1(z)
	 use DATA01
	 use DATA04
	 implicit none
   integer j,num
	 real(8) h,xjh,dxjh,jj
	 complex(8) DE_nu1,z,total

	 num=226
	 h=4.d-2
	 total=(0.d0,0.d0)

	 !式(4.39)参照
	 do j=-num/2,num/2,1
	  jj=dfloat(j)
	  xjh=dexp(pi*dsinh(jj*h)/2.d0)
	  dxjh=(pi*dcosh(jj*h)/2.d0)*dexp(pi*dsinh(jj*h)/2.d0)
	  total=total+((((xjh)**x_shisu_nu1)*dexp(-xjh)*dxjh)/((xjh+z)**z_shisu))
	 end do

	 DE_nu1=total*h

   return
   end
