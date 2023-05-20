program main
  implicit none
  integer N,M
  parameter(N=3,M=2)
  real a(N,N),b(N,M)
  integer i,j
  integer eflag
  

!元の連立方程式
  write(*,*) 'A='
  do i=1,N
    write(*,*) (a(i,j),j=1,N)
  end do
  write(*,*) 'B='
  do i=1,N
    write(*,*) (b(i,j),j=1,N)
  end do

  eflag=0
  
  call gauss_elim(a,N,b,M,eflag)
  
  stop
end program main

subroutine gauss_elim(a,n,b,m,eflag)
implicit none
  integer n,m
  real a(n,n),b(n,m)
  integer eflag
  
  real TINY
  parameter(TINY=1.0e-20)
  
  integer i,j,k
  real r

!Gaussian elimination
do j=1,n-1
  call ppivot(a,n,b,m,j)
  if(ABS(a(i,j))<TINY)then
    eflag=1
    return
  end if
  
!forward elimination
 do i=j+1,n
   r=a(i,j)/a(j,j)
   do k=j+1,n
     a(i,k)=a(i,k)-a(j,k)*r
   end do
   do k=1,m
     b(i,k)=b(i,k)-b(j,k)*r
   end do
  end do
end do

 if(ABS(a(n,n))<TINY)then
    eflag=1
    return
  end if

!Back substitution
  do j=1,m
    do i=n,1,-1
     r=b(i,j)
     do k=i+1,n
       r=r-a(i,k)*b(k,j)
     end do
        b(i,j)=r/a(i,i)
     end do
   end do
   return 
end subroutine gauss_elim

subroutine ppivot(a,n,b,m,j)
  implicit none
  integer n,m,j
  real a(n,n),b(n,m)
  
  integer i,k,imax
  real amax,r
!find the maximum element
  amax=ABS(a(j,j))
  imax=j
  do i=j+1,n
    r=ABS(a(i,j))
    if(r>amax)then
      amax=r
      imax=i
    end if
  end do
!exchange rows if necessary
if(j==imax)return
 do k=1,n
  r=a(imax,k)
  a(imax,k)=a(j,k)
  a(j,k)=r
 end do
 do k=1,m
  r=b(imax,k)
  b(imax,k)=b(j,k)
  b(j,k)=r
 end do
 return
end subroutine ppivot