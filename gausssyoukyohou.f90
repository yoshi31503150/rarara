subroutine gauss(n,a,b)
implicit none
!Gaussの消去法によりn元連立方程式
!Ax=bの解を求める
!入力
!n:次元
!a:係数行列Aの要素（n x n）
!b:ベクトルbの要素(n)
!出力
!a:行列Aを上三角行列に変換したもの
!（左下の非対角要素は意味がない）
!b:解　x
!input/output
  integer n
  real a(n,n),b(n)
!local
  integer i,j,k
!begin:
  do i=1,n
    do j=i+1,n
     a(j,i)=a(j,i)/a(i,i)        !Pivot:c=a_ji/a_ii
     do k=i+1,n
       a(j,k)=a(j,k)-a(j,i)*a(i,k)
      end do
        b(j)=b(j)-a(j,i)*b(i)
     end do
   end do
   
   do i=n ,1,-1
     do j=i+1,n
       b(i)=b(i)-a(i,j)+b(j)
      end do
      b(i)=b(i)/a(i,i)           !Ax=bの解x
    end do
    return 
    end