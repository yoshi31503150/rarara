field=0.d0
      do n=-1799,1799,1
        the1=dble(n)*0.1d0
        theta=the1*(pi/180.d0)
        ab1=dcmplx(-ka* cos(theta),0.d0)
        ab11=dcmplx(-kk* cos(theta),0.d0)
  !---
        if(n .gt. -1800 .and. n .lt.(-1800+int(10.d0*the0))) then
          sigma=abs(( ((exp(-I*ab1)*Uminus(ab11))+(exp(I*ab1)         &
            *Uplus(ab11)))*(-0.5d0*i)                                 &
            -0.5d0*((exp(-I*ab1)*Vminus(ab11))                        &
            +(exp(I*ab1)*Vplus(ab11)))*kk*sin(abs(theta))))
          db=sigma
        else if(n.eq.(-1800+int(10.d0*the0))) then
          sigma=0.d0
        else if(n .gt. (-1800+int(10.d0*the0)) .and. n .lt. 0) then
          sigma=abs(( ((exp(-I*ab1)*Uminus(ab11))+(exp(I*ab1)         &
            *Uplus(ab11)))*(-0.5d0*i)                                 &
            -0.5d0*((exp(-I*ab1)*Vminus(ab11))                        &
            +(exp(I*ab1)*Vplus(ab11)))*kk*sin(abs(theta))))
          db=sigma
!            db=20.d0*dlog10(sigma)
        else if(n.eq.0) then
          sigma=0.d0
        else if(n .gt.0 .and. n .lt. (1800-int(10.d0*the0))) then
          sigma=abs(( ((exp(-I*ab1)*Uminus(ab11))+(exp(I*ab1)         &
            *Uplus(ab11)))*(-0.5d0*i)                                 &
            +0.5d0*((exp(-I*ab1)*Vminus(ab11))                        &
            +(exp(I*ab1)*Vplus(ab11)))*kk*sin(theta)))
          db=sigma
        else if(n .eq. (1800-int(10.d0*the0))) then
          sigma=0.d0
        else
          sigma=abs(( ((exp(-I*ab1)*Uminus(ab11))+(exp(I*ab1)         &
            *Uplus(ab11)))*(-0.5d0*i)                                 &
            +0.5d0*((exp(-I*ab1)*Vminus(ab11))                        &
            +(exp(I*ab1)*Vplus(ab11)))*kk*sin(theta)))
          db=sigma
        end if
  !-----------------------------------
        write(*,*) n, db
        field(n)=db
!          write(6,*) db
!----------------------------------
      end do
      write(*,*)'finished'
      write(*,*)'db?(Y=else,N=0)'
      read(*,*) YYY1
      if(YYY1 == 0)then
        do n=-1799,1799,1
          write(10,*) field(n)
        end do
      else
        mmax=maxval(field)
        write(*,*) 'MMAX=',mmax
        do n=-1799,1799,1
          write(10,*) 20.d0*log10(field(n)/mmax)
        end do
      endif
      CLOSE(10)
      write(*,*)'2a, b', ipka, ipkb
      write(*,*)'epr, mur', epsilonr, mur
      STOP
      END PROGRAM CHECK