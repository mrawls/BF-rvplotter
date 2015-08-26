c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..
c..
c..
      subroutine convlv(data,n,respns,m,isign,ans)
c..
c..
c.......................................................................
c..
C      implicit none
      include             'convspec.h'
      integer             i,n,m,isign,no2
      double precision    data(n),respns(n)
      double complex      fft(maxcon),ans(n)
c..
c..
      do 11 i=1,(m-1)/2
       respns(n+1-i) = respns(m+1-i)
 11   continue
      do 12 i=(m+3)/2,n-(m-1)/2
       respns(i) = 0.0d0
 12   continue
c..
      call twofft(data,respns,fft,ans,n)
c..
      no2 = n/2
      do 13 i=1,no2+1
       if (isign.eq.1) then
        ans(i) = fft(i)*ans(i)/no2
       else if (isign.eq.-1) then
        if (abs(ans(i)).eq.0.0) pause 'death in convlv'
         ans(i) = fft(i)/ans(i)/no2
       else
        pause 'no meaning for isign'
       endif
 13   continue
      ans(1) = dcmplx(dble(ans(1)),dble(ans(no2+1)))
      call realft(ans,no2,-1)
c..
      return
      end
c..
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..
c..
c..
      subroutine twofft(data1,data2,fft1,fft2,n)
c..
c..
c.......................................................................
c..
      implicit none
      integer             j,n,n2
      double precision    data1(n),data2(n)
      double complex      fft1(n),fft2(n),h1,h2,c1,c2
c..
c..
      c1 = dcmplx(0.5,0.0)
      c2 = dcmplx(0.0,-0.5)
      do 11 j=1,n
        fft1(j) = dcmplx(data1(j),data2(j))
 11   continue
c..
      call four1(fft1,n,1)
c..
      fft1(n+1) = fft1(1)
      n2        = n+2
      do 12 j=1,n/2+1
       h1         = c1 * (fft1(j) + conjg(fft1(n2-j)))
       h2         = c2 * (fft1(j) - conjg(fft1(n2-j)))
       fft1(j)    = h1
       fft1(n2-j) = conjg(h1)
       fft2(j)    = h2
       fft2(n2-j) = conjg(h2)
 12   continue
c..
      return
      end
c..
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..
c..
c..
      subroutine realft(data,n,isign)
c..
c..
c.......................................................................
c..
      implicit none
      integer             n,isign,i,n2p3,i1,i2,i3,i4
      double precision    wr,wi,wpr,wpi,wtemp,wrs,wis,
     @                    h1r,h1i,h2r,h2i,theta,c1,c2
      double precision    data(2*n)
c..
c..
      theta = 6.28318530717959d0/2.0d0/dble(n)
      wr = 1.0d0
      wi = 0.0d0
      c1 = 0.5
c..
      if (isign.eq.1) then
       c2 = -0.5
       call four1(data,n,+1)
       data(2*n+1) = data(1)
       data(2*n+2) = data(2)
      else
       c2          = 0.5
       theta       = -theta
       data(2*n+1) = data(2)
       data(2*n+2) = 0.0
       data(2)     = 0.0
      endif
c..
      wpr  = -2.0d0 * (sin(0.5d0*theta)**2)
      wpi  = sin(theta)
      n2p3 = 2*n+3
      do 11 i=1,n/2+1
       i1       = 2*i-1
       i2       = i1+1
       i3       = n2p3-i2
       i4       = i3+1
       wrs      = dble(wr)
       wis      = dble(wi)
       h1r      = c1*(data(i1)+data(i3))
       h1i      = c1*(data(i2)-data(i4))
       h2r      = -c2*(data(i2)+data(i4))
       h2i      = c2*(data(i1)-data(i3))
       data(i1) = h1r+wrs*h2r-wis*h2i
       data(i2) = h1i+wrs*h2i+wis*h2r
       data(i3) = h1r-wrs*h2r+wis*h2i
       data(i4) = -h1i+wrs*h2i+wis*h2r
       wtemp    = wr
       wr       = wr*wpr - wi*wpi+wr
       wi       = wi*wpr + wtemp*wpi+wi
 11   continue
c..
      if (isign.eq.1) then
       data(2) = data(2*n+1)
      else
       call four1(data,n,-1)
      endif
c..
      return
      end
c..
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..
c..
c..
      subroutine four1(data,nn,isign)
c..
c..compute the fft
c..isign = 1 foward transform ; isign = -1 inverse transform
c.......................................................................
c..
      implicit none
      integer              nn,i,j,n,isign,m,mmax,istep
      double precision     wr,wi,wpr,wpi,wtemp,theta,tempr,tempi
      double precision     data(2*nn)
c..
c..
      n = 2*nn
      j = 1
      do 11 i=1,n,2
       if (j.gt.i) then
        tempr     = data(j)
        tempi     = data(j+1)
        data(j)   = data(i)
        data(j+1) = data(i+1)
        data(i)   = tempr
        data(i+1) = tempi
       endif
       m = n/2
 1     if ((m.ge.2).and.(j.gt.m)) then
        j = j - m
        m = m/2
        goto 1
       endif
       j = j + m
 11   continue
c..
      mmax = 2
 2    if (n.gt.mmax) then
       istep = 2*mmax
       theta = 6.28318530717959d0/(isign*mmax)
       wpr   = -2.d0*sin(0.5d0*theta)**2
       wpi   = sin(theta)
       wr    = 1.d0
       wi    = 0.d0
       do 13 m=1,mmax,2
        do 12 i=m,n,istep
         j=i+mmax
         tempr     = dble(wr)*data(j)   - dble(wi)*data(j+1)
         tempi     = dble(wr)*data(j+1) + dble(wi)*data(j)
         data(j)   = data(i)-tempr
         data(j+1) = data(i+1)-tempi
         data(i)   = data(i)+tempr
         data(i+1) = data(i+1)+tempi
 12     continue
        wtemp = wr
        wr    = wr*wpr - wi*wpi + wr
        wi    = wi*wpr + wtemp*wpi + wi
 13    continue
       mmax = istep
       go to 2
      endif
c..
      return
      end
c.......................................................................
c..eof
