      subroutine spline(x,y,n,y2)   

      implicit none
c.. 
c..given arrays x y of length n containing a tabulated function i.e y=f(x), 
c..with the x monotonically increasing and given values for yp1 and ypn, the
c..first derivative at the pints 1 and n respectively, this routine returns 
c..the array y2 of length n which contains the second derivatives of the    
c..at the tabulated points x.   
c..if yp1 and/or ypn are set larger than 1e30, the routine sets the boundary
c..condtion for a natural spline (one with zero second derivative) at the   
c..boundary.
c..this routine is only called once for any given x and y.  
c.. 
c..declare  
      integer            n,i,k,nmax 
      parameter          (nmax=10000) 
      double precision   x(nmax),y(nmax),y2(nmax),yp1,ypn,
     &                   u(nmax),sig,p,qn,un
      parameter          (yp1=1.0e33, ypn=1.0e33)
c.. 
c.. 
c..the lower boundary condition is set to be either natural 
      if (yp1 .gt. .99E30) then 
        y2(1)=0.0   
        u(1)=0.0
c..or it has a specified first derivative   
      else  
        y2(1)=-0.5  
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      end if
c.. 
c..this is the decomposition loop of the tridiagonal algorithm. y2 and u
c..are used for temporary storage of the decomposition factors. 
      do 11 i=2,n-1 
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))   
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     +      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p   
11    continue  
c.. 
c..the upper boundary condition is set to be either natural 
      if (ypn .gt. .99E30) then 
        qn=0.0  
        un=0.0  
c.. 
c..or it has a specified first derivative   
      else  
        qn=0.5  
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      end if
c.. 
c..this is the backsubstitution loop of the tridiagonal algorithm   
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)  
      do 12 k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue  
c.. 
c..bye bye  
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine splint(xa,ya,y2a,n,x,y)

      implicit none
c.. 
c..given the arrays xa ya of length n, which tablulate a monotonic function,
c..and given the array y2a, which is the output of spline (above), and given
c..a value of x this routine returns a cubic spline interpolated value of y.
c.. 
c..declare  

      integer             n,klo,khi,k,nmax   
      parameter          (nmax=800000) 
      double precision    xa(nmax),ya(nmax),y2a(nmax),x,y,h,a,b  
c.. 
c.. 
c..find the right place in the table by bisection. this is optimal if   
c..the sequential calls to this routine are at random values of x.  
c..if the sequential calls are in order and closely spaced, one might store 
c..the values of klo and khi and test if they remain appropriate on next call   
      klo=1 
      khi=n 
1     if (khi-klo .gt. 1) then  
       k=(khi+klo)/2
       if (xa(k) .gt. x) then   
        khi=k   
       else 
        klo=k   
       end if   
       goto 1   
      end if    
c.. 
c..klo and khi now bracket the input value of x 
c..the xa's must be distinct
      h=xa(khi)-xa(klo) 
      if (h .eq. 0.0d0) then
        write(6,99) khi,xa(khi),klo,xa(klo)
   99   format(1x,'  khi=',i4,'  xa(khi)=',1pe13.5,'  klo=',
     1         i4,'  xa(klo)=',1pe13.5)
        stop 'bad xa input in routine splint' 
      end if
c.. 
c..evaluate the cubic spline
      a=(xa(khi)-x)/h   
      b=(x-xa(klo))/h   
      y=a*ya(klo)+b*ya(khi)+
     +      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
c.. 
c..bye bye  
      return    
      end   
