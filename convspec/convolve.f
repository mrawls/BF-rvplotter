c..............................................................................
c..............................................................................
c  
c  
c  
      SUBROUTINE          convolve
c  
c  
c  this routine sets up the call to the Num Recipes convolver
c  convolution resolution is higher than data resolution so 
c  that it data are "smoother"
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include             'convspec.h'
      integer             i,k,nmax2,ndat,pixi,mfft
      parameter           (nmax2 = maxcon*2)
      double precision    xa(maxvec),ya(maxvec),x,y,y2a(maxvec)
      double precision    ans(nmax2+2),phiov(maxcon)
c
c
c
c
       ndat    = ndata
       mfft    = nfft
      
c
c
c  stuff the "smoothing" data array, one only need keep the index
c  as the abscissa; then obtain the spline coefficients
c
        do 27 pixi=1,ndata
         xa(pixi)  = real(pixi)
         ya(pixi)  = wrkflx(pixi)
         y2a(pixi) = 0.0d0
 27     continue
      call spline(xa,ya,ndat,y2a)
c  
c
c  stuff the end points and the interpolation points
c
      convdata(1)       = ya(1)
      convdata(ncondat) = ya(ndat)
      do 29 i=2,ncondat-1
       x = 1.0d0 + real(i-1)/resfac
       call splint(xa,ya,y2a,ndat,x,y)
       convdata(i) = y
 29   continue
c
c  pad the data array with unity, the continuum flux level
c  this assumes that no features are on the edges of the 
c  data array
c
      do 25 i=ncondat+1,nfft
       convdata(i) = 1.0d0
 25   continue
c
c
c  the FFT's do violence to the response function-
c  so re-load it into phi(v)
c
      do 26 i=1,nresponse
        phiov(i) = response(i)
 26   continue
c
c
c  we are "go" for the FFT's
c
      CALL convlv(convdata,mfft,phiov,nresponse,+1,ans)
c
c
c
c  now, restuff the wrkflx array and bail; pick off every resfac
c  element and stuff into idx element of wrkflx
c
        do 33 pixi=1,ndata
         k = 1 + (pixi-1)*int(resfac)
         wrkflx(pixi) = ans(k)
 33     continue
c
      return
      end
c
c
c..............................................................................
