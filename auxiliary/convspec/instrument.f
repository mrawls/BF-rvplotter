c.......................................................................
c
c  routine instrument  stuffs the instrumental profile array; this
c                      routine called once
c  routine phi         computes the relative instrumental response 
c                      (assumes a Gaussian)
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..
c..
c..
      subroutine          instrument(m,wcen,flag)
c..
c..
c..given the instrumental profile sigma in velocity units
c..
c..these profiles are loaded in wrap-around order for the convolution
c..zero spatial information is the first index
c
c  flag = 0 ; call for communicating instrumetal parameters
c  flag = 1 ; call for setting up convolution
c
c.......................................................................
      include             'convspec.h'
      include             'const.dek'
      integer             i,m,iplus,iminus,np2,flag
      parameter           (np2 = 17)
      double precision    dv,xdv,norm,phi,pixpres,wcen
      integer             pwrsof2(np2)
c
c
      data pwrsof2 /16,32,64,128,256,512,1024,2048,4096,8192,
     &              16384,32768,65536,131072,262144,524288,
     &              1048576 /
c
c
c
c
c  compute the instrumental resolution in velocity units [km/s]
c  sigma in km/s given by FWHM/2.35

       profile = ckms/(2.35*R_fac)
c
c
c  dv is the velocity sampling of the pixels km/s/pixel...
c  the number of pixels per resolution element =  profile/dv
c

      dv = dwave/wcen * ckms
      pixpres = 2.35 * profile/dv
      hdv = dv/resfac
      nresponse  = 2*int(conwindo*profile/hdv) + 1
c
c
c  now stuff the response function in wrap around order
c
      response(1) = phi(0.0d0,profile)
      norm        = response(1)
      do 11 i=1,int(nresponse/2)
       xdv              = real(i)*hdv 
       iplus            = i + 1
       iminus           = nresponse - (i-1)
       response(iplus)  = phi(xdv,profile)
       response(iminus) = response(iplus)
       norm = norm + response(iplus) + response(iminus)
 11   continue
c  
c  for the convolution integral, the response function
c  needs to be normalized or flux is not conserved...
c  unit width is used because data are discretized by
c  pixel indices in the convolution
      do 13 i=1,nresponse
       response(i) = response(i)/norm
 13   continue
c  
c
c  compute the length of the convolution functions..
c  NCONDAT and NFFT; do a table search on the powers 
c  of 2; this ensures the smallest padding of the
c  CONVDATA array in routine convolve
c
      ncondat = int(resfac)*(m-1) - 1
      nfft   = ncondat + (nresponse-1)/2 + 1 
      do 15 i=1,np2
       if (nfft.le.pwrsof2(i)) then
        nfft = pwrsof2(i)
        GOTO 17
       end if
 15   continue

c
c  this is seen only if we exceed the powers of 2 table
      write(STDOUT,*) ' NRESPONSE = ',nresponse
      write(STDOUT,*) ' M         = ',m
      write(STDOUT,*) ' NCONDAT   = ',ncondat
      write(STDOUT,*) ' NFFT      = ',nfft
      write(STDOUT,*) ' NFFT MAX  = ',maxcon
      stop ' ERROR:(instrument): NFFT not defined? too big?'
c
 17   if (nfft.gt.maxcon) then
       stop ' ERROR:(instrument): nfft > maxcon'
      end if
c
c
      if (flag.eq.1) then 
        write(STDOUT,'(a,f8.0)') 
     @  ' Spectrograph R     [lam/dlam] =',R_fac
        write(STDOUT,'(a,f8.3)') 
     @  ' Instrumental Sigma     [km/s] =',profile
       write(STDOUT,'(a,f8.3)') 
     @  ' Pixel Resolution       [km/s] =',dv
       write(STDOUT,'(a,f8.3)') 
     @ ' Pix / Res Elem         [FWHM] =',pixpres
       write(STDOUT,'(a,f8.3)') 
     @ ' Response Widths       [sigma] =',conwindo
       write(STDOUT,'(a,f8.3)') 
     @ ' Convolution Res       [1/pix] =',resfac
       write(STDOUT,'(a,i5)')   
     @ ' Response Length         [pix] =',nresponse
       write(STDOUT,'(a,i7)')   
     @ ' Convolution Pixels      [pix] =',nfft
      end if
c
c  
      return
c
      end
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  
c  
c  
      double precision function phi(dv,width)
c  
c  
c    given the sigma in velocity units this routine computes the 
c    relative Gaussian value of the instrumental profile.  
c    called by routine instrument iteratively
c.......................................................................
c  
C      implicit none
      include             'convspec.h'
      double precision     dv,z,width
c
c   dv is value at which the instrumental response is to be evaluated
c
      z   = dv/width
      phi = dexp(-(z**2)/2.0d0)
c
      return
      end
c
c  eof
c.......................................................................
