c
c
      implicit none
c
c  integer parameters
c
c  maxspecies --- max number of species to be fit
c  maxions    --- max number of ion transitions to be fit
c  maxpix     --- max number of pixels in a given ion transition
c  maxlines   --- max number of Voigt profiles to be fit
c  maxcoeffs  --- max number of coefficients to be fit
c
c
c
c
      integer            maxpix,maxlines,maxvec
      parameter          (maxlines =  4000       ,  
     @                    maxpix   =  1048576    ,
     @                    maxvec   =  5.0*maxpix )

      integer            SCREEN,STDIN,STDOUT
      parameter          (SCREEN = 6, STDIN = 5, STDOUT = SCREEN)
c
c
c
c  atomic constants and general acounting
c
      integer            lines,ndata
      double precision   con1,con2,lambda0,wave_max,wave_min,dwave,zabs
      character*80       ion_name
c
c
      COMMON/atomblck/con1(maxlines),con2(maxlines),lambda0(maxlines),
     &                wave_max,wave_min,dwave,zabs
      COMMON/atomchar/ion_name(maxlines)
      COMMON/accntint/lines,ndata
c
c
c  data arrays and parameter arrays
c
      double precision   lambda,wrkflx
      double precision   zline,nline,bline
c
      COMMON/fitblck/lambda(maxpix),wrkflx(maxpix)
      COMMON/parblck/zline(maxlines),nline(maxlines),bline(maxlines)
c
c
c  convolution arrays 
c  ***** IMPORTANT NOTE FOR CONVOLUTION *****
c  NOTE:  there is a DATA statement in instruments.f that contains
c  a list of powers of 2.  It goes up to MAXCON.  
c  Parameter MAXCON must be a power of 2.
c
      logical           convolving
      integer           maxcon
      parameter         (MAXCON = 1048576)
      integer           nresponse,ncondat,nfft
      double precision  response,convdata
      double precision  R_fac,slit,conwindo,resfac,profile,hdv
c
      COMMON/convlblck/convolving
      COMMON/convrblck/response(maxpix),convdata(maxcon)
      COMMON/conviblck/nresponse,ncondat,nfft
      COMMON/instblck/R_fac,slit,conwindo,resfac,profile,hdv
c
c
c............................................... end minfit.h .................
