c..............................................................................
c
c
      PROGRAM    cnv_spec
c
c
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'convspec.h' 
      integer            m,flag
      double precision   wcen
      character*80       inputspec,paramlist
c
c
c
      flag  = 1
      lines = 0

c     read in the input spectrum

c     grab the instrument resolution and convolution parameters

      CALL getarg(1,inputspec)  
      CALL getarg(2,paramlist)  

      CALL getparams(paramlist)            ! load the input deck
      CALL initspectrum(inputspec)         ! initializes the continuum

c     for the convolution approximate the sampling at the center of the
c     spectrum

      m = ndata
      wcen = 0.5*(lambda(ndata)+lambda(1)) ! this is an approximation

      CALL instrument(m,wcen,flag)         

      if (convolving) CALL convolve        ! FFT convolution with ISF

      CALL output(inputspec)               ! write output file (root name 
                                           ! is same as paramlist) ".out"
c
c  we are done
c

      STOP

      END
c
c
c..............................................................................
c  
c  
c  
      SUBROUTINE        initspectrum(input)
c  
c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'convspec.h'
      integer           pixi,maxdata
      parameter         (maxdata=1048576)
      double precision  dw1,dwn,dum
      character*80      input

      OPEN(unit=1,file=input)

      ndata = 0
      DO 11 pixi=1,maxdata
        READ(1,*,end=12) lambda(pixi),dum,wrkflx(pixi)
        ndata = ndata + 1
 11   CONTINUE

      IF (ndata.gt.maxpix) then
        WRITE(SCREEN,*) 'NDATA > MAXPIX'
        STOP
      END IF

 12   CLOSE(unit=1)

      dw1   = 0.5*(lambda(2)-lambda(1))
      dwn   = 0.5*(lambda(ndata)-lambda(ndata-1))
      dwave = 0.5*(dw1+dwn)
  
      RETURN
      END
c
c..............................................................................
c
c
      SUBROUTINE         output(input)
c
c

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'convspec.h'
      include           'const.dek'
      integer           pixi
      double precision  wave,vel,flux,noise
      character*80      out_file,input
c
c
c
c
c
      CALL fappend(input,'cnv.out',out_file)

      OPEN(unit=1,file=out_file,status='unknown')
c
       do 05 pixi=1,ndata
        wave  = lambda(pixi)
        vel   = 0.00
        flux  = wrkflx(pixi)
        noise = 1.0d0/30.0d0
        write(1,'(4f12.6)') wave,vel,flux,noise
 05    continue
c
       CLOSE(unit=1)
 03   continue

c
c
c

      WRITE(6,'(a,a)') 'OUTPUT FILE is : ',out_file

c
c
c
      return
      end
c
c  
c..............................................................................

      include 'convolve.f'
      include 'fft.f'
      include 'instrument.f'
      include 'getparms.f'
      include 'voigt.f'
      include 'spline.f'
      include 'strappend.f'

c /END/

