c------------------------------------------------------------------------------
c
c
c
	SUBROUTINE        getparams(paramlist)
c
c
c
c    R_FAC is a double precision real input variable.  R_FAC is the 
c    instrument resolution in (lambda/delta lambda) units.  For example, 
c    HIRES has several R factors, one of which is R=45000., which gives 
c    a velocity resolution of 6.6 km/s/arcsecond.  The code computes the 
c    instrument resolution as a Gaussian sigma, where the FWHM is given by 
c    c/R_FAC times SLIT (see below), where c is the speed of light.  If 
c    R_FAC=0 and CONFLAG=1, then the  the instrument resolution as a Gaussian 
c    sigma is taken directly from the value of SLIT.
c
c
c    CONWINDO is a double precision real input variable.  CONWINDO is the 
c    width in units of Gaussian sigma of the instrumental profile (response) 
c    over which the convolution is performed (to be thought of as the half 
c    length of the response function in sigma units).  Setting CWIDTH to 3.0 
c    is usually sufficient.
c
c
c    RESFAC is a double precision real input variable.  RESFAC is the 
c    sampling rate for the convolution.  Since instruments usually have a 
c    near critical sampling rate, convolution may result in artifacts such 
c    as flux non-conservation.  RESFAC is used to increase the "sampling 
c    rate" by adding pixels in the convolution integrals only.  The values 
c    of the WRKFLX (a smooth noiseless model function) in these added pixels 
c    are obtained by interpolation of a cubic spline.  RESFAC should be an 
c    odd whole number.  Setting RESFAC to 3.0 is usually sufficient.
c  
c
c
C***************** input parameters for specsynth *****************************
C 45000.  R_fac - the spectrograph resolution R=lambda/(Delta lambda)
C 3.0     conwindo -  # of inst. profile sigma to include in convolution 
C 3       resfac - (odd integer) sample rate factor for convolution arrays
C******************************************************************************
c
	include           'convspec.h'
        integer           conflag
        character*80      mf_path,par_file,paramlist
c
c
c
c  check current directory for parameter file, if DNE goto error trap, which
c  grabs the default file which lives where the user has set their MFDIR
c  environment variable
c

        convolving = .false.

c
c  did we enter the parameter file on the command line, or
c  if not, use the default
c

        IF (paramlist.eq.' ') then 
          STOP 'no input parameter file!'
        END IF

	OPEN(unit=3,file=paramlist,err=999,status='old')

 01	read(3,*)                 ! by pass the header
        read(3,*) R_fac
        read(3,*) conwindo
        read(3,*) resfac        ! must be odd integer
c
	CLOSE(unit=3)

        convolving = .true.

	WRITE(6,*) 'SPECSYNTH: input parameters are:'                

	RETURN
c
c
c
c
c

c  error opening local parameter file; check for the default settings
c  get the Unix environment variable MFDIR, which contains the path
c  to where the file lives;  try to open it; on fail, abort 

 999	write(STDOUT,*) ' cannot find/open parameter file'
        STOP 
c
	end
c

