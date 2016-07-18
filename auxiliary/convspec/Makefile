### I HAVE NO IDEA HOW TO MAKEFILE
# typing 'make' would make an executable but it would crash as it finished running
# here is what I did to get it to work on my Mac OS X 10.10.4
# gfortran -O -Bstatic -c fft.f convolve.f getparms.f convspec.f instrument.f voigt.f
# gfortran -O -Bstatic convspec.f -o convspec
###
FC=gfortran
FFLAGS=-O -Bstatic
TLIB=./tool.a
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; $(FC) $(FFLAGS) -c $*.F
#
default:: convspec
#

convspec: fft.f convolve.f getparms.f convspec.f instrument.f voigt.f convspec.h 
	$(FC) -w convspec.f -o convspec $(TLIB) 
#
