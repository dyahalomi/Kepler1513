NESTLIBDIR = /Users/Danielyahalomi/Documents/Columbia/Research/Kipping/MultiNest/MultiNest_v2.17/
LIBS = -L/Users/dkipping/Storage1/Work/Documents/Transit_Work/CODES/swiftlib16/ -lswift -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
LIBS = -L/Users/Danielyahalomi/Documents/Columbia/Research/Kipping/SWIFT/ -lswift -L$(NESTLIBDIR) -lnest3
FFLAGS = -O2
FC=gfortran

OBJFILES = params.o ntrans_rmvs3.o transit.o like.o nestwrap.o main.o

all: zwift

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90 
 
zwift: $(OBJFILES)
	$(FC) -o ../zwift $(OBJFILES) \
	$(FFLAGS) $(LIBS) -llapack  

clean:
	rm -f *.o *.mod ../zwift
