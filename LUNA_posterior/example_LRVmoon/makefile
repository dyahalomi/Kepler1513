NESTLIBDIR = ..
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
FC=mpif90
OBJFILES = params.o jasmine.o mandel.o luna.o transit.o radial.o like.o nestwrap.o main.o

all: LRVmoon 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
LRVmoon: $(OBJFILES)
	$(FC) -o ../LRVmoon $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../LRVmoon
