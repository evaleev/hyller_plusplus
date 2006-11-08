# site-specific variables loaded first
PSIROOT = /home/evaleev/Development/psi3/recent/x86-linux-gcc41
GSLPATH = /usr/local

# CODE = $(shell basename `pwd`)
CODE = hyller++
vpath %.a $(PSIROOT)/lib

#defines for AIX go here
CFLAGS = -g
CXXFLAGS = -g

LIBS = -lPSI_ciomr -lPSI_ipv1 -lgsl -lgslcblas -llapack -lf77blas -latlas -lm
FLIBS = -L/usr/lib/gcc/i386-redhat-linux/4.1.1 -L/usr/lib/gcc/i386-redhat-linux/4.1.1/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
LIBSPATH = -L$(PSIROOT)/lib -L$(GSLPATH)/lib -L/usr/local/LAPACK -L/usr/local/ATLAS/lib/Linux_IntelCore
INCLUDES = -I$(PSIROOT)/include -I$(GSLPATH)/include/gsl -I.

CXXSRC = main.cc matrix.cc polynom.cc misc.cc hylleraas.cc slaterhylleraas.cc orbital.cc \
projector.cc determinant.cc csf.cc except.cc gsh_basissets.cc clock.cc fock.cc
CXXOBJ = $(CXXSRC:%.cc=%.o)

$(CODE): $(CXXOBJ) 
	$(CXX) $^ $(LDFLAGS) $(LIBSPATH) $(LIBS) $(FLIBS) -o $(CODE)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	/bin/rm -f $(CODE) *.o

realclean::
	/bin/rm -f $(CODE) *.o core

depend:
	makedepend -- $(INCLUDES)  -- $(CSRC) -- $(CXXSRC)

