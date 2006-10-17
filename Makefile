# site-specific variables loaded first
PSIROOT = /home/evaleev/Development/psi3/recent/x86-linux-gcc41
GSLPATH = /usr/local

# CODE = $(shell basename `pwd`)
CODE = hyller++
vpath %.a $(PSIROOT)/lib

#defines for AIX go here
CFLAGS = -g
CXXFLAGS = -g

LIBS = -lPSI_ciomr -lPSI_ipv1 -lgsl -lgslcblas -lm
LIBSPATH = -L$(PSIROOT)/lib -L$(GSLPATH)/lib
INCLUDES = -I$(PSIROOT)/include -I$(GSLPATH)/include/gsl -I.

CXXSRC = main.cc matrix.cc polynom.cc misc.cc hylleraas.cc slaterhylleraas.cc orbital.cc \
projector.cc determinant.cc csf.cc except.cc gsh_basissets.cc clock.cc
CXXOBJ = $(CXXSRC:%.cc=%.o)

$(CODE): $(CXXOBJ) 
	$(CXX) $^ $(LDFLAGS) $(LIBSPATH) $(LIBS) -o $(CODE)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	/bin/rm -f $(CODE) *.o

realclean::
	/bin/rm -f $(CODE) *.o core

depend:
	makedepend -- $(INCLUDES)  -- $(CSRC) -- $(CXXSRC)

