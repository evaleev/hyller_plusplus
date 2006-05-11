# site-specific variables loaded first
PSIROOT = /Users/evaleev/Development/QuantumChemistry/psi3/recent/ppc-osx-gcc34

# CODE = $(shell basename `pwd`)
CODE = hyller++
vpath %.a $(PSIROOT)/lib

#defines for AIX go here
CFLAGS = -g
CXXFLAGS = -g

LIBS = -lPSI_ciomr -lPSI_ipv1 -lm
LIBSPATH = -L$(PSIROOT)/lib
INCLUDES = -I$(PSIROOT)/include -I.

CXXSRC = main.cc matrix.cc polynom.cc misc.cc hylleraas.cc orbital.cc \
projector.cc determinant.cc csf.cc except.cc
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

