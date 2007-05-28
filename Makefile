#
# The MPQC section
#
SCCONFIG = /Users/evaleev/Development/mpqc/x86-linux-gcc40/bin/sc-config
CXX := $(shell $(SCCONFIG) --cxx)
CXXFLAGS := $(shell $(SCCONFIG) --cxxflags)
CPPFLAGS := $(shell $(SCCONFIG) --cppflags)
MPQCLIBS := $(shell $(SCCONFIG) --libs)
MPQCLIBDIR  := $(shell $(SCCONFIG) --libdir)
LTLINK := $(shell $(SCCONFIG) --ltlink)
LTLINKBINOPTS := $(shell $(SCCONFIG) --ltlinkbinopts)

# site-specific variables loaded first
#PSIROOT = /home/evaleev/Development/psi3/recent/x86-linux-gcc41
PSIROOT = /Users/evaleev/Development/Psi3/x86-linux-gcc40
#GSLPATH = /usr/local
GSLPATH = /opt/local
BOOSTPATH = /opt/local
BOOSTINCLUDE = $(BOOSTPATH)/include
BOOSTLIB = $(BOOSTPATH)/lib

CODE = hyller++
vpath %.a $(PSIROOT)/lib

#BLAS =  -L/usr/local/LAPACK -L/usr/local/ATLAS/lib/Linux_IntelCore -llapack -lf77blas -latlas 
BLAS = /System/Library/Frameworks/vecLib.framework/vecLib
#LIBS = -lPSI_ciomr -lPSI_ipv1 -lgsl -lgslcblas $(BLAS) -lm
LIBS = $(MPQCLIBS) -lPSI_ciomr -lPSI_ipv1 -lgsl $(BLAS) -lm
#FLIBS = -L/usr/lib/gcc/i386-redhat-linux/4.1.1 -L/usr/lib/gcc/i386-redhat-linux/4.1.1/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
FLIBS =
LIBDIR = $(MPQCLIBDIR) -L$(PSIROOT)/lib -L$(GSLPATH)/lib
INCLUDES = -I$(PSIROOT)/include -I$(GSLPATH)/include/gsl -I$(BOOSTINCLUDE) -I.

CXXSRC = main.cc matrix.cc polynom.cc misc.cc hylleraas.cc slaterhylleraas.cc orbital.cc \
projector.cc determinant.cc csf.cc except.cc gsh_basissets.cc clock.cc fock.cc integrate.cc \
product_ansatz.cc
CXXOBJ = $(CXXSRC:%.cc=%.o)

$(CODE): $(CXXOBJ) 
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS) $(FLIBS) $(LTLINKBINOPTS)
#	$(CXX) $^ $(LDFLAGS) $(LIBSPATH) $(LIBS) $(FLIBS) -o $(CODE)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c $<

clean:
	/bin/rm -f $(CODE) *.o

realclean::
	/bin/rm -f $(CODE) *.o core

depend:
	makedepend -- $(INCLUDES)  -- $(CSRC) -- $(CXXSRC)

