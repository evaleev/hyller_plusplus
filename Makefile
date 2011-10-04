#
# The MPQC section
#
SCCONFIG = /Users/evaleev/Development/workspace/install/mpqc-gcc/bin/sc-config
CXX := $(shell $(SCCONFIG) --cxx)
CXXFLAGS := $(shell $(SCCONFIG) --cxxflags)
CPPFLAGS := $(shell $(SCCONFIG) --cppflags)
CPPFLAGS += -I.
MPQCLIBS := $(shell $(SCCONFIG) --libs)
MPQCLIBDIR  := $(shell $(SCCONFIG) --libdir)
LTLINK := $(shell $(SCCONFIG) --ltlink)
LTLINKBINOPTS := $(shell $(SCCONFIG) --ltlinkbinopts)
# The suffix generated by the -M compiler option
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M
CXXDEPEND = $(CXX)
OBJSUF = o

# site-specific variables loaded first
GSLPATH = /usr/local
BOOSTPATH = /Users/evaleev/Development/boost/install
BOOSTINCLUDE = $(BOOSTPATH)/include
BOOSTLIB = $(BOOSTPATH)/lib

CODE = hyller++
vpath %.a $(PSIROOT)/lib

#BLAS =  -L/usr/local/LAPACK -L/usr/local/ATLAS/lib/Linux_IntelCore -llapack -lf77blas -latlas 
BLAS = /System/Library/Frameworks/vecLib.framework/vecLib
#LIBS = -lPSI_ciomr -lPSI_ipv1 -lgsl -lgslcblas $(BLAS) -lm
LIBS = $(MPQCLIBS) -lPSI_ciomr -lPSI_ipv1 $(GSLPATH)/lib/libgsl.dylib $(BLAS) -lm
#FLIBS = -L/usr/lib/gcc/i386-redhat-linux/4.1.1 -L/usr/lib/gcc/i386-redhat-linux/4.1.1/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
FLIBS =
LIBDIR = $(MPQCLIBDIR)
INCLUDES = -I$(GSLPATH)/include/gsl -I$(BOOSTINCLUDE) -I.
CPPFLAGS += $(INCLUDES)

CXXSRC = main.cc matrix.cc polynom.cc misc.cc hylleraas.cc slaterhylleraas.cc orbital.cc \
projector.cc determinant.cc csf.cc except.cc gsh_basissets.cc clock.cc fock.cc
CXXOBJ = $(CXXSRC:%.cc=%.o)

$(CODE): $(CXXOBJ)
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS) $(FLIBS) $(LTLINKBINOPTS)
#	$(CXX) $^ $(LDFLAGS) $(LIBSPATH) $(LIBS) $(FLIBS) -o $(CODE)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c $<

clean:
	/bin/rm -f $(CODE) *.o *.d

realclean:: clean
	/bin/rm -f core

depend:
	makedepend -- $(INCLUDES)  -- $(CSRC) -- $(CXXSRC)

#
# dependencies
#

ifneq ($(CXXDEPENDSUF),none)
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< > /dev/null
	sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)
	/bin/rm -f $(*F).$(CXXDEPENDSUF)
else
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
endif

ifneq ($(DODEPEND),no)
include $(CXXOBJ:%.o=%.d)
endif
