# Makefile for SimulatePointSourceInFP
#
#
# 
#


ARCH := $(shell uname)

# linux flags
ifeq ($(ARCH),Linux)
CXX           = g++ 
CXXFLAGS      =  -g -O3 -Wall -fPIC -fno-strict-aliasing  -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -D_LARGEFILE64_SOURCE
LD            = g++
LDFLAGS       = -g
# LDFLAGS       =  -pg -O
SOFLAGS       = -shared
endif

# Apple OS X flags
ifeq ($(ARCH),Darwin)
CXX           = g++ 
CXXFLAGS      = -g -O3 -Wall -fPIC  -fno-strict-aliasing
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

OutPutOpt     = -o

# Root

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLDFLAGS  := $(shell root-config --ldflags)

CLLFLAGS += $(ROOTLDFLAGS)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTGLIBS)
LIBS         += -lMinuit

INCLUDEFLAGS  = -I. -I./inc/
CXXFLAGS     += $(INCLUDEFLAGS)

ALLFLAGS = $(CXXFLAGS) $(CPPFLAGS) -Wall


# rule for any compiling any .cpp file
%.o : %.cpp  
	@printf "Compiling $< ... "
	@g++ $(ALLFLAGS) -c $<
	@echo "Done"

SimulatePointSourceInFP: SimulatePointSourceInFP.o GPilot.o GRootWriter.o GDefinition.o
	        $(LD)  $(CLLFLAGS) $(LIBS)  $(LDFLAGS) $^ $(OutPutOpt) $@ -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS
	        @echo "$@ done"

clean:	
	rm -f *.o

.SUFFIXES: .o


depend:
	$(DEPEND)

Makefile.depend:
	$(DEPEND)

DEPEND=echo > Makefile.depend0 &&\
	makedepend -s "\#DEPEND LIST DONT DELETE" -- $(INCLUDEFLAGS) \
		-Y --  *.cpp *.c \
		-f Makefile.depend0 > /dev/null 2>&1 &&\
	sed "s/^[a-zA-Z0-9]*\///" Makefile.depend0 > Makefile.depend &&\
	rm -f Makefile.depend0

include Makefile.depend
