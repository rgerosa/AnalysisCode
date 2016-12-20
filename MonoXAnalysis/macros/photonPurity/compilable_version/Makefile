ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

BINFILES =  MakeFitPurityNewId.cc 

PROGRAMS = $(patsubst %.cc,%,$(BINFILES))


# --- External configuration ---------------------------------
#CC         = g++ -Wno-deprecated 
CC         = g++ -Wall -pedantic 
#CC         = g++ -Wno-deprecated -pedantic
CCFLAGS    =  -g -O3 -std=c++11 -std=gnu++11
MFLAGS     = -MM
INCLUDES   =
WORKDIR    = tmp/
LIBDIR     = $(WORKDIR)
OBJDIR=$(WORKDIR)/objects/
# -------------------------------------------------------------

#ROOFITVER = 5.25.02-cms6
#ROOFITDIR = $(ROOTSYS)/../../roofit/$(ROOFITVER)

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')
ROOFIT_LIBS := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIB= | sed 's|LIB=||')
ROOFIT_LIBS += $(shell cd $(CMSSW_BASE); scram tool info roofit | grep LIB= | sed 's|LIB=||') 


INCLUDES += -I.  -I.. -I../include -I$(ROOTSYS)/include/  -I$(ROOFIT_INCLUDE)/
ROOTSYS  ?= ERROR_RootSysIsNotDefined

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -lRooStats

# CC files excluding the binaries
#CCFILES=(FitPurityNewId.C CMS_lumi.C)
CCFILES=$(filter-out $(BINFILES),$(wildcard *.C))

# List of all object files to build
OLIST=$(patsubst %.C,$(OBJDIR)/%.o,$(CCFILES))


# Implicit rule to compile all classes
$(OBJDIR)/%.o : %.C
	      @echo "Compiling $<"
	      @mkdir -p $(OBJDIR)
	      @$(CC) $(CCFLAGS) $(ROOTLIBS) $(EXTRALIBS) $(INCLUDES) -c $< -o $@  


$(PROGRAMS) : $(OLIST)
	    @echo "======================================"
	    @echo "Source file is  $(BINFILES)"	
	    @echo "Executable will be created in  directory \"$(WORKDIR)\""	
	    @echo "======================================"
	    @echo "Linking $@"
	    @$(CC) $(CCFLAGS)  $(INCLUDES) $(OLIST) \
	    $(ROOTLIBS) $(EXTRALIBS) -o $(WORKDIR)/$@   $(patsubst %,%.cc,$@)

default : Finalyzer

all : Finalyzer

clean:
	rm -Rf $(WORKDIR)/*
	@#rm -f $(LIBFILE)
	@rm -Rf *.o

veryclean:
	rm -Rf $(WORKDIR)

