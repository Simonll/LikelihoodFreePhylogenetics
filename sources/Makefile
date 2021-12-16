CC=g++
CPPFLAGS=-c -g -std=c++11 -O2 -Wall -pedantic -fopenmp -fno-var-tracking -Wfatal-errors
LDFLAGS=-fopenmp
SRCS=GlobalParameters.cpp TaxonSet.cpp Tree.cpp  Random.cpp  SequenceAlignment.cpp CodonSequenceAlignment.cpp\
  StateSpace.cpp CodonStateSpace.cpp  LocalParameters.cpp SummaryStatistics.cpp SiteInterSubMatrix.cpp SiteInterSubMatrixCodeMLM7M8.cpp SiteInterSubMatrixCABC2018.cpp SiteInterSubMatrixCodeMLFMUTSEL.cpp EvolHistStatistics.cpp\
  PriorSampler.cpp AncestralSequence.cpp TreeSimulator.cpp BranchSpecificParameters.cpp  LocalData.cpp\
  Posterior.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

#OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=../data/LFP
PROGSDIR=../data
ALL= LFP codemlM7M8
#all: $(SOURCES) $(EXECUTABLE)
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

EXEEXT=
uname := $(shell uname -s 2>/dev/null)
ifeq (MINGW,$(findstring MINGW,$(uname)))
	EXEEXT=.exe
	PROGS=$(addsuffix $(EXEEXT), $(ALL))
endif

.PHONY: all clean
all: $(PROGS)

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$

%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

LFP$(EXEEXT): $(PROGSDIR)/LFP$(EXEEXT)
$(PROGSDIR)/LFP$(EXEEXT): main.o $(OBJS)
	$(CC) main.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codemlM8M7$(EXEEXT): $(PROGSDIR)/codemlM7M8$(EXEEXT)
$(PROGSDIR)/codemlM7M8$(EXEEXT): codemlM7M8.o $(OBJS)
	$(CC) codemlM7M8.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@


#$(EXECUTABLE): $(OBJECTS)
#	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

#.cpp.o:
#	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm *.d
	-rm *.o
	-rm -f $(PROGS)
