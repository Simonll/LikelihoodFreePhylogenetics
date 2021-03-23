CC=g++
CFLAGS=-c -g -std=c++17 -O2 -Wall -pedantic -fopenmp -fno-var-tracking -Wfatal-errors
LDFLAGS=-fopenmp
SOURCES=main.cpp GlobalParameters.cpp TaxonSet.cpp Tree.cpp  Random.cpp  SequenceAlignment.cpp CodonSequenceAlignment.cpp\
  StateSpace.cpp CodonStateSpace.cpp  LocalParameters.cpp SummaryStatistics.cpp SiteInterSubMatrix.cpp EvolHistStatistics.cpp\
  PriorSampler.cpp AncestralSequence.cpp TreeSimulator.cpp BranchSpecificParameters.cpp  LocalData.cpp\
  Posterior.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../data/LFP
PROGSDIR=../data

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o
