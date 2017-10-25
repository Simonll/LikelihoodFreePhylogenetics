#ifndef TREESIMULATOR_H
#define TREESIMULATOR_H

#include <iostream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

#define isnan std::isnan
#define isinf std::isinf
#define string std::string
#define ostream std::ostream
#define ofstream std::ofstream
#define istream std::istream
#define ifstream std::ifstream
#define cin std::cin
#define cerr std::cerr
#define cout std::cout
#define setw std::setw
#define ostringstream std::ostringstream
#define istringstream std::istringstream
#define IOS_APPEND std::ios_base::app
#define APPEND std::ios_base::app
#define OUT std::ios_base::out

#include "Random.h"
#include "Tree.h"
#include "SequenceAlignment.h"
#include "CodonStateSpace.h"
#include "CodonSequenceAlignment.h"
#include "LocalParameters.h"
#include "EvolHistStatistics.h"
#include "SiteInterSubMatrix.h"
#include "AncestralSequence.h"


class TreeSimulator
{
public:




    //Simulating along a tree
    int*** CurrentAncestralCodonSequence;
    int** CurrentNodeCodonSequence;
    int** CurrentNodeNucSequence;
    int** CurrentLeafNodeNucSequence;
    int** CurrentLeafNodeCodonSequences;

    //parameters
    LocalParameters* lparam;

    //Evolutionary Statistics
    EvolHistStatistics* treeEvoStats;
    EvolHistStatistics* rootBranchEvoStats;


    //Substitution Matrix
    SiteInterSubMatrix* submatrix;

    //Ancestral Sequence
    AncestralSequence* ancestralseq;


    //Simulation functions
    void UpdateSubMatrixTreeSim(int NnodeIndex, int site_codon);
    void RegisterSubTreeSim(int NodeIndex, int site_nuc, int nucTo);

    void ComputeRecursiveSimulation(Link* from);
    void resetSimulator();
    void resetEvoStatVectors();
    void WriteAncestral();

    //Setters
    void SetAncestralSequence();


    //Getters
    void GetNewSimulatedCodonAlignment();
    void GetSampleAncestralCodonSequence(int FromNodeIndex,int interval);


    TreeSimulator(LocalParameters* lparam,SiteInterSubMatrix* submatrix,AncestralSequence* ancestralseq);
    ~TreeSimulator();

protected:

private:
};

#endif // TREESIMULATOR_H



