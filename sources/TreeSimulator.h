/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/
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
    void ComputeSeqProb();
    void resetSimulator();
    void resetSimulatorSeq();
    void resetEvoStatVectors();
    void WriteAncestral();

    //Setters
    void SetAncestralSequence();

    //Getters
    void GetNewSimulatedCodonAlignment();
    void GetNewProbSeq();
    void GetSampleAncestralCodonSequence(int FromNodeIndex,int interval);

    TreeSimulator(LocalParameters* lparam,SiteInterSubMatrix* submatrix,AncestralSequence* ancestralseq);
    TreeSimulator(LocalParameters* lparam,SiteInterSubMatrix* submatrix);
    ~TreeSimulator();

protected:

private:
};

#endif // TREESIMULATOR_H



