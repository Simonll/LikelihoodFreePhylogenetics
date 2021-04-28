/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version. LikelihoodFreePhylogenetics is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details. You should have received a copy of the GNU
General Public License along with LikelihoodFreePhylogenetics. If not, see
<http://www.gnu.org/licenses/>.
*/
#ifndef SOURCES_TREESIMULATOR_H_
#define SOURCES_TREESIMULATOR_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AncestralSequence.h"
#include "CodonSequenceAlignment.h"
#include "CodonStateSpace.h"
#include "EvolHistStatistics.h"
#include "LocalParameters.h"
#include "Random.h"
#include "SequenceAlignment.h"
#include "SiteInterSubMatrix.h"
#include "Tree.h"

class TreeSimulator {
 public:
  // Simulating along a tree
  int*** CurrentAncestralCodonSequence;
  int** CurrentNodeCodonSequence;
  int** CurrentNodeNucSequence;
  int** CurrentLeafNodeNucSequence;
  int** CurrentLeafNodeCodonSequences;

  // parameters
  LocalParameters* lparam;

  // Evolutionary Statistics
  EvolHistStatistics* treeEvoStats;
  EvolHistStatistics* rootBranchEvoStats;

  // Substitution Matrix
  SiteInterSubMatrix* submatrix;

  // Ancestral Sequence
  AncestralSequence* ancestralseq;

  // Simulation functions
  void UpdateSubMatrixTreeSim(int NnodeIndex, int site_codon);
  void RegisterSubTreeSim(int NodeIndex, int site_nuc, int nucTo);
  void ComputeRecursiveSimulation(Link* from);
  void ComputeSeqProb();
  void resetSimulator();
  void resetSimulatorSeq();
  void resetEvoStatVectors();
  void WriteAncestral();

  // Setters
  void SetAncestralSequence();
  void SetAncestralCodonSequence(int FromNodeIndex, int interval);

  // Getters
  void GenerateCodonAlignment();
  void GetNewProbSeq();

  TreeSimulator(LocalParameters* lparam, SiteInterSubMatrix* submatrix,
                AncestralSequence* ancestralseq);
  TreeSimulator(LocalParameters* lparam, SiteInterSubMatrix* submatrix);
  ~TreeSimulator();
};

#endif  // SOURCES_TREESIMULATOR_H_
