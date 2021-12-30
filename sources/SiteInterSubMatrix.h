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
#ifndef SOURCES_SITEINTERSUBMATRIX_H_
#define SOURCES_SITEINTERSUBMATRIX_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "CodonSequenceAlignment.h"
#include "CodonStateSpace.h"
#include "EvolHistStatistics.h"
#include "LocalParameters.h"
#include "Random.h"
#include "SequenceAlignment.h"
#include "Tree.h"

class SiteInterSubMatrix {
 public:
  // SiteInterSubMatrix
  double*** submatrixTreeSim;
  double*** mutmatrixTreeSim;
  double* TotalSubRate;
  double* TotalMutRate;
  double* PartialSubRate;
  double* PartialMutRate;

  // parameters
  LocalParameters* lparam;

  // Constructor
  explicit SiteInterSubMatrix(LocalParameters* lparam);
  SiteInterSubMatrix(LocalParameters* lparam, std::string s);
  virtual ~SiteInterSubMatrix();

  std::tuple<double, double, double> ComputeCore(double MutRate, double SubRate,
                                                 double S, int* nucposFrom,
                                                 int* nucposTo, int codonPos,
                                                 int NodeIndex, int site_nuc,
                                                 int site_codon_i,
                                                 int** CurrentNodeNucSequence);

  // Getters
  double GetSubRate(int NodeIndex, int site_codon);
  double GetMutRate(int NodeIndex, int site_codon);
  std::tuple<int, int> getStartEndCodons(int site_codon);

  // Setters
  void setSubMatrix();
  void setSubMatrixFromLeaves();
  void resetSubMatrix();
  void resetSubMatrixFromLeaves();
  void transfertTotalRate(int sourceNodeIndex, int sinkNodeIndex);
  void findCodonContext(int NodeIndex, int site_nuc, int nucFrom, int nucTo,
                        int& pos1From, int& pos2From, int& pos3From,
                        int& pos1To, int& pos2To, int& pos3To,
                        int** CurrentNodeNucSequence);
  void ComputePartialRates(int NodeIndex, int site_codon,
                           int** CurrentNodeNucSequence);
  void UpdateSubMatrixTreeSim(int NnodeIndex, int site_codon,
                              int** CurrentNodeNucSequence);
  void UpdateSubMatrixFromLeaves(int NnodeIndex, int** CurrentNodeNucSequence);
  int testGpTcontext(int NnodeIndex, int site, int nucFrom, int nucTo,
                     int** CurrentNodeNucSequence);
  int testCpGcontext(int NnodeIndex, int site, int nucFrom, int nucTo,
                     int** CurrentNodeNucSequence);
  int testTpAcontext(int NnodeIndex, int site, int nucFrom, int nucTo,
                     int** CurrentNodeNucSequence);
  int testGCPref(int innucFrom, int innucTo);
  int testContextDinuc(int NodeIndex, int site_nuc, int* context, int nucTo,
                       int** CurrentNodeNucSequence);

  void transfertNodeMatrix(int sourceNodeIndex, int sinkNodeIndex,
                           int site_nuc);

  // writers
  void WriteSubMatrix(ostream& mutation_os, ostream& selection_os, int nucTo);

  ////
  // Getters TotalRates
  ////

  double GetTotalMutRate(int NodeIndex) { return TotalMutRate[NodeIndex]; }
  double GetTotalSubRate(int NodeIndex) { return TotalSubRate[NodeIndex]; }

  ////
  // Getters PartialRates
  ////
  double GetPartialMutRate(int NodeIndex) { return PartialMutRate[NodeIndex]; }
  double GetPartialSubRate(int NodeIndex) { return PartialSubRate[NodeIndex]; }
  double GetSubRate(int NodeIndex, int site_nuc, int nucTo) {
    return submatrixTreeSim[NodeIndex][site_nuc][nucTo];
  }
  double GetMutRate(int NodeIndex, int site_nuc, int nucTo) {
    return mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
  }
  bool CheckStop(int pos1, int pos2, int pos3) {
    if (lparam->codonstatespace->CheckStop(pos1 = pos1, pos2 = pos2,
                                           pos3 = pos3)) {
      std::cerr << pos1 << " " << pos2 << " " << pos3 << " "
                << "\n";
      exit(0);
    }
    return true;
  }

  double ComputeFixationFactor(double S, double SubRate);
};

#endif  // SOURCES_SITEINTERSUBMATRIX_H_
