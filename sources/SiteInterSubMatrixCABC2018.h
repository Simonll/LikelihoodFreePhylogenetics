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
#ifndef SOURCES_SITEINTERSUBMATRIXCABC2018_H_
#define SOURCES_SITEINTERSUBMATRIXCABC2018_H_

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
#include "SiteInterSubMatrix.h"
#include "Tree.h"

class SiteInterSubMatrixCABC2018 : public SiteInterSubMatrix {
 public:
  double*** selmatrixTreeSim;
  double* TotalSubRateNonSyn;
  double* TotalMutRateNonSyn;
  double* TotalSubRateSyn;
  double* TotalMutRateSyn;
  double* PartialSubRateNonSyn;
  double* PartialMutRateNonSyn;
  double* PartialSubRateSyn;
  double* PartialMutRateSyn;

  void setSubMatrix() override;
  void resetSubMatrix() override;
  void init() override;
  explicit SiteInterSubMatrixCABC2018(LocalParameters* lparam)
      : SiteInterSubMatrix(lparam) {}
  virtual ~SiteInterSubMatrixCABC2018();

  std::tuple<double, double, double> ComputeCore(int* nucposFrom, int* nucposTo,
                                                 int codonPos, int NodeIndex,
                                                 int site_nuc, int site_codon_i,
                                                 int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesCpG(int NodeIndex, int site_codon,
                                         int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesNonSynCpG(int NodeIndex, int site_codon,
                                               int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesSynCpG(int NodeIndex, int site_codon,
                                            int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesNonSyn(int NodeIndex, int site_codon,
                                            int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesSyn(int NodeIndex, int site_codon,
                                         int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesWeakStrong(int NodeIndex, int site_codon,
                                                int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesStrongWeak(int NodeIndex, int site_codon,
                                                int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesWeakWeak(int NodeIndex, int site_codon,
                                              int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesStrongStrong(int NodeIndex, int site_codon,
                                                  int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesTransition(int NodeIndex, int site_codon,
                                                int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesTransversion(int NodeIndex, int site_codon,
                                                  int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesConsPol(int NodeIndex, int site_codon,
                                             int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesRadPol(int NodeIndex, int site_codon,
                                            int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesConsVol(int NodeIndex, int site_codon,
                                             int** CurrentNodeNucSequence);

  std::tuple<double, double> GetRatesRadVol(int NodeIndex, int site_codon,
                                            int** CurrentNodeNucSequence);

  bool isWeakStrong(int nucFrom, int nucTo);
  bool isTransition(int nucFrom, int nucTo);
  bool isStrongWeak(int nucFrom, int nucTo);
  bool isStrongStrong(int nucFrom, int nucTo);
  bool isWeakWeak(int nucFrom, int nucTo);
};

#endif  // SOURCES_SITEINTERSUBMATRIXCABC2018_H_
