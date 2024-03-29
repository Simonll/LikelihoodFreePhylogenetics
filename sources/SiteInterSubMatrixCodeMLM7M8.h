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
#ifndef SOURCES_SITEINTERSUBMATRIXCODEMLM7M8_H_
#define SOURCES_SITEINTERSUBMATRIXCODEMLM7M8_H_

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

class SiteInterSubMatrixCodeMLM7M8 : public SiteInterSubMatrix {
 public:
  using SiteInterSubMatrix::SiteInterSubMatrix;
  std::tuple<double, double, double> ComputeCore(
      int* nucposFrom, int* nucposTo, int codonPos, int NodeIndex, int site_nuc,
      int site_codon_i, int** CurrentNodeNucSequence) override;
};

#endif  // SOURCES_SITEINTERSUBMATRIXCODEMLM7M8_H_
