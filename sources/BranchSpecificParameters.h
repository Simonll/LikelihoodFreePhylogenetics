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
#ifndef SOURCES_BRANCHSPECIFICPARAMETERS_H_
#define SOURCES_BRANCHSPECIFICPARAMETERS_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class BranchSpecificParameters {
 public:
  static const int Nnucp = 4;
  static const int Nnucrr = 6;
  static const int Ndinuc = 16;
  double *nucp, *nucrr;  // 6 param
  double** nucrrnr;      // 12 param
  double** gtnr;

  BranchSpecificParameters();
  virtual ~BranchSpecificParameters();

  void SetLocalParaemters(double* nucp, double* nucrr);
};

#endif  // SOURCES_BRANCHSPECIFICPARAMETERS_H_
