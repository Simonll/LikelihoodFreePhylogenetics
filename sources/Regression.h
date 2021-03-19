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
#ifndef REGRESSION_H
#define REGRESSION_H

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

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

#include "GlobalParameters.h"
#include "Posterior.h"
#include "linalg.h"

class Regression {
 public:
  GlobalParameters* gparam;
  Posterior* post;

  // local regression adjustment
  //(X_T dot W dot X)^{-1} dot X_T dot W dot Y
  std::vector<std::vector<double>> Y_HAT;
  std::vector<std::vector<double>> Yres;   // Matrix of residual = Y - Y_HAT
  std::vector<std::vector<double>> B_HAT;  // Matrix of regression coefficient

  Regression(GlobalParameters* gparam, Posterior* post);
  virtual ~Regression();

  // Getters
  double GetRSquare(int i);
  double GetAjustedRSquare(int i);
  double GetRSquarePred(int i);
  double GetS(int i);
  double** GetXWX_t();
  double** GetX_tWY();

  void ComputeMultipleRegression();
  void ComputeX_TWY();

 protected:
 private:
};

#endif  // REGRESSION_H
