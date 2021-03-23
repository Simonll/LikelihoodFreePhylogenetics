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
#ifndef SOURCES_PRIORSAMPLER_H_
#define SOURCES_PRIORSAMPLER_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "BiologicalSequences.h"
#include "LocalParameters.h"
#include "Random.h"

class PriorSampler {
 public:
  LocalParameters* lparam;

  explicit PriorSampler(LocalParameters* lparam);
  virtual ~PriorSampler();

  void sample();

 protected:
  double Unif(double inf, double sup) {
    return lparam->rnd->Uniform() * (sup - inf) + inf;
  }

  double log10Unif() { return pow(10, lparam->rnd->Uniform() * 2 - 1); }

  double logNUnif(int n) { return pow(n, lparam->rnd->Uniform() * 2 - 1); }

  double log50Unif() { return pow(50, lparam->rnd->Uniform() * 2 - 1); }

  double logUnifTruncated() {
    double a = pow(100, lparam->rnd->Uniform() * 2 - 1);
    // 1 < x <= 4
    if (a <= 20.0 && a > 5.0) {
      a /= 5;

      // 4 < x <= 20
    } else if (a > 20) {
      a /= 5;
    }
    return a;
  }
};

#endif  // SOURCES_PRIORSAMPLER_H_
