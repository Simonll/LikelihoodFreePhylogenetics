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
#ifndef SOURCES_POSTERIORSAMPLER_H_
#define SOURCES_POSTERIORSAMPLER_H_

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

#include "LocalParameters.h"
#include "Random.h"

class PosteriorSampler {
 public:
  LocalParameters* lparam;
  PosteriorSampler();
  virtual ~PosteriorSampler();

 protected:
  void sample();
  std::tuple<double, double> sampleFromAjustedDensity(int param, double inf,
                                                      double sup);

 private:
};

#endif  // SOURCES_POSTERIORSAMPLER_H_
