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
#ifndef SOURCES_POSTERIOR_H_
#define SOURCES_POSTERIOR_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "GlobalParameters.h"
#include "Random.h"

class Posterior {
 public:
  explicit Posterior(GlobalParameters* gparam);
  virtual ~Posterior();

  static const int chainIDGetter = 0;
  static const int paramGetter = 1;
  static const int summariesGetter = 2;
  static const int accsummariesGetter = 3;
  static const int evoancstatsGetter = 4;
  static const int evostatsGetter = 5;
  static const int ssevostatsGetter = 6;

  double TOOSMALL;
  double TOOLARGE;
  double TOOLARGENEGATIVE;

  int verbose;
  int NSummaries;
  int NParam;
  int NEvoStats;
  int NSiteSpecificEvoStats;
  int NAccessorySummaries;
  int NEvoAncStats;

  string* listParam;
  string* listSummaries;
  string* listEvoStats;
  string* listSiteSpecificEvoStats;

  std::map<string, int> mapUsedParam;
  std::map<string, int> mapUsedSummaries;
  std::map<string, int> mapUsedAccessorySummaries;
  std::map<string, int> mapUsedEvoStats;
  std::map<string, int> mapUsedSiteSpecificEvoStats;
  std::map<string, int> mapUsedEvoAncStats;

  int NusedEvoStats;
  int NusedSiteSpecificEvoStats;
  int NusedEvoAncStats;
  int NusedParam;
  int NusedSummaries;
  int NusedAccessorySummaries;
  int Ngenes;

  string localcontrolfile, output, model;

  std::vector<std::tuple<int, std::vector<double>, std::vector<double>,
                         std::vector<double>, std::vector<double>,
                         std::vector<double>, std::vector<double>>>
      population_t;
  std::vector<std::vector<double>> posterior;
  double* empVar;
  double* empMean;

  Random* rnd;
  int randomseed;

  int Niter, Nrun, threshold, Nthread, Nsite_codon;
  bool sorted;
  // writter
  void writePosterior(ofstream& os);
  void writeSimu(ofstream& os);
  void writePosterior(ofstream& os, int Nsimu);
  void writePosteriorPredictiveStatistics(
      ofstream& os, std::vector<double> realDataSummaries);
  void writeHeader(ofstream& os);
  // readers
  void readPosterior(string posteriorfile);
  void readPosterior(ifstream& is);

  // Getters
  int PosteriorGetSize();
  void GetEmpVar();

  // Setters
  void registerSimulation(int chainID, std::vector<double> param,
                          std::vector<double> summaries,
                          std::vector<double> accsummaries,
                          std::vector<double> ancevostat,
                          std::vector<double> evostat,
                          std::vector<double> ssevostat);

  void SetNsite(int i);

  int GetSize() { return population_t.size(); }
};

#endif  // SOURCES_POSTERIOR_H_
