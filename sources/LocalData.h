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

#ifndef SOURCES_LOCALDATA_H_
#define SOURCES_LOCALDATA_H_

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

#include "CodonSequenceAlignment.h"
#include "CodonStateSpace.h"
#include "GlobalParameters.h"
#include "SequenceAlignment.h"
#include "StringStreamUtils.h"
#include "Tree.h"

class LocalData {
 public:
  static const int Nnucp = 4;
  static const int Nnucrr = 6;
  static const int Ndinuc = 16;
  static const int Nstate_aa = 20;

  // to be set by global param
  double TOOSMALL;
  double TOOLARGE;
  double TOOLARGENEGATIVE;

  int NSummaries;
  int NParam;
  int NMapStats;

  string* listParam;
  string* listSummaries;
  string* listMapStats;

  std::map<string, int> mapUsedParam;
  std::map<string, int> mapUsedSummaries;
  std::map<string, int> mapUsedMapStats;
  std::map<string, int> mapUsedMapAncStats;

  int NusedMapStats;
  int NusedMapAncStats;
  int NusedParam;
  int NusedSummaries;
  int Ngenes;

  std::vector<string> listGenes;

  string localcontrolfile, output, model;
  std::vector<double> summariesRealData;

  void writeRealDataSummaries(ofstream& os, bool headers = true);
  void toFasta(ofstream& os, int** currentNodeleafCodonSequence);
  void toAli(ofstream& os, int** currentNodeleafCodonSequence);

  // GlobalParameters* gparam;
  string data;
  string tree;
  FileSequenceAlignment* dnadata;
  CodonStateSpace* codonstatespace;
  CodonSequenceAlignment* codondata;
  Tree* refTree;
  const TaxonSet* taxonset;
  bool iscodon;
  string code;
  int Nsite_codon, Nsite_nuc, Ntaxa, Nnode, Nstate_codon;

  void readLocalData(int k);
  explicit LocalData(GlobalParameters* gparam);
  virtual ~LocalData();
};

#endif  // SOURCES_LOCALDATA_H_
