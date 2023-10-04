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
#ifndef SOURCES_LOCALPARAMETERS_H_
#define SOURCES_LOCALPARAMETERS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "BranchSpecificParameters.h"
#include "CodonSequenceAlignment.h"
#include "CodonStateSpace.h"
#include "GlobalParameters.h"
#include "Random.h"
#include "SequenceAlignment.h"
#include "StringStreamUtils.h"
#include "Tree.h"

class LocalParameters {
 public:
  static const int Nnucp = 4;
  static const int Nnucrr = 6;
  static const int Ndinuc = 16;
  static const int Nstate_aa = 20;
  static constexpr double TOOSMALL = 1e-30;
  static constexpr double TOOLARGE = 500;
  static constexpr double TOOLARGENEGATIVE = -500;

  // to be set by global param

  int NSummaries;
  int NParam;
  int NEvoStats;
  int NSiteSpecificEvoStats;

  string *listParam;
  string *listSummaries;
  string *listEvoStats;
  string *listSiteSpecificEvoStats;
  string *listSpecies;

  std::map<string, int> mapUsedParam;
  std::map<string, int> mapUsedSummaries;
  std::map<string, int> mapUsedAncSummaries;
  std::map<string, int> mapUsedAccessorySummaries;
  std::map<string, int> mapUsedEvoStats;
  std::map<string, int> mapUsedEvoAncStats;
  std::map<string, int> mapUsedSiteSpecificEvoStats;

  int NusedEvoStats;
  int NusedSiteSpecificEvoStats;
  int NusedEvoAncStats;
  int NusedParam;
  int NusedSummaries;
  int NusedAncSummaries;
  int NusedAccessorySummaries;
  int Ngenes;

  string data, localcontrolfile, output, model;

  BranchSpecificParameters **bparam;

  FileSequenceAlignment *dnadata;
  CodonStateSpace *codonstatespace;
  CodonSequenceAlignment *codondata;
  Tree *refTree;
  const TaxonSet *taxonset;
  Link *outgroupLink;
  Link *newlink;
  Link *newnext;
  Branch *branchToInGroup;
  Branch *branchToOutGroup;
  Node *newnode;
  Random *rnd;

  bool iscodon, isdata;
  std::map<int, int> gtrMap;
  int gtr1NodeIndex;
  int gtr2NodeIndex;
  int Ninterval;
  int sampleAncSeq;
  double rootlength, percentFromOutGroup, branchLengthBetweenInAndOutGroup;

  int Nsite_codon, Nsite_nuc, Ntaxa, Nstate_codon, startPoint, endPoint,
      everyPoint, tofasta, Nrep, N_profile;
  string controlfile, code, chain, posteriorfile, taxa_a, taxa_b, taxa_gtr1_a,
      taxa_gtr1_b, taxa_gtr1_c, taxa_gtr2_a, taxa_gtr2_b, taxa_gtr2_c,
      transformation;

  bool getrate;
  bool getrate1;
  bool getrate2;

  double *site_omega;
  double omega;
  double *nucp, *nucrr, *nucp1, *nucrr1, *nucp2, *nucrr2;  // 6 param
  double **nucrrnr, **nucrrnr1, **nucrrnr2;                // 12 param
  double **gtnr, **gtnr1, **gtnr2;
  double **ssaaprofiles;
  double **sscodonprofiles;
  double *codonprofile;
  int *alloc;

  // param value
  double lambda_TBL, lambda_omega, lambda_CpG, lambda_TpA, MutationNormFactor,
      MutationNormFactor1, MutationNormFactor2, fitCpG, fitTpA, lambda_tvCpG,
      lambda_tvTpA, lambda_tstvCpG, lambda_tstvTpA, fitGC, lambda_R;
  double *muBranch;
  double *AAadj;
  double Aadj, Cadj, Dadj, Eadj, Fadj, Gadj, Hadj, Iadj, Kadj, Ladj, Madj, Nadj,
      Padj, Qadj, Radj, Sadj, Tadj, Vadj, Wadj, Yadj;

  // swhitch fix or free param
  int fixNsite, fixomega, fixlambda_omega, fixlambda_TBL, fixlambda_CpG,
      fixlambda_TpA, fixgtr, fixgtr1, fixgtr2, fixgtnr, fixstat, fixts, fixtr,
      fixrr, fixkappa, fixhky, randomseed, verbose, rooted, fixroot, fixss,
      fixfitCpG, fixlambda_tvCpG, fixlambda_tvTpA, fixlambda_tstvCpG,
      fixlambda_tstvTpA, fixfitTpA, fixfitGC, fixlambda_R, fixAAadj,
      fixCODONadj;
  int MCMCpointID;

  string lambda_TBL_prior, lambda_CpG_prior, lambda_TpA_prior,
      lambda_omega_prior;

  std::vector<double> summariesRealData;
  std::vector<double> accessorysummariesRealData;
  std::vector<double> summariesSimulatedData;
  std::vector<double> accessorysummariesSimulatedData;
  std::vector<double> summariesAncestralData;
  std::vector<double> evostats;
  std::vector<double> ancevostats;
  std::vector<double> sitespecificevostats;

  // Constructor
  explicit LocalParameters(GlobalParameters *gparam);
  virtual ~LocalParameters();

  // Writers
  void writeRealDataSummaries(ofstream &os, bool headers = true);
  void writeAncestralDataSummaries(ofstream &os, bool headers);

  void tobstats(ofstream &os);
  void tobstats(ofstream &os, const Link *from);
  void toSsstats(ofstream &os);
  void toFasta(ofstream &os, int **currentNodeleafCodonSequence);
  void toAli(ofstream &os, int **currentNodeleafCodonSequence);

  // Readers
  void readChainCodonMutSelSBDP();
  void readChainCodonMutSelSBDP(int pt_i);
  int readParametersMutSelC(int pt_i);
  int readBayescodeParametersMutSelAAC(int pt_i);
  int readBayescodeParametersMutSelC(int pt_i);

  void readChainCodonMutSelFinite();
  void readChainCodonMutSelFinite(int pt_i);

  void readFMutSelCodeML();
  void readLocalInstructions();
  int readParametersCodemlM7M8(int pt_i);

  // Setters
  void SetCurrentParametersFromPosterior(
      std::vector<std::vector<double>> posterior, int it);
  void SetTree();
  void SetBranchesLengthsBetweenInAndOutGroup();
  void SetRootBetweenInAndOutGroup();
  void SetRootLCA();
  void SetTreeStuff();
  void SetTreeStuffRecursively(Link *from, int notNodeIndex, int gtrIndex);
  void Setgtr2gtnr();

  // Getters
  int GetPointID();
  std::vector<double> GetCurrentParameters();
  std::vector<double> GetCurrentSummaries();
  std::vector<double> GetCurrentAccessorySummaries();
  std::vector<double> GetCurrentEvoStats();
  std::vector<double> GetCurrentAncEvoStats();
  std::vector<double> GetCurrentSiteSpecificEvoStats();
  void GetGTR1();
  void GetGTR2();
  void incrementStartPoint() { this->startPoint++; }

  int get_number_of_samplesCodemlM7M8(std::ifstream &is) {
    std::string line = "";
    int k = 0;
    while (std::getline(is, line)) {
      std::getline(is, line);
      std::getline(is, line);
      std::getline(is, line);
      k++;
    }
    return k;
  }

  int get_number_of_samples(std::ifstream &is) {
    std::string line;
    int k = 0;
    std::getline(is, line);
    while (getline(is, line)) {
      while (getline(is, line)) {
        if (line.empty()) continue;
        if (line[0] == '(') break;
      }
      k++;
    }
    return k;
  }

  int GetNsiteCodon() { return Nsite_codon; }

  double GetGTNR(int i, int j) {
    return this
        ->nucrrnr[i][j];  // (2*3);//((this->nucrrnr[i][j])/GetRateGTNR());
  }
  double GetGTR(int i, int j) {
    return ((this->nucp[j] * this->nucrrnr[i][j]) / GetRate());
  }

  double GetGTRCodeML(int i, int j) {
    return (this->nucp[j] * this->nucrrnr[i][j]);
  }

  int GetNucRRIndex(int i, int j) {
    return (i < j) ? (2 * 4 - i - 1) * i / 2 + j - i - 1
                   : (2 * 4 - j - 1) * j / 2 + i - j - 1;
  }

  double GetRateGTNR() {
    if (getrate) {
      return MutationNormFactor;
    }
    double norm = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        if (i != j) {
          norm += this->nucrrnr[i][j];
        }
      }
    }
    getrate = true;
    MutationNormFactor = 2 * (norm * 3);
    return MutationNormFactor;
  }

  double GetRate() {
    if (getrate) {
      return MutationNormFactor;
    }
    double norm = 0.0;
    for (int i = 0; i < 4 - 1; i++) {
      for (int j = i + 1; j < 4; j++) {
        norm += this->nucp[i] * this->nucp[j] * this->nucrrnr[i][j];
      }
    }
    // 2 for the symetry of the matrix??, and 3 for the number of codon positons
    getrate = true;
    MutationNormFactor = 2 * (norm * 3);
    return MutationNormFactor;
  }
};

#endif  // SOURCES_LOCALPARAMETERS_H_
