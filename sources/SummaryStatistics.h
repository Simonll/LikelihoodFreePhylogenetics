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
#ifndef SOURCES_SUMMARYSTATISTICS_H_
#define SOURCES_SUMMARYSTATISTICS_H_

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

#include "BiologicalSequences.h"
#include "CodonSequenceAlignment.h"
#include "CodonStateSpace.h"
#include "LocalData.h"
#include "LocalParameters.h"
#include "Random.h"
#include "SequenceAlignment.h"
#include "StringStreamUtils.h"
#include "Tree.h"

class SummaryStatistics {
 public:
  // Parameters
  LocalParameters* lparam;
  LocalData* ldata;
  typedef double (SummaryStatistics::*funcpt)(
      CodonSequenceAlignment* codondata);
  std::map<string, funcpt> GetSummariesMap;

  // Summaries
  double** dinuc_usage;
  double** dinuc12_usage;
  double** dinuc23_usage;
  double** dinuc31_usage;
  double** dicodon_usage;
  double** diaa_usage;
  double* RSCU;
  double* relativeCodonFrequency;
  double* codon_usage_wonR;
  double* relativeAAFrequency;
  double* nuc_usage;
  double* nuc1_usage;
  double* nuc2_usage;
  double* nuc3_usage;
  double* nuc_meandiff;
  double* nuc1_meandiff;
  double* nuc2_meandiff;
  double* nuc3_meandiff;
  double* codon_meandiff;
  double* aa_meandiff;
  double* aa_wonR_meandiff;
  double* CGNAGR;
  double* nuc_pairwiseSyn;
  double* nuc_pairwiseSyn10;
  double* nuc_pairwiseSyn30;
  double* nuc_pairwiseSyn50;
  double* nuc_pairwiseSyn70;
  double* nuc_pairwiseSyn90;
  double* nuc_pairwise;
  double* nuc_pairwise10;
  double* nuc_pairwise30;
  double* nuc_pairwise50;
  double* nuc_pairwise70;
  double* nuc_pairwise90;
  double* nuc1_pairwise;
  double* nuc1_pairwise10;
  double* nuc1_pairwise30;
  double* nuc1_pairwise50;
  double* nuc1_pairwise70;
  double* nuc1_pairwise90;
  double* nuc2_pairwise;
  double* nuc2_pairwise10;
  double* nuc2_pairwise30;
  double* nuc2_pairwise50;
  double* nuc2_pairwise70;
  double* nuc2_pairwise90;
  double* nuc3_pairwise;
  double* nuc3_pairwise10;
  double* nuc3_pairwise30;
  double* nuc3_pairwise50;
  double* nuc3_pairwise70;
  double* nuc3_pairwise90;
  double* aa_pairwise;
  double* aa_pairwise10;
  double* aa_pairwise30;
  double* aa_pairwise50;
  double* aa_pairwise70;
  double* aa_pairwise90;
  int* dinucCpG_pairwise;

  double nuc_site_comphet;
  double nuc1_site_comphet;
  double nuc2_site_comphet;
  double nuc3_site_comphet;
  double nuc_taxa_comphet;
  double nuc1_taxa_comphet;
  double nuc2_taxa_comphet;
  double nuc3_taxa_comphet;
  double codon_site_comphet;
  double codon_taxa_comphet;
  double aa_site_comphet;
  double aa_taxa_comphet;
  double RSCUentropy;
  double GC, GC1, GC2, GC3;

  bool RSCU_bool;
  bool relativeCodonFrequency_bool;
  bool codon_wonR_bool;
  bool dinuc_bool;
  bool dinuc12_bool;
  bool dinuc23_bool;
  bool dinuc31_bool;
  bool relativeAAFrequency_bool;
  bool dicodon_bool;
  bool diaa_bool;
  bool nuc_bool;
  bool nuc1_bool;
  bool nuc2_bool;
  bool nuc3_bool;
  bool nuc_meandiff_bool;
  bool nuc1_meandiff_bool;
  bool nuc2_meandiff_bool;
  bool nuc3_meandiff_bool;
  bool codon_meandiff_bool;
  bool aa_meandiff_bool;
  bool nuc_pairwise_bool;
  bool nuc_pairwise_bool10;
  bool nuc_pairwise_bool30;
  bool nuc_pairwise_bool50;
  bool nuc_pairwise_bool70;
  bool nuc_pairwise_bool90;
  bool nuc_pairwiseSyn_bool;
  bool nuc_pairwiseSyn_bool10;
  bool nuc_pairwiseSyn_bool30;
  bool nuc_pairwiseSyn_bool50;
  bool nuc_pairwiseSyn_bool70;
  bool nuc_pairwiseSyn_bool90;
  bool nuc1_pairwise_bool;
  bool nuc1_pairwise_bool10;
  bool nuc1_pairwise_bool30;
  bool nuc1_pairwise_bool50;
  bool nuc1_pairwise_bool70;
  bool nuc1_pairwise_bool90;
  bool nuc2_pairwise_bool;
  bool nuc2_pairwise_bool10;
  bool nuc2_pairwise_bool30;
  bool nuc2_pairwise_bool50;
  bool nuc2_pairwise_bool70;
  bool nuc2_pairwise_bool90;
  bool nuc3_pairwise_bool;
  bool nuc3_pairwise_bool10;
  bool nuc3_pairwise_bool30;
  bool nuc3_pairwise_bool50;
  bool nuc3_pairwise_bool70;
  bool nuc3_pairwise_bool90;
  bool aa_pairwise_bool;
  bool aa_pairwise_bool10;
  bool aa_pairwise_bool30;
  bool aa_pairwise_bool50;
  bool aa_pairwise_bool70;
  bool aa_pairwise_bool90;
  bool dinucCpG_pairwise_bool;
  bool nuc_site_comphet_bool;
  bool nuc1_site_comphet_bool;
  bool nuc2_site_comphet_bool;
  bool nuc3_site_comphet_bool;
  bool nuc_taxa_comphet_bool;
  bool nuc1_taxa_comphet_bool;
  bool nuc2_taxa_comphet_bool;
  bool nuc3_taxa_comphet_bool;
  bool codon_site_comphet_bool;
  bool codon_taxa_comphet_bool;
  bool aa_site_comphet_bool;
  bool aa_taxa_comphet_bool;
  bool RSCUentropy_bool;
  bool GC_bool, GC1_bool, GC2_bool, GC3_bool;

  // Constructors
  explicit SummaryStatistics(LocalParameters* lparam);
  explicit SummaryStatistics(LocalData* ldata);

  // Destructors
  virtual ~SummaryStatistics();

  void MapFunctions();
  void computeSummariesFromData();
  void computeSummaries();
  void computeSummaries(int** CurrentNodeLeafCodonSequence);
  void computeSummariesAncestralSequence(int** CurrentAncestralCodonSequence);
  double transformSummaryStatistics(double s);
  void GetRealStat();
  //        std::vector<double> ReadRealStat(string inrealstat) {
  //            std::vector<double> cur_real_statistics;
  //            ostringstream buffer;
  //            buffer << inrealstat;
  //            ifstream is(buffer.str());
  //            if (!is)       {
  //                std::cerr << "error: did not find " << buffer.str() << "\n";
  //                exit(1);
  //            }
  //
  //            ReadLine(is); // Read headers
  //            double cur_dl;
  //            while(is.eof()) {
  //                is >> cur_dl;
  //                cur_real_statistics.push_back(cur_dl);
  //
  //            }
  //
  //            is.close();
  //            return cur_real_statistics;
  //
  //        }

  /////////////////
  // compositional heterogenity
  /////////////////
  double Getnuc_site_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc_site_comphet_bool) {
      cur_v = codondata->nuc_site_comphet();
      nuc_site_comphet_bool = true;
    }
    return cur_v;
  }
  double Getnuc1_site_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc1_site_comphet_bool) {
      cur_v = codondata->nuc1_site_comphet();
      nuc1_site_comphet_bool = true;
    }
    return cur_v;
  }
  double Getnuc2_site_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc2_site_comphet_bool) {
      cur_v = codondata->nuc2_site_comphet();
      nuc2_site_comphet_bool = true;
    }
    return cur_v;
  }
  double Getnuc3_site_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc3_site_comphet_bool) {
      cur_v = codondata->nuc3_site_comphet();
      nuc3_site_comphet_bool = true;
    }
    return cur_v;
  }
  double Getnuc_taxa_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc_taxa_comphet_bool) {
      cur_v = codondata->nuc_taxa_comphet();
      nuc_taxa_comphet_bool = true;
    }
    return cur_v;
  }
  double Getnuc1_taxa_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc1_taxa_comphet_bool) {
      cur_v = codondata->nuc1_taxa_comphet();
      nuc1_taxa_comphet_bool = true;
    }
    return cur_v;
  }
  double Getnuc2_taxa_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc2_taxa_comphet_bool) {
      cur_v = codondata->nuc2_taxa_comphet();
      nuc2_taxa_comphet_bool = true;
    }
    return cur_v;
  }
  double Getnuc3_taxa_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!nuc3_taxa_comphet_bool) {
      cur_v = codondata->nuc3_taxa_comphet();
      nuc3_taxa_comphet_bool = true;
    }
    return cur_v;
  }
  double Getcodon_site_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!codon_site_comphet_bool) {
      cur_v = codondata->codon_site_comphet();
      codon_site_comphet_bool = true;
    }
    return cur_v;
  }
  double Getcodon_taxa_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!codon_taxa_comphet_bool) {
      cur_v = codondata->codon_taxa_comphet();
      codon_taxa_comphet_bool = true;
    }
    return cur_v;
  }
  double Getaa_site_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!aa_site_comphet_bool) {
      cur_v = codondata->aa_site_comphet();
      aa_site_comphet_bool = true;
    }
    return cur_v;
  }
  double Getaa_taxa_comphet(CodonSequenceAlignment* codondata) {
    double cur_v = 0.0;
    if (!aa_taxa_comphet_bool) {
      cur_v = codondata->aa_taxa_comphet();
      aa_taxa_comphet_bool = true;
    }
    return cur_v;
  }

  double GetDinucentropy(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    double H = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        H -= dinuc_usage[i][j] * log2(dinuc_usage[i][j]);
      }
    }
    return H;
  }

  double GetDinuc12entropy(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    double H = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        H -= dinuc12_usage[i][j] * log2(dinuc12_usage[i][j]);
      }
    }
    return H;
  }

  double GetDinuc23entropy(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    double H = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        H -= dinuc23_usage[i][j] * log2(dinuc23_usage[i][j]);
      }
    }
    return H;
  }

  double GetDinuc31entropy(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    double H = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        H -= dinuc31_usage[i][j] * log2(dinuc31_usage[i][j]);
      }
    }
    return H;
  }

  double GetRSCUentropy(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    double H = 0.0;
    for (int i = 0; i < codondata->GetCodonStateSpace()->GetNstate(); i++) {
      H -= RSCU[i] * log2(RSCU[i]);
    }
    return H;
  }

  double GetRelativeCodonFrequencyEntropy(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    double H = 0.0;
    for (int i = 0; i < codondata->GetCodonStateSpace()->GetNstate(); i++) {
      H -= relativeCodonFrequency[i] * log2(relativeCodonFrequency[i]);
    }
    return H;
  }

  double GetAAentropy(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    double H = 0.0;
    for (int i = 0; i < Naa; i++) {
      H -= relativeAAFrequency[i] * log2(relativeAAFrequency[i]);
    }
    return H;
  }

  double GetCGNonAGR(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    double CGN = RSCU[codondata->GetCodonStateSpace()->GetState("CGA")] +
                 RSCU[codondata->GetCodonStateSpace()->GetState("CGC")] +
                 RSCU[codondata->GetCodonStateSpace()->GetState("CGG")] +
                 RSCU[codondata->GetCodonStateSpace()->GetState("CGT")];
    double AGR = RSCU[codondata->GetCodonStateSpace()->GetState("AGA")] +
                 RSCU[codondata->GetCodonStateSpace()->GetState("AGG")];

    return CGN / AGR;
  }

  double GetCHQWonR(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }

    int C = 1;
    int H = 6;
    int Q = 13;
    int W = 18;
    int R = 14;

    double CHQW = relativeAAFrequency[C] + relativeAAFrequency[H] +
                  relativeAAFrequency[Q] + relativeAAFrequency[W];

    return CHQW / relativeAAFrequency[R];
  }

  double GetLMVonAPST(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }

    int A = 0;
    int L = 9;
    int M = 10;
    int P = 12;
    int S = 15;
    int T = 16;
    int V = 17;

    double MLV = relativeAAFrequency[M] + relativeAAFrequency[L] +
                 relativeAAFrequency[V];

    double APST = relativeAAFrequency[A] + relativeAAFrequency[P] +
                  relativeAAFrequency[T] + relativeAAFrequency[S];

    return MLV / APST;
  }

  /////////////////
  // codon_usage
  /////////////////

  double GetGGG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GGG")]);
  }
  double GetGGA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GGA")]);
  }
  double GetGGC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GGC")]);
  }
  double GetGGT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GGT")]);
  }
  double GetGAG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GAG")]);
  }
  double GetGAA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GAA")]);
  }
  double GetGAC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GAC")]);
  }
  double GetGAT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GAT")]);
  }
  double GetGCG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GCG")]);
  }
  double GetGCA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GCA")]);
  }
  double GetGCC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GCC")]);
  }
  double GetGCT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GCT")]);
  }
  double GetGTG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GTG")]);
  }
  double GetGTA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GTA")]);
  }
  double GetGTC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GTC")]);
  }
  double GetGTT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("GTT")]);
  }
  double GetAGG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AGG")]);
  }
  double GetAGA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AGA")]);
  }
  double GetAGC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AGC")]);
  }
  double GetAGT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AGT")]);
  }
  double GetAAG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AAG")]);
  }
  double GetAAA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AAA")]);
  }
  double GetAAC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AAC")]);
  }
  double GetAAT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("AAT")]);
  }
  double GetACG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ACG")]);
  }
  double GetACA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ACA")]);
  }
  double GetACC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ACC")]);
  }
  double GetACT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ACT")]);
  }
  double GetATG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ATG")]);
  }
  double GetATA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ATA")]);
  }
  double GetATC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ATC")]);
  }
  double GetATT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("ATT")]);
  }
  double GetCGG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CGG")]);
  }
  double GetCGA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CGA")]);
  }
  double GetCGC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CGC")]);
  }
  double GetCGT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CGT")]);
  }
  double GetCAG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CAG")]);
  }
  double GetCAA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CAA")]);
  }
  double GetCAC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CAC")]);
  }
  double GetCAT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CAT")]);
  }
  double GetCCG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CCG")]);
  }
  double GetCCA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CCA")]);
  }
  double GetCCC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CCC")]);
  }
  double GetCCT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CCT")]);
  }
  double GetCTG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CTG")]);
  }
  double GetCTA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CTA")]);
  }
  double GetCTC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CTC")]);
  }
  double GetCTT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("CTT")]);
  }
  double GetTGG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TGG")]);
  }
  double GetTGA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return 0.0;
  }
  double GetTGC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TGC")]);
  }
  double GetTGT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TGT")]);
  }
  double GetTAG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return 0.0;
  }

  double GetTAA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return 0.0;
  }

  double GetTAC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TAC")]);
  }

  double GetTAT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TAT")]);
  }
  double GetTCG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TCG")]);
  }

  double GetTCA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TCA")]);
  }
  double GetTCC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TCC")]);
  }
  double GetTCT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TCT")]);
  }

  double GetTTG(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TTG")]);
  }

  double GetTTA(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TTA")]);
  }
  double GetTTC(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TTC")]);
  }
  double GetTTT(CodonSequenceAlignment* codondata) {
    if (!RSCU_bool) {
      codondata->RSCU(RSCU);
      RSCU_bool = true;
    }
    return static_cast<double>(
        RSCU[codondata->GetCodonStateSpace()->GetState("TTT")]);
  }

  double GetfGGG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GGG")]);
  }
  double GetfGGA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GGA")]);
  }
  double GetfGGC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GGC")]);
  }
  double GetfGGT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GGT")]);
  }
  double GetfGAG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GAG")]);
  }
  double GetfGAA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GAA")]);
  }
  double GetfGAC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GAC")]);
  }
  double GetfGAT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GAT")]);
  }
  double GetfGCG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GCG")]);
  }
  double GetfGCA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GCA")]);
  }
  double GetfGCC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GCC")]);
  }
  double GetfGCT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GCT")]);
  }
  double GetfGTG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GTG")]);
  }
  double GetfGTA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GTA")]);
  }
  double GetfGTC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GTC")]);
  }
  double GetfGTT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "GTT")]);
  }
  double GetfAGG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AGG")]);
  }
  double GetfAGA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AGA")]);
  }
  double GetfAGC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AGC")]);
  }
  double GetfAGT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AGT")]);
  }
  double GetfAAG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AAG")]);
  }
  double GetfAAA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AAA")]);
  }
  double GetfAAC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AAC")]);
  }
  double GetfAAT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "AAT")]);
  }
  double GetfACG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ACG")]);
  }
  double GetfACA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ACA")]);
  }
  double GetfACC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ACC")]);
  }
  double GetfACT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ACT")]);
  }
  double GetfATG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ATG")]);
  }
  double GetfATA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ATA")]);
  }
  double GetfATC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ATC")]);
  }
  double GetfATT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "ATT")]);
  }
  double GetfCGG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CGG")]);
  }
  double GetfCGA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CGA")]);
  }
  double GetfCGC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CGC")]);
  }
  double GetfCGT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CGT")]);
  }
  double GetfCAG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CAG")]);
  }
  double GetfCAA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CAA")]);
  }
  double GetfCAC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CAC")]);
  }
  double GetfCAT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CAT")]);
  }
  double GetfCCG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CCG")]);
  }
  double GetfCCA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CCA")]);
  }
  double GetfCCC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CCC")]);
  }
  double GetfCCT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CCT")]);
  }
  double GetfCTG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CTG")]);
  }
  double GetfCTA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CTA")]);
  }
  double GetfCTC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CTC")]);
  }
  double GetfCTT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "CTT")]);
  }
  double GetfTGG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TGG")]);
  }
  double GetfTGA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return 0.0;
  }
  double GetfTGC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TGC")]);
  }
  double GetfTGT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TGT")]);
  }
  double GetfTAG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return 0.0;
  }

  double GetfTAA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return 0.0;
  }

  double GetfTAC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TAC")]);
  }

  double GetfTAT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TAT")]);
  }
  double GetfTCG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TCG")]);
  }

  double GetfTCA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TCA")]);
  }
  double GetfTCC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TCC")]);
  }
  double GetfTCT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TCT")]);
  }

  double GetfTTG(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TTG")]);
  }

  double GetfTTA(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TTA")]);
  }
  double GetfTTC(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TTC")]);
  }
  double GetfTTT(CodonSequenceAlignment* codondata) {
    if (!relativeCodonFrequency_bool) {
      codondata->relativeCodonFrequency(relativeCodonFrequency);
      relativeCodonFrequency_bool = true;
    }
    return static_cast<double>(
        relativeCodonFrequency[codondata->GetCodonStateSpace()->GetState(
            "TTT")]);
  }

  /////////////////
  // aa_usage
  /////////////////

  double GetA(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[0]);
  }
  double GetC(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[1]);
  }
  double GetD(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[2]);
  }
  double GetE(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[3]);
  }
  double GetF(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[4]);
  }
  double GetG(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[5]);
  }
  double GetH(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[6]);
  }
  double GetI(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[7]);
  }
  double GetK(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[8]);
  }
  double GetL(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[9]);
  }
  double GetM(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[10]);
  }
  double GetN(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[11]);
  }
  double GetP(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[12]);
  }
  double GetQ(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[13]);
  }
  double GetR(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[14]);
  }
  double GetS(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[15]);
  }
  double GetT(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[16]);
  }
  double GetV(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[17]);
  }
  double GetW(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[18]);
  }
  double GetY(CodonSequenceAlignment* codondata) {
    if (!relativeAAFrequency_bool) {
      codondata->relativeAAFrequency(relativeAAFrequency);
      relativeAAFrequency_bool = true;
    }
    return static_cast<double>(relativeAAFrequency[19]);
  }
  /////////////////
  // dinuc31
  /////////////////
  double GetDinuc31AA(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[0][0]);
  }

  double GetDinuc31AC(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[0][1]);
  }

  double GetDinuc31AG(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[0][2]);
  }
  double GetDinuc31AT(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[0][3]);
  }
  double GetDinuc31CA(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[1][0]);
  }

  double GetDinuc31CC(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[1][1]);
  }
  double GetDinuc31CG(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[1][2]);
  }
  double GetDinuc31CT(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[1][3]);
  }
  double GetDinuc31GA(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[2][0]);
  }

  double GetDinuc31GC(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[2][1]);
  }
  double GetDinuc31GG(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[2][2]);
  }
  double GetDinuc31GT(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[2][3]);
  }
  double GetDinuc31TA(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[3][0]);
  }

  double GetDinuc31TC(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[3][1]);
  }
  double GetDinuc31TG(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[3][2]);
  }
  double GetDinuc31TT(CodonSequenceAlignment* codondata) {
    if (!dinuc31_bool) {
      codondata->dinuc31_usage(dinuc31_usage);
      dinuc31_bool = true;
    }
    return static_cast<double>(dinuc31_usage[3][3]);
  }
  /////////////////
  // dinuc23
  /////////////////
  double GetDinuc23AA(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[0][0]);
  }

  double GetDinuc23AC(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[0][1]);
  }
  double GetDinuc23AG(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[0][2]);
  }
  double GetDinuc23AT(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[0][3]);
  }
  double GetDinuc23CA(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[1][0]);
  }

  double GetDinuc23CC(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[1][1]);
  }
  double GetDinuc23CG(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[1][2]);
  }
  double GetDinuc23CT(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[1][3]);
  }
  double GetDinuc23GA(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[2][0]);
  }

  double GetDinuc23GC(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[2][1]);
  }
  double GetDinuc23GG(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[2][2]);
  }
  double GetDinuc23GT(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[2][3]);
  }
  double GetDinuc23TA(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[3][0]);
  }

  double GetDinuc23TC(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[3][1]);
  }
  double GetDinuc23TG(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[3][2]);
  }
  double GetDinuc23TT(CodonSequenceAlignment* codondata) {
    if (!dinuc23_bool) {
      codondata->dinuc23_usage(dinuc23_usage);
      dinuc23_bool = true;
    }
    return static_cast<double>(dinuc23_usage[3][3]);
  }

  /////////////////
  // dinuc12
  /////////////////
  double GetDinuc12AA(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[0][0]);
  }

  double GetDinuc12AC(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[0][1]);
  }
  double GetDinuc12AG(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[0][2]);
  }
  double GetDinuc12AT(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[0][3]);
  }
  double GetDinuc12CA(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[1][0]);
  }

  double GetDinuc12CC(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[1][1]);
  }
  double GetDinuc12CG(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[1][2]);
  }
  double GetDinuc12CT(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[1][3]);
  }
  double GetDinuc12GA(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[2][0]);
  }

  double GetDinuc12GC(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[2][1]);
  }
  double GetDinuc12GG(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[2][2]);
  }
  double GetDinuc12GT(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[2][3]);
  }
  double GetDinuc12TA(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[3][0]);
  }

  double GetDinuc12TC(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[3][1]);
  }
  double GetDinuc12TG(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[3][2]);
  }
  double GetDinuc12TT(CodonSequenceAlignment* codondata) {
    if (!dinuc12_bool) {
      codondata->dinuc12_usage(dinuc12_usage);
      dinuc12_bool = true;
    }
    return static_cast<double>(dinuc12_usage[3][3]);
  }
  /////////////////
  // dinuc
  /////////////////
  double GetDinucAA(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[0][0];
  }

  double GetDinucAC(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[0][1];
  }
  double GetDinucAG(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[0][2];
  }
  double GetDinucAT(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[0][3];
  }
  double GetDinucCA(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[1][0];
  }

  double GetDinucCC(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[1][1];
  }
  double GetDinucCG(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[1][2];
  }
  double GetDinucCT(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[1][3];
  }
  double GetDinucGA(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[2][0];
  }

  double GetDinucGC(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[2][1];
  }
  double GetDinucGG(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[2][2];
  }
  double GetDinucGT(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[2][3];
  }
  double GetDinucTA(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[3][0];
  }

  double GetDinucTC(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[3][1];
  }
  double GetDinucTG(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[3][2];
  }
  double GetDinucTT(CodonSequenceAlignment* codondata) {
    if (!dinuc_bool) {
      codondata->dinuc_usage(dinuc_usage);
      dinuc_bool = true;
    }
    return dinuc_usage[3][3];
  }

  /////////////////
  // nuc1
  /////////////////
  double GetNuc1A(CodonSequenceAlignment* codondata) {
    if (!nuc1_bool) {
      codondata->nuc1_usage(nuc1_usage);
      nuc1_bool = true;
    }
    return nuc1_usage[0];
  }

  double GetNuc1C(CodonSequenceAlignment* codondata) {
    if (!nuc1_bool) {
      codondata->nuc1_usage(nuc1_usage);
      nuc1_bool = true;
    }
    return nuc1_usage[1];
  }

  double GetNuc1G(CodonSequenceAlignment* codondata) {
    if (!nuc1_bool) {
      codondata->nuc1_usage(nuc1_usage);
      nuc1_bool = true;
    }
    return nuc1_usage[2];
  }

  double GetNuc1T(CodonSequenceAlignment* codondata) {
    if (!nuc1_bool) {
      codondata->nuc1_usage(nuc1_usage);
      nuc1_bool = true;
    }
    return nuc1_usage[3];
  }
  /////////////////
  // nuc2
  /////////////////
  double GetNuc2A(CodonSequenceAlignment* codondata) {
    if (!nuc2_bool) {
      codondata->nuc2_usage(nuc2_usage);
      nuc2_bool = true;
    }
    return nuc2_usage[0];
  }

  double GetNuc2C(CodonSequenceAlignment* codondata) {
    if (!nuc2_bool) {
      codondata->nuc2_usage(nuc2_usage);
      nuc2_bool = true;
    }
    return nuc2_usage[1];
  }

  double GetNuc2G(CodonSequenceAlignment* codondata) {
    if (!nuc2_bool) {
      codondata->nuc2_usage(nuc2_usage);
      nuc2_bool = true;
    }
    return nuc2_usage[2];
  }

  double GetNuc2T(CodonSequenceAlignment* codondata) {
    if (!nuc2_bool) {
      codondata->nuc2_usage(nuc2_usage);
      nuc2_bool = true;
    }
    return nuc2_usage[3];
  }
  /////////////////
  // GC
  /////////////////
  double GetGC(CodonSequenceAlignment* codondata) {
    if (!GC_bool) {
      codondata->nuc_usage(nuc_usage);
      nuc_bool = true;
    }

    return nuc_usage[1] + nuc_usage[2];
  }

  double GetGC1(CodonSequenceAlignment* codondata) {
    if (!GC1_bool) {
      codondata->nuc1_usage(nuc1_usage);
      nuc1_bool = true;
    }

    return nuc1_usage[1] + nuc1_usage[2];
  }

  double GetGC2(CodonSequenceAlignment* codondata) {
    if (!GC2_bool) {
      codondata->nuc2_usage(nuc2_usage);
      nuc2_bool = true;
    }

    return nuc2_usage[1] + nuc2_usage[2];
  }

  double GetGC3(CodonSequenceAlignment* codondata) {
    if (!GC3_bool) {
      codondata->nuc3_usage(nuc3_usage);
      nuc3_bool = true;
    }

    return nuc3_usage[1] + nuc3_usage[2];
  }

  /////////////////
  // nuc3
  /////////////////
  double GetNuc3A(CodonSequenceAlignment* codondata) {
    if (!nuc3_bool) {
      codondata->nuc3_usage(nuc3_usage);
      nuc3_bool = true;
    }
    return nuc3_usage[0];
  }

  double GetNuc3C(CodonSequenceAlignment* codondata) {
    if (!nuc3_bool) {
      codondata->nuc3_usage(nuc3_usage);
      nuc3_bool = true;
    }
    return nuc3_usage[1];
  }

  double GetNuc3G(CodonSequenceAlignment* codondata) {
    if (!nuc3_bool) {
      codondata->nuc3_usage(nuc3_usage);
      nuc3_bool = true;
    }
    return nuc3_usage[2];
  }

  double GetNuc3T(CodonSequenceAlignment* codondata) {
    if (!nuc3_bool) {
      codondata->nuc3_usage(nuc3_usage);
      nuc3_bool = true;
    }
    return nuc3_usage[3];
  }
  /////////////////
  // nuc
  /////////////////
  double GetNucA(CodonSequenceAlignment* codondata) {
    if (!nuc_bool) {
      codondata->nuc_usage(nuc_usage);
      nuc_bool = true;
    }
    return nuc_usage[0];
  }

  double GetNucC(CodonSequenceAlignment* codondata) {
    if (!nuc_bool) {
      codondata->nuc_usage(nuc_usage);
      nuc_bool = true;
    }
    return nuc_usage[1];
  }

  double GetNucG(CodonSequenceAlignment* codondata) {
    if (!nuc_bool) {
      codondata->nuc_usage(nuc_usage);
      nuc_bool = true;
    }
    return nuc_usage[2];
  }

  double GetNucT(CodonSequenceAlignment* codondata) {
    if (!nuc_bool) {
      codondata->nuc_usage(nuc_usage);
      nuc_bool = true;
    }
    return nuc_usage[3];
  }

  /////////////////
  // nuc3_meandiff
  /////////////////
  double GetNuc3mean(CodonSequenceAlignment* codondata) {
    if (!nuc3_meandiff_bool) {
      codondata->nuc3_meandiff(nuc3_meandiff);
      nuc3_meandiff_bool = true;
    }
    return nuc3_meandiff[0];
  }

  double GetNuc3var(CodonSequenceAlignment* codondata) {
    if (!nuc3_meandiff_bool) {
      codondata->nuc3_meandiff(nuc3_meandiff);
      nuc3_meandiff_bool = true;
    }
    return nuc3_meandiff[3];
  }
  /////////////////
  // nuc2_meandiff
  /////////////////
  double GetNuc2mean(CodonSequenceAlignment* codondata) {
    if (!nuc2_meandiff_bool) {
      codondata->nuc2_meandiff(nuc2_meandiff);
      nuc2_meandiff_bool = true;
    }
    return nuc2_meandiff[0];
  }

  double GetNuc2var(CodonSequenceAlignment* codondata) {
    if (!nuc2_meandiff_bool) {
      codondata->nuc2_meandiff(nuc2_meandiff);
      nuc2_meandiff_bool = true;
    }
    return nuc2_meandiff[2];
  }
  /////////////////
  // nuc1_meandiff
  /////////////////
  double GetNuc1mean(CodonSequenceAlignment* codondata) {
    if (!nuc1_meandiff_bool) {
      codondata->nuc1_meandiff(nuc1_meandiff);
      nuc1_meandiff_bool = true;
    }
    return nuc1_meandiff[0];
  }

  double GetNuc1var(CodonSequenceAlignment* codondata) {
    if (!nuc1_meandiff_bool) {
      codondata->nuc1_meandiff(nuc1_meandiff);
      nuc1_meandiff_bool = true;
    }
    return nuc1_meandiff[1];
  }
  /////////////////
  // nuc_meandiff
  /////////////////
  double GetNucmean(CodonSequenceAlignment* codondata) {
    if (!nuc_meandiff_bool) {
      codondata->nuc_meandiff(nuc_meandiff);
      nuc_meandiff_bool = true;
    }
    return nuc_meandiff[0];
  }

  double GetNucvar(CodonSequenceAlignment* codondata) {
    if (!nuc_meandiff_bool) {
      codondata->nuc_meandiff(nuc_meandiff);
      nuc_meandiff_bool = true;
    }
    return nuc_meandiff[1];
  }

  /////////////////
  // codon_meandiff
  /////////////////
  double GetCodonmean(CodonSequenceAlignment* codondata) {
    if (!codon_meandiff_bool) {
      codondata->codon_meandiff(codon_meandiff);
      codon_meandiff_bool = true;
    }
    return codon_meandiff[0];
  }

  double GetCodonvar(CodonSequenceAlignment* codondata) {
    if (!codon_meandiff_bool) {
      codondata->codon_meandiff(codon_meandiff);
      codon_meandiff_bool = true;
    }
    return codon_meandiff[1];
  }
  /////////////////
  // aa_meandiff
  /////////////////
  double GetAAmean(CodonSequenceAlignment* codondata) {
    if (!aa_meandiff_bool) {
      codondata->aa_meandiff(aa_meandiff);
      aa_meandiff_bool = true;
    }
    return aa_meandiff[0];
  }

  double GetAAvar(CodonSequenceAlignment* codondata) {
    if (!aa_meandiff_bool) {
      codondata->aa_meandiff(aa_meandiff);
      aa_meandiff_bool = true;
    }
    return aa_meandiff[1];
  }

  /////////////////
  // dinucCpG_pairwise
  /////////////////
  double GetdinucCpG_TpG(CodonSequenceAlignment* codondata) {
    if (!dinucCpG_pairwise_bool) {
      codondata->dinucCpG_pairwise(dinucCpG_pairwise);
      dinucCpG_pairwise_bool = true;
    }
    return static_cast<double>(dinucCpG_pairwise[0]);
  }

  double GetdinucCpG_CpA(CodonSequenceAlignment* codondata) {
    if (!dinucCpG_pairwise_bool) {
      codondata->dinucCpG_pairwise(dinucCpG_pairwise);
      dinucCpG_pairwise_bool = true;
    }
    return static_cast<double>(dinucCpG_pairwise[1]);
  }

  double GetdinucApG_TpG(CodonSequenceAlignment* codondata) {
    if (!dinucCpG_pairwise_bool) {
      codondata->dinucCpG_pairwise(dinucCpG_pairwise);
      dinucCpG_pairwise_bool = true;
    }
    return static_cast<double>(dinucCpG_pairwise[2]);
  }

  /////////////////
  // aa_pairwise
  /////////////////
  double GetpwAA(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool) {
      codondata->aa_pairwise(aa_pairwise);
      aa_pairwise_bool = true;
    }

    return aa_pairwise[0];
  }

  double GetpwAA10(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool10) {
      codondata->aa_pairwise10(aa_pairwise10);
      aa_pairwise_bool10 = true;
    }

    return aa_pairwise10[0];
  }

  double GetpwAA30(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool30) {
      codondata->aa_pairwise30(aa_pairwise30);
      aa_pairwise_bool30 = true;
    }

    return aa_pairwise30[0];
  }

  double GetpwAA50(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool50) {
      codondata->aa_pairwise50(aa_pairwise50);
      aa_pairwise_bool50 = true;
    }

    return aa_pairwise50[0];
  }

  double GetpwAA70(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool70) {
      codondata->aa_pairwise70(aa_pairwise70);
      aa_pairwise_bool70 = true;
    }

    return aa_pairwise70[0];
  }

  double GetpwAA90(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool90) {
      codondata->aa_pairwise90(aa_pairwise90);
      aa_pairwise_bool90 = true;
    }

    return aa_pairwise90[0];
  }

  double GetK80aa(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool) {
      codondata->aa_pairwise(aa_pairwise);
      aa_pairwise_bool = true;
    }

    return aa_pairwise[1];
  }

  double GetK80aa10(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool10) {
      codondata->aa_pairwise10(aa_pairwise10);
      aa_pairwise_bool10 = true;
    }

    return aa_pairwise10[1];
  }

  double GetK80aa30(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool30) {
      codondata->aa_pairwise30(aa_pairwise30);
      aa_pairwise_bool30 = true;
    }

    return aa_pairwise30[1];
  }

  double GetK80aa50(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool50) {
      codondata->aa_pairwise50(aa_pairwise50);
      aa_pairwise_bool50 = true;
    }

    return aa_pairwise50[1];
  }

  double GetK80aa70(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool70) {
      codondata->aa_pairwise70(aa_pairwise70);
      aa_pairwise_bool70 = true;
    }

    return aa_pairwise70[1];
  }

  double GetK80aa90(CodonSequenceAlignment* codondata) {
    if (!aa_pairwise_bool90) {
      codondata->aa_pairwise90(aa_pairwise90);
      aa_pairwise_bool90 = true;
    }

    return aa_pairwise90[1];
  }

  /////////////////
  // nuc3_pairwise
  /////////////////

  double Getpw3GT(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool) {
      codondata->nuc3_pairwise(nuc3_pairwise);
      nuc3_pairwise_bool = true;
    }
    return static_cast<double>(nuc3_pairwise[5]);
  }

  double Getpw3CT(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool) {
      codondata->nuc3_pairwise(nuc3_pairwise);
      nuc3_pairwise_bool = true;
    }
    return static_cast<double>(nuc3_pairwise[4]);
  }

  double Getpw3CG(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool) {
      codondata->nuc3_pairwise(nuc3_pairwise);
      nuc3_pairwise_bool = true;
    }
    return static_cast<double>(nuc3_pairwise[3]);
  }

  double Getpw3AT(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool) {
      codondata->nuc3_pairwise(nuc3_pairwise);
      nuc3_pairwise_bool = true;
    }
    return static_cast<double>(nuc3_pairwise[2]);
  }

  double Getpw3AG(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool) {
      codondata->nuc3_pairwise(nuc3_pairwise);
      nuc3_pairwise_bool = true;
    }
    return static_cast<double>(nuc3_pairwise[1]);
  }

  double Getpw3AC(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool) {
      codondata->nuc3_pairwise(nuc3_pairwise);
      nuc3_pairwise_bool = true;
    }
    return static_cast<double>(nuc3_pairwise[0]);
  }

  double GetK80nuc3(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool) {
      codondata->nuc3_pairwise10(nuc3_pairwise);
      nuc3_pairwise_bool = true;
    }
    return nuc3_pairwise[6];
  }

  double GetK80nuc310(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool10) {
      codondata->nuc3_pairwise10(nuc3_pairwise10);
      nuc3_pairwise_bool10 = true;
    }
    return nuc3_pairwise10[6];
  }

  double GetK80nuc330(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool30) {
      codondata->nuc3_pairwise30(nuc3_pairwise30);
      nuc3_pairwise_bool30 = true;
    }
    return nuc3_pairwise30[6];
  }

  double GetK80nuc350(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool50) {
      codondata->nuc3_pairwise50(nuc3_pairwise50);
      nuc3_pairwise_bool50 = true;
    }
    return nuc3_pairwise50[6];
  }

  double GetK80nuc370(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool70) {
      codondata->nuc3_pairwise70(nuc3_pairwise70);
      nuc3_pairwise_bool70 = true;
    }
    return nuc3_pairwise70[6];
  }

  double GetK80nuc390(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool90) {
      codondata->nuc3_pairwise90(nuc3_pairwise90);
      nuc3_pairwise_bool90 = true;
    }
    return nuc3_pairwise90[6];
  }

  // 10
  double Getpw3GT10(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool10) {
      codondata->nuc3_pairwise10(nuc3_pairwise10);
      nuc3_pairwise_bool10 = true;
    }
    return nuc3_pairwise10[5];
  }

  double Getpw3CT10(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool10) {
      codondata->nuc3_pairwise10(nuc3_pairwise10);
      nuc3_pairwise_bool10 = true;
    }
    return nuc3_pairwise10[4];
  }

  double Getpw3CG10(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool10) {
      codondata->nuc3_pairwise10(nuc3_pairwise10);
      nuc3_pairwise_bool10 = true;
    }
    return nuc3_pairwise10[3];
  }

  double Getpw3AT10(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool10) {
      codondata->nuc3_pairwise10(nuc3_pairwise10);
      nuc3_pairwise_bool10 = true;
    }
    return nuc3_pairwise10[2];
  }

  double Getpw3AG10(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool10) {
      codondata->nuc3_pairwise10(nuc3_pairwise10);
      nuc3_pairwise_bool10 = true;
    }
    return nuc3_pairwise10[1];
  }

  double Getpw3AC10(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool10) {
      codondata->nuc3_pairwise10(nuc3_pairwise10);
      nuc3_pairwise_bool10 = true;
    }
    return nuc3_pairwise10[0];
  }

  // 30
  double Getpw3GT30(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool30) {
      codondata->nuc3_pairwise30(nuc3_pairwise30);
      nuc3_pairwise_bool30 = true;
    }
    return nuc3_pairwise30[5];
  }

  double Getpw3CT30(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool30) {
      codondata->nuc3_pairwise30(nuc3_pairwise30);
      nuc3_pairwise_bool30 = true;
    }
    return nuc3_pairwise30[4];
  }

  double Getpw3CG30(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool30) {
      codondata->nuc3_pairwise30(nuc3_pairwise30);
      nuc3_pairwise_bool30 = true;
    }
    return nuc3_pairwise30[3];
  }

  double Getpw3AT30(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool30) {
      codondata->nuc3_pairwise30(nuc3_pairwise30);
      nuc3_pairwise_bool30 = true;
    }
    return nuc3_pairwise30[2];
  }

  double Getpw3AG30(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool30) {
      codondata->nuc3_pairwise30(nuc3_pairwise30);
      nuc3_pairwise_bool30 = true;
    }
    return nuc3_pairwise30[1];
  }

  double Getpw3AC30(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool30) {
      codondata->nuc3_pairwise30(nuc3_pairwise30);
      nuc3_pairwise_bool30 = true;
    }
    return nuc3_pairwise30[0];
  }
  // 50

  double Getpw3GT50(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool50) {
      codondata->nuc3_pairwise50(nuc3_pairwise50);
      nuc3_pairwise_bool50 = true;
    }
    return nuc3_pairwise50[5];
  }

  double Getpw3CT50(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool50) {
      codondata->nuc3_pairwise50(nuc3_pairwise50);
      nuc3_pairwise_bool50 = true;
    }
    return nuc3_pairwise50[4];
  }

  double Getpw3CG50(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool50) {
      codondata->nuc3_pairwise50(nuc3_pairwise50);
      nuc3_pairwise_bool50 = true;
    }
    return nuc3_pairwise50[3];
  }

  double Getpw3AT50(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool50) {
      codondata->nuc3_pairwise50(nuc3_pairwise50);
      nuc3_pairwise_bool50 = true;
    }
    return nuc3_pairwise50[2];
  }

  double Getpw3AG50(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool50) {
      codondata->nuc3_pairwise50(nuc3_pairwise50);
      nuc3_pairwise_bool50 = true;
    }
    return nuc3_pairwise50[1];
  }

  double Getpw3AC50(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool50) {
      codondata->nuc3_pairwise50(nuc3_pairwise50);
      nuc3_pairwise_bool50 = true;
    }
    return nuc3_pairwise50[0];
  }

  // 70

  double Getpw3GT70(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool70) {
      codondata->nuc3_pairwise70(nuc3_pairwise70);
      nuc3_pairwise_bool70 = true;
    }
    return nuc3_pairwise70[5];
  }

  double Getpw3CT70(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool70) {
      codondata->nuc3_pairwise70(nuc3_pairwise70);
      nuc3_pairwise_bool70 = true;
    }
    return nuc3_pairwise70[4];
  }

  double Getpw3CG70(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool70) {
      codondata->nuc3_pairwise70(nuc3_pairwise70);
      nuc3_pairwise_bool70 = true;
    }
    return nuc3_pairwise70[3];
  }

  double Getpw3AT70(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool70) {
      codondata->nuc3_pairwise70(nuc3_pairwise70);
      nuc3_pairwise_bool70 = true;
    }
    return nuc3_pairwise70[2];
  }

  double Getpw3AG70(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool70) {
      codondata->nuc3_pairwise70(nuc3_pairwise70);
      nuc3_pairwise_bool70 = true;
    }
    return nuc3_pairwise70[1];
  }

  double Getpw3AC70(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool70) {
      codondata->nuc3_pairwise70(nuc3_pairwise70);
      nuc3_pairwise_bool70 = true;
    }
    return nuc3_pairwise70[0];
  }

  // 90
  double Getpw3GT90(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool90) {
      codondata->nuc3_pairwise90(nuc3_pairwise90);
      nuc3_pairwise_bool90 = true;
    }
    return nuc3_pairwise90[5];
  }

  double Getpw3CT90(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool90) {
      codondata->nuc3_pairwise90(nuc3_pairwise90);
      nuc3_pairwise_bool90 = true;
    }
    return nuc3_pairwise90[4];
  }

  double Getpw3CG90(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool90) {
      codondata->nuc3_pairwise90(nuc3_pairwise90);
      nuc3_pairwise_bool90 = true;
    }
    return nuc3_pairwise90[3];
  }

  double Getpw3AT90(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool90) {
      codondata->nuc3_pairwise90(nuc3_pairwise90);
      nuc3_pairwise_bool90 = true;
    }
    return nuc3_pairwise90[2];
  }

  double Getpw3AG90(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool90) {
      codondata->nuc3_pairwise90(nuc3_pairwise90);
      nuc3_pairwise_bool90 = true;
    }
    return nuc3_pairwise90[1];
  }

  double Getpw3AC90(CodonSequenceAlignment* codondata) {
    if (!nuc3_pairwise_bool90) {
      codondata->nuc3_pairwise90(nuc3_pairwise90);
      nuc3_pairwise_bool90 = true;
    }
    return nuc3_pairwise90[0];
  }

  /////////////////
  // nuc2_pairwise
  /////////////////

  double Getpw2GT(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool) {
      codondata->nuc2_pairwise(nuc2_pairwise);
      nuc2_pairwise_bool = true;
    }
    return static_cast<double>(nuc2_pairwise[5]);
  }

  double Getpw2CT(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool) {
      codondata->nuc2_pairwise(nuc2_pairwise);
      nuc2_pairwise_bool = true;
    }
    return static_cast<double>(nuc2_pairwise[4]);
  }

  double Getpw2CG(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool) {
      codondata->nuc2_pairwise(nuc2_pairwise);
      nuc2_pairwise_bool = true;
    }
    return static_cast<double>(nuc2_pairwise[3]);
  }

  double Getpw2AT(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool) {
      codondata->nuc2_pairwise(nuc2_pairwise);
      nuc2_pairwise_bool = true;
    }
    return static_cast<double>(nuc2_pairwise[2]);
  }

  double Getpw2AG(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool) {
      codondata->nuc2_pairwise(nuc2_pairwise);
      nuc2_pairwise_bool = true;
    }
    return static_cast<double>(nuc2_pairwise[1]);
  }

  double Getpw2AC(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool) {
      codondata->nuc2_pairwise(nuc2_pairwise);
      nuc2_pairwise_bool = true;
    }
    return static_cast<double>(nuc2_pairwise[0]);
  }

  double GetK80nuc2(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool) {
      codondata->nuc2_pairwise10(nuc2_pairwise);
      nuc2_pairwise_bool = true;
    }
    return nuc2_pairwise[6];
  }

  double GetK80nuc210(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool10) {
      codondata->nuc2_pairwise10(nuc2_pairwise10);
      nuc2_pairwise_bool10 = true;
    }
    return nuc2_pairwise10[6];
  }

  double GetK80nuc230(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool30) {
      codondata->nuc2_pairwise30(nuc2_pairwise30);
      nuc2_pairwise_bool30 = true;
    }
    return nuc2_pairwise30[6];
  }

  double GetK80nuc250(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool50) {
      codondata->nuc2_pairwise50(nuc2_pairwise50);
      nuc2_pairwise_bool50 = true;
    }
    return nuc2_pairwise50[6];
  }

  double GetK80nuc270(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool70) {
      codondata->nuc2_pairwise70(nuc2_pairwise70);
      nuc2_pairwise_bool70 = true;
    }
    return nuc2_pairwise70[6];
  }

  double GetK80nuc290(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool90) {
      codondata->nuc2_pairwise90(nuc2_pairwise90);
      nuc2_pairwise_bool90 = true;
    }
    return nuc2_pairwise90[6];
  }

  // 10
  double Getpw2GT10(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool10) {
      codondata->nuc2_pairwise10(nuc2_pairwise10);
      nuc2_pairwise_bool10 = true;
    }
    return nuc2_pairwise10[5];
  }

  double Getpw2CT10(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool10) {
      codondata->nuc2_pairwise10(nuc2_pairwise10);
      nuc2_pairwise_bool10 = true;
    }
    return nuc2_pairwise10[4];
  }

  double Getpw2CG10(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool10) {
      codondata->nuc2_pairwise10(nuc2_pairwise10);
      nuc2_pairwise_bool10 = true;
    }
    return nuc2_pairwise10[3];
  }

  double Getpw2AT10(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool10) {
      codondata->nuc2_pairwise10(nuc2_pairwise10);
      nuc2_pairwise_bool10 = true;
    }
    return nuc2_pairwise10[2];
  }

  double Getpw2AG10(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool10) {
      codondata->nuc2_pairwise10(nuc2_pairwise10);
      nuc2_pairwise_bool10 = true;
    }
    return nuc2_pairwise10[1];
  }

  double Getpw2AC10(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool10) {
      codondata->nuc2_pairwise10(nuc2_pairwise10);
      nuc2_pairwise_bool10 = true;
    }
    return nuc2_pairwise10[0];
  }

  // 30
  double Getpw2GT30(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool30) {
      codondata->nuc2_pairwise30(nuc2_pairwise30);
      nuc2_pairwise_bool30 = true;
    }
    return nuc2_pairwise30[5];
  }

  double Getpw2CT30(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool30) {
      codondata->nuc2_pairwise30(nuc2_pairwise30);
      nuc2_pairwise_bool30 = true;
    }
    return nuc2_pairwise30[4];
  }

  double Getpw2CG30(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool30) {
      codondata->nuc2_pairwise30(nuc2_pairwise30);
      nuc2_pairwise_bool30 = true;
    }
    return nuc2_pairwise30[3];
  }

  double Getpw2AT30(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool30) {
      codondata->nuc2_pairwise30(nuc2_pairwise30);
      nuc2_pairwise_bool30 = true;
    }
    return nuc2_pairwise30[2];
  }

  double Getpw2AG30(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool30) {
      codondata->nuc2_pairwise30(nuc2_pairwise30);
      nuc2_pairwise_bool30 = true;
    }
    return nuc2_pairwise30[1];
  }

  double Getpw2AC30(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool30) {
      codondata->nuc2_pairwise30(nuc2_pairwise30);
      nuc2_pairwise_bool30 = true;
    }
    return nuc2_pairwise30[0];
  }
  // 50

  double Getpw2GT50(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool50) {
      codondata->nuc2_pairwise50(nuc2_pairwise50);
      nuc2_pairwise_bool50 = true;
    }
    return nuc2_pairwise50[5];
  }

  double Getpw2CT50(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool50) {
      codondata->nuc2_pairwise50(nuc2_pairwise50);
      nuc2_pairwise_bool50 = true;
    }
    return nuc2_pairwise50[4];
  }

  double Getpw2CG50(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool50) {
      codondata->nuc2_pairwise50(nuc2_pairwise50);
      nuc2_pairwise_bool50 = true;
    }
    return nuc2_pairwise50[3];
  }

  double Getpw2AT50(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool50) {
      codondata->nuc2_pairwise50(nuc2_pairwise50);
      nuc2_pairwise_bool50 = true;
    }
    return nuc2_pairwise50[2];
  }

  double Getpw2AG50(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool50) {
      codondata->nuc2_pairwise50(nuc2_pairwise50);
      nuc2_pairwise_bool50 = true;
    }
    return nuc2_pairwise50[1];
  }

  double Getpw2AC50(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool50) {
      codondata->nuc2_pairwise50(nuc2_pairwise50);
      nuc2_pairwise_bool50 = true;
    }
    return nuc2_pairwise50[0];
  }

  // 70

  double Getpw2GT70(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool70) {
      codondata->nuc2_pairwise70(nuc2_pairwise70);
      nuc2_pairwise_bool70 = true;
    }
    return nuc2_pairwise70[5];
  }

  double Getpw2CT70(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool70) {
      codondata->nuc2_pairwise70(nuc2_pairwise70);
      nuc2_pairwise_bool70 = true;
    }
    return nuc2_pairwise70[4];
  }

  double Getpw2CG70(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool70) {
      codondata->nuc2_pairwise70(nuc2_pairwise70);
      nuc2_pairwise_bool70 = true;
    }
    return nuc2_pairwise70[3];
  }

  double Getpw2AT70(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool70) {
      codondata->nuc2_pairwise70(nuc2_pairwise70);
      nuc2_pairwise_bool70 = true;
    }
    return nuc2_pairwise70[2];
  }

  double Getpw2AG70(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool70) {
      codondata->nuc2_pairwise70(nuc2_pairwise70);
      nuc2_pairwise_bool70 = true;
    }
    return nuc2_pairwise70[1];
  }

  double Getpw2AC70(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool70) {
      codondata->nuc2_pairwise70(nuc2_pairwise70);
      nuc2_pairwise_bool70 = true;
    }
    return nuc2_pairwise70[0];
  }

  // 90
  double Getpw2GT90(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool90) {
      codondata->nuc2_pairwise90(nuc2_pairwise90);
      nuc2_pairwise_bool90 = true;
    }
    return nuc2_pairwise90[5];
  }

  double Getpw2CT90(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool90) {
      codondata->nuc2_pairwise90(nuc2_pairwise90);
      nuc2_pairwise_bool90 = true;
    }
    return nuc2_pairwise90[4];
  }

  double Getpw2CG90(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool90) {
      codondata->nuc2_pairwise90(nuc2_pairwise90);
      nuc2_pairwise_bool90 = true;
    }
    return nuc2_pairwise90[3];
  }

  double Getpw2AT90(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool90) {
      codondata->nuc2_pairwise90(nuc2_pairwise90);
      nuc2_pairwise_bool90 = true;
    }
    return nuc2_pairwise90[2];
  }

  double Getpw2AG90(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool90) {
      codondata->nuc2_pairwise90(nuc2_pairwise90);
      nuc2_pairwise_bool90 = true;
    }
    return nuc2_pairwise90[1];
  }

  double Getpw2AC90(CodonSequenceAlignment* codondata) {
    if (!nuc2_pairwise_bool90) {
      codondata->nuc2_pairwise90(nuc2_pairwise90);
      nuc2_pairwise_bool90 = true;
    }
    return nuc2_pairwise90[0];
  }

  /////////////////
  // nuc1_pairwise
  /////////////////

  double Getpw1GT(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool) {
      codondata->nuc1_pairwise(nuc1_pairwise);
      nuc1_pairwise_bool = true;
    }
    return static_cast<double>(nuc1_pairwise[5]);
  }

  double Getpw1CT(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool) {
      codondata->nuc1_pairwise(nuc1_pairwise);
      nuc1_pairwise_bool = true;
    }
    return static_cast<double>(nuc1_pairwise[4]);
  }

  double Getpw1CG(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool) {
      codondata->nuc1_pairwise(nuc1_pairwise);
      nuc1_pairwise_bool = true;
    }
    return static_cast<double>(nuc1_pairwise[3]);
  }

  double Getpw1AT(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool) {
      codondata->nuc1_pairwise(nuc1_pairwise);
      nuc1_pairwise_bool = true;
    }
    return static_cast<double>(nuc1_pairwise[2]);
  }

  double Getpw1AG(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool) {
      codondata->nuc1_pairwise(nuc1_pairwise);
      nuc1_pairwise_bool = true;
    }
    return static_cast<double>(nuc1_pairwise[1]);
  }

  double Getpw1AC(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool) {
      codondata->nuc1_pairwise(nuc1_pairwise);
      nuc1_pairwise_bool = true;
    }
    return static_cast<double>(nuc1_pairwise[0]);
  }

  double GetK80nuc1(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool) {
      codondata->nuc1_pairwise10(nuc1_pairwise);
      nuc1_pairwise_bool = true;
    }
    return nuc1_pairwise[6];
  }

  double GetK80nuc110(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool10) {
      codondata->nuc1_pairwise10(nuc1_pairwise10);
      nuc1_pairwise_bool10 = true;
    }
    return nuc1_pairwise10[6];
  }

  double GetK80nuc130(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool30) {
      codondata->nuc1_pairwise30(nuc1_pairwise30);
      nuc1_pairwise_bool30 = true;
    }
    return nuc1_pairwise30[6];
  }

  double GetK80nuc150(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool50) {
      codondata->nuc1_pairwise50(nuc1_pairwise50);
      nuc1_pairwise_bool50 = true;
    }
    return nuc1_pairwise50[6];
  }

  double GetK80nuc170(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool70) {
      codondata->nuc1_pairwise70(nuc1_pairwise70);
      nuc1_pairwise_bool70 = true;
    }
    return nuc1_pairwise70[6];
  }

  double GetK80nuc190(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool90) {
      codondata->nuc1_pairwise90(nuc1_pairwise90);
      nuc1_pairwise_bool90 = true;
    }
    return nuc1_pairwise90[6];
  }

  // 10
  double Getpw1GT10(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool10) {
      codondata->nuc1_pairwise10(nuc1_pairwise10);
      nuc1_pairwise_bool10 = true;
    }
    return nuc1_pairwise10[5];
  }

  double Getpw1CT10(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool10) {
      codondata->nuc1_pairwise10(nuc1_pairwise10);
      nuc1_pairwise_bool10 = true;
    }
    return nuc1_pairwise10[4];
  }

  double Getpw1CG10(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool10) {
      codondata->nuc1_pairwise10(nuc1_pairwise10);
      nuc1_pairwise_bool10 = true;
    }
    return nuc1_pairwise10[3];
  }

  double Getpw1AT10(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool10) {
      codondata->nuc1_pairwise10(nuc1_pairwise10);
      nuc1_pairwise_bool10 = true;
    }
    return nuc1_pairwise10[2];
  }

  double Getpw1AG10(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool10) {
      codondata->nuc1_pairwise10(nuc1_pairwise10);
      nuc1_pairwise_bool10 = true;
    }
    return nuc1_pairwise10[1];
  }

  double Getpw1AC10(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool10) {
      codondata->nuc1_pairwise10(nuc1_pairwise10);
      nuc1_pairwise_bool10 = true;
    }
    return nuc1_pairwise10[0];
  }

  // 30
  double Getpw1GT30(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool30) {
      codondata->nuc1_pairwise30(nuc1_pairwise30);
      nuc1_pairwise_bool30 = true;
    }
    return nuc1_pairwise30[5];
  }

  double Getpw1CT30(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool30) {
      codondata->nuc1_pairwise30(nuc1_pairwise30);
      nuc1_pairwise_bool30 = true;
    }
    return nuc1_pairwise30[4];
  }

  double Getpw1CG30(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool30) {
      codondata->nuc1_pairwise30(nuc1_pairwise30);
      nuc1_pairwise_bool30 = true;
    }
    return nuc1_pairwise30[3];
  }

  double Getpw1AT30(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool30) {
      codondata->nuc1_pairwise30(nuc1_pairwise30);
      nuc1_pairwise_bool30 = true;
    }
    return nuc1_pairwise30[2];
  }

  double Getpw1AG30(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool30) {
      codondata->nuc1_pairwise30(nuc1_pairwise30);
      nuc1_pairwise_bool30 = true;
    }
    return nuc1_pairwise30[1];
  }

  double Getpw1AC30(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool30) {
      codondata->nuc1_pairwise30(nuc1_pairwise30);
      nuc1_pairwise_bool30 = true;
    }
    return nuc1_pairwise30[0];
  }
  // 50

  double Getpw1GT50(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool50) {
      codondata->nuc1_pairwise50(nuc1_pairwise50);
      nuc1_pairwise_bool50 = true;
    }
    return nuc1_pairwise50[5];
  }

  double Getpw1CT50(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool50) {
      codondata->nuc1_pairwise50(nuc1_pairwise50);
      nuc1_pairwise_bool50 = true;
    }
    return nuc1_pairwise50[4];
  }

  double Getpw1CG50(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool50) {
      codondata->nuc1_pairwise50(nuc1_pairwise50);
      nuc1_pairwise_bool50 = true;
    }
    return nuc1_pairwise50[3];
  }

  double Getpw1AT50(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool50) {
      codondata->nuc1_pairwise50(nuc1_pairwise50);
      nuc1_pairwise_bool50 = true;
    }
    return nuc1_pairwise50[2];
  }

  double Getpw1AG50(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool50) {
      codondata->nuc1_pairwise50(nuc1_pairwise50);
      nuc1_pairwise_bool50 = true;
    }
    return nuc1_pairwise50[1];
  }

  double Getpw1AC50(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool50) {
      codondata->nuc1_pairwise50(nuc1_pairwise50);
      nuc1_pairwise_bool50 = true;
    }
    return nuc1_pairwise50[0];
  }

  // 70

  double Getpw1GT70(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool70) {
      codondata->nuc1_pairwise70(nuc1_pairwise70);
      nuc1_pairwise_bool70 = true;
    }
    return nuc1_pairwise70[5];
  }

  double Getpw1CT70(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool70) {
      codondata->nuc1_pairwise70(nuc1_pairwise70);
      nuc1_pairwise_bool70 = true;
    }
    return nuc1_pairwise70[4];
  }

  double Getpw1CG70(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool70) {
      codondata->nuc1_pairwise70(nuc1_pairwise70);
      nuc1_pairwise_bool70 = true;
    }
    return nuc1_pairwise70[3];
  }

  double Getpw1AT70(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool70) {
      codondata->nuc1_pairwise70(nuc1_pairwise70);
      nuc1_pairwise_bool70 = true;
    }
    return nuc1_pairwise70[2];
  }

  double Getpw1AG70(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool70) {
      codondata->nuc1_pairwise70(nuc1_pairwise70);
      nuc1_pairwise_bool70 = true;
    }
    return nuc1_pairwise70[1];
  }

  double Getpw1AC70(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool70) {
      codondata->nuc1_pairwise70(nuc1_pairwise70);
      nuc1_pairwise_bool70 = true;
    }
    return nuc1_pairwise70[0];
  }

  // 90
  double Getpw1GT90(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool90) {
      codondata->nuc1_pairwise90(nuc1_pairwise90);
      nuc1_pairwise_bool90 = true;
    }
    return nuc1_pairwise90[5];
  }

  double Getpw1CT90(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool90) {
      codondata->nuc1_pairwise90(nuc1_pairwise90);
      nuc1_pairwise_bool90 = true;
    }
    return nuc1_pairwise90[4];
  }

  double Getpw1CG90(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool90) {
      codondata->nuc1_pairwise90(nuc1_pairwise90);
      nuc1_pairwise_bool90 = true;
    }
    return nuc1_pairwise90[3];
  }

  double Getpw1AT90(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool90) {
      codondata->nuc1_pairwise90(nuc1_pairwise90);
      nuc1_pairwise_bool90 = true;
    }
    return nuc1_pairwise90[2];
  }

  double Getpw1AG90(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool90) {
      codondata->nuc1_pairwise90(nuc1_pairwise90);
      nuc1_pairwise_bool90 = true;
    }
    return nuc1_pairwise90[1];
  }

  double Getpw1AC90(CodonSequenceAlignment* codondata) {
    if (!nuc1_pairwise_bool90) {
      codondata->nuc1_pairwise90(nuc1_pairwise90);
      nuc1_pairwise_bool90 = true;
    }
    return nuc1_pairwise90[0];
  }

  /////////////////
  // nuc_pairwise
  /////////////////

  double GetpwGT(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool) {
      codondata->nuc_pairwise(nuc_pairwise);
      nuc_pairwise_bool = true;
    }
    return static_cast<double>(nuc_pairwise[5]);
  }

  double GetpwCT(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool) {
      codondata->nuc_pairwise(nuc_pairwise);
      nuc_pairwise_bool = true;
    }
    return static_cast<double>(nuc_pairwise[4]);
  }

  double GetpwCG(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool) {
      codondata->nuc_pairwise(nuc_pairwise);
      nuc_pairwise_bool = true;
    }
    return static_cast<double>(nuc_pairwise[3]);
  }

  double GetpwAT(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool) {
      codondata->nuc_pairwise(nuc_pairwise);
      nuc_pairwise_bool = true;
    }
    return static_cast<double>(nuc_pairwise[2]);
  }

  double GetpwAG(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool) {
      codondata->nuc_pairwise(nuc_pairwise);
      nuc_pairwise_bool = true;
    }
    return static_cast<double>(nuc_pairwise[1]);
  }

  double GetpwAC(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool) {
      codondata->nuc_pairwise(nuc_pairwise);
      nuc_pairwise_bool = true;
    }
    return static_cast<double>(nuc_pairwise[0]);
  }

  double GetK80nuc(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool) {
      codondata->nuc_pairwise10(nuc_pairwise);
      nuc_pairwise_bool = true;
    }
    return nuc_pairwise[6];
  }

  double GetK80nuc10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool10) {
      codondata->nuc_pairwise10(nuc_pairwise10);
      nuc_pairwise_bool10 = true;
    }
    return nuc_pairwise10[6];
  }

  double GetK80nuc30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool30) {
      codondata->nuc_pairwise30(nuc_pairwise30);
      nuc_pairwise_bool30 = true;
    }
    return nuc_pairwise30[6];
  }

  double GetK80nuc50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool50) {
      codondata->nuc_pairwise50(nuc_pairwise50);
      nuc_pairwise_bool50 = true;
    }
    return nuc_pairwise50[6];
  }

  double GetK80nuc70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool70) {
      codondata->nuc_pairwise70(nuc_pairwise70);
      nuc_pairwise_bool70 = true;
    }
    return nuc_pairwise70[6];
  }

  double GetK80nuc90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool90) {
      codondata->nuc_pairwise90(nuc_pairwise90);
      nuc_pairwise_bool90 = true;
    }
    return nuc_pairwise90[6];
  }

  // 10
  double GetpwGT10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool10) {
      codondata->nuc_pairwise10(nuc_pairwise10);
      nuc_pairwise_bool10 = true;
    }
    return nuc_pairwise10[5];
  }

  double GetpwCT10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool10) {
      codondata->nuc_pairwise10(nuc_pairwise10);
      nuc_pairwise_bool10 = true;
    }
    return nuc_pairwise10[4];
  }

  double GetpwCG10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool10) {
      codondata->nuc_pairwise10(nuc_pairwise10);
      nuc_pairwise_bool10 = true;
    }
    return nuc_pairwise10[3];
  }

  double GetpwAT10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool10) {
      codondata->nuc_pairwise10(nuc_pairwise10);
      nuc_pairwise_bool10 = true;
    }
    return nuc_pairwise10[2];
  }

  double GetpwAG10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool10) {
      codondata->nuc_pairwise10(nuc_pairwise10);
      nuc_pairwise_bool10 = true;
    }
    return nuc_pairwise10[1];
  }

  double GetpwAC10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool10) {
      codondata->nuc_pairwise10(nuc_pairwise10);
      nuc_pairwise_bool10 = true;
    }
    return nuc_pairwise10[0];
  }

  // 30
  double GetpwGT30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool30) {
      codondata->nuc_pairwise30(nuc_pairwise30);
      nuc_pairwise_bool30 = true;
    }
    return nuc_pairwise30[5];
  }

  double GetpwCT30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool30) {
      codondata->nuc_pairwise30(nuc_pairwise30);
      nuc_pairwise_bool30 = true;
    }
    return nuc_pairwise30[4];
  }

  double GetpwCG30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool30) {
      codondata->nuc_pairwise30(nuc_pairwise30);
      nuc_pairwise_bool30 = true;
    }
    return nuc_pairwise30[3];
  }

  double GetpwAT30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool30) {
      codondata->nuc_pairwise30(nuc_pairwise30);
      nuc_pairwise_bool30 = true;
    }
    return nuc_pairwise30[2];
  }

  double GetpwAG30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool30) {
      codondata->nuc_pairwise30(nuc_pairwise30);
      nuc_pairwise_bool30 = true;
    }
    return nuc_pairwise30[1];
  }

  double GetpwAC30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool30) {
      codondata->nuc_pairwise30(nuc_pairwise30);
      nuc_pairwise_bool30 = true;
    }
    return nuc_pairwise30[0];
  }
  // 50

  double GetpwGT50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool50) {
      codondata->nuc_pairwise50(nuc_pairwise50);
      nuc_pairwise_bool50 = true;
    }
    return nuc_pairwise50[5];
  }

  double GetpwCT50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool50) {
      codondata->nuc_pairwise50(nuc_pairwise50);
      nuc_pairwise_bool50 = true;
    }
    return nuc_pairwise50[4];
  }

  double GetpwCG50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool50) {
      codondata->nuc_pairwise50(nuc_pairwise50);
      nuc_pairwise_bool50 = true;
    }
    return nuc_pairwise50[3];
  }

  double GetpwAT50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool50) {
      codondata->nuc_pairwise50(nuc_pairwise50);
      nuc_pairwise_bool50 = true;
    }
    return nuc_pairwise50[2];
  }

  double GetpwAG50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool50) {
      codondata->nuc_pairwise50(nuc_pairwise50);
      nuc_pairwise_bool50 = true;
    }
    return nuc_pairwise50[1];
  }

  double GetpwAC50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool50) {
      codondata->nuc_pairwise50(nuc_pairwise50);
      nuc_pairwise_bool50 = true;
    }
    return nuc_pairwise50[0];
  }

  // 70

  double GetpwGT70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool70) {
      codondata->nuc_pairwise70(nuc_pairwise70);
      nuc_pairwise_bool70 = true;
    }
    return nuc_pairwise70[5];
  }

  double GetpwCT70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool70) {
      codondata->nuc_pairwise70(nuc_pairwise70);
      nuc_pairwise_bool70 = true;
    }
    return nuc_pairwise70[4];
  }

  double GetpwCG70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool70) {
      codondata->nuc_pairwise70(nuc_pairwise70);
      nuc_pairwise_bool70 = true;
    }
    return nuc_pairwise70[3];
  }

  double GetpwAT70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool70) {
      codondata->nuc_pairwise70(nuc_pairwise70);
      nuc_pairwise_bool70 = true;
    }
    return nuc_pairwise70[2];
  }

  double GetpwAG70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool70) {
      codondata->nuc_pairwise70(nuc_pairwise70);
      nuc_pairwise_bool70 = true;
    }
    return nuc_pairwise70[1];
  }

  double GetpwAC70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool70) {
      codondata->nuc_pairwise70(nuc_pairwise70);
      nuc_pairwise_bool70 = true;
    }
    return nuc_pairwise70[0];
  }

  // 90
  double GetpwGT90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool90) {
      codondata->nuc_pairwise90(nuc_pairwise90);
      nuc_pairwise_bool90 = true;
    }
    return nuc_pairwise90[5];
  }

  double GetpwCT90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool90) {
      codondata->nuc_pairwise90(nuc_pairwise90);
      nuc_pairwise_bool90 = true;
    }
    return nuc_pairwise90[4];
  }

  double GetpwCG90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool90) {
      codondata->nuc_pairwise90(nuc_pairwise90);
      nuc_pairwise_bool90 = true;
    }
    return nuc_pairwise90[3];
  }

  double GetpwAT90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool90) {
      codondata->nuc_pairwise90(nuc_pairwise90);
      nuc_pairwise_bool90 = true;
    }
    return nuc_pairwise90[2];
  }

  double GetpwAG90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool90) {
      codondata->nuc_pairwise90(nuc_pairwise90);
      nuc_pairwise_bool90 = true;
    }
    return nuc_pairwise90[1];
  }

  double GetpwAC90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwise_bool90) {
      codondata->nuc_pairwise90(nuc_pairwise90);
      nuc_pairwise_bool90 = true;
    }
    return nuc_pairwise90[0];
  }
  ////
  //
  ////

  double GetpwSynGT(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool) {
      codondata->nuc_pairwiseSyn(nuc_pairwiseSyn);
      nuc_pairwiseSyn_bool = true;
    }
    return static_cast<double>(nuc_pairwiseSyn[5]);
  }

  double GetpwSynCT(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool) {
      codondata->nuc_pairwiseSyn(nuc_pairwiseSyn);
      nuc_pairwiseSyn_bool = true;
    }
    return static_cast<double>(nuc_pairwiseSyn[4]);
  }

  double GetpwSynCG(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool) {
      codondata->nuc_pairwiseSyn(nuc_pairwiseSyn);
      nuc_pairwiseSyn_bool = true;
    }
    return static_cast<double>(nuc_pairwiseSyn[3]);
  }

  double GetpwSynAT(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool) {
      codondata->nuc_pairwiseSyn(nuc_pairwiseSyn);
      nuc_pairwiseSyn_bool = true;
    }
    return static_cast<double>(nuc_pairwiseSyn[2]);
  }

  double GetpwSynAG(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool) {
      codondata->nuc_pairwiseSyn(nuc_pairwiseSyn);
      nuc_pairwiseSyn_bool = true;
    }
    return static_cast<double>(nuc_pairwiseSyn[1]);
  }

  double GetpwSynAC(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool) {
      codondata->nuc_pairwiseSyn(nuc_pairwiseSyn);
      nuc_pairwiseSyn_bool = true;
    }
    return static_cast<double>(nuc_pairwiseSyn[0]);
  }

  double GetK80nucSyn(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn);
      nuc_pairwiseSyn_bool = true;
    }
    return nuc_pairwiseSyn[6];
  }

  double GetK80nucSyn10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool10) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn10);
      nuc_pairwiseSyn_bool10 = true;
    }
    return nuc_pairwiseSyn10[6];
  }

  double GetK80nucSyn30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool30) {
      codondata->nuc_pairwiseSyn30(nuc_pairwiseSyn30);
      nuc_pairwiseSyn_bool30 = true;
    }
    return nuc_pairwiseSyn30[6];
  }

  double GetK80nucSyn50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool50) {
      codondata->nuc_pairwiseSyn50(nuc_pairwiseSyn50);
      nuc_pairwiseSyn_bool50 = true;
    }
    return nuc_pairwiseSyn50[6];
  }

  double GetK80nucSyn70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool70) {
      codondata->nuc_pairwiseSyn70(nuc_pairwiseSyn70);
      nuc_pairwiseSyn_bool70 = true;
    }
    return nuc_pairwiseSyn70[6];
  }

  double GetK80nucSyn90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool90) {
      codondata->nuc_pairwiseSyn90(nuc_pairwiseSyn90);
      nuc_pairwiseSyn_bool90 = true;
    }
    return nuc_pairwiseSyn90[6];
  }

  // 10
  double GetpwSynGT10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool10) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn10);
      nuc_pairwiseSyn_bool10 = true;
    }
    return nuc_pairwiseSyn10[5];
  }

  double GetpwSynCT10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool10) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn10);
      nuc_pairwiseSyn_bool10 = true;
    }
    return nuc_pairwiseSyn10[4];
  }

  double GetpwSynCG10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool10) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn10);
      nuc_pairwiseSyn_bool10 = true;
    }
    return nuc_pairwiseSyn10[3];
  }

  double GetpwSynAT10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool10) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn10);
      nuc_pairwiseSyn_bool10 = true;
    }
    return nuc_pairwiseSyn10[2];
  }

  double GetpwSynAG10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool10) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn10);
      nuc_pairwiseSyn_bool10 = true;
    }
    return nuc_pairwiseSyn10[1];
  }

  double GetpwSynAC10(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool10) {
      codondata->nuc_pairwiseSyn10(nuc_pairwiseSyn10);
      nuc_pairwiseSyn_bool10 = true;
    }
    return nuc_pairwiseSyn10[0];
  }

  // 30
  double GetpwSynGT30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool30) {
      codondata->nuc_pairwiseSyn30(nuc_pairwiseSyn30);
      nuc_pairwiseSyn_bool30 = true;
    }
    return nuc_pairwiseSyn30[5];
  }

  double GetpwSynCT30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool30) {
      codondata->nuc_pairwiseSyn30(nuc_pairwiseSyn30);
      nuc_pairwiseSyn_bool30 = true;
    }
    return nuc_pairwiseSyn30[4];
  }

  double GetpwSynCG30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool30) {
      codondata->nuc_pairwiseSyn30(nuc_pairwiseSyn30);
      nuc_pairwiseSyn_bool30 = true;
    }
    return nuc_pairwiseSyn30[3];
  }

  double GetpwSynAT30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool30) {
      codondata->nuc_pairwiseSyn30(nuc_pairwiseSyn30);
      nuc_pairwiseSyn_bool30 = true;
    }
    return nuc_pairwiseSyn30[2];
  }

  double GetpwSynAG30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool30) {
      codondata->nuc_pairwiseSyn30(nuc_pairwiseSyn30);
      nuc_pairwiseSyn_bool30 = true;
    }
    return nuc_pairwiseSyn30[1];
  }

  double GetpwSynAC30(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool30) {
      codondata->nuc_pairwiseSyn30(nuc_pairwiseSyn30);
      nuc_pairwiseSyn_bool30 = true;
    }
    return nuc_pairwiseSyn30[0];
  }
  // 50

  double GetpwSynGT50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool50) {
      codondata->nuc_pairwiseSyn50(nuc_pairwiseSyn50);
      nuc_pairwiseSyn_bool50 = true;
    }
    return nuc_pairwiseSyn50[5];
  }

  double GetpwSynCT50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool50) {
      codondata->nuc_pairwiseSyn50(nuc_pairwiseSyn50);
      nuc_pairwiseSyn_bool50 = true;
    }
    return nuc_pairwiseSyn50[4];
  }

  double GetpwSynCG50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool50) {
      codondata->nuc_pairwiseSyn50(nuc_pairwiseSyn50);
      nuc_pairwiseSyn_bool50 = true;
    }
    return nuc_pairwiseSyn50[3];
  }

  double GetpwSynAT50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool50) {
      codondata->nuc_pairwiseSyn50(nuc_pairwiseSyn50);
      nuc_pairwiseSyn_bool50 = true;
    }
    return nuc_pairwiseSyn50[2];
  }

  double GetpwSynAG50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool50) {
      codondata->nuc_pairwiseSyn50(nuc_pairwiseSyn50);
      nuc_pairwiseSyn_bool50 = true;
    }
    return nuc_pairwiseSyn50[1];
  }

  double GetpwSynAC50(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool50) {
      codondata->nuc_pairwiseSyn50(nuc_pairwiseSyn50);
      nuc_pairwiseSyn_bool50 = true;
    }
    return nuc_pairwiseSyn50[0];
  }

  // 70

  double GetpwSynGT70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool70) {
      codondata->nuc_pairwiseSyn70(nuc_pairwiseSyn70);
      nuc_pairwiseSyn_bool70 = true;
    }
    return nuc_pairwiseSyn70[5];
  }

  double GetpwSynCT70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool70) {
      codondata->nuc_pairwiseSyn70(nuc_pairwiseSyn70);
      nuc_pairwiseSyn_bool70 = true;
    }
    return nuc_pairwiseSyn70[4];
  }

  double GetpwSynCG70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool70) {
      codondata->nuc_pairwiseSyn70(nuc_pairwiseSyn70);
      nuc_pairwiseSyn_bool70 = true;
    }
    return nuc_pairwiseSyn70[3];
  }

  double GetpwSynAT70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool70) {
      codondata->nuc_pairwiseSyn70(nuc_pairwiseSyn70);
      nuc_pairwiseSyn_bool70 = true;
    }
    return nuc_pairwiseSyn70[2];
  }

  double GetpwSynAG70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool70) {
      codondata->nuc_pairwiseSyn70(nuc_pairwiseSyn70);
      nuc_pairwiseSyn_bool70 = true;
    }
    return nuc_pairwiseSyn70[1];
  }

  double GetpwSynAC70(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool70) {
      codondata->nuc_pairwiseSyn70(nuc_pairwiseSyn70);
      nuc_pairwiseSyn_bool70 = true;
    }
    return nuc_pairwiseSyn70[0];
  }

  // 90
  double GetpwSynGT90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool90) {
      codondata->nuc_pairwiseSyn90(nuc_pairwiseSyn90);
      nuc_pairwiseSyn_bool90 = true;
    }
    return nuc_pairwiseSyn90[5];
  }

  double GetpwSynCT90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool90) {
      codondata->nuc_pairwiseSyn90(nuc_pairwiseSyn90);
      nuc_pairwiseSyn_bool90 = true;
    }
    return nuc_pairwiseSyn90[4];
  }

  double GetpwSynCG90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool90) {
      codondata->nuc_pairwiseSyn90(nuc_pairwiseSyn90);
      nuc_pairwiseSyn_bool90 = true;
    }
    return nuc_pairwiseSyn90[3];
  }

  double GetpwSynAT90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool90) {
      codondata->nuc_pairwiseSyn90(nuc_pairwiseSyn90);
      nuc_pairwiseSyn_bool90 = true;
    }
    return nuc_pairwiseSyn90[2];
  }

  double GetpwSynAG90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool90) {
      codondata->nuc_pairwiseSyn90(nuc_pairwiseSyn90);
      nuc_pairwiseSyn_bool90 = true;
    }
    return nuc_pairwiseSyn90[1];
  }

  double GetpwSynAC90(CodonSequenceAlignment* codondata) {
    if (!nuc_pairwiseSyn_bool90) {
      codondata->nuc_pairwiseSyn90(nuc_pairwiseSyn90);
      nuc_pairwiseSyn_bool90 = true;
    }
    return nuc_pairwiseSyn90[0];
  }

  //
  //
  //
  /*    double GetDIAA_AA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][0];
     }
     double GetDIAA_AC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][1];
     }
     double GetDIAA_AD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][2];
     }
     double GetDIAA_AE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][3];
     }
     double GetDIAA_AF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][4];
     }
     double GetDIAA_AG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][5];
     }
     double GetDIAA_AH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][6];
     }
     double GetDIAA_AI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][7];
     }
     double GetDIAA_AK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][8];
     }
     double GetDIAA_AL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][9];
     }
     double GetDIAA_AM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][10];
     }
     double GetDIAA_AN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][11];
     }
     double GetDIAA_AP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][12];
     }
     double GetDIAA_AQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][13];
     }
     double GetDIAA_AR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][14];
     }
     double GetDIAA_AS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][15];
     }
     double GetDIAA_AT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][16];
     }
     double GetDIAA_AV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][17];
     }
     double GetDIAA_AW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][18];
     }
     double GetDIAA_AY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[0][19];
     }
     double GetDIAA_CA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][0];
     }
     double GetDIAA_CC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][1];
     }
     double GetDIAA_CD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][2];
     }
     double GetDIAA_CE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][3];
     }
     double GetDIAA_CF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][4];
     }
     double GetDIAA_CG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][5];
     }
     double GetDIAA_CH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][6];
     }
     double GetDIAA_CI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][7];
     }
     double GetDIAA_CK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][8];
     }
     double GetDIAA_CL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][9];
     }
     double GetDIAA_CM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][10];
     }
     double GetDIAA_CN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][11];
     }
     double GetDIAA_CP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][12];
     }
     double GetDIAA_CQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][13];
     }
     double GetDIAA_CR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][14];
     }
     double GetDIAA_CS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][15];
     }
     double GetDIAA_CT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][16];
     }
     double GetDIAA_CV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][17];
     }
     double GetDIAA_CW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][18];
     }
     double GetDIAA_CY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[1][19];
     }
     double GetDIAA_DA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][0];
     }
     double GetDIAA_DC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][1];
     }
     double GetDIAA_DD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][2];
     }
     double GetDIAA_DE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][3];
     }
     double GetDIAA_DF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][4];
     }
     double GetDIAA_DG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][5];
     }
     double GetDIAA_DH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][6];
     }
     double GetDIAA_DI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][7];
     }
     double GetDIAA_DK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][8];
     }
     double GetDIAA_DL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][9];
     }
     double GetDIAA_DM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][10];
     }
     double GetDIAA_DN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][11];
     }
     double GetDIAA_DP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][12];
     }
     double GetDIAA_DQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][13];
     }
     double GetDIAA_DR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][14];
     }
     double GetDIAA_DS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][15];
     }
     double GetDIAA_DT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][16];
     }
     double GetDIAA_DV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][17];
     }
     double GetDIAA_DW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][18];
     }
     double GetDIAA_DY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[2][19];
     }
     double GetDIAA_EA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][0];
     }
     double GetDIAA_EC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][1];
     }
     double GetDIAA_ED(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][2];
     }
     double GetDIAA_EE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][3];
     }
     double GetDIAA_EF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][4];
     }
     double GetDIAA_EG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][5];
     }
     double GetDIAA_EH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][6];
     }
     double GetDIAA_EI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][7];
     }
     double GetDIAA_EK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][8];
     }
     double GetDIAA_EL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][9];
     }
     double GetDIAA_EM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][10];
     }
     double GetDIAA_EN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][11];
     }
     double GetDIAA_EP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][12];
     }
     double GetDIAA_EQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][13];
     }
     double GetDIAA_ER(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][14];
     }
     double GetDIAA_ES(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][15];
     }
     double GetDIAA_ET(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][16];
     }
     double GetDIAA_EV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][17];
     }
     double GetDIAA_EW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][18];
     }
     double GetDIAA_EY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[3][19];
     }
     double GetDIAA_FA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][0];
     }
     double GetDIAA_FC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][1];
     }
     double GetDIAA_FD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][2];
     }
     double GetDIAA_FE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][3];
     }
     double GetDIAA_FF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][4];
     }
     double GetDIAA_FG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][5];
     }
     double GetDIAA_FH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][6];
     }
     double GetDIAA_FI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][7];
     }
     double GetDIAA_FK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][8];
     }
     double GetDIAA_FL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][9];
     }
     double GetDIAA_FM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][10];
     }
     double GetDIAA_FN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][11];
     }
     double GetDIAA_FP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][12];
     }
     double GetDIAA_FQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][13];
     }
     double GetDIAA_FR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][14];
     }
     double GetDIAA_FS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][15];
     }
     double GetDIAA_FT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][16];
     }
     double GetDIAA_FV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][17];
     }
     double GetDIAA_FW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][18];
     }
     double GetDIAA_FY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[4][19];
     }
     double GetDIAA_GA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][0];
     }
     double GetDIAA_GC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][1];
     }
     double GetDIAA_GD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][2];
     }
     double GetDIAA_GE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][3];
     }
     double GetDIAA_GF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][4];
     }
     double GetDIAA_GG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][5];
     }
     double GetDIAA_GH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][6];
     }
     double GetDIAA_GI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][7];
     }
     double GetDIAA_GK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][8];
     }
     double GetDIAA_GL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][9];
     }
     double GetDIAA_GM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][10];
     }
     double GetDIAA_GN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][11];
     }
     double GetDIAA_GP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][12];
     }
     double GetDIAA_GQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][13];
     }
     double GetDIAA_GR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][14];
     }
     double GetDIAA_GS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][15];
     }
     double GetDIAA_GT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][16];
     }
     double GetDIAA_GV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][17];
     }
     double GetDIAA_GW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][18];
     }
     double GetDIAA_GY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[5][19];
     }
     double GetDIAA_HA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][0];
     }
     double GetDIAA_HC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][1];
     }
     double GetDIAA_HD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][2];
     }
     double GetDIAA_HE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][3];
     }
     double GetDIAA_HF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][4];
     }
     double GetDIAA_HG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][5];
     }
     double GetDIAA_HH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][6];
     }
     double GetDIAA_HI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][7];
     }
     double GetDIAA_HK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][8];
     }
     double GetDIAA_HL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][9];
     }
     double GetDIAA_HM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][10];
     }
     double GetDIAA_HN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][11];
     }
     double GetDIAA_HP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][12];
     }
     double GetDIAA_HQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][13];
     }
     double GetDIAA_HR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][14];
     }
     double GetDIAA_HS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][15];
     }
     double GetDIAA_HT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][16];
     }
     double GetDIAA_HV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][17];
     }
     double GetDIAA_HW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][18];
     }
     double GetDIAA_HY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[6][19];
     }
     double GetDIAA_IA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][0];
     }
     double GetDIAA_IC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][1];
     }
     double GetDIAA_ID(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][2];
     }
     double GetDIAA_IE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][3];
     }
     double GetDIAA_IF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][4];
     }
     double GetDIAA_IG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][5];
     }
     double GetDIAA_IH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][6];
     }
     double GetDIAA_II(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][7];
     }
     double GetDIAA_IK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][8];
     }
     double GetDIAA_IL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][9];
     }
     double GetDIAA_IM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][10];
     }
     double GetDIAA_IN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][11];
     }
     double GetDIAA_IP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][12];
     }
     double GetDIAA_IQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][13];
     }
     double GetDIAA_IR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][14];
     }
     double GetDIAA_IS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][15];
     }
     double GetDIAA_IT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][16];
     }
     double GetDIAA_IV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][17];
     }
     double GetDIAA_IW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][18];
     }
     double GetDIAA_IY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[7][19];
     }
     double GetDIAA_KA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][0];
     }
     double GetDIAA_KC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][1];
     }
     double GetDIAA_KD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][2];
     }
     double GetDIAA_KE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][3];
     }
     double GetDIAA_KF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][4];
     }
     double GetDIAA_KG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][5];
     }
     double GetDIAA_KH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][6];
     }
     double GetDIAA_KI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][7];
     }
     double GetDIAA_KK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][8];
     }
     double GetDIAA_KL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][9];
     }
     double GetDIAA_KM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][10];
     }
     double GetDIAA_KN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][11];
     }
     double GetDIAA_KP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][12];
     }
     double GetDIAA_KQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][13];
     }
     double GetDIAA_KR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][14];
     }
     double GetDIAA_KS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][15];
     }
     double GetDIAA_KT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][16];
     }
     double GetDIAA_KV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][17];
     }
     double GetDIAA_KW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][18];
     }
     double GetDIAA_KY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[8][19];
     }
     double GetDIAA_LA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][0];
     }
     double GetDIAA_LC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][1];
     }
     double GetDIAA_LD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][2];
     }
     double GetDIAA_LE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][3];
     }
     double GetDIAA_LF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][4];
     }
     double GetDIAA_LG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][5];
     }
     double GetDIAA_LH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][6];
     }
     double GetDIAA_LI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][7];
     }
     double GetDIAA_LK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][8];
     }
     double GetDIAA_LL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][9];
     }
     double GetDIAA_LM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][10];
     }
     double GetDIAA_LN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][11];
     }
     double GetDIAA_LP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][12];
     }
     double GetDIAA_LQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][13];
     }
     double GetDIAA_LR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][14];
     }
     double GetDIAA_LS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][15];
     }
     double GetDIAA_LT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][16];
     }
     double GetDIAA_LV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][17];
     }
     double GetDIAA_LW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][18];
     }
     double GetDIAA_LY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[9][19];
     }
     double GetDIAA_MA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][0];
     }
     double GetDIAA_MC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][1];
     }
     double GetDIAA_MD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][2];
     }
     double GetDIAA_ME(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][3];
     }
     double GetDIAA_MF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][4];
     }
     double GetDIAA_MG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][5];
     }
     double GetDIAA_MH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][6];
     }
     double GetDIAA_MI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][7];
     }
     double GetDIAA_MK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][8];
     }
     double GetDIAA_ML(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][9];
     }
     double GetDIAA_MM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][10];
     }
     double GetDIAA_MN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][11];
     }
     double GetDIAA_MP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][12];
     }
     double GetDIAA_MQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][13];
     }
     double GetDIAA_MR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][14];
     }
     double GetDIAA_MS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][15];
     }
     double GetDIAA_MT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][16];
     }
     double GetDIAA_MV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][17];
     }
     double GetDIAA_MW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][18];
     }
     double GetDIAA_MY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[10][19];
     }
     double GetDIAA_NA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][0];
     }
     double GetDIAA_NC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][1];
     }
     double GetDIAA_ND(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][2];
     }
     double GetDIAA_NE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][3];
     }
     double GetDIAA_NF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][4];
     }
     double GetDIAA_NG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][5];
     }
     double GetDIAA_NH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][6];
     }
     double GetDIAA_NI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][7];
     }
     double GetDIAA_NK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][8];
     }
     double GetDIAA_NL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][9];
     }
     double GetDIAA_NM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][10];
     }
     double GetDIAA_NN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][11];
     }
     double GetDIAA_NP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][12];
     }
     double GetDIAA_NQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][13];
     }
     double GetDIAA_NR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][14];
     }
     double GetDIAA_NS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][15];
     }
     double GetDIAA_NT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][16];
     }
     double GetDIAA_NV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][17];
     }
     double GetDIAA_NW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][18];
     }
     double GetDIAA_NY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[11][19];
     }
     double GetDIAA_PA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][0];
     }
     double GetDIAA_PC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][1];
     }
     double GetDIAA_PD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][2];
     }
     double GetDIAA_PE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][3];
     }
     double GetDIAA_PF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][4];
     }
     double GetDIAA_PG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][5];
     }
     double GetDIAA_PH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][6];
     }
     double GetDIAA_PI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][7];
     }
     double GetDIAA_PK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][8];
     }
     double GetDIAA_PL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][9];
     }
     double GetDIAA_PM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][10];
     }
     double GetDIAA_PN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][11];
     }
     double GetDIAA_PP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][12];
     }
     double GetDIAA_PQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][13];
     }
     double GetDIAA_PR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][14];
     }
     double GetDIAA_PS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][15];
     }
     double GetDIAA_PT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][16];
     }
     double GetDIAA_PV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][17];
     }
     double GetDIAA_PW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][18];
     }
     double GetDIAA_PY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[12][19];
     }
     double GetDIAA_QA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][0];
     }
     double GetDIAA_QC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][1];
     }
     double GetDIAA_QD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][2];
     }
     double GetDIAA_QE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][3];
     }
     double GetDIAA_QF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][4];
     }
     double GetDIAA_QG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][5];
     }
     double GetDIAA_QH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][6];
     }
     double GetDIAA_QI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][7];
     }
     double GetDIAA_QK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][8];
     }
     double GetDIAA_QL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][9];
     }
     double GetDIAA_QM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][10];
     }
     double GetDIAA_QN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][11];
     }
     double GetDIAA_QP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][12];
     }
     double GetDIAA_QQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][13];
     }
     double GetDIAA_QR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][14];
     }
     double GetDIAA_QS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][15];
     }
     double GetDIAA_QT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][16];
     }
     double GetDIAA_QV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][17];
     }
     double GetDIAA_QW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][18];
     }
     double GetDIAA_QY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[13][19];
     }
     double GetDIAA_RA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][0];
     }
     double GetDIAA_RC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][1];
     }
     double GetDIAA_RD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][2];
     }
     double GetDIAA_RE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][3];
     }
     double GetDIAA_RF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][4];
     }
     double GetDIAA_RG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][5];
     }
     double GetDIAA_RH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][6];
     }
     double GetDIAA_RI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][7];
     }
     double GetDIAA_RK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][8];
     }
     double GetDIAA_RL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][9];
     }
     double GetDIAA_RM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][10];
     }
     double GetDIAA_RN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][11];
     }
     double GetDIAA_RP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][12];
     }
     double GetDIAA_RQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][13];
     }
     double GetDIAA_RR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][14];
     }
     double GetDIAA_RS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][15];
     }
     double GetDIAA_RT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][16];
     }
     double GetDIAA_RV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][17];
     }
     double GetDIAA_RW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][18];
     }
     double GetDIAA_RY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[14][19];
     }
     double GetDIAA_SA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][0];
     }
     double GetDIAA_SC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][1];
     }
     double GetDIAA_SD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][2];
     }
     double GetDIAA_SE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][3];
     }
     double GetDIAA_SF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][4];
     }
     double GetDIAA_SG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][5];
     }
     double GetDIAA_SH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][6];
     }
     double GetDIAA_SI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][7];
     }
     double GetDIAA_SK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][8];
     }
     double GetDIAA_SL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][9];
     }
     double GetDIAA_SM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][10];
     }
     double GetDIAA_SN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][11];
     }
     double GetDIAA_SP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][12];
     }
     double GetDIAA_SQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][13];
     }
     double GetDIAA_SR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][14];
     }
     double GetDIAA_SS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][15];
     }
     double GetDIAA_ST(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][16];
     }
     double GetDIAA_SV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][17];
     }
     double GetDIAA_SW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][18];
     }
     double GetDIAA_SY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[15][19];
     }
     double GetDIAA_TA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][0];
     }
     double GetDIAA_TC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][1];
     }
     double GetDIAA_TD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][2];
     }
     double GetDIAA_TE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][3];
     }
     double GetDIAA_TF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][4];
     }
     double GetDIAA_TG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][5];
     }
     double GetDIAA_TH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][6];
     }
     double GetDIAA_TI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][7];
     }
     double GetDIAA_TK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][8];
     }
     double GetDIAA_TL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][9];
     }
     double GetDIAA_TM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][10];
     }
     double GetDIAA_TN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][11];
     }
     double GetDIAA_TP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][12];
     }
     double GetDIAA_TQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][13];
     }
     double GetDIAA_TR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][14];
     }
     double GetDIAA_TS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][15];
     }
     double GetDIAA_TT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][16];
     }
     double GetDIAA_TV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][17];
     }
     double GetDIAA_TW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][18];
     }
     double GetDIAA_TY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[16][19];
     }
     double GetDIAA_VA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][0];
     }
     double GetDIAA_VC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][1];
     }
     double GetDIAA_VD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][2];
     }
     double GetDIAA_VE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][3];
     }
     double GetDIAA_VF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][4];
     }
     double GetDIAA_VG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][5];
     }
     double GetDIAA_VH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][6];
     }
     double GetDIAA_VI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][7];
     }
     double GetDIAA_VK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][8];
     }
     double GetDIAA_VL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][9];
     }
     double GetDIAA_VM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][10];
     }
     double GetDIAA_VN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][11];
     }
     double GetDIAA_VP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][12];
     }
     double GetDIAA_VQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][13];
     }
     double GetDIAA_VR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][14];
     }
     double GetDIAA_VS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][15];
     }
     double GetDIAA_VT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][16];
     }
     double GetDIAA_VV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][17];
     }
     double GetDIAA_VW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][18];
     }
     double GetDIAA_VY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[17][19];
     }
     double GetDIAA_WA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][0];
     }
     double GetDIAA_WC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][1];
     }
     double GetDIAA_WD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][2];
     }
     double GetDIAA_WE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][3];
     }
     double GetDIAA_WF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][4];
     }
     double GetDIAA_WG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][5];
     }
     double GetDIAA_WH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][6];
     }
     double GetDIAA_WI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][7];
     }
     double GetDIAA_WK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][8];
     }
     double GetDIAA_WL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][9];
     }
     double GetDIAA_WM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][10];
     }
     double GetDIAA_WN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][11];
     }
     double GetDIAA_WP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][12];
     }
     double GetDIAA_WQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][13];
     }
     double GetDIAA_WR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][14];
     }
     double GetDIAA_WS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][15];
     }
     double GetDIAA_WT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][16];
     }
     double GetDIAA_WV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][17];
     }
     double GetDIAA_WW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][18];
     }
     double GetDIAA_WY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[18][19];
     }
     double GetDIAA_YA(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][0];
     }
     double GetDIAA_YC(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][1];
     }
     double GetDIAA_YD(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][2];
     }
     double GetDIAA_YE(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][3];
     }
     double GetDIAA_YF(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][4];
     }
     double GetDIAA_YG(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][5];
     }
     double GetDIAA_YH(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][6];
     }
     double GetDIAA_YI(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][7];
     }
     double GetDIAA_YK(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][8];
     }
     double GetDIAA_YL(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][9];
     }
     double GetDIAA_YM(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][10];
     }
     double GetDIAA_YN(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][11];
     }
     double GetDIAA_YP(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][12];
     }
     double GetDIAA_YQ(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][13];
     }
     double GetDIAA_YR(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][14];
     }
     double GetDIAA_YS(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][15];
     }
     double GetDIAA_YT(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][16];
     }
     double GetDIAA_YV(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][17];
     }
     double GetDIAA_YW(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][18];
     }
     double GetDIAA_YY(CodonSequenceAlignment* codondata)
     {
         if(!diaa_bool)
         {
             codondata->diaa_usage(diaa_usage);
             diaa_bool = true;
         }
         return diaa_usage[19][19];
     }


  */
};

#endif  // SOURCES_SUMMARYSTATISTICS_H_
