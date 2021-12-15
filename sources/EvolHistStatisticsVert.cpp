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

#include "EvolHistStatisticsVert.h"

EvolHistStatistics::EvolHistStatistics(LocalParameters* lparam) {
  this->lparam = lparam;

  Nsub = 0.0;
  Nsynsub = 0.0;

  MutRate = new double*[2];
  SubRate = new double*[2];
  for (int i = 0; i < 2; i++) {
    MutRate[i] = new double[3];
    SubRate[i] = new double[3];
  }

  gtnr_stat = new double**[4];
  gtnrSyn_stat = new double**[4];
  gtnrNSyn_stat = new double**[4];

  for (int i = 0; i < 4; i++) {
    gtnr_stat[i] = new double*[lparam->Nnucp];
    gtnrSyn_stat[i] = new double*[lparam->Nnucp];
    gtnrNSyn_stat[i] = new double*[lparam->Nnucp];
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < lparam->Nnucp; j++) {
      gtnr_stat[i][j] = new double[lparam->Nnucp];
      gtnrSyn_stat[i][j] = new double[lparam->Nnucp];
      gtnrNSyn_stat[i][j] = new double[lparam->Nnucp];
    }
  }

  dinuc_stat = new double**[4];
  dinucSyn_stat = new double**[4];
  dinucNSyn_stat = new double**[4];
  for (int i = 0; i < 4; i++) {
    dinuc_stat[i] = new double*[lparam->Ndinuc];
    dinucSyn_stat[i] = new double*[lparam->Ndinuc];
    dinucNSyn_stat[i] = new double*[lparam->Ndinuc];
  }
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < lparam->Ndinuc; j++) {
      dinuc_stat[i][j] = new double[lparam->Ndinuc];
      dinucSyn_stat[i][j] = new double[lparam->Ndinuc];
      dinucNSyn_stat[i][j] = new double[lparam->Ndinuc];
    }
  }
  branch_stat = new double*[lparam->refTree->GetNnode()];
  for (int node = 0; node < lparam->refTree->GetNnode(); node++) {
    branch_stat[node] = new double[18];
  }

  GetEvoStatMap["gtnrAA"] = &EvolHistStatistics::GetGTNR_AA;
  GetEvoStatMap["gtnrAC"] = &EvolHistStatistics::GetGTNR_AC;
  GetEvoStatMap["gtnrAG"] = &EvolHistStatistics::GetGTNR_AG;
  GetEvoStatMap["gtnrAT"] = &EvolHistStatistics::GetGTNR_AT;
  GetEvoStatMap["gtnrCA"] = &EvolHistStatistics::GetGTNR_CA;
  GetEvoStatMap["gtnrCC"] = &EvolHistStatistics::GetGTNR_CC;
  GetEvoStatMap["gtnrCG"] = &EvolHistStatistics::GetGTNR_CG;
  GetEvoStatMap["gtnrCT"] = &EvolHistStatistics::GetGTNR_CT;
  GetEvoStatMap["gtnrGA"] = &EvolHistStatistics::GetGTNR_GA;
  GetEvoStatMap["gtnrGC"] = &EvolHistStatistics::GetGTNR_GC;
  GetEvoStatMap["gtnrGG"] = &EvolHistStatistics::GetGTNR_GG;
  GetEvoStatMap["gtnrGT"] = &EvolHistStatistics::GetGTNR_GT;
  GetEvoStatMap["gtnrTA"] = &EvolHistStatistics::GetGTNR_TA;
  GetEvoStatMap["gtnrTC"] = &EvolHistStatistics::GetGTNR_TC;
  GetEvoStatMap["gtnrTG"] = &EvolHistStatistics::GetGTNR_TG;
  GetEvoStatMap["gtnrTT"] = &EvolHistStatistics::GetGTNR_TT;

  GetEvoStatMap["gtnrSAA"] = &EvolHistStatistics::GetGTNRSyn_AA;
  GetEvoStatMap["gtnrSAC"] = &EvolHistStatistics::GetGTNRSyn_AC;
  GetEvoStatMap["gtnrSAG"] = &EvolHistStatistics::GetGTNRSyn_AG;
  GetEvoStatMap["gtnrSAT"] = &EvolHistStatistics::GetGTNRSyn_AT;
  GetEvoStatMap["gtnrSCA"] = &EvolHistStatistics::GetGTNRSyn_CA;
  GetEvoStatMap["gtnrSCC"] = &EvolHistStatistics::GetGTNRSyn_CC;
  GetEvoStatMap["gtnrSCG"] = &EvolHistStatistics::GetGTNRSyn_CG;
  GetEvoStatMap["gtnrSCT"] = &EvolHistStatistics::GetGTNRSyn_CT;
  GetEvoStatMap["gtnrSGA"] = &EvolHistStatistics::GetGTNRSyn_GA;
  GetEvoStatMap["gtnrSGC"] = &EvolHistStatistics::GetGTNRSyn_GC;
  GetEvoStatMap["gtnrSGG"] = &EvolHistStatistics::GetGTNRSyn_GG;
  GetEvoStatMap["gtnrSGT"] = &EvolHistStatistics::GetGTNRSyn_GT;
  GetEvoStatMap["gtnrSTA"] = &EvolHistStatistics::GetGTNRSyn_TA;
  GetEvoStatMap["gtnrSTC"] = &EvolHistStatistics::GetGTNRSyn_TC;
  GetEvoStatMap["gtnrSTG"] = &EvolHistStatistics::GetGTNRSyn_TG;
  GetEvoStatMap["gtnrSTT"] = &EvolHistStatistics::GetGTNRSyn_TT;

  GetEvoStatMap["gtnrNSAA"] = &EvolHistStatistics::GetGTNRNSyn_AA;
  GetEvoStatMap["gtnrNSAC"] = &EvolHistStatistics::GetGTNRNSyn_AC;
  GetEvoStatMap["gtnrNSAG"] = &EvolHistStatistics::GetGTNRNSyn_AG;
  GetEvoStatMap["gtnrNSAT"] = &EvolHistStatistics::GetGTNRNSyn_AT;
  GetEvoStatMap["gtnrNSCA"] = &EvolHistStatistics::GetGTNRNSyn_CA;
  GetEvoStatMap["gtnrNSCC"] = &EvolHistStatistics::GetGTNRNSyn_CC;
  GetEvoStatMap["gtnrNSCG"] = &EvolHistStatistics::GetGTNRNSyn_CG;
  GetEvoStatMap["gtnrNSCT"] = &EvolHistStatistics::GetGTNRNSyn_CT;
  GetEvoStatMap["gtnrNSGA"] = &EvolHistStatistics::GetGTNRNSyn_GA;
  GetEvoStatMap["gtnrNSGC"] = &EvolHistStatistics::GetGTNRNSyn_GC;
  GetEvoStatMap["gtnrNSGG"] = &EvolHistStatistics::GetGTNRNSyn_GG;
  GetEvoStatMap["gtnrNSGT"] = &EvolHistStatistics::GetGTNRNSyn_GT;
  GetEvoStatMap["gtnrNSTA"] = &EvolHistStatistics::GetGTNRNSyn_TA;
  GetEvoStatMap["gtnrNSTC"] = &EvolHistStatistics::GetGTNRNSyn_TC;
  GetEvoStatMap["gtnrNSTG"] = &EvolHistStatistics::GetGTNRNSyn_TG;
  GetEvoStatMap["gtnrNSTT"] = &EvolHistStatistics::GetGTNRNSyn_TT;

  GetEvoStatMap["dinucCGCA"] = &EvolHistStatistics::GetDinuc_CGCA;
  GetEvoStatMap["dinucCGTG"] = &EvolHistStatistics::GetDinuc_CGTG;
  GetEvoStatMap["dinuc12CGCA"] = &EvolHistStatistics::GetDinuc12_CGCA;
  GetEvoStatMap["dinuc12CGTC"] = &EvolHistStatistics::GetDinuc12_CGTC;
  GetEvoStatMap["dinuc12CGTG"] = &EvolHistStatistics::GetDinuc12_CGTG;
  GetEvoStatMap["dinuc23CGCA"] = &EvolHistStatistics::GetDinuc23_CGCA;
  GetEvoStatMap["dinuc23CGTG"] = &EvolHistStatistics::GetDinuc23_CGTG;
  GetEvoStatMap["dinuc31CGCA"] = &EvolHistStatistics::GetDinuc31_CGCA;
  GetEvoStatMap["dinuc31CGTG"] = &EvolHistStatistics::GetDinuc31_CGTG;
  GetEvoStatMap["dinucNSCGCA"] = &EvolHistStatistics::GetDinucNSyn_CGCA;
  GetEvoStatMap["dinucNSCGTG"] = &EvolHistStatistics::GetDinucNSyn_CGTG;

  GetEvoStatMap["dinuc12NSCGCA"] = &EvolHistStatistics::GetDinucNSyn12_CGCA;

  GetEvoStatMap["dinuc12NSCGTG"] = &EvolHistStatistics::GetDinucNSyn12_CGTG;

  GetEvoStatMap["dinuc23NSCGCA"] = &EvolHistStatistics::GetDinucNSyn23_CGCA;

  GetEvoStatMap["dinuc23NSCGTG"] = &EvolHistStatistics::GetDinucNSyn23_CGTG;

  GetEvoStatMap["dinuc31NSCGCA"] = &EvolHistStatistics::GetDinucNSyn31_CGCA;

  GetEvoStatMap["dinuc31NSCGTG"] = &EvolHistStatistics::GetDinucNSyn31_CGTG;

  GetEvoStatMap["dinucSCGCA"] = &EvolHistStatistics::GetDinucSyn_CGCA;

  GetEvoStatMap["dinucSCGTG"] = &EvolHistStatistics::GetDinucSyn_CGTG;

  GetEvoStatMap["dinuc12SCGCA"] = &EvolHistStatistics::GetDinucSyn12_CGCA;

  GetEvoStatMap["dinuc12SCGTG"] = &EvolHistStatistics::GetDinucSyn12_CGTG;

  GetEvoStatMap["dinuc23SCGCA"] = &EvolHistStatistics::GetDinucSyn23_CGCA;

  GetEvoStatMap["dinuc23SCGTG"] = &EvolHistStatistics::GetDinucSyn23_CGTG;

  GetEvoStatMap["dinuc31SCGCA"] = &EvolHistStatistics::GetDinucSyn31_CGCA;

  GetEvoStatMap["dinuc31SCGTG"] = &EvolHistStatistics::GetDinucSyn31_CGTG;

  GetEvoStatMap["Nsub"] = &EvolHistStatistics::GetNsub;
  GetEvoStatMap["Nsynsub"] = &EvolHistStatistics::GetNSynsub;

  GetEvoStatMap["MutRateStart"] = &EvolHistStatistics::GetMutRateStart;
  GetEvoStatMap["SubRateStart"] = &EvolHistStatistics::GetSubRateStart;
  GetEvoStatMap["MutRateNonSynStart"] =
      &EvolHistStatistics::GetMutRateNonSynStart;
  GetEvoStatMap["SubRateNonSynStart"] =
      &EvolHistStatistics::GetSubRateNonSynStart;
  GetEvoStatMap["MutRateSynStart"] = &EvolHistStatistics::GetMutRateSynStart;
  GetEvoStatMap["SubRateSynStart"] = &EvolHistStatistics::GetSubRateSynStart;
  GetEvoStatMap["MutRateEnd"] = &EvolHistStatistics::GetMutRateEnd;
  GetEvoStatMap["SubRateEnd"] = &EvolHistStatistics::GetSubRateEnd;
  GetEvoStatMap["MutRateNonSynEnd"] = &EvolHistStatistics::GetMutRateNonSynEnd;
  GetEvoStatMap["SubRateNonSynEnd"] = &EvolHistStatistics::GetSubRateNonSynEnd;
  GetEvoStatMap["MutRateSynEnd"] = &EvolHistStatistics::GetMutRateSynEnd;
  GetEvoStatMap["SubRateSynEnd"] = &EvolHistStatistics::GetSubRateSynEnd;
}

EvolHistStatistics::~EvolHistStatistics() {}

void EvolHistStatistics::GetEvoAncStats() {
  if (lparam->NusedEvoAncStats > 0) {
    std::string* arrStats = new std::string[lparam->NusedEvoAncStats];
    for (int map_i = 0; map_i < lparam->NEvoStats; map_i++) {
      auto it = lparam->mapUsedEvoAncStats.find(lparam->listEvoStats[map_i]);
      if (it != lparam->mapUsedEvoAncStats.end()) {
        if (it->second != -1) {
          arrStats[it->second] = it->first;
        }
      }
    }
    for (int map_i = 0; map_i < lparam->NusedEvoAncStats; map_i++) {
      funcpt f = GetEvoStatMap[arrStats[map_i]];
      double s = static_cast<double>((this->*f)());
      lparam->ancevostats.push_back(s);
    }
    delete[] arrStats;
  }
}

void EvolHistStatistics::GetEvoStats() {
  if (lparam->NusedEvoStats > 0) {
    string* arrStats = new string[lparam->NusedEvoStats];
    for (int map_i = 0; map_i < lparam->NEvoStats; map_i++) {
      auto it = lparam->mapUsedEvoStats.find(lparam->listEvoStats[map_i]);
      if (it != lparam->mapUsedEvoStats.end()) {
        if (it->second != -1) {
          arrStats[it->second] = it->first;
        }
      }
    }
    for (int map_i = 0; map_i < lparam->NusedEvoStats; map_i++) {
      funcpt f = GetEvoStatMap[arrStats[map_i]];
      double s = static_cast<double>((this->*f)());

      lparam->evostats.push_back(s);
    }
    delete[] arrStats;
  }
}

void EvolHistStatistics::resetEvoStats() {
  int verbose = lparam->verbose;
  for (int h = 0; h < 4; h++) {
    for (int i = 0; i < lparam->Nnucp; i++) {
      for (int j = 0; j < lparam->Nnucp; j++) {
        gtnr_stat[h][i][j] = 0;
        gtnrSyn_stat[h][i][j] = 0;
        gtnrNSyn_stat[h][i][j] = 0;
      }
    }
  }
  for (int h = 0; h < 4; h++) {
    for (int i = 0; i < lparam->Ndinuc; i++) {
      for (int j = 0; j < lparam->Ndinuc; j++) {
        dinuc_stat[h][i][j] = 0;
        dinucSyn_stat[h][i][j] = 0;
        dinucNSyn_stat[h][i][j] = 0;
      }
    }
  }
  Nsynsub = 0;
  Nsub = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      MutRate[i][j] = 0.0;
      SubRate[i][j] = 0.0;
    }
  }
}
