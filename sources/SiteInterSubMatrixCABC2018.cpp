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

#include "SiteInterSubMatrixCABC2018.h"

#include <string>

SiteInterSubMatrixCABC2018::~SiteInterSubMatrixCABC2018() {
  // dtor
}

void SiteInterSubMatrixCABC2018::init() { setSubMatrix(); }

void SiteInterSubMatrixCABC2018::initFromLeaves() { setSubMatrixFromLeaves(); }

std::tuple<double, double, double> SiteInterSubMatrixCABC2018::ComputeCore(
    int* nucposFrom, int* nucposTo, int codonPos, int NodeIndex, int site_nuc,
    int site_codon_i, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  double S = 0.0;
  int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
      nucposFrom[0], nucposFrom[1], nucposFrom[2]);
  int codonTo = lparam->codonstatespace->GetCodonFromDNA(
      nucposTo[0], nucposTo[1], nucposTo[2]);
  if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                          nucposTo[2])) {
    MutRate = lparam->gtnr[nucposFrom[codonPos]][nucposTo[codonPos]];
    int CpGcont = testCpGcontext(NodeIndex, site_nuc, nucposFrom[codonPos],
                                 nucposTo[codonPos], CurrentNodeNucSequence);

    if (CpGcont == 1 || CpGcont == 2) {
      // tsCpG
      MutRate *= lparam->lambda_CpG;
      // CpG>TpG
    }

    if (MutRate < lparam->TOOSMALL) {
      MutRate = lparam->TOOSMALL;
    }

    MutRate *= lparam->lambda_TBL;

    if (!lparam->codonstatespace->Synonymous(codonFrom, codonTo)) {
      int aaTo = lparam->codonstatespace->Translation(codonTo);
      int aaFrom = lparam->codonstatespace->Translation(codonFrom);
      S = log(
          (lparam->ssaaprofiles[lparam->alloc[site_codon_i]][aaTo] /
           lparam->ssaaprofiles[lparam->alloc[site_codon_i]][aaFrom]) *
          (lparam->codonprofile[codonTo] / lparam->codonprofile[codonFrom]));

      SubRate = MutRate * lparam->lambda_omega * lparam->omega;
    } else {
      S = log(lparam->codonprofile[codonTo] / lparam->codonprofile[codonFrom]);
      SubRate = MutRate;
    }
    SubRate = ComputeFixationFactor(S, SubRate);
  }

  return {MutRate, S, SubRate};
}

void SiteInterSubMatrixCABC2018::setSubMatrix() {
  submatrixTreeSim = new double**[lparam->refTree->GetNnode()];
  mutmatrixTreeSim = new double**[lparam->refTree->GetNnode()];
  selmatrixTreeSim = new double**[lparam->refTree->GetNnode()];
  for (int node = 0; node < lparam->refTree->GetNnode(); node++) {
    submatrixTreeSim[node] = new double*[lparam->Nsite_nuc];
    mutmatrixTreeSim[node] = new double*[lparam->Nsite_nuc];
    selmatrixTreeSim[node] = new double*[lparam->Nsite_nuc];
    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++) {
      submatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
      mutmatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
      selmatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
    }
  }

  TotalSubRate = new double[lparam->refTree->GetNnode()];
  TotalMutRate = new double[lparam->refTree->GetNnode()];
  // TotalSubRateNonSyn = new double[lparam->refTree->GetNnode()];
  // TotalMutRateNonSyn = new double[lparam->refTree->GetNnode()];
  // TotalSubRateSyn = new double[lparam->refTree->GetNnode()];
  // TotalMutRateSyn = new double[lparam->refTree->GetNnode()];

  PartialSubRate = new double[lparam->refTree->GetNnode()];
  PartialMutRate = new double[lparam->refTree->GetNnode()];
  // PartialSubRateNonSyn = new double[lparam->refTree->GetNnode()];
  // PartialMutRateNonSyn = new double[lparam->refTree->GetNnode()];
  // PartialSubRateSyn = new double[lparam->refTree->GetNnode()];
  // PartialMutRateSyn = new double[lparam->refTree->GetNnode()];
}

void SiteInterSubMatrixCABC2018::resetSubMatrix() {
  for (int node = 0; node < lparam->refTree->GetNnode(); node++) {
    TotalMutRate[node] = 0.0;
    TotalSubRate[node] = 0.0;
    // TotalMutRateNonSyn[node] = 0.0;
    // TotalSubRateNonSyn[node] = 0.0;
    // TotalMutRateSyn[node] = 0.0;
    // TotalSubRateSyn[node] = 0.0;

    PartialMutRate[node] = 0.0;
    PartialSubRate[node] = 0.0;
    // PartialMutRateNonSyn[node] = 0.0;
    // PartialSubRateNonSyn[node] = 0.0;
    // PartialMutRateSyn[node] = 0.0;
    // PartialSubRateSyn[node] = 0.0;
    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        submatrixTreeSim[node][site_nuc][nuc] = 0.0;
        mutmatrixTreeSim[node][site_nuc][nuc] = 0.0;
        selmatrixTreeSim[node][site_nuc][nuc] = 0.0;
      }
    }
  }
}

void SiteInterSubMatrixCABC2018::setSubMatrixFromLeaves() {
  submatrixTreeSim = new double**[lparam->Ntaxa];
  mutmatrixTreeSim = new double**[lparam->Ntaxa];
  selmatrixTreeSim = new double**[lparam->Ntaxa];

  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    submatrixTreeSim[taxa_i] = new double*[lparam->Nsite_nuc];
    mutmatrixTreeSim[taxa_i] = new double*[lparam->Nsite_nuc];
    selmatrixTreeSim[taxa_i] = new double*[lparam->Nsite_nuc];
    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++) {
      submatrixTreeSim[taxa_i][site_nuc] = new double[lparam->Nnucp];
      mutmatrixTreeSim[taxa_i][site_nuc] = new double[lparam->Nnucp];
      selmatrixTreeSim[taxa_i][site_nuc] = new double[lparam->Nnucp];
    }
  }
  TotalSubRate = new double[lparam->Ntaxa];
  TotalMutRate = new double[lparam->Ntaxa];
  PartialSubRate = new double[lparam->Ntaxa];
  PartialMutRate = new double[lparam->Ntaxa];
}

void SiteInterSubMatrixCABC2018::resetSubMatrixFromLeaves() {
  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    TotalMutRate[taxa_i] = 0.0;
    TotalSubRate[taxa_i] = 0.0;
    // TotalMutRateNonSyn[taxa_i] = 0.0;
    // TotalSubRateNonSyn[taxa_i] = 0.0;
    // TotalMutRateSyn[taxa_i] = 0.0;
    // TotalSubRateSyn[taxa_i] = 0.0;

    PartialMutRate[taxa_i] = 0.0;
    PartialSubRate[taxa_i] = 0.0;
    // PartialMutRateNonSyn[taxa_i] = 0.0;
    // PartialSubRateNonSyn[taxa_i] = 0.0;
    // PartialMutRateSyn[taxa_i] = 0.0;
    // PartialSubRateSyn[taxa_i] = 0.0;
    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        submatrixTreeSim[taxa_i][site_nuc][nuc] = 0.0;
        mutmatrixTreeSim[taxa_i][site_nuc][nuc] = 0.0;
        selmatrixTreeSim[taxa_i][site_nuc][nuc] = 0.0;
      }
    }
  }
}

std::tuple<double, double> SiteInterSubMatrixCABC2018::GetRatesCpG(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  std::tie(site_codon_start, site_codon_end) = getStartEndCodons(site_codon);
  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
                nucposFrom[0], nucposFrom[1], nucposFrom[2]);
            int codonTo = lparam->codonstatespace->GetCodonFromDNA(
                nucposTo[0], nucposTo[1], nucposTo[2]);
            int CpGcont =
                testCpGcontext(NodeIndex, site_nuc, nucposFrom[codonPos],
                               nucposTo[codonPos], CurrentNodeNucSequence);

            if (CpGcont == 1 || CpGcont == 2) {
              SubRate += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
              MutRate += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
  return std::make_tuple(MutRate, SubRate);
}

std::tuple<double, double> SiteInterSubMatrixCABC2018::GetRatesNonSyn(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  std::tie(site_codon_start, site_codon_end) = getStartEndCodons(site_codon);
  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
                nucposFrom[0], nucposFrom[1], nucposFrom[2]);
            int codonTo = lparam->codonstatespace->GetCodonFromDNA(
                nucposTo[0], nucposTo[1], nucposTo[2]);
            if (!lparam->codonstatespace->Synonymous(codonFrom, codonTo)) {
              SubRate += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
              MutRate += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
  return std::make_tuple(MutRate, SubRate);
}

std::tuple<double, double> SiteInterSubMatrixCABC2018::GetRatesSyn(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  std::tie(site_codon_start, site_codon_end) = getStartEndCodons(site_codon);

  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
                nucposFrom[0], nucposFrom[1], nucposFrom[2]);
            int codonTo = lparam->codonstatespace->GetCodonFromDNA(
                nucposTo[0], nucposTo[1], nucposTo[2]);
            if (lparam->codonstatespace->Synonymous(codonFrom, codonTo)) {
              SubRate += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
              MutRate += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
  return std::make_tuple(MutRate, SubRate);
}

std::tuple<double, double> SiteInterSubMatrixCABC2018::GetRatesStrongWeak(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  std::tie(site_codon_start, site_codon_end) = getStartEndCodons(site_codon);

  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
                nucposFrom[0], nucposFrom[1], nucposFrom[2]);
            int codonTo = lparam->codonstatespace->GetCodonFromDNA(
                nucposTo[0], nucposTo[1], nucposTo[2]);
            if (!isStrongWeak(nucposFrom[codonPos], nucTo)) {
              SubRate += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
              MutRate += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
  return std::make_tuple(MutRate, SubRate);
}

std::tuple<double, double> SiteInterSubMatrixCABC2018::GetRatesWeakStrong(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  std::tie(site_codon_start, site_codon_end) = getStartEndCodons(site_codon);

  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
                nucposFrom[0], nucposFrom[1], nucposFrom[2]);
            int codonTo = lparam->codonstatespace->GetCodonFromDNA(
                nucposTo[0], nucposTo[1], nucposTo[2]);
            if (isWeakStrong(nucposFrom[codonPos], nucTo)) {
              SubRate += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
              MutRate += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
  return std::make_tuple(MutRate, SubRate);
}

std::tuple<double, double> SiteInterSubMatrixCABC2018::GetRatesTransition(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  std::tie(site_codon_start, site_codon_end) = getStartEndCodons(site_codon);

  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
                nucposFrom[0], nucposFrom[1], nucposFrom[2]);
            int codonTo = lparam->codonstatespace->GetCodonFromDNA(
                nucposTo[0], nucposTo[1], nucposTo[2]);
            if (isTransition(nucposFrom[codonPos], nucTo)) {
              SubRate += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
              MutRate += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
  return std::make_tuple(MutRate, SubRate);
}

std::tuple<double, double> SiteInterSubMatrixCABC2018::GetRatesTransversion(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  std::tie(site_codon_start, site_codon_end) = getStartEndCodons(site_codon);

  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentNodeNucSequence[NodeIndex][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
                nucposFrom[0], nucposFrom[1], nucposFrom[2]);
            int codonTo = lparam->codonstatespace->GetCodonFromDNA(
                nucposTo[0], nucposTo[1], nucposTo[2]);
            if (!isTransition(nucposFrom[codonPos], nucTo)) {
              SubRate += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
              MutRate += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
  return std::make_tuple(MutRate, SubRate);
}

bool SiteInterSubMatrixCABC2018::isTransition(int nucFrom, int nucTo) {
  if ((nucFrom + nucTo) % 2 == 0) {
    return true;
  } else {
    return false;
  }
}

bool SiteInterSubMatrixCABC2018::isStrongWeak(int nucFrom, int nucTo) {
  if ((nucTo == 1 || nucTo == 2) && (nucFrom == 0 || nucFrom == 3)) {
    return true;
  } else {
    return false;
  }
}

bool SiteInterSubMatrixCABC2018::isWeakStrong(int nucFrom, int nucTo) {
  if ((nucFrom == 0 || nucFrom == 3) && (nucTo == 1 || nucTo == 2)) {
    return true;
  } else {
    return false;
  }
}

void SiteInterSubMatrixCABC2018::UpdateSubMatrixFromLeaves(
    int taxa, int** CurrentLeafNodeNucSequences) {
  double deltaTotalSubRate = 0.0;
  double deltaTotalMutRate = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  int* nucposFrom = new int[3];
  int* nucposTo = new int[3];
  for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
       site_codon_i++) {
    int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      nucposFrom[codonPos] =
          CurrentLeafNodeNucSequences[taxa][site_nuc_start + codonPos];
      nucposTo[codonPos] =
          CurrentLeafNodeNucSequences[taxa][site_nuc_start + codonPos];
    }
    for (int codonPos = 0; codonPos < 3; codonPos++) {
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        double S = 0.0;
        double MutRate = 0.0;
        double SubRate = 0.0;
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          std::tie(MutRate, S, SubRate) =
              ComputeCore(nucposFrom, nucposTo, codonPos, taxa, site_nuc,
                          site_codon_i, CurrentLeafNodeNucSequences);
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
        deltaTotalSubRate += SubRate;
        deltaTotalMutRate += MutRate;
        mutmatrixTreeSim[taxa][site_nuc][nucTo] = MutRate;
        submatrixTreeSim[taxa][site_nuc][nucTo] = SubRate;
        selmatrixTreeSim[taxa][site_nuc][nucTo] = S;
      }
    }
  }
  TotalSubRate[taxa] = deltaTotalSubRate;
  TotalMutRate[taxa] = deltaTotalMutRate;
  delete[] nucposFrom;
  delete[] nucposTo;
}