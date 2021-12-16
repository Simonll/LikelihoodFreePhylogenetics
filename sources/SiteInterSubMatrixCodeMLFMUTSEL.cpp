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

#include "SiteInterSubMatrixCodeMLFMUTSEL.h"

#include <string>

void SiteInterSubMatrixCodeMLFMUTSEL::UpdateSubMatrixTreeSim(
    int NodeIndex, int site_codon, int** CurrentNodeNucSequence) {
  // int verbose = lparam->verbose;
  bool whole = false;
  if (site_codon == -1) {
    whole = true;
  }
  double deltaTotalSubRate = 0;
  double deltaTotalMutRate = 0;

  // if -1, loop over all sites
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  // else loop over site_codon-1 to site_codon+1: 3 codons to take CpG into
  // account, which is in fact the worst case.
  if (!whole) {
    if (site_codon < lparam->Nsite_codon - 2) {
      site_codon_end = site_codon + 2;
    }
    if (site_codon > 0) {
      site_codon_start = site_codon - 1;
    }

    deltaTotalSubRate -= GetPartialSubRate(NodeIndex);
    deltaTotalMutRate -= GetPartialMutRate(NodeIndex);
  }

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
      // for each nucleotide codon postions [0,2] we will be
      // computing adjacent nucleotide.
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        double S = 0.0;
        double MutRate = 0.0;
        double SubRate = 0.0;
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          std::tie(MutRate, S, SubRate) = ComputeCore(
              MutRate, SubRate, S, nucposFrom, nucposTo, codonPos, NodeIndex,
              site_nuc, site_codon_i, CurrentNodeNucSequence);
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
        deltaTotalSubRate += SubRate;
        deltaTotalMutRate += MutRate;
        mutmatrixTreeSim[NodeIndex][site_nuc][nucTo] = MutRate;
        submatrixTreeSim[NodeIndex][site_nuc][nucTo] = SubRate;
      }
    }
  }
  if (whole) {
    TotalSubRate[NodeIndex] = deltaTotalSubRate;
    TotalMutRate[NodeIndex] = deltaTotalMutRate;
  } else {
    TotalSubRate[NodeIndex] += deltaTotalSubRate;
    TotalMutRate[NodeIndex] += deltaTotalMutRate;
  }
  delete[] nucposFrom;
  delete[] nucposTo;
}

void SiteInterSubMatrixCodeMLFMUTSEL::UpdateSubMatrixSeq(
    int taxa, int** CurrentLeafNodeNucSequences) {
  // int verbose = lparam->verbose;

  double deltaTotalSubRate = 0;
  double deltaTotalMutRate = 0;
  // if -1, loop over all sites
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  // else loop over site_codon-1 to site_codon+1: 3 codons to take CpG into
  // account, which is in fact the worst case.

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
      // for each nucleotide codon postions [0,2] we will be
      // computing adjacent nucleotide.
      int site_nuc = site_nuc_start + codonPos;
      for (int nucTo = 0; nucTo < 4; nucTo++) {
        double S = 0.0;
        double MutRate = 0.0;
        double SubRate = 0.0;
        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          std::tie(MutRate, S, SubRate) = ComputeCore(
              MutRate, SubRate, S, nucposFrom, nucposTo, codonPos, taxa,
              site_nuc, site_codon_i, CurrentLeafNodeNucSequences);
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
        ////
        // IF nucposFrom[codonPos] == nucTo set prob to 0
        ////
        deltaTotalSubRate += SubRate;
        deltaTotalMutRate += MutRate;
        mutmatrixTreeSim[taxa][site_nuc][nucTo] = MutRate;
        selmatrixTreeSim[taxa][site_nuc][nucTo] = S;
        submatrixTreeSim[taxa][site_nuc][nucTo] = SubRate;
      }
    }
  }

  TotalSubRate[taxa] = deltaTotalSubRate;
  TotalMutRate[taxa] = deltaTotalMutRate;
  delete[] nucposFrom;
  delete[] nucposTo;
}

std::tuple<double, double, double> SiteInterSubMatrixCodeMLFMUTSEL::ComputeCore(
    double MutRate, double SubRate, double S, int* nucposFrom, int* nucposTo,
    int codonPos, int NodeIndex, int site_nuc, int site_codon_i,
    int** CurrentNodeNucSequence) {
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

      SubRate =
          MutRate * lparam->lambda_omega * lparam->omega_site[site_codon_i];
    } else {
      S = log(lparam->codonprofile[codonTo] / lparam->codonprofile[codonFrom]);
      SubRate = MutRate;
    }

    SubRate = ComputeFixationFactor(S, SubRate);
  }

  return {MutRate, S, SubRate};
}
