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

std::tuple<double, double, double> SiteInterSubMatrixCodeMLFMUTSEL::ComputeCore(
    int* nucposFrom, int* nucposTo, int codonPos, int NodeIndex, int site_nuc,
    int site_codon_i, int** CurrentNodeNucSequence) {
  double MutRate = 0.0;
  double SubRate = 0.0;
  double S = 0.0;
  int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
      nucposFrom[0], nucposFrom[1], nucposFrom[2]);
  int codonTo = lparam->codonstatespace->GetCodonFromDNA(
      nucposTo[0], nucposTo[1], nucposTo[2]);
  if (!lparam->codonstatespace->isStop(nucposTo[0], nucposTo[1], nucposTo[2])) {
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
