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
#include "SummaryStatistics.h"

SummaryStatistics::SummaryStatistics(LocalParameters* lparam) {
  // speudo count are used to ensure the possibitlity of using log
  // standardization

  this->lparam = lparam;

  // Statistics containers
  nuc_usage = new double[lparam->Nnucp];
  nuc1_usage = new double[lparam->Nnucp];
  nuc2_usage = new double[lparam->Nnucp];
  nuc3_usage = new double[lparam->Nnucp];
  relativeAAFrequency = new double[lparam->Nstate_aa];
  RSCU = new double[lparam->Nstate_codon];
  relativeCodonFrequency = new double[lparam->Nstate_codon];

  ////
  // Adding pseudo-code
  ////
  for (int nuc = 0; nuc < lparam->Nnucp; nuc++) {
    nuc_usage[nuc] = 1.0;
    nuc1_usage[nuc] = 1.0;
    nuc2_usage[nuc] = 1.0;
    nuc3_usage[nuc] = 1.0;
  }

  for (int aa = 0; aa < lparam->Nstate_aa; aa++) {
    relativeAAFrequency[aa] = 1.0;
  }

  for (int codon = 0; codon < lparam->Nstate_codon; codon++) {
    relativeCodonFrequency[codon] = 1.0;
    RSCU[codon] = 1.0;
  }

  dinuc_usage = new double*[lparam->Nnucp];
  dinuc12_usage = new double*[lparam->Nnucp];
  dinuc23_usage = new double*[lparam->Nnucp];
  dinuc31_usage = new double*[lparam->Nnucp];

  for (int nuc = 0; nuc < lparam->Nnucp; nuc++) {
    dinuc_usage[nuc] = new double[lparam->Nnucp];
    dinuc12_usage[nuc] = new double[lparam->Nnucp];
    dinuc23_usage[nuc] = new double[lparam->Nnucp];
    dinuc31_usage[nuc] = new double[lparam->Nnucp];
  }

  for (int nuc1 = 0; nuc1 < lparam->Nnucp; nuc1++) {
    for (int nuc2 = 0; nuc2 < lparam->Nnucp; nuc2++) {
      dinuc_usage[nuc1][nuc2] = 1.0;
      dinuc12_usage[nuc1][nuc2] = 1.0;
      dinuc23_usage[nuc1][nuc2] = 1.0;
      dinuc31_usage[nuc1][nuc2] = 1.0;
    }
  }

  /* diaa_usage = new double*[lparam->Nstate_aa];

  for (int aa = 0 ; aa < lparam->Nstate_aa; aa++)
  {
      diaa_usage[aa] = new double [lparam->Nstate_aa];

  }

  for (int aa1 = 0 ; aa1 < lparam->Nstate_aa; aa1++)
  {
      for (int aa2 = 0 ; aa2 < lparam->Nstate_aa; aa2++)
      {
          diaa_usage[aa1][aa2] = 1.0;

      }
  }

  dicodon_usage = new double* [lparam->Nstate_codon];

  for(int codon = 0 ; codon < lparam->Nstate_codon; codon++)
  {
      dicodon_usage[codon] = new double [lparam->Nstate_codon];

  }

  for(int codon1 = 0 ; codon1 < lparam->Nstate_codon; codon1++)
  {
      for(int codon2 = 0 ; codon2 < lparam->Nstate_codon; codon2++)
      {
          dicodon_usage[codon1][codon2] = 1.0;

      }
  } */

  nuc_meandiff = new double[2];
  nuc1_meandiff = new double[2];
  nuc2_meandiff = new double[2];
  nuc3_meandiff = new double[2];
  codon_meandiff = new double[2];
  aa_meandiff = new double[2];

  for (int i = 0; i < 2; i++) {
    nuc_meandiff[i] = 0.0;
    nuc1_meandiff[i] = 0.0;
    nuc2_meandiff[i] = 0.0;
    nuc3_meandiff[i] = 0.0;
    codon_meandiff[i] = 0.0;
    aa_meandiff[i] = 0.0;
  }

  nuc_pairwise = new double[7];
  nuc_pairwise10 = new double[7];
  nuc_pairwise30 = new double[7];
  nuc_pairwise50 = new double[7];
  nuc_pairwise70 = new double[7];
  nuc_pairwise90 = new double[7];
  nuc_pairwiseSyn = new double[7];
  nuc_pairwiseSyn10 = new double[7];
  nuc_pairwiseSyn30 = new double[7];
  nuc_pairwiseSyn50 = new double[7];
  nuc_pairwiseSyn70 = new double[7];
  nuc_pairwiseSyn90 = new double[7];
  nuc1_pairwise = new double[7];
  nuc1_pairwise10 = new double[7];
  nuc1_pairwise30 = new double[7];
  nuc1_pairwise50 = new double[7];
  nuc1_pairwise70 = new double[7];
  nuc1_pairwise90 = new double[7];
  nuc2_pairwise = new double[7];
  nuc2_pairwise10 = new double[7];
  nuc2_pairwise30 = new double[7];
  nuc2_pairwise50 = new double[7];
  nuc2_pairwise70 = new double[7];
  nuc2_pairwise90 = new double[7];
  nuc3_pairwise = new double[7];
  nuc3_pairwise10 = new double[7];
  nuc3_pairwise30 = new double[7];
  nuc3_pairwise50 = new double[7];
  nuc3_pairwise70 = new double[7];
  nuc3_pairwise90 = new double[7];

  for (int i = 0; i < 7; i++) {
    nuc_pairwise[i] = 0;
    nuc_pairwise10[i] = 0;
    nuc_pairwise30[i] = 0;
    nuc_pairwise50[i] = 0;
    nuc_pairwise70[i] = 0;
    nuc_pairwise90[i] = 0;
    nuc_pairwiseSyn[i] = 0;
    nuc_pairwiseSyn10[i] = 0;
    nuc_pairwiseSyn30[i] = 0;
    nuc_pairwiseSyn50[i] = 0;
    nuc_pairwiseSyn70[i] = 0;
    nuc_pairwiseSyn90[i] = 0;
    nuc1_pairwise[i] = 0;
    nuc1_pairwise10[i] = 0;
    nuc1_pairwise30[i] = 0;
    nuc1_pairwise50[i] = 0;
    nuc1_pairwise70[i] = 0;
    nuc1_pairwise90[i] = 0;
    nuc2_pairwise[i] = 0;
    nuc2_pairwise10[i] = 0;
    nuc2_pairwise30[i] = 0;
    nuc2_pairwise50[i] = 0;
    nuc2_pairwise70[i] = 0;
    nuc2_pairwise90[i] = 0;
    nuc3_pairwise[i] = 0;
    nuc3_pairwise10[i] = 0;
    nuc3_pairwise30[i] = 0;
    nuc3_pairwise50[i] = 0;
    nuc3_pairwise70[i] = 0;
    nuc3_pairwise90[i] = 0;
  }

  aa_pairwise = new double[2];
  aa_pairwise10 = new double[2];
  aa_pairwise30 = new double[2];
  aa_pairwise50 = new double[2];
  aa_pairwise70 = new double[2];
  aa_pairwise90 = new double[2];

  for (int i = 0; i < 2; i++) {
    aa_pairwise[i] = 0;
    aa_pairwise10[i] = 0;
    aa_pairwise30[i] = 0;
    aa_pairwise50[i] = 0;
    aa_pairwise70[i] = 0;
    aa_pairwise90[i] = 0;
  }

  dinucCpG_pairwise = new int[4];

  for (int i = 0; i < 4; i++) {
    dinucCpG_pairwise[i] = 0;
  }

  RSCU_bool = false;
  relativeCodonFrequency_bool = false;
  dinuc_bool = false;
  dinuc12_bool = false;
  dinuc23_bool = false;
  dinuc31_bool = false;
  relativeAAFrequency_bool = false;
  dicodon_bool = false;
  diaa_bool = false;
  nuc_bool = false;
  nuc1_bool = false;
  nuc2_bool = false;
  nuc3_bool = false;
  nuc_meandiff_bool = false;
  nuc1_meandiff_bool = false;
  nuc2_meandiff_bool = false;
  nuc3_meandiff_bool = false;
  codon_meandiff_bool = false;
  aa_meandiff_bool = false;
  nuc_pairwise_bool = false;
  nuc_pairwise_bool10 = false;
  nuc_pairwise_bool30 = false;
  nuc_pairwise_bool50 = false;
  nuc_pairwise_bool70 = false;
  nuc_pairwise_bool90 = false;
  nuc_pairwiseSyn_bool = false;
  nuc_pairwiseSyn_bool10 = false;
  nuc_pairwiseSyn_bool30 = false;
  nuc_pairwiseSyn_bool50 = false;
  nuc_pairwiseSyn_bool70 = false;
  nuc_pairwiseSyn_bool90 = false;
  nuc1_pairwise_bool = false;
  nuc1_pairwise_bool10 = false;
  nuc1_pairwise_bool30 = false;
  nuc1_pairwise_bool50 = false;
  nuc1_pairwise_bool70 = false;
  nuc1_pairwise_bool90 = false;
  nuc2_pairwise_bool = false;
  nuc2_pairwise_bool10 = false;
  nuc2_pairwise_bool30 = false;
  nuc2_pairwise_bool50 = false;
  nuc2_pairwise_bool70 = false;
  nuc2_pairwise_bool90 = false;
  nuc3_pairwise_bool = false;
  nuc3_pairwise_bool10 = false;
  nuc3_pairwise_bool30 = false;
  nuc3_pairwise_bool50 = false;
  nuc3_pairwise_bool70 = false;
  nuc3_pairwise_bool90 = false;
  aa_pairwise_bool = false;
  aa_pairwise_bool10 = false;
  aa_pairwise_bool30 = false;
  aa_pairwise_bool50 = false;
  aa_pairwise_bool70 = false;
  aa_pairwise_bool90 = false;
  dinucCpG_pairwise_bool = false;
  nuc_site_comphet_bool = false;
  nuc1_site_comphet_bool = false;
  nuc2_site_comphet_bool = false;
  nuc3_site_comphet_bool = false;
  nuc_taxa_comphet_bool = false;
  nuc1_taxa_comphet_bool = false;
  nuc2_taxa_comphet_bool = false;
  nuc3_taxa_comphet_bool = false;
  codon_site_comphet_bool = false;
  codon_taxa_comphet_bool = false;
  aa_site_comphet_bool = false;
  aa_taxa_comphet_bool = false;
  RSCUentropy_bool = false;
  GC_bool = false;
  GC1_bool = false;
  GC2_bool = false;
  GC3_bool = false;

  MapFunctions();
}

SummaryStatistics::SummaryStatistics(LocalData* ldata) {
  // speudo count are used to ensure the possibitlity of using log
  // standardization

  this->ldata = ldata;

  // Statistics containers
  nuc_usage = new double[ldata->Nnucp];
  nuc1_usage = new double[ldata->Nnucp];
  nuc2_usage = new double[ldata->Nnucp];
  nuc3_usage = new double[ldata->Nnucp];
  relativeAAFrequency = new double[ldata->Nstate_aa];
  RSCU = new double[ldata->Nstate_codon];
  relativeCodonFrequency = new double[ldata->Nstate_codon];

  for (int nuc = 0; nuc < ldata->Nnucp; nuc++) {
    nuc_usage[nuc] = 1.0;
    nuc1_usage[nuc] = 1.0;
    nuc2_usage[nuc] = 1.0;
    nuc3_usage[nuc] = 1.0;
  }

  for (int aa = 0; aa < ldata->Nstate_aa; aa++) {
    relativeAAFrequency[aa] = 1.0;
  }

  for (int codon = 0; codon < ldata->Nstate_codon; codon++) {
    RSCU[codon] = 1.0;
  }

  dinuc_usage = new double*[ldata->Nnucp];
  dinuc12_usage = new double*[ldata->Nnucp];
  dinuc23_usage = new double*[ldata->Nnucp];
  dinuc31_usage = new double*[ldata->Nnucp];

  for (int nuc = 0; nuc < ldata->Nnucp; nuc++) {
    dinuc_usage[nuc] = new double[ldata->Nnucp];
    dinuc12_usage[nuc] = new double[ldata->Nnucp];
    dinuc23_usage[nuc] = new double[ldata->Nnucp];
    dinuc31_usage[nuc] = new double[ldata->Nnucp];
  }

  for (int nuc1 = 0; nuc1 < ldata->Nnucp; nuc1++) {
    for (int nuc2 = 0; nuc2 < ldata->Nnucp; nuc2++) {
      dinuc_usage[nuc1][nuc2] = 1.0;
      dinuc12_usage[nuc1][nuc2] = 1.0;
      dinuc23_usage[nuc1][nuc2] = 1.0;
      dinuc31_usage[nuc1][nuc2] = 1.0;
    }
  }

  /* diaa_usage = new double*[ldata->Nstate_aa];

  for (int aa = 0 ; aa < ldata->Nstate_aa; aa++)
  {
      diaa_usage[aa] = new double [ldata->Nstate_aa];

  }

  for (int aa1 = 0 ; aa1 < ldata->Nstate_aa; aa1++)
  {
      for (int aa2 = 0 ; aa2 < ldata->Nstate_aa; aa2++)
      {
          diaa_usage[aa1][aa2] = 1.0;

      }
  }

  dicodon_usage = new double* [ldata->Nstate_codon];

  for(int codon = 0 ; codon < ldata->Nstate_codon; codon++)
  {
      dicodon_usage[codon] = new double [ldata->Nstate_codon];

  }

  for(int codon1 = 0 ; codon1 < ldata->Nstate_codon; codon1++)
  {
      for(int codon2 = 0 ; codon2 < ldata->Nstate_codon; codon2++)
      {
          dicodon_usage[codon1][codon2] = 1.0;

      }
  } */

  nuc_meandiff = new double[2];
  nuc1_meandiff = new double[2];
  nuc2_meandiff = new double[2];
  nuc3_meandiff = new double[2];
  codon_meandiff = new double[2];
  aa_meandiff = new double[2];
  for (int i = 0; i < 2; i++) {
    nuc_meandiff[i] = 0.0;
    nuc1_meandiff[i] = 0.0;
    nuc2_meandiff[i] = 0.0;
    nuc3_meandiff[i] = 0.0;
    codon_meandiff[i] = 0.0;
    aa_meandiff[i] = 0.0;
  }

  nuc_pairwise = new double[7];
  nuc_pairwise10 = new double[7];
  nuc_pairwise30 = new double[7];
  nuc_pairwise50 = new double[7];
  nuc_pairwise70 = new double[7];
  nuc_pairwise90 = new double[7];
  nuc_pairwiseSyn = new double[7];
  nuc_pairwiseSyn10 = new double[7];
  nuc_pairwiseSyn30 = new double[7];
  nuc_pairwiseSyn50 = new double[7];
  nuc_pairwiseSyn70 = new double[7];
  nuc_pairwiseSyn90 = new double[7];
  nuc1_pairwise = new double[7];
  nuc1_pairwise10 = new double[7];
  nuc1_pairwise30 = new double[7];
  nuc1_pairwise50 = new double[7];
  nuc1_pairwise70 = new double[7];
  nuc1_pairwise90 = new double[7];
  nuc2_pairwise = new double[7];
  nuc2_pairwise10 = new double[7];
  nuc2_pairwise30 = new double[7];
  nuc2_pairwise50 = new double[7];
  nuc2_pairwise70 = new double[7];
  nuc2_pairwise90 = new double[7];
  nuc3_pairwise = new double[7];
  nuc3_pairwise10 = new double[7];
  nuc3_pairwise30 = new double[7];
  nuc3_pairwise50 = new double[7];
  nuc3_pairwise70 = new double[7];
  nuc3_pairwise90 = new double[7];

  for (int i = 0; i < 7; i++) {
    nuc_pairwise[i] = 0;
    nuc_pairwise10[i] = 0;
    nuc_pairwise30[i] = 0;
    nuc_pairwise50[i] = 0;
    nuc_pairwise70[i] = 0;
    nuc_pairwise90[i] = 0;
    nuc_pairwiseSyn[i] = 0;
    nuc_pairwiseSyn10[i] = 0;
    nuc_pairwiseSyn30[i] = 0;
    nuc_pairwiseSyn50[i] = 0;
    nuc_pairwiseSyn70[i] = 0;
    nuc_pairwiseSyn90[i] = 0;
    nuc1_pairwise[i] = 0;
    nuc1_pairwise10[i] = 0;
    nuc1_pairwise30[i] = 0;
    nuc1_pairwise50[i] = 0;
    nuc1_pairwise70[i] = 0;
    nuc1_pairwise90[i] = 0;
    nuc2_pairwise[i] = 0;
    nuc2_pairwise10[i] = 0;
    nuc2_pairwise30[i] = 0;
    nuc2_pairwise50[i] = 0;
    nuc2_pairwise70[i] = 0;
    nuc2_pairwise90[i] = 0;
    nuc3_pairwise[i] = 0;
    nuc3_pairwise10[i] = 0;
    nuc3_pairwise30[i] = 0;
    nuc3_pairwise50[i] = 0;
    nuc3_pairwise70[i] = 0;
    nuc3_pairwise90[i] = 0;
  }

  aa_pairwise = new double[2];
  aa_pairwise10 = new double[2];
  aa_pairwise30 = new double[2];
  aa_pairwise50 = new double[2];
  aa_pairwise70 = new double[2];
  aa_pairwise90 = new double[2];

  for (int i = 0; i < 2; i++) {
    aa_pairwise[i] = 0;
    aa_pairwise10[i] = 0;
    aa_pairwise30[i] = 0;
    aa_pairwise50[i] = 0;
    aa_pairwise70[i] = 0;
    aa_pairwise90[i] = 0;
  }

  dinucCpG_pairwise = new int[4];

  for (int i = 0; i < 4; i++) {
    dinucCpG_pairwise[i] = 0;
  }

  // ss.computeSummaries();

  RSCU_bool = false;
  relativeCodonFrequency_bool = false;
  dinuc_bool = false;
  dinuc12_bool = false;
  dinuc23_bool = false;
  dinuc31_bool = false;
  relativeAAFrequency_bool = false;
  dicodon_bool = false;
  diaa_bool = false;
  nuc_bool = false;
  nuc1_bool = false;
  nuc2_bool = false;
  nuc3_bool = false;
  nuc_meandiff_bool = false;
  nuc1_meandiff_bool = false;
  nuc2_meandiff_bool = false;
  nuc3_meandiff_bool = false;
  codon_meandiff_bool = false;
  aa_meandiff_bool = false;
  nuc_pairwise_bool = false;
  nuc_pairwise_bool10 = false;
  nuc_pairwise_bool30 = false;
  nuc_pairwise_bool50 = false;
  nuc_pairwise_bool70 = false;
  nuc_pairwise_bool90 = false;
  nuc_pairwiseSyn_bool = false;
  nuc_pairwiseSyn_bool10 = false;
  nuc_pairwiseSyn_bool30 = false;
  nuc_pairwiseSyn_bool50 = false;
  nuc_pairwiseSyn_bool70 = false;
  nuc_pairwiseSyn_bool90 = false;
  nuc1_pairwise_bool = false;
  nuc1_pairwise_bool10 = false;
  nuc1_pairwise_bool30 = false;
  nuc1_pairwise_bool50 = false;
  nuc1_pairwise_bool70 = false;
  nuc1_pairwise_bool90 = false;
  nuc2_pairwise_bool = false;
  nuc2_pairwise_bool10 = false;
  nuc2_pairwise_bool30 = false;
  nuc2_pairwise_bool50 = false;
  nuc2_pairwise_bool70 = false;
  nuc2_pairwise_bool90 = false;
  nuc3_pairwise_bool = false;
  nuc3_pairwise_bool10 = false;
  nuc3_pairwise_bool30 = false;
  nuc3_pairwise_bool50 = false;
  nuc3_pairwise_bool70 = false;
  nuc3_pairwise_bool90 = false;
  aa_pairwise_bool = false;
  aa_pairwise_bool10 = false;
  aa_pairwise_bool30 = false;
  aa_pairwise_bool50 = false;
  aa_pairwise_bool70 = false;
  aa_pairwise_bool90 = false;
  dinucCpG_pairwise_bool = false;
  nuc_site_comphet_bool = false;
  nuc1_site_comphet_bool = false;
  nuc2_site_comphet_bool = false;
  nuc3_site_comphet_bool = false;
  nuc_taxa_comphet_bool = false;
  nuc1_taxa_comphet_bool = false;
  nuc2_taxa_comphet_bool = false;
  nuc3_taxa_comphet_bool = false;
  codon_site_comphet_bool = false;
  codon_taxa_comphet_bool = false;
  aa_site_comphet_bool = false;
  aa_taxa_comphet_bool = false;
  RSCUentropy_bool = false;
  GC_bool = false;
  GC1_bool = false;
  GC2_bool = false;
  GC3_bool = false;

  MapFunctions();
}

void SummaryStatistics::MapFunctions() {
  GetSummariesMap["K80nuc"] = &SummaryStatistics::GetK80nuc;
  GetSummariesMap["pwAC"] = &SummaryStatistics::GetpwAC;
  GetSummariesMap["pwAG"] = &SummaryStatistics::GetpwAG;
  GetSummariesMap["pwAT"] = &SummaryStatistics::GetpwAT;
  GetSummariesMap["pwCG"] = &SummaryStatistics::GetpwCG;
  GetSummariesMap["pwCT"] = &SummaryStatistics::GetpwCT;
  GetSummariesMap["pwGT"] = &SummaryStatistics::GetpwGT;

  GetSummariesMap["K80nuc10"] = &SummaryStatistics::GetK80nuc10;
  GetSummariesMap["pwAC10"] = &SummaryStatistics::GetpwAC10;
  GetSummariesMap["pwAG10"] = &SummaryStatistics::GetpwAG10;
  GetSummariesMap["pwAT10"] = &SummaryStatistics::GetpwAT10;
  GetSummariesMap["pwCG10"] = &SummaryStatistics::GetpwCG10;
  GetSummariesMap["pwCT10"] = &SummaryStatistics::GetpwCT10;
  GetSummariesMap["pwGT10"] = &SummaryStatistics::GetpwGT10;

  GetSummariesMap["K80nuc30"] = &SummaryStatistics::GetK80nuc30;
  GetSummariesMap["pwAC30"] = &SummaryStatistics::GetpwAC30;
  GetSummariesMap["pwAG30"] = &SummaryStatistics::GetpwAG30;
  GetSummariesMap["pwAT30"] = &SummaryStatistics::GetpwAT30;
  GetSummariesMap["pwCG30"] = &SummaryStatistics::GetpwCG30;
  GetSummariesMap["pwCT30"] = &SummaryStatistics::GetpwCT30;
  GetSummariesMap["pwGT30"] = &SummaryStatistics::GetpwGT30;

  GetSummariesMap["K80nuc50"] = &SummaryStatistics::GetK80nuc50;
  GetSummariesMap["pwAC50"] = &SummaryStatistics::GetpwAC50;
  GetSummariesMap["pwAG50"] = &SummaryStatistics::GetpwAG50;
  GetSummariesMap["pwAT50"] = &SummaryStatistics::GetpwAT50;
  GetSummariesMap["pwCG50"] = &SummaryStatistics::GetpwCG50;
  GetSummariesMap["pwCT50"] = &SummaryStatistics::GetpwCT50;
  GetSummariesMap["pwGT50"] = &SummaryStatistics::GetpwGT50;

  GetSummariesMap["K80nuc70"] = &SummaryStatistics::GetK80nuc70;
  GetSummariesMap["pwAC70"] = &SummaryStatistics::GetpwAC70;
  GetSummariesMap["pwAG70"] = &SummaryStatistics::GetpwAG70;
  GetSummariesMap["pwAT70"] = &SummaryStatistics::GetpwAT70;
  GetSummariesMap["pwCG70"] = &SummaryStatistics::GetpwCG70;
  GetSummariesMap["pwCT70"] = &SummaryStatistics::GetpwCT70;
  GetSummariesMap["pwGT70"] = &SummaryStatistics::GetpwGT70;

  GetSummariesMap["K80nuc90"] = &SummaryStatistics::GetK80nuc90;
  GetSummariesMap["pwAC90"] = &SummaryStatistics::GetpwAC90;
  GetSummariesMap["pwAG90"] = &SummaryStatistics::GetpwAG90;
  GetSummariesMap["pwAT90"] = &SummaryStatistics::GetpwAT90;
  GetSummariesMap["pwCG90"] = &SummaryStatistics::GetpwCG90;
  GetSummariesMap["pwCT90"] = &SummaryStatistics::GetpwCT90;
  GetSummariesMap["pwGT90"] = &SummaryStatistics::GetpwGT90;

  GetSummariesMap["K80nucSyn"] = &SummaryStatistics::GetK80nucSyn;
  GetSummariesMap["pwSynAC"] = &SummaryStatistics::GetpwSynAC;
  GetSummariesMap["pwSynAG"] = &SummaryStatistics::GetpwSynAG;
  GetSummariesMap["pwSynAT"] = &SummaryStatistics::GetpwSynAT;
  GetSummariesMap["pwSynCG"] = &SummaryStatistics::GetpwSynCG;
  GetSummariesMap["pwSynCT"] = &SummaryStatistics::GetpwSynCT;
  GetSummariesMap["pwSynGT"] = &SummaryStatistics::GetpwSynGT;

  GetSummariesMap["K80nucSyn10"] = &SummaryStatistics::GetK80nucSyn10;
  GetSummariesMap["pwSynAC10"] = &SummaryStatistics::GetpwSynAC10;
  GetSummariesMap["pwSynAG10"] = &SummaryStatistics::GetpwSynAG10;
  GetSummariesMap["pwSynAT10"] = &SummaryStatistics::GetpwSynAT10;
  GetSummariesMap["pwSynCG10"] = &SummaryStatistics::GetpwSynCG10;
  GetSummariesMap["pwSynCT10"] = &SummaryStatistics::GetpwSynCT10;
  GetSummariesMap["pwSynGT10"] = &SummaryStatistics::GetpwSynGT10;

  GetSummariesMap["K80nucSyn30"] = &SummaryStatistics::GetK80nucSyn30;
  GetSummariesMap["pwSynAC30"] = &SummaryStatistics::GetpwSynAC30;
  GetSummariesMap["pwSynAG30"] = &SummaryStatistics::GetpwSynAG30;
  GetSummariesMap["pwSynAT30"] = &SummaryStatistics::GetpwSynAT30;
  GetSummariesMap["pwSynCG30"] = &SummaryStatistics::GetpwSynCG30;
  GetSummariesMap["pwSynCT30"] = &SummaryStatistics::GetpwSynCT30;
  GetSummariesMap["pwSynGT30"] = &SummaryStatistics::GetpwSynGT30;

  GetSummariesMap["K80nucSyn50"] = &SummaryStatistics::GetK80nucSyn50;
  GetSummariesMap["pwSynAC50"] = &SummaryStatistics::GetpwSynAC50;
  GetSummariesMap["pwSynAG50"] = &SummaryStatistics::GetpwSynAG50;
  GetSummariesMap["pwSynAT50"] = &SummaryStatistics::GetpwSynAT50;
  GetSummariesMap["pwSynCG50"] = &SummaryStatistics::GetpwSynCG50;
  GetSummariesMap["pwSynCT50"] = &SummaryStatistics::GetpwSynCT50;
  GetSummariesMap["pwSynGT50"] = &SummaryStatistics::GetpwSynGT50;

  GetSummariesMap["K80nucSyn70"] = &SummaryStatistics::GetK80nucSyn70;
  GetSummariesMap["pwSynAC70"] = &SummaryStatistics::GetpwSynAC70;
  GetSummariesMap["pwSynAG70"] = &SummaryStatistics::GetpwSynAG70;
  GetSummariesMap["pwSynAT70"] = &SummaryStatistics::GetpwSynAT70;
  GetSummariesMap["pwSynCG70"] = &SummaryStatistics::GetpwSynCG70;
  GetSummariesMap["pwSynCT70"] = &SummaryStatistics::GetpwSynCT70;
  GetSummariesMap["pwSynGT70"] = &SummaryStatistics::GetpwSynGT70;

  GetSummariesMap["K80nucSyn90"] = &SummaryStatistics::GetK80nucSyn90;
  GetSummariesMap["pwSynAC90"] = &SummaryStatistics::GetpwSynAC90;
  GetSummariesMap["pwSynAG90"] = &SummaryStatistics::GetpwSynAG90;
  GetSummariesMap["pwSynAT90"] = &SummaryStatistics::GetpwSynAT90;
  GetSummariesMap["pwSynCG90"] = &SummaryStatistics::GetpwSynCG90;
  GetSummariesMap["pwSynCT90"] = &SummaryStatistics::GetpwSynCT90;
  GetSummariesMap["pwSynGT90"] = &SummaryStatistics::GetpwSynGT90;

  GetSummariesMap["K80nuc1"] = &SummaryStatistics::GetK80nuc1;
  GetSummariesMap["pw1AC"] = &SummaryStatistics::Getpw1AC;
  GetSummariesMap["pw1AG"] = &SummaryStatistics::Getpw1AG;
  GetSummariesMap["pw1AT"] = &SummaryStatistics::Getpw1AT;
  GetSummariesMap["pw1CG"] = &SummaryStatistics::Getpw1CG;
  GetSummariesMap["pw1CT"] = &SummaryStatistics::Getpw1CT;
  GetSummariesMap["pw1GT"] = &SummaryStatistics::Getpw1GT;

  GetSummariesMap["K80nuc110"] = &SummaryStatistics::GetK80nuc110;
  GetSummariesMap["pw1AC10"] = &SummaryStatistics::Getpw1AC10;
  GetSummariesMap["pw1AG10"] = &SummaryStatistics::Getpw1AG10;
  GetSummariesMap["pw1AT10"] = &SummaryStatistics::Getpw1AT10;
  GetSummariesMap["pw1CG10"] = &SummaryStatistics::Getpw1CG10;
  GetSummariesMap["pw1CT10"] = &SummaryStatistics::Getpw1CT10;
  GetSummariesMap["pw1GT10"] = &SummaryStatistics::Getpw1GT10;

  GetSummariesMap["K80nuc130"] = &SummaryStatistics::GetK80nuc130;
  GetSummariesMap["pw1AC30"] = &SummaryStatistics::Getpw1AC30;
  GetSummariesMap["pw1AG30"] = &SummaryStatistics::Getpw1AG30;
  GetSummariesMap["pw1AT30"] = &SummaryStatistics::Getpw1AT30;
  GetSummariesMap["pw1CG30"] = &SummaryStatistics::Getpw1CG30;
  GetSummariesMap["pw1CT30"] = &SummaryStatistics::Getpw1CT30;
  GetSummariesMap["pw1GT30"] = &SummaryStatistics::Getpw1GT30;

  GetSummariesMap["K80nuc150"] = &SummaryStatistics::GetK80nuc150;
  GetSummariesMap["pw1AC50"] = &SummaryStatistics::Getpw1AC50;
  GetSummariesMap["pw1AG50"] = &SummaryStatistics::Getpw1AG50;
  GetSummariesMap["pw1AT50"] = &SummaryStatistics::Getpw1AT50;
  GetSummariesMap["pw1CG50"] = &SummaryStatistics::Getpw1CG50;
  GetSummariesMap["pw1CT50"] = &SummaryStatistics::Getpw1CT50;
  GetSummariesMap["pw1GT50"] = &SummaryStatistics::Getpw1GT50;

  GetSummariesMap["K80nuc170"] = &SummaryStatistics::GetK80nuc170;
  GetSummariesMap["pw1AC70"] = &SummaryStatistics::Getpw1AC70;
  GetSummariesMap["pw1AG70"] = &SummaryStatistics::Getpw1AG70;
  GetSummariesMap["pw1AT70"] = &SummaryStatistics::Getpw1AT70;
  GetSummariesMap["pw1CG70"] = &SummaryStatistics::Getpw1CG70;
  GetSummariesMap["pw1CT70"] = &SummaryStatistics::Getpw1CT70;
  GetSummariesMap["pw1GT70"] = &SummaryStatistics::Getpw1GT70;

  GetSummariesMap["K80nuc190"] = &SummaryStatistics::GetK80nuc190;
  GetSummariesMap["pw1AC90"] = &SummaryStatistics::Getpw1AC90;
  GetSummariesMap["pw1AG90"] = &SummaryStatistics::Getpw1AG90;
  GetSummariesMap["pw1AT90"] = &SummaryStatistics::Getpw1AT90;
  GetSummariesMap["pw1CG90"] = &SummaryStatistics::Getpw1CG90;
  GetSummariesMap["pw1CT90"] = &SummaryStatistics::Getpw1CT90;
  GetSummariesMap["pw1GT90"] = &SummaryStatistics::Getpw1GT90;

  GetSummariesMap["K80nuc2"] = &SummaryStatistics::GetK80nuc2;
  GetSummariesMap["pw2AC"] = &SummaryStatistics::Getpw2AC;
  GetSummariesMap["pw2AG"] = &SummaryStatistics::Getpw2AG;
  GetSummariesMap["pw2AT"] = &SummaryStatistics::Getpw2AT;
  GetSummariesMap["pw2CG"] = &SummaryStatistics::Getpw2CG;
  GetSummariesMap["pw2CT"] = &SummaryStatistics::Getpw2CT;
  GetSummariesMap["pw2GT"] = &SummaryStatistics::Getpw2GT;

  GetSummariesMap["K80nuc210"] = &SummaryStatistics::GetK80nuc210;
  GetSummariesMap["pw2AC10"] = &SummaryStatistics::Getpw2AC10;
  GetSummariesMap["pw2AG10"] = &SummaryStatistics::Getpw2AG10;
  GetSummariesMap["pw2AT10"] = &SummaryStatistics::Getpw2AT10;
  GetSummariesMap["pw2CG10"] = &SummaryStatistics::Getpw2CG10;
  GetSummariesMap["pw2CT10"] = &SummaryStatistics::Getpw2CT10;
  GetSummariesMap["pw2GT10"] = &SummaryStatistics::Getpw2GT10;

  GetSummariesMap["K80nuc230"] = &SummaryStatistics::GetK80nuc230;
  GetSummariesMap["pw2AC30"] = &SummaryStatistics::Getpw2AC30;
  GetSummariesMap["pw2AG30"] = &SummaryStatistics::Getpw2AG30;
  GetSummariesMap["pw2AT30"] = &SummaryStatistics::Getpw2AT30;
  GetSummariesMap["pw2CG30"] = &SummaryStatistics::Getpw2CG30;
  GetSummariesMap["pw2CT30"] = &SummaryStatistics::Getpw2CT30;
  GetSummariesMap["pw2GT30"] = &SummaryStatistics::Getpw2GT30;

  GetSummariesMap["K80nuc250"] = &SummaryStatistics::GetK80nuc250;
  GetSummariesMap["pw2AC50"] = &SummaryStatistics::Getpw2AC50;
  GetSummariesMap["pw2AG50"] = &SummaryStatistics::Getpw2AG50;
  GetSummariesMap["pw2AT50"] = &SummaryStatistics::Getpw2AT50;
  GetSummariesMap["pw2CG50"] = &SummaryStatistics::Getpw2CG50;
  GetSummariesMap["pw2CT50"] = &SummaryStatistics::Getpw2CT50;
  GetSummariesMap["pw2GT50"] = &SummaryStatistics::Getpw2GT50;

  GetSummariesMap["K80nuc270"] = &SummaryStatistics::GetK80nuc270;
  GetSummariesMap["pw2AC70"] = &SummaryStatistics::Getpw2AC70;
  GetSummariesMap["pw2AG70"] = &SummaryStatistics::Getpw2AG70;
  GetSummariesMap["pw2AT70"] = &SummaryStatistics::Getpw2AT70;
  GetSummariesMap["pw2CG70"] = &SummaryStatistics::Getpw2CG70;
  GetSummariesMap["pw2CT70"] = &SummaryStatistics::Getpw2CT70;
  GetSummariesMap["pw2GT70"] = &SummaryStatistics::Getpw2GT70;

  GetSummariesMap["K80nuc290"] = &SummaryStatistics::GetK80nuc290;
  GetSummariesMap["pw2AC90"] = &SummaryStatistics::Getpw2AC90;
  GetSummariesMap["pw2AG90"] = &SummaryStatistics::Getpw2AG90;
  GetSummariesMap["pw2AT90"] = &SummaryStatistics::Getpw2AT90;
  GetSummariesMap["pw2CG90"] = &SummaryStatistics::Getpw2CG90;
  GetSummariesMap["pw2CT90"] = &SummaryStatistics::Getpw2CT90;
  GetSummariesMap["pw2GT90"] = &SummaryStatistics::Getpw2GT90;

  GetSummariesMap["K80nuc3"] = &SummaryStatistics::GetK80nuc3;
  GetSummariesMap["pw3AC"] = &SummaryStatistics::Getpw3AC;
  GetSummariesMap["pw3AG"] = &SummaryStatistics::Getpw3AG;
  GetSummariesMap["pw3AT"] = &SummaryStatistics::Getpw3AT;
  GetSummariesMap["pw3CG"] = &SummaryStatistics::Getpw3CG;
  GetSummariesMap["pw3CT"] = &SummaryStatistics::Getpw3CT;
  GetSummariesMap["pw3GT"] = &SummaryStatistics::Getpw3GT;

  GetSummariesMap["K80nuc310"] = &SummaryStatistics::GetK80nuc310;
  GetSummariesMap["pw3AC10"] = &SummaryStatistics::Getpw3AC10;
  GetSummariesMap["pw3AG10"] = &SummaryStatistics::Getpw3AG10;
  GetSummariesMap["pw3AT10"] = &SummaryStatistics::Getpw3AT10;
  GetSummariesMap["pw3CG10"] = &SummaryStatistics::Getpw3CG10;
  GetSummariesMap["pw3CT10"] = &SummaryStatistics::Getpw3CT10;
  GetSummariesMap["pw3GT10"] = &SummaryStatistics::Getpw3GT10;

  GetSummariesMap["K80nuc330"] = &SummaryStatistics::GetK80nuc330;
  GetSummariesMap["pw3AC30"] = &SummaryStatistics::Getpw3AC30;
  GetSummariesMap["pw3AG30"] = &SummaryStatistics::Getpw3AG30;
  GetSummariesMap["pw3AT30"] = &SummaryStatistics::Getpw3AT30;
  GetSummariesMap["pw3CG30"] = &SummaryStatistics::Getpw3CG30;
  GetSummariesMap["pw3CT30"] = &SummaryStatistics::Getpw3CT30;
  GetSummariesMap["pw3GT30"] = &SummaryStatistics::Getpw3GT30;

  GetSummariesMap["K80nuc350"] = &SummaryStatistics::GetK80nuc350;
  GetSummariesMap["pw3AC50"] = &SummaryStatistics::Getpw3AC50;
  GetSummariesMap["pw3AG50"] = &SummaryStatistics::Getpw3AG50;
  GetSummariesMap["pw3AT50"] = &SummaryStatistics::Getpw3AT50;
  GetSummariesMap["pw3CG50"] = &SummaryStatistics::Getpw3CG50;
  GetSummariesMap["pw3CT50"] = &SummaryStatistics::Getpw3CT50;
  GetSummariesMap["pw3GT50"] = &SummaryStatistics::Getpw3GT50;

  GetSummariesMap["K80nuc370"] = &SummaryStatistics::GetK80nuc370;
  GetSummariesMap["pw3AC70"] = &SummaryStatistics::Getpw3AC70;
  GetSummariesMap["pw3AG70"] = &SummaryStatistics::Getpw3AG70;
  GetSummariesMap["pw3AT70"] = &SummaryStatistics::Getpw3AT70;
  GetSummariesMap["pw3CG70"] = &SummaryStatistics::Getpw3CG70;
  GetSummariesMap["pw3CT70"] = &SummaryStatistics::Getpw3CT70;
  GetSummariesMap["pw3GT70"] = &SummaryStatistics::Getpw3GT70;

  GetSummariesMap["K80nuc390"] = &SummaryStatistics::GetK80nuc390;
  GetSummariesMap["pw3AC90"] = &SummaryStatistics::Getpw3AC90;
  GetSummariesMap["pw3AG90"] = &SummaryStatistics::Getpw3AG90;
  GetSummariesMap["pw3AT90"] = &SummaryStatistics::Getpw3AT90;
  GetSummariesMap["pw3CG90"] = &SummaryStatistics::Getpw3CG90;
  GetSummariesMap["pw3CT90"] = &SummaryStatistics::Getpw3CT90;
  GetSummariesMap["pw3GT90"] = &SummaryStatistics::Getpw3GT90;

  GetSummariesMap["pwAA"] = &SummaryStatistics::GetpwAA;
  GetSummariesMap["pwAA10"] = &SummaryStatistics::GetpwAA10;
  GetSummariesMap["pwAA30"] = &SummaryStatistics::GetpwAA30;
  GetSummariesMap["pwAA50"] = &SummaryStatistics::GetpwAA50;
  GetSummariesMap["pwAA70"] = &SummaryStatistics::GetpwAA70;
  GetSummariesMap["pwAA90"] = &SummaryStatistics::GetpwAA90;

  GetSummariesMap["K80aa"] = &SummaryStatistics::GetK80aa;
  GetSummariesMap["K80aa10"] = &SummaryStatistics::GetK80aa10;
  GetSummariesMap["K80aa30"] = &SummaryStatistics::GetK80aa30;
  GetSummariesMap["K80aa50"] = &SummaryStatistics::GetK80aa50;
  GetSummariesMap["K80aa70"] = &SummaryStatistics::GetK80aa70;
  GetSummariesMap["K80aa90"] = &SummaryStatistics::GetK80aa90;

  GetSummariesMap["pwCpGTpG"] = &SummaryStatistics::GetdinucCpG_TpG;
  GetSummariesMap["pwCpGCpA"] = &SummaryStatistics::GetdinucCpG_CpA;
  GetSummariesMap["pwApGTpG"] = &SummaryStatistics::GetdinucApG_TpG;

  GetSummariesMap["mAA"] = &SummaryStatistics::GetAAmean;
  GetSummariesMap["vAA"] = &SummaryStatistics::GetAAvar;

  GetSummariesMap["mCodon"] = &SummaryStatistics::GetCodonmean;
  GetSummariesMap["vCodon"] = &SummaryStatistics::GetCodonvar;

  GetSummariesMap["mNuc"] = &SummaryStatistics::GetNucmean;
  GetSummariesMap["vNuc"] = &SummaryStatistics::GetNucvar;

  GetSummariesMap["mNuc1"] = &SummaryStatistics::GetNuc1mean;
  GetSummariesMap["vNuc1"] = &SummaryStatistics::GetNuc1var;

  GetSummariesMap["mNuc2"] = &SummaryStatistics::GetNuc2mean;
  GetSummariesMap["vNuc2"] = &SummaryStatistics::GetNuc2var;

  GetSummariesMap["mNuc3"] = &SummaryStatistics::GetNuc2mean;
  GetSummariesMap["vNuc3"] = &SummaryStatistics::GetNuc2var;

  GetSummariesMap["nucA"] = &SummaryStatistics::GetNucA;
  GetSummariesMap["nucC"] = &SummaryStatistics::GetNucC;
  GetSummariesMap["nucG"] = &SummaryStatistics::GetNucG;
  GetSummariesMap["nucT"] = &SummaryStatistics::GetNucT;

  GetSummariesMap["nuc1A"] = &SummaryStatistics::GetNuc1A;
  GetSummariesMap["nuc1C"] = &SummaryStatistics::GetNuc1C;
  GetSummariesMap["nuc1G"] = &SummaryStatistics::GetNuc1G;
  GetSummariesMap["nuc1T"] = &SummaryStatistics::GetNuc1T;

  GetSummariesMap["nuc2A"] = &SummaryStatistics::GetNuc2A;
  GetSummariesMap["nuc2C"] = &SummaryStatistics::GetNuc2C;
  GetSummariesMap["nuc2G"] = &SummaryStatistics::GetNuc2G;
  GetSummariesMap["nuc2T"] = &SummaryStatistics::GetNuc2T;

  GetSummariesMap["nuc3A"] = &SummaryStatistics::GetNuc3A;
  GetSummariesMap["nuc3C"] = &SummaryStatistics::GetNuc3C;
  GetSummariesMap["nuc3G"] = &SummaryStatistics::GetNuc3G;
  GetSummariesMap["nuc3T"] = &SummaryStatistics::GetNuc3T;

  GetSummariesMap["dinucAA"] = &SummaryStatistics::GetDinucAA;
  GetSummariesMap["dinucAC"] = &SummaryStatistics::GetDinucAC;
  GetSummariesMap["dinucAG"] = &SummaryStatistics::GetDinucAG;
  GetSummariesMap["dinucAT"] = &SummaryStatistics::GetDinucAT;

  GetSummariesMap["dinucCA"] = &SummaryStatistics::GetDinucCA;
  GetSummariesMap["dinucCC"] = &SummaryStatistics::GetDinucCC;
  GetSummariesMap["dinucCG"] = &SummaryStatistics::GetDinucCG;
  GetSummariesMap["dinucCT"] = &SummaryStatistics::GetDinucCT;

  GetSummariesMap["dinucGA"] = &SummaryStatistics::GetDinucGA;
  GetSummariesMap["dinucGC"] = &SummaryStatistics::GetDinucGC;
  GetSummariesMap["dinucGG"] = &SummaryStatistics::GetDinucGG;
  GetSummariesMap["dinucGT"] = &SummaryStatistics::GetDinucGT;

  GetSummariesMap["dinucTA"] = &SummaryStatistics::GetDinucTA;
  GetSummariesMap["dinucTC"] = &SummaryStatistics::GetDinucTC;
  GetSummariesMap["dinucTG"] = &SummaryStatistics::GetDinucTG;
  GetSummariesMap["dinucTT"] = &SummaryStatistics::GetDinucTT;

  GetSummariesMap["dinuc12AA"] = &SummaryStatistics::GetDinuc12AA;
  GetSummariesMap["dinuc12AC"] = &SummaryStatistics::GetDinuc12AC;
  GetSummariesMap["dinuc12AG"] = &SummaryStatistics::GetDinuc12AG;
  GetSummariesMap["dinuc12AT"] = &SummaryStatistics::GetDinuc12AT;

  GetSummariesMap["dinuc12CA"] = &SummaryStatistics::GetDinuc12CA;
  GetSummariesMap["dinuc12CC"] = &SummaryStatistics::GetDinuc12CC;
  GetSummariesMap["dinuc12CG"] = &SummaryStatistics::GetDinuc12CG;
  GetSummariesMap["dinuc12CT"] = &SummaryStatistics::GetDinuc12CT;

  GetSummariesMap["dinuc12GA"] = &SummaryStatistics::GetDinuc12GA;
  GetSummariesMap["dinuc12GC"] = &SummaryStatistics::GetDinuc12GC;
  GetSummariesMap["dinuc12GG"] = &SummaryStatistics::GetDinuc12GG;
  GetSummariesMap["dinuc12GT"] = &SummaryStatistics::GetDinuc12GT;

  GetSummariesMap["dinuc12TA"] = &SummaryStatistics::GetDinuc12TA;
  GetSummariesMap["dinuc12TC"] = &SummaryStatistics::GetDinuc12TC;
  GetSummariesMap["dinuc12TG"] = &SummaryStatistics::GetDinuc12TG;
  GetSummariesMap["dinuc12TT"] = &SummaryStatistics::GetDinuc12TT;

  GetSummariesMap["dinuc23AA"] = &SummaryStatistics::GetDinuc23AA;
  GetSummariesMap["dinuc23AC"] = &SummaryStatistics::GetDinuc23AC;
  GetSummariesMap["dinuc23AG"] = &SummaryStatistics::GetDinuc23AG;
  GetSummariesMap["dinuc23AT"] = &SummaryStatistics::GetDinuc23AT;

  GetSummariesMap["dinuc23CA"] = &SummaryStatistics::GetDinuc23CA;
  GetSummariesMap["dinuc23CC"] = &SummaryStatistics::GetDinuc23CC;
  GetSummariesMap["dinuc23CG"] = &SummaryStatistics::GetDinuc23CG;
  GetSummariesMap["dinuc23CT"] = &SummaryStatistics::GetDinuc23CT;

  GetSummariesMap["dinuc23GA"] = &SummaryStatistics::GetDinuc23GA;
  GetSummariesMap["dinuc23GC"] = &SummaryStatistics::GetDinuc23GC;
  GetSummariesMap["dinuc23GG"] = &SummaryStatistics::GetDinuc23GG;
  GetSummariesMap["dinuc23GT"] = &SummaryStatistics::GetDinuc23GT;

  GetSummariesMap["dinuc23TA"] = &SummaryStatistics::GetDinuc23TA;
  GetSummariesMap["dinuc23TC"] = &SummaryStatistics::GetDinuc23TC;
  GetSummariesMap["dinuc23TG"] = &SummaryStatistics::GetDinuc23TG;
  GetSummariesMap["dinuc23TT"] = &SummaryStatistics::GetDinuc23TT;

  GetSummariesMap["dinuc31AA"] = &SummaryStatistics::GetDinuc31AA;
  GetSummariesMap["dinuc31AC"] = &SummaryStatistics::GetDinuc31AC;
  GetSummariesMap["dinuc31AG"] = &SummaryStatistics::GetDinuc31AG;
  GetSummariesMap["dinuc31AT"] = &SummaryStatistics::GetDinuc31AT;

  GetSummariesMap["dinuc31CA"] = &SummaryStatistics::GetDinuc31CA;
  GetSummariesMap["dinuc31CC"] = &SummaryStatistics::GetDinuc31CC;
  GetSummariesMap["dinuc31CG"] = &SummaryStatistics::GetDinuc31CG;
  GetSummariesMap["dinuc31CT"] = &SummaryStatistics::GetDinuc31CT;

  GetSummariesMap["dinuc31GA"] = &SummaryStatistics::GetDinuc31GA;
  GetSummariesMap["dinuc31GC"] = &SummaryStatistics::GetDinuc31GC;
  GetSummariesMap["dinuc31GG"] = &SummaryStatistics::GetDinuc31GG;
  GetSummariesMap["dinuc31GT"] = &SummaryStatistics::GetDinuc31GT;

  GetSummariesMap["dinuc31TA"] = &SummaryStatistics::GetDinuc31TA;
  GetSummariesMap["dinuc31TC"] = &SummaryStatistics::GetDinuc31TC;
  GetSummariesMap["dinuc31TG"] = &SummaryStatistics::GetDinuc31TG;
  GetSummariesMap["dinuc31TT"] = &SummaryStatistics::GetDinuc31TT;

  GetSummariesMap["A"] = &SummaryStatistics::GetA;
  GetSummariesMap["C"] = &SummaryStatistics::GetC;
  GetSummariesMap["D"] = &SummaryStatistics::GetD;
  GetSummariesMap["E"] = &SummaryStatistics::GetE;
  GetSummariesMap["F"] = &SummaryStatistics::GetF;
  GetSummariesMap["G"] = &SummaryStatistics::GetG;
  GetSummariesMap["H"] = &SummaryStatistics::GetH;
  GetSummariesMap["I"] = &SummaryStatistics::GetI;
  GetSummariesMap["K"] = &SummaryStatistics::GetK;
  GetSummariesMap["L"] = &SummaryStatistics::GetL;
  GetSummariesMap["M"] = &SummaryStatistics::GetM;
  GetSummariesMap["N"] = &SummaryStatistics::GetN;
  GetSummariesMap["P"] = &SummaryStatistics::GetP;
  GetSummariesMap["Q"] = &SummaryStatistics::GetQ;
  GetSummariesMap["R"] = &SummaryStatistics::GetR;
  GetSummariesMap["S"] = &SummaryStatistics::GetS;
  GetSummariesMap["T"] = &SummaryStatistics::GetT;
  GetSummariesMap["V"] = &SummaryStatistics::GetV;
  GetSummariesMap["W"] = &SummaryStatistics::GetW;
  GetSummariesMap["Y"] = &SummaryStatistics::GetY;

  GetSummariesMap["TTT"] = &SummaryStatistics::GetTTT;
  GetSummariesMap["TTC"] = &SummaryStatistics::GetTTC;
  GetSummariesMap["TTA"] = &SummaryStatistics::GetTTA;
  GetSummariesMap["TTG"] = &SummaryStatistics::GetTTG;
  GetSummariesMap["TCT"] = &SummaryStatistics::GetTCT;
  GetSummariesMap["TCC"] = &SummaryStatistics::GetTCC;
  GetSummariesMap["TCA"] = &SummaryStatistics::GetTCA;
  GetSummariesMap["TCG"] = &SummaryStatistics::GetTCG;
  GetSummariesMap["TAT"] = &SummaryStatistics::GetTAT;
  GetSummariesMap["TAC"] = &SummaryStatistics::GetTAC;
  GetSummariesMap["TAA"] = &SummaryStatistics::GetTAA;
  GetSummariesMap["TAG"] = &SummaryStatistics::GetTAG;
  GetSummariesMap["TGT"] = &SummaryStatistics::GetTGT;
  GetSummariesMap["TGC"] = &SummaryStatistics::GetTGC;
  GetSummariesMap["TGA"] = &SummaryStatistics::GetTGA;
  GetSummariesMap["TGG"] = &SummaryStatistics::GetTGG;
  GetSummariesMap["CTT"] = &SummaryStatistics::GetCTT;
  GetSummariesMap["CTC"] = &SummaryStatistics::GetCTC;
  GetSummariesMap["CTA"] = &SummaryStatistics::GetCTA;
  GetSummariesMap["CTG"] = &SummaryStatistics::GetCTG;
  GetSummariesMap["CCT"] = &SummaryStatistics::GetCCT;
  GetSummariesMap["CCC"] = &SummaryStatistics::GetCCC;
  GetSummariesMap["CCA"] = &SummaryStatistics::GetCCA;
  GetSummariesMap["CCG"] = &SummaryStatistics::GetCCG;
  GetSummariesMap["CAT"] = &SummaryStatistics::GetCAT;
  GetSummariesMap["CAC"] = &SummaryStatistics::GetCAC;
  GetSummariesMap["CAA"] = &SummaryStatistics::GetCAA;
  GetSummariesMap["CAG"] = &SummaryStatistics::GetCAG;
  GetSummariesMap["CGT"] = &SummaryStatistics::GetCGT;
  GetSummariesMap["CGC"] = &SummaryStatistics::GetCGC;
  GetSummariesMap["CGA"] = &SummaryStatistics::GetCGA;
  GetSummariesMap["CGG"] = &SummaryStatistics::GetCGG;
  GetSummariesMap["ATT"] = &SummaryStatistics::GetATT;
  GetSummariesMap["ATC"] = &SummaryStatistics::GetATC;
  GetSummariesMap["ATA"] = &SummaryStatistics::GetATA;
  GetSummariesMap["ATG"] = &SummaryStatistics::GetATG;
  GetSummariesMap["ACT"] = &SummaryStatistics::GetACT;
  GetSummariesMap["ACC"] = &SummaryStatistics::GetACC;
  GetSummariesMap["ACA"] = &SummaryStatistics::GetACA;
  GetSummariesMap["ACG"] = &SummaryStatistics::GetACG;
  GetSummariesMap["AAT"] = &SummaryStatistics::GetAAT;
  GetSummariesMap["AAC"] = &SummaryStatistics::GetAAC;
  GetSummariesMap["AAA"] = &SummaryStatistics::GetAAA;
  GetSummariesMap["AAG"] = &SummaryStatistics::GetAAG;
  GetSummariesMap["AGT"] = &SummaryStatistics::GetAGT;
  GetSummariesMap["AGC"] = &SummaryStatistics::GetAGC;
  GetSummariesMap["AGA"] = &SummaryStatistics::GetAGA;
  GetSummariesMap["AGG"] = &SummaryStatistics::GetAGG;
  GetSummariesMap["GTT"] = &SummaryStatistics::GetGTT;
  GetSummariesMap["GTC"] = &SummaryStatistics::GetGTC;
  GetSummariesMap["GTA"] = &SummaryStatistics::GetGTA;
  GetSummariesMap["GTG"] = &SummaryStatistics::GetGTG;
  GetSummariesMap["GCT"] = &SummaryStatistics::GetGCT;
  GetSummariesMap["GCC"] = &SummaryStatistics::GetGCC;
  GetSummariesMap["GCA"] = &SummaryStatistics::GetGCA;
  GetSummariesMap["GCG"] = &SummaryStatistics::GetGCG;
  GetSummariesMap["GAT"] = &SummaryStatistics::GetGAT;
  GetSummariesMap["GAC"] = &SummaryStatistics::GetGAC;
  GetSummariesMap["GAA"] = &SummaryStatistics::GetGAA;
  GetSummariesMap["GAG"] = &SummaryStatistics::GetGAG;
  GetSummariesMap["GGT"] = &SummaryStatistics::GetGGT;
  GetSummariesMap["GGC"] = &SummaryStatistics::GetGGC;
  GetSummariesMap["GGA"] = &SummaryStatistics::GetGGA;
  GetSummariesMap["GGG"] = &SummaryStatistics::GetGGG;

  GetSummariesMap["fTTT"] = &SummaryStatistics::GetfTTT;
  GetSummariesMap["fTTC"] = &SummaryStatistics::GetfTTC;
  GetSummariesMap["fTTA"] = &SummaryStatistics::GetfTTA;
  GetSummariesMap["fTTG"] = &SummaryStatistics::GetfTTG;
  GetSummariesMap["fTCT"] = &SummaryStatistics::GetfTCT;
  GetSummariesMap["fTCC"] = &SummaryStatistics::GetfTCC;
  GetSummariesMap["fTCA"] = &SummaryStatistics::GetfTCA;
  GetSummariesMap["fTCG"] = &SummaryStatistics::GetfTCG;
  GetSummariesMap["fTAT"] = &SummaryStatistics::GetfTAT;
  GetSummariesMap["fTAC"] = &SummaryStatistics::GetfTAC;
  GetSummariesMap["fTAA"] = &SummaryStatistics::GetfTAA;
  GetSummariesMap["fTAG"] = &SummaryStatistics::GetfTAG;
  GetSummariesMap["fTGT"] = &SummaryStatistics::GetfTGT;
  GetSummariesMap["fTGC"] = &SummaryStatistics::GetfTGC;
  GetSummariesMap["fTGA"] = &SummaryStatistics::GetfTGA;
  GetSummariesMap["fTGG"] = &SummaryStatistics::GetfTGG;
  GetSummariesMap["fCTT"] = &SummaryStatistics::GetfCTT;
  GetSummariesMap["fCTC"] = &SummaryStatistics::GetfCTC;
  GetSummariesMap["fCTA"] = &SummaryStatistics::GetfCTA;
  GetSummariesMap["fCTG"] = &SummaryStatistics::GetfCTG;
  GetSummariesMap["fCCT"] = &SummaryStatistics::GetfCCT;
  GetSummariesMap["fCCC"] = &SummaryStatistics::GetfCCC;
  GetSummariesMap["fCCA"] = &SummaryStatistics::GetfCCA;
  GetSummariesMap["fCCG"] = &SummaryStatistics::GetfCCG;
  GetSummariesMap["fCAT"] = &SummaryStatistics::GetfCAT;
  GetSummariesMap["fCAC"] = &SummaryStatistics::GetfCAC;
  GetSummariesMap["fCAA"] = &SummaryStatistics::GetfCAA;
  GetSummariesMap["fCAG"] = &SummaryStatistics::GetfCAG;
  GetSummariesMap["fCGT"] = &SummaryStatistics::GetfCGT;
  GetSummariesMap["fCGC"] = &SummaryStatistics::GetfCGC;
  GetSummariesMap["fCGA"] = &SummaryStatistics::GetfCGA;
  GetSummariesMap["fCGG"] = &SummaryStatistics::GetfCGG;
  GetSummariesMap["fATT"] = &SummaryStatistics::GetfATT;
  GetSummariesMap["fATC"] = &SummaryStatistics::GetfATC;
  GetSummariesMap["fATA"] = &SummaryStatistics::GetfATA;
  GetSummariesMap["fATG"] = &SummaryStatistics::GetfATG;
  GetSummariesMap["fACT"] = &SummaryStatistics::GetfACT;
  GetSummariesMap["fACC"] = &SummaryStatistics::GetfACC;
  GetSummariesMap["fACA"] = &SummaryStatistics::GetfACA;
  GetSummariesMap["fACG"] = &SummaryStatistics::GetfACG;
  GetSummariesMap["fAAT"] = &SummaryStatistics::GetfAAT;
  GetSummariesMap["fAAC"] = &SummaryStatistics::GetfAAC;
  GetSummariesMap["fAAA"] = &SummaryStatistics::GetfAAA;
  GetSummariesMap["fAAG"] = &SummaryStatistics::GetfAAG;
  GetSummariesMap["fAGT"] = &SummaryStatistics::GetfAGT;
  GetSummariesMap["fAGC"] = &SummaryStatistics::GetfAGC;
  GetSummariesMap["fAGA"] = &SummaryStatistics::GetfAGA;
  GetSummariesMap["fAGG"] = &SummaryStatistics::GetfAGG;
  GetSummariesMap["fGTT"] = &SummaryStatistics::GetfGTT;
  GetSummariesMap["fGTC"] = &SummaryStatistics::GetfGTC;
  GetSummariesMap["fGTA"] = &SummaryStatistics::GetfGTA;
  GetSummariesMap["fGTG"] = &SummaryStatistics::GetfGTG;
  GetSummariesMap["fGCT"] = &SummaryStatistics::GetfGCT;
  GetSummariesMap["fGCC"] = &SummaryStatistics::GetfGCC;
  GetSummariesMap["fGCA"] = &SummaryStatistics::GetfGCA;
  GetSummariesMap["fGCG"] = &SummaryStatistics::GetfGCG;
  GetSummariesMap["fGAT"] = &SummaryStatistics::GetfGAT;
  GetSummariesMap["fGAC"] = &SummaryStatistics::GetfGAC;
  GetSummariesMap["fGAA"] = &SummaryStatistics::GetfGAA;
  GetSummariesMap["fGAG"] = &SummaryStatistics::GetfGAG;
  GetSummariesMap["fGGT"] = &SummaryStatistics::GetfGGT;
  GetSummariesMap["fGGC"] = &SummaryStatistics::GetfGGC;
  GetSummariesMap["fGGA"] = &SummaryStatistics::GetfGGA;
  GetSummariesMap["fGGG"] = &SummaryStatistics::GetfGGG;

  GetSummariesMap["nucsitecomphet"] = &SummaryStatistics::Getnuc_site_comphet;
  GetSummariesMap["nuc1sitecomphet"] = &SummaryStatistics::Getnuc1_site_comphet;
  GetSummariesMap["nuc2sitecomphet"] = &SummaryStatistics::Getnuc2_site_comphet;
  GetSummariesMap["nuc3sitecomphet"] = &SummaryStatistics::Getnuc3_site_comphet;
  GetSummariesMap["nuc_taxacomphet"] = &SummaryStatistics::Getnuc_taxa_comphet;
  GetSummariesMap["nuc1taxacomphet"] = &SummaryStatistics::Getnuc1_taxa_comphet;
  GetSummariesMap["nuc2taxacomphet"] = &SummaryStatistics::Getnuc2_taxa_comphet;
  GetSummariesMap["nuc3taxacomphet"] = &SummaryStatistics::Getnuc3_taxa_comphet;
  GetSummariesMap["codonsitecomphet"] =
      &SummaryStatistics::Getcodon_site_comphet;
  GetSummariesMap["codontaxacomphet"] =
      &SummaryStatistics::Getcodon_taxa_comphet;
  GetSummariesMap["aasitecomphet"] = &SummaryStatistics::Getaa_site_comphet;
  GetSummariesMap["aataxacomphet"] = &SummaryStatistics::Getaa_taxa_comphet;

  /* GetSummariesMap["diaaAA"]=&SummaryStatistics::GetDIAA_AA;
  GetSummariesMap["diaaAC"]=&SummaryStatistics::GetDIAA_AC;
  GetSummariesMap["diaaAD"]=&SummaryStatistics::GetDIAA_AD;
  GetSummariesMap["diaaAE"]=&SummaryStatistics::GetDIAA_AE;
  GetSummariesMap["diaaAF"]=&SummaryStatistics::GetDIAA_AF;
  GetSummariesMap["diaaAG"]=&SummaryStatistics::GetDIAA_AG;
  GetSummariesMap["diaaAH"]=&SummaryStatistics::GetDIAA_AH;
  GetSummariesMap["diaaAI"]=&SummaryStatistics::GetDIAA_AI;
  GetSummariesMap["diaaAK"]=&SummaryStatistics::GetDIAA_AK;
  GetSummariesMap["diaaAL"]=&SummaryStatistics::GetDIAA_AL;
  GetSummariesMap["diaaAM"]=&SummaryStatistics::GetDIAA_AM;
  GetSummariesMap["diaaAN"]=&SummaryStatistics::GetDIAA_AN;
  GetSummariesMap["diaaAP"]=&SummaryStatistics::GetDIAA_AP;
  GetSummariesMap["diaaAQ"]=&SummaryStatistics::GetDIAA_AQ;
  GetSummariesMap["diaaAR"]=&SummaryStatistics::GetDIAA_AR;
  GetSummariesMap["diaaAS"]=&SummaryStatistics::GetDIAA_AS;
  GetSummariesMap["diaaAT"]=&SummaryStatistics::GetDIAA_AT;
  GetSummariesMap["diaaAV"]=&SummaryStatistics::GetDIAA_AV;
  GetSummariesMap["diaaAW"]=&SummaryStatistics::GetDIAA_AW;
  GetSummariesMap["diaaAY"]=&SummaryStatistics::GetDIAA_AY;
  GetSummariesMap["diaaCA"]=&SummaryStatistics::GetDIAA_CA;
  GetSummariesMap["diaaCC"]=&SummaryStatistics::GetDIAA_CC;
  GetSummariesMap["diaaCD"]=&SummaryStatistics::GetDIAA_CD;
  GetSummariesMap["diaaCE"]=&SummaryStatistics::GetDIAA_CE;
  GetSummariesMap["diaaCF"]=&SummaryStatistics::GetDIAA_CF;
  GetSummariesMap["diaaCG"]=&SummaryStatistics::GetDIAA_CG;
  GetSummariesMap["diaaCH"]=&SummaryStatistics::GetDIAA_CH;
  GetSummariesMap["diaaCI"]=&SummaryStatistics::GetDIAA_CI;
  GetSummariesMap["diaaCK"]=&SummaryStatistics::GetDIAA_CK;
  GetSummariesMap["diaaCL"]=&SummaryStatistics::GetDIAA_CL;
  GetSummariesMap["diaaCM"]=&SummaryStatistics::GetDIAA_CM;
  GetSummariesMap["diaaCN"]=&SummaryStatistics::GetDIAA_CN;
  GetSummariesMap["diaaCP"]=&SummaryStatistics::GetDIAA_CP;
  GetSummariesMap["diaaCQ"]=&SummaryStatistics::GetDIAA_CQ;
  GetSummariesMap["diaaCR"]=&SummaryStatistics::GetDIAA_CR;
  GetSummariesMap["diaaCS"]=&SummaryStatistics::GetDIAA_CS;
  GetSummariesMap["diaaCT"]=&SummaryStatistics::GetDIAA_CT;
  GetSummariesMap["diaaCV"]=&SummaryStatistics::GetDIAA_CV;
  GetSummariesMap["diaaCW"]=&SummaryStatistics::GetDIAA_CW;
  GetSummariesMap["diaaCY"]=&SummaryStatistics::GetDIAA_CY;
  GetSummariesMap["diaaDA"]=&SummaryStatistics::GetDIAA_DA;
  GetSummariesMap["diaaDC"]=&SummaryStatistics::GetDIAA_DC;
  GetSummariesMap["diaaDD"]=&SummaryStatistics::GetDIAA_DD;
  GetSummariesMap["diaaDE"]=&SummaryStatistics::GetDIAA_DE;
  GetSummariesMap["diaaDF"]=&SummaryStatistics::GetDIAA_DF;
  GetSummariesMap["diaaDG"]=&SummaryStatistics::GetDIAA_DG;
  GetSummariesMap["diaaDH"]=&SummaryStatistics::GetDIAA_DH;
  GetSummariesMap["diaaDI"]=&SummaryStatistics::GetDIAA_DI;
  GetSummariesMap["diaaDK"]=&SummaryStatistics::GetDIAA_DK;
  GetSummariesMap["diaaDL"]=&SummaryStatistics::GetDIAA_DL;
  GetSummariesMap["diaaDM"]=&SummaryStatistics::GetDIAA_DM;
  GetSummariesMap["diaaDN"]=&SummaryStatistics::GetDIAA_DN;
  GetSummariesMap["diaaDP"]=&SummaryStatistics::GetDIAA_DP;
  GetSummariesMap["diaaDQ"]=&SummaryStatistics::GetDIAA_DQ;
  GetSummariesMap["diaaDR"]=&SummaryStatistics::GetDIAA_DR;
  GetSummariesMap["diaaDS"]=&SummaryStatistics::GetDIAA_DS;
  GetSummariesMap["diaaDT"]=&SummaryStatistics::GetDIAA_DT;
  GetSummariesMap["diaaDV"]=&SummaryStatistics::GetDIAA_DV;
  GetSummariesMap["diaaDW"]=&SummaryStatistics::GetDIAA_DW;
  GetSummariesMap["diaaDY"]=&SummaryStatistics::GetDIAA_DY;
  GetSummariesMap["diaaEA"]=&SummaryStatistics::GetDIAA_EA;
  GetSummariesMap["diaaEC"]=&SummaryStatistics::GetDIAA_EC;
  GetSummariesMap["diaaED"]=&SummaryStatistics::GetDIAA_ED;
  GetSummariesMap["diaaEE"]=&SummaryStatistics::GetDIAA_EE;
  GetSummariesMap["diaaEF"]=&SummaryStatistics::GetDIAA_EF;
  GetSummariesMap["diaaEG"]=&SummaryStatistics::GetDIAA_EG;
  GetSummariesMap["diaaEH"]=&SummaryStatistics::GetDIAA_EH;
  GetSummariesMap["diaaEI"]=&SummaryStatistics::GetDIAA_EI;
  GetSummariesMap["diaaEK"]=&SummaryStatistics::GetDIAA_EK;
  GetSummariesMap["diaaEL"]=&SummaryStatistics::GetDIAA_EL;
  GetSummariesMap["diaaEM"]=&SummaryStatistics::GetDIAA_EM;
  GetSummariesMap["diaaEN"]=&SummaryStatistics::GetDIAA_EN;
  GetSummariesMap["diaaEP"]=&SummaryStatistics::GetDIAA_EP;
  GetSummariesMap["diaaEQ"]=&SummaryStatistics::GetDIAA_EQ;
  GetSummariesMap["diaaER"]=&SummaryStatistics::GetDIAA_ER;
  GetSummariesMap["diaaES"]=&SummaryStatistics::GetDIAA_ES;
  GetSummariesMap["diaaET"]=&SummaryStatistics::GetDIAA_ET;
  GetSummariesMap["diaaEV"]=&SummaryStatistics::GetDIAA_EV;
  GetSummariesMap["diaaEW"]=&SummaryStatistics::GetDIAA_EW;
  GetSummariesMap["diaaEY"]=&SummaryStatistics::GetDIAA_EY;
  GetSummariesMap["diaaFA"]=&SummaryStatistics::GetDIAA_FA;
  GetSummariesMap["diaaFC"]=&SummaryStatistics::GetDIAA_FC;
  GetSummariesMap["diaaFD"]=&SummaryStatistics::GetDIAA_FD;
  GetSummariesMap["diaaFE"]=&SummaryStatistics::GetDIAA_FE;
  GetSummariesMap["diaaFF"]=&SummaryStatistics::GetDIAA_FF;
  GetSummariesMap["diaaFG"]=&SummaryStatistics::GetDIAA_FG;
  GetSummariesMap["diaaFH"]=&SummaryStatistics::GetDIAA_FH;
  GetSummariesMap["diaaFI"]=&SummaryStatistics::GetDIAA_FI;
  GetSummariesMap["diaaFK"]=&SummaryStatistics::GetDIAA_FK;
  GetSummariesMap["diaaFL"]=&SummaryStatistics::GetDIAA_FL;
  GetSummariesMap["diaaFM"]=&SummaryStatistics::GetDIAA_FM;
  GetSummariesMap["diaaFN"]=&SummaryStatistics::GetDIAA_FN;
  GetSummariesMap["diaaFP"]=&SummaryStatistics::GetDIAA_FP;
  GetSummariesMap["diaaFQ"]=&SummaryStatistics::GetDIAA_FQ;
  GetSummariesMap["diaaFR"]=&SummaryStatistics::GetDIAA_FR;
  GetSummariesMap["diaaFS"]=&SummaryStatistics::GetDIAA_FS;
  GetSummariesMap["diaaFT"]=&SummaryStatistics::GetDIAA_FT;
  GetSummariesMap["diaaFV"]=&SummaryStatistics::GetDIAA_FV;
  GetSummariesMap["diaaFW"]=&SummaryStatistics::GetDIAA_FW;
  GetSummariesMap["diaaFY"]=&SummaryStatistics::GetDIAA_FY;
  GetSummariesMap["diaaGA"]=&SummaryStatistics::GetDIAA_GA;
  GetSummariesMap["diaaGC"]=&SummaryStatistics::GetDIAA_GC;
  GetSummariesMap["diaaGD"]=&SummaryStatistics::GetDIAA_GD;
  GetSummariesMap["diaaGE"]=&SummaryStatistics::GetDIAA_GE;
  GetSummariesMap["diaaGF"]=&SummaryStatistics::GetDIAA_GF;
  GetSummariesMap["diaaGG"]=&SummaryStatistics::GetDIAA_GG;
  GetSummariesMap["diaaGH"]=&SummaryStatistics::GetDIAA_GH;
  GetSummariesMap["diaaGI"]=&SummaryStatistics::GetDIAA_GI;
  GetSummariesMap["diaaGK"]=&SummaryStatistics::GetDIAA_GK;
  GetSummariesMap["diaaGL"]=&SummaryStatistics::GetDIAA_GL;
  GetSummariesMap["diaaGM"]=&SummaryStatistics::GetDIAA_GM;
  GetSummariesMap["diaaGN"]=&SummaryStatistics::GetDIAA_GN;
  GetSummariesMap["diaaGP"]=&SummaryStatistics::GetDIAA_GP;
  GetSummariesMap["diaaGQ"]=&SummaryStatistics::GetDIAA_GQ;
  GetSummariesMap["diaaGR"]=&SummaryStatistics::GetDIAA_GR;
  GetSummariesMap["diaaGS"]=&SummaryStatistics::GetDIAA_GS;
  GetSummariesMap["diaaGT"]=&SummaryStatistics::GetDIAA_GT;
  GetSummariesMap["diaaGV"]=&SummaryStatistics::GetDIAA_GV;
  GetSummariesMap["diaaGW"]=&SummaryStatistics::GetDIAA_GW;
  GetSummariesMap["diaaGY"]=&SummaryStatistics::GetDIAA_GY;
  GetSummariesMap["diaaHA"]=&SummaryStatistics::GetDIAA_HA;
  GetSummariesMap["diaaHC"]=&SummaryStatistics::GetDIAA_HC;
  GetSummariesMap["diaaHD"]=&SummaryStatistics::GetDIAA_HD;
  GetSummariesMap["diaaHE"]=&SummaryStatistics::GetDIAA_HE;
  GetSummariesMap["diaaHF"]=&SummaryStatistics::GetDIAA_HF;
  GetSummariesMap["diaaHG"]=&SummaryStatistics::GetDIAA_HG;
  GetSummariesMap["diaaHH"]=&SummaryStatistics::GetDIAA_HH;
  GetSummariesMap["diaaHI"]=&SummaryStatistics::GetDIAA_HI;
  GetSummariesMap["diaaHK"]=&SummaryStatistics::GetDIAA_HK;
  GetSummariesMap["diaaHL"]=&SummaryStatistics::GetDIAA_HL;
  GetSummariesMap["diaaHM"]=&SummaryStatistics::GetDIAA_HM;
  GetSummariesMap["diaaHN"]=&SummaryStatistics::GetDIAA_HN;
  GetSummariesMap["diaaHP"]=&SummaryStatistics::GetDIAA_HP;
  GetSummariesMap["diaaHQ"]=&SummaryStatistics::GetDIAA_HQ;
  GetSummariesMap["diaaHR"]=&SummaryStatistics::GetDIAA_HR;
  GetSummariesMap["diaaHS"]=&SummaryStatistics::GetDIAA_HS;
  GetSummariesMap["diaaHT"]=&SummaryStatistics::GetDIAA_HT;
  GetSummariesMap["diaaHV"]=&SummaryStatistics::GetDIAA_HV;
  GetSummariesMap["diaaHW"]=&SummaryStatistics::GetDIAA_HW;
  GetSummariesMap["diaaHY"]=&SummaryStatistics::GetDIAA_HY;
  GetSummariesMap["diaaIA"]=&SummaryStatistics::GetDIAA_IA;
  GetSummariesMap["diaaIC"]=&SummaryStatistics::GetDIAA_IC;
  GetSummariesMap["diaaID"]=&SummaryStatistics::GetDIAA_ID;
  GetSummariesMap["diaaIE"]=&SummaryStatistics::GetDIAA_IE;
  GetSummariesMap["diaaIF"]=&SummaryStatistics::GetDIAA_IF;
  GetSummariesMap["diaaIG"]=&SummaryStatistics::GetDIAA_IG;
  GetSummariesMap["diaaIH"]=&SummaryStatistics::GetDIAA_IH;
  GetSummariesMap["diaaII"]=&SummaryStatistics::GetDIAA_II;
  GetSummariesMap["diaaIK"]=&SummaryStatistics::GetDIAA_IK;
  GetSummariesMap["diaaIL"]=&SummaryStatistics::GetDIAA_IL;
  GetSummariesMap["diaaIM"]=&SummaryStatistics::GetDIAA_IM;
  GetSummariesMap["diaaIN"]=&SummaryStatistics::GetDIAA_IN;
  GetSummariesMap["diaaIP"]=&SummaryStatistics::GetDIAA_IP;
  GetSummariesMap["diaaIQ"]=&SummaryStatistics::GetDIAA_IQ;
  GetSummariesMap["diaaIR"]=&SummaryStatistics::GetDIAA_IR;
  GetSummariesMap["diaaIS"]=&SummaryStatistics::GetDIAA_IS;
  GetSummariesMap["diaaIT"]=&SummaryStatistics::GetDIAA_IT;
  GetSummariesMap["diaaIV"]=&SummaryStatistics::GetDIAA_IV;
  GetSummariesMap["diaaIW"]=&SummaryStatistics::GetDIAA_IW;
  GetSummariesMap["diaaIY"]=&SummaryStatistics::GetDIAA_IY;
  GetSummariesMap["diaaKA"]=&SummaryStatistics::GetDIAA_KA;
  GetSummariesMap["diaaKC"]=&SummaryStatistics::GetDIAA_KC;
  GetSummariesMap["diaaKD"]=&SummaryStatistics::GetDIAA_KD;
  GetSummariesMap["diaaKE"]=&SummaryStatistics::GetDIAA_KE;
  GetSummariesMap["diaaKF"]=&SummaryStatistics::GetDIAA_KF;
  GetSummariesMap["diaaKG"]=&SummaryStatistics::GetDIAA_KG;
  GetSummariesMap["diaaKH"]=&SummaryStatistics::GetDIAA_KH;
  GetSummariesMap["diaaKI"]=&SummaryStatistics::GetDIAA_KI;
  GetSummariesMap["diaaKK"]=&SummaryStatistics::GetDIAA_KK;
  GetSummariesMap["diaaKL"]=&SummaryStatistics::GetDIAA_KL;
  GetSummariesMap["diaaKM"]=&SummaryStatistics::GetDIAA_KM;
  GetSummariesMap["diaaKN"]=&SummaryStatistics::GetDIAA_KN;
  GetSummariesMap["diaaKP"]=&SummaryStatistics::GetDIAA_KP;
  GetSummariesMap["diaaKQ"]=&SummaryStatistics::GetDIAA_KQ;
  GetSummariesMap["diaaKR"]=&SummaryStatistics::GetDIAA_KR;
  GetSummariesMap["diaaKS"]=&SummaryStatistics::GetDIAA_KS;
  GetSummariesMap["diaaKT"]=&SummaryStatistics::GetDIAA_KT;
  GetSummariesMap["diaaKV"]=&SummaryStatistics::GetDIAA_KV;
  GetSummariesMap["diaaKW"]=&SummaryStatistics::GetDIAA_KW;
  GetSummariesMap["diaaKY"]=&SummaryStatistics::GetDIAA_KY;
  GetSummariesMap["diaaLA"]=&SummaryStatistics::GetDIAA_LA;
  GetSummariesMap["diaaLC"]=&SummaryStatistics::GetDIAA_LC;
  GetSummariesMap["diaaLD"]=&SummaryStatistics::GetDIAA_LD;
  GetSummariesMap["diaaLE"]=&SummaryStatistics::GetDIAA_LE;
  GetSummariesMap["diaaLF"]=&SummaryStatistics::GetDIAA_LF;
  GetSummariesMap["diaaLG"]=&SummaryStatistics::GetDIAA_LG;
  GetSummariesMap["diaaLH"]=&SummaryStatistics::GetDIAA_LH;
  GetSummariesMap["diaaLI"]=&SummaryStatistics::GetDIAA_LI;
  GetSummariesMap["diaaLK"]=&SummaryStatistics::GetDIAA_LK;
  GetSummariesMap["diaaLL"]=&SummaryStatistics::GetDIAA_LL;
  GetSummariesMap["diaaLM"]=&SummaryStatistics::GetDIAA_LM;
  GetSummariesMap["diaaLN"]=&SummaryStatistics::GetDIAA_LN;
  GetSummariesMap["diaaLP"]=&SummaryStatistics::GetDIAA_LP;
  GetSummariesMap["diaaLQ"]=&SummaryStatistics::GetDIAA_LQ;
  GetSummariesMap["diaaLR"]=&SummaryStatistics::GetDIAA_LR;
  GetSummariesMap["diaaLS"]=&SummaryStatistics::GetDIAA_LS;
  GetSummariesMap["diaaLT"]=&SummaryStatistics::GetDIAA_LT;
  GetSummariesMap["diaaLV"]=&SummaryStatistics::GetDIAA_LV;
  GetSummariesMap["diaaLW"]=&SummaryStatistics::GetDIAA_LW;
  GetSummariesMap["diaaLY"]=&SummaryStatistics::GetDIAA_LY;
  GetSummariesMap["diaaMA"]=&SummaryStatistics::GetDIAA_MA;
  GetSummariesMap["diaaMC"]=&SummaryStatistics::GetDIAA_MC;
  GetSummariesMap["diaaMD"]=&SummaryStatistics::GetDIAA_MD;
  GetSummariesMap["diaaME"]=&SummaryStatistics::GetDIAA_ME;
  GetSummariesMap["diaaMF"]=&SummaryStatistics::GetDIAA_MF;
  GetSummariesMap["diaaMG"]=&SummaryStatistics::GetDIAA_MG;
  GetSummariesMap["diaaMH"]=&SummaryStatistics::GetDIAA_MH;
  GetSummariesMap["diaaMI"]=&SummaryStatistics::GetDIAA_MI;
  GetSummariesMap["diaaMK"]=&SummaryStatistics::GetDIAA_MK;
  GetSummariesMap["diaaML"]=&SummaryStatistics::GetDIAA_ML;
  GetSummariesMap["diaaMM"]=&SummaryStatistics::GetDIAA_MM;
  GetSummariesMap["diaaMN"]=&SummaryStatistics::GetDIAA_MN;
  GetSummariesMap["diaaMP"]=&SummaryStatistics::GetDIAA_MP;
  GetSummariesMap["diaaMQ"]=&SummaryStatistics::GetDIAA_MQ;
  GetSummariesMap["diaaMR"]=&SummaryStatistics::GetDIAA_MR;
  GetSummariesMap["diaaMS"]=&SummaryStatistics::GetDIAA_MS;
  GetSummariesMap["diaaMT"]=&SummaryStatistics::GetDIAA_MT;
  GetSummariesMap["diaaMV"]=&SummaryStatistics::GetDIAA_MV;
  GetSummariesMap["diaaMW"]=&SummaryStatistics::GetDIAA_MW;
  GetSummariesMap["diaaMY"]=&SummaryStatistics::GetDIAA_MY;
  GetSummariesMap["diaaNA"]=&SummaryStatistics::GetDIAA_NA;
  GetSummariesMap["diaaNC"]=&SummaryStatistics::GetDIAA_NC;
  GetSummariesMap["diaaND"]=&SummaryStatistics::GetDIAA_ND;
  GetSummariesMap["diaaNE"]=&SummaryStatistics::GetDIAA_NE;
  GetSummariesMap["diaaNF"]=&SummaryStatistics::GetDIAA_NF;
  GetSummariesMap["diaaNG"]=&SummaryStatistics::GetDIAA_NG;
  GetSummariesMap["diaaNH"]=&SummaryStatistics::GetDIAA_NH;
  GetSummariesMap["diaaNI"]=&SummaryStatistics::GetDIAA_NI;
  GetSummariesMap["diaaNK"]=&SummaryStatistics::GetDIAA_NK;
  GetSummariesMap["diaaNL"]=&SummaryStatistics::GetDIAA_NL;
  GetSummariesMap["diaaNM"]=&SummaryStatistics::GetDIAA_NM;
  GetSummariesMap["diaaNN"]=&SummaryStatistics::GetDIAA_NN;
  GetSummariesMap["diaaNP"]=&SummaryStatistics::GetDIAA_NP;
  GetSummariesMap["diaaNQ"]=&SummaryStatistics::GetDIAA_NQ;
  GetSummariesMap["diaaNR"]=&SummaryStatistics::GetDIAA_NR;
  GetSummariesMap["diaaNS"]=&SummaryStatistics::GetDIAA_NS;
  GetSummariesMap["diaaNT"]=&SummaryStatistics::GetDIAA_NT;
  GetSummariesMap["diaaNV"]=&SummaryStatistics::GetDIAA_NV;
  GetSummariesMap["diaaNW"]=&SummaryStatistics::GetDIAA_NW;
  GetSummariesMap["diaaNY"]=&SummaryStatistics::GetDIAA_NY;
  GetSummariesMap["diaaPA"]=&SummaryStatistics::GetDIAA_PA;
  GetSummariesMap["diaaPC"]=&SummaryStatistics::GetDIAA_PC;
  GetSummariesMap["diaaPD"]=&SummaryStatistics::GetDIAA_PD;
  GetSummariesMap["diaaPE"]=&SummaryStatistics::GetDIAA_PE;
  GetSummariesMap["diaaPF"]=&SummaryStatistics::GetDIAA_PF;
  GetSummariesMap["diaaPG"]=&SummaryStatistics::GetDIAA_PG;
  GetSummariesMap["diaaPH"]=&SummaryStatistics::GetDIAA_PH;
  GetSummariesMap["diaaPI"]=&SummaryStatistics::GetDIAA_PI;
  GetSummariesMap["diaaPK"]=&SummaryStatistics::GetDIAA_PK;
  GetSummariesMap["diaaPL"]=&SummaryStatistics::GetDIAA_PL;
  GetSummariesMap["diaaPM"]=&SummaryStatistics::GetDIAA_PM;
  GetSummariesMap["diaaPN"]=&SummaryStatistics::GetDIAA_PN;
  GetSummariesMap["diaaPP"]=&SummaryStatistics::GetDIAA_PP;
  GetSummariesMap["diaaPQ"]=&SummaryStatistics::GetDIAA_PQ;
  GetSummariesMap["diaaPR"]=&SummaryStatistics::GetDIAA_PR;
  GetSummariesMap["diaaPS"]=&SummaryStatistics::GetDIAA_PS;
  GetSummariesMap["diaaPT"]=&SummaryStatistics::GetDIAA_PT;
  GetSummariesMap["diaaPV"]=&SummaryStatistics::GetDIAA_PV;
  GetSummariesMap["diaaPW"]=&SummaryStatistics::GetDIAA_PW;
  GetSummariesMap["diaaPY"]=&SummaryStatistics::GetDIAA_PY;
  GetSummariesMap["diaaQA"]=&SummaryStatistics::GetDIAA_QA;
  GetSummariesMap["diaaQC"]=&SummaryStatistics::GetDIAA_QC;
  GetSummariesMap["diaaQD"]=&SummaryStatistics::GetDIAA_QD;
  GetSummariesMap["diaaQE"]=&SummaryStatistics::GetDIAA_QE;
  GetSummariesMap["diaaQF"]=&SummaryStatistics::GetDIAA_QF;
  GetSummariesMap["diaaQG"]=&SummaryStatistics::GetDIAA_QG;
  GetSummariesMap["diaaQH"]=&SummaryStatistics::GetDIAA_QH;
  GetSummariesMap["diaaQI"]=&SummaryStatistics::GetDIAA_QI;
  GetSummariesMap["diaaQK"]=&SummaryStatistics::GetDIAA_QK;
  GetSummariesMap["diaaQL"]=&SummaryStatistics::GetDIAA_QL;
  GetSummariesMap["diaaQM"]=&SummaryStatistics::GetDIAA_QM;
  GetSummariesMap["diaaQN"]=&SummaryStatistics::GetDIAA_QN;
  GetSummariesMap["diaaQP"]=&SummaryStatistics::GetDIAA_QP;
  GetSummariesMap["diaaQQ"]=&SummaryStatistics::GetDIAA_QQ;
  GetSummariesMap["diaaQR"]=&SummaryStatistics::GetDIAA_QR;
  GetSummariesMap["diaaQS"]=&SummaryStatistics::GetDIAA_QS;
  GetSummariesMap["diaaQT"]=&SummaryStatistics::GetDIAA_QT;
  GetSummariesMap["diaaQV"]=&SummaryStatistics::GetDIAA_QV;
  GetSummariesMap["diaaQW"]=&SummaryStatistics::GetDIAA_QW;
  GetSummariesMap["diaaQY"]=&SummaryStatistics::GetDIAA_QY;
  GetSummariesMap["diaaRA"]=&SummaryStatistics::GetDIAA_RA;
  GetSummariesMap["diaaRC"]=&SummaryStatistics::GetDIAA_RC;
  GetSummariesMap["diaaRD"]=&SummaryStatistics::GetDIAA_RD;
  GetSummariesMap["diaaRE"]=&SummaryStatistics::GetDIAA_RE;
  GetSummariesMap["diaaRF"]=&SummaryStatistics::GetDIAA_RF;
  GetSummariesMap["diaaRG"]=&SummaryStatistics::GetDIAA_RG;
  GetSummariesMap["diaaRH"]=&SummaryStatistics::GetDIAA_RH;
  GetSummariesMap["diaaRI"]=&SummaryStatistics::GetDIAA_RI;
  GetSummariesMap["diaaRK"]=&SummaryStatistics::GetDIAA_RK;
  GetSummariesMap["diaaRL"]=&SummaryStatistics::GetDIAA_RL;
  GetSummariesMap["diaaRM"]=&SummaryStatistics::GetDIAA_RM;
  GetSummariesMap["diaaRN"]=&SummaryStatistics::GetDIAA_RN;
  GetSummariesMap["diaaRP"]=&SummaryStatistics::GetDIAA_RP;
  GetSummariesMap["diaaRQ"]=&SummaryStatistics::GetDIAA_RQ;
  GetSummariesMap["diaaRR"]=&SummaryStatistics::GetDIAA_RR;
  GetSummariesMap["diaaRS"]=&SummaryStatistics::GetDIAA_RS;
  GetSummariesMap["diaaRT"]=&SummaryStatistics::GetDIAA_RT;
  GetSummariesMap["diaaRV"]=&SummaryStatistics::GetDIAA_RV;
  GetSummariesMap["diaaRW"]=&SummaryStatistics::GetDIAA_RW;
  GetSummariesMap["diaaRY"]=&SummaryStatistics::GetDIAA_RY;
  GetSummariesMap["diaaSA"]=&SummaryStatistics::GetDIAA_SA;
  GetSummariesMap["diaaSC"]=&SummaryStatistics::GetDIAA_SC;
  GetSummariesMap["diaaSD"]=&SummaryStatistics::GetDIAA_SD;
  GetSummariesMap["diaaSE"]=&SummaryStatistics::GetDIAA_SE;
  GetSummariesMap["diaaSF"]=&SummaryStatistics::GetDIAA_SF;
  GetSummariesMap["diaaSG"]=&SummaryStatistics::GetDIAA_SG;
  GetSummariesMap["diaaSH"]=&SummaryStatistics::GetDIAA_SH;
  GetSummariesMap["diaaSI"]=&SummaryStatistics::GetDIAA_SI;
  GetSummariesMap["diaaSK"]=&SummaryStatistics::GetDIAA_SK;
  GetSummariesMap["diaaSL"]=&SummaryStatistics::GetDIAA_SL;
  GetSummariesMap["diaaSM"]=&SummaryStatistics::GetDIAA_SM;
  GetSummariesMap["diaaSN"]=&SummaryStatistics::GetDIAA_SN;
  GetSummariesMap["diaaSP"]=&SummaryStatistics::GetDIAA_SP;
  GetSummariesMap["diaaSQ"]=&SummaryStatistics::GetDIAA_SQ;
  GetSummariesMap["diaaSR"]=&SummaryStatistics::GetDIAA_SR;
  GetSummariesMap["diaaSS"]=&SummaryStatistics::GetDIAA_SS;
  GetSummariesMap["diaaST"]=&SummaryStatistics::GetDIAA_ST;
  GetSummariesMap["diaaSV"]=&SummaryStatistics::GetDIAA_SV;
  GetSummariesMap["diaaSW"]=&SummaryStatistics::GetDIAA_SW;
  GetSummariesMap["diaaSY"]=&SummaryStatistics::GetDIAA_SY;
  GetSummariesMap["diaaTA"]=&SummaryStatistics::GetDIAA_TA;
  GetSummariesMap["diaaTC"]=&SummaryStatistics::GetDIAA_TC;
  GetSummariesMap["diaaTD"]=&SummaryStatistics::GetDIAA_TD;
  GetSummariesMap["diaaTE"]=&SummaryStatistics::GetDIAA_TE;
  GetSummariesMap["diaaTF"]=&SummaryStatistics::GetDIAA_TF;
  GetSummariesMap["diaaTG"]=&SummaryStatistics::GetDIAA_TG;
  GetSummariesMap["diaaTH"]=&SummaryStatistics::GetDIAA_TH;
  GetSummariesMap["diaaTI"]=&SummaryStatistics::GetDIAA_TI;
  GetSummariesMap["diaaTK"]=&SummaryStatistics::GetDIAA_TK;
  GetSummariesMap["diaaTL"]=&SummaryStatistics::GetDIAA_TL;
  GetSummariesMap["diaaTM"]=&SummaryStatistics::GetDIAA_TM;
  GetSummariesMap["diaaTN"]=&SummaryStatistics::GetDIAA_TN;
  GetSummariesMap["diaaTP"]=&SummaryStatistics::GetDIAA_TP;
  GetSummariesMap["diaaTQ"]=&SummaryStatistics::GetDIAA_TQ;
  GetSummariesMap["diaaTR"]=&SummaryStatistics::GetDIAA_TR;
  GetSummariesMap["diaaTS"]=&SummaryStatistics::GetDIAA_TS;
  GetSummariesMap["diaaTT"]=&SummaryStatistics::GetDIAA_TT;
  GetSummariesMap["diaaTV"]=&SummaryStatistics::GetDIAA_TV;
  GetSummariesMap["diaaTW"]=&SummaryStatistics::GetDIAA_TW;
  GetSummariesMap["diaaTY"]=&SummaryStatistics::GetDIAA_TY;
  GetSummariesMap["diaaVA"]=&SummaryStatistics::GetDIAA_VA;
  GetSummariesMap["diaaVC"]=&SummaryStatistics::GetDIAA_VC;
  GetSummariesMap["diaaVD"]=&SummaryStatistics::GetDIAA_VD;
  GetSummariesMap["diaaVE"]=&SummaryStatistics::GetDIAA_VE;
  GetSummariesMap["diaaVF"]=&SummaryStatistics::GetDIAA_VF;
  GetSummariesMap["diaaVG"]=&SummaryStatistics::GetDIAA_VG;
  GetSummariesMap["diaaVH"]=&SummaryStatistics::GetDIAA_VH;
  GetSummariesMap["diaaVI"]=&SummaryStatistics::GetDIAA_VI;
  GetSummariesMap["diaaVK"]=&SummaryStatistics::GetDIAA_VK;
  GetSummariesMap["diaaVL"]=&SummaryStatistics::GetDIAA_VL;
  GetSummariesMap["diaaVM"]=&SummaryStatistics::GetDIAA_VM;
  GetSummariesMap["diaaVN"]=&SummaryStatistics::GetDIAA_VN;
  GetSummariesMap["diaaVP"]=&SummaryStatistics::GetDIAA_VP;
  GetSummariesMap["diaaVQ"]=&SummaryStatistics::GetDIAA_VQ;
  GetSummariesMap["diaaVR"]=&SummaryStatistics::GetDIAA_VR;
  GetSummariesMap["diaaVS"]=&SummaryStatistics::GetDIAA_VS;
  GetSummariesMap["diaaVT"]=&SummaryStatistics::GetDIAA_VT;
  GetSummariesMap["diaaVV"]=&SummaryStatistics::GetDIAA_VV;
  GetSummariesMap["diaaVW"]=&SummaryStatistics::GetDIAA_VW;
  GetSummariesMap["diaaVY"]=&SummaryStatistics::GetDIAA_VY;
  GetSummariesMap["diaaWA"]=&SummaryStatistics::GetDIAA_WA;
  GetSummariesMap["diaaWC"]=&SummaryStatistics::GetDIAA_WC;
  GetSummariesMap["diaaWD"]=&SummaryStatistics::GetDIAA_WD;
  GetSummariesMap["diaaWE"]=&SummaryStatistics::GetDIAA_WE;
  GetSummariesMap["diaaWF"]=&SummaryStatistics::GetDIAA_WF;
  GetSummariesMap["diaaWG"]=&SummaryStatistics::GetDIAA_WG;
  GetSummariesMap["diaaWH"]=&SummaryStatistics::GetDIAA_WH;
  GetSummariesMap["diaaWI"]=&SummaryStatistics::GetDIAA_WI;
  GetSummariesMap["diaaWK"]=&SummaryStatistics::GetDIAA_WK;
  GetSummariesMap["diaaWL"]=&SummaryStatistics::GetDIAA_WL;
  GetSummariesMap["diaaWM"]=&SummaryStatistics::GetDIAA_WM;
  GetSummariesMap["diaaWN"]=&SummaryStatistics::GetDIAA_WN;
  GetSummariesMap["diaaWP"]=&SummaryStatistics::GetDIAA_WP;
  GetSummariesMap["diaaWQ"]=&SummaryStatistics::GetDIAA_WQ;
  GetSummariesMap["diaaWR"]=&SummaryStatistics::GetDIAA_WR;
  GetSummariesMap["diaaWS"]=&SummaryStatistics::GetDIAA_WS;
  GetSummariesMap["diaaWT"]=&SummaryStatistics::GetDIAA_WT;
  GetSummariesMap["diaaWV"]=&SummaryStatistics::GetDIAA_WV;
  GetSummariesMap["diaaWW"]=&SummaryStatistics::GetDIAA_WW;
  GetSummariesMap["diaaWY"]=&SummaryStatistics::GetDIAA_WY;
  GetSummariesMap["diaaYA"]=&SummaryStatistics::GetDIAA_YA;
  GetSummariesMap["diaaYC"]=&SummaryStatistics::GetDIAA_YC;
  GetSummariesMap["diaaYD"]=&SummaryStatistics::GetDIAA_YD;
  GetSummariesMap["diaaYE"]=&SummaryStatistics::GetDIAA_YE;
  GetSummariesMap["diaaYF"]=&SummaryStatistics::GetDIAA_YF;
  GetSummariesMap["diaaYG"]=&SummaryStatistics::GetDIAA_YG;
  GetSummariesMap["diaaYH"]=&SummaryStatistics::GetDIAA_YH;
  GetSummariesMap["diaaYI"]=&SummaryStatistics::GetDIAA_YI;
  GetSummariesMap["diaaYK"]=&SummaryStatistics::GetDIAA_YK;
  GetSummariesMap["diaaYL"]=&SummaryStatistics::GetDIAA_YL;
  GetSummariesMap["diaaYM"]=&SummaryStatistics::GetDIAA_YM;
  GetSummariesMap["diaaYN"]=&SummaryStatistics::GetDIAA_YN;
  GetSummariesMap["diaaYP"]=&SummaryStatistics::GetDIAA_YP;
  GetSummariesMap["diaaYQ"]=&SummaryStatistics::GetDIAA_YQ;
  GetSummariesMap["diaaYR"]=&SummaryStatistics::GetDIAA_YR;
  GetSummariesMap["diaaYS"]=&SummaryStatistics::GetDIAA_YS;
  GetSummariesMap["diaaYT"]=&SummaryStatistics::GetDIAA_YT;
  GetSummariesMap["diaaYV"]=&SummaryStatistics::GetDIAA_YV;
  GetSummariesMap["diaaYW"]=&SummaryStatistics::GetDIAA_YW;
  GetSummariesMap["diaaYY"]=&SummaryStatistics::GetDIAA_YY; */
  GetSummariesMap["RSCUentropy"] = &SummaryStatistics::GetRSCUentropy;
  GetSummariesMap["GC"] = &SummaryStatistics::GetGC;
  GetSummariesMap["GC1"] = &SummaryStatistics::GetGC1;
  GetSummariesMap["GC2"] = &SummaryStatistics::GetGC2;
  GetSummariesMap["GC3"] = &SummaryStatistics::GetGC3;
  GetSummariesMap["Codonfentropy"] =
      &SummaryStatistics::GetRelativeCodonFrequencyEntropy;
  GetSummariesMap["AAentropy"] = &SummaryStatistics::GetAAentropy;
  GetSummariesMap["Dinucentropy"] = &SummaryStatistics::GetDinucentropy;
  GetSummariesMap["Dinuc12entropy"] = &SummaryStatistics::GetDinuc12entropy;
  GetSummariesMap["Dinuc23entropy"] = &SummaryStatistics::GetDinuc23entropy;
  GetSummariesMap["Dinuc31entropy"] = &SummaryStatistics::GetDinuc31entropy;
}

SummaryStatistics::~SummaryStatistics() {
  for (int nuc = 0; nuc < lparam->Nnucp; nuc++) {
    delete[] dinuc_usage[nuc];
    delete[] dinuc12_usage[nuc];
    delete[] dinuc23_usage[nuc];
    delete[] dinuc31_usage[nuc];
  }

  delete[] dinuc_usage;
  delete[] dinuc12_usage;
  delete[] dinuc23_usage;
  delete[] dinuc31_usage;
  delete[] dicodon_usage;
  delete[] diaa_usage;
  delete[] relativeAAFrequency;
  delete[] nuc_usage;
  delete[] nuc1_usage;
  delete[] nuc2_usage;
  delete[] nuc3_usage;
  delete[] nuc_meandiff;
  delete[] nuc1_meandiff;
  delete[] nuc2_meandiff;
  delete[] nuc3_meandiff;
  delete[] codon_meandiff;
  delete[] aa_meandiff;
  delete[] nuc_pairwise;
  delete[] nuc_pairwise10;
  delete[] nuc_pairwise30;
  delete[] nuc_pairwise50;
  delete[] nuc_pairwise70;
  delete[] nuc_pairwise90;
  delete[] nuc1_pairwise;
  delete[] nuc1_pairwise10;
  delete[] nuc1_pairwise30;
  delete[] nuc1_pairwise50;
  delete[] nuc1_pairwise70;
  delete[] nuc1_pairwise90;
  delete[] nuc2_pairwise;
  delete[] nuc2_pairwise10;
  delete[] nuc2_pairwise30;
  delete[] nuc2_pairwise50;
  delete[] nuc2_pairwise70;
  delete[] nuc2_pairwise90;
  delete[] nuc3_pairwise;
  delete[] nuc3_pairwise10;
  delete[] nuc3_pairwise30;
  delete[] nuc3_pairwise50;
  delete[] nuc3_pairwise70;
  delete[] nuc3_pairwise90;
  delete[] aa_pairwise;
  delete[] aa_pairwise10;
  delete[] aa_pairwise30;
  delete[] aa_pairwise50;
  delete[] aa_pairwise70;
  delete[] aa_pairwise90;
}

// void SummaryStatistics::computeSummaries(CodonSequenceAlignment* codondata) {

//}

void SummaryStatistics::computeSummariesAncestralSequence(
    int** CurrentAncestralCodonSequence) {
  lparam->summariesAncestralData.clear();
  lparam->summariesAncestralData.shrink_to_fit();

  RSCU_bool = false;
  relativeCodonFrequency_bool = false;
  dinuc_bool = false;
  dinuc12_bool = false;
  dinuc23_bool = false;
  dinuc31_bool = false;
  relativeAAFrequency_bool = false;
  dicodon_bool = false;
  diaa_bool = false;
  nuc_bool = false;
  nuc1_bool = false;
  nuc2_bool = false;
  nuc3_bool = false;
  nuc_meandiff_bool = false;
  nuc1_meandiff_bool = false;
  nuc2_meandiff_bool = false;
  nuc3_meandiff_bool = false;
  codon_meandiff_bool = false;
  aa_meandiff_bool = false;
  nuc_pairwise_bool = false;
  nuc_pairwise_bool10 = false;
  nuc_pairwise_bool30 = false;
  nuc_pairwise_bool50 = false;
  nuc_pairwise_bool90 = false;
  nuc_pairwiseSyn_bool = false;
  nuc_pairwiseSyn_bool10 = false;
  nuc_pairwiseSyn_bool30 = false;
  nuc_pairwiseSyn_bool50 = false;
  nuc_pairwiseSyn_bool70 = false;
  nuc_pairwiseSyn_bool90 = false;
  nuc1_pairwise_bool = false;
  nuc2_pairwise_bool = false;
  nuc3_pairwise_bool = false;
  nuc3_pairwise_bool10 = false;
  nuc3_pairwise_bool30 = false;
  nuc3_pairwise_bool50 = false;
  aa_pairwise_bool = false;
  aa_pairwise_bool10 = false;
  aa_pairwise_bool30 = false;
  aa_pairwise_bool50 = false;
  aa_pairwise_bool90 = false;
  dinucCpG_pairwise_bool = false;
  nuc_site_comphet_bool = false;
  nuc1_site_comphet_bool = false;
  nuc2_site_comphet_bool = false;
  nuc3_site_comphet_bool = false;
  nuc_taxa_comphet_bool = false;
  nuc1_taxa_comphet_bool = false;
  nuc2_taxa_comphet_bool = false;
  nuc3_taxa_comphet_bool = false;
  codon_site_comphet_bool = false;
  codon_taxa_comphet_bool = false;
  aa_site_comphet_bool = false;
  aa_taxa_comphet_bool = false;
  RSCUentropy_bool = false;
  GC_bool = false;
  GC1_bool = false;
  GC2_bool = false;
  GC3_bool = false;

  CodonSequenceAlignment* simulation = new CodonSequenceAlignment(
      lparam->codondata, 1, CurrentAncestralCodonSequence);

  string* arrSummaries = new string[lparam->NusedAncSummaries];
  for (int i_summary = 0; i_summary < lparam->NSummaries; i_summary++) {
    auto it =
        lparam->mapUsedAncSummaries.find(lparam->listSummaries[i_summary]);
    if (it != lparam->mapUsedAncSummaries.end() && it->second != -1) {
      arrSummaries[it->second] = it->first;
    }
  }

  for (int i_summary = 0; i_summary < lparam->NusedAncSummaries; i_summary++) {
    auto it = GetSummariesMap.find(arrSummaries[i_summary]);
    if (it != GetSummariesMap.end()) {
      funcpt f = GetSummariesMap[arrSummaries[i_summary]];
      double s = (this->*f)(simulation);

      if (s < lparam->TOOSMALL || isinf(s)) {
        s = lparam->TOOSMALL;
      }

      s = transformSummaryStatistics(s);

      lparam->summariesAncestralData.push_back(s);
    }
  }

  delete[] arrSummaries;
}

void SummaryStatistics::computeSummaries(int** CurrentNodeLeafCodonSequence) {
  int verbose = lparam->verbose;
  if (verbose) {
    std::cerr << "SummaryStatistics::computeSummaries\n";
  }

  int** CurrentLeafNodeCodonSequences_ = new int*[lparam->Ntaxa];
  for (int taxa = 0; taxa < lparam->Ntaxa; taxa++) {
    CurrentLeafNodeCodonSequences_[taxa] = new int[lparam->Nsite_codon];
  }

  for (int taxa = 0; taxa < lparam->Ntaxa; taxa++) {
    for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
      int state = lparam->codondata->GetState(taxa, site_codon);
      if (state != unknown) {
        CurrentLeafNodeCodonSequences_[taxa][site_codon] =
            CurrentNodeLeafCodonSequence[taxa][site_codon];

      } else {
        CurrentLeafNodeCodonSequences_[taxa][site_codon] = -1;
      }
    }
  }

  lparam->summariesSimulatedData.clear();
  lparam->summariesSimulatedData.shrink_to_fit();

  lparam->accessorysummariesSimulatedData.clear();
  lparam->accessorysummariesSimulatedData.shrink_to_fit();
  if (verbose) {
    std::cerr << "clear\n";
  }
  RSCU_bool = false;
  relativeCodonFrequency_bool = false;
  dinuc_bool = false;
  dinuc12_bool = false;
  dinuc23_bool = false;
  dinuc31_bool = false;
  relativeAAFrequency_bool = false;
  dicodon_bool = false;
  diaa_bool = false;
  nuc_bool = false;
  nuc1_bool = false;
  nuc2_bool = false;
  nuc3_bool = false;
  nuc_meandiff_bool = false;
  nuc1_meandiff_bool = false;
  nuc2_meandiff_bool = false;
  nuc3_meandiff_bool = false;
  codon_meandiff_bool = false;
  aa_meandiff_bool = false;
  nuc_pairwise_bool = false;
  nuc_pairwise_bool10 = false;
  nuc_pairwise_bool30 = false;
  nuc_pairwise_bool50 = false;
  nuc_pairwise_bool70 = false;
  nuc_pairwise_bool90 = false;
  nuc_pairwiseSyn_bool = false;
  nuc_pairwiseSyn_bool10 = false;
  nuc_pairwiseSyn_bool30 = false;
  nuc_pairwiseSyn_bool50 = false;
  nuc_pairwiseSyn_bool70 = false;
  nuc_pairwiseSyn_bool90 = false;
  nuc1_pairwise_bool = false;
  nuc1_pairwise_bool10 = false;
  nuc1_pairwise_bool30 = false;
  nuc1_pairwise_bool50 = false;
  nuc1_pairwise_bool70 = false;
  nuc1_pairwise_bool90 = false;
  nuc2_pairwise_bool = false;
  nuc2_pairwise_bool10 = false;
  nuc2_pairwise_bool30 = false;
  nuc2_pairwise_bool50 = false;
  nuc2_pairwise_bool70 = false;
  nuc2_pairwise_bool90 = false;
  nuc3_pairwise_bool = false;
  nuc3_pairwise_bool10 = false;
  nuc3_pairwise_bool30 = false;
  nuc3_pairwise_bool50 = false;
  nuc3_pairwise_bool70 = false;
  nuc3_pairwise_bool90 = false;
  aa_pairwise_bool = false;
  aa_pairwise_bool10 = false;
  aa_pairwise_bool30 = false;
  aa_pairwise_bool50 = false;
  aa_pairwise_bool70 = false;
  aa_pairwise_bool90 = false;
  dinucCpG_pairwise_bool = false;
  nuc_site_comphet_bool = false;
  nuc1_site_comphet_bool = false;
  nuc2_site_comphet_bool = false;
  nuc3_site_comphet_bool = false;
  nuc_taxa_comphet_bool = false;
  nuc1_taxa_comphet_bool = false;
  nuc2_taxa_comphet_bool = false;
  nuc3_taxa_comphet_bool = false;
  codon_site_comphet_bool = false;
  codon_taxa_comphet_bool = false;
  aa_site_comphet_bool = false;
  aa_taxa_comphet_bool = false;
  RSCUentropy_bool = false;
  GC_bool = false;
  GC1_bool = false;
  GC2_bool = false;
  GC3_bool = false;

  CodonSequenceAlignment* simulation = new CodonSequenceAlignment(
      lparam->codondata, CurrentLeafNodeCodonSequences_);

  if (verbose) {
    std::cerr << "computeSummaries(int** CurrentNodeLeafCodonSequence)\n";
  }
  string* arrSummaries = new string[lparam->NusedSummaries];
  for (int i_summary = 0; i_summary < lparam->NSummaries; i_summary++) {
    auto it = lparam->mapUsedSummaries.find(lparam->listSummaries[i_summary]);
    if (it != lparam->mapUsedSummaries.end() && it->second != -1) {
      arrSummaries[it->second] = it->first;
    }
  }

  if (verbose) {
    std::cerr << "computeSummaries(int** CurrentNodeLeafCodonSequence)2\n";
  }
  for (int i_summary = 0; i_summary < lparam->NusedSummaries; i_summary++) {
    auto it = GetSummariesMap.find(arrSummaries[i_summary]);
    if (it != GetSummariesMap.end()) {
      funcpt f = GetSummariesMap[arrSummaries[i_summary]];
      double s = (this->*f)(simulation);

      if (s < lparam->TOOSMALL || isinf(s)) {
        s = lparam->TOOSMALL;
      }

      s = transformSummaryStatistics(s);

      lparam->summariesSimulatedData.push_back(s);
    }
  }
  if (verbose) {
    std::cerr << "computeSummaries(int** CurrentNodeLeafCodonSequence)3\n";
  }
  string* arrAccSummaries = new string[lparam->NusedAccessorySummaries];
  for (int i_summary = 0; i_summary < lparam->NSummaries; i_summary++) {
    auto it = lparam->mapUsedAccessorySummaries.find(
        lparam->listSummaries[i_summary]);
    if (it != lparam->mapUsedAccessorySummaries.end() && it->second != -1) {
      arrAccSummaries[it->second] = it->first;
    }
  }

  if (verbose) {
    std::cerr << "computeSummaries(int** CurrentNodeLeafCodonSequence)4\n";
  }
  for (int i_summary = 0; i_summary < lparam->NusedAccessorySummaries;
       i_summary++) {
    auto it = GetSummariesMap.find(arrAccSummaries[i_summary]);
    if (it != GetSummariesMap.end()) {
      funcpt f = GetSummariesMap[arrAccSummaries[i_summary]];
      double s = (this->*f)(simulation);

      if (s < lparam->TOOSMALL || isinf(s)) {
        s = lparam->TOOSMALL;
      }

      s = transformSummaryStatistics(s);

      lparam->accessorysummariesSimulatedData.push_back(s);
    }
  }

  delete simulation;
  delete[] arrAccSummaries;
  delete[] arrSummaries;
  for (int taxa = 0; taxa < lparam->Ntaxa; taxa++) {
    delete[] CurrentLeafNodeCodonSequences_[taxa];
  }
  delete[] CurrentLeafNodeCodonSequences_;
}

void SummaryStatistics::computeSummaries() {
  lparam->summariesRealData.clear();
  lparam->summariesRealData.shrink_to_fit();

  lparam->accessorysummariesRealData.clear();
  lparam->accessorysummariesRealData.shrink_to_fit();

  RSCU_bool = false;
  relativeCodonFrequency_bool = false;
  dinuc_bool = false;
  dinuc12_bool = false;
  dinuc23_bool = false;
  dinuc31_bool = false;
  relativeAAFrequency_bool = false;
  dicodon_bool = false;
  diaa_bool = false;
  nuc_bool = false;
  nuc1_bool = false;
  nuc2_bool = false;
  nuc3_bool = false;
  nuc_meandiff_bool = false;
  nuc1_meandiff_bool = false;
  nuc2_meandiff_bool = false;
  nuc3_meandiff_bool = false;
  codon_meandiff_bool = false;
  aa_meandiff_bool = false;
  nuc_pairwise_bool = false;
  nuc_pairwise_bool10 = false;
  nuc_pairwise_bool30 = false;
  nuc_pairwise_bool50 = false;
  nuc_pairwise_bool70 = false;
  nuc_pairwise_bool90 = false;
  nuc_pairwiseSyn_bool = false;
  nuc_pairwiseSyn_bool10 = false;
  nuc_pairwiseSyn_bool30 = false;
  nuc_pairwiseSyn_bool50 = false;
  nuc_pairwiseSyn_bool70 = false;
  nuc_pairwiseSyn_bool90 = false;
  nuc1_pairwise_bool = false;
  nuc1_pairwise_bool10 = false;
  nuc1_pairwise_bool30 = false;
  nuc1_pairwise_bool50 = false;
  nuc1_pairwise_bool70 = false;
  nuc1_pairwise_bool90 = false;
  nuc2_pairwise_bool = false;
  nuc2_pairwise_bool10 = false;
  nuc2_pairwise_bool30 = false;
  nuc2_pairwise_bool50 = false;
  nuc2_pairwise_bool70 = false;
  nuc2_pairwise_bool90 = false;
  nuc3_pairwise_bool = false;
  nuc3_pairwise_bool10 = false;
  nuc3_pairwise_bool30 = false;
  nuc3_pairwise_bool50 = false;
  nuc3_pairwise_bool70 = false;
  nuc3_pairwise_bool90 = false;
  aa_pairwise_bool = false;
  aa_pairwise_bool10 = false;
  aa_pairwise_bool30 = false;
  aa_pairwise_bool50 = false;
  aa_pairwise_bool70 = false;
  aa_pairwise_bool90 = false;
  dinucCpG_pairwise_bool = false;
  nuc_site_comphet_bool = false;
  nuc1_site_comphet_bool = false;
  nuc2_site_comphet_bool = false;
  nuc3_site_comphet_bool = false;
  nuc_taxa_comphet_bool = false;
  nuc1_taxa_comphet_bool = false;
  nuc2_taxa_comphet_bool = false;
  nuc3_taxa_comphet_bool = false;
  codon_site_comphet_bool = false;
  codon_taxa_comphet_bool = false;
  aa_site_comphet_bool = false;
  aa_taxa_comphet_bool = false;
  RSCUentropy_bool = false;
  GC_bool = false;
  GC1_bool = false;
  GC2_bool = false;
  GC3_bool = false;

  string* arrSummaries = new string[lparam->NusedSummaries];
  for (int i_summary = 0; i_summary < lparam->NSummaries; i_summary++) {
    auto it = lparam->mapUsedSummaries.find(lparam->listSummaries[i_summary]);
    if (it != lparam->mapUsedSummaries.end() && it->second != -1) {
      arrSummaries[it->second] = it->first;
    }
  }

  for (int i_summary = 0; i_summary < lparam->NusedSummaries; i_summary++) {
    auto it = GetSummariesMap.find(arrSummaries[i_summary]);
    if (it != GetSummariesMap.end()) {
      funcpt f = GetSummariesMap[arrSummaries[i_summary]];

      double s = (this->*f)(lparam->codondata);

      if (s < lparam->TOOSMALL || isinf(s)) {
        s = lparam->TOOSMALL;
      }

      s = transformSummaryStatistics(s);

      lparam->summariesRealData.push_back(s);
    }
  }

  string* arrAccSummaries = new string[lparam->NusedAccessorySummaries];
  for (int i_summary = 0; i_summary < lparam->NSummaries; i_summary++) {
    auto it = lparam->mapUsedAccessorySummaries.find(
        lparam->listSummaries[i_summary]);
    if (it != lparam->mapUsedAccessorySummaries.end() && it->second != -1) {
      arrAccSummaries[it->second] = it->first;
    }
  }

  for (int i_summary = 0; i_summary < lparam->NusedAccessorySummaries;
       i_summary++) {
    auto it = GetSummariesMap.find(arrAccSummaries[i_summary]);
    if (it != GetSummariesMap.end()) {
      funcpt f = GetSummariesMap[arrAccSummaries[i_summary]];

      double s = (this->*f)(lparam->codondata);

      if (s < lparam->TOOSMALL || isinf(s)) {
        s = lparam->TOOSMALL;
      }

      s = transformSummaryStatistics(s);

      lparam->accessorysummariesRealData.push_back(s);
    }
  }
  delete[] arrAccSummaries;
  delete[] arrSummaries;
}

void SummaryStatistics::computeSummariesFromData() {
  ldata->summariesRealData.clear();
  ldata->summariesRealData.shrink_to_fit();

  RSCU_bool = false;
  relativeCodonFrequency_bool = false;
  dinuc_bool = false;
  dinuc12_bool = false;
  dinuc23_bool = false;
  dinuc31_bool = false;
  relativeAAFrequency_bool = false;
  dicodon_bool = false;
  diaa_bool = false;
  nuc_bool = false;
  nuc1_bool = false;
  nuc2_bool = false;
  nuc3_bool = false;
  nuc_meandiff_bool = false;
  nuc1_meandiff_bool = false;
  nuc2_meandiff_bool = false;
  nuc3_meandiff_bool = false;
  codon_meandiff_bool = false;
  aa_meandiff_bool = false;
  nuc_pairwise_bool = false;
  nuc_pairwise_bool10 = false;
  nuc_pairwise_bool30 = false;
  nuc_pairwise_bool50 = false;
  nuc_pairwise_bool70 = false;
  nuc_pairwise_bool90 = false;
  nuc_pairwiseSyn_bool = false;
  nuc_pairwiseSyn_bool10 = false;
  nuc_pairwiseSyn_bool30 = false;
  nuc_pairwiseSyn_bool50 = false;
  nuc_pairwiseSyn_bool70 = false;
  nuc_pairwiseSyn_bool90 = false;
  nuc1_pairwise_bool = false;
  nuc1_pairwise_bool10 = false;
  nuc1_pairwise_bool30 = false;
  nuc1_pairwise_bool50 = false;
  nuc1_pairwise_bool70 = false;
  nuc1_pairwise_bool90 = false;
  nuc2_pairwise_bool = false;
  nuc2_pairwise_bool10 = false;
  nuc2_pairwise_bool30 = false;
  nuc2_pairwise_bool50 = false;
  nuc2_pairwise_bool70 = false;
  nuc2_pairwise_bool90 = false;
  nuc3_pairwise_bool = false;
  nuc3_pairwise_bool10 = false;
  nuc3_pairwise_bool30 = false;
  nuc3_pairwise_bool50 = false;
  nuc3_pairwise_bool70 = false;
  nuc3_pairwise_bool90 = false;
  aa_pairwise_bool = false;
  aa_pairwise_bool10 = false;
  aa_pairwise_bool30 = false;
  aa_pairwise_bool50 = false;
  aa_pairwise_bool70 = false;
  aa_pairwise_bool90 = false;
  dinucCpG_pairwise_bool = false;
  nuc_site_comphet_bool = false;
  nuc1_site_comphet_bool = false;
  nuc2_site_comphet_bool = false;
  nuc3_site_comphet_bool = false;
  nuc_taxa_comphet_bool = false;
  nuc1_taxa_comphet_bool = false;
  nuc2_taxa_comphet_bool = false;
  nuc3_taxa_comphet_bool = false;
  codon_site_comphet_bool = false;
  codon_taxa_comphet_bool = false;
  aa_site_comphet_bool = false;
  aa_taxa_comphet_bool = false;
  RSCUentropy_bool = false;
  GC_bool = false;
  GC1_bool = false;
  GC2_bool = false;
  GC3_bool = false;

  string* arrSummaries = new string[ldata->NusedSummaries];
  for (int i_summary = 0; i_summary < ldata->NSummaries; i_summary++) {
    auto it = ldata->mapUsedSummaries.find(ldata->listSummaries[i_summary]);
    if (it != ldata->mapUsedSummaries.end() && it->second != -1) {
      arrSummaries[it->second] = it->first;
    }
  }

  for (int i_summary = 0; i_summary < ldata->NusedSummaries; i_summary++) {
    auto it = GetSummariesMap.find(arrSummaries[i_summary]);
    if (it != GetSummariesMap.end()) {
      funcpt f = GetSummariesMap[arrSummaries[i_summary]];

      double s = (this->*f)(ldata->codondata);
      // double s = log2(it.*second(simulation));
      if (s < ldata->TOOSMALL || isinf(s)) {
        s = ldata->TOOSMALL;
      }

      s = transformSummaryStatistics(s);

      ldata->summariesRealData.push_back(s);
    }
  }

  delete[] arrSummaries;
}

double SummaryStatistics::transformSummaryStatistics(double s) {
  if (lparam->transformation == "log2") {
    s = log2(s);
  } else if (lparam->transformation == "log10") {
    s = log10(s);
  } else if (lparam->transformation == "ln") {
    s = log(s);
  }
  return s;
}
