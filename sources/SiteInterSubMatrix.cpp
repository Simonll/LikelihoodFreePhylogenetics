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

#include "SiteInterSubMatrix.h"

#include <string>

SiteInterSubMatrix::SiteInterSubMatrix(LocalParameters* lparam) {
  this->lparam = lparam;

  // init SiteInterSubMatrix containers

  submatrixTreeSim = new double**[lparam->refTree->GetNnode()];
  mutmatrixTreeSim = new double**[lparam->refTree->GetNnode()];
  // selmatrixTreeSim = new double**[lparam->refTree->GetNnode()];
  for (int node = 0; node < lparam->refTree->GetNnode(); node++) {
    submatrixTreeSim[node] = new double*[lparam->Nsite_nuc];
    mutmatrixTreeSim[node] = new double*[lparam->Nsite_nuc];
    // selmatrixTreeSim[node] = new double*[lparam->Nsite_nuc];
    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++) {
      submatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
      mutmatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
      // selmatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
    }
  }

  TotalSubRate = new double[lparam->refTree->GetNnode()];
  TotalMutRate = new double[lparam->refTree->GetNnode()];
  /* TotalSubRateNonSyn = new double[lparam->refTree->GetNnode()];
  TotalMutRateNonSyn = new double[lparam->refTree->GetNnode()];
  TotalSubRateSyn = new double[lparam->refTree->GetNnode()];
  TotalMutRateSyn = new double[lparam->refTree->GetNnode()]; */

  PartialSubRate = new double[lparam->refTree->GetNnode()];
  PartialMutRate = new double[lparam->refTree->GetNnode()];
  /* PartialSubRateNonSyn = new double[lparam->refTree->GetNnode()];
  PartialMutRateNonSyn = new double[lparam->refTree->GetNnode()];
  PartialSubRateSyn = new double[lparam->refTree->GetNnode()];
  PartialMutRateSyn = new double[lparam->refTree->GetNnode()]; */
}

SiteInterSubMatrix::SiteInterSubMatrix(LocalParameters* lparam,
                                       std::string seq) {
  this->lparam = lparam;

  // init SiteInterSubMatrix containers

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
  /* TotalSubRateNonSyn = new double[lparam->refTree->GetNnode()];
  TotalMutRateNonSyn = new double[lparam->refTree->GetNnode()];
  TotalSubRateSyn = new double[lparam->refTree->GetNnode()];
  TotalMutRateSyn = new double[lparam->refTree->GetNnode()]; */

  PartialSubRate = new double[lparam->Ntaxa];
  PartialMutRate = new double[lparam->Ntaxa];
  /* PartialSubRateNonSyn = new double[lparam->refTree->GetNnode()];
  PartialMutRateNonSyn = new double[lparam->refTree->GetNnode()];
  PartialSubRateSyn = new double[lparam->refTree->GetNnode()];
  PartialMutRateSyn = new double[lparam->refTree->GetNnode()]; */
}

SiteInterSubMatrix::~SiteInterSubMatrix() {
  // dtor
}

void SiteInterSubMatrix::resetSubMatrix() {
  for (int node = 0; node < lparam->refTree->GetNnode(); node++) {
    TotalMutRate[node] = 0.0;
    TotalSubRate[node] = 0.0;
    /* TotalMutRateNonSyn[node] = 0.0;
    TotalSubRateNonSyn[node] = 0.0;
    TotalMutRateSyn[node] = 0.0;
    TotalSubRateSyn[node] = 0.0; */

    PartialMutRate[node] = 0.0;
    PartialSubRate[node] = 0.0;
    /*  PartialMutRateNonSyn[node] = 0.0;
     PartialSubRateNonSyn[node] = 0.0;
     PartialMutRateSyn[node] = 0.0;
     PartialSubRateSyn[node] = 0.0; */
    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        submatrixTreeSim[node][site_nuc][nuc] = 0.0;
        mutmatrixTreeSim[node][site_nuc][nuc] = 0.0;
        // selmatrixTreeSim[node][site_nuc][nuc] = 0.0;
      }
    }
  }
}

void SiteInterSubMatrix::resetSubMatrixSeq() {
  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    TotalMutRate[taxa_i] = 0.0;
    TotalSubRate[taxa_i] = 0.0;
    /* TotalMutRateNonSyn[node] = 0.0;
    TotalSubRateNonSyn[node] = 0.0;
    TotalMutRateSyn[node] = 0.0;
    TotalSubRateSyn[node] = 0.0; */

    PartialMutRate[taxa_i] = 0.0;
    PartialSubRate[taxa_i] = 0.0;
    /*  PartialMutRateNonSyn[node] = 0.0;
     PartialSubRateNonSyn[node] = 0.0;
     PartialMutRateSyn[node] = 0.0;
     PartialSubRateSyn[node] = 0.0; */
    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        submatrixTreeSim[taxa_i][site_nuc][nuc] = 0.0;
        mutmatrixTreeSim[taxa_i][site_nuc][nuc] = 0.0;
        selmatrixTreeSim[taxa_i][site_nuc][nuc] = 0.0;
      }
    }
  }
}

int SiteInterSubMatrix::testGpTcontext(int inNnodeIndex, int insite,
                                       int innucFrom, int innucTo,
                                       int** CurrentNodeNucSequence) {
  // if insite is one before end seq and  cur_insite = G and cur_insite+1 = T
  // and G goes to T then GpTHyp
  if (insite < lparam->Nsite_codon * 3 - 1 &&
      (CurrentNodeNucSequence[inNnodeIndex][insite] == 2 &&
       CurrentNodeNucSequence[inNnodeIndex][insite + 1] == 3) &&
      innucFrom == 2 && innucTo == 3) {
    return 1;  // TpT
  } else {
    return -1;
  }
}

int SiteInterSubMatrix::testGCPref(int innucFrom, int innucTo) {
  // from A|T
  if (innucFrom == 0 || innucFrom == 3) {
    // To C|G
    if (innucTo == 1 || innucTo == 2) {
      return -1;
    } else {  // To A|T
      return 0;
    }

  } else {
    // From C|G to A|T
    if (innucTo == 0 || innucTo == 3) {
      return 1;
    } else {  // To C|G
      return 0;
    }
  }
}

int SiteInterSubMatrix::testCpGcontext(int inNnodeIndex, int insite,
                                       int innucFrom, int innucTo,
                                       int** CurrentNodeNucSequence) {
  int X_ = 1;
  int Y_ = 2;

  // from XpY to NpY
  if (insite < lparam->Nsite_codon * 3 - 1 &&
      (CurrentNodeNucSequence[inNnodeIndex][insite] == X_ &&
       CurrentNodeNucSequence[inNnodeIndex][insite + 1] == Y_)) {
    if (((innucFrom + innucTo) % 2) == 0) {
      // we are leaving XY, coordinate 0,1
      // trought ts
      return 1;
    } else {
      // we are leaving XY, coordinate 0,1
      // trought tv
      return 9;
    }
  } else if (insite > 0 &&
             (CurrentNodeNucSequence[inNnodeIndex][insite - 1] == X_ &&
              CurrentNodeNucSequence[inNnodeIndex][insite] == Y_)) {
    if (((innucFrom + innucTo) % 2) == 0) {  // from XpY to XpN
      // we are leaving XY, coordinate -1,0
      // trought ts
      return 2;
    } else {
      // we are leaving XY, coordinate -1,0
      // trought tv
      return 9;
    }
  } else if (insite < lparam->Nsite_codon * 3 - 1 &&
             (CurrentNodeNucSequence[inNnodeIndex][insite] != X_ &&
              CurrentNodeNucSequence[inNnodeIndex][insite + 1] == Y_) &&
             (innucFrom != X_ && innucTo == X_)) {
    if (((innucFrom + innucTo) % 2) == 0) {  // from {N}pY to XpY

      // we are landing on XY, coordinate 0,1
      // trought  ts
      return -1;
    } else {
      // we are landing on XY, coordinate 0,1
      // trought  ts
      return -2;
    }
  } else if (insite > 0 &&
             (CurrentNodeNucSequence[inNnodeIndex][insite - 1] == X_ &&
              CurrentNodeNucSequence[inNnodeIndex][insite] != Y_) &&
             (innucFrom != Y_ && innucTo == Y_)) {  // from Xp{N} to XpY
    if (((innucFrom + innucTo) % 2) == 0) {
      // we are landing on XY, coordinate -1,0
      // trought ts
      return -3;
    } else {
      // we are landing on XY, coordinate -1,0
      // trought  tv
      return -4;
    }
  } else {  // not XpY context
    return 0;
  }
}

int SiteInterSubMatrix::testTpAcontext(int inNnodeIndex, int insite,
                                       int innucFrom, int innucTo,
                                       int** CurrentNodeNucSequence) {
  int X_ = 3;
  int Y_ = 0;

  // from XpY to NpY
  if (insite < lparam->Nsite_codon * 3 - 1 &&
      (CurrentNodeNucSequence[inNnodeIndex][insite] == X_ &&
       CurrentNodeNucSequence[inNnodeIndex][insite + 1] == Y_)) {
    if (((innucFrom + innucTo) % 2) == 0) {
      // we are leaving XY, coordinate 0,1
      // trought ts
      return 1;
    } else {
      // we are leaving XY, coordinate 0,1
      // trought tv
      return 9;
    }
  } else if (insite > 0 &&
             (CurrentNodeNucSequence[inNnodeIndex][insite - 1] == X_ &&
              CurrentNodeNucSequence[inNnodeIndex][insite] ==
                  Y_)) {  // from XpY to XpN
    if (((innucFrom + innucTo) % 2) == 0) {
      // we are leaving XY, coordinate -1,0
      // trought ts
      return 2;
    } else {
      // we are leaving XY, coordinate -1,0
      // trought tv
      return 9;
    }
  } else if (insite < lparam->Nsite_codon * 3 - 1 &&
             (CurrentNodeNucSequence[inNnodeIndex][insite] != X_ &&
              CurrentNodeNucSequence[inNnodeIndex][insite + 1] == Y_) &&
             (innucFrom != X_ && innucTo == X_)) {  // from {N}pY to XpY
    if (((innucFrom + innucTo) % 2) == 0) {
      // we are landing on XY, coordinate 0,1
      // trought  ts
      return -1;
    } else {
      // we are landing on XY, coordinate 0,1
      // trought  ts
      return -2;
    }
  } else if (insite > 0 &&
             (CurrentNodeNucSequence[inNnodeIndex][insite - 1] == X_ &&
              CurrentNodeNucSequence[inNnodeIndex][insite] != Y_) &&
             (innucFrom != Y_ && innucTo == Y_)) {  // from Xp{N} to XpY
    if (((innucFrom + innucTo) % 2) == 0) {
      // we are landing on XY, coordinate -1,0
      // trought ts
      return -3;
    } else {
      // we are landing on XY, coordinate -1,0
      // trought  tv
      return -4;
    }
  } else {  // not XpY context
    return 0;
  }
}

void SiteInterSubMatrix::UpdateSubMatrixTreeSim(int NodeIndex, int site_codon,
                                                int** CurrentNodeNucSequence) {
  // int verbose = lparam->verbose;
  bool whole = false;
  if (site_codon == -1) {
    whole = true;
  }
  double deltaTotalSubRate = 0;
  double deltaTotalMutRate = 0;
  //    double deltaTotalSubRateNonSyn = 0;
  //    double deltaTotalMutRateNonSyn = 0;
  //    double deltaTotalSubRateSyn = 0;
  //    double deltaTotalMutRateSyn = 0;

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
    //        deltaTotalSubRateNonSyn -= GetPartialSubRateNonSyn(NodeIndex);
    //        deltaTotalMutRateNonSyn -= GetPartialMutRateNonSyn(NodeIndex);
    //        deltaTotalSubRateSyn -= GetPartialSubRateSyn(NodeIndex);
    //        deltaTotalMutRateSyn -= GetPartialMutRateSyn(NodeIndex);

    //        if(verbose)
    //        {
    //            std::cerr << "SiteInterSubMatrix::UpdateSubMatrixTreeSim " <<
    //            deltaTotalSubRateNonSyn << " " << deltaTotalSubRateSyn <<
    //            "\n";
    //        }
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
        //                double  MutRateNonSyn = 0.0;
        //                double  SubRateNonSyn = 0.0;
        //                double  MutRateSyn = 0.0;
        //                double  SubRateSyn = 0.0;

        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
              nucposFrom[0], nucposFrom[1], nucposFrom[2]);
          int codonTo = lparam->codonstatespace->GetCodonFromDNA(
              nucposTo[0], nucposTo[1], nucposTo[2]);
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            double m;

            auto it = lparam->gtrMap.find(NodeIndex);
            if (it == lparam->gtrMap.end()) {
              std::cerr << "Error when looking to the gtr along the tree\n";
              exit(0);
            }
            if (it->second == 1) {
              // std::cerr << "GTR " << 1 << "\n";
              m = lparam->gtnr1[nucposFrom[codonPos]][nucposTo[codonPos]];
            } else if (it->second == 2) {
              // std::cerr << "GTR " << 2 << "\n";
              m = lparam->gtnr2[nucposFrom[codonPos]][nucposTo[codonPos]];
            } else {
              m = lparam->gtnr[nucposFrom[codonPos]][nucposTo[codonPos]];
            }

            int CpGcont =
                testCpGcontext(NodeIndex, site_nuc, nucposFrom[codonPos],
                               nucposTo[codonPos], CurrentNodeNucSequence);
            int TpAcont =
                testTpAcontext(NodeIndex, site_nuc, nucposFrom[codonPos],
                               nucposTo[codonPos], CurrentNodeNucSequence);
            int GCPref = testGCPref(nucposFrom[codonPos], nucposTo[codonPos]);
            if (CpGcont == 1 || CpGcont == 2) {
              // tsCpG
              m *= lparam->lambda_CpG;
              m *= lparam->lambda_tstvCpG;
              // CpG>TpG
            } else if (CpGcont > 2) {
              // tvCpG
              m *= lparam->lambda_tvCpG;
              m *= lparam->lambda_tstvCpG;
            }

            if (TpAcont == 1 || TpAcont == 2) {
              // tsTpA
              m *= lparam->lambda_TpA;
              m *= lparam->lambda_tstvTpA;
            } else if (TpAcont > 2) {
              // tvTpA
              m *= lparam->lambda_tvTpA;
              m *= lparam->lambda_tstvTpA;
            }

            if (m < lparam->TOOSMALL) {
              MutRate = lparam->TOOSMALL;
            } else {
              MutRate = m;
            }

            MutRate *= lparam->lambda_TBL;

            if (!lparam->codonstatespace->Synonymous(codonFrom, codonTo)) {
              int aaTo = lparam->codonstatespace->Translation(codonTo);
              int aaFrom = lparam->codonstatespace->Translation(codonFrom);

              S = log(
                  std::pow(
                      lparam->ssaaprofiles[lparam->alloc[site_codon_i]][aaTo],
                      lparam->AAadj[aaTo]) /
                  std::pow(
                      lparam->ssaaprofiles[lparam->alloc[site_codon_i]][aaFrom],
                      lparam->AAadj[aaFrom]) *
                  lparam->codonprofile[codonTo] /
                  lparam->codonprofile[codonFrom]);

              /* if(aaTo != 14 && aaFrom == 14)
              {
                  S += log((1.0-lparam->lambda_R) / lparam->lambda_R);
              }
              else if (aaTo == 14 && aaFrom != 14)
              {
                  S += log(lparam->lambda_R / (1.0 - lparam->lambda_R));
              } */

              if (aaTo != 14 && aaFrom == 14) {
                S += log((1.0 - lparam->lambda_R) / lparam->lambda_R);
              } else if (aaTo == 14 && aaFrom != 14) {
                S += log(lparam->lambda_R / (1.0 - lparam->lambda_R));
              }

              SubRate = MutRate * lparam->lambda_omega * lparam->omega *
                        lparam->omega_site[site_codon_i];
              // MutRateNonSyn = MutRate;
              // SubRateNonSyn = SubRate;
            } else {
              S = log(lparam->codonprofile[codonTo] /
                      lparam->codonprofile[codonFrom]);
              SubRate = MutRate;
              // MutRateSyn = MutRate;
              // SubRateSyn = MutRate;
            }

            if (GCPref > 0) {
              S += log((1.0 - lparam->fitGC) / lparam->fitGC);
            } else if (GCPref < 0) {
              S += log(lparam->fitGC / (1.0 - lparam->fitGC));
            }
            // leave from CpG
            if (CpGcont > 0) {
              S += log((1.0 - lparam->fitCpG) / lparam->fitCpG);
            } else if (CpGcont < 0) {  // land on CpG
              S += log(lparam->fitCpG / (1.0 - lparam->fitCpG));
            }

            if (TpAcont > 0) {
              S += log((1.0 - lparam->fitTpA) / lparam->fitTpA);
            } else if (TpAcont < 0) {
              S += log(lparam->fitTpA / (1.0 - lparam->fitTpA));
            }

            // NUMERICAL SECURITY
            if (fabs(S) < lparam->TOOSMALL) {
              SubRate /= (1.0 - (S / 2));
            } else if (S > lparam->TOOLARGE) {
              SubRate *= S;

            } else if (S < lparam->TOOLARGENEGATIVE) {
              SubRate = 0.0;
            } else {
              SubRate *= (S / (1.0 - exp(-S)));
            }

            //                        if
            //                        (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
            //                        {
            //
            //                            SubRateNonSyn = SubRate;
            //                            MutRateNonSyn = MutRate;
            //                        }
            //                        else
            //                        {
            //
            //                            SubRateSyn = SubRate;
            //                            MutRateSyn = MutRate;
            //                        }

            // NUMERICAL SECURITY linked to substitution process
            if (SubRate < 0) {
              std::cerr << "negative entry in matrix\n";
              std::cerr << "S: " << S << "\n";
              exit(1);
            }
            if (isinf(SubRate)) {
              std::cerr << "isinf\n";
              std::cerr << site_nuc << "\t" << nucTo << "\t" << MutRate << "\t"
                        << S << "\n";
              exit(1);
            }
            if (isnan(SubRate)) {
              std::cerr << "isnan\n";
              std::cerr << site_nuc << "\t" << nucTo << "\t" << MutRate << "\t"
                        << S << "\n";
              exit(1);
            }
            if (lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                   nucposTo[2])) {
              std::cerr << nucposTo[0] << " " << nucposTo[1] << " "
                        << nucposTo[2] << " " << MutRate << " " << SubRate
                        << " " << submatrixTreeSim[0][site_nuc][nucTo] << "\n";
              exit(0);
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
        ////
        // IF nucposFrom[codonPos] == nucTo set prob to 0
        ////
        deltaTotalSubRate += SubRate;
        deltaTotalMutRate += MutRate;
        //                deltaTotalSubRateNonSyn += SubRateNonSyn;
        //                deltaTotalMutRateNonSyn += MutRateNonSyn;
        //                deltaTotalSubRateSyn += SubRateSyn;
        //                deltaTotalMutRateSyn += MutRateSyn;

        mutmatrixTreeSim[NodeIndex][site_nuc][nucTo] = MutRate;
        // selmatrixTreeSim[NodeIndex][site_nuc][nucTo] = S;
        submatrixTreeSim[NodeIndex][site_nuc][nucTo] = SubRate;

        //              std::cerr << "MutRate" <<
        //              submatrix_mut[inNodeIndex][site_nuc][nucTo] << "\n";
        /* if (lparam->opt == 3)
        {
            submatrixTreeSim[NodeIndex][site_nuc][nucTo] = MutRate;
        } */
      }
    }
  }

  if (whole) {
    TotalSubRate[NodeIndex] = deltaTotalSubRate;
    TotalMutRate[NodeIndex] = deltaTotalMutRate;
    //        TotalSubRateNonSyn[NodeIndex] = deltaTotalSubRateNonSyn;
    //        TotalMutRateNonSyn[NodeIndex] = deltaTotalMutRateNonSyn;
    //        TotalSubRateSyn[NodeIndex] = deltaTotalSubRateSyn;
    //        TotalMutRateSyn[NodeIndex] = deltaTotalMutRateSyn;
  } else {
    TotalSubRate[NodeIndex] += deltaTotalSubRate;
    TotalMutRate[NodeIndex] += deltaTotalMutRate;
    //        TotalSubRateNonSyn[NodeIndex] += deltaTotalSubRateNonSyn;
    //        TotalMutRateNonSyn[NodeIndex] += deltaTotalMutRateNonSyn;
    //        TotalSubRateSyn[NodeIndex] += deltaTotalSubRateSyn;
    //        TotalMutRateSyn[NodeIndex] += deltaTotalMutRateSyn;
  }

  /* if (lparam->opt == 3)
  {
      TotalSubRate[NodeIndex] = TotalMutRate[NodeIndex];
  } */
  delete[] nucposFrom;
  delete[] nucposTo;
}

void SiteInterSubMatrix::UpdateSubMatrixSeq(int taxa,
                                            int** CurrentLeafNodeNucSequences) {
  // int verbose = lparam->verbose;

  double deltaTotalSubRate = 0;
  double deltaTotalMutRate = 0;
  //    double deltaTotalSubRateNonSyn = 0;
  //    double deltaTotalMutRateNonSyn = 0;
  //    double deltaTotalSubRateSyn = 0;
  //    double deltaTotalMutRateSyn = 0;

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
        //                double  MutRateNonSyn = 0.0;
        //                double  SubRateNonSyn = 0.0;
        //                double  MutRateSyn = 0.0;
        //                double  SubRateSyn = 0.0;

        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          int codonFrom = lparam->codonstatespace->GetCodonFromDNA(
              nucposFrom[0], nucposFrom[1], nucposFrom[2]);
          int codonTo = lparam->codonstatespace->GetCodonFromDNA(
              nucposTo[0], nucposTo[1], nucposTo[2]);
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            double m;

            m = lparam->gtnr[nucposFrom[codonPos]][nucposTo[codonPos]];

            int CpGcont =
                testCpGcontext(taxa, site_nuc, nucposFrom[codonPos],
                               nucposTo[codonPos], CurrentLeafNodeNucSequences);
            int TpAcont =
                testTpAcontext(taxa, site_nuc, nucposFrom[codonPos],
                               nucposTo[codonPos], CurrentLeafNodeNucSequences);
            int GCPref = testGCPref(nucposFrom[codonPos], nucposTo[codonPos]);
            if (CpGcont == 1 || CpGcont == 2) {
              // tsCpG
              m *= lparam->lambda_CpG;
              m *= lparam->lambda_tstvCpG;
              // CpG>TpG
            } else if (CpGcont > 2) {
              // tvCpG
              m *= lparam->lambda_tvCpG;
              m *= lparam->lambda_tstvCpG;
            }

            if (TpAcont == 1 || TpAcont == 2) {
              // tsTpA
              m *= lparam->lambda_TpA;
              m *= lparam->lambda_tstvTpA;
            } else if (TpAcont > 2) {
              // tvTpA
              m *= lparam->lambda_tvTpA;
              m *= lparam->lambda_tstvTpA;
            }

            if (m < lparam->TOOSMALL) {
              MutRate = lparam->TOOSMALL;
            } else {
              MutRate = m;
            }

            MutRate *= lparam->lambda_TBL;

            if (!lparam->codonstatespace->Synonymous(codonFrom, codonTo)) {
              int aaTo = lparam->codonstatespace->Translation(codonTo);
              int aaFrom = lparam->codonstatespace->Translation(codonFrom);

              S = log(
                  std::pow(
                      lparam->ssaaprofiles[lparam->alloc[site_codon_i]][aaTo],
                      lparam->AAadj[aaTo]) /
                  std::pow(
                      lparam->ssaaprofiles[lparam->alloc[site_codon_i]][aaFrom],
                      lparam->AAadj[aaFrom]) *
                  lparam->codonprofile[codonTo] /
                  lparam->codonprofile[codonFrom]);

              /* if(aaTo != 14 && aaFrom == 14)
              {
                  S += log((1.0-lparam->lambda_R) / lparam->lambda_R);
              }
              else if (aaTo == 14 && aaFrom != 14)
              {
                  S += log(lparam->lambda_R / (1.0 - lparam->lambda_R));
              } */

              if (aaTo != 14 && aaFrom == 14) {
                S += log((1.0 - lparam->lambda_R) / lparam->lambda_R);
              } else if (aaTo == 14 && aaFrom != 14) {
                S += log(lparam->lambda_R / (1.0 - lparam->lambda_R));
              }

              SubRate = MutRate * lparam->lambda_omega * lparam->omega *
                        lparam->omega_site[site_codon_i];
              // MutRateNonSyn = MutRate;
              // SubRateNonSyn = SubRate;
            } else {
              S = log(lparam->codonprofile[codonTo] /
                      lparam->codonprofile[codonFrom]);
              SubRate = MutRate;
              // MutRateSyn = MutRate;
              // SubRateSyn = MutRate;
            }

            if (GCPref > 0) {
              S += log((1.0 - lparam->fitGC) / lparam->fitGC);
            } else if (GCPref < 0) {
              S += log(lparam->fitGC / (1.0 - lparam->fitGC));
            }
            // leave from CpG
            if (CpGcont > 0) {
              S += log((1.0 - lparam->fitCpG) / lparam->fitCpG);
            } else if (CpGcont < 0) {  // land on CpG
              S += log(lparam->fitCpG / (1.0 - lparam->fitCpG));
            }

            if (TpAcont > 0) {
              S += log((1.0 - lparam->fitTpA) / lparam->fitTpA);
            } else if (TpAcont < 0) {
              S += log(lparam->fitTpA / (1.0 - lparam->fitTpA));
            }

            // NUMERICAL SECURITY
            if (fabs(S) < lparam->TOOSMALL) {
              SubRate /= (1.0 - (S / 2));
            } else if (S > lparam->TOOLARGE) {
              SubRate *= S;

            } else if (S < lparam->TOOLARGENEGATIVE) {
              SubRate = 0.0;
            } else {
              SubRate *= (S / (1.0 - exp(-S)));
            }

            //                        if
            //                        (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
            //                        {
            //
            //                            SubRateNonSyn = SubRate;
            //                            MutRateNonSyn = MutRate;
            //                        }
            //                        else
            //                        {
            //
            //                            SubRateSyn = SubRate;
            //                            MutRateSyn = MutRate;
            //                        }

            // NUMERICAL SECURITY linked to substitution process
            if (SubRate < 0) {
              std::cerr << "negative entry in matrix\n";
              std::cerr << "S: " << S << "\n";
              exit(1);
            }
            if (isinf(SubRate)) {
              std::cerr << "isinf\n";
              std::cerr << site_nuc << "\t" << nucTo << "\t" << MutRate << "\t"
                        << S << "\n";
              exit(1);
            }
            if (isnan(SubRate)) {
              std::cerr << "isnan\n";
              std::cerr << site_nuc << "\t" << nucTo << "\t" << MutRate << "\t"
                        << S << "\n";
              exit(1);
            }
            if (lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                   nucposTo[2])) {
              std::cerr << nucposTo[0] << " " << nucposTo[1] << " "
                        << nucposTo[2] << " " << MutRate << " " << SubRate
                        << " " << submatrixTreeSim[0][site_nuc][nucTo] << "\n";
              exit(0);
            }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
        ////
        // IF nucposFrom[codonPos] == nucTo set prob to 0
        ////
        deltaTotalSubRate += SubRate;
        deltaTotalMutRate += MutRate;
        //                deltaTotalSubRateNonSyn += SubRateNonSyn;
        //                deltaTotalMutRateNonSyn += MutRateNonSyn;
        //                deltaTotalSubRateSyn += SubRateSyn;
        //                deltaTotalMutRateSyn += MutRateSyn;

        mutmatrixTreeSim[taxa][site_nuc][nucTo] = MutRate;
        selmatrixTreeSim[taxa][site_nuc][nucTo] = S;
        submatrixTreeSim[taxa][site_nuc][nucTo] = SubRate;

        //              std::cerr << "MutRate" <<
        //              submatrix_mut[inNodeIndex][site_nuc][nucTo] << "\n";
        /* if (lparam->opt == 3)
        {
            submatrixTreeSim[NodeIndex][site_nuc][nucTo] = MutRate;
        } */
      }
    }
  }

  TotalSubRate[taxa] = deltaTotalSubRate;
  TotalMutRate[taxa] = deltaTotalMutRate;
  //        TotalSubRateNonSyn[NodeIndex] = deltaTotalSubRateNonSyn;
  //        TotalMutRateNonSyn[NodeIndex] = deltaTotalMutRateNonSyn;
  //        TotalSubRateSyn[NodeIndex] = deltaTotalSubRateSyn;
  //        TotalMutRateSyn[NodeIndex] = deltaTotalMutRateSyn;

  /* if (lparam->opt == 3)
  {
      TotalSubRate[NodeIndex] = TotalMutRate[NodeIndex];
  } */
  delete[] nucposFrom;
  delete[] nucposTo;
}

void SiteInterSubMatrix::transfertTotalRate(int sourceNodeIndex,
                                            int sinkNodeIndex) {
  TotalSubRate[sinkNodeIndex] = TotalSubRate[sourceNodeIndex];
  TotalMutRate[sinkNodeIndex] = TotalMutRate[sourceNodeIndex];
  //    TotalSubRateNonSyn[sinkNodeIndex] = TotalSubRateNonSyn[sourceNodeIndex]
  //    ; TotalMutRateNonSyn[sinkNodeIndex] =
  //    TotalMutRateNonSyn[sourceNodeIndex] ; TotalSubRateSyn[sinkNodeIndex] =
  //    TotalSubRateSyn[sourceNodeIndex] ; TotalMutRateSyn[sinkNodeIndex] =
  //    TotalMutRateSyn[sourceNodeIndex] ;
}

void SiteInterSubMatrix::transfertNodeMatrix(int sourceNodeIndex,
                                             int sinkNodeIndex, int site_nuc) {
  for (int i = 0; i < 4; i++) {
    submatrixTreeSim[sinkNodeIndex][site_nuc][i] =
        submatrixTreeSim[sourceNodeIndex][site_nuc][i];
    mutmatrixTreeSim[sinkNodeIndex][site_nuc][i] =
        mutmatrixTreeSim[sourceNodeIndex][site_nuc][i];
    // selmatrixTreeSim[sinkNodeIndex][site_nuc][i] =
    // selmatrixTreeSim[sourceNodeIndex][site_nuc][i];
  }
}

int SiteInterSubMatrix::testContextDinuc(int NodeIndex, int site_nuc,
                                         int* context, int nucTo,
                                         int** CurrentNodeNucSequence) {
  int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
  int MutPos = -1;

  // Find dinucleotide mutating position for the actual partern to be tested
  if (context[0] == context[2] && context[1] == nucFrom &&
      context[3] == nucTo) {
    MutPos = 1;
  } else if (context[0] == nucFrom && context[1] == context[3] &&
             context[2] == nucTo) {
    MutPos = 0;
  }
  // Test if the mutating patern is recovered from the actual sequence
  if (MutPos == 0 && site_nuc < lparam->Nsite_nuc - 1) {
    if (CurrentNodeNucSequence[NodeIndex][site_nuc] == context[MutPos] &&
        CurrentNodeNucSequence[NodeIndex][site_nuc + 1] == context[1]) {
      return 1;
    }
  } else if (MutPos == 1 && site_nuc > 0) {
    if (CurrentNodeNucSequence[NodeIndex][site_nuc] == context[MutPos] &&
        CurrentNodeNucSequence[NodeIndex][site_nuc - 1] == context[0]) {
      return 1;
    }
  }
  return -1;
}

//
// double SiteInterSubMatrix::GetMutRateSyn(int NodeIndex,int**
// CurrentNodeNucSequence)
//{
//    int nucFrom, pos1From, pos2From, pos3From, pos1To, pos2To, pos3To,
//    codonFrom, codonTo; double sum  = 0; for (int site_nuc = 0; site_nuc <
//    lparam->Nsite_nuc; site_nuc++)
//    {
//        nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
//        for (int nucTo = 0; nucTo < 4; nucTo++)
//        {
//            findCodonContext(NodeIndex, site_nuc, nucFrom, nucTo, pos1From,
//            pos2From, pos3From, pos1To, pos2To,
//            pos3To,CurrentNodeNucSequence);
//            if(!lparam->codonstatespace->CheckStop(pos1To, pos2To, pos3To))
//            {
//                codonFrom = lparam->codonstatespace->GetCodonFromDNA(pos1From,
//                pos2From, pos3From); codonTo =
//                lparam->codonstatespace->GetCodonFromDNA(pos1To, pos2To,
//                pos3To) ; if
//                (lparam->codonstatespace->Synonymous(codonFrom,codonTo))
//                {
//                    sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
//                }
//            }
//        }
//    }
//    return sum;
//}
//
// double SiteInterSubMatrix::GetSubRateCpG(int NodeIndex,int**
// CurrentNodeNucSequence)
//{
//    double sum  = 0.0;
//    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++)
//    {
//        int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
//        for (int nucTo = 0; nucTo < 4; nucTo++)
//        {
//            if(testCpGcontext(NodeIndex, site_nuc, nucFrom,
//            nucTo,CurrentNodeNucSequence) > 0)
//            {
//                sum += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
//            }
//        }
//    }
//    return sum;
//}
//
// double SiteInterSubMatrix::GetMutRateCpG(int NodeIndex,int**
// CurrentNodeNucSequence)
//{
//    double sum  = 0.0;
//    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++)
//    {
//        int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
//        for (int nucTo = 0; nucTo < 4; nucTo++)
//        {
//            if(testCpGcontext(NodeIndex, site_nuc, nucFrom,
//            nucTo,CurrentNodeNucSequence)> 0)
//            {
//                sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
//            }
//        }
//    }
//    return sum;
//}

double SiteInterSubMatrix::GetMutRate(int NodeIndex, int site_codon) {
  double sum = 0.0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  if (site_codon > -1) {
    if (site_codon < lparam->Nsite_codon - 1) {
      site_codon_end = site_codon + 2;
    }
    if (site_codon > 0) {
      site_codon_start = site_codon - 1;
    }
  }
  int site_nuc_start = site_codon_start * 3;
  int site_nuc_end = site_codon_end * 3;
  for (int site_nuc = site_nuc_start; site_nuc < site_nuc_end; site_nuc++) {
    for (int nucTo = 0; nucTo < 4; nucTo++) {
      sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
    }
  }
  // TotalMut[NodeIndex] = sum;
  return sum;
}

double SiteInterSubMatrix::GetSubRate(int NodeIndex, int site_codon) {
  double sum = 0;
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  if (site_codon > -1) {
    if (site_codon < lparam->Nsite_codon - 1) {
      site_codon_end = site_codon + 2;
    }
    if (site_codon > 0) {
      site_codon_start = site_codon - 1;
    }
  }
  int site_nuc_start = site_codon_start * 3;
  int site_nuc_end = site_codon_end * 3;
  for (int site_nuc = site_nuc_start; site_nuc < site_nuc_end; site_nuc++) {
    for (int nucTo = 0; nucTo < 4; nucTo++) {
      sum += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
    }
  }
  // TotalSub[NodeIndex] = sum;
  return sum;
}

void SiteInterSubMatrix::ComputePartialRates(int NodeIndex, int site_codon,
                                             int** CurrentNodeNucSequence) {
  int site_codon_start = 0;
  int site_codon_end = lparam->Nsite_codon;
  if (site_codon > -1) {
    if (site_codon < lparam->Nsite_codon - 2) {
      site_codon_end = site_codon + 2;
    }
    if (site_codon > 0) {
      site_codon_start = site_codon - 1;
    }
  }
  PartialSubRate[NodeIndex] = 0.0;
  PartialMutRate[NodeIndex] = 0.0;
  //    PartialSubRateNonSyn[NodeIndex] = 0.0;
  //    PartialMutRateNonSyn[NodeIndex] = 0.0;
  //    PartialSubRateSyn[NodeIndex] = 0.0;
  //    PartialMutRateSyn[NodeIndex] = 0.0;
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
        double MutRate = 0.0;
        double SubRate = 0.0;
        //                double  MutRateNonSyn = 0.0;
        //                double  SubRateNonSyn = 0.0;
        //                double  MutRateSyn = 0.0;
        //                double  SubRateSyn = 0.0;

        if (nucposFrom[codonPos] != nucTo) {
          nucposTo[codonPos] = nucTo;
          // int codonFrom =
          // lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0],
          // nucposFrom[1], nucposFrom[2]); int codonTo =
          // lparam->codonstatespace->GetCodonFromDNA(nucposTo[0], nucposTo[1],
          // nucposTo[2]);
          if (!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1],
                                                  nucposTo[2])) {
            SubRate = submatrixTreeSim[NodeIndex][site_nuc][nucTo];
            MutRate = mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];

            //                        if
            //                        (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
            //                        {
            //
            //                            SubRateNonSyn =
            //                            submatrixTreeSim[NodeIndex][site_nuc][nucTo];
            //                            MutRateNonSyn =
            //                            mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            //                        }
            //                        else
            //                        {
            //
            //                            SubRateSyn =
            //                            submatrixTreeSim[NodeIndex][site_nuc][nucTo];
            //                            MutRateSyn =
            //                            mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
            //                        }
          }
          nucposTo[codonPos] = nucposFrom[codonPos];
        }
        ////
        // IF nucposFrom[codonPos] == nucTo set prob to 0
        ////

        PartialSubRate[NodeIndex] += SubRate;
        PartialMutRate[NodeIndex] += MutRate;
        //                PartialSubRateNonSyn[NodeIndex] += SubRateNonSyn;
        //                PartialMutRateNonSyn[NodeIndex] += MutRateNonSyn;
        //                PartialSubRateSyn[NodeIndex] += SubRateSyn;
        //                PartialMutRateSyn[NodeIndex] += MutRateSyn;
      }
    }
  }
  delete[] nucposFrom;
  delete[] nucposTo;
}

void SiteInterSubMatrix::WriteSubMatrix(ostream& mutation_os,
                                        ostream& selection_os, int nucTo) {
  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    mutation_os << lparam->taxonset->GetTaxon(taxa_i) << "\t";
    selection_os << lparam->taxonset->GetTaxon(taxa_i) << "\t";

    for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
      mutation_os << mutmatrixTreeSim[taxa_i][site_codon][nucTo] << "\t";
      selection_os << selmatrixTreeSim[taxa_i][site_codon][nucTo] << "\t";
    }
    mutation_os << "\n";
    selection_os << "\n";
  }
}

//
// double SiteInterSubMatrix::GetSubRateNonSyn(int NodeIndex, int
// site_codon,int** CurrentNodeNucSequence)
//{
//
//    int verbose = lparam->verbose;
//    double sum  = 0.0;
//    int site_codon_start = 0 ;
//    int site_codon_end = lparam->Nsite_codon;
//
//    if(verbose)
//    {
//        std::cerr << "SiteInterSubMatrix::GetSubRateNonSyn " << NodeIndex << "
//        " << site_codon << "\n";
//    }
//
//    if (site_codon > -1)
//    {
//        if (site_codon < lparam->Nsite_codon-1)
//        {
//            site_codon_end = site_codon + 2;
//        }
//        if (site_codon > 0)
//        {
//            site_codon_start =  site_codon - 1;
//        }
//    }
//
//
//    int* nucposFrom = new int[3];
//    int* nucposTo = new int[3];
//
//    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ;
//    site_codon_i++)
//    {
//
//        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
//        for (int codonPos = 0; codonPos < 3; codonPos++)
//        {
//            nucposFrom[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//            nucposTo[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//        }
//
//        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each
//        nucleotide codon postions [0,2] we will be computing adjacent
//        nucleotide.
//        {
//            int site_nuc =  site_nuc_start+codonPos;
//            for (int nucTo = 0; nucTo < 4; nucTo++)
//            {
//                if (nucposFrom[codonPos] != nucTo)
//                {
//                    nucposTo[codonPos] = nucTo;
//
//                    if(lparam->codonstatespace->CheckStop(nucposTo[0],
//                    nucposTo[1], nucposTo[2]))
//                    {
//
//                        int codonFrom =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0],
//                        nucposFrom[1], nucposFrom[2]); int codonTo =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposTo[0],
//                        nucposTo[1], nucposTo[2]);
//
//                        if
//                        (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
//                        {
//
//                            sum +=
//                            submatrixTreeSim[NodeIndex][site_nuc][nucTo];
//
//                        }
//
//                    }
//                }
//
//            }
//
//        }
//    }
//
//
//
//    delete [] nucposFrom;
//    delete [] nucposTo;
//
//
//    return sum;
//}
//
// double SiteInterSubMatrix::GetSubRateSyn(int NodeIndex, int site_codon,int**
// CurrentNodeNucSequence)
//{
//    double sum  = 0.0;
//    int site_codon_start = 0 ;
//    int site_codon_end = lparam->Nsite_codon;
//
//
//    if (site_codon > -1)
//    {
//        if (site_codon < lparam->Nsite_codon-1)
//        {
//            site_codon_end = site_codon + 2;
//        }
//        if (site_codon > 0)
//        {
//            site_codon_start =  site_codon - 1;
//        }
//    }
//
//
//    int* nucposFrom = new int[3];
//    int* nucposTo = new int[3];
//
//    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ;
//    site_codon_i++)
//    {
//
//        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
//        for (int codonPos = 0; codonPos < 3; codonPos++)
//        {
//            nucposFrom[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//            nucposTo[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//        }
//
//        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each
//        nucleotide codon postions [0,2] we will be computing adjacent
//        nucleotide.
//        {
//            int site_nuc =  site_nuc_start+codonPos;
//            for (int nucTo = 0; nucTo < 4; nucTo++)
//            {
//                if (nucposFrom[codonPos] != nucTo)
//                {
//                    nucposTo[codonPos] = nucTo;
//
//                    if(lparam->codonstatespace->CheckStop(nucposTo[0],
//                    nucposTo[1], nucposTo[2]))
//                    {
//
//                        int codonFrom =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0],
//                        nucposFrom[1], nucposFrom[2]); int codonTo =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposTo[0],
//                        nucposTo[1], nucposTo[2]);
//
//                        if
//                        (lparam->codonstatespace->Synonymous(codonFrom,codonTo))
//                        {
//
//                            sum +=
//                            submatrixTreeSim[NodeIndex][site_nuc][nucTo];
//
//                        }
//
//
//
//                    }
//
//                }
//
//            }
//        }
//
//    }
//
//    delete [] nucposFrom;
//    delete [] nucposTo;
//
//
//    return sum;
//}
//
// double SiteInterSubMatrix::GetMutRateNonSyn(int NodeIndex, int
// site_codon,int** CurrentNodeNucSequence)
//{
//    double sum  = 0.0;
//    int site_codon_start = 0 ;
//    int site_codon_end = lparam->Nsite_codon;
//
//
//    if (site_codon > -1)
//    {
//        if (site_codon < lparam->Nsite_codon-1)
//        {
//            site_codon_end = site_codon + 2;
//        }
//        if (site_codon > 0)
//        {
//            site_codon_start =  site_codon - 1;
//        }
//    }
//
//
//    int* nucposFrom = new int[3];
//    int* nucposTo = new int[3];
//
//    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ;
//    site_codon_i++)
//    {
//
//        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
//        for (int codonPos = 0; codonPos < 3; codonPos++)
//        {
//            nucposFrom[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//            nucposTo[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//        }
//
//        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each
//        nucleotide codon postions [0,2] we will be computing adjacent
//        nucleotide.
//        {
//            int site_nuc =  site_nuc_start+codonPos;
//            for (int nucTo = 0; nucTo < 4; nucTo++)
//            {
//                if (nucposFrom[codonPos] != nucTo)
//                {
//                    nucposTo[codonPos] = nucTo;
//
//                    if(lparam->codonstatespace->CheckStop(nucposTo[0],
//                    nucposTo[1], nucposTo[2]))
//                    {
//
//                        int codonFrom =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0],
//                        nucposFrom[1], nucposFrom[2]); int codonTo =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposTo[0],
//                        nucposTo[1], nucposTo[2]);
//
//                        if
//                        (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
//                        {
//
//                            sum +=
//                            mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
//
//                        }
//
//                    }
//
//                }
//
//            }
//        }
//
//    }
//
//    delete [] nucposFrom;
//    delete [] nucposTo;
//
//
//    return sum;
//}
//
// double SiteInterSubMatrix::GetMutRateSyn(int NodeIndex, int site_codon,int**
// CurrentNodeNucSequence)
//{
//    double sum  = 0.0;
//    int site_codon_start = 0 ;
//    int site_codon_end = lparam->Nsite_codon;
//
//
//    if (site_codon > -1)
//    {
//        if (site_codon < lparam->Nsite_codon-1)
//        {
//            site_codon_end = site_codon + 2;
//        }
//        if (site_codon > 0)
//        {
//            site_codon_start =  site_codon - 1;
//        }
//    }
//
//
//
//    int* nucposFrom = new int[3];
//    int* nucposTo = new int[3];
//
//    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ;
//    site_codon_i++)
//    {
//
//        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
//        for (int codonPos = 0; codonPos < 3; codonPos++)
//        {
//            nucposFrom[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//            nucposTo[codonPos] =
//            CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
//        }
//
//        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each
//        nucleotide codon postions [0,2] we will be computing adjacent
//        nucleotide.
//        {
//            int site_nuc =  site_nuc_start+codonPos;
//            for (int nucTo = 0; nucTo < 4; nucTo++)
//            {
//                if (nucposFrom[codonPos] != nucTo)
//                {
//                    nucposTo[codonPos] = nucTo;
//
//                    if(lparam->codonstatespace->CheckStop(nucposTo[0],
//                    nucposTo[1], nucposTo[2]))
//                    {
//
//                        int codonFrom =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0],
//                        nucposFrom[1], nucposFrom[2]); int codonTo =
//                        lparam->codonstatespace->GetCodonFromDNA(nucposTo[0],
//                        nucposTo[1], nucposTo[2]);
//
//                        if
//                        (lparam->codonstatespace->Synonymous(codonFrom,codonTo))
//                        {
//
//                            sum +=
//                            mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
//
//                        }
//
//                    }
//
//                }
//
//            }
//        }
//
//    }
//
//    delete [] nucposFrom;
//    delete [] nucposTo;
//
//
//    return sum;
//}
