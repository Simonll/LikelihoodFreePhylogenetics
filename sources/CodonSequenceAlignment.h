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
/****

Adapted from: PhyloBayes MPI code

****/

/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel
Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. PhyloBayes is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a
copy of the GNU General Public License along with PhyloBayes. If not, see
<http://www.gnu.org/licenses/>.

**********************/

#ifndef CODONSEQUENCEALIGNMENT_H
#define CODONSEQUENCEALIGNMENT_H

//#include "ContinuousData.h"
#include <algorithm>
#include <cmath>
#include <vector>

#include "CodonStateSpace.h"
#include "SequenceAlignment.h"

class CodonSequenceAlignment : public SequenceAlignment {
 public:
  CodonSequenceAlignment(CodonSequenceAlignment* from)
      : SequenceAlignment((SequenceAlignment*)from) {}
  CodonSequenceAlignment(CodonSequenceAlignment* from, int** indata);
  CodonSequenceAlignment(CodonSequenceAlignment* from, int Ntaxa, int** indata);
  CodonSequenceAlignment(SequenceAlignment* from, bool force_stops = false,
                         GeneticCodeType type = Universal);

  ~CodonSequenceAlignment() {}

  int** GetData();

  void DeleteAAConstantSites() {
    int i = 0;
    int j = 0;
    int Eliminated = 0;
    while (i < Nsite) {
      int k = 0;
      while ((k < Ntaxa) && (Data[k][i] == unknown)) k++;
      if (k < Ntaxa) {
        int a = GetCodonStateSpace()->Translation(Data[k][i]);
        k++;
        while ((k < Ntaxa) &&
               ((Data[k][i] == unknown) ||
                (GetCodonStateSpace()->Translation(Data[k][i]) == a)))
          k++;
        if (k == Ntaxa) {
          Eliminated++;
        } else {
          for (int k = 0; k < Ntaxa; k++) {
            Data[k][j] = Data[k][i];
          }
          j++;
        }
      }
      i++;
    }
    Nsite -= Eliminated;
    cout << "number of positions eliminated : " << Eliminated << '\n';
  }

  virtual double GetTotalDiversity(int sitemin, int sitemax) {
    double total = 0;
    int obs[Naa];
    for (int i = sitemin; i < sitemax; i++) {
      for (int k = 0; k < Naa; k++) {
        obs[k] = 0;
      }
      for (int j = 0; j < Ntaxa; j++) {
        int state = GetState(j, i);
        if (state != unknown) {
          int aa = GetCodonStateSpace()->Translation(state);
          obs[aa]++;
        }
      }
      int div = 0;
      for (int k = 0; k < Naa; k++) {
        if (obs[k]) {
          div++;
        }
      }
      total += div;
    }
    return total;
  }

  double Nucleotide123CompositionalHeterogeneity() {
    double** taxfreq12 = new double*[Ntaxa];
    double** taxfreq3 = new double*[Ntaxa];
    for (int j = 0; j < Ntaxa; j++) {
      taxfreq12[j] = new double[Nnuc];
      taxfreq3[j] = new double[Nnuc];
      for (int k = 0; k < Nnuc; k++) {
        taxfreq12[j][k] = 0.0;
        taxfreq3[j][k] = 0.0;
      }
    }

    for (int i = 0; i < Nsite; i++) {
      for (int j = 0; j < Ntaxa; j++) {
        int state = GetState(j, i);
        if (state != unknown) {
          taxfreq12[j][GetCodonStateSpace()->GetCodonPosition(0, state)]++;
          taxfreq12[j][GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          taxfreq3[j][GetCodonStateSpace()->GetCodonPosition(2, state)]++;
        }
      }
    }

    for (int j = 0; j < Ntaxa; j++) {
      double total = 0;
      for (int k = 0; k < Nnuc; k++) {
        total += taxfreq12[j][k];
      }
      for (int k = 0; k < Nnuc; k++) {
        taxfreq12[j][k] /= total;
      }
    }

    for (int j = 0; j < Ntaxa; j++) {
      double total = 0;
      for (int k = 0; k < Nnuc; k++) {
        total += taxfreq3[j][k];
      }
      for (int k = 0; k < Nnuc; k++) {
        taxfreq3[j][k] /= total;
      }
    }

    // compute max distance
    double maxdist = 0;
    for (int j = 0; j < Ntaxa; j++) {
      double dist = 0;
      for (int k = 0; k < Nnuc; k++) {
        double tmp = (taxfreq12[j][k] - taxfreq3[j][k]);
        dist += tmp * tmp;
      }
      if (maxdist < dist) {
        maxdist = dist;
      }
    }

    for (int j = 0; j < Ntaxa; j++) {
      delete[] taxfreq12[j];
      delete[] taxfreq3[j];
    }
    delete[] taxfreq12;
    delete[] taxfreq3;

    return maxdist;
  }

  double NucleotideCompositionalHeterogeneity(ostream* os, int pos = -1,
                                              double** comp = 0) {
    double** taxfreq = 0;
    if (comp) {
      taxfreq = comp;
    } else {
      taxfreq = new double*[Ntaxa];
      for (int j = 0; j < Ntaxa; j++) {
        taxfreq[j] = new double[Nnuc];
        for (int k = 0; k < Nnuc; k++) {
          taxfreq[j][k] = 0;
        }
      }
    }

    for (int j = 0; j < Ntaxa; j++) {
      for (int k = 0; k < Nnuc; k++) {
        taxfreq[j][k] = 0;
      }
    }

    for (int i = 0; i < Nsite; i++) {
      for (int j = 0; j < Ntaxa; j++) {
        int state = GetState(j, i);
        if (state != unknown) {
          if (pos == -1) {
            taxfreq[j][GetCodonStateSpace()->GetCodonPosition(0, state)]++;
            taxfreq[j][GetCodonStateSpace()->GetCodonPosition(1, state)]++;
            taxfreq[j][GetCodonStateSpace()->GetCodonPosition(2, state)]++;
          } else if (pos > 2) {
            cerr << "error in "
                    "CodonSequenceAlignment::"
                    "NucleotideCompositionHeterogeneity : "
                 << pos << '\n';
            exit(1);
          } else {
            taxfreq[j][GetCodonStateSpace()->GetCodonPosition(pos, state)]++;
          }
        }
      }
    }

    // make global freqs out of tax-specific freqs
    double* globalfreq = new double[Nnuc];
    for (int k = 0; k < Nnuc; k++) {
      globalfreq[k] = 0;
      for (int j = 0; j < Ntaxa; j++) {
        globalfreq[k] += taxfreq[j][k];
      }
    }

    // normalise
    double total = 0;
    for (int k = 0; k < Nnuc; k++) {
      total += globalfreq[k];
    }
    for (int k = 0; k < Nnuc; k++) {
      globalfreq[k] /= total;
    }
    for (int j = 0; j < Ntaxa; j++) {
      double total = 0;
      for (int k = 0; k < Nnuc; k++) {
        total += taxfreq[j][k];
      }
      for (int k = 0; k < Nnuc; k++) {
        taxfreq[j][k] /= total;
        if (os) {
          (*os) << taxfreq[j][k] << '\t';
        }
      }
      if (os) {
        (*os) << '\n';
      }
    }
    if (os) {
      (*os) << '\n';
    }

    // compute max distance
    double maxdist = 0;
    for (int j = 0; j < Ntaxa; j++) {
      double dist = 0;
      for (int k = 0; k < Nnuc; k++) {
        double tmp = (taxfreq[j][k] - globalfreq[k]);
        dist += tmp * tmp;
      }
      if (maxdist < dist) {
        maxdist = dist;
      }
    }

    delete[] globalfreq;
    if (!comp) {
      for (int j = 0; j < Ntaxa; j++) {
        delete[] taxfreq[j];
      }
      delete[] taxfreq;
    }

    return maxdist;
  }

  double AminoAcidCompositionalHeterogeneity(ostream* os) {
    double** taxfreq = new double*[Ntaxa];
    for (int j = 0; j < Ntaxa; j++) {
      taxfreq[j] = new double[Naa];
      for (int k = 0; k < Naa; k++) {
        taxfreq[j][k] = 0;
      }
    }

    for (int i = 0; i < Nsite; i++) {
      for (int j = 0; j < Ntaxa; j++) {
        int state = GetState(j, i);
        if (state != unknown) {
          taxfreq[j][GetCodonStateSpace()->Translation(state)]++;
        }
      }
    }

    // make global freqs out of tax-specific freqs
    double* globalfreq = new double[Naa];
    for (int k = 0; k < Naa; k++) {
      globalfreq[k] = 0;
      for (int j = 0; j < Ntaxa; j++) {
        globalfreq[k] += taxfreq[j][k];
      }
    }

    // normalise
    double total = 0;
    for (int k = 0; k < Naa; k++) {
      total += globalfreq[k];
    }
    for (int k = 0; k < Naa; k++) {
      globalfreq[k] /= total;
    }
    for (int j = 0; j < Ntaxa; j++) {
      double total = 0;
      for (int k = 0; k < Naa; k++) {
        total += taxfreq[j][k];
      }
      for (int k = 0; k < Naa; k++) {
        taxfreq[j][k] /= total;
        if (os) {
          (*os) << taxfreq[j][k] << '\t';
        }
      }
      if (os) {
        (*os) << '\n';
      }
    }
    if (os) {
      (*os) << '\n';
    }

    // compute max distance
    double maxdist = 0;
    for (int j = 0; j < Ntaxa; j++) {
      double dist = 0;
      for (int k = 0; k < Naa; k++) {
        double tmp = (taxfreq[j][k] - globalfreq[k]);
        dist += tmp * tmp;
      }
      if (maxdist < dist) {
        maxdist = dist;
      }
    }

    delete[] globalfreq;
    for (int j = 0; j < Ntaxa; j++) {
      delete[] taxfreq[j];
    }
    delete[] taxfreq;

    return maxdist;
  }

  //-->SLL
  void nuc_meandiff(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int Nsite_nuc = Nsite * 3;
    int** in = new int*[Nsite_nuc];

    for (int site_nuc = 0; site_nuc < Nsite_nuc; site_nuc++) {
      in[site_nuc] = new int[Nnuc];
    }

    for (int site_nuc = 0; site_nuc < Nsite_nuc; site_nuc++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        in[site_nuc][nuc] = 0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon) ;
        if (state != unknown) {
          for (int pos = 0; pos < 3; pos++) {
            in[site_codon * 3 + pos]
              [GetCodonStateSpace()->GetCodonPosition(pos, state)] = 1;
          }
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_nuc = 0; site_nuc < Nsite_nuc; site_nuc++) {
      int Ndiff = 0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        if (in[site_nuc][nuc] != 0) {
          Ndiff++;
        }
      }
      mean += static_cast<double>(Ndiff);
      var += static_cast<double>(Ndiff * Ndiff);
    }

    mean /= static_cast<double>(Nsite_nuc);
    var /= static_cast<double>(Nsite_nuc);
    var = var - (mean * mean);

    stat_container[0] = mean;
    stat_container[1] = var;

    for (int site_nuc = 0; site_nuc < Nsite_nuc; site_nuc++) {
      delete[] in[site_nuc];
    }
    delete[] in;
  }

  void nuc1_meandiff(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int** in = new int*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = new int[Nnuc];
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        in[site_codon][nuc] = 0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          in[site_codon][GetCodonStateSpace()->GetCodonPosition(0, state)] = 1;
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      int Ndiff = 0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        if (in[site_codon][nuc] != 0) {
          Ndiff++;
        }
      }
      mean += static_cast<double>(Ndiff);
      var += static_cast<double>(Ndiff * Ndiff);
    }
    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);
    stat_container[0] = mean;
    stat_container[1] = var;

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] in[site_codon];
    }

    delete[] in;
  }

  void nuc2_meandiff(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int** in = new int*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = new int[Nnuc];
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        in[site_codon][nuc] = 0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          in[site_codon][GetCodonStateSpace()->GetCodonPosition(1, state)] = 1;
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      int Ndiff = 0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        if (in[site_codon][nuc] != 0) {
          Ndiff++;
        }
      }

      mean += static_cast<double>(Ndiff);
      var += static_cast<double>(Ndiff * Ndiff);
    }
    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);
    stat_container[0] = mean;
    stat_container[1] = var;

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] in[site_codon];
    }
    delete[] in;
  }

  void nuc3_meandiff(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int** in = new int*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = new int[Nnuc];
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        in[site_codon][nuc] = 0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          in[site_codon][GetCodonStateSpace()->GetCodonPosition(2, state)] = 1;
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      int Ndiff = 0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        if (in[site_codon][nuc] != 0) {
          Ndiff++;
        }
      }
      mean += static_cast<double>(Ndiff);
      var += static_cast<double>(Ndiff * Ndiff);
    }
    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);
    stat_container[0] = mean;
    stat_container[1] = var;

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] in[site_codon];
    }
    delete[] in;
  }

  void nuc_usage(double* stat_container) {
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] = 1.0;
    }
    int total = 4;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          for (int pos = 0; pos < 3; pos++) {
            stat_container[GetCodonStateSpace()->GetCodonPosition(pos,
                                                                  state)]++;
            total++;
          }
        }
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] /= total;
    }
  }

  void nuc1_usage(double* stat_container) {
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] = 1.0;
    }
    int total = 4;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState (taxa, site_codon);
        if (state != unknown) {
          stat_container[GetCodonStateSpace()->GetCodonPosition(0, state)]++;
          total++;
        }
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] /= total;
    }
  }

  void nuc2_usage(double* stat_container) {
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] = 1.0;
    }
    int total = 4;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState (taxa, site_codon);
        if (state != unknown) {
          stat_container[GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          total++;
        }
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] /= total;
    }
  }

  void nuc3_usage(double* stat_container) {
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] = 1.0;
    }
    int total = 4;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState (taxa, site_codon);
        if (state != unknown) {
          stat_container[GetCodonStateSpace()->GetCodonPosition(2, state)]++;
          total++;
        }
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      stat_container[nuc] /= total;
    }
  }

  void dinuc_usage(double** stat_container) {
    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < Nnuc; nuc2++) {
        stat_container[nuc1][nuc2] = 1.0;
      }
    }

    int tot = 16;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        int next_state = -1;
        if (site_codon < Nsite - 1) {
          next_state =
              Data[taxa][site_codon + 1];  // GetState(taxa,site_codon+1);
        }
        if (state != unknown && next_state != unknown) {
          stat_container[GetCodonStateSpace()->GetCodonPosition(0, state)]
                        [GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          stat_container[GetCodonStateSpace()->GetCodonPosition(1, state)]
                        [GetCodonStateSpace()->GetCodonPosition(2, state)]++;
          stat_container[GetCodonStateSpace()->GetCodonPosition(2, state)]
                        [GetCodonStateSpace()->GetCodonPosition(0,
                                                                next_state)]++;
          tot += 3;
        } else if (state != unknown && next_state == unknown) {
          stat_container[GetCodonStateSpace()->GetCodonPosition(0, state)]
                        [GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          stat_container[GetCodonStateSpace()->GetCodonPosition(1, state)]
                        [GetCodonStateSpace()->GetCodonPosition(2, state)]++;
          tot += 2;
        }
      }
    }

    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < Nnuc; nuc2++) {
        stat_container[nuc1][nuc2] /= tot;
      }
    }
  }

  void dinuc12_usage(double** stat_container) {
    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < 4; nuc2++) {
        stat_container[nuc1][nuc2] = 1.0;
      }
    }
    int tot = 16;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          stat_container[GetCodonStateSpace()->GetCodonPosition(0, state)]
                        [GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          tot++;
        }
      }
    }
    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < Nnuc; nuc2++) {
        stat_container[nuc1][nuc2] /= tot;
      }
    }
  }

  void dinuc23_usage(double** stat_container) {
    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < 4; nuc2++) {
        stat_container[nuc1][nuc2] = 1.0;
      }
    }
    int tot = 16;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          stat_container[GetCodonStateSpace()->GetCodonPosition(1, state)]
                        [GetCodonStateSpace()->GetCodonPosition(2, state)]++;
          tot++;
        }
      }
    }
    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < 4; nuc2++) {
        stat_container[nuc1][nuc2] /= tot;
      }
    }
  }

  void dinuc31_usage(double** stat_container) {
    int tot = 16;
    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < Nnuc; nuc2++) {
        stat_container[nuc1][nuc2] = 1.0;
      }
    }
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (site_codon < Nsite - 1) {
          int next_state =
              Data[taxa][site_codon + 1];  // GetState(taxa,site_codon+1);
          if (state != unknown && next_state != unknown) {
            stat_container[GetCodonStateSpace()->GetCodonPosition(2, state)]
                          [GetCodonStateSpace()->GetCodonPosition(
                              0, next_state)]++;
            tot++;
          }
        }
      }
    }

    for (int nuc1 = 0; nuc1 < Nnuc; nuc1++) {
      for (int nuc2 = 0; nuc2 < 4; nuc2++) {
        stat_container[nuc1][nuc2] /= tot;
      }
    }
  }

  // CpG+meandiff
  void CpG1_mean(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int* in = new int[Nsite];

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = 0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          int posa1 = GetCodonStateSpace()->GetCodonPosition(0, state);
          int posa2 = GetCodonStateSpace()->GetCodonPosition(1, state);
          if (posa1 == 1 && posa2 == 2) {
            in[site_codon]++;
          }
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      mean += static_cast<double>(in[site_codon]);
      var += static_cast<double>(in[site_codon] * in[site_codon]);
    }

    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);
    stat_container[0] = mean;
    stat_container[1] = var;

    delete[] in;
  }

  void CpG2_mean(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int* in = new int[Nsite];

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = 0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          int posa1 = GetCodonStateSpace()->GetCodonPosition(1, state);
          int posa2 = GetCodonStateSpace()->GetCodonPosition(2, state);
          if (posa1 == 1 && posa2 == 2) {
            in[site_codon]++;
          }
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      mean += static_cast<double>(in[site_codon]);
      var += static_cast<double>(in[site_codon] * in[site_codon]);
    }

    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);
    stat_container[0] = mean;
    stat_container[1] = var;

    delete[] in;
  }

  void CpG3_mean(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int* in = new int[Nsite];

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = 0;
    }

    for (int site_codon = 0; site_codon < Nsite - 1; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        int next_state =
            Data[taxa][site_codon + 1];  // GetState(taxa, site_codon);
        if (state != unknown) {
          int posa1 = GetCodonStateSpace()->GetCodonPosition(2, state);
          int posa2 = GetCodonStateSpace()->GetCodonPosition(0, next_state);
          if (posa1 == 1 && posa2 == 2) {
            in[site_codon]++;
          }
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      mean += static_cast<double>(in[site_codon]);
      var += static_cast<double>(in[site_codon] * in[site_codon]);
    }

    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);
    stat_container[0] = mean;
    stat_container[1] = var;

    delete[] in;
  }

  // CODON stat_containerS
  void codon_meandiff(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int** in = new int*[Nsite];
    int** inaa = new int*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = new int[GetCodonStateSpace()->GetNstate()];
      inaa[site_codon] = new int[Naa];
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        in[site_codon][codon] = 0;
      }
      for (int aa = 0; aa < Naa; aa++) {
        inaa[site_codon][aa] = 0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          in[site_codon][state] = 1;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int aa = 0; aa < Naa; aa++) {
        int NDiffsynPerAA = 0;
        for (int codon = 0; codon < GetCodonStateSpace()->GetNstate();
             codon++) {
          if (aa == GetCodonStateSpace()->Translation(codon) &&
              in[site_codon][codon] == 1) {
            NDiffsynPerAA++;
          }
        }
        if (NDiffsynPerAA > 0) {
          inaa[site_codon][aa] = NDiffsynPerAA - 1;
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      int Ndiff = 0;
      for (int aa = 0; aa < Naa; aa++) {
        Ndiff += inaa[site_codon][aa];
      }
      mean += static_cast<double>(Ndiff);
      var += static_cast<double>(Ndiff * Ndiff);
    }

    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);
    stat_container[0] = mean;
    stat_container[1] = var;

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] in[site_codon];
    }
    delete[] in;
  }

  void relativeCodonFrequency(double* stat_container) {
    int tot = GetCodonStateSpace()->GetNstate();
    for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
      stat_container[codon] = 1.0;
    }
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          stat_container[state]++;
          tot++;
        }
      }
    }
    for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
      stat_container[codon] /= tot;
    }
  }

  void CGNAGR(double* stat_container) {
    stat_container[0] = 1.0;
    stat_container[1] = 1.0;
    stat_container[2] = 1.0;
    stat_container[3] = 1.0;

    //    int Nsite_constant = 0;
    //    double meanRatio = 0.0;
    //    double varRatio = 0.0;

    int CGN_const = 0;
    int AGR_const = 0;
    int CGN_var = 0;
    int AGR_var = 0;

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      bool constant_site_R = true;
      int CGN_cur = 0;
      int AGR_cur = 0;
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        int CGT = 28;
        int CGC = 29;
        int CGA = 30;
        int CGG = 31;
        int AGA = 46;
        int AGG = 47;

        if (state != unknown) {
          if (state == CGT || state == CGC || state == CGA || state == CGG) {
            CGN_cur++;

          } else if (state == AGA || state == AGG) {
            AGR_cur++;

          } else {
            constant_site_R = false;
          }
        }
      }
      if (constant_site_R) {
        CGN_const += CGN_cur;
        AGR_const += AGR_cur;
        //            Nsite_constant++;
        //            meanRatio += (double) (CGN/AGR);
        //            varRatio += (double) (CGN/AGR) * (double) (CGN/AGR);
      } else {
        CGN_var += CGN_cur;
        AGR_var += AGR_cur;
      }
    }

    //    meanRatio /= Nsite_constant;
    //    varRatio = varRatio - (meanRatio * meanRatio);

    stat_container[0] += static_cast<double>(CGN_const);
    stat_container[1] += static_cast<double>(AGR_const);
    stat_container[2] += static_cast<double>(CGN_var);
    stat_container[3] += static_cast<double>(AGR_var);
  }

  /*    void CodonUsagePerAAwonR(double* stat_container)
     {
         //int tot = GetCodonStateSpace()->GetNstate();
         for ( int codon  = 0 ; codon < GetCodonStateSpace()->GetNstate(); codon
     ++ )
         {
             stat_container[codon] = 1.0 ;
         }


     //    int CCT = GetCodonStateSpace()->GetState("CCT"); // P
     //    int CCC = GetCodonStateSpace()->GetState("CCC"); // P
     //    int CCA = GetCodonStateSpace()->GetState("CCA"); // P
     //    int CCG = GetCodonStateSpace()->GetState("CCG"); // P

     //    int CTT = GetCodonStateSpace()->GetState("CTT"); // L
     //    int CTC = GetCodonStateSpace()->GetState("CTC"); // L
     //    int CTA = GetCodonStateSpace()->GetState("CTA"); // L
         int CTG = GetCodonStateSpace()->GetState("CTG"); // L

         int TGG = GetCodonStateSpace()->GetState("TGG"); // W

         int TCT = GetCodonStateSpace()->GetState("TCT"); // C
         int TGC = GetCodonStateSpace()->GetState("TGC"); // C

     //    int CAT = GetCodonStateSpace()->GetState("CAT"); // H
         int CAC = GetCodonStateSpace()->GetState("CAC"); // H

         int CAA = GetCodonStateSpace()->GetState("CAA"); // Q
         int CAG = GetCodonStateSpace()->GetState("CAG"); // Q

     //    int AGT = GetCodonStateSpace()->GetState("AGT"); // S
     //    int AGC = GetCodonStateSpace()->GetState("AGC"); // S

         int CGT = GetCodonStateSpace()->GetState("CGT"); // R
         int CGC = GetCodonStateSpace()->GetState("CGC"); // R
         int CGA = GetCodonStateSpace()->GetState("CGA"); // R
         int CGG = GetCodonStateSpace()->GetState("CGG"); // R

         int AGA = GetCodonStateSpace()->GetState("AGA"); // R
         int AGG = GetCodonStateSpace()->GetState("AGG"); // R

         for (int site_codon = 0 ; site_codon < Nsite; site_codon ++ )
         {
             bool test = true;
             int* cur = new int[GetCodonStateSpace()->GetNstate()];
             for (int codon = 0; codon <GetCodonStateSpace()->GetNstate();
     codon++)
             {
                 cur[codon] = 1;
             }
             for(int taxa = 0; taxa < Ntaxa; taxa++)
             {

                 int state = Data[taxa][site_codon]; //GetState(taxa,
     site_codon);


                 if (state != unknown)
                 {

                     if (
                         //state != CCT && state != CCC && state != CCA && state
     != CCG &&
                         //state != CTT && state != CTC && state != CTA && state
     != CAT && state != TGG && state != TCT && state != TGC && state != CTG &&
     state != CAC && state != CAA && state != CAG &&
                         //state != AGT && state != AGC
                         state != CGT && state != CGC && state != CGA && state
     != CGG && state != AGA && state != AGG
                     )
                     {
                         cur[state]++;

                     }
                     else
                     {
                         test = false;
                     }
                 }
             }
             if (test)
             {

                 for (int codon = 0; codon <Naa; codon++)
                 {
                     stat_container[codon] += cur[codon];
                 }
             }
             delete [] cur;
         }




         for (int aa = 0 ; aa < Naa ; aa++)
         {
             int aa_tot = 0;
             for (int codon = 0 ; codon < GetCodonStateSpace()->GetNstate();
     codon++)
             {
                 if (aa == GetCodonStateSpace()->Translation(codon))
                 {
                     aa_tot += stat_container[codon];
                 }
             }
             for (int codon = 0 ; codon < GetCodonStateSpace()->GetNstate();
     codon++)
             {
                 if (aa == GetCodonStateSpace()->Translation(codon))
                 {
                     stat_container[codon]/=aa_tot;
                 }
             }

         }

     }
  */

  double RSCUEntropy() {
    double* stat_container = new double[GetCodonStateSpace()->GetNstate()];

    int tot = GetCodonStateSpace()->GetNstate();
    for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
      stat_container[codon] = 1.0;
    }
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          stat_container[state]++;
          tot++;
        }
      }
    }
    for (int aa = 0; aa < Naa; aa++) {
      int aa_tot = 0;
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          aa_tot += stat_container[codon];
        }
      }
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          stat_container[codon] /= aa_tot;
        }
      }
    }

    double entropy = 0.0;
    for (int codon_i = 0; codon_i < GetCodonStateSpace()->GetNstate();
         codon_i++) {
      entropy -= stat_container[codon_i] * log2(stat_container[codon_i]);
    }

    delete[] stat_container;

    return entropy;
  }

  void RSCU(double* stat_container) {
    int tot = GetCodonStateSpace()->GetNstate();
    for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
      stat_container[codon] = 1.0;
    }
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          stat_container[state]++;
          tot++;
        }
      }
    }
    for (int aa = 0; aa < Naa; aa++) {
      int aa_tot = 0;
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          aa_tot += stat_container[codon];
        }
      }
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          stat_container[codon] /= aa_tot;
        }
      }
    }
  }
  /*

      void dicodon_usage(double** stat_container)
      {
          int tot =
     GetCodonStateSpace()->GetNstate()*GetCodonStateSpace()->GetNstate() ; for
     (int codon1 = 0 ; codon1 < GetCodonStateSpace()->GetNstate() ; codon1 ++ )
          {
              for (int codon2 = 0 ; codon2 < GetCodonStateSpace()->GetNstate();
     codon2 ++ )
              {
                  stat_container[codon1][codon2]  = 1.0 ;
              }
          }
          for (int taxa = 0; taxa < Ntaxa; taxa++)
          {
              for (int site_codon = 0; site_codon < Nsite; site_codon++)
              {
                  int state = Data[taxa][site_codon];
     //GetState(taxa,site_codon); if(site_codon < Nsite-1)
                  {
                      int next_state = Data[taxa][site_codon+1];
     //GetState(taxa,site_codon+1); if(state != unknown && next_state !=
     unknown)
                      {
                          stat_container[state][next_state]++;
                          tot++ ;
                      }
                  }
              }
          }

          for (int codon1 = 0 ; codon1 < GetCodonStateSpace()->GetNstate() ;
     codon1 ++ )
          {
              for (int codon2 = 0 ; codon2 < GetCodonStateSpace()->GetNstate();
     codon2 ++ )
              {
                  stat_container[codon1][codon2] /=  tot ;
              }
          }

      }
   */
  void aa_meandiff(double* stat_container) {
    stat_container[0] = 0.0;
    stat_container[1] = 0.0;
    int** in = new int*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      in[site_codon] = new int[Naa];
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int aa = 0; aa < Naa; aa++) {
        in[site_codon][aa] = 0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa, site_codon);
        if (state != unknown) {
          in[site_codon][GetCodonStateSpace()->Translation(state)] = 1;
        }
      }
    }

    double mean = 0.0;
    double var = 0.0;

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      int Ndiff = 0;
      for (int aa = 0; aa < Naa; aa++) {
        if (in[site_codon][aa] != 0) {
          Ndiff++;
        }
      }
      mean += static_cast<double>(Ndiff);
      var += static_cast<double>(Ndiff * Ndiff);
    }
    mean /= static_cast<double>(Nsite);
    var /= static_cast<double>(Nsite);
    var = var - (mean * mean);

    stat_container[0] = mean;
    stat_container[1] = var;

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] in[site_codon];
    }
    delete[] in;
  }

  void relativeAAFrequency(double* stat_container) {
    int tot = 20;
    for (int aa = 0; aa < Naa; aa++) {
      stat_container[aa] = 1.0;
    }
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          stat_container[GetCodonStateSpace()->Translation(state)]++;
          tot++;
        }
      }
    }
    for (int aa = 0; aa < Naa; aa++) {
      stat_container[aa] /= tot;
    }
  }

  /* void diaa_usage(double** stat_container)
  {
      int tot = 400 ;
      for (int aa1 = 0 ; aa1 < Naa ; aa1 ++ )
      {
          for (int aa2 = 0 ; aa2 < Naa; aa2 ++ )
          {
              stat_container[aa1][aa2]  = 1.0 ;
          }
      }
      for (int taxa = 0; taxa < Ntaxa; taxa++)
      {
          for (int site_codon = 0; site_codon < Nsite; site_codon++)
          {
              int state = Data[taxa][site_codon]; //GetState(taxa,site_codon);
              if(site_codon < Nsite-1)
              {
                  int next_state = Data[taxa][site_codon+1];
  //GetState(taxa,site_codon+1); if(state != unknown && next_state != unknown)
                  {
                      stat_container[GetCodonStateSpace()->Translation(state)][GetCodonStateSpace()->Translation(next_state)]++;
                      tot++ ;
                  }
              }
          }
      }

      for (int aa1 = 0 ; aa1 < Naa ; aa1 ++ )
      {
          for (int aa2 = 0 ; aa2 < Naa; aa2 ++ )
          {
              stat_container[aa1][aa2] /=  tot ;
          }
      }

  } */

  void dNdS_pairwise(int* invec) {
    for (int i = 0; i < 2; i++) {
      invec[i] = 0;
    }
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);

          if (state_seq1 != unknown && state_seq2 != unknown) {
            if (!GetCodonStateSpace()->Synonymous(state_seq1, state_seq2)) {
              invec[0]++;
            } else {
              invec[1]++;
            }
          }
        }
      }
    }
  }

  // pairwise stat_containers
  void dinucCpG_pairwise(int* invec) {
    for (int i = 0; i < 4; i++) {
      invec[i] = 1;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);

          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          int next_state_seq1 = -1;
          int next_state_seq2 = -1;
          if (site_codon < Nsite - 1) {
            next_state_seq1 = Data[taxa1][site_codon + 1];
            next_state_seq1 = Data[taxa1][site_codon + 1];
          }

          if (state_seq1 != unknown && state_seq2 != unknown) {
            for (int pos = 0; pos < 2; pos++) {
              int posa1 =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq1);
              int posa2 =
                  GetCodonStateSpace()->GetCodonPosition(pos + 1, state_seq1);
              int posb1 =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq2);
              int posb2 =
                  GetCodonStateSpace()->GetCodonPosition(pos + 1, state_seq2);

              if ((posa1 == 1 && posb1 == 2 && posa2 == 3 && posb2 == 2) ||
                  (posa1 == 3 && posb1 == 2 && posa2 == 1 &&
                   posb2 == 2))  // CpG<->TpG
              {
                invec[0]++;
                invec[3]++;
              } else if ((posa1 == 1 && posb1 == 2 && posa2 == 1 &&
                          posb2 == 0) ||
                         (posa1 == 1 && posb1 == 0 && posa2 == 1 &&
                          posb2 == 2))  // CpG<->CpA
              {
                invec[1]++;
                invec[3]++;
              } else if ((posa1 == 0 && posb1 == 2 && posa2 == 3 &&
                          posb2 == 2) ||
                         (posa1 == 3 && posb1 == 2 && posa2 == 0 &&
                          posb2 == 2))  // ApG<->TpG
              {
                invec[2]++;
                invec[3]++;
              }
            }
          }

          if (state_seq1 != unknown && state_seq2 != unknown &&
              next_state_seq1 != unknown && next_state_seq1 != unknown) {
            int posa1 = GetCodonStateSpace()->GetCodonPosition(2, state_seq1);
            int posa2 =
                GetCodonStateSpace()->GetCodonPosition(0, next_state_seq1);
            int posb1 = GetCodonStateSpace()->GetCodonPosition(2, state_seq2);
            int posb2 =
                GetCodonStateSpace()->GetCodonPosition(0, next_state_seq2);

            if ((posa1 == 1 && posb1 == 2 && posa2 == 3 && posb2 == 2) ||
                (posa1 == 3 && posb1 == 2 && posa2 == 1 &&
                 posb2 == 2))  // CpG<->TpG
            {
              invec[0]++;
              invec[2]++;
            } else if ((posa1 == 1 && posb1 == 2 && posa2 == 1 && posb2 == 0) ||
                       (posa1 == 1 && posb1 == 0 && posa2 == 1 &&
                        posb2 == 2))  // CpG<->CpA
            {
              invec[1]++;
              invec[2]++;
            }
          }
        }
      }
    }
  }

  // pairwise stat_containers
  /*     void nuc_pairwise(int* invec)
      {
          for (int i = 0 ; i < 7; i ++ )
          {
              invec[i] = 1 ;
          }

          for (int site_codon = 0; site_codon < Nsite; site_codon++)
          {
              for (int taxa1 = 0; taxa1 < Ntaxa-1; taxa1++)
              {
                  for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++)
                  {
                      int state_seq1 = Data[taxa1][site_codon];
     //GetState(taxa1,site_codon); int state_seq2 = Data[taxa2][site_codon];
     //GetState(taxa2,site_codon);

                      if(state_seq1 != unknown && state_seq2 != unknown)
                      {
                          for(int pos = 0 ; pos < 3 ; pos++ )
                          {

                              int posa =
     GetCodonStateSpace()->GetCodonPosition(pos,state_seq1); int posb =
     GetCodonStateSpace()->GetCodonPosition(pos,state_seq2);

                              if ((posa == 0 && posb == 1 ) || (posa == 1 &&
     posb == 0))  // ac_ca
                              {
                                  invec[0]++;
                                  //invec[6]++;
                              }
                              else if ((posa == 0 && posb == 2) || (posa == 2 &&
     posb == 0) )    //ag_ga
                              {
                                  invec[1]++;
                                  //invec[6]++;
                              }
                              else if ((posa == 0 && posb == 3 ) || (posa == 3
     && posb == 0))    //at_ta
                              {
                                  invec[2]++;
                                  //invec[6]++;
                              }
                              else if ((posa == 1 && posb == 2 ) || (posa == 2
     && posb == 1))    //cg_gc
                              {
                                  invec[3]++;
                                  //invec[6]++;
                              }
                              else if ( (posa == 1 && posb == 3) || (posa == 3
     && posb == 1))    //ct_tc
                              {
                                  invec[4]++;
                                  //invec[6]++;
                              }
                              else if ((posa == 3 && posb == 2 ) || (posa == 2
     && posb == 3))    //tg_gt
                              {
                                  invec[5]++;
                                  //invec[6]++;

                              }

                          }
                      }
                  }

              }

          }

      }
      void nuc1_pairwise(int* invec)
      {
          for (int i = 0 ; i < 7; i ++ )
          {
              invec[i] = 0 ;
          }
          for (int site_codon = 0; site_codon < Nsite; site_codon++)
          {
              for (int taxa1 = 0; taxa1 < Ntaxa-1; taxa1++)
              {
                  for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++)
                  {
                      int state_seq1 = Data[taxa1][site_codon];
     //GetState(taxa1,site_codon); int state_seq2 = Data[taxa2][site_codon];
     //GetState(taxa2,site_codon); if(state_seq1 != unknown && state_seq2 !=
     unknown)
                      {
                          int posa =
     GetCodonStateSpace()->GetCodonPosition(0,state_seq1); int posb =
     GetCodonStateSpace()->GetCodonPosition(0,state_seq2); if ((posa == 0 &&
     posb == 1 ) || (posa == 1 && posb == 0))  // ac_ca
                          {
                              invec[0]++;
                              //invec[6]++;
                          }
                          else if ((posa == 0 && posb == 2) || (posa == 2 &&
     posb == 0) )    //ag_ga
                          {
                              invec[1]++;
                              //invec[6]++;
                          }
                          else if ((posa == 0 && posb == 3 ) || (posa == 3 &&
     posb == 0))    //at_ta
                          {
                              invec[2]++;
                              //invec[6]++;
                          }
                          else if ((posa == 1 && posb == 2 ) || (posa == 2 &&
     posb == 1))    //cg_gc
                          {
                              invec[3]++;
                              //invec[6]++;
                          }
                          else if ( (posa == 1 && posb == 3) || (posa == 3 &&
     posb == 1))    //ct_tc
                          {
                              invec[4]++;
                              //invec[6]++;
                          }
                          else if ((posa == 3 && posb == 2 ) || (posa == 2 &&
     posb == 3))    //tg_gt
                          {
                              invec[5]++;
                              //invec[6]++;
                          }
                      }
                  }
              }
          }
      }
      void nuc2_pairwise(int* invec)
      {
          for (int i = 0 ; i < 7; i ++ )
          {
              invec[i] = 0;
          }
          for (int site_codon = 0; site_codon < Nsite; site_codon++)
          {
              for (int taxa1 = 0; taxa1 < Ntaxa-1; taxa1++)
              {
                  for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++)
                  {
                      int state_seq1 = Data[taxa1][site_codon];
     //GetState(taxa1,site_codon); int state_seq2 = Data[taxa2][site_codon];
     //GetState(taxa2,site_codon); if(state_seq1 != unknown && state_seq2 !=
     unknown)
                      {
                          int posa =
     GetCodonStateSpace()->GetCodonPosition(1,state_seq1); int posb =
     GetCodonStateSpace()->GetCodonPosition(1,state_seq2); if ((posa == 0 &&
     posb == 1 ) || (posa == 1 && posb == 0))  // ac_ca
                          {
                              invec[0]++;
                              //invec[6]++;
                          }
                          else if ((posa == 0 && posb == 2) || (posa == 2 &&
     posb == 0) )    //ag_ga
                          {
                              invec[1]++;
                              //invec[6]++;
                          }
                          else if ((posa == 0 && posb == 3 ) || (posa == 3 &&
     posb == 0))    //at_ta
                          {
                              invec[2]++;
                              //invec[6]++;
                          }
                          else if ((posa == 1 && posb == 2 ) || (posa == 2 &&
     posb == 1))    //cg_gc
                          {
                              invec[3]++;
                              //invec[6]++;
                          }
                          else if ( (posa == 1 && posb == 3) || (posa == 3 &&
     posb == 1))    //ct_tc
                          {
                              invec[4]++;
                              //invec[6]++;
                          }
                          else if ((posa == 3 && posb == 2 ) || (posa == 2 &&
     posb == 3))    //tg_gt
                          {
                              invec[5]++;
                              //invec[6]++;
                          }
                      }
                  }
              }
          }
      }
   */
  /* void nuc3_pairwise(int* invec)
  {
      for (int i = 0 ; i < 7; i ++ )
      {
          invec[i] = 0 ;
      }
      int cons = 1;

      for (int site_codon = 0; site_codon < Nsite; site_codon++)
      {
          for (int taxa1 = 0; taxa1 < Ntaxa-1; taxa1++)
          {
              for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++)
              {
                  int state_seq1 = Data[taxa1][site_codon];
  //GetState(taxa1,site_codon); int state_seq2 = Data[taxa2][site_codon];
  //GetState(taxa2,site_codon); if(state_seq1 != unknown && state_seq2 !=
  unknown)
                  {
                      int posa =
  GetCodonStateSpace()->GetCodonPosition(2,state_seq1); int posb =
  GetCodonStateSpace()->GetCodonPosition(2,state_seq2); if ((posa == 0 && posb
  == 1 ) || (posa == 1 && posb == 0))  // ca_ac
                      {
                          invec[0]++;
                          //invec[6]++;
                      }
                      else if ((posa == 0 && posb == 2) || (posa == 2 && posb ==
  0) )    //ag_ga
                      {
                          invec[1]++;
                          //invec[6]++;
                      }
                      else if ((posa == 0 && posb == 3 ) || (posa == 3 && posb
  == 0))    //at_ta
                      {
                          invec[2]++;
                          //invec[6]++;
                      }
                      else if ((posa == 1 && posb == 2 ) || (posa == 2 && posb
  == 1))    //cg_gc
                      {
                          invec[3]++;
                          //invec[6]++;
                      }
                      else if ( (posa == 1 && posb == 3) || (posa == 3 && posb
  == 1))    //ct_tc
                      {
                          invec[4]++;
                          //invec[6]++;
                      }
                      else if ((posa == 3 && posb == 2 ) || (posa == 2 && posb
  == 3))    //gt_tg
                      {
                          invec[5]++;
                          //invec[6]++;
                      }
                      else if (posa == posb){
                          cons++;

                      }
                  }
              }
          }
      }
      double all = 0 ;
      for (int i =0; i <6; i++)
      {
          all += (double) invec[i];
      }
      double p = (invec[1] + invec[3]) / all;
      double q = (invec[0]+invec[2]+invec[4]+invec[5]) / all;
      double d_K80 = -0.5*log((1-2*p-q)*sqrt(1-2*q));
      invec[6]

  } */

  void nuc_pairwise(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            for (int i = 0; i < 3; i++) {
              int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
              int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              } else if (posa == posb) {
                cons++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size();

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwiseSyn(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown &&
              GetCodonStateSpace()->Synonymous(state_seq1, state_seq2)) {
            int pos = GetCodonStateSpace()->GetDifferingPosition(state_seq1,
                                                                 state_seq2);
            if (pos == -1) {
              cons++;
            } else if (pos != 3) {
              int posa =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq1);
              int posb =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size();

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwiseSyn10(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown &&
              GetCodonStateSpace()->Synonymous(state_seq1, state_seq2)) {
            int pos = GetCodonStateSpace()->GetDifferingPosition(state_seq1,
                                                                 state_seq2);
            if (pos == -1) {
              cons++;
            } else if (pos != 3) {
              int posa =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq1);
              int posb =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.1;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwiseSyn30(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown &&
              GetCodonStateSpace()->Synonymous(state_seq1, state_seq2)) {
            int pos = GetCodonStateSpace()->GetDifferingPosition(state_seq1,
                                                                 state_seq2);
            if (pos == -1) {
              cons++;
            } else if (pos != 3) {
              int posa =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq1);
              int posb =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.3;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwiseSyn50(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown &&
              GetCodonStateSpace()->Synonymous(state_seq1, state_seq2)) {
            int pos = GetCodonStateSpace()->GetDifferingPosition(state_seq1,
                                                                 state_seq2);
            if (pos == -1) {
              cons++;
            } else if (pos != 3) {
              int posa =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq1);
              int posb =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.5;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwiseSyn70(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown &&
              GetCodonStateSpace()->Synonymous(state_seq1, state_seq2)) {
            int pos = GetCodonStateSpace()->GetDifferingPosition(state_seq1,
                                                                 state_seq2);
            if (pos == -1) {
              cons++;
            } else if (pos != 3) {
              int posa =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq1);
              int posb =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.7;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwiseSyn90(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown &&
              GetCodonStateSpace()->Synonymous(state_seq1, state_seq2)) {
            int pos = GetCodonStateSpace()->GetDifferingPosition(state_seq1,
                                                                 state_seq2);
            if (pos == -1) {
              cons++;
            } else if (pos != 3) {
              int posa =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq1);
              int posb =
                  GetCodonStateSpace()->GetCodonPosition(pos, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.9;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc1_pairwise(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 0;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size();

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc1_pairwise10(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 0;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.1;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc1_pairwise30(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 0;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.3;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc1_pairwise50(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 0;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.5;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc1_pairwise70(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 0;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.7;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc1_pairwise90(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 0;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.9;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc2_pairwise(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 1;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size();

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc2_pairwise10(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 1;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.1;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc2_pairwise30(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 1;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.3;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc2_pairwise50(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 1;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.5;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc2_pairwise70(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 1;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.7;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc2_pairwise90(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 1;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.9;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc3_pairwise(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 2;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size();

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc3_pairwise10(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 2;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double t1 = 1 - 2 * p - q;
        double t2 = 1 - 2 * q;
        if (t1 > 0 && t2 > 0) {
          double d_K80 = -0.5 * log(t1 * sqrt(t2));
          // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
          vec_ac.push_back(ac);
          vec_ag.push_back(ag);
          vec_at.push_back(at);
          vec_cg.push_back(cg);
          vec_ct.push_back(ct);
          vec_gt.push_back(gt);
          vec_d.push_back(d_K80);
        }
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.1;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc3_pairwise30(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 2;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.3;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc3_pairwise50(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 2;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.5;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc3_pairwise70(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 2;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.7;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc3_pairwise90(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            int i = 2;
            int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
            int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
            if ((posa == 0 && posb == 1) || (posa == 1 && posb == 0))  // ca_ac
            {
              ac++;
            } else if ((posa == 0 && posb == 2) ||
                       (posa == 2 && posb == 0))  // ag_ga
            {
              ag++;
            } else if ((posa == 0 && posb == 3) ||
                       (posa == 3 && posb == 0))  // at_ta
            {
              at++;
            } else if ((posa == 1 && posb == 2) ||
                       (posa == 2 && posb == 1))  // cg_gc
            {
              cg++;
            } else if ((posa == 1 && posb == 3) ||
                       (posa == 3 && posb == 1))  // ct_tc
            {
              ct++;
            } else if ((posa == 3 && posb == 2) ||
                       (posa == 2 && posb == 3))  // gt_tg
            {
              gt++;
            } else if (posa == posb) {
              cons++;
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.9;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwise10(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            for (int i = 0; i < 3; i++) {
              int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
              int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              } else if (posa == posb) {
                cons++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.10;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwise30(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            for (int i = 0; i < 3; i++) {
              int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
              int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              } else if (posa == posb) {
                cons++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.30;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwise50(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            for (int i = 0; i < 3; i++) {
              int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
              int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              } else if (posa == posb) {
                cons++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.50;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwise70(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            for (int i = 0; i < 3; i++) {
              int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
              int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              } else if (posa == posb) {
                cons++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.70;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }

  void nuc_pairwise90(double* invec) {
    for (int i = 0; i < 7; i++) {
      invec[i] = 0;
    }

    std::vector<double> vec_ac;
    std::vector<double> vec_ag;
    std::vector<double> vec_at;
    std::vector<double> vec_cg;
    std::vector<double> vec_ct;
    std::vector<double> vec_gt;
    std::vector<double> vec_d;

    for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
      for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
        double ac, ag, at, cg, ct, gt, cons;
        ac = ag = at = cg = ct = gt = cons = 0;
        for (int site_codon = 0; site_codon < Nsite; site_codon++) {
          int state_seq1 =
              Data[taxa1][site_codon];  // GetState(taxa1,site_codon);
          int state_seq2 =
              Data[taxa2][site_codon];  // GetState(taxa2,site_codon);
          if (state_seq1 != unknown && state_seq2 != unknown) {
            for (int i = 0; i < 3; i++) {
              int posa = GetCodonStateSpace()->GetCodonPosition(i, state_seq1);
              int posb = GetCodonStateSpace()->GetCodonPosition(i, state_seq2);
              if ((posa == 0 && posb == 1) ||
                  (posa == 1 && posb == 0))  // ca_ac
              {
                ac++;
              } else if ((posa == 0 && posb == 2) ||
                         (posa == 2 && posb == 0))  // ag_ga
              {
                ag++;
              } else if ((posa == 0 && posb == 3) ||
                         (posa == 3 && posb == 0))  // at_ta
              {
                at++;
              } else if ((posa == 1 && posb == 2) ||
                         (posa == 2 && posb == 1))  // cg_gc
              {
                cg++;
              } else if ((posa == 1 && posb == 3) ||
                         (posa == 3 && posb == 1))  // ct_tc
              {
                ct++;
              } else if ((posa == 3 && posb == 2) ||
                         (posa == 2 && posb == 3))  // gt_tg
              {
                gt++;
              } else if (posa == posb) {
                cons++;
              }
            }
          }
        }
        double all = ac + ag + at + cg + ct + gt + cons;
        double p = (ag + ct) / all;
        double q = (ac + at + cg + gt) / all;
        double d_K80 = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q));
        // double d_kimura  = -0.5*log(1-2*p-q)-0.25*log(1-2*q);
        vec_ac.push_back(ac);
        vec_ag.push_back(ag);
        vec_at.push_back(at);
        vec_cg.push_back(cg);
        vec_ct.push_back(ct);
        vec_gt.push_back(gt);
        vec_d.push_back(d_K80);
      }
    }

    std::sort(vec_ac.begin(), vec_ac.end());
    std::sort(vec_ag.begin(), vec_ag.end());
    std::sort(vec_at.begin(), vec_at.end());
    std::sort(vec_cg.begin(), vec_cg.end());
    std::sort(vec_ct.begin(), vec_ct.end());
    std::sort(vec_gt.begin(), vec_gt.end());
    std::sort(vec_d.begin(), vec_d.end());

    int size = (int)vec_ac.size() * 0.90;

    double sum_ac, sum_ag, sum_at, sum_cg, sum_ct, sum_gt, sum_all, sum_d_K80;
    sum_ac = sum_ag = sum_at = sum_cg = sum_ct = sum_gt = sum_all = sum_d_K80 =
        0;

    for (int i = 0; i < size; i++) {
      sum_ac += vec_ac[i];
      sum_ag += vec_ag[i];
      sum_at += vec_at[i];
      sum_cg += vec_cg[i];
      sum_ct += vec_ct[i];
      sum_gt += vec_gt[i];
      sum_d_K80 += vec_d[i];
    }
    invec[0] = sum_ac;
    invec[1] = sum_ag;
    invec[2] = sum_at;
    invec[3] = sum_cg;
    invec[4] = sum_ct;
    invec[5] = sum_gt;
    invec[6] = sum_d_K80;
  }
  void aa_pairwise10(double* invec) {
    for (int i = 0; i < 2; i++) {
      invec[i] = 0;
    }
    std::vector<double> vec_aa;
    std::vector<double> vec_K80;

    if (Ntaxa > 1) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          // added a pseudocount
          double aa = 1;
          double cons = 1;
          for (int site_codon = 0; site_codon < Nsite; site_codon++) {
            int state_seq1 = GetCodonStateSpace()->Translation(
                Data[taxa1][site_codon]);  // GetState(taxa1,site_codon);
            int state_seq2 = GetCodonStateSpace()->Translation(
                Data[taxa2][site_codon]);  // GetState(taxa2,site_codon);
            if (state_seq1 != unknown && state_seq2 != unknown) {
              if (state_seq1 != state_seq2) {
                aa++;
              } else {
                cons++;
              }
            }
          }
          if (aa == 0) {
            double d_K80 = 0;
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          } else {
            double freqVar = aa / (aa + cons);
            double d_K80 = -log(1 - freqVar - 0.2 * (freqVar * freqVar));
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          }
        }
      }
    }
    std::sort(vec_aa.begin(), vec_aa.end());
    int size = (int)vec_aa.size() * 0.1;
    double sum_aa = 0;
    double sum_K80 = 0;
    for (int i = 0; i < size; i++) {
      sum_aa += vec_aa[i];
      sum_K80 += vec_K80[i];
    }
    invec[0] = sum_aa;
    invec[1] = sum_K80;
  }

  void aa_pairwise30(double* invec) {
    for (int i = 0; i < 2; i++) {
      invec[i] = 0;
    }
    std::vector<double> vec_aa;
    std::vector<double> vec_K80;

    if (Ntaxa > 1) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          // added a pseudocount
          double aa = 1;
          double cons = 1;
          for (int site_codon = 0; site_codon < Nsite; site_codon++) {
            int state_seq1 = GetCodonStateSpace()->Translation(
                Data[taxa1][site_codon]);  // GetState(taxa1,site_codon);
            int state_seq2 = GetCodonStateSpace()->Translation(
                Data[taxa2][site_codon]);  // GetState(taxa2,site_codon);
            if (state_seq1 != unknown && state_seq2 != unknown) {
              if (state_seq1 != state_seq2) {
                aa++;
              } else {
                cons++;
              }
            }
          }
          if (aa == 0) {
            double d_K80 = 0;
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          } else {
            double freqVar = aa / (aa + cons);
            double d_K80 = -log(1 - freqVar - 0.2 * (freqVar * freqVar));
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          }
        }
      }
    }
    std::sort(vec_aa.begin(), vec_aa.end());
    int size = (int)vec_aa.size() * 0.3;
    double sum_aa = 0;
    double sum_K80 = 0;
    for (int i = 0; i < size; i++) {
      sum_aa += vec_aa[i];
      sum_K80 += vec_K80[i];
    }
    invec[0] = sum_aa;
    invec[1] = sum_K80;
  }
  void aa_pairwise50(double* invec) {
    for (int i = 0; i < 2; i++) {
      invec[i] = 0;
    }
    std::vector<double> vec_aa;
    std::vector<double> vec_K80;

    if (Ntaxa > 1) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          // added a pseudocount
          double aa = 1;
          double cons = 1;
          for (int site_codon = 0; site_codon < Nsite; site_codon++) {
            int state_seq1 = GetCodonStateSpace()->Translation(
                Data[taxa1][site_codon]);  // GetState(taxa1,site_codon);
            int state_seq2 = GetCodonStateSpace()->Translation(
                Data[taxa2][site_codon]);  // GetState(taxa2,site_codon);
            if (state_seq1 != unknown && state_seq2 != unknown) {
              if (state_seq1 != state_seq2) {
                aa++;
              } else {
                cons++;
              }
            }
          }
          if (aa == 0) {
            double d_K80 = 0;
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          } else {
            double freqVar = aa / (aa + cons);
            double d_K80 = -log(1 - freqVar - 0.2 * (freqVar * freqVar));
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          }
        }
      }
    }
    std::sort(vec_aa.begin(), vec_aa.end());
    int size = (int)vec_aa.size() * 0.5;
    double sum_aa = 0;
    double sum_K80 = 0;
    for (int i = 0; i < size; i++) {
      sum_aa += vec_aa[i];
      sum_K80 += vec_K80[i];
    }
    invec[0] = sum_aa;
    invec[1] = sum_K80;
  }

  void aa_pairwise70(double* invec) {
    for (int i = 0; i < 2; i++) {
      invec[i] = 0;
    }
    std::vector<double> vec_aa;
    std::vector<double> vec_K80;

    if (Ntaxa > 1) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          // added a pseudocount
          double aa = 1;
          double cons = 1;
          for (int site_codon = 0; site_codon < Nsite; site_codon++) {
            int state_seq1 = GetCodonStateSpace()->Translation(
                Data[taxa1][site_codon]);  // GetState(taxa1,site_codon);
            int state_seq2 = GetCodonStateSpace()->Translation(
                Data[taxa2][site_codon]);  // GetState(taxa2,site_codon);
            if (state_seq1 != unknown && state_seq2 != unknown) {
              if (state_seq1 != state_seq2) {
                aa++;
              } else {
                cons++;
              }
            }
          }
          if (aa == 0) {
            double d_K80 = 0;
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          } else {
            double freqVar = aa / (aa + cons);
            double d_K80 = -log(1 - freqVar - 0.2 * (freqVar * freqVar));
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          }
        }
      }
    }
    std::sort(vec_aa.begin(), vec_aa.end());
    int size = (int)vec_aa.size() * 0.7;
    double sum_aa = 0;
    double sum_K80 = 0;
    for (int i = 0; i < size; i++) {
      sum_aa += vec_aa[i];
      sum_K80 += vec_K80[i];
    }
    invec[0] = sum_aa;
    invec[1] = sum_K80;
  }
  void aa_pairwise90(double* invec) {
    for (int i = 0; i < 2; i++) {
      invec[i] = 0;
    }
    std::vector<double> vec_aa;
    std::vector<double> vec_K80;

    if (Ntaxa > 1) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          // added a pseudocount
          double aa = 1;
          double cons = 1;
          for (int site_codon = 0; site_codon < Nsite; site_codon++) {
            int state_seq1 = GetCodonStateSpace()->Translation(
                Data[taxa1][site_codon]);  // GetState(taxa1,site_codon);
            int state_seq2 = GetCodonStateSpace()->Translation(
                Data[taxa2][site_codon]);  // GetState(taxa2,site_codon);
            if (state_seq1 != unknown && state_seq2 != unknown) {
              if (state_seq1 != state_seq2) {
                aa++;
              } else {
                cons++;
              }
            }
          }
          if (aa == 0) {
            double d_K80 = 0;
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          } else {
            double freqVar = aa / (aa + cons);
            double d_K80 = -log(1 - freqVar - 0.2 * (freqVar * freqVar));
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          }
        }
      }
    }
    std::sort(vec_aa.begin(), vec_aa.end());
    int size = (int)vec_aa.size() * 0.9;
    double sum_aa = 0;
    double sum_K80 = 0;
    for (int i = 0; i < size; i++) {
      sum_aa += vec_aa[i];
      sum_K80 += vec_K80[i];
    }
    invec[0] = sum_aa;
    invec[1] = sum_K80;
  }

  void aa_pairwise(double* invec) {
    for (int i = 0; i < 2; i++) {
      invec[i] = 0;
    }
    std::vector<double> vec_aa;
    std::vector<double> vec_K80;

    if (Ntaxa > 1) {
      for (int taxa1 = 0; taxa1 < Ntaxa - 1; taxa1++) {
        for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++) {
          // added a pseudocount
          double aa = 1;
          double cons = 1;
          for (int site_codon = 0; site_codon < Nsite; site_codon++) {
            int state_seq1 = GetCodonStateSpace()->Translation(
                Data[taxa1][site_codon]);  // GetState(taxa1,site_codon);
            int state_seq2 = GetCodonStateSpace()->Translation(
                Data[taxa2][site_codon]);  // GetState(taxa2,site_codon);
            if (state_seq1 != unknown && state_seq2 != unknown) {
              if (state_seq1 != state_seq2) {
                aa++;
              } else {
                cons++;
              }
            }
          }
          if (aa == 0) {
            double d_K80 = 0;
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          } else {
            double freqVar = aa / (aa + cons);
            double d_K80 = -log(1 - freqVar - 0.2 * (freqVar * freqVar));
            vec_K80.push_back(d_K80);
            vec_aa.push_back(aa);
          }
        }
      }
    }
    std::sort(vec_aa.begin(), vec_aa.end());
    int size = (int)vec_aa.size();
    double sum_aa = 0;
    double sum_K80 = 0;
    for (int i = 0; i < size; i++) {
      sum_aa += vec_aa[i];
      sum_K80 += vec_K80[i];
    }
    invec[0] = sum_aa;
    invec[1] = sum_K80;
  }

  /* double aa_pairwise()
  {
      double sum = 0;
      if (Ntaxa > 1)
      {
          int ** m = new int* [Naa] ;
          for(int aa = 0 ; aa < Naa; aa ++ )
          {
              m[aa] = new int [Naa];
          }

          for(int aa1 = 0 ; aa1 < Naa; aa1 ++ )
          {
              for(int aa2 = 0 ; aa2 < Naa; aa2 ++ )
              {
                  m[aa1][aa2] = 0 ;
              }
          }
          for (int site_codon = 0; site_codon < Nsite; site_codon++)
          {
              for (int taxa1 = 0; taxa1 < Ntaxa-1; taxa1++)
              {
                  for (int taxa2 = taxa1 + 1; taxa2 < Ntaxa; taxa2++)
                  {
                      int state_seq1 =
  GetCodonStateSpace()->Translation(Data[taxa1][site_codon]);
  //GetState(taxa1,site_codon); int state_seq2 =
  GetCodonStateSpace()->Translation(Data[taxa2][site_codon]);
  //GetState(taxa2,site_codon); if(state_seq1 != unknown && state_seq2 !=
  unknown && state_seq1 != state_seq2)
                      {
                          m[state_seq1][state_seq2]++;
                      }
                  }
              }
          }

          for(int aa1 = 0 ; aa1 < Naa-1; aa1 ++ )
          {
              for(int aa2 = aa1+1 ; aa2 < Naa; aa2 ++ )
              {
                  sum += m[aa1][aa2]+m[aa2][aa1];
              }
          }

          for(int aa = 0 ; aa < Naa; aa ++ )
          {
              delete [] m[aa];
          }
          delete [] m;
      }
      return sum;
  }
*/

  // Compositional heterogeneity

  double dinuc_site_comphet() {
    double** sitefreq = new double*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      sitefreq[site_codon] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(0, state)]++;
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(2, state)]++;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += sitefreq[site_codon][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += sitefreq[site_codon][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Nsite;
    }

    double dist = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (sitefreq[site_codon][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }
    for (int j = 0; j < Nsite; j++) {
      delete[] sitefreq[j];
    }
    delete[] sitefreq;
    delete[] globalfreq;
    return dist / Nsite;
  }

  double nuc_site_comphet() {
    double** sitefreq = new double*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      sitefreq[site_codon] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] = 1.0;
      }
    }
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(0, state)]++;
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(2, state)]++;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += sitefreq[site_codon][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += sitefreq[site_codon][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Nsite;
    }

    double dist = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (sitefreq[site_codon][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete sitefreq[site_codon];
    }

    delete[] sitefreq;
    delete[] globalfreq;

    dist /= Nsite;
    return dist;
  }

  double nuc1_site_comphet() {
    double** sitefreq = new double*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      sitefreq[site_codon] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(0, state)]++;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += sitefreq[site_codon][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += sitefreq[site_codon][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Nsite;
    }

    double dist = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (sitefreq[site_codon][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] sitefreq[site_codon];
    }
    delete[] sitefreq;
    delete[] globalfreq;

    dist /= Nsite;
    return dist;
  }

  double nuc2_site_comphet() {
    double** sitefreq = new double*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      sitefreq[site_codon] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(1, state)]++;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += sitefreq[site_codon][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += sitefreq[site_codon][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Nsite;
    }

    double dist = 0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (sitefreq[site_codon][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] sitefreq[site_codon];
    }
    delete[] sitefreq;
    delete[] globalfreq;
    dist /= Nsite;
    return dist;
  }

  double nuc3_site_comphet() {
    double** sitefreq = new double*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      sitefreq[site_codon] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          sitefreq[site_codon]
                  [GetCodonStateSpace()->GetCodonPosition(2, state)]++;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += sitefreq[site_codon][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        sitefreq[site_codon][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += sitefreq[site_codon][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Nsite;
    }

    double dist = 0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (sitefreq[site_codon][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] sitefreq[site_codon];
    }
    delete[] sitefreq;
    delete[] globalfreq;
    dist /= Nsite;
    return dist;
  }

  double nuc_taxa_comphet() {
    double** taxafreq = new double*[Ntaxa];
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      taxafreq[taxa] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          taxafreq[taxa][GetCodonStateSpace()->GetCodonPosition(0, state)]++;
          taxafreq[taxa][GetCodonStateSpace()->GetCodonPosition(1, state)]++;
          taxafreq[taxa][GetCodonStateSpace()->GetCodonPosition(2, state)]++;
        }
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += taxafreq[taxa][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += taxafreq[taxa][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Ntaxa;
    }

    // compute max distance
    double dist = 0.0;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (taxafreq[taxa][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      delete[] taxafreq[taxa];
    }
    delete[] taxafreq;
    delete[] globalfreq;

    dist /= Ntaxa;
    return dist;
  }

  double nuc1_taxa_comphet() {
    double** taxafreq = new double*[Ntaxa];
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      taxafreq[taxa] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          taxafreq[taxa][GetCodonStateSpace()->GetCodonPosition(0, state)]++;
        }
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += taxafreq[taxa][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += taxafreq[taxa][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Ntaxa;
    }

    // compute max distance
    double dist = 0.0;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (taxafreq[taxa][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      delete[] taxafreq[taxa];
    }
    delete[] taxafreq;
    delete[] globalfreq;

    dist /= Ntaxa;
    return dist;
  }
  double nuc2_taxa_comphet() {
    double** taxafreq = new double*[Ntaxa];
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      taxafreq[taxa] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          taxafreq[taxa][GetCodonStateSpace()->GetCodonPosition(1, state)]++;
        }
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += taxafreq[taxa][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += taxafreq[taxa][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Ntaxa;
    }

    // compute max distance
    double dist = 0.0;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (taxafreq[taxa][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      delete[] taxafreq[taxa];
    }
    delete[] taxafreq;
    delete[] globalfreq;

    dist /= Ntaxa;
    return dist;
  }
  double nuc3_taxa_comphet() {
    double** taxafreq = new double*[Ntaxa];
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      taxafreq[taxa] = new double[Nnuc];
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          taxafreq[taxa][GetCodonStateSpace()->GetCodonPosition(2, state)]++;
        }
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      double total = 0.0;
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        total += taxafreq[taxa][nuc];
      }
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        taxafreq[taxa][nuc] /= total;
      }
    }

    double* globalfreq = new double[Nnuc];
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] = 0.0;
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        globalfreq[nuc] += taxafreq[taxa][nuc];
      }
    }
    for (int nuc = 0; nuc < Nnuc; nuc++) {
      globalfreq[nuc] /= Ntaxa;
    }

    // compute max distance
    double dist = 0.0;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int nuc = 0; nuc < Nnuc; nuc++) {
        double tmp = (taxafreq[taxa][nuc] - globalfreq[nuc]);
        dist += tmp * tmp;
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      delete[] taxafreq[taxa];
    }
    delete[] taxafreq;
    delete[] globalfreq;

    dist /= Ntaxa;
    return dist;
  }

  double codon_site_comphet() {
    double** sitefreq = new double*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      sitefreq[site_codon] = new double[GetCodonStateSpace()->GetNstate()];
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        sitefreq[site_codon][codon] = 1.0;
      }
    }

    // make global freqs out of tax-specific freqs
    double* globalfreq = new double[GetCodonStateSpace()->GetNstate()];
    for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
      globalfreq[codon] = 1.0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          sitefreq[site_codon][state]++;
          globalfreq[state]++;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int aa = 0; aa < Naa; aa++) {
        double aa_tot = 0;
        for (int codon = 0; codon < GetCodonStateSpace()->GetNstate();
             codon++) {
          if (aa == GetCodonStateSpace()->Translation(codon)) {
            aa_tot += sitefreq[site_codon][codon];
          }
        }
        for (int codon = 0; codon < GetCodonStateSpace()->GetNstate();
             codon++) {
          if (aa == GetCodonStateSpace()->Translation(codon)) {
            sitefreq[site_codon][codon] /= aa_tot;
          }
        }
      }
    }

    for (int aa = 0; aa < Naa; aa++) {
      double aa_tot = 0;
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          aa_tot += globalfreq[codon];
        }
      }
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          globalfreq[codon] /= aa_tot;
        }
      }
    }

    // compute sum of squared
    double dist = 0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        double tmp = (sitefreq[site_codon][codon] - globalfreq[codon]);
        dist += tmp * tmp;
      }
    }

    delete[] globalfreq;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] sitefreq[site_codon];
    }
    delete[] sitefreq;

    dist /= Nsite;
    return dist;
  }

  double codon_taxa_comphet() {
    double** taxafreq = new double*[Ntaxa];
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      taxafreq[taxa] = new double[GetCodonStateSpace()->GetNstate()];
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        taxafreq[taxa][codon] = 1.0;
      }
    }

    // make global freqs out of tax-specific freqs
    double* globalfreq = new double[GetCodonStateSpace()->GetNstate()];
    for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
      globalfreq[codon] = 1.0;
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          taxafreq[taxa][state]++;
          globalfreq[state]++;
        }
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int aa = 0; aa < Naa; aa++) {
        double aa_tot = 0;
        for (int codon = 0; codon < GetCodonStateSpace()->GetNstate();
             codon++) {
          if (aa == GetCodonStateSpace()->Translation(codon)) {
            aa_tot += taxafreq[taxa][codon];
          }
        }
        for (int codon = 0; codon < GetCodonStateSpace()->GetNstate();
             codon++) {
          if (aa == GetCodonStateSpace()->Translation(codon)) {
            taxafreq[taxa][codon] /= aa_tot;
          }
        }
      }
    }

    for (int aa = 0; aa < Naa; aa++) {
      double aa_tot = 0;
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          aa_tot += globalfreq[codon];
        }
      }
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        if (aa == GetCodonStateSpace()->Translation(codon)) {
          globalfreq[codon] /= aa_tot;
        }
      }
    }

    // compute sum of squared
    double dist = 0;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int codon = 0; codon < GetCodonStateSpace()->GetNstate(); codon++) {
        double tmp = (taxafreq[taxa][codon] - globalfreq[codon]);
        dist += tmp * tmp;
      }
    }

    delete[] globalfreq;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      delete[] taxafreq[taxa];
    }
    delete[] taxafreq;

    dist /= Ntaxa;
    return dist;
  }

  double aa_site_comphet() {
    double** sitefreq = new double*[Nsite];
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      sitefreq[site_codon] = new double[Naa];
      for (int aa = 0; aa < Naa; aa++) {
        sitefreq[site_codon][aa] = 1.0;
      }
    }
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          sitefreq[site_codon][GetCodonStateSpace()->Translation(state)]++;
        }
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      double total = 0.0;
      for (int aa = 0; aa < Naa; aa++) {
        total += sitefreq[site_codon][aa];
      }
      for (int aa = 0; aa < Naa; aa++) {
        sitefreq[site_codon][aa] /= total;
      }
    }

    // make global freqs out of tax-specific freqs
    double* globalfreq = new double[Naa];
    for (int aa = 0; aa < Naa; aa++) {
      globalfreq[aa] = 0.0;
      for (int site_codon = 0; site_codon < Nsite; site_codon++) {
        globalfreq[aa] += sitefreq[site_codon][aa];
      }
    }
    for (int aa = 0; aa < Naa; aa++) {
      globalfreq[aa] /= Nsite;
    }

    // compute sum of squared
    double dist = 0.0;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int aa = 0; aa < Naa; aa++) {
        double tmp = (sitefreq[site_codon][aa] - globalfreq[aa]);
        dist += tmp * tmp;
      }
    }

    delete[] globalfreq;
    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      delete[] sitefreq[site_codon];
    }
    delete[] sitefreq;

    dist /= Nsite;
    return dist;
  }

  double aa_taxa_comphet() {
    double** taxafreq = new double*[Ntaxa];
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      taxafreq[taxa] = new double[Naa];
      for (int aa = 0; aa < Naa; aa++) {
        taxafreq[taxa][aa] = 1.0;
      }
    }

    for (int site_codon = 0; site_codon < Nsite; site_codon++) {
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        int state = Data[taxa][site_codon];  // GetState(taxa,site_codon);
        if (state != unknown) {
          taxafreq[taxa][GetCodonStateSpace()->Translation(state)]++;
        }
      }
    }

    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      double total = 0.0;
      for (int aa = 0; aa < Naa; aa++) {
        total += taxafreq[taxa][aa];
      }
      for (int aa = 0; aa < Naa; aa++) {
        taxafreq[taxa][aa] /= total;
      }
    }

    // make global freqs out of tax-specific freqs
    double* globalfreq = new double[Naa];
    for (int aa = 0; aa < Naa; aa++) {
      globalfreq[aa] = 0.0;
      for (int taxa = 0; taxa < Ntaxa; taxa++) {
        globalfreq[aa] += taxafreq[taxa][aa];
      }
    }
    for (int aa = 0; aa < Naa; aa++) {
      globalfreq[aa] /= Ntaxa;
    }

    // compute sum of squared
    double dist = 0;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      for (int aa = 0; aa < Naa; aa++) {
        double tmp = (taxafreq[taxa][aa] - globalfreq[aa]);
        dist += tmp * tmp;
      }
    }

    delete[] globalfreq;
    for (int taxa = 0; taxa < Ntaxa; taxa++) {
      delete[] taxafreq[taxa];
    }
    delete[] taxafreq;

    dist /= Ntaxa;
    return dist;
  }

  //<--SLL

  CodonStateSpace* GetCodonStateSpace() {
    // return static_cast<CodonStateSpace*>(statespace);
    return (CodonStateSpace*)(statespace);
  }

  void ToStream(ostream& os);

 private:
  SequenceAlignment* DNAsource;
};

/*
class GCContinuousData : public ContinuousData {

        public:

        GCContinuousData(CodonSequenceAlignment* from, int pos)	{
                taxset = from->GetTaxonSet();
                double** freq = new double*[taxset->GetNtaxa()];
                for (int i=0; i<taxset->GetNtaxa(); i++)	{
                        freq[i] = new double[Nnuc];
                }
                from->NucleotideCompositionalHeterogeneity(0,pos,freq);
                Data = new double*[taxset->GetNtaxa()];
                Nsite = 1;
                for (int i=0; i<taxset->GetNtaxa(); i++)	{
                        Data[i] = new double[1];
                        double tmp = freq[i][1] + freq[i][2];
                        Data[i][0] = tmp;
                        // Data[i][0] = log(tmp / (1-tmp));
                        cerr << taxset->GetTaxon(i) << '\t' << tmp << '\t' <<
Data[i][0] << '\n';
                }
        }
};
*/

#endif
