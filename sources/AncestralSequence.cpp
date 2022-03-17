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
#include "AncestralSequence.h"

AncestralSequence::AncestralSequence(LocalParameters* lparam) {
  this->lparam = lparam;
  this->CurrentStationaryCodonSequence = new double*[lparam->Nsite_codon];
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    this->CurrentStationaryCodonSequence[site_codon] =
        new double[lparam->Nstate_codon];
  }
  this->CurrentCodonSequence = new int[lparam->Nsite_codon];
  this->CurrentAncestralCodonSequence = new int[lparam->Nsite_codon];
  this->CurrentNucSequence = new int[lparam->Nsite_codon * 3];
  this->CurrentAncestralNucSequence = new int[lparam->Nsite_codon * 3];
}

AncestralSequence::~AncestralSequence() {
  // dtor
}

void AncestralSequence::WriteStationaryCodonSequence() {
  // Write stationary
  ofstream stationary_os((lparam->output + ".stat").c_str(),
                         std::ios_base::app);
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    for (int state = 0; state < lparam->Nstate_codon; state++) {
      if (state < lparam->Nstate_codon - 1) {
        stationary_os << this->CurrentStationaryCodonSequence[site_codon][state]
                      << "\t";
      } else {
        stationary_os << this->CurrentStationaryCodonSequence[site_codon][state]
                      << "\n";
      }
    }
  }
  stationary_os << "\n";
  stationary_os.close();
}

void AncestralSequence::ComputeStationaryCodon() {
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    double Z = 0.0;
    for (int state = 0; state < lparam->Nstate_codon; state++) {
      this->CurrentStationaryCodonSequence[site_codon][state] =
          lparam->nucp[lparam->codonstatespace->GetCodonPosition(0, state)] *
          lparam->nucp[lparam->codonstatespace->GetCodonPosition(1, state)] *
          lparam->nucp[lparam->codonstatespace->GetCodonPosition(2, state)] *
          lparam->ssaaprofiles[lparam->alloc[site_codon]]
                              [lparam->codonstatespace->Translation(state)] *
          lparam->codonprofile[state];

      Z += CurrentStationaryCodonSequence[site_codon][state];
    }
    for (int state = 0; state < lparam->Nstate_codon; state++) {
      CurrentStationaryCodonSequence[site_codon][state] /= Z;
    }
  }
}

void AncestralSequence::SampleAncestralCodonSequenceFromStationaryCodon(
    int site_codon) {
  double u = lparam->rnd->Uniform();
  int state = 0;
  double testcummul = this->CurrentStationaryCodonSequence[site_codon][state];
  while (testcummul < u) {
    state++;
    testcummul += this->CurrentStationaryCodonSequence[site_codon][state];
  }
  if (state >= lparam->Nstate_codon) {
    std::cerr << "failed to draw stationary\n";
    exit(0);
  }
  this->CurrentCodonSequence[site_codon] = state;
  for (int j = 0; j < 3; j++) {
    this->CurrentNucSequence[site_codon * 3 + j] =
        lparam->codonstatespace->GetCodonPosition(j, state);
  }
}

void AncestralSequence::SampleAncestralCodonSequenceFromStationaryCodon() {
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    SampleAncestralCodonSequenceFromStationaryCodon(site_codon);
  }
}

void AncestralSequence::SampleAncestralCodonSequenceFromCodonData() {
  int** cur_data = lparam->codondata->GetData();
  double u = lparam->rnd->Uniform();
  int choosen_taxa = static_cast<int>((lparam->Ntaxa - 1) * u);
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    int codon_state = cur_data[choosen_taxa][site_codon];

    if (codon_state == unknown) {
      SampleAncestralCodonSequenceFromStationaryCodon(site_codon);
      codon_state = GetCurrentAncestralCodonSequence(site_codon);
    }

    this->CurrentCodonSequence[site_codon] = codon_state;
    for (int j = 0; j < 3; j++) {
      this->CurrentNucSequence[site_codon * 3 + j] =
          lparam->codonstatespace->GetCodonPosition(j, codon_state);
    }
  }

  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    delete[] cur_data[taxa_i];
  }
  delete[] cur_data;
}

int AncestralSequence::GetCurrentAncestralCodonSequence(int site_codon) {
  return this->CurrentCodonSequence[site_codon];
}

int AncestralSequence::GetCurrentAncestralNucSequence(int site_nuc) {
  return this->CurrentNucSequence[site_nuc];
}
