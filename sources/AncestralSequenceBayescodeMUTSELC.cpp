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
#include "AncestralSequenceBayescodeMUTSELC.h"

void AncestralSequenceBayescodeMUTSELC::ComputeStationaryCodon() {
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    double Z = 0.0;
    for (int state = 0; state < lparam->Nstate_codon; state++) {
      this->CurrentStationaryCodonSequence[site_codon][state] =
          lparam->nucp[lparam->codonstatespace->GetCodonPosition(0, state)] *
          lparam->nucp[lparam->codonstatespace->GetCodonPosition(1, state)] *
          lparam->nucp[lparam->codonstatespace->GetCodonPosition(2, state)] *
          lparam->sscodonprofiles[lparam->alloc[site_codon]][state];
        
      Z += CurrentStationaryCodonSequence[site_codon][state];
    }
    for (int state = 0; state < lparam->Nstate_codon; state++) {
      CurrentStationaryCodonSequence[site_codon][state] /= Z;
    }
  }
}