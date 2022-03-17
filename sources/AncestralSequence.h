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
#ifndef SOURCES_ANCESTRALSEQUENCE_H_
#define SOURCES_ANCESTRALSEQUENCE_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CodonSequenceAlignment.h"
#include "CodonStateSpace.h"
#include "LocalParameters.h"
#include "Random.h"
#include "SequenceAlignment.h"
#include "Tree.h"

class AncestralSequence {
 public:
  double** CurrentStationaryCodonSequence;
  int* CurrentCodonSequence;
  int* CurrentNucSequence;
  int* CurrentAncestralCodonSequence;
  int* CurrentAncestralNucSequence;
  int choosen_taxa;
  LocalParameters* lparam;

  explicit AncestralSequence(LocalParameters* lparam);
  ~AncestralSequence();

  void WriteStationaryCodonSequence();
  void ComputeStationaryCodon();
  void SampleAncestralCodonSequenceFromLeaves();
  void SampleAncestralCodonSequenceFromStationary();
  void SampleAncestralCodonSequenceFromStationaryCodon(int site_codon);
  int GetCurrentAncestralCodonSequence(int site_codon);
  int GetCurrentAncestralNucSequence(int site_nuc);

 protected:
 private:
};

#endif  // SOURCES_ANCESTRALSEQUENCE_H_
