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
#include "TreeSimulator.h"

TreeSimulator::TreeSimulator(LocalParameters* lparam,
                             SiteInterSubMatrix* submatrix,
                             AncestralSequence* ancestralseq) {
  this->lparam = lparam;
  this->submatrix = submatrix;
  this->ancestralseq = ancestralseq;
  setSimulator();
}

TreeSimulator::TreeSimulator(LocalParameters* lparam,
                             SiteInterSubMatrix* submatrix) {
  this->lparam = lparam;
  this->submatrix = submatrix;
  setSimulatorFromLeaves();
}

TreeSimulator::~TreeSimulator() {
  // dtor
}

void TreeSimulator::setSimulatorFromLeaves() {
  CurrentLeafNodeCodonSequences = new int*[lparam->Ntaxa];
  CurrentLeafNodeNucSequence = new int*[lparam->Ntaxa];

  for (int taxa = 0; taxa < lparam->Ntaxa; taxa++) {
    CurrentLeafNodeCodonSequences[taxa] = new int[lparam->Nsite_codon];
    CurrentLeafNodeNucSequence[taxa] = new int[lparam->Nsite_nuc];
  }
}
void TreeSimulator::setSimulator() {
  this->treeEvoStats = new EvolHistStatistics(this->lparam);
  this->rootBranchEvoStats = new EvolHistStatistics(this->lparam);
  CurrentNodeCodonSequence = new int*[lparam->refTree->GetNnode()];
  CurrentNodeNucSequence = new int*[lparam->refTree->GetNnode()];

  for (int node = 0; node < lparam->refTree->GetNnode(); node++) {
    CurrentNodeCodonSequence[node] = new int[lparam->Nsite_codon];
    CurrentNodeNucSequence[node] = new int[lparam->Nsite_codon * 3];
  }

  CurrentLeafNodeCodonSequences = new int*[lparam->Ntaxa];
  CurrentLeafNodeNucSequence = new int*[lparam->Ntaxa];

  for (int taxa = 0; taxa < lparam->Ntaxa; taxa++) {
    CurrentLeafNodeCodonSequences[taxa] = new int[lparam->Nsite_codon];
    CurrentLeafNodeNucSequence[taxa] = new int[lparam->Nsite_nuc];
  }

  CurrentAncestralCodonSequence = new int**[lparam->Ninterval];
  for (int point_i = 0; point_i < lparam->Ninterval; point_i++) {
    CurrentAncestralCodonSequence[point_i] = new int*[1];
    CurrentAncestralCodonSequence[point_i][0] = new int[lparam->Nsite_codon];
  }
}

void TreeSimulator::resetSimulator() {
  for (int node = 0; node < lparam->refTree->GetNnode(); node++) {
    for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
      CurrentNodeCodonSequence[node][site_codon] = -2;
      for (int j = 0; j < 3; j++) {
        CurrentNodeNucSequence[node][site_codon * 3 + j] = -2;
      }
    }
  }

  // reset nodeleaf sequences
  for (int taxa = 0; taxa < lparam->Ntaxa; taxa++) {
    for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
      CurrentLeafNodeCodonSequences[taxa][site_codon] = -2;
      for (int j = 0; j < 3; j++) {
        CurrentLeafNodeNucSequence[taxa][site_codon * 3 + j] = -2;
      }
    }
  }

  // reset ancestral sequences
  for (int point_i = 0; point_i < lparam->Ninterval; point_i++) {
    for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
      CurrentAncestralCodonSequence[point_i][0][site_codon] = -2;
    }
  }
}

void TreeSimulator::resetSimulatorFromLeaves() {
  int** cur_data = lparam->codondata->GetData();
  // reset nodeleaf sequences
  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
      int codon_state = cur_data[taxa_i][site_codon];
      if (codon_state != unknown) {
        CurrentLeafNodeCodonSequences[taxa_i][site_codon] = codon_state;
        for (int j = 0; j < 3; j++) {
          CurrentLeafNodeNucSequence[taxa_i][site_codon * 3 + j] =
              lparam->codonstatespace->GetCodonPosition(j, codon_state);
        }
      } else {
        double u = lparam->rnd->Uniform();
        std::vector<int> cur_vec;
        for (int taxa_j = 0; taxa_j < lparam->Ntaxa; taxa_j++) {
          codon_state = cur_data[taxa_j][site_codon];
          if (codon_state != unknown) {
            cur_vec.push_back(codon_state);
          }
        }
        int p = static_cast<int>(cur_vec.size() * u);
        codon_state = cur_vec[p];
        CurrentLeafNodeCodonSequences[taxa_i][site_codon] = codon_state;
        for (int j = 0; j < 3; j++) {
          CurrentLeafNodeNucSequence[taxa_i][site_codon * 3 + j] =
              lparam->codonstatespace->GetCodonPosition(j, codon_state);
        }
      }
    }
  }
  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    delete[] cur_data[taxa_i];
  }
  delete[] cur_data;
}

void TreeSimulator::GenerateFromLeaves() {
  submatrix->resetSubMatrixFromLeaves();
  resetSimulatorFromLeaves();
  // launch recursive simulation on a phylogenetic tree
  computeFromLeaves();
}

void TreeSimulator::GenerateCodonAlignment() {
  submatrix->resetSubMatrix();
  rootBranchEvoStats->resetEvoStats();
  treeEvoStats->resetEvoStats();
  resetSimulator();
  ancestralseq->ComputeStationaryCodon();

  if (this->lparam->rootlength == -1) {
    ancestralseq->SampleAncestralCodonSequenceFromCodonData();
  } else {
    ancestralseq->SampleAncestralCodonSequenceFromStationaryCodon();
  }

  SetAncestralSequence();
  // launch recursive simulation on a phylogenetic tree
  ComputeRecursiveSimulation(lparam->refTree->GetRoot());
  // register mappingstats

  resetEvoStatVectors();
  rootBranchEvoStats->GetEvoAncStats();
  treeEvoStats->GetEvoStats();
  treeEvoStats->GetSiteSpecificEvoStats();
}

void TreeSimulator::resetEvoStatVectors() {
  lparam->ancevostats.clear();
  lparam->ancevostats.shrink_to_fit();

  lparam->evostats.clear();
  lparam->evostats.shrink_to_fit();

  lparam->sitespecificevostats.clear();
  lparam->sitespecificevostats.shrink_to_fit();
}

void TreeSimulator::SetAncestralSequence() {
  int NodeIndex = lparam->refTree->GetRoot()->GetNode()->GetIndex();
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    this->CurrentNodeCodonSequence[NodeIndex][site_codon] =
        this->ancestralseq->GetCurrentAncestralCodonSequence(site_codon);
    for (int j = 0; j < 3; j++) {
      this->CurrentNodeNucSequence[NodeIndex][site_codon * 3 + j] =
          this->ancestralseq->GetCurrentAncestralNucSequence(site_codon * 3 +
                                                             j);
    }
  }
}

void TreeSimulator::RegisterSubTreeSim(int NodeIndex, int site_nuc, int nucTo) {
  int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
  int CodonPos = site_nuc % 3;
  // 3 codons length context

  int* nucposFrom = new int[9];
  int* nucposTo = new int[9];

  int* codonFrom = new int[3];
  int* codonTo = new int[3];

  int site_codon = static_cast<int>(site_nuc / 3);
  int site_codon_start = site_codon - 1;
  int site_codon_end = site_codon + 2;

  ///
  /// GET 3 ADJACENT CODONS FOR CODON INTERFACE STATS
  {
    int codon_count = 0;
    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end;
         site_codon_i++) {
      // int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
      for (int codonPos = 0; codonPos < 3; codonPos++) {
        if (site_codon_i >= 0 && site_codon_i < lparam->Nsite_codon) {
          nucposFrom[codon_count * 3 + codonPos] =
              CurrentNodeNucSequence[NodeIndex][site_codon_i * 3 + codonPos];
          // nucposTo[codon_count*3+codonPos] =
          // CurrentNodeNucSequence[NodeIndex][site_codon_i*3+codonPos];
          nucposTo[codon_count * 3 + codonPos] =
              nucposFrom[codon_count * 3 + codonPos];

        } else {
          nucposFrom[codon_count * 3 + codonPos] = -1;
          nucposTo[codon_count * 3 + codonPos] = -1;
        }
      }
      if (site_codon_i >= 0 && site_codon_i < lparam->Nsite_codon) {
        codonTo[codon_count] = lparam->codonstatespace->GetCodonFromDNA(
            nucposTo[codon_count * 3], nucposTo[codon_count * 3 + 1],
            nucposTo[codon_count * 3 + 2]);
        codonFrom[codon_count] = codonTo[codon_count];
      } else {
        codonTo[codon_count] = -1;
        codonFrom[codon_count] = -1;
      }
      codon_count++;
    }
  }

  int site_nuc_To = site_nuc % 3 + 3;
  nucposTo[site_nuc_To] = nucTo;
  codonTo[1] = lparam->codonstatespace->GetCodonFromDNA(
      nucposTo[3], nucposTo[4], nucposTo[5]);

  if (lparam->codonstatespace->CheckStop(nucposTo[3], nucposTo[4],
                                         nucposTo[5])) {
    std::cerr << "error while registring a substitution\n";
    std::cerr << nucTo << "\n";
    std::cerr << nucposTo[3] << " " << nucposTo[4] << " " << nucposTo[5] << " "
              << codonTo[1] << "\n";
    std::cerr << nucposFrom[3] << " " << nucposFrom[4] << " " << nucposFrom[5]
              << " " << codonFrom[1] << "\n";
    exit(0);
  }

  // Mapping statitics
  // Statistics on the ancestral sequence
  if (NodeIndex == lparam->refTree->GetRoot()->GetNode()->GetIndex()) {
    rootBranchEvoStats->Nsub++;
    // rootBranchEvoStats->ssNsub[site_codon]++;
    if (lparam->codonstatespace->Synonymous(codonFrom[1], codonTo[1])) {
      rootBranchEvoStats->Nsynsub++;
    }

    if (CodonPos == 0) {
      rootBranchEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
      rootBranchEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;
      rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                    [rootBranchEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To],
                                        nucposTo[site_nuc_To + 1])]++;  //
      rootBranchEvoStats
          ->dinuc_stat[CodonPos][rootBranchEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                      [rootBranchEvoStats->GetDinucContext(
                          nucposTo[site_nuc_To],
                          nucposTo[site_nuc_To + 1])]++;  //
      if (site_nuc > 0) {
        rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                      [rootBranchEvoStats->GetDinucContext(
                                          nucposTo[site_nuc_To - 1],
                                          nucposTo[site_nuc_To])]++;  //
        rootBranchEvoStats->dinuc_stat[2][rootBranchEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                      [rootBranchEvoStats->GetDinucContext(
                                          nucposTo[site_nuc_To - 1],
                                          nucposTo[site_nuc_To])]++;  //
      }

    } else if (CodonPos == 1) {
      rootBranchEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
      rootBranchEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;

      rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                    [rootBranchEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To],
                                        nucposTo[site_nuc_To + 1])]++;  //
      rootBranchEvoStats
          ->dinuc_stat[CodonPos][rootBranchEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                      [rootBranchEvoStats->GetDinucContext(
                          nucposTo[site_nuc_To],
                          nucposTo[site_nuc_To + 1])]++;  //

      rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [rootBranchEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;  //
      rootBranchEvoStats->dinuc_stat[0][rootBranchEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [rootBranchEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;  //

    } else if (CodonPos == 2) {
      rootBranchEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
      rootBranchEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;

      rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [rootBranchEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;  //
      rootBranchEvoStats->dinuc_stat[1][rootBranchEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [rootBranchEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;  //

      if (site_codon < lparam->Nsite_codon - 2) {
        rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                      [rootBranchEvoStats->GetDinucContext(
                                          nucposTo[site_nuc_To],
                                          nucposTo[site_nuc_To + 1])]++;  //
        rootBranchEvoStats
            ->dinuc_stat[CodonPos][rootBranchEvoStats->GetDinucContext(
                nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                        [rootBranchEvoStats->GetDinucContext(
                            nucposTo[site_nuc_To],
                            nucposTo[site_nuc_To + 1])]++;  //
      }
    }

  } else {
    // Non-root events
    treeEvoStats->Nsub++;
    treeEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
    treeEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;

    if (lparam->codonstatespace->Synonymous(codonFrom[1], codonTo[1])) {
      treeEvoStats->Nsynsub++;
      treeEvoStats->gtnrSyn_stat[3][nucFrom][nucTo]++;
      treeEvoStats->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;
    } else {
      treeEvoStats->gtnrNSyn_stat[3][nucFrom][nucTo]++;
      treeEvoStats->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;
    }
    if (CodonPos == 0) {
      treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                              [treeEvoStats->GetDinucContext(
                                  nucposTo[site_nuc_To],
                                  nucposTo[site_nuc_To + 1])]++;  //
      treeEvoStats->dinuc_stat[CodonPos][treeEvoStats->GetDinucContext(
          nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                              [treeEvoStats->GetDinucContext(
                                  nucposTo[site_nuc_To],
                                  nucposTo[site_nuc_To + 1])]++;  //

      // SECOND FRAM FOR DINUC
      if (site_nuc > 0) {
        treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                [treeEvoStats->GetDinucContext(
                                    nucposTo[site_nuc_To - 1],
                                    nucposTo[site_nuc_To])]++;
        treeEvoStats->dinuc_stat[2][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                [treeEvoStats->GetDinucContext(
                                    nucposTo[site_nuc_To - 1],
                                    nucposTo[site_nuc_To])]++;
      }

      if (lparam->codonstatespace->Synonymous(codonFrom[1], codonTo[1])) {
        treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To],
                                       nucposTo[site_nuc_To + 1])]++;
        treeEvoStats->dinucSyn_stat[CodonPos][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To],
                                       nucposTo[site_nuc_To + 1])]++;

        if (site_nuc > 0) {
          treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                     [treeEvoStats->GetDinucContext(
                                         nucposTo[site_nuc_To - 1],
                                         nucposTo[site_nuc_To])]++;
          treeEvoStats->dinucSyn_stat[2][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                     [treeEvoStats->GetDinucContext(
                                         nucposTo[site_nuc_To - 1],
                                         nucposTo[site_nuc_To])]++;
        }

      } else {
        treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To],
                                        nucposTo[site_nuc_To + 1])]++;
        treeEvoStats->dinucNSyn_stat[CodonPos][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To],
                                        nucposTo[site_nuc_To + 1])]++;

        if (site_nuc > 0) {
          treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                      [treeEvoStats->GetDinucContext(
                                          nucposTo[site_nuc_To - 1],
                                          nucposTo[site_nuc_To])]++;
          treeEvoStats->dinucNSyn_stat[2][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                      [treeEvoStats->GetDinucContext(
                                          nucposTo[site_nuc_To - 1],
                                          nucposTo[site_nuc_To])]++;
        }
      }

    } else if (CodonPos == 1) {
      treeEvoStats
          ->dinuc_stat[3][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                      [treeEvoStats->GetDinucContext(
                          nucposTo[site_nuc_To], nucposTo[site_nuc_To + 1])]++;
      treeEvoStats
          ->dinuc_stat[CodonPos][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                      [treeEvoStats->GetDinucContext(
                          nucposTo[site_nuc_To], nucposTo[site_nuc_To + 1])]++;
      treeEvoStats
          ->dinuc_stat[3][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                      [treeEvoStats->GetDinucContext(nucposTo[site_nuc_To - 1],
                                                     nucposTo[site_nuc_To])]++;
      treeEvoStats
          ->dinuc_stat[0][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                      [treeEvoStats->GetDinucContext(nucposTo[site_nuc_To - 1],
                                                     nucposTo[site_nuc_To])]++;
      if (lparam->codonstatespace->Synonymous(codonFrom[1], codonTo[1])) {
        treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To],
                                       nucposTo[site_nuc_To + 1])]++;
        treeEvoStats->dinucSyn_stat[CodonPos][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To],
                                       nucposTo[site_nuc_To + 1])]++;
        treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To - 1],
                                       nucposTo[site_nuc_To])]++;
        treeEvoStats->dinucSyn_stat[0][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To - 1],
                                       nucposTo[site_nuc_To])]++;
      } else {
        treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To],
                                        nucposTo[site_nuc_To + 1])]++;
        treeEvoStats->dinucNSyn_stat[CodonPos][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To],
                                        nucposTo[site_nuc_To + 1])]++;

        treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;
        treeEvoStats->dinucNSyn_stat[0][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;
      }
    } else if (CodonPos == 2) {
      treeEvoStats
          ->dinuc_stat[3][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                      [treeEvoStats->GetDinucContext(nucposTo[site_nuc_To - 1],
                                                     nucposTo[site_nuc_To])]++;
      treeEvoStats
          ->dinuc_stat[1][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                      [treeEvoStats->GetDinucContext(nucposTo[site_nuc_To - 1],
                                                     nucposTo[site_nuc_To])]++;

      if (site_codon < lparam->Nsite_codon - 2) {
        treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                [treeEvoStats->GetDinucContext(
                                    nucposTo[site_nuc_To],
                                    nucposTo[site_nuc_To + 1])]++;
        treeEvoStats->dinuc_stat[CodonPos][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                [treeEvoStats->GetDinucContext(
                                    nucposTo[site_nuc_To],
                                    nucposTo[site_nuc_To + 1])]++;
      }

      if (lparam->codonstatespace->Synonymous(codonFrom[1], codonTo[1])) {
        treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To - 1],
                                       nucposTo[site_nuc_To])]++;
        treeEvoStats->dinucSyn_stat[1][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                   [treeEvoStats->GetDinucContext(
                                       nucposTo[site_nuc_To - 1],
                                       nucposTo[site_nuc_To])]++;

        if (site_codon < lparam->Nsite_codon - 2) {
          treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                     [treeEvoStats->GetDinucContext(
                                         nucposTo[site_nuc_To],
                                         nucposTo[site_nuc_To + 1])]++;
          treeEvoStats->dinucSyn_stat[CodonPos][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                     [treeEvoStats->GetDinucContext(
                                         nucposTo[site_nuc_To],
                                         nucposTo[site_nuc_To + 1])]++;
        }

      } else {
        treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;
        treeEvoStats->dinucNSyn_stat[1][treeEvoStats->GetDinucContext(
            nucposFrom[site_nuc_To - 1], nucposFrom[site_nuc_To])]
                                    [treeEvoStats->GetDinucContext(
                                        nucposTo[site_nuc_To - 1],
                                        nucposTo[site_nuc_To])]++;

        if (site_codon < lparam->Nsite_codon - 2) {
          treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                      [treeEvoStats->GetDinucContext(
                                          nucposTo[site_nuc_To],
                                          nucposTo[site_nuc_To + 1])]++;
          treeEvoStats->dinucNSyn_stat[CodonPos][treeEvoStats->GetDinucContext(
              nucposFrom[site_nuc_To], nucposFrom[site_nuc_To + 1])]
                                      [treeEvoStats->GetDinucContext(
                                          nucposTo[site_nuc_To],
                                          nucposTo[site_nuc_To + 1])]++;
        }
      }
    }
  }

  // The evolving sequence is updated here

  CurrentNodeNucSequence[NodeIndex][site_nuc] = nucTo;
  CurrentNodeCodonSequence[NodeIndex][site_codon] = codonTo[1];

  delete[] nucposFrom;
  delete[] nucposTo;
  delete[] codonFrom;
  delete[] codonTo;
}

void TreeSimulator::ComputeRecursiveSimulation(Link* from) {
  int FromNodeIndex = from->GetNode()->GetIndex();

  if (from->isRoot()) {
    submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, -1,
                                      CurrentNodeNucSequence);

    double rate = 0.0;
    double time = 0.0;
    double blength = 0.0;

    rate = submatrix->GetTotalSubRate(FromNodeIndex);

    if (lparam->model == "FMutSelSimu") {
      rate = submatrix->GetTotalMutRate(FromNodeIndex);
    }

    //        rootBranchEvoStats->MutRate[0][0] =
    //        submatrix->GetTotalMutRate(FromNodeIndex);
    //        rootBranchEvoStats->MutRate[0][1] =
    //        submatrix->GetTotalMutRateNonSyn(FromNodeIndex);
    //        rootBranchEvoStats->MutRate[0][2] =
    //        submatrix->GetTotalMutRateSyn(FromNodeIndex);
    //        rootBranchEvoStats->SubRate[0][0] =
    //        submatrix->GetTotalSubRate(FromNodeIndex);
    //        rootBranchEvoStats->SubRate[0][1] =
    //        submatrix->GetTotalSubRateNonSyn(FromNodeIndex);
    //        rootBranchEvoStats->SubRate[0][2] =
    //        submatrix->GetTotalSubRateSyn(FromNodeIndex);

    time = (lparam->rnd->sExpo()) / rate;

    blength = lparam->rootlength;

    double IntervalLength = 1.0;
    int interval = 0;
    SetAncestralCodonSequence(FromNodeIndex, interval);
    interval++;
    while (time < blength) {
      double u =
          lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);
      int site_nuc = 0;
      int nucTo = 0;
      double testcummul = submatrix->GetSubRate(FromNodeIndex, site_nuc, nucTo);

      int k = 0;
      while (testcummul < u) {
        k++;
        site_nuc = static_cast<int>(k / 4);
        nucTo = k % 4;

        if (site_nuc >= lparam->Nsite_nuc) {
          std::cerr << "error " << site_nuc << " > " << lparam->Nsite_nuc
                    << "testcummul " << testcummul << "\n";
          exit(0);
        }

        testcummul += submatrix->GetSubRate(FromNodeIndex, site_nuc, nucTo);
      }

      int site_codon = static_cast<int>(site_nuc / 3);

      submatrix->ComputePartialRates(FromNodeIndex, site_codon,
                                     CurrentNodeNucSequence);
      RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo);

      submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, site_codon,
                                        CurrentNodeNucSequence);

      rate = submatrix->GetTotalSubRate(FromNodeIndex);
      if (lparam->model == "FMutSelSimu") {
        rate = submatrix->GetTotalMutRate(FromNodeIndex);
      }

      time += (lparam->rnd->sExpo()) / rate;

      if (time > IntervalLength * interval && interval < lparam->Ninterval) {
        SetAncestralCodonSequence(FromNodeIndex, interval);
        interval++;
      }
    }

    //        rootBranchEvoStats->MutRate[1][0] =
    //        submatrix->GetTotalMutRate(FromNodeIndex);
    //        rootBranchEvoStats->MutRate[1][1] =
    //        submatrix->GetTotalMutRateNonSyn(FromNodeIndex);
    //        rootBranchEvoStats->MutRate[1][2] =
    //        submatrix->GetTotalMutRateSyn(FromNodeIndex);
    //        rootBranchEvoStats->SubRate[1][0] =
    //        submatrix->GetTotalSubRate(FromNodeIndex);
    //        rootBranchEvoStats->SubRate[1][1] =
    //        submatrix->GetTotalSubRateNonSyn(FromNodeIndex);
    //        rootBranchEvoStats->SubRate[1][2] =
    //        submatrix->GetTotalSubRateSyn(FromNodeIndex);

    for (Link* link = from->Next(); link != from; link = link->Next()) {
      int OutNodeIndex = link->Out()->GetNode()->GetIndex();
      submatrix->transfertTotalRate(FromNodeIndex, OutNodeIndex);

      for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
        CurrentNodeCodonSequence[OutNodeIndex][site_codon] =
            CurrentNodeCodonSequence[FromNodeIndex][site_codon];
        for (int j = 0; j < 3; j++) {
          CurrentNodeNucSequence[OutNodeIndex][site_codon * 3 + j] =
              CurrentNodeNucSequence[FromNodeIndex][site_codon * 3 + j];
          submatrix->transfertNodeMatrix(FromNodeIndex, OutNodeIndex,
                                         site_codon * 3 + j);
        }
      }

      ComputeRecursiveSimulation(link->Out());
    }
  } else if (!from->isLeaf()) {
    double rate = 0.0;
    double time = 0.0;
    double blength = 0.0;

    rate = submatrix->GetTotalSubRate(FromNodeIndex);

    if (lparam->model == "FMutSelSimu") {
      rate = submatrix->GetTotalMutRate(FromNodeIndex);
    }

    time = (lparam->rnd->sExpo()) / rate;

    blength = atof(from->GetBranch()->GetName().c_str());

    while (time < blength) {
      double u =
          lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);
      int site_nuc = 0;
      int nucTo = 0;
      double testcummul = submatrix->GetSubRate(FromNodeIndex, site_nuc, nucTo);
      int k = 0;
      while (testcummul < u) {
        k++;
        site_nuc = static_cast<int>(k / 4);
        nucTo = k % 4;
        if (site_nuc >= lparam->Nsite_nuc) {
          std::cerr << "error " << site_nuc << " >= " << lparam->Nsite_nuc
                    << "\n";
          exit(0);
        }
        testcummul += submatrix->GetSubRate(FromNodeIndex, site_nuc, nucTo);
      }

      int site_codon = static_cast<int>(site_nuc / 3);

      submatrix->ComputePartialRates(FromNodeIndex, site_codon,
                                     CurrentNodeNucSequence);
      RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo);

      submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, site_codon,
                                        CurrentNodeNucSequence);
      rate = submatrix->GetTotalSubRate(FromNodeIndex);

      if (lparam->model == "FMutSelSimu") {
        rate = submatrix->GetTotalMutRate(FromNodeIndex);
      }

      time += (lparam->rnd->sExpo()) / rate;
    }

    for (Link* link = from->Next(); link != from; link = link->Next()) {
      int OutNodeIndex = link->Out()->GetNode()->GetIndex();
      submatrix->transfertTotalRate(FromNodeIndex, OutNodeIndex);
      for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
        CurrentNodeCodonSequence[OutNodeIndex][site_codon] =
            CurrentNodeCodonSequence[FromNodeIndex][site_codon];
        for (int j = 0; j < 3; j++) {
          CurrentNodeNucSequence[OutNodeIndex][site_codon * 3 + j] =
              CurrentNodeNucSequence[FromNodeIndex][site_codon * 3 + j];
          submatrix->transfertNodeMatrix(FromNodeIndex, OutNodeIndex,
                                         site_codon * 3 + j);
        }
      }

      ComputeRecursiveSimulation(link->Out());
    }
  } else if (from->isLeaf()) {
    double rate = 0.0;
    double time = 0.0;
    double blength = 0.0;
    submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, -1,
                                      CurrentNodeNucSequence);
    rate = submatrix->GetTotalSubRate(FromNodeIndex);
    if (lparam->model == "FMutSelSimu") {
      rate = submatrix->GetTotalMutRate(FromNodeIndex);
    }

    time = (lparam->rnd->sExpo()) / rate;

    blength = atof(from->GetBranch()->GetName().c_str());

    while (time < blength) {
      double u =
          lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);

      int site_nuc = 0;

      int nucTo = 0;
      double testcummul = submatrix->GetSubRate(FromNodeIndex, site_nuc, nucTo);
      int k = 0;
      while (testcummul < u) {
        k++;
        site_nuc = static_cast<int>(k / 4);
        nucTo = k % 4;
        if (site_nuc >= lparam->Nsite_nuc) {
          std::cerr << "error " << site_nuc << " >= " << lparam->Nsite_nuc
                    << "\n";
          exit(0);
        }
        testcummul += submatrix->GetSubRate(FromNodeIndex, site_nuc, nucTo);
      }
      int site_codon = static_cast<int>(site_nuc / 3);

      submatrix->ComputePartialRates(FromNodeIndex, site_codon,
                                     CurrentNodeNucSequence);
      RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo);
      submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, site_codon,
                                        CurrentNodeNucSequence);
      rate = submatrix->GetTotalSubRate(FromNodeIndex);

      if (lparam->model == "FMutSelSimu") {
        rate = submatrix->GetTotalMutRate(FromNodeIndex);
      }

      time += (lparam->rnd->sExpo()) / rate;
    }

    for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
      CurrentLeafNodeCodonSequences[FromNodeIndex][site_codon] =
          CurrentNodeCodonSequence[FromNodeIndex][site_codon];
      for (int j = 0; j < 3; j++) {
        CurrentLeafNodeNucSequence[FromNodeIndex][site_codon * 3 + j] =
            CurrentNodeNucSequence[FromNodeIndex][site_codon * 3 + j];
      }
    }
  }
}

void TreeSimulator::computeFromLeaves() {
  for (int taxa_i = 0; taxa_i < lparam->Ntaxa; taxa_i++) {
    submatrix->UpdateSubMatrixFromLeaves(taxa_i, CurrentLeafNodeCodonSequences);
  }
}

void TreeSimulator::SetAncestralCodonSequence(int FromNodeIndex, int interval) {
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    CurrentAncestralCodonSequence[interval][0][site_codon] =
        CurrentNodeCodonSequence[FromNodeIndex][site_codon];
  }
}

void TreeSimulator::WriteAncestral() {
  // Write stationary
  ofstream ancestral_os((lparam->output + ".ancestral").c_str(),
                        std::ios_base::app);
  for (int site_codon = 0; site_codon < lparam->Nsite_codon; site_codon++) {
    ancestral_os << lparam->codonstatespace->GetState(
        CurrentNodeCodonSequence[0][site_codon]);
  }
  ancestral_os << "\n";
  ancestral_os.close();
}
