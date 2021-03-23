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
#include "LocalData.h"

LocalData::LocalData(GlobalParameters* gparam) {
  // this->gparam = gparam;
  this->localcontrolfile = gparam->localcontrolfile;
  this->output = gparam->output;
  this->model = gparam->model;

  this->TOOSMALL = gparam->TOOSMALL;
  this->TOOLARGE = gparam->TOOLARGE;
  this->TOOLARGENEGATIVE = gparam->TOOLARGENEGATIVE;

  this->NSummaries = gparam->NSummaries;
  this->NParam = gparam->NParam;
  this->NMapStats = gparam->NEvoStats;

  this->listParam = new string[this->NParam];
  for (int param_i = 0; param_i < this->NParam; param_i++) {
    this->listParam[param_i] = gparam->listParam[param_i];
  }

  this->listSummaries = new string[this->NSummaries];
  for (int summary_i = 0; summary_i < this->NSummaries; summary_i++) {
    this->listSummaries[summary_i] = gparam->listSummaries[summary_i];
  }

  this->listMapStats = new string[this->NMapStats];
  for (int MapStats_i = 0; MapStats_i < this->NMapStats; MapStats_i++) {
    this->listMapStats[MapStats_i] = gparam->listEvoStats[MapStats_i];
  }

  this->NusedMapStats = gparam->NusedEvoStats;
  this->NusedMapAncStats = gparam->NusedEvoAncStats;
  this->NusedParam = gparam->NusedParam;
  this->NusedSummaries = gparam->NusedSummaries;
  this->Ngenes = gparam->Ngenes;

  this->mapUsedParam.insert(gparam->mapUsedParam.begin(),
                            gparam->mapUsedParam.end());
  this->mapUsedSummaries.insert(gparam->mapUsedSummaries.begin(),
                                gparam->mapUsedSummaries.end());
  this->mapUsedMapStats.insert(gparam->mapUsedEvoStats.begin(),
                               gparam->mapUsedEvoStats.end());
  this->mapUsedMapAncStats.insert(gparam->mapUsedEvoAncStats.begin(),
                                  gparam->mapUsedEvoAncStats.end());

  for (int gene_i = 0; gene_i < this->Ngenes; gene_i++) {
    this->listGenes.push_back(gparam->listGenes[gene_i]);
  }
}

LocalData::~LocalData() {
  // dtor
}

void LocalData::readLocalData(int k) {
  std::stringstream ss;
  ss << this->listGenes[0] << this->listGenes[k];
  this->data = ss.str();
  std::cerr << "data " << this->data << "\n";
  /*
          std::stringstream ss;
          ss << this->gparam->listTrees[0] << gparam->listTrees[k];
          this->tree = ss.str();
  */

  this->dnadata = new FileSequenceAlignment((data).c_str(), 0, 0);
  this->iscodon = true;
  this->code = "Universal";
  if (iscodon) {
    if (code == "Universal") {
      std::cerr << "Universal\n";
      this->codonstatespace = new CodonStateSpace(Universal);
      this->codondata = new CodonSequenceAlignment(dnadata, true, Universal);
      this->taxonset = this->codondata->GetTaxonSet();
    } else if (code == "MtMam") {
      std::cerr << "MtMam\n";
      this->codonstatespace = new CodonStateSpace(MtMam);
      this->codondata = new CodonSequenceAlignment(dnadata, true, MtMam);
      this->taxonset = this->codondata->GetTaxonSet();
    } else if (code == "MtInv") {
      std::cerr << "MtInv\n";
      this->codonstatespace = new CodonStateSpace(MtInv);
      this->codondata = new CodonSequenceAlignment(dnadata, true, MtInv);
      this->taxonset = this->codondata->GetTaxonSet();
    } else {
      std::cerr << "wrong genetic code\n";
    }
  }

  this->Nsite_codon = this->codondata->GetNsite();
  this->Ntaxa = this->codondata->GetNtaxa();
  this->Nstate_codon = this->codondata->GetNstate();
  this->Nsite_nuc = this->Nsite_codon * 3;
  this->Nnode = 0;  // Defined when the tree is red from chain file

  std::cerr << "Nsite_codon " << this->Nsite_codon << "\n";
  std::cerr << "Nsite_nuc " << this->Nsite_nuc << "\n";
  std::cerr << "Ntaxa " << this->Ntaxa << "\n";
  std::cerr << "Nstate_codon " << this->Nstate_codon << "\n";
  std::cerr << "Nstate_aa " << this->Nstate_aa << "\n";
}

void LocalData::toFasta(ofstream& os, int** curent_nodeleaf_sequence_codon) {
  for (int taxa = 0; taxa < Ntaxa; taxa++) {
    os << ">" << codondata->taxset->GetTaxon(taxa) << "\n";
    for (int site_codon = 0; site_codon < Nsite_codon; site_codon++) {
      os << codonstatespace->GetState(
          curent_nodeleaf_sequence_codon[taxa][site_codon]);
    }
    os << "\n";
  }
}

void LocalData::toAli(ofstream& os, int** curent_nodeleaf_sequence_codon) {
  os << Ntaxa << "\t" << (3 * Nsite_codon) << "\n";
  for (int taxa = 0; taxa < Ntaxa; taxa++) {
    os << codondata->taxset->GetTaxon(taxa) << ' ';
    for (int site_codon = 0; site_codon < Nsite_codon; site_codon++) {
      os << codonstatespace->GetState(
          curent_nodeleaf_sequence_codon[taxa][site_codon]);
    }
    os << "\n";
  }
  os << "\n";
}

void LocalData::writeRealDataSummaries(ofstream& os, bool headers) {
  string* arrSummaries = new string[NusedSummaries];
  for (int summary_i = 0; summary_i < NSummaries; summary_i++) {
    auto it = mapUsedSummaries.find(listSummaries[summary_i]);
    if (it != mapUsedSummaries.end()) {
      if (it->second != -1) {
        arrSummaries[it->second] = it->first;
      }
    }
  }

  if (headers) {
    for (int summary_i = 0; summary_i < NusedSummaries; summary_i++) {
      if (summary_i < NusedSummaries - 1) {
        os << arrSummaries[summary_i] << "\t";

      } else {
        os << arrSummaries[summary_i] << "\n";
      }
    }
  }

  for (int summary_i = 0; summary_i < NusedSummaries; summary_i++) {
    if (summary_i < NusedSummaries - 1) {
      os << summariesRealData[summary_i] << "\t";

    } else {
      os << summariesRealData[summary_i] << "\n";
    }
  }

  delete[] arrSummaries;
}
