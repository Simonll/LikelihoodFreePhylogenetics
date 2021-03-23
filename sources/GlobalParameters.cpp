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
#include "GlobalParameters.h"

GlobalParameters::GlobalParameters(string model, string controlfile) {
  this->model = model;
  this->controlfile = controlfile;
  this->Nrun = 1;
  this->Ncon = 1000;
  this->Nsimu = 1;
  this->Nrep = 1;
  this->Nthread = 1;
  this->Niter = 0;
  this->threshold = 10000;
  this->verbose = 0;

  this->chainPointStart = 1;
  this->chainPointEnd = 101;
  this->chainPointEvery = 1;

  this->NusedEvoStats = 0;
  this->NusedSiteSpecificEvoStats = 0;
  this->NusedEvoAncStats = 0;
  this->NusedParam = 0;
  this->NusedSummaries = 0;
  this->NusedAccessorySummaries = 0;
  this->Ngenes = 0;

  this->distance = "Euclidian";
  this->transformation = "log2";
  this->OutPartialDistance = 0;

  std::cerr << "Constructing mapUsedParam\n";
  for (unsigned int i = 0; i < this->NParam; i++) {
    mapUsedParam[listParam[i]] = -1;
  }

  std::cerr << "Constructing mapUsedAncSummaries\n";
  for (unsigned int i = 0; i < this->NSummaries; i++) {
    mapUsedAncSummaries[listSummaries[i]] = -1;
  }

  std::cerr << "Constructing mapUsedSummaries\n";
  for (unsigned int i = 0; i < this->NSummaries; i++) {
    mapUsedSummaries[listSummaries[i]] = -1;
  }

  std::cerr << "Constructing mapUsedAccessorySummaries\n";
  for (unsigned int i = 0; i < this->NSummaries; i++) {
    mapUsedAccessorySummaries[listSummaries[i]] = -1;
  }

  std::cerr << "Constructing mapUsedEvoStats\n";
  for (unsigned int i = 0; i < this->NEvoStats; i++) {
    mapUsedEvoStats[listEvoStats[i]] = -1;
  }

  std::cerr << "Constructing mapUsedSiteSpecificEvoStats\n";
  for (unsigned int i = 0; i < this->NSiteSpecificEvoStats; i++) {
    mapUsedSiteSpecificEvoStats[listSiteSpecificEvoStats[i]] = -1;
  }

  std::cerr << "Constructing mapUsedEvoAncStats\n";
  for (unsigned int i = 0; i < this->NEvoStats; i++) {
    mapUsedEvoAncStats[listEvoStats[i]] = -1;
  }

  std::cerr << "Reading instructions\n";
  readInstructions();
}
GlobalParameters::GlobalParameters() {}
GlobalParameters::~GlobalParameters() {
  // dtor
}

void GlobalParameters::readInstructions() {
  ifstream is(this->controlfile.c_str());
  if (!is) {
    std::cerr << "error: did not find " << this->controlfile << "\n";
    exit(1);
  }

  string line = "";
  while (std::getline(is, line)) {
    // std::cerr << line << "\n";
    if (!line.empty() && line.substr(0, 10) == "#SUMMARIES") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        auto it = mapUsedSummaries.find(w);
        if (it != mapUsedSummaries.end()) {
          it->second = k;
          k++;

        } else if (w != "#SUMMARIES") {
          std::cerr << "Undefined summary " << w << "\n";
          exit(0);
        }
      }
      NusedSummaries = k;
      std::cerr << "#SUMMARIES " << NusedSummaries << "\n";

    } else if (!line.empty() && line.substr(0, 13) == "#ACCSUMMARIES") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        auto it = mapUsedAccessorySummaries.find(w);
        if (it != mapUsedAccessorySummaries.end()) {
          it->second = k;
          k++;

        } else if (w != "#ACCSUMMARIES") {
          std::cerr << "Undefined summary " << w << "\n";
          exit(0);
        }
      }
      NusedAccessorySummaries = k;
      std::cerr << "#ACCSUMMARIES " << NusedAccessorySummaries << "\n";

    } else if (!line.empty() && line.substr(0, 13) == "#ANCSUMMARIES") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        auto it = mapUsedAncSummaries.find(w);
        if (it != mapUsedAncSummaries.end()) {
          it->second = k;
          k++;

        } else if (w != "#ANCSUMMARIES") {
          std::cerr << "Undefined ancestral summary " << w << "\n";
          exit(0);
        }
      }
      NusedAncSummaries = k;
      std::cerr << "#ANCSUMMARIES " << NusedAncSummaries << "\n";

    } else if (!line.empty() && line.substr(0, 6) == "#PARAM") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        auto it = mapUsedParam.find(w);
        if (it != mapUsedParam.end()) {
          it->second = k;
          std::cerr << "UsedParam " << k << " " << w << "\n";
          k++;

        } else if (w != "#PARAM") {
          std::cerr << "Undefined parameter " << w << "\n";
          exit(0);
        }
      }
      NusedParam = k;
      std::cerr << "#PARAM " << NusedParam << "\n";

    } else if (!line.empty() && line.substr(0, 6) == "#SSMAP") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        auto it = mapUsedSiteSpecificEvoStats.find(w);
        if (it != mapUsedSiteSpecificEvoStats.end()) {
          it->second = k;
          k++;

        } else if (w != "#SSMAP") {
          std::cerr << "Undefined SiteSpecificEvoStats " << w << "\n";
          exit(0);
        }
      }
      NusedSiteSpecificEvoStats = k;
      std::cerr << "#SSMAP " << NusedSiteSpecificEvoStats << "\n";

    } else if (!line.empty() && line.substr(0, 4) == "#MAP") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        auto it = mapUsedEvoStats.find(w);
        if (it != mapUsedEvoStats.end()) {
          it->second = k;
          k++;

        } else if (w != "#MAP") {
          std::cerr << "Undefined EvoStats parameter " << w << "\n";
          exit(0);
        }
      }
      NusedEvoStats = k;
      std::cerr << "#MAP " << NusedEvoStats << "\n";

    } else if (!line.empty() && line.substr(0, 13) == "#ANCESTRALMAP") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        auto it = mapUsedEvoAncStats.find(w);
        if (it != mapUsedEvoAncStats.end()) {
          it->second = k;
          k++;

        } else if (w != "#ANCESTRALMAP") {
          std::cerr << "Undefined EvoAncStats parameter" << w << "\n";
          exit(0);
        }
      }
      this->NusedEvoAncStats = k;
      std::cerr << "#ANCESTRALMAP " << NusedEvoAncStats << "\n";

    } else if (!line.empty() && line.substr(0, 6) == "#GENES") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        if (w != "#GENES") {
          k++;
          this->listGenes.push_back(w);
        }
      }
      this->Ngenes = k;
      std::cerr << "#GENES\t" << this->Ngenes << "\n";
      for (int i = 0; i < static_cast<int>(this->listGenes.size()); i++) {
        std::cerr << this->listGenes[i] << "\t";
      }
      std::cerr << "\n";

    } else if (!line.empty() && line.substr(0, 6) == "#TRANS") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        if (w != "#TRANS") {
          k++;
          this->transformation = w;
        }
      }

      std::cerr << "#TRANS\t" << this->transformation << "\n";

    } else if (!line.empty() && line.substr(0, 8) == "#OUTDIST") {
      this->OutPartialDistance = 1;

      std::cerr << "#OUTDIST\t" << this->distance << "\n";

    } else if (!line.empty() && line.substr(0, 11) == "#SPEUDODATA") {
      istringstream iss(line);
      string w;
      int k = 0;
      iss >> w;
      iss >> w;
      this->Ntaxa = atoi(w.c_str());
      iss >> w;
      this->Nsite_codon = atoi(w.c_str());
      while (iss >> w) {
        if (w != "#SPEUDODATA") {
          k++;
          this->listSpecies.push_back(w);
        }
      }
      if (this->Ntaxa != static_cast<int>(this->listSpecies.size())) {
        std::cerr << "Error the number of taxa does not match the list species "
                     "size\n";
        exit(0);
      }

      std::cerr << "#SPEUDODATA " << this->Ntaxa << " " << this->Nsite_codon
                << "\n";
      for (auto i : this->listSpecies) {
        std::cerr << i << " ";
      }
      std::cerr << "\n";

    } else if (!line.empty() && line.substr(0, 5) == "#DIST") {
      istringstream iss(line);
      string w;
      int k = 0;
      while (iss >> w) {
        if (w != "#DIST") {
          k++;
          this->distance = w;
        }
      }

      std::cerr << "#DIST\t" << this->distance << "\n";

    } else if (!line.empty() && line.substr(0, 7) == "#CHAINS") {
      istringstream iss(line);
      string w;
      while (iss >> w) {
        listChains.push_back(w);
        iss >> w;
        listPoints.push_back(atoi(w.c_str()));
      }

    } else if (!line.empty() && line.substr(0, 15) == "#HYPERMUT-DINUC") {
      std::cerr << "### HYPERMUT-DINUC ###\n";
      istringstream iss(line);
      string w;
      while (iss >> w) {
        std::cerr << w << "\n";
      }

    } else if (!line.empty() && line.substr(0, 5) == "#RUN") {
      istringstream iss(line);
      string w;
      iss >> w;
      iss >> w;
      this->Nrep = atoi(w.c_str());
      iss >> w;
      this->Nsimu = atoi(w.c_str());
      iss >> w;
      this->Ncon = atoi(w.c_str());
      std::cerr << "#Nrep " << this->Nrep << "#Nsimu" << this->Nsimu << "#Ncon"
                << this->Ncon << "\n";

    } else if (!line.empty() && line.substr(0, 5) == "#NRUN") {
      istringstream iss(line);
      string w;
      iss >> w;
      iss >> w;
      this->Nrun = atoi(w.c_str());
      iss >> w;
      this->threshold = atoi(w.c_str());
      std::cerr << "#Nrun " << this->Nrun << " threshold " << this->threshold
                << "\n";

    } else if (!line.empty() && line.substr(0, 9) == "#NTHREADS") {
      istringstream iss(line);
      string w;
      iss >> w;
      iss >> w;
      this->Nthread = atoi(w.c_str());
      std::cerr << "#Nthread " << this->Nthread << "\n";

    } else if (!line.empty() && line.substr(0, 7) == "#OUTPUT") {
      istringstream iss(line);
      string w;
      iss >> w;
      iss >> w;
      this->output = w;
      std::cerr << "#OUTPUT " << this->output << "\n";

    } else if (!line.empty() && line.substr(0, 8) == "#VERBOSE") {
      this->verbose = 1;
      std::cerr << "#VERBOSE " << this->verbose << "\n";

    } else if (!line.empty() && line.substr(0, 9) == "#SAMPLING") {
      istringstream iss(line);
      string w;
      iss >> w;
      iss >> w;
      this->chainPointStart = atoi(w.c_str());
      iss >> w;
      this->chainPointEvery = atoi(w.c_str());
      iss >> w;
      this->chainPointEnd = atoi(w.c_str());
      std::cerr << "#SAMPLING " << this->chainPointStart << " "
                << this->chainPointEvery << " " << this->chainPointEnd << "\n";

    } else if (!line.empty() && line.substr(0, 10) == "#OLDPARAMS") {
      std::cerr << "### OLDPARAMS ### \n";
      localcontrolfile = line;

    } else if (!line.empty() && line.substr(0, 11) == "#LOCALPARAM") {
      // localcontrolfile should be a list so different localcontrolfile could
      // be tested
      std::cerr << "### LOCALPARAM ### \n";
      localcontrolfile = line;
    }
  }
  is.close();
}
