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

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "AncestralSequence.h"
#include "BiologicalSequences.h"
#include "GlobalParameters.h"
#include "LocalParameters.h"
#include "Posterior.h"
#include "PriorSampler.h"
#include "SiteInterSubMatrixCABC2018.h"
#include "SummaryStatistics.h"
#include "TreeSimulator.h"

int main(int argc, char* argv[]) {
  // Comments

  // program options
  std::string model = "";
  std::string controlfile = "";

  try {
    if (argc < 2) {
      throw(0);
    }
    int i = 1;
    while (i < argc) {
      std::string s = argv[i];
      if (s == "-v" || s == "--version") {
        throw(0);
      } else if (s == "-m") {
        i++;
        model = argv[i];
        i++;
        controlfile = argv[i];
      }
      i++;
    }  // end while
  }    // end try
  catch (...) {
    std::cerr << "\n";
    std::cerr << "version 1.0\n";
    std::cerr << "###########################\n";
    std::cerr << "-m < show | CodonMutSelFinite | CodonMutSelSBDP \n";
    std::cerr << "###########################\n";
    std::cerr << "#SUMMARIES\n";
    std::cerr << "#ANCSUMMARIES\n";
    std::cerr << "#ACCSUMMARIES\n";
    std::cerr << "#PARAM\n";
    std::cerr << "#SSMAP\n";
    std::cerr << "#MAP\n";
    std::cerr << "#ANCESTRALMAP\n";
    std::cerr << "#DIST\n";
    std::cerr << "#TRANS\n";
    std::cerr << "#SPEUDODATA\n";
    std::cerr << "#NRUN\n";
    std::cerr << "#SAMPLING\n";
    std::cerr << "#LOCALPARAM\n";
    std::cerr << "###########################\n";
    exit(1);
  }
  if (model == "CodonMutSelFinite" || model == "CodonMutSelSBDP") {
    // the chain pointS are extract from the posterior file according to chainID
    std::cerr << model << "\n";

    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    Posterior* post = new Posterior(gparam);
    LocalParameters* lparam = new LocalParameters(gparam);

    if (model == "CodonMutSelSBDPSeq") {
      lparam->readChainCodonMutSelSBDP();
    } else if (model == "CodonMutSelFinite") {
      lparam->readChainCodonMutSelFinite();
    }
    std::string s =
        "compute MutRate and SubRate from leaves given posterior values";
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam, s);
    std::cerr << lparam->Nsite_codon << "\n";
    TreeSimulator* simulator = new TreeSimulator(lparam, submatrix);
    post->readPosterior(lparam->posteriorfile);
    std::cerr << "The simulation process started\n";

    ofstream rates_os((gparam->output + ".rates").c_str(), std::ios_base::out);
    rates_os << "chainID"
             << "\n"
             << "taxaID"
             << "\t"
             << "MutRate"
             << "\t"
             << "SubRate"
             << "\n";

    if (!post->posterior.empty()) {
      int it = 0;
      while (it < post->threshold) {
        int chainID =
            static_cast<int>(lparam->rnd->Uniform() * post->posterior.size());
        lparam->SetCurrentParametersFromPosterior(post->posterior, chainID);

        if (model == "CodonMutSelSBDPPPred") {
          lparam->readChainCodonMutSelSBDP(lparam->GetPointID());

        } else if (model == "CodonMutSelFinitePPred") {
          lparam->readChainCodonMutSelFinite(lparam->GetPointID());
        }

        int rep = 0;
        while (rep < post->Nrun) {
          simulator->GenerateFromLeaves();
          double u = lparam->rnd->Uniform();
          int taxaID = static_cast<int>(lparam->Ntaxa * u);
          double MutRate = 0.0;
          double SubRate = 0.0;
          for (int site_codon = 0; site_codon < lparam->Nsite_codon;
               site_codon++) {
            MutRate += submatrix->GetMutRate(taxaID, site_codon);
            SubRate += submatrix->GetSubRate(taxaID, site_codon);
          }
          ofstream rates_os((gparam->output + ".rates").c_str(),
                            std::ios_base::app);
          rates_os << chainID << "\t" << lparam->taxonset->GetTaxon(taxaID)
                   << "\t" << MutRate << "\t" << SubRate << "\n";
          rep++;
          std::cerr << ".";
        }
        it++;
      }
    }
  }
}
