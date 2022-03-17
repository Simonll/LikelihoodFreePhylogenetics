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
  int start = 0;
  int every = 0;
  int until = 0;
  std::string output = "";
  std::string code = "Universal";
  std::string taxa_a = "";
  std::string taxa_b = "";
  int rootlength = 100;
  std::string phylip = "";
  std::string mcmc = "";
  std::string abc = "";
  std::string seqtype = "";
  int Nrep = 0;
  int Nrun = 0;

  try {
    if (argc < 2) {
      throw(0);
    }
    int i = 1;
    while (i < argc) {
      std::string s = argv[i];
      if (s == "-v" || s == "--version") {
        throw(0);
      } else if (s == "-model") {
        i++;
        model = argv[i];
      } else if (s == "-x") {
        i++;
        start = atoi(argv[i]);
        i++;
        every = atoi(argv[i]);
        i++;
        until = atoi(argv[i]);
      } else if (s == "-output") {
        i++;
        output = argv[i];
      } else if (s == "-code") {
        i++;
        code = argv[i];
      } else if (s == "-root") {
        i++;
        taxa_a = argv[i];
        i++;
        taxa_b = argv[i];
      } else if (s == "-rep") {
        i++;
        Nrep = atoi(argv[i]);
      } else if (s == "-run") {
        i++;
        Nrun = atoi(argv[i]);
      } else if (s == "-mcmc") {
        i++;
        mcmc = argv[i];
      } else if (s == "-abc") {
        i++;
        abc = argv[i];
      } else if (s == "-d") {
        i++;
        phylip = argv[i];
      } else if (s == "-seqtype") {
        i++;
        seqtype = argv[i];
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
    std::cerr << "--version\n"
              << "-model <CodonMutSelFinite|CodonMutSelSBDP>\n"
              << "-x <start> <every> <until>\n"
              << "-output <>\n"
              << "-code <>\n"
              << "-root <taxa_a> <taxa_b>\n"
              << "-rep <int>\n"
              << "-run <int>\n"
              << "-mcmc <chain>\n"
              << "-abc <post>\n"
              << "-d <phylip>\n"
              << "-seqtype <stationary|data>\n";
    exit(1);
  }
  if (model == "CodonMutSelFinite" || model == "CodonMutSelSBDP") {
    // the chain pointS are extract from the posterior file according to chainID
    std::cerr << model << "\n";

    GlobalParameters* gparam = new GlobalParameters();
    gparam->model = model;
    gparam->chainPointStart = start;
    gparam->chainPointEvery = every;
    gparam->chainPointEnd = until;
    gparam->Nrun = Nrun;
    gparam->Nrep = Nrep;
    gparam->output = output;

    std::cerr << "global parameters registred"
              << "\n";
    LocalParameters* lparam = new LocalParameters(gparam);
    lparam->taxa_a = taxa_a;
    lparam->taxa_b = taxa_b;
    lparam->posteriorfile = abc;
    lparam->chain = mcmc;
    lparam->code = code;
    lparam->data = phylip;
    lparam->rootlength = rootlength;
    std::cerr << "local parameters registred"
              << "\n";

    if (model == "CodonMutSelSBDP") {
      lparam->readChainCodonMutSelSBDP();
    } else if (model == "CodonMutSelFinite") {
      lparam->readChainCodonMutSelFinite();
    }

    SiteInterSubMatrixCABC2018* submatrix =
        new SiteInterSubMatrixCABC2018(lparam);
    submatrix->initFromLeaves();
    std::cerr << lparam->Nsite_codon << "\n";
    TreeSimulator* simulator = new TreeSimulator(lparam, submatrix);

    Posterior* post = new Posterior(gparam);
    post->readPosterior(lparam->posteriorfile);
    std::cerr << "The simulation process started\n";

    ofstream rates_os((gparam->output + ".rates").c_str(), std::ios_base::out);
    submatrix->writeHeaderFromLeaves(rates_os);
    rates_os.close();
    if (!post->posterior.empty()) {
      int it = 0;
      while (it < gparam->Nrun) {
        int pointID = static_cast<int>(
            lparam->rnd->Uniform() * post->posterior.size() - 1);
        lparam->SetCurrentParametersFromPosterior(post->posterior, pointID);
        if (model == "CodonMutSelSBDP") {
          lparam->readChainCodonMutSelSBDP(lparam->GetPointID());

        } else if (model == "CodonMutSelFinite") {
          lparam->readChainCodonMutSelFinite(lparam->GetPointID());
        }

        int rep = 0;
        while (rep < gparam->Nrep) {
          simulator->GenerateFromLeaves();

          double u = lparam->rnd->Uniform();
          int taxaID = static_cast<int>((lparam->Ntaxa - 1) * u);
          double MutRate = 0.0;
          double SubRate = 0.0;
          std::tie(MutRate, SubRate) = submatrix->GetRates(
              taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          double MutRateNonSyn = 0.0;
          double SubRateNonSyn = 0.0;
          std::tie(MutRateNonSyn, SubRateNonSyn) = submatrix->GetRatesNonSyn(
              taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          double MutRateSyn = 0.0;
          double SubRateSyn = 0.0;
          std::tie(MutRateSyn, SubRateSyn) = submatrix->GetRatesSyn(
              taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          double MutRateCpG = 0.0;
          double SubRateCpG = 0.0;
          std::tie(MutRateCpG, SubRateCpG) = submatrix->GetRatesCpG(
              taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          double MutRateWeakStrong = 0.0;
          double SubRateWeakStrong = 0.0;
          std::tie(MutRateWeakStrong, SubRateWeakStrong) =
              submatrix->GetRatesWeakStrong(
                  taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          double MutRateStrongWeak = 0.0;
          double SubRateStrongWeak = 0.0;
          std::tie(MutRateStrongWeak, SubRateStrongWeak) =
              submatrix->GetRatesStrongWeak(
                  taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          double MutRateTransition = 0.0;
          double SubRateTransition = 0.0;
          std::tie(MutRateTransition, SubRateTransition) =
              submatrix->GetRatesTransition(
                  taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          double MutRateTransversion = 0.0;
          double SubRateTransversion = 0.0;
          std::tie(MutRateTransversion, SubRateTransversion) =
              submatrix->GetRatesTransversion(
                  taxaID, -1, simulator->CurrentLeafNodeNucSequence);

          ofstream rates_os((gparam->output + ".rates").c_str(),
                            std::ios_base::app);
          rates_os << pointID << "\t" << lparam->taxonset->GetTaxon(taxaID)
                   << "\t" << MutRate << "\t" << SubRate << "\t"
                   << MutRateNonSyn << "\t" << SubRateNonSyn << "\t"
                   << MutRateSyn << "\t" << SubRateSyn << "\t" << MutRateCpG
                   << "\t" << SubRateCpG << "\t" << MutRateWeakStrong << "\t"
                   << SubRateWeakStrong << "\t" << MutRateStrongWeak << "\t"
                   << SubRateStrongWeak << "\t" << MutRateTransition << "\t"
                   << SubRateTransition << "\t" << MutRateTransversion << "\t"
                   << SubRateTransversion << "\n";

          rates_os.close();
          rep++;
          std::cerr << ".";
        }
      }
      it++;
    }
  }
}
