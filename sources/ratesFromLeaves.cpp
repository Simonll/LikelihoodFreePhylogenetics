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

void writeHeaderFromLeaves(ofstream& os) {
  os << "chainID"
     << "\t"
     << "taxaID"
     << "\t"
     << "MutRate"
     << "\t"
     << "SubRate"
     << "\t"
     << "MutRateNonSyn"
     << "\t"
     << "SubRateNonSyn"
     << "\t"
     << "MutRateSyn"
     << "\t"
     << "SubRateSyn"
     << "\t"
     << "MutRateCpG"
     << "\t"
     << "SubRateCpG"
     << "\t"
     << "MutRateWeakStrong"
     << "\t"
     << "SubRateWeakStrong"
     << "\t"
     << "MutRateStrongWeak"
     << "\t"
     << "SubRateStrongWeak"
     << "\t"
     << "MutRateTransition"
     << "\t"
     << "SubRateTransition"
     << "\t"
     << "MutRateTransversion"
     << "\t"
     << "SubRateTransversion"
     << "\t"
     << "MutRateConsPol"
     << "\t"
     << "SubRateConsPol"
     << "\t"
     << "MutRateRadPol"
     << "\t"
     << "SubRateRadPol"
     << "\n";
}

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
  std::string controlfile = "";
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
      } else if (s == "-m") {
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
      } else if (s == "-conf") {
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
              << "-seqtype <stationary|data>\n"
              << "-conf <configfile>";
    exit(1);
  }
  if (model == "CodonMutSelFinite" || model == "CodonMutSelSBDP") {
    // the chain pointS are extract from the posterior file according to chainID
    std::cerr << model << "\n";

    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    // gparam->chainPointStart = start;
    // gparam->chainPointEvery = every;
    // gparam->chainPointEnd = until;
    // gparam->Nrun = Nrun;
    // gparam->Nrep = Nrep;
    // gparam->output = output;
    std::cerr << "global parameters registred"
              << "\n";
    LocalParameters* lparam = new LocalParameters(gparam);
    // lparam->taxa_a = taxa_a;
    // lparam->taxa_b = taxa_b;
    // lparam->posteriorfile = abc;
    // lparam->chain = mcmc;
    // lparam->code = code;
    // lparam->data = phylip;
    // lparam->rootlength = rootlength;
    std::cerr << "local parameters registred"
              << "\n";

    if (model == "CodonMutSelSBDP") {
      lparam->readChainCodonMutSelSBDP();
    } else if (model == "CodonMutSelFinite") {
      lparam->readChainCodonMutSelFinite();
    }

    SiteInterSubMatrixCABC2018* submatrix =
        new SiteInterSubMatrixCABC2018(lparam);
    submatrix->init();
    AncestralSequence* ancestraseq = new AncestralSequence(lparam);
    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);

    Posterior* post = new Posterior(gparam);
    post->readPosterior(lparam->posteriorfile);
    std::cerr << "The simulation process started\n";

    SummaryStatistics* ss = new SummaryStatistics(lparam);

    ofstream rates_os((gparam->output + ".rates").c_str(), std::ios_base::out);
    writeHeaderFromLeaves(rates_os);
    rates_os.close();

    ofstream ancestral_ss_os((gparam->output + ".anc").c_str(),
                             std::ios_base::out);
    ancestral_ss_os.close();
    bool writeHeaderAnc = true;

    if (!post->posterior.empty()) {
      int it = 0;
      while (it < gparam->Nsimu) {
        int pointID = static_cast<int>(
            lparam->rnd->Uniform() * post->posterior.size() - 1);
        lparam->SetCurrentParametersFromPosterior(post->posterior, pointID);
        if (model == "CodonMutSelSBDP") {
          lparam->readChainCodonMutSelSBDP(lparam->GetPointID());

        } else if (model == "CodonMutSelFinite") {
          lparam->readChainCodonMutSelFinite(lparam->GetPointID());
        }
        int NodeIndex = lparam->refTree->GetRoot()->GetNode()->GetIndex();

        double MutRate = 0.0;
        double SubRate = 0.0;
        double MutRateNonSyn = 0.0;
        double SubRateNonSyn = 0.0;
        double MutRateSyn = 0.0;
        double SubRateSyn = 0.0;
        double MutRateSyn = 0.0;
        double SubRateSyn = 0.0;
        double MutRateCpG = 0.0;
        double SubRateCpG = 0.0;
        double MutRateWeakStrong = 0.0;
        double SubRateWeakStrong = 0.0;
        double MutRateStrongWeak = 0.0;
        double SubRateStrongWeak = 0.0;
        double MutRateTransition = 0.0;
        double SubRateTransition = 0.0;
        double MutRateTransversion = 0.0;
        double SubRateTransversion = 0.0;
        double MutRateConsPol = 0.0;
        double SubRateConsPol = 0.0;
        double MutRateRadPol = 0.0;
        double SubRateRadPol = 0.0;
        int rep = 0;
        while (rep < gparam->Nrep) {
          simulator->run_jump_chain_over_seq(seqtype);
          double MutRate_ = 0.0;
          double SubRate_ = 0.0;
          std::tie(MutRate_, SubRate_) = submatrix->GetRates(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);

          MutRate += MutRate_;
          SubRate += SubRate_;

          double MutRateNonSyn_ = 0.0;
          double SubRateNonSyn_ = 0.0;
          std::tie(MutRateNonSyn_, SubRateNonSyn_) = submatrix->GetRatesNonSyn(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateNonSyn + = MutRateNonSyn_;
          SubRateNonSyn += SubRateNonSyn_;
          double MutRateSyn_ = 0.0;
          double SubRateSyn_ = 0.0;
          std::tie(MutRateSyn_, SubRateSyn_) = submatrix->GetRatesSyn(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateSyn += MutRateSyn_;
          SubRateSyn += SubRateSyn_;

          double MutRateCpG_ = 0.0;
          double SubRateCpG_ = 0.0;
          std::tie(MutRateCpG_, SubRateCpG_) = submatrix->GetRatesCpG(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateCpG += MutRateCpG_;
          SubRateCpG + = SubRateCpG_;
          double MutRateWeakStrong_ = 0.0;
          double SubRateWeakStrong_ = 0.0;
          std::tie(MutRateWeakStrong_, SubRateWeakStrong_) =
              submatrix->GetRatesWeakStrong(NodeIndex, -1,
                                            simulator->CurrentNodeNucSequence);
          MutRateWeakStrong + = MutRateWeakStrong_;
          SubRateWeakStrong += SubRateWeakStrong_;
          double MutRateStrongWeak_ = 0.0;
          double SubRateStrongWeak_ = 0.0;
          std::tie(MutRateStrongWeak_, SubRateStrongWeak_) =
              submatrix->GetRatesStrongWeak(NodeIndex, -1,
                                            simulator->CurrentNodeNucSequence);
          MutRateStrongWeak += MutRateStrongWeak_;
          SubRateStrongWeak += SubRateStrongWeak_;

          double MutRateTransition_ = 0.0;
          double SubRateTransition_ = 0.0;
          std::tie(MutRateTransition_, SubRateTransition_) =
              submatrix->GetRatesTransition(NodeIndex, -1,
                                            simulator->CurrentNodeNucSequence);
          MutRateTransition += MutRateTransition_;
          SubRateTransition += SubRateTransition_;

          double MutRateTransversion_ = 0.0;
          double SubRateTransversion_ = 0.0;
          std::tie(MutRateTransversion_, SubRateTransversion_) =
              submatrix->GetRatesTransversion(
                  NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateTransversion += MutRateTransversion;
          SubRateTransversion += SubRateTransversion;

          double MutRateConsPol_ = 0.0;
          double SubRateConsPol_ = 0.0;
          std::tie(MutRateConsPol_, SubRateConsPol_) =
              submatrix->GetRatesConsPol(NodeIndex, -1,
                                         simulator->CurrentNodeNucSequence);
          MutRateConsPol += MutRateConsPol_;
          SubRateConsPol += SubRateConsPol_;

          double MutRateRadPol_ = 0.0;
          double SubRateRadPol_ = 0.0;
          std::tie(MutRateRadPol_, SubRateRadPol_) = submatrix->GetRatesRadPol(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateRadPol += MutRateRadPol_;
          SubRateRadPol += SubRateRadPol_;
          int** ancestralCodonSequence_ = new int*[1];
          ancestralCodonSequence_[0] = new int[lparam->Nsite_codon];
          for (int site_codon = 0; site_codon < lparam->Nsite_codon;
               site_codon++) {
            ancestralCodonSequence_[0][site_codon] =
                simulator->CurrentNodeCodonSequence[NodeIndex][site_codon];
          }

          ss->computeSummariesAncestralSequence(ancestralCodonSequence_);
          ofstream ancestral_ss_os((gparam->output + ".anc").c_str(),
                                   std::ios_base::app);

          lparam->writeAncestralDataSummaries(ancestral_ss_os, writeHeaderAnc);
          ancestral_ss_os.close();
          writeHeaderAnc = false;

          it++;
          rep++;
          std::cerr << ".";
        }
        MutRate /= rep;
        SubRate /= rep;
        MutRateNonSyn /= rep;
        SubRateNonSyn /= rep;
        MutRateSyn /= rep;
        SubRateSyn /= rep;
        MutRateSyn /= rep;
        SubRateSyn /= rep;
        MutRateCpG /= rep;
        SubRateCpG /= rep;
        MutRateWeakStrong /= rep;
        SubRateWeakStrong /= rep;
        MutRateStrongWeak /= rep;
        SubRateStrongWeak /= rep;
        MutRateTransition /= rep;
        SubRateTransition /= rep;
        MutRateTransversion /= rep;
        SubRateTransversion /= rep;
        MutRateConsPol /= rep;
        SubRateConsPol /= rep;
        MutRateRadPol /= rep;
        SubRateRadPol /= rep;
        ofstream rates_os((gparam->output + ".rates").c_str(),
                          std::ios_base::app);
        rates_os << pointID << "\t"
                 << ((seqtype == "stationary") ? "NA"
                                               : lparam->taxonset->GetTaxon(
                                                     ancestraseq->choosen_taxa))
                 << "\t" << MutRate << "\t" << SubRate << "\t" << MutRateNonSyn
                 << "\t" << SubRateNonSyn << "\t" << MutRateSyn << "\t"
                 << SubRateSyn << "\t" << MutRateCpG << "\t" << SubRateCpG
                 << "\t" << MutRateWeakStrong << "\t" << SubRateWeakStrong
                 << "\t" << MutRateStrongWeak << "\t" << SubRateStrongWeak
                 << "\t" << MutRateTransition << "\t" << SubRateTransition
                 << "\t" << MutRateTransversion << "\t" << SubRateTransversion
                 << "\t" << MutRateConsPol << "\t" << SubRateConsPol << "\t"
                 << MutRateRadPol << "\t" << SubRateRadPol << "\n";

        rates_os.close();
      }
    }
  }
}
