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
     << "MutRateNonSynTs"
     << "\t"
     << "SubRateNonSynTs"
     << "\t"
     << "MutRateNonSynTr"
     << "\t"
     << "SubRateNonSynTr"
     << "\t"
     << "MutRateSyn"
     << "\t"
     << "SubRateSyn"
     << "\t"
     << "MutRateCpGTs"
     << "\t"
     << "SubRateCpGTs"
     << "\t"
     << "MutRateNonSynCpGTs"
     << "\t"
     << "SubRateNonSynCpGTs"
     << "\t"
     << "MutRateSynCpGTs"
     << "\t"
     << "SubRateSynCpGTs"
     << "\t"
     << "MutRateWeakStrong"
     << "\t"
     << "SubRateWeakStrong"
     << "\t"
     << "MutRateStrongWeak"
     << "\t"
     << "SubRateStrongWeak"
     << "\t"
     << "MutRateStrongStrong"
     << "\t"
     << "SubRateStrongStrong"
     << "\t"
     << "MutRateWeakWeak"
     << "\t"
     << "SubRateWeakWeak"
     << "\t"
     << "MutRateTs"
     << "\t"
     << "SubRateTs"
     << "\t"
     << "MutRateTr"
     << "\t"
     << "SubRateTr"
     << "\t"
     << "MutRateConsPol"
     << "\t"
     << "SubRateConsPol"
     << "\t"
     << "MutRateRadPol"
     << "\t"
     << "SubRateRadPol"
     << "\t"
     << "MutRateConsVol"
     << "\t"
     << "SubRateConsVol"
     << "\t"
     << "MutRateRadVol"
     << "\t"
     << "SubRateRadVol"
     << "\t"
     << "MutRateConsPolTs"
     << "\t"
     << "SubRateConsPolTs"
     << "\t"
     << "MutRateRadPolTs"
     << "\t"
     << "SubRateRadPolTs"
     << "\t"
     << "MutRateConsVolTs"
     << "\t"
     << "SubRateConsVolTs"
     << "\t"
     << "MutRateRadVolTs"
     << "\t"
     << "SubRateRadVolTs"
     << "\t"
     << "MutRateConsPolTr"
     << "\t"
     << "SubRateConsPolTr"
     << "\t"
     << "MutRateRadPolTr"
     << "\t"
     << "SubRateRadPolTr"
     << "\t"
     << "MutRateConsVolTr"
     << "\t"
     << "SubRateConsVolTr"
     << "\t"
     << "MutRateRadVolTr"
     << "\t"
     << "SubRateRadVolTr"
     << "\t"
     << "MutRateConsPolCpGTs"
     << "\t"
     << "SubRateConsPolCpGTs"
     << "\t"
     << "MutRateRadPolCpGTs"
     << "\t"
     << "SubRateRadPolCpGTs"
     << "\t"
     << "MutRateConsVolCpGTs"
     << "\t"
     << "SubRateConsVolCpGTs"
     << "\t"
     << "MutRateRadVolCpGTs"
     << "\t"
     << "SubRateRadVolCpGTs"
     << "\t"
     << "MutRateCGNTs"
     << "\t"
     << "SubRateCGNTs"
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

        double MutRateNonSynTs = 0.0;
        double SubRateNonSynTs = 0.0;
        double MutRateNonSynTr = 0.0;
        double SubRateNonSynTr = 0.0;

        double MutRateCpGTs = 0.0;
        double SubRateCpGTs = 0.0;
        double MutRateNonSynCpGTs = 0.0;
        double SubRateNonSynCpGTs = 0.0;
        double MutRateSynCpGTs = 0.0;
        double SubRateSynCpGTs = 0.0;

        double MutRateWeakStrong = 0.0;
        double SubRateWeakStrong = 0.0;
        double MutRateStrongWeak = 0.0;
        double SubRateStrongWeak = 0.0;
        double MutRateWeakWeak = 0.0;
        double SubRateWeakWeak = 0.0;
        double MutRateStrongStrong = 0.0;
        double SubRateStrongStrong = 0.0;

        double MutRateTs = 0.0;
        double SubRateTs = 0.0;
        double MutRateTr = 0.0;
        double SubRateTr = 0.0;

        double MutRateConsPol = 0.0;
        double SubRateConsPol = 0.0;
        double MutRateRadPol = 0.0;
        double SubRateRadPol = 0.0;

        double MutRateConsVol = 0.0;
        double SubRateConsVol = 0.0;
        double MutRateRadVol = 0.0;
        double SubRateRadVol = 0.0;

        double MutRateConsPolTs = 0.0;
        double SubRateConsPolTs = 0.0;
        double MutRateRadPolTs = 0.0;
        double SubRateRadPolTs = 0.0;

        double MutRateConsVolTs = 0.0;
        double SubRateConsVolTs = 0.0;
        double MutRateRadVolTs = 0.0;
        double SubRateRadVolTs = 0.0;

        double MutRateConsPolTr = 0.0;
        double SubRateConsPolTr = 0.0;
        double MutRateRadPolTr = 0.0;
        double SubRateRadPolTr = 0.0;

        double MutRateConsVolTr = 0.0;
        double SubRateConsVolTr = 0.0;
        double MutRateRadVolTr = 0.0;
        double SubRateRadVolTr = 0.0;

        double MutRateConsPolCpGTs = 0.0;
        double SubRateConsPolCpGTs = 0.0;
        double MutRateRadPolCpGTs = 0.0;
        double SubRateRadPolCpGTs = 0.0;

        double MutRateConsVolCpGTs = 0.0;
        double SubRateConsVolCpGTs = 0.0;
        double MutRateRadVolCpGTs = 0.0;
        double SubRateRadVolCpGTs = 0.0;

        double MutRateCGNTs = 0.0;
        double SubRateCGNTs = 0.0;

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
          MutRateNonSyn += MutRateNonSyn_;
          SubRateNonSyn += SubRateNonSyn_;

          double MutRateNonSynTs_ = 0.0;
          double SubRateNonSynTs_ = 0.0;
          std::tie(MutRateNonSynTs_, SubRateNonSynTs_) =
              submatrix->GetRatesNonSynTs(NodeIndex, -1,
                                          simulator->CurrentNodeNucSequence);
          MutRateNonSynTs += MutRateNonSynTs_;
          SubRateNonSynTs += SubRateNonSynTs_;

          double MutRateNonSynTr_ = 0.0;
          double SubRateNonSynTr_ = 0.0;
          std::tie(MutRateNonSyn_, SubRateNonSyn_) =
              submatrix->GetRatesNonSynTr(NodeIndex, -1,
                                          simulator->CurrentNodeNucSequence);
          MutRateNonSynTr += MutRateNonSynTr_;
          SubRateNonSynTr += SubRateNonSynTr_;

          double MutRateSyn_ = 0.0;
          double SubRateSyn_ = 0.0;
          std::tie(MutRateSyn_, SubRateSyn_) = submatrix->GetRatesSyn(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateSyn += MutRateSyn_;
          SubRateSyn += SubRateSyn_;

          double MutRateCpGTs_ = 0.0;
          double SubRateCpGTs_ = 0.0;
          std::tie(MutRateCpGTs_, SubRateCpGTs_) = submatrix->GetRatesCpGTs(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateCpGTs += MutRateCpGTs_;
          SubRateCpGTs += SubRateCpGTs_;

          double MutRateSynCpGTs_ = 0.0;
          double SubRateSynCpGTs_ = 0.0;
          std::tie(MutRateSynCpGTs_, SubRateSynCpGTs_) =
              submatrix->GetRatesSynCpGTs(NodeIndex, -1,
                                          simulator->CurrentNodeNucSequence);
          MutRateSynCpGTs += MutRateSynCpGTs_;
          SubRateSynCpGTs += SubRateSynCpGTs_;

          double MutRateNonSynCpGTs_ = 0.0;
          double SubRateNonSynCpGTs_ = 0.0;
          std::tie(MutRateNonSynCpGTs_, SubRateNonSynCpGTs_) =
              submatrix->GetRatesNonSynCpGTs(NodeIndex, -1,
                                             simulator->CurrentNodeNucSequence);
          MutRateNonSynCpGTs += MutRateNonSynCpGTs_;
          SubRateNonSynCpGTs += SubRateNonSynCpGTs_;

          double MutRateWeakStrong_ = 0.0;
          double SubRateWeakStrong_ = 0.0;
          std::tie(MutRateWeakStrong_, SubRateWeakStrong_) =
              submatrix->GetRatesWeakStrong(NodeIndex, -1,
                                            simulator->CurrentNodeNucSequence);
          MutRateWeakStrong += MutRateWeakStrong_;
          SubRateWeakStrong += SubRateWeakStrong_;

          double MutRateStrongWeak_ = 0.0;
          double SubRateStrongWeak_ = 0.0;
          std::tie(MutRateStrongWeak_, SubRateStrongWeak_) =
              submatrix->GetRatesStrongWeak(NodeIndex, -1,
                                            simulator->CurrentNodeNucSequence);
          MutRateStrongWeak += MutRateStrongWeak_;
          SubRateStrongWeak += SubRateStrongWeak_;

          double MutRateTs_ = 0.0;
          double SubRateTs_ = 0.0;
          std::tie(MutRateTs_, SubRateTs_) = submatrix->GetRatesTransition(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateTs += MutRateTs_;
          SubRateTs += SubRateTs_;

          double MutRateTr_ = 0.0;
          double SubRateTr_ = 0.0;
          std::tie(MutRateTr_, SubRateTr_) = submatrix->GetRatesTransversion(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateTr += MutRateTr_;
          SubRateTr += SubRateTr_;

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

          double MutRateConsVol_ = 0.0;
          double SubRateConsVol_ = 0.0;
          std::tie(MutRateConsVol_, SubRateConsVol_) =
              submatrix->GetRatesConsVol(NodeIndex, -1,
                                         simulator->CurrentNodeNucSequence);
          MutRateConsVol += MutRateConsVol_;
          SubRateConsVol += SubRateConsVol_;

          double MutRateRadVol_ = 0.0;
          double SubRateRadVol_ = 0.0;
          std::tie(MutRateRadVol_, SubRateRadVol_) = submatrix->GetRatesRadVol(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateRadVol += MutRateRadVol_;
          SubRateRadVol += SubRateRadVol_;

          double MutRateConsPolTs_ = 0.0;
          double SubRateConsPolTs_ = 0.0;
          std::tie(MutRateConsPolTs_, SubRateConsPolTs_) =
              submatrix->GetRatesConsPolTs(NodeIndex, -1,
                                           simulator->CurrentNodeNucSequence);
          MutRateConsPolTs += MutRateConsPolTs_;
          SubRateConsPolTs += SubRateConsPolTs_;

          double MutRateRadPolTs_ = 0.0;
          double SubRateRadPolTs_ = 0.0;
          std::tie(MutRateRadPolTs_, SubRateRadPolTs_) =
              submatrix->GetRatesRadPolTs(NodeIndex, -1,
                                          simulator->CurrentNodeNucSequence);
          MutRateRadPolTs += MutRateRadPolTs_;
          SubRateRadPolTs += SubRateRadPolTs_;

          double MutRateConsVolTs_ = 0.0;
          double SubRateConsVolTs_ = 0.0;
          std::tie(MutRateConsVolTs_, SubRateConsVolTs_) =
              submatrix->GetRatesConsVolTs(NodeIndex, -1,
                                           simulator->CurrentNodeNucSequence);
          MutRateConsVolTs += MutRateConsVolTs_;
          SubRateConsVolTs += SubRateConsVolTs_;

          double MutRateRadVolTs_ = 0.0;
          double SubRateRadVolTs_ = 0.0;
          std::tie(MutRateRadVolTs_, SubRateRadVolTs_) =
              submatrix->GetRatesRadVolTs(NodeIndex, -1,
                                          simulator->CurrentNodeNucSequence);
          MutRateRadVolTs += MutRateRadVolTs_;
          SubRateRadVolTs += SubRateRadVolTs_;

          double MutRateConsPolTr_ = 0.0;
          double SubRateConsPolTr_ = 0.0;
          std::tie(MutRateConsPolTr_, SubRateConsPolTr_) =
              submatrix->GetRatesConsPolTr(NodeIndex, -1,
                                           simulator->CurrentNodeNucSequence);
          MutRateConsPolTr += MutRateConsPolTr_;
          SubRateConsPolTr += SubRateConsPolTr_;

          double MutRateRadPolTr_ = 0.0;
          double SubRateRadPolTr_ = 0.0;
          std::tie(MutRateRadPolTr_, SubRateRadPolTr_) =
              submatrix->GetRatesRadPolTr(NodeIndex, -1,
                                          simulator->CurrentNodeNucSequence);
          MutRateRadPolTr += MutRateRadPolTr_;
          SubRateRadPolTr += SubRateRadPolTr_;

          double MutRateConsVolTr_ = 0.0;
          double SubRateConsVolTr_ = 0.0;
          std::tie(MutRateConsVolTr_, SubRateConsVolTr_) =
              submatrix->GetRatesConsVolTr(NodeIndex, -1,
                                           simulator->CurrentNodeNucSequence);
          MutRateConsVolTr += MutRateConsVolTr_;
          SubRateConsVolTr += SubRateConsVolTr_;

          double MutRateRadVolTr_ = 0.0;
          double SubRateRadVolTr_ = 0.0;
          std::tie(MutRateRadVolTr_, SubRateRadVolTr_) =
              submatrix->GetRatesRadVolTr(NodeIndex, -1,
                                          simulator->CurrentNodeNucSequence);
          MutRateRadVolTr += MutRateRadVolTr_;
          SubRateRadVolTr += SubRateRadVolTr_;

          double MutRateConsPolCpGTs_ = 0.0;
          double SubRateConsPolCpGTs_ = 0.0;
          std::tie(MutRateConsPolCpGTs_, SubRateConsPolCpGTs_) =
              submatrix->GetRatesConsPolCpGTs(
                  NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateConsPolCpGTs += MutRateConsPolCpGTs_;
          SubRateConsPolCpGTs += SubRateConsPolCpGTs_;

          double MutRateRadPolCpGTs_ = 0.0;
          double SubRateRadPolCpGTs_ = 0.0;
          std::tie(MutRateRadPolCpGTs_, SubRateRadPolCpGTs_) =
              submatrix->GetRatesRadPolCpGTs(NodeIndex, -1,
                                             simulator->CurrentNodeNucSequence);
          MutRateRadPolCpGTs += MutRateRadPolCpGTs_;
          SubRateRadPolCpGTs += SubRateRadPolCpGTs_;

          double MutRateConsVolCpGTs_ = 0.0;
          double SubRateConsVolCpGTs_ = 0.0;
          std::tie(MutRateConsVolCpGTs_, SubRateConsVolCpGTs_) =
              submatrix->GetRatesConsVolCpGTs(
                  NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateConsVolCpGTs += MutRateConsVolCpGTs_;
          SubRateConsVolCpGTs += SubRateConsVolCpGTs_;

          double MutRateCGNTs_ = 0.0;
          double SubRateCGNTs_ = 0.0;
          std::tie(MutRateCGNTs_, SubRateCGNTs_) = submatrix->GetRatesCGNTs(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateCGNTs += MutRateCGNTs_;
          SubRateCGNTs += SubRateCGNTs_;

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

        MutRateNonSynTs /= rep;
        SubRateNonSynTs /= rep;
        MutRateNonSynTr /= rep;
        SubRateNonSynTr /= rep;

        MutRateCpGTs /= rep;
        SubRateCpGTs /= rep;
        MutRateNonSynCpGTs /= rep;
        SubRateNonSynCpGTs /= rep;
        MutRateSynCpGTs /= rep;
        SubRateSynCpGTs /= rep;

        MutRateWeakStrong /= rep;
        SubRateWeakStrong /= rep;
        MutRateStrongWeak /= rep;
        SubRateStrongWeak /= rep;
        MutRateWeakWeak /= rep;
        SubRateWeakWeak /= rep;
        MutRateStrongStrong /= rep;
        SubRateStrongStrong /= rep;

        MutRateTs /= rep;
        SubRateTs /= rep;
        MutRateTr /= rep;
        SubRateTr /= rep;

        MutRateConsPol /= rep;
        SubRateConsPol /= rep;
        MutRateRadPol /= rep;
        SubRateRadPol /= rep;

        MutRateConsVol /= rep;
        SubRateConsVol /= rep;
        MutRateRadVol /= rep;
        SubRateRadVol /= rep;

        MutRateConsPolTs /= rep;
        SubRateConsPolTs /= rep;
        MutRateRadPolTs /= rep;
        SubRateRadPolTs /= rep;

        MutRateConsVolTs /= rep;
        SubRateConsVolTs /= rep;
        MutRateRadVolTs /= rep;
        SubRateRadVolTs /= rep;

        MutRateConsPolTr /= rep;
        SubRateConsPolTr /= rep;
        MutRateRadPolTr /= rep;
        SubRateRadPolTr /= rep;

        MutRateConsVolTr /= rep;
        SubRateConsVolTr /= rep;
        MutRateRadVolTr /= rep;
        SubRateRadVolTr /= rep;

        MutRateConsPolCpGTs /= rep;
        SubRateConsPolCpGTs /= rep;
        MutRateRadPolCpGTs /= rep;
        SubRateRadPolCpGTs /= rep;

        MutRateConsVolCpGTs /= rep;
        SubRateConsVolCpGTs /= rep;
        MutRateRadVolCpGTs /= rep;
        SubRateRadVolCpGTs /= rep;

        MutRateCGNTs /= rep;
        SubRateCGNTs /= rep;
        ofstream rates_os((gparam->output + ".rates").c_str(),
                          std::ios_base::app);
        rates_os << pointID << "\t"
                 << ((seqtype == "stationary") ? "NA"
                                               : lparam->taxonset->GetTaxon(
                                                     ancestraseq->choosen_taxa))
                 << "\t" << MutRate << "\t" << SubRate << "\t" << MutRateNonSyn
                 << "\t" << SubRateNonSyn << "\t" << MutRateNonSynTs << "\t"
                 << SubRateNonSynTs << "\t" << MutRateNonSynTr << "\t"
                 << SubRateNonSynTr << "\t" << MutRateSyn << "\t" << SubRateSyn
                 << "\t" << MutRateCpGTs << "\t" << SubRateCpGTs << "\t"
                 << MutRateNonSynCpGTs << "\t" << SubRateNonSynCpGTs << "\t"
                 << MutRateSynCpGTs << "\t" << SubRateSynCpGTs << "\t"
                 << MutRateWeakStrong << "\t" << SubRateWeakStrong << "\t"
                 << MutRateStrongWeak << "\t" << SubRateStrongWeak << "\t"
                 << MutRateStrongStrong << "\t" << SubRateStrongStrong << "\t"
                 << MutRateWeakWeak << "\t" << SubRateWeakWeak << "\t"
                 << MutRateTs << "\t" << SubRateTs << "\t" << MutRateTr << "\t"
                 << SubRateTr << "\t" << MutRateConsPol << "\t"
                 << SubRateConsPol << "\t" << MutRateRadPol << "\t"
                 << SubRateRadPol << "\t" << MutRateConsVol << "\t"
                 << SubRateConsVol << "\t" << MutRateRadVol << "\t"
                 << SubRateRadVol << "\t" << MutRateConsPolTs << "\t"
                 << SubRateConsPolTs << "\t" << MutRateRadPolTs << "\t"
                 << SubRateRadPolTs << "\t" << MutRateConsVolTs << "\t"
                 << SubRateConsVolTs << "\t" << MutRateRadVolTs << "\t"
                 << SubRateRadVolTs << "\t" << MutRateConsPolTr << "\t"
                 << SubRateConsPolTr << "\t" << MutRateRadPolTr << "\t"
                 << SubRateRadPolTr << "\t" << MutRateConsVolTr << "\t"
                 << SubRateConsVolTr << "\t" << MutRateRadVolTr << "\t"
                 << SubRateRadVolTr << "\t" << MutRateConsPolCpGTs << "\t"
                 << SubRateConsPolCpGTs << "\t" << MutRateRadPolCpGTs << "\t"
                 << SubRateRadPolCpGTs << "\t" << MutRateConsVolCpGTs << "\t"
                 << SubRateConsVolCpGTs << "\t" << MutRateRadVolCpGTs << "\t"
                 << SubRateRadVolCpGTs << "\t" << MutRateCGNTs << "\t"
                 << SubRateCGNTs << "\n";
        rates_os.close();
      }
    }
  }
}
