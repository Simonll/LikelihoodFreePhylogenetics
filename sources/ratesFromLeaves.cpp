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
     << "MutRateNonSynCpG"
     << "\t"
     << "SubRateNonSynCpG"
     << "\t"
     << "MutRateSynCpG"
     << "\t"
     << "SubRateSynCpG"
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
     << "MutRateConsPolCpG"
     << "\t"
     << "SubRateConsPolCpG"
     << "\t"
     << "MutRateRadPolCpG"
     << "\t"
     << "SubRateRadPolCpG"
     << "\t"
     << "MutRateConsVolCpG"
     << "\t"
     << "SubRateConsVolCpG"
     << "\t"
     << "MutRateRadVolCpG"
     << "\t"
     << "SubRateRadVolCpG"
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
        double MutRateCpG = 0.0;
        double SubRateCpG = 0.0;
        double MutRateNonSynCpG = 0.0;
        double SubRateNonSynCpG = 0.0;
        double MutRateSynCpG = 0.0;
        double SubRateSynCpG = 0.0;
        double MutRateWeakStrong = 0.0;
        double SubRateWeakStrong = 0.0;
        double MutRateStrongWeak = 0.0;
        double SubRateStrongWeak = 0.0;
        double MutRateWeakWeak = 0.0;
        double SubRateWeakWeak = 0.0;
        double MutRateStrongStrong = 0.0;
        double SubRateStrongStrong = 0.0;
        double MutRateTransition = 0.0;
        double SubRateTransition = 0.0;
        double MutRateTransversion = 0.0;
        double SubRateTransversion = 0.0;
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

        double MutRateConsPolCpG = 0.0;
        double SubRateConsPolCpG = 0.0;
        double MutRateRadPolCpG = 0.0;
        double SubRateRadPolCpG = 0.0;
        double MutRateConsVolCpG = 0.0;
        double SubRateConsVolCpG = 0.0;
        double MutRateRadVolCpG = 0.0;
        double SubRateRadVolCpG = 0.0;

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
          SubRateCpG += SubRateCpG_;

          double MutRateSynCpG_ = 0.0;
          double SubRateSynCpG_ = 0.0;
          std::tie(MutRateSynCpG_, SubRateSynCpG_) = submatrix->GetRatesSynCpG(
              NodeIndex, -1, simulator->CurrentNodeNucSequence);
          MutRateSynCpG += MutRateSynCpG_;
          SubRateSynCpG += SubRateSynCpG_;

          double MutRateNonSynCpG_ = 0.0;
          double SubRateNonSynCpG_ = 0.0;
          std::tie(MutRateNonSynCpG_, SubRateNonSynCpG_) =
              submatrix->GetRatesNonSynCpG(NodeIndex, -1,
                                           simulator->CurrentNodeNucSequence);
          MutRateNonSynCpG += MutRateNonSynCpG_;
          SubRateNonSynCpG += SubRateNonSynCpG_;

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
          MutRateTransversion += MutRateTransversion_;
          SubRateTransversion += SubRateTransversion_;

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

          double MutRateConsPolCpG_ = 0.0;
          double SubRateConsPolCpG_ = 0.0;
          std::tie(MutRateConsPolCpG_, SubRateConsPolCpG_) =
              submatrix->GetRatesConsPolCpG(NodeIndex, -1,
                                            simulator->CurrentNodeNucSequence);
          MutRateConsPolCpG += MutRateConsPolCpG_;
          SubRateConsPolCpG += SubRateConsPolCpG_;

          double MutRateRadPolCpG_ = 0.0;
          double SubRateRadPolCpG_ = 0.0;
          std::tie(MutRateRadPolCpG_, SubRateRadPolCpG_) =
              submatrix->GetRatesRadPolCpG(NodeIndex, -1,
                                           simulator->CurrentNodeNucSequence);
          MutRateRadPolCpG += MutRateRadPolCpG_;
          SubRateRadPolCpG += SubRateRadPolCpG_;

          double MutRateConsVolCpG_ = 0.0;
          double SubRateConsVolCpG_ = 0.0;
          std::tie(MutRateConsVolCpG_, SubRateConsVolCpG_) =
              submatrix->GetRatesConsVolCpG(NodeIndex, -1,
                                            simulator->CurrentNodeNucSequence);
          MutRateConsVolCpG += MutRateConsVolCpG_;
          SubRateConsVolCpG += SubRateConsVolCpG_;

          double MutRateRadVolCpG_ = 0.0;
          double SubRateRadVolCpG_ = 0.0;
          std::tie(MutRateRadVolCpG_, SubRateRadVolCpG_) =
              submatrix->GetRatesRadVolCpG(NodeIndex, -1,
                                           simulator->CurrentNodeNucSequence);
          MutRateRadVolCpG += MutRateRadVolCpG_;
          SubRateRadVolCpG += SubRateRadVolCpG_;

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
        MutRateNonSynCpG /= rep;
        SubRateNonSynCpG /= rep;
        MutRateSynCpG /= rep;
        SubRateSynCpG /= rep;
        MutRateWeakStrong /= rep;
        SubRateWeakStrong /= rep;
        MutRateStrongWeak /= rep;
        SubRateStrongWeak /= rep;
        MutRateStrongStrong /= rep;
        SubRateStrongStrong /= rep;
        MutRateWeakWeak /= rep;
        SubRateWeakWeak /= rep;
        MutRateTransition /= rep;
        SubRateTransition /= rep;
        MutRateTransversion /= rep;
        SubRateTransversion /= rep;
        MutRateConsPol /= rep;
        SubRateConsPol /= rep;
        MutRateRadPol /= rep;
        SubRateRadPol /= rep;
        MutRateConsVol /= rep;
        SubRateConsVol /= rep;
        MutRateRadVol /= rep;
        SubRateRadVol /= rep;
        ofstream rates_os((gparam->output + ".rates").c_str(),
                          std::ios_base::app);
        rates_os << pointID << "\t"
                 << ((seqtype == "stationary") ? "NA"
                                               : lparam->taxonset->GetTaxon(
                                                     ancestraseq->choosen_taxa))
                 << "\t" << MutRate << "\t" << SubRate << "\t" << MutRateNonSyn
                 << "\t" << SubRateNonSyn << "\t" << MutRateSyn << "\t"
                 << SubRateSyn << "\t" << MutRateCpG << "\t" << SubRateCpG
                 << "\t" << MutRateNonSynCpG << "\t" << SubRateNonSynCpG << "\t"
                 << MutRateSynCpG << "\t" << SubRateSynCpG << "\t"
                 << MutRateWeakStrong << "\t" << SubRateWeakStrong << "\t"
                 << MutRateStrongWeak << "\t" << SubRateStrongWeak << "\t"
                 << MutRateStrongStrong << "\t" << SubRateStrongStrong << "\t"
                 << MutRateWeakWeak << "\t" << SubRateWeakWeak << "\t"
                 << MutRateTransition << "\t" << SubRateTransition << "\t"
                 << MutRateTransversion << "\t" << SubRateTransversion << "\t"
                 << MutRateConsPol << "\t" << SubRateConsPol << "\t"
                 << MutRateRadPol << "\t" << SubRateRadPol << "\t"
                 << MutRateConsVol << "\t" << SubRateConsVol << "\t"
                 << MutRateRadVol << "\t" << SubRateRadVol << "\t"
                 << MutRateConsPolTs << "\t" << SubRateConsPolTs << "\t"
                 << MutRateRadPolTs << "\t" << SubRateRadPolTs << "\t"
                 << MutRateConsVolTs << "\t" << SubRateConsVolTs << "\t"
                 << MutRateRadVolTs << "\t" << SubRateRadVolTs << "\t"
                 << MutRateConsPolTr << "\t" << SubRateConsPolTr << "\t"
                 << MutRateRadPolTr << "\t" << SubRateRadPolTr << "\t"
                 << MutRateConsVolTr << "\t" << SubRateConsVolTr << "\t"
                 << MutRateRadVolTr << "\t" << SubRateRadVolTr << "\t"
                 << MutRateConsPolCpG << "\t" << SubRateConsPolCpG << "\t"
                 << MutRateRadPolCpG << "\t" << SubRateRadPolCpG << "\t"
                 << MutRateConsVolCpG << "\t" << SubRateConsVolCpG << "\t"
                 << MutRateRadVolCpG << "\t" << SubRateRadVolCpG << "\n";

        rates_os.close();
      }
    }
  }
}
