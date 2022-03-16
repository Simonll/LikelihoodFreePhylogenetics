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

#include <omp.h>

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
#include "SiteInterSubMatrix.h"
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
    std::cerr
        << "-m < stats | show | CodonMutSelFiniteABC | CodonMutSelSBDPABC | "
           "CodonDegMutSelFiniteABC | CodonDegMutSelSBDPABC | "
           "CodonMutSelFinite | CodonMutSelSBDP | CodonMutSelFinitePPred | "
           "CodonMutSelSBDPPPred | CodonMutSelSBDPSeq > <controlfile>\n";
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

  std::cerr << "models\n";
  if (model == "stats") {
    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    // Posterior* post = new Posterior(gparam);

    LocalData* ldata = new LocalData(gparam);
    ldata->readLocalData(1);
    SummaryStatistics* ss = new SummaryStatistics(ldata);
    ofstream realDataSummaries_os((ldata->output + ".stats").c_str());
    int k = 1;

    while (k < static_cast<int>(gparam->listGenes.size())) {
      std::cerr << "#########\n";
      std::cerr << k << "\n";
      std::cerr << "#########\n";

      ldata->readLocalData(k);

      ss->computeSummariesFromData();

      ldata->writeRealDataSummaries(realDataSummaries_os, k == 1);
      k++;
    }
    realDataSummaries_os.close();
  } else if (model == "show") {
    std::cerr << "show program content\n";
    GlobalParameters* gparam = new GlobalParameters();
    std::cerr << "#PARAMETERS\n";
    for (auto i : gparam->listParam) std::cerr << i << "\t";
    std::cerr << "\n";
    std::cerr << "#SUMMARY\n";
    for (auto i : gparam->listSummaries) std::cerr << i << "\t";
    std::cerr << "\n";
    std::cerr << "#EVOLSTATS\n";
    for (auto i : gparam->listEvoStats) std::cerr << i << "\t";
    std::cerr << "\n";
  } else if (model == "CodonMutSelSBDP") {
    std::cerr << "CodonMutSelSBDP"
              << "\n";
    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    Posterior* post = new Posterior(gparam);

    LocalParameters* lparam = new LocalParameters(gparam);
    lparam->readChainCodonMutSelSBDP();

    SummaryStatistics* ss = new SummaryStatistics(lparam);
    ss->computeSummaries();

    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();

    PriorSampler* sampler = new PriorSampler(lparam);
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
    submatrix->init();
    AncestralSequence* ancestraseq = new AncestralSequence(lparam);
    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);

    std::cerr << "Start simulating \n";
    while (lparam->startPoint < lparam->endPoint) {
      lparam->readChainCodonMutSelSBDP();

      int iter = 0;
      while (iter < post->Nrun) {
        sampler->sample();
        simulator->GenerateCodonAlignment();

        ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

        if (lparam->sampleAncSeq) {
          for (int sample_i = 0; sample_i < lparam->Ninterval; sample_i++) {
            ss->computeSummariesAncestralSequence(
                simulator->CurrentAncestralCodonSequence[sample_i]);

            if (post->Niter == 0) {
              ofstream AncestralDataSummaries_os(
                  (gparam->output + ".ancestral").c_str(), std::ios_base::out);
              bool headers = true;
              lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,
                                                  headers);
              AncestralDataSummaries_os.close();
            } else {
              ofstream AncestralDataSummaries_os(
                  (gparam->output + ".ancestral").c_str(), std::ios_base::app);
              bool headers = false;
              lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,
                                                  headers);
              AncestralDataSummaries_os.close();
            }
          }
        }

        post->registerSimulation(
            lparam->startPoint, lparam->GetCurrentParameters(),
            lparam->GetCurrentSummaries(),
            lparam->GetCurrentAccessorySummaries(),
            lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
            lparam->GetCurrentSiteSpecificEvoStats());

        if (lparam->tofasta) {
          ostringstream oss;
          oss << gparam->output << "-pt" << lparam->startPoint << "-" << iter
              << ".phylip";
          std::string output = oss.str();
          ofstream ali_os((output).c_str(), std::ios_base::out);
          lparam->toAli(ali_os, simulator->CurrentLeafNodeCodonSequences);
          ali_os.close();
        }
        iter++;
      }

      int rep = 1;
      while (rep < lparam->everyPoint) {
        rep++;
        lparam->incrementStartPoint();
      }
      lparam->incrementStartPoint();
      std::cerr << ".";
    }
    std::cerr << "End of the simulation process\n";

    ofstream dist_os((gparam->output + ".simu").c_str(), std::ios_base::out);
    post->writeHeader(dist_os);
    post->writePosterior(dist_os);
    dist_os.close();

    ofstream ppp_os((gparam->output + ".ppp").c_str(), std::ios_base::out);
    post->writePosteriorPredictiveStatistics(ppp_os, lparam->summariesRealData);
    ppp_os.close();

    exit(0);

  } else if (model == "CodonMutSelFinite") {
    std::cerr << "CodonMutSelFinite"
              << "\n";
    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    Posterior* post = new Posterior(gparam);
    LocalParameters* lparam = new LocalParameters(gparam);
    lparam->readChainCodonMutSelFinite();

    SummaryStatistics* ss = new SummaryStatistics(lparam);
    ss->computeSummaries();

    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();

    PriorSampler* prior = new PriorSampler(lparam);
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
    submatrix->init();
    AncestralSequence* ancestraseq = new AncestralSequence(lparam);
    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);

    std::cerr << "Start simulating \n";
    std::cerr << "Nsimu to be generated: " << gparam->Nsimu << " with Nrep "
              << gparam->Nrep << "\n";

    while (post->Niter < gparam->Nsimu) {
      int k =
          static_cast<int>(lparam->rnd->Uniform() * (gparam->chainPointEnd -
                                                     gparam->chainPointStart) +
                           gparam->chainPointStart);
      lparam->readChainCodonMutSelFinite(k);
      for (int i = 0; i < gparam->Nrep; i++) {
        prior->sample();
        simulator->GenerateCodonAlignment();
        ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);
        post->registerSimulation(
            k, lparam->GetCurrentParameters(), lparam->GetCurrentSummaries(),
            lparam->GetCurrentAccessorySummaries(),
            lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
            lparam->GetCurrentSiteSpecificEvoStats());

        if (lparam->tofasta) {
          ostringstream oss;
          oss << gparam->output << "-" << post->Niter << "_" << i << ".fasta";
          std::string output = oss.str();
          ofstream fasta_os((output).c_str(), std::ios_base::out);
          lparam->toFasta(fasta_os, simulator->CurrentLeafNodeCodonSequences);
          fasta_os.close();
        }
        std::cerr << ".";
      }
    }
    std::cerr << "End of the simulation process\n";

    ofstream dist_os((gparam->output + ".simu").c_str(), std::ios_base::out);
    post->writeHeader(dist_os);
    post->writeSimu(dist_os);
    dist_os.close();

    ofstream ppp_os((gparam->output + ".ppp").c_str(), std::ios_base::out);
    post->writePosteriorPredictiveStatistics(ppp_os, lparam->summariesRealData);
    ppp_os.close();

    exit(0);

  } else if (model == "CodonMutSelFinitePPPred" ||
             model == "CodonMutSelSBDPPPred") {
    // the chain points are extract from the posterior file according to chainID
    std::cerr << model << "\n";

    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    Posterior* post = new Posterior(gparam);
    LocalParameters* lparam = new LocalParameters(gparam);

    if (model == "CodonMutSelSBDPPPred") {
      lparam->readChainCodonMutSelSBDP();

    } else if (model == "CodonMutSelFinitePPred") {
      lparam->readChainCodonMutSelFinite();
    }
    SummaryStatistics* ss = new SummaryStatistics(lparam);
    ss->computeSummaries();
    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
    submatrix->init();
    AncestralSequence* ancestraseq = new AncestralSequence(lparam);
    std::cerr << lparam->Nsite_codon << "\n";
    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);
    post->readPosterior(lparam->posteriorfile);
    std::cerr << "The simulation process started\n";
    if (!post->posterior.empty()) {
      int it = 0;
      while (it < post->threshold) {
        int point = static_cast<int>(
            lparam->rnd->Uniform() * post->posterior.size() - 1);
        lparam->SetCurrentParametersFromPosterior(post->posterior, point);

        if (model == "CodonMutSelSBDPPPred") {
          lparam->readChainCodonMutSelSBDP(lparam->GetPointID());

        } else if (model == "CodonMutSelFinitePPred") {
          lparam->readChainCodonMutSelFinite(lparam->GetPointID());
        }
        int rep = 0;
        while (rep < post->Nrun) {
          rep++;

          simulator->GenerateCodonAlignment();

          ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

          post->registerSimulation(
              lparam->GetPointID(), lparam->GetCurrentParameters(),
              lparam->GetCurrentSummaries(),
              lparam->GetCurrentAccessorySummaries(),
              lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
              lparam->GetCurrentSiteSpecificEvoStats());
          it++;
        }

        std::cerr << ".";
      }

      std::cerr << "End of the simulation process\n";

      std::ofstream dist_os((gparam->output + ".simu").c_str(),
                            std::ios_base::out);
      post->writeHeader(dist_os);
      post->writePosterior(dist_os);
      dist_os.close();

      std::ofstream ppp_os((gparam->output + ".ppp").c_str(),
                           std::ios_base::out);
      post->writePosteriorPredictiveStatistics(ppp_os,
                                               lparam->summariesRealData);
      ppp_os.close();
    }
  }
}
