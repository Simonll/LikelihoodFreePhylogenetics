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
    std::cerr << "-m < show | M7 | M8 | M8a > < controlfile >\n";
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
  if (model == "M7" || model == "M8" || model == "M8a") {
    cerr << "simulating under " << model << "\n";

    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    LocalParameters* lparam = new LocalParameters(gparam);
    Posterior* post = new Posterior(gparam);
    post->SetNsite(lparam->Nsite_codon);
    SummaryStatistics* ss = new SummaryStatistics(lparam);

    lparam->readParametersCodemlM7M8(0);

    ss->computeSummaries();
    PriorSampler* prior = new PriorSampler(lparam);

    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);

    AncestralSequence* ancestraseq = new AncestralSequence(lparam);

    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);
    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();
    cerr << "starting to simulate\n";
    int k = 0;
    while (post->Niter < post->Nrun) {
      lparam->readParametersCodemlM7M8(k);
      prior->sample();
      simulator->GetNewSimulatedCodonAlignment();
      ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);
      post->registerNewSimulation(
          1, lparam->GetCurrentParameters(), lparam->GetCurrentSummaries(),
          lparam->GetCurrentAccessorySummaries(),
          lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
          lparam->GetCurrentSiteSpecificEvoStats(),
          lparam->GetCurrentDistances(), lparam->GetCurrentWeights());

      if (lparam->tofasta) {
        ostringstream oss;
        oss << gparam->output << "-" << post->Niter << ".fasta";
        std::string output = oss.str();
        ofstream fasta_os((output).c_str(), std::ios_base::out);
        lparam->toFasta(fasta_os, simulator->CurrentLeafNodeCodonSequences);
        fasta_os.close();
        std::cerr << ".";
      }
      k++;
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
  }
}
