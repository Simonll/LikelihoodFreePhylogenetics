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

// SimuMutSelSBDP: posterior-predictive / ABC simulator driven by a PhyloBayes
// MSAAGTRW chain (AACODONMUTSELSBDP). Same layout as BayescodeSimuMutSelAA, but
// each parameter point is read directly from <chain>.chain through
// LocalParameters::readChainCodonMutSelSBDP(k) instead of a BayesCode .pvalues.
//
// Chain points visited: start, start+every, ... < until, taken from
// -start/-every/-until in #LOCALPARAM, or from #SAMPLING when that local range
// is empty. #NRUN <Nrun> <threshold>: Nrun replicates per point, at most
// threshold simulations in total.

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

int main(int argc, char *argv[])
{
  // program options
  std::string model = "";
  std::string controlfile = "";

  try
  {
    if (argc < 2)
    {
      throw(0);
    }
    int i = 1;
    while (i < argc)
    {
      std::string s = argv[i];
      if (s == "-v" || s == "--version")
      {
        throw(0);
      }
      else if (s == "-m")
      {
        if (i + 2 >= argc)
        {
          throw(0);
        }
        i++;
        model = argv[i];
        i++;
        controlfile = argv[i];
      }
      i++;
    } // end while
  } // end try
  catch (...)
  {
    std::cerr << "\n";
    std::cerr << "version 1.0\n";
    std::cerr << "###########################\n";
    std::cerr << "-m < CodonMutSelSBDP > < controlfile >\n";
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
    std::cerr << "#NRUN <Nrun per point> <threshold>\n";
    std::cerr << "#SAMPLING <start> <every> <until>\n";
    std::cerr << "#LOCALPARAM\n";
    std::cerr << "###########################\n";
    exit(1);
  }
  if (model == "CodonMutSelSBDP")
  {
    cerr << "simulating under " << model << "\n";

    GlobalParameters *gparam = new GlobalParameters(model, controlfile);
    LocalParameters *lparam = new LocalParameters(gparam);
    Posterior *post = new Posterior(gparam);
    post->SetNsite(lparam->Nsite_codon);

    // Chain range: #LOCALPARAM -start/-every/-until, else #SAMPLING.
    if (lparam->startPoint >= lparam->endPoint)
    {
      lparam->startPoint = gparam->chainPointStart;
      lparam->everyPoint = gparam->chainPointEvery;
      lparam->endPoint = gparam->chainPointEnd;
    }
    if (lparam->everyPoint < 1)
    {
      lparam->everyPoint = 1;
    }
    const int start = lparam->startPoint;
    const int every = lparam->everyPoint;
    const int until = lparam->endPoint;
    std::cerr << "chain range: start " << start << " every " << every
              << " until " << until << "\n";

    lparam->readChainCodonMutSelSBDP(start);
    SummaryStatistics *ss = new SummaryStatistics(lparam);
    ss->computeSummaries();
    PriorSampler *prior = new PriorSampler(lparam);
    SiteInterSubMatrix *submatrix = new SiteInterSubMatrix(lparam);
    submatrix->init();
    AncestralSequence *ancestraseq = new AncestralSequence(lparam);
    TreeSimulator *simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);
    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();
    std::cerr << "starting to simulate\n";
    std::cerr << "Nrun per point: " << gparam->Nrun << " threshold "
              << gparam->threshold << "\n";

    for (int k = start; k < until && post->Niter < gparam->threshold;
         k += every)
    {
      lparam->readChainCodonMutSelSBDP(k);
      for (int i = 0; i < gparam->Nrun && post->Niter < gparam->threshold;
           i++)
      {
        prior->sample();
        simulator->run_jump_chain_over_tree();
        ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);
        post->registerSimulation(
            k, lparam->GetCurrentParameters(), lparam->GetCurrentSummaries(),
            lparam->GetCurrentAccessorySummaries(),
            lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
            lparam->GetCurrentSiteSpecificEvoStats());

        if (lparam->tofasta)
        {
          // For generating the .fasta file
          std::ostringstream oss_fasta;
          oss_fasta << gparam->output << "-pt" << k << "_" << i << ".fasta";
          std::ofstream fasta_os(oss_fasta.str().c_str(), std::ios_base::out);
          lparam->toFastaWithMissingData(
              fasta_os, simulator->CurrentLeafNodeCodonSequences);
          fasta_os.close();

          // For generating the .tre file
          std::ostringstream oss_tre;
          oss_tre << gparam->output << "-pt" << k << "_" << i << ".tre";
          std::ofstream tre_os(oss_tre.str().c_str(), std::ios_base::out);
          lparam->refTree->ToStream(tre_os);
          tre_os.close();
        }
      }
    }
    std::cerr << "\nEnd of the simulation process\n";
    std::cerr << "Nsimu generated: " << post->Niter << "\n";
    std::cerr << "writing files " << (gparam->output + ".simu") << " and "
              << (gparam->output + ".ppp") << "\n";
    ofstream dist_os((gparam->output + ".simu").c_str(), std::ios_base::out);
    post->writeHeader(dist_os);
    post->writePosterior(dist_os);
    dist_os.close();

    ofstream ppp_os((gparam->output + ".ppp").c_str(), std::ios_base::out);
    post->writePosteriorPredictiveStatistics(ppp_os, lparam->summariesRealData);
    ppp_os.close();
    exit(0);
  }
  else
  {
    std::cerr << "This model {" << model << "} is not available yet.\n";
    exit(1);
  }
}
