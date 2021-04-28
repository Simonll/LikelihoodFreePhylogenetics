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
    std::cerr << "#DIST\n";
    for (auto i : gparam->listDistances) std::cerr << i << "\t";
    std::cerr << "\n";
    std::cerr << "#TRANS\n";
    for (auto i : gparam->listTransformtations) std::cerr << i << "\t";
    std::cerr << "\n";
  } else if (model == "CodonMutSelFiniteABC") {
    std::cerr << "CodonMutSelFiniteCpG\n";
    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    int Npoint =
        static_cast<int>(gparam->chainPointEnd - gparam->chainPointStart) /
        gparam->chainPointEvery;
    Posterior* post3 = new Posterior(gparam);
    std::cerr << gparam->chainPointStart << " " << gparam->chainPointEnd << " "
              << gparam->chainPointEvery << " " << Npoint << "\n";

    LocalParameters** lparam = new LocalParameters*[Npoint];
    SummaryStatistics** ss = new SummaryStatistics*[Npoint];
    PriorSampler** sampler = new PriorSampler*[Npoint];
    SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[Npoint];
    AncestralSequence** ancestraseq = new AncestralSequence*[Npoint];
    TreeSimulator** simulator = new TreeSimulator*[Npoint];

    omp_set_dynamic(0);
    omp_set_num_threads(gparam->Nthread);
    int pt_i;
    std::cerr << "Npoint\t" << Npoint << "\n";
#pragma omp parallel for
    for (pt_i = gparam->chainPointStart; pt_i < gparam->chainPointEnd;
         pt_i += gparam->chainPointEvery) {
      int l = static_cast<int>(pt_i - gparam->chainPointStart) /
              gparam->chainPointEvery;
      lparam[l] = new LocalParameters(gparam);
      lparam[l]->readChainCodonMutSelFinite(pt_i);
      ss[l] = new SummaryStatistics(lparam[l]);
      ss[l]->computeSummaries();
      sampler[l] = new PriorSampler(lparam[l]);
      submatrix[l] = new SiteInterSubMatrix(lparam[l]);
      ancestraseq[l] = new AncestralSequence(lparam[l]);
      simulator[l] = new TreeSimulator(lparam[l], submatrix[l], ancestraseq[l]);

      if (l == 0) {
        post3->SetNsite(lparam[l]->Nsite_codon);

        std::ostringstream ost1;
        ost1 << gparam->output << ".inputparam";
        ofstream lparam_os(ost1.str());
        lparam[l]->writeParam(lparam_os);
        lparam_os.close();

        std::ostringstream ost2;
        ost2 << gparam->output << ".realdata";
        ofstream realDataSummaries_os(ost2.str());
        lparam[l]->writeRealDataSummaries(realDataSummaries_os);
        realDataSummaries_os.close();
      }
    }

    ifstream monitor_is((gparam->output + ".monitor-tmp").c_str());
    if (monitor_is) {
      post3->readMonitor(monitor_is);
      monitor_is.close();

      ifstream posterior_is((gparam->output + ".post-tmp").c_str());
      if (!posterior_is) {
        std::cerr << "error: did not find posteriorfile"
                  << "\n";
        exit(1);
      }

      post3->readPosterior(posterior_is);
      posterior_is.close();
      std::cerr << post3->Niter << " on " << post3->Nrun << "\n";
    }

    std::cerr << "The simulation process started\n";

    while (post3->Niter < post3->Nrun) {
      std::cerr << ".";
      omp_set_dynamic(0);
      omp_set_num_threads(gparam->Nthread);
#pragma omp parallel
      {
#pragma omp for
        for (int l = 0; l < Npoint; l++) {
          sampler[l]->sample();

          simulator[l]->GenerateCodonAlignment();

          ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);

#pragma omp critical
          {
            if (post3->Niter < post3->Nrun) {
              post3->registerNewSimulation(
                  lparam[l]->MCMCpointID, lparam[l]->GetCurrentParameters(),
                  lparam[l]->GetCurrentSummaries(),
                  lparam[l]->GetCurrentAccessorySummaries(),
                  lparam[l]->GetCurrentAncEvoStats(),
                  lparam[l]->GetCurrentEvoStats(),
                  lparam[l]->GetCurrentSiteSpecificEvoStats(),
                  lparam[l]->GetCurrentDistances(),
                  lparam[l]->GetCurrentWeights());

              if (lparam[l]->tofasta) {
                ostringstream oss;
                oss << gparam->output << "-" << post3->Niter << ".fasta";
                std::string output = oss.str();
                ofstream fasta_os((output).c_str(), std::ios_base::out);
                lparam[l]->toFasta(fasta_os,
                                   simulator[l]->CurrentLeafNodeCodonSequences);
                fasta_os.close();
              }
            }
          }
        }
      }
    }

    ofstream dist_os1((gparam->output + ".simu").c_str(), std::ios_base::out);
    post3->writeHeader(dist_os1);
    post3->writePosterior(dist_os1);
    dist_os1.close();

    ofstream monitor_os1((gparam->output + ".monitor").c_str(),
                         std::ios_base::out);
    post3->writeMonitorPosterior(monitor_os1);
    monitor_os1.close();
    std::cerr << "End of the simulation process\n";
    exit(0);
  } else if (model == "CodonMutSelSBDPABC") {
    std::cerr << "CodonMutSelSBDPCpG\n";

    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    std::cerr << gparam->chainPointStart << " " << gparam->chainPointEnd << " "
              << gparam->chainPointEvery << "\n";
    int Npoint =
        static_cast<int>(gparam->chainPointEnd - gparam->chainPointStart) /
        gparam->chainPointEvery;
    Posterior* post3 = new Posterior(gparam);
    LocalParameters** lparam = new LocalParameters*[Npoint];
    std::cerr << "Npoint\t" << Npoint << "\n";
    SummaryStatistics** ss = new SummaryStatistics*[Npoint];
    PriorSampler** sampler = new PriorSampler*[Npoint];
    SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[Npoint];
    AncestralSequence** ancestraseq = new AncestralSequence*[Npoint];
    TreeSimulator** simulator = new TreeSimulator*[Npoint];

    omp_set_dynamic(0);
    omp_set_num_threads(gparam->Nthread);
    int pt_i;
#pragma omp parallel for
    for (pt_i = gparam->chainPointStart; pt_i < gparam->chainPointEnd;
         pt_i += gparam->chainPointEvery) {
      int l = static_cast<int>(pt_i - gparam->chainPointStart) /
              gparam->chainPointEvery;

      lparam[l] = new LocalParameters(gparam);

      lparam[l]->readChainCodonMutSelSBDP(pt_i);

      ss[l] = new SummaryStatistics(lparam[l]);

      ss[l]->computeSummaries();

      sampler[l] = new PriorSampler(lparam[l]);

      submatrix[l] = new SiteInterSubMatrix(lparam[l]);

      ancestraseq[l] = new AncestralSequence(lparam[l]);

      simulator[l] = new TreeSimulator(lparam[l], submatrix[l], ancestraseq[l]);

      if (l == 0) {
        post3->SetNsite(lparam[l]->Nsite_codon);

        std::ostringstream ost1;
        ost1 << gparam->output << ".inputparam";
        ofstream lparam_os(ost1.str());
        lparam[l]->writeParam(lparam_os);
        lparam_os.close();

        std::ostringstream ost2;
        ost2 << gparam->output << ".realdata";
        ofstream realDataSummaries_os(ost2.str());
        lparam[l]->writeRealDataSummaries(realDataSummaries_os);
        realDataSummaries_os.close();
      }
    }

    ifstream monitor_is((gparam->output + ".monitor-tmp").c_str());
    if (monitor_is) {
      post3->readMonitor(monitor_is);
      monitor_is.close();

      ifstream posterior_is((gparam->output + ".post-tmp").c_str());
      if (!posterior_is) {
        std::cerr << "error: did not find posteriorfile"
                  << "\n";
        exit(1);
      }

      post3->readPosterior(posterior_is);
      posterior_is.close();
      std::cerr << post3->Niter << " on " << post3->Nrun << "\n";
    }

    std::cerr << "The simulation process started\n";
    while (post3->Niter < post3->Nrun) {
      std::cerr << ".";
      omp_set_dynamic(0);
      omp_set_num_threads(gparam->Nthread);
#pragma omp parallel
      {
#pragma omp for
        for (int l = 0; l < Npoint; l++) {
          sampler[l]->sample();
          simulator[l]->GenerateCodonAlignment();
          ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);

#pragma omp critical
          {
            post3->registerNewSimulation(
                lparam[l]->MCMCpointID, lparam[l]->GetCurrentParameters(),
                lparam[l]->GetCurrentSummaries(),
                lparam[l]->GetCurrentAccessorySummaries(),
                lparam[l]->GetCurrentAncEvoStats(),
                lparam[l]->GetCurrentEvoStats(),
                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                lparam[l]->GetCurrentDistances(),
                lparam[l]->GetCurrentWeights());

            if (post3->Niter == post3->Nrun) {
              ofstream dist_os1((gparam->output + ".simu").c_str(),
                                std::ios_base::out);
              post3->writeHeader(dist_os1);
              post3->writePosterior(dist_os1);
              dist_os1.close();

              ofstream monitor_os1((gparam->output + ".monitor").c_str(),
                                   std::ios_base::out);
              post3->writeMonitorPosterior(monitor_os1);
              monitor_os1.close();
              std::cerr << "End of the simulation process\n";
              exit(0);
            }
          }
        }
      }
    }
  } else if (model == "CodonMutSelSBDPABC-v2") {
    std::cerr << "CodonMutSelSBDPCpG-v2\n";

    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    std::cerr << gparam->chainPointStart << " " << gparam->chainPointEnd << " "
              << gparam->chainPointEvery << "\n";
    int Npoint =
        static_cast<int>(gparam->chainPointEnd - gparam->chainPointStart) /
        gparam->chainPointEvery;
    Posterior* postMaster = new Posterior(gparam);
    Posterior** postSlave = new Posterior*[Npoint];
    LocalParameters** lparam = new LocalParameters*[Npoint];

    SummaryStatistics** ss = new SummaryStatistics*[Npoint];

    PriorSampler** sampler = new PriorSampler*[Npoint];

    SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[Npoint];

    AncestralSequence** ancestraseq = new AncestralSequence*[Npoint];

    TreeSimulator** simulator = new TreeSimulator*[Npoint];

    omp_set_dynamic(0);
    omp_set_num_threads(gparam->Nthread);
    int pt_i;
#pragma omp parallel for
    for (pt_i = gparam->chainPointStart; pt_i < gparam->chainPointEnd;
         pt_i += gparam->chainPointEvery) {
      int l = static_cast<int>(pt_i - gparam->chainPointStart) /
              gparam->chainPointEvery;

      postSlave[l] = new Posterior(gparam);

      lparam[l] = new LocalParameters(gparam);

      lparam[l]->readChainCodonMutSelSBDP(pt_i);

      ss[l] = new SummaryStatistics(lparam[l]);

      ss[l]->computeSummaries();

      sampler[l] = new PriorSampler(lparam[l]);

      submatrix[l] = new SiteInterSubMatrix(lparam[l]);

      ancestraseq[l] = new AncestralSequence(lparam[l]);

      simulator[l] = new TreeSimulator(lparam[l], submatrix[l], ancestraseq[l]);

      postSlave[l]->SetNsite(lparam[l]->Nsite_codon);

      if (l == 0) {
        std::ostringstream ost1;
        ost1 << gparam->output << ".inputparam";
        ofstream lparam_os(ost1.str());
        lparam[l]->writeParam(lparam_os);
        lparam_os.close();

        std::ostringstream ost2;
        ost2 << gparam->output << ".realdata";
        ofstream realDataSummaries_os(ost2.str());
        lparam[l]->writeRealDataSummaries(realDataSummaries_os);
        realDataSummaries_os.close();
      }
    }

    std::cerr << "The simulation process started\n";
    std::cerr << postMaster->Niter << " on " << postMaster->Nrun << "\n";

    ifstream monitor_is((gparam->output + "-1M.monitor").c_str());
    ifstream monitor_is_100K((gparam->output + "-100K.monitor").c_str());
    if (monitor_is) {
      postMaster->readMonitor(monitor_is);
      monitor_is.close();

      ifstream posterior_is((gparam->output + "-1M.post").c_str());
      if (!posterior_is) {
        std::cerr << "error: did not find posteriorfile"
                  << "\n";
        exit(1);
      }

      postMaster->readPosterior(posterior_is);
      posterior_is.close();

    } else if (monitor_is_100K) {
      postMaster->readMonitor(monitor_is_100K);
      monitor_is_100K.close();

      ifstream posterior_is((gparam->output + "-100K.post").c_str());
      if (!posterior_is) {
        std::cerr << "error: did not find posteriorfile"
                  << "\n";
        exit(1);
      }

      postMaster->readPosterior(posterior_is);
      posterior_is.close();
    }

    int runTodo =
        static_cast<int>((postMaster->Nrun - postMaster->Niter) / Npoint);
    int run = 0;

    omp_set_dynamic(0);
    omp_set_num_threads(gparam->Nthread);

#pragma omp parallel for private(run) collapse(2)
    for (int l = 0; l < Npoint; l++) {
      for (run = 0; run < runTodo; run++) {
        sampler[l]->sample();

        simulator[l]->GenerateCodonAlignment();

        for (int interval_i = 0; interval_i < 11; interval_i++) {
          ss[l]->computeSummariesAncestralSequence(
              simulator[l]->CurrentAncestralCodonSequence[interval_i]);
        }

        ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);

        postSlave[l]->slaveRegisterNewSimulation(
            lparam[l]->MCMCpointID, lparam[l]->GetCurrentParameters(),
            lparam[l]->GetCurrentSummaries(),
            lparam[l]->GetCurrentAccessorySummaries(),
            lparam[l]->GetCurrentAncEvoStats(), lparam[l]->GetCurrentEvoStats(),
            lparam[l]->GetCurrentSiteSpecificEvoStats(),
            lparam[l]->GetCurrentDistances(), lparam[l]->GetCurrentWeights());

        std::cerr << ".";
      }
    }

    for (int l = 0; l < Npoint; l++) {
      postMaster->slaveToMaster(postSlave[l]->population_t);
    }

    postMaster->sortPopulation();

    ofstream dist_os2((gparam->output + "-1M.post").c_str(),
                      std::ios_base::out);
    postMaster->writeHeader(dist_os2);
    postMaster->writePosterior(dist_os2);
    dist_os2.close();
    std::cerr << "End of the simulation process\n";
    exit(0);

  } else if (model == "FMutSelSimu") {
    GlobalParameters* gparam = new GlobalParameters(model, controlfile);

    LocalParameters* lparam = new LocalParameters(gparam);
    lparam->readFMutSelCodeML();

    Posterior* post = new Posterior(gparam);
    post->SetNsite(lparam->Nsite_codon);

    SummaryStatistics* ss = new SummaryStatistics(lparam);

    ss->computeSummaries();

    PriorSampler* sampler = new PriorSampler(lparam);

    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);

    AncestralSequence* ancestraseq = new AncestralSequence(lparam);

    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);

    ofstream lparam_os((gparam->output + ".inputparam").c_str());
    lparam->writeParam(lparam_os);
    lparam_os.close();

    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();

    while (post->Niter < post->Nrun) {
      sampler->sample();

      simulator->GenerateCodonAlignment();

      ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

      ss->computeSummariesAncestralSequence(
          simulator->CurrentAncestralCodonSequence[10]);

      if (post->Niter == 0) {
        ofstream AncestralDataSummaries_os(
            (gparam->output + ".ancestral").c_str(), std::ios_base::out);
        bool headers = true;
        lparam->writeAncestralDataSummaries(AncestralDataSummaries_os, headers);
        AncestralDataSummaries_os.close();
      } else {
        ofstream AncestralDataSummaries_os(
            (gparam->output + ".ancestral").c_str(), std::ios_base::app);
        bool headers = false;
        lparam->writeAncestralDataSummaries(AncestralDataSummaries_os, headers);
        AncestralDataSummaries_os.close();
      }

      post->registerNewSimulation(
          1, lparam->GetCurrentParameters(), lparam->GetCurrentSummaries(),
          lparam->GetCurrentAccessorySummaries(),
          lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
          lparam->GetCurrentSiteSpecificEvoStats(),
          lparam->GetCurrentDistances(), lparam->GetCurrentWeights());

      if (lparam->tofasta) {
        std::ostringstream oss;
        oss << gparam->output << "-" << post->Niter << ".phylip";
        std::string output = oss.str();
        std::ofstream ali_os((output).c_str(), std::ios_base::out);
        lparam->toAli(ali_os, simulator->CurrentLeafNodeCodonSequences);
        ali_os.close();
      }
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

  } else if (model == "CodonMutSelSBDP") {
    std::cerr << "CodonMutSelSBDP"
              << "\n";
    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    Posterior* post = new Posterior(gparam);

    LocalParameters* lparam = new LocalParameters(gparam);
    lparam->readChainCodonMutSelSBDP();

    ofstream lparam_os((gparam->output + ".inputparam").c_str());
    lparam->writeParam(lparam_os);
    lparam_os.close();

    SummaryStatistics* ss = new SummaryStatistics(lparam);
    ss->computeSummaries();

    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();

    PriorSampler* sampler = new PriorSampler(lparam);
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
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

        post->registerNewSimulation(
            lparam->startPoint, lparam->GetCurrentParameters(),
            lparam->GetCurrentSummaries(),
            lparam->GetCurrentAccessorySummaries(),
            lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
            lparam->GetCurrentSiteSpecificEvoStats(),
            lparam->GetCurrentDistances(), lparam->GetCurrentWeights());

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

    ofstream lparam_os((gparam->output + ".inputparam").c_str());
    lparam->writeParam(lparam_os);
    lparam_os.close();

    SummaryStatistics* ss = new SummaryStatistics(lparam);
    ss->computeSummaries();

    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();

    PriorSampler* sampler = new PriorSampler(lparam);
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
    AncestralSequence* ancestraseq = new AncestralSequence(lparam);
    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);

    std::cerr << "Start simulating \n";
    while (lparam->startPoint < lparam->endPoint) {
      lparam->readChainCodonMutSelFinite();
      int iter = 0;
      while (iter < post->Nrun) {
        sampler->sample();

        simulator->GenerateCodonAlignment();

        ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

        post->registerNewSimulation(
            lparam->startPoint, lparam->GetCurrentParameters(),
            lparam->GetCurrentSummaries(),
            lparam->GetCurrentAccessorySummaries(),
            lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
            lparam->GetCurrentSiteSpecificEvoStats(),
            lparam->GetCurrentDistances(), lparam->GetCurrentWeights());

        if (lparam->tofasta) {
          ostringstream oss;
          oss << gparam->output << "-pt" << lparam->startPoint << "-" << iter
              << ".phylip";
          std::string output = oss.str();
          ofstream fasta_os((output).c_str(), std::ios_base::out);
          lparam->toFasta(fasta_os, simulator->CurrentLeafNodeCodonSequences);
          fasta_os.close();
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
    cerr << "writing .simu file"
         << "\n";
    ofstream dist_os((gparam->output + ".simu").c_str(), std::ios_base::out);

    post->writeHeader(dist_os);
    post->writePosterior(dist_os);
    dist_os.close();

    cerr << "writing .ppp file"
         << "\n";
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

    AncestralSequence* ancestraseq = new AncestralSequence(lparam);

    std::cerr << lparam->Nsite_codon << "\n";
    TreeSimulator* simulator =
        new TreeSimulator(lparam, submatrix, ancestraseq);

    post->readPosterior(lparam->posteriorfile);

    std::cerr << "The simulation process started\n";
    if (!post->posterior.empty()) {
      int it = 0;
      while (it < post->threshold) {
        int point =
            static_cast<int>(lparam->rnd->Uniform() * post->posterior.size());
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

          post->registerNewSimulation(
              lparam->GetPointID(), lparam->GetCurrentParameters(),
              lparam->GetCurrentSummaries(),
              lparam->GetCurrentAccessorySummaries(),
              lparam->GetCurrentAncEvoStats(), lparam->GetCurrentEvoStats(),
              lparam->GetCurrentSiteSpecificEvoStats(),
              lparam->GetCurrentDistances(), lparam->GetCurrentWeights());
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

  } else if (model == "CodonMutSelFiniteSeq" || model == "CodonMutSelSBDPSeq") {
    // the chain pointS are extract from the posterior file according to chainID
    std::cerr << model << "\n";

    GlobalParameters* gparam = new GlobalParameters(model, controlfile);
    Posterior* post = new Posterior(gparam);
    LocalParameters* lparam = new LocalParameters(gparam);

    if (model == "CodonMutSelSBDPSeq") {
      lparam->readChainCodonMutSelSBDP();

    } else if (model == "CodonMutSelFiniteSeq") {
      lparam->readChainCodonMutSelFinite();
    }

    SummaryStatistics* ss = new SummaryStatistics(lparam);
    ss->computeSummaries();

    ofstream realDataSummaries_os((gparam->output + ".realdata").c_str());
    lparam->writeRealDataSummaries(realDataSummaries_os);
    realDataSummaries_os.close();
    std::string s = "seq";
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam, s);

    std::cerr << lparam->Nsite_codon << "\n";
    TreeSimulator* simulator = new TreeSimulator(lparam, submatrix);

    post->readPosterior(lparam->posteriorfile);

    std::cerr << "The simulation process started\n";
    if (!post->posterior.empty()) {
      int it = 0;
      while (it < post->threshold) {
        int point =
            static_cast<int>(lparam->rnd->Uniform() * post->posterior.size());
        lparam->SetCurrentParametersFromPosterior(post->posterior, point);

        if (model == "CodonMutSelSBDPPPred") {
          lparam->readChainCodonMutSelSBDP(lparam->GetPointID());

        } else if (model == "CodonMutSelFinitePPred") {
          lparam->readChainCodonMutSelFinite(lparam->GetPointID());
        }

        int rep = 0;

        while (rep < post->Nrun) {
          rep++;

          simulator->GetNewProbSeq();

          if (it == 0) {
            ofstream mutmatrix_A_os((gparam->output + ".mutmatrix_A").c_str(),
                                    std::ios_base::out);
            ofstream selmatrix_A_os((gparam->output + ".selmatrix_A").c_str(),
                                    std::ios_base::out);
            submatrix->WriteSubMatrix(mutmatrix_A_os, selmatrix_A_os, 0);
            mutmatrix_A_os.close();
            selmatrix_A_os.close();

            ofstream mutmatrix_C_os((gparam->output + ".mutmatrix_C").c_str(),
                                    std::ios_base::out);
            ofstream selmatrix_C_os((gparam->output + ".selmatrix_C").c_str(),
                                    std::ios_base::out);
            submatrix->WriteSubMatrix(mutmatrix_C_os, selmatrix_C_os, 1);
            mutmatrix_C_os.close();
            selmatrix_C_os.close();

            ofstream mutmatrix_G_os((gparam->output + ".mutmatrix_G").c_str(),
                                    std::ios_base::out);
            ofstream selmatrix_G_os((gparam->output + ".selmatrix_G").c_str(),
                                    std::ios_base::out);
            submatrix->WriteSubMatrix(mutmatrix_G_os, selmatrix_G_os, 2);
            mutmatrix_G_os.close();
            selmatrix_G_os.close();

            ofstream mutmatrix_T_os((gparam->output + ".mutmatrix_T").c_str(),
                                    std::ios_base::out);
            ofstream selmatrix_T_os((gparam->output + ".selmatrix_T").c_str(),
                                    std::ios_base::out);
            submatrix->WriteSubMatrix(mutmatrix_T_os, selmatrix_T_os, 0);
            mutmatrix_T_os.close();
            selmatrix_T_os.close();
          } else {
            ofstream mutmatrix_A_os((gparam->output + ".mutmatrix_A").c_str(),
                                    std::ios_base::app);
            ofstream selmatrix_A_os((gparam->output + ".selmatrix_A").c_str(),
                                    std::ios_base::app);
            submatrix->WriteSubMatrix(mutmatrix_A_os, selmatrix_A_os, 0);
            mutmatrix_A_os.close();
            selmatrix_A_os.close();

            ofstream mutmatrix_C_os((gparam->output + ".mutmatrix_C").c_str(),
                                    std::ios_base::app);
            ofstream selmatrix_C_os((gparam->output + ".selmatrix_C").c_str(),
                                    std::ios_base::app);
            submatrix->WriteSubMatrix(mutmatrix_C_os, selmatrix_C_os, 1);
            mutmatrix_C_os.close();
            selmatrix_C_os.close();

            ofstream mutmatrix_G_os((gparam->output + ".mutmatrix_G").c_str(),
                                    std::ios_base::app);
            ofstream selmatrix_G_os((gparam->output + ".selmatrix_G").c_str(),
                                    std::ios_base::app);
            submatrix->WriteSubMatrix(mutmatrix_G_os, selmatrix_G_os, 2);
            mutmatrix_G_os.close();
            selmatrix_G_os.close();

            ofstream mutmatrix_T_os((gparam->output + ".mutmatrix_T").c_str(),
                                    std::ios_base::app);
            ofstream selmatrix_T_os((gparam->output + ".selmatrix_T").c_str(),
                                    std::ios_base::app);
            submatrix->WriteSubMatrix(mutmatrix_T_os, selmatrix_T_os, 3);
            mutmatrix_T_os.close();
            selmatrix_T_os.close();
          }

          it++;
        }
        std::cerr << ".";
      }
    }
  }
}
