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
    std::cerr << "-m < show | CodonMutSelFiniteSeq | CodonMutSelSBDPSeq \n";
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
  if (model == "CodonMutSelFiniteSeq" || model == "CodonMutSelSBDPSeq") {
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
    SiteInterSubMatrix* submatrix = new SiteInterSubMatrixCABC2018(lparam, s);

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
