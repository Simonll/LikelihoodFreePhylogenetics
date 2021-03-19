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
#include "PosteriorSampler.h"

PosteriorSampler::PosteriorSampler() { this->lparam = lparam; }

PosteriorSampler::~PosteriorSampler() {
  // dtor
}

// std::tuple<double,double> PosteriorSampler::sampleFromAjustedDensity(int
// param,double b_inf, double b_sup) {
//
//    double weight_sum = 0.0 ;
//    for (unsigned int simu_i =0 ; simu_i <
//    lparam->gparam->population_t.size(); simu_i++) {
//        weight_sum +=
//        std::get<4>(lparam->gparam->population_t[simu_i])[param];
//
//    }
//    double u = b_inf -1;
//    while (u < b_inf || u > b_sup) {
//        int simu_i = 0 ;
//
//        double testcummul = weight_sum * lparam->rnd->Uniform();
//
//        double test =
//        std::get<4>(lparam->gparam->population_t[simu_i])[param]; while (test
//        < testcummul && simu_i < lparam->gparam->population_t.size()) {
//            simu_i++;
//            test += std::get<4>(lparam->gparam->population_t[simu_i])[param];
//        }
//        if (simu_i >=  lparam->gparam->population_t.size()) {
//            cerr << "error when sampling particules\n";
//            exit(0);
//        }
//        u = std::get<0>(lparam->gparam->population_t[simu_i])[param]  + 2.0 *
//        (lparam->gparam->empVar[param]) * lparam->rnd->sNormal();
//    }
//
//
//
//
//
//    return  std::make_tuple(u,0.0) ;
//
//}

void PosteriorSampler::sample() {
  //            for (unsigned int i = 0 ; i <
  //            lparam->gparam->listUsedParam.size(); i++){
  //                if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "root") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i, 0.0, 1.0);
  //                    lparam->percentFromOutGroup = std::get<0>(cur_tup);
  //                    lparam->SetBranchesLengthsBetweenInAndOutGroup();
  //
  //
  //                } else
  //                if(lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "lambda") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 100); lparam->lambdaCpG
  //                    = std::get<0>(cur_tup);
  //
  //
  //                } else
  //                if(lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "lambdaCpG") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 100); lparam->lambdaCpG
  //                    = std::get<0>(cur_tup);
  //
  //
  //                } else
  //                if(lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "lambdaTpA") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 100); lparam->lambdaTpA
  //                    = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "lambdaTG") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 100); lparam->lambdaTG
  //                    = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "lambdaCA") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 100); lparam->lambdaCA
  //                    = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "mu") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.1, 10); lparam->mu =
  //                    std::get<0>(cur_tup); for(int node = 0; node <
  //                    lparam->refTree->GetNnode(); node++){
  //                        lparam->mu_branch[node] = lparam->mu;
  //                    }
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "muomega") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10); lparam->muomega =
  //                    std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucsA") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10); lparam->nucp[0] =
  //                    std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucsC") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10); lparam->nucp[1] =
  //                    std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucsG") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10); lparam->nucp[2] =
  //                    std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucsT") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10); lparam->nucp[3] =
  //                    std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrAC") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[0][1] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrAG") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[0][2] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrAT") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[0][3] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrCA") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[1][0] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrCG") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[1][2] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrCT") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[1][3] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrGA") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[2][0] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrGC") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[2][1] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrGT") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[2][3] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrTA") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[3][0] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrTC") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[3][1] = std::get<0>(cur_tup);
  //
  //
  //                } else if
  //                (lparam->gparam->listParam[lparam->gparam->listUsedParam[i]]
  //                == "nucrrTG") {
  //                    std::tuple<double,double> cur_tup =
  //                    sampleFromAjustedDensity(i,0.01, 10);
  //                    lparam->nucrrnr[3][2] = std::get<0>(cur_tup);
  //
  //
  //                }
  //            }
  //
  //            double sum =  0.0 ;
  //            for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++){
  //                sum += lparam->nucp[nuc1];
  //            }
  //            for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++){
  //                lparam->nucp[nuc1]/=sum;
  //            }
  //
  //            sum = 0.0 ;
  //            for (int nuc1=0; nuc1<4-1; nuc1++) {
  //                for (int nuc2=nuc1+1; nuc2<4; nuc2++) {
  //                   sum += lparam->nucrrnr[nuc1][nuc2];
  //                }
  //            }
  //
  //            for (int nuc1=0; nuc1<4; nuc1++) {
  //                for (int nuc2=0; nuc2<4; nuc2++) {
  //                   lparam->nucrrnr[nuc1][nuc2] /= sum;
  //                }
  //            }
  //
  //            lparam->getrate = false;
  //
  //            lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
  //            lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
  //            lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
  //            lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
  //            //gtnr[1][1] = 0.0; //cc
  //            lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
  //            lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
  //            lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
  //            lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
  //            //gtnr[2][2] = 0.0; //gg
  //            lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
  //            lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
  //            lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
  //            lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
  //            //gtnr[3][3] = 0.0; //tt
}
