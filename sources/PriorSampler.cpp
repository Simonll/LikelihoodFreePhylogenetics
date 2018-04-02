/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/
#include "PriorSampler.h"


PriorSampler::PriorSampler(LocalParameters* lparam)
{
    this->lparam = lparam;
}

PriorSampler::~PriorSampler()
{
    //dtor
}


void PriorSampler::sample()
{
/*     if (lparam->fixwR_CHQW != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freewR_CHQW\n";
        }
        lparam->wR_CHQW = lparam->rnd->Uniform();
    } */

    if (lparam->fixfitCpG !=1)
    {
        if (lparam->verbose)
        {
            cerr << "freefitCpG\n";
        }
        lparam->fitCpG = lparam->rnd->Uniform();
    }

    if (lparam->fixfitTpA !=1)
    {
        if (lparam->verbose)
        {
            cerr << "freefitTpA\n";
        }
        lparam->fitTpA = lparam->rnd->Uniform();
    }

    if (lparam->fixroot == 0)
    {
        lparam->percentFromOutGroup = lparam->rnd->Uniform();
        lparam->SetBranchesLengthsBetweenInAndOutGroup();
    }

    if(lparam->fixlambda_omega != 1)
    {
        if (lparam->lambda_omega_prior == "log2Unif")
        {
            lparam->lambda_omega = logNUnif(2);
        }
        else if (lparam->lambda_omega_prior == "log3Unif")
        {
            lparam->lambda_omega = logNUnif(3);
        }
    }

    if(lparam->fixlambda_TBL != 1)
    {
        if (lparam->lambda_TBL_prior == "log2Unif")
        {
            lparam->lambda_TBL = logNUnif(2);
        }
        else if (lparam->lambda_TBL_prior == "log3Unif")
        {
            lparam->lambda_TBL = logNUnif(3);
        }


//                    for(int node = 0; node < lparam->refTree->GetNnode(); node++){
//                        lparam->muBranch[node] = lparam->lambda_TBL;
//                    }
    }

    if(lparam->fixlambda_CpG != 1)
    {
        if (lparam->lambda_CpG_prior == "log10Unif")
        {
            lparam->lambda_CpG = logNUnif(10);
        }
        else if (lparam->lambda_CpG_prior == "log20Unif")
        {
            lparam->lambda_CpG = logNUnif(20);
        }
        else if (lparam->lambda_CpG_prior == "log50Unif")
        {
            lparam->lambda_CpG = logNUnif(50);
        }
        else if (lparam->lambda_CpG_prior == "log100Unif")
        {
            lparam->lambda_CpG = logNUnif(100);
        }
    }

/*     if(lparam->fixlambda_CpG_GpG != 1)
    {
        if (lparam->lambda_CpG_GpG_prior == "log10Unif")
        {
            lparam->lambda_CpG_GpG = logNUnif(10);
        }
        else if (lparam->lambda_CpG_GpG_prior == "log20Unif")
        {
            lparam->lambda_CpG_GpG = logNUnif(20);
        }
        else if (lparam->lambda_CpG_GpG_prior == "log50Unif")
        {
            lparam->lambda_CpG_GpG = logNUnif(50);
        }
        else if (lparam->lambda_CpG_GpG_prior == "log100Unif")
        {
            lparam->lambda_CpG_GpG = logNUnif(100);
        }
    }
 */

    if(lparam->fixlambda_TpA != 1)
    {
        if (lparam->lambda_TpA_prior == "log10Unif")
        {
            lparam->lambda_TpA = logNUnif(10);

        }
        else if (lparam->lambda_TpA_prior == "log20Unif")
        {
            lparam->lambda_TpA = logNUnif(20);

        }
        else if (lparam->lambda_TpA_prior == "log50Unif")
        {
            lparam->lambda_TpA = logNUnif(50);

        }
        else if (lparam->lambda_TpA_prior == "log100Unif")
        {
            lparam->lambda_TpA = logNUnif(100);
        }
    }

    if(lparam->fixlambda_tvTpA != 1)
    {
        lparam->lambda_tvTpA = logNUnif(10);
    }

    if(lparam->fixlambda_tstvTpA != 1)
    {
        lparam->lambda_tstvTpA = logNUnif(10);
    }

    if(lparam->fixlambda_tvCpG != 1)
    {
        lparam->lambda_tvCpG = logNUnif(10);
    }

    if(lparam->fixlambda_tstvCpG != 1)
    {
        lparam->lambda_tstvCpG = logNUnif(10);
    }

/*     if(lparam->fixlambda_GpT != 1)
    {

        if (lparam->lambda_GpT_prior == "log10Unif")
        {
            lparam->lambda_GpT = logNUnif(10);

        }
        else if (lparam->lambda_GpT_prior == "log20Unif")
        {
            lparam->lambda_GpT = logNUnif(20);

        }
        else if (lparam->lambda_GpT_prior == "log50Unif")
        {
            lparam->lambda_GpT = logNUnif(50);

        }
        else if (lparam->lambda_GpT_prior == "log100Unif")
        {
            lparam->lambda_GpT = logNUnif(100);
        }


    } */

    if (lparam->fixrr != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freenucrr\n";
        }
        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = lparam->rnd->Gamma(1,2);  //ac
        lparam->nucrrnr[0][2] = lparam->rnd->Gamma(2,2);  //ag
        lparam->nucrrnr[0][3] = lparam->rnd->Gamma(1,2);  //at
        lparam->nucrrnr[1][0] = lparam->nucrrnr[0][1]; //ca
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = lparam->rnd->Gamma(1,2);  //cg
        lparam->nucrrnr[1][3] = lparam->rnd->Gamma(2,2); ; //ct
        lparam->nucrrnr[2][0] = lparam->nucrrnr[0][2]; //ga
        lparam->nucrrnr[2][1] = lparam->nucrrnr[1][2]; //gc
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = lparam->rnd->Gamma(1,2);  //gt
        lparam->nucrrnr[3][0] = lparam->nucrrnr[0][3]; //ta
        lparam->nucrrnr[3][1] = lparam->nucrrnr[1][3]; //tc
        lparam->nucrrnr[3][2] = lparam->nucrrnr[2][3]; //tg
        //nucrrnr[3][3] = 0.0; //tt


        double sum = 0.0;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }

    else if (lparam->fixkappa != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freekappa\n";
        }
        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = 1.0; //ac
        lparam->nucrrnr[0][2] = Unif(1.0,10.0); //ag
        lparam->nucrrnr[0][3] = 1.0; //at
        lparam->nucrrnr[1][0] = 1.0; //ca
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = 1.0; //cg
        lparam->nucrrnr[1][3] = lparam->nucrrnr[0][2]; //ct
        lparam->nucrrnr[2][0] = lparam->nucrrnr[0][2]; //ga
        lparam->nucrrnr[2][1] = 1.0; //gc
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = 1.0; //gt
        lparam->nucrrnr[3][0] = 1.0; //ta
        lparam->nucrrnr[3][1] = lparam->nucrrnr[0][2]; //tc
        lparam->nucrrnr[3][2] = 1.0; //tg
        //nucrrnr[3][3] = 0.0; //tt

        double sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;


        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }
    else if (lparam->fixhky != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freehky\n";
        }

        lparam->nucp[0] = lparam->rnd->sExpo();
        lparam->nucp[1] = lparam->rnd->sExpo();
        lparam->nucp[2] = lparam->rnd->sExpo();
        lparam->nucp[3] = lparam->rnd->sExpo();

        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = 1.0; //ac
        lparam->nucrrnr[0][2] = Unif(1.0,10.0); //ag
        lparam->nucrrnr[0][3] = 1.0; //at
        lparam->nucrrnr[1][0] = 1.0; //ca
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = 1.0; //cg
        lparam->nucrrnr[1][3] = lparam->nucrrnr[0][2]; //ct
        lparam->nucrrnr[2][0] = lparam->nucrrnr[0][2]; //ga
        lparam->nucrrnr[2][1] = 1.0; //gc
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = 1.0; //gt
        lparam->nucrrnr[3][0] = 1.0; //ta
        lparam->nucrrnr[3][1] = lparam->nucrrnr[0][2]; //tc
        lparam->nucrrnr[3][2] = 1.0; //tg
        //nucrrnr[3][3] = 0.0; //tt

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            sum += lparam->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            lparam->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt
    }
    else if (lparam->fixstat != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freestat\n";
        }
        lparam->nucp[0] = lparam->rnd->sExpo();
        lparam->nucp[1] = lparam->rnd->sExpo();
        lparam->nucp[2] = lparam->rnd->sExpo();
        lparam->nucp[3] = lparam->rnd->sExpo();

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            sum += lparam->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            lparam->nucp[nuc1]/=sum;
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }

    else if (lparam->fixtr != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freetr\n";
        }

        lparam->nucp[0] = lparam->rnd->sExpo();
        lparam->nucp[1] = lparam->rnd->sExpo();
        lparam->nucp[2] = lparam->rnd->sExpo();
        lparam->nucp[3] = lparam->rnd->sExpo();
        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = lparam->rnd->Gamma(1,2);  //ac
        //lparam->nucrrnr[0][2] = lparam->rnd->Gamma(2,2);  //ag
        lparam->nucrrnr[0][3] = lparam->rnd->Gamma(1,2);  //at
        lparam->nucrrnr[1][0] = lparam->nucrrnr[0][1];    //ca
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = lparam->rnd->Gamma(1,2);  //cg
        //lparam->nucrrnr[1][3] = lparam->rnd->Gamma(2,2);  //ct
        //lparam->nucrrnr[2][0] = lparam->nucrrnr[0][2];    //ga
        lparam->nucrrnr[2][1] = lparam->nucrrnr[1][2];    //gc
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = lparam->rnd->Gamma(1,2);  //gt
        lparam->nucrrnr[3][0] = lparam->nucrrnr[0][3];    //ta
        //lparam->nucrrnr[3][1] = lparam->nucrrnr[1][3];    //tc
        lparam->nucrrnr[3][2] = lparam->nucrrnr[2][3];    //tg
        //nucrrnr[3][3] = 0.0; //tt

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            sum += lparam->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            lparam->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }

    else  if (lparam->fixts != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freets\n";
        }

        lparam->nucp[0] = lparam->rnd->sExpo();
        lparam->nucp[1] = lparam->rnd->sExpo();
        lparam->nucp[2] = lparam->rnd->sExpo();
        lparam->nucp[3] = lparam->rnd->sExpo();
        //nucrrnr[0][0] = 0.0; //AA
        //lparam->nucrrnr[0][1] = lparam->rnd->Gamma(1,2);  //ac
        lparam->nucrrnr[0][2] = lparam->rnd->Gamma(2,2);  //ag
        //lparam->nucrrnr[0][3] = lparam->rnd->Gamma(1,2);  //at
        //lparam->nucrrnr[1][0] = lparam->nucrrnr[0][1];    //ca
        //nucrrnr[1][1] = 0.0; //cc
        //lparam->nucrrnr[1][2] = lparam->rnd->Gamma(1,2);  //cg
        lparam->nucrrnr[1][3] = lparam->rnd->Gamma(2,2);  //ct
        lparam->nucrrnr[2][0] = lparam->nucrrnr[0][2];    //ga
        //lparam->nucrrnr[2][1] = lparam->nucrrnr[1][2];    //gc
        //nucrrnr[2][2] = 0.0; //gg
        //lparam->nucrrnr[2][3] = lparam->rnd->Gamma(1,2);  //gt
        //lparam->nucrrnr[3][0] = lparam->nucrrnr[0][3];    //ta
        lparam->nucrrnr[3][1] = lparam->nucrrnr[1][3];    //tc
        //lparam->nucrrnr[3][2] = lparam->nucrrnr[2][3];    //tg
        //nucrrnr[3][3] = 0.0; //tt

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            sum += lparam->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            lparam->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }
    else if (lparam->fixgtr != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freegtr\n";
        }

        lparam->nucp[0] = lparam->rnd->sExpo();
        lparam->nucp[1] = lparam->rnd->sExpo();
        lparam->nucp[2] = lparam->rnd->sExpo();
        lparam->nucp[3] = lparam->rnd->sExpo();
        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = lparam->rnd->Gamma(1,2);  //ac
        lparam->nucrrnr[0][2] = lparam->rnd->Gamma(2,2);  //ag
        lparam->nucrrnr[0][3] = lparam->rnd->Gamma(1,2);  //at
        lparam->nucrrnr[1][0] = lparam->nucrrnr[0][1];    //ca
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = lparam->rnd->Gamma(1,2);  //cg
        lparam->nucrrnr[1][3] = lparam->rnd->Gamma(2,2);  //ct
        lparam->nucrrnr[2][0] = lparam->nucrrnr[0][2];    //ga
        lparam->nucrrnr[2][1] = lparam->nucrrnr[1][2];    //gc
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = lparam->rnd->Gamma(1,2);  //gt
        lparam->nucrrnr[3][0] = lparam->nucrrnr[0][3];    //ta
        lparam->nucrrnr[3][1] = lparam->nucrrnr[1][3];    //tc
        lparam->nucrrnr[3][2] = lparam->nucrrnr[2][3];    //tg
        //nucrrnr[3][3] = 0.0; //tt

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            sum += lparam->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            lparam->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }

    else if (lparam->fixgtr2 == 0)
    {
        if (lparam->verbose)
        {
            cerr << "freegtr2\n";
        }

        lparam->nucp[0] = lparam->rnd->sExpo();
        lparam->nucp[1] = lparam->rnd->sExpo();
        lparam->nucp[2] = lparam->rnd->sExpo();
        lparam->nucp[3] = lparam->rnd->sExpo();
        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = lparam->rnd->Gamma(1,2);  //ac
        lparam->nucrrnr[0][2] = lparam->rnd->Gamma(2,2);  //ag
        lparam->nucrrnr[0][3] = lparam->rnd->Gamma(1,2);  //at
        lparam->nucrrnr[1][0] = lparam->nucrrnr[0][1];    //ca
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = lparam->rnd->Gamma(1,2);  //cg
        lparam->nucrrnr[1][3] = lparam->rnd->Gamma(2,2);  //ct
        lparam->nucrrnr[2][0] = lparam->nucrrnr[0][2];    //ga
        lparam->nucrrnr[2][1] = lparam->nucrrnr[1][2];    //gc
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = lparam->rnd->Gamma(1,2);  //gt
        lparam->nucrrnr[3][0] = lparam->nucrrnr[0][3];    //ta
        lparam->nucrrnr[3][1] = lparam->nucrrnr[1][3];    //tc
        lparam->nucrrnr[3][2] = lparam->nucrrnr[2][3];    //tg
        //nucrrnr[3][3] = 0.0; //tt

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            sum += lparam->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            lparam->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }

    else if (lparam->fixss != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freess\n";
        }

        lparam->nucp[0] = lparam->rnd->sExpo();
        lparam->nucp[1] = lparam->rnd->sExpo();
        lparam->nucp[2] = lparam->rnd->sExpo();
        lparam->nucp[3] = lparam->rnd->sExpo();
        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = lparam->rnd->Gamma(1,2);  //at2gc
        lparam->nucrrnr[0][2] = lparam->rnd->Gamma(2,2);  //at2gc
        lparam->nucrrnr[0][3] = lparam->rnd->Gamma(1,2);  //at2ta
        lparam->nucrrnr[1][0] = lparam->rnd->Gamma(1,2);  //cg2at
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = lparam->rnd->Gamma(1,2);  //cg2gc
        lparam->nucrrnr[1][3] = lparam->rnd->Gamma(2,2);  //cg2ta
        lparam->nucrrnr[2][0] = lparam->nucrrnr[1][3];     //gc2at
        lparam->nucrrnr[2][1] = lparam->nucrrnr[1][2];    //gc2cg
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = lparam->nucrrnr[0][1];    //gc2ta
        lparam->nucrrnr[3][0] = lparam->nucrrnr[0][3];    //ta2at
        lparam->nucrrnr[3][1] = lparam->nucrrnr[0][2];    //ta2cg
        lparam->nucrrnr[3][2] = lparam->nucrrnr[0][1];    //ta2gc
        //nucrrnr[3][3] = 0.0; //tt

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            sum += lparam->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            lparam->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }

    else if (lparam->fixgtnr != 1)
    {
        if (lparam->verbose)
        {
            cerr << "freegtnr\n";
        }
   /*      lparam->nucp[0] = 1.0;
        lparam->nucp[1] = 1.0;
        lparam->nucp[2] = 1.0;
        lparam->nucp[3] = 1.0; 
    */
        //nucrrnr[0][0] = 0.0; //AA
        lparam->nucrrnr[0][1] = lparam->rnd->Gamma(1,2);  //ac
        lparam->nucrrnr[0][2] = lparam->rnd->Gamma(2,2);  //ag
        lparam->nucrrnr[0][3] = lparam->rnd->Gamma(1,2);  //at
        lparam->nucrrnr[1][0] = lparam->rnd->Gamma(1,2);  //ca
        //nucrrnr[1][1] = 0.0; //cc
        lparam->nucrrnr[1][2] = lparam->rnd->Gamma(1,2);  //cg
        lparam->nucrrnr[1][3] = lparam->rnd->Gamma(2,2);  //ct
        lparam->nucrrnr[2][0] = lparam->rnd->Gamma(2,2);  //ga
        lparam->nucrrnr[2][1] = lparam->rnd->Gamma(1,2);  //gc
        //nucrrnr[2][2] = 0.0; //gg
        lparam->nucrrnr[2][3] = lparam->rnd->Gamma(1,2);  //gt
        lparam->nucrrnr[3][0] = lparam->rnd->Gamma(1,2);  //ta
        lparam->nucrrnr[3][1] = lparam->rnd->Gamma(2,2);  //tc
        lparam->nucrrnr[3][2] = lparam->rnd->Gamma(1,2);  //tg
        //nucrrnr[3][3] = 0.0; //tt

        double sum =  0.0 ;
        // for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        // {
        //     sum += lparam->nucp[nuc1];
        // }
        // for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        // {
        //     lparam->nucp[nuc1]/=sum;
        // }

        sum = 0.0 ;
        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            for (int nuc2 = 0 ; nuc2 < lparam->Nnucp; nuc2++)
            {
                sum += lparam->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1 = 0 ; nuc1 < lparam->Nnucp; nuc1++)
        {
            for (int nuc2 = 0 ; nuc2 < lparam->Nnucp; nuc2++)
            {
                lparam->nucrrnr[nuc1][nuc2] /= sum;
            }
        }

        lparam->getrate = false;

        lparam->gtnr[0][1] = lparam->GetGTNR(0,1); //ac
        lparam->gtnr[0][2] = lparam->GetGTNR(0,2); //ag
        lparam->gtnr[0][3] = lparam->GetGTNR(0,3); //at
        lparam->gtnr[1][0] = lparam->GetGTNR(1,0); //ca
        //gtnr[1][1] = 0.0; //cc
        lparam->gtnr[1][2] = lparam->GetGTNR(1,2); //cg
        lparam->gtnr[1][3] = lparam->GetGTNR(1,3); //ct
        lparam->gtnr[2][0] = lparam->GetGTNR(2,0); //ga
        lparam->gtnr[2][1] = lparam->GetGTNR(2,1); //gc
        //gtnr[2][2] = 0.0; //gg
        lparam->gtnr[2][3] = lparam->GetGTNR(2,3); //gt
        lparam->gtnr[3][0] = lparam->GetGTNR(3,0); //ta
        lparam->gtnr[3][1] = lparam->GetGTNR(3,1); //tc
        lparam->gtnr[3][2] = lparam->GetGTNR(3,2); //tg
        //gtnr[3][3] = 0.0; //tt

    }
}


