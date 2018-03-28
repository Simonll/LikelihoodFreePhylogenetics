/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/
#include "Regression.h"

Regression::Regression(GlobalParameters* gparam, Posterior* post)
{
    this->gparam = gparam;
    this->post = post;


    std::vector<std::vector<double>> B_HAT(gparam->NusedSummaries+1,vector<double>(gparam->NusedParam,0.0));


}

Regression::~Regression()
{
    //dtor
}

void Regression::ComputeMultipleRegression()
{

    std::vector<std::vector<double>> X      = post->GetPartialDistances();
    std::vector<double> W                   = post->GetWeights();
    std::vector<std::vector<double>> X_T    = post->GetPartialDistancesT();
    std::vector<std::vector<double>> Y      = post->GetTheta();
    std::vector<std::vector<double>> X_TW;
    std::vector<std::vector<double>> X_TWX;
    std::vector<std::vector<double>> X_TWXinv;
    std::vector<std::vector<double>> X_TWY;


    for (int summaries_i = 0 ; summaries_i < gparam->NusedSummaries+1; summaries_i++)
    {
        std::vector<double> X_TWi;
        for (int simu_i = 0 ; simu_i < gparam->threshold; simu_i++)
        {
            X_TWi.push_back(X_TW[summaries_i][simu_i] * W[simu_i]);
        }
        X_TW.push_back(X_TWi);
    }


    for (int summaries_i = 0 ; summaries_i < gparam->NusedSummaries+1; summaries_i++)
    {
        std::vector<double>X_TWXi;
        for (int summaries_j = 0 ; summaries_j < gparam->NusedSummaries+1; summaries_j++)
        {
            double sum  = 0.0;
            for (int simu_i = 0 ; simu_i < gparam->threshold; simu_i++)
            {
                sum += X_TW[summaries_i][simu_i] * X[simu_i][summaries_j];
            }
            X_TWXi.push_back(sum);
        }
        X_TWX.push_back(X_TWXi);

    }


    LinAlg::Gauss(X_TWX, gparam->NusedSummaries+1, X_TWXinv);




    for (int summaries_i = 0 ; summaries_i < gparam->NusedSummaries+1; summaries_i++)
    {
        std::vector<double>X_TWYi;
        for (int param_i = 0 ; param_i < gparam->NusedParam; param_i++)
        {
            double sum  = 0.0;
            for (int simu_i = 0 ; simu_i < gparam->threshold; simu_i++)
            {
                sum += X_TW[summaries_i][simu_i] * Y[simu_i][param_i];
            }
            X_TWYi.push_back(sum);
        }
        X_TWY.push_back(X_TWYi);

    }


    for (int summaries_i = 0 ; summaries_i < gparam->NusedSummaries+1; summaries_i++)
    {
        for (int summaries_j = 0 ; summaries_j < gparam->NusedSummaries+1; summaries_j++)
        {
            for (int param_i = 0 ; param_i < gparam->NusedParam; param_i++)
            {

                B_HAT[summaries_i][param_i]+= X_TWXinv[summaries_i][summaries_j] * X_TWY[summaries_j][param_i];

            }
        }
    }






}


void Regression::ComputeX_TWY()
{




}
