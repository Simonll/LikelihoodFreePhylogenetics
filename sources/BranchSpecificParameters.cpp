/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/
#include "BranchSpecificParameters.h"

BranchSpecificParameters::BranchSpecificParameters()
{
    nucrrnr = new double*[Nnucp];
    gtnr = new double*[Nnucp];
    for (int i = 0; i < Nnucp; i++)
    {
        gtnr[i] = new double[Nnucp];
        nucrrnr[i] = new double[Nnucp];
    }
    for (int i = 0; i < Nnucp; i++)
    {
        for (int j = 0; j < Nnucp; j++)
        {
            gtnr[i][j] = 0.0;
            nucrrnr[i][j] = 0.0;
        }
    }

    nucp = new double [Nnucp];
    for (int i = 0; i < Nnucp; i++)
    {
        nucp[i] = 0.0;
    }

    nucrr = new double [Nnucrr];
    for (int i= 0; i < Nnucrr; i++)
    {
        nucrr[i] = 0.0;
    }


}

BranchSpecificParameters::~BranchSpecificParameters()
{
    //dtor
}

void BranchSpecificParameters::SetLocalParaemters(double* nucp, double* nucrr)
{

    for (int i = 0; i < Nnucp; i++)
    {
        this->nucp[i] = nucp[i];
    }


    for (int i= 0; i < Nnucrr; i++)
    {
        this->nucrr[i] = nucrr[i];
    }

}
