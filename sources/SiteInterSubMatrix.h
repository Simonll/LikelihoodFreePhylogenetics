/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef SITEINTERSUBMATRIX_H
#define SITEINTERSUBMATRIX_H

#include <iostream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

#define isnan std::isnan
#define isinf std::isinf
#define string std::string
#define ostream std::ostream
#define ofstream std::ofstream
#define istream std::istream
#define ifstream std::ifstream
#define cin std::cin
#define cerr std::cerr
#define cout std::cout
#define setw std::setw
#define ostringstream std::ostringstream
#define istringstream std::istringstream
#define IOS_APPEND std::ios_base::app
#define APPEND std::ios_base::app
#define OUT std::ios_base::out

#include "Random.h"
#include "Tree.h"
#include "SequenceAlignment.h"
#include "CodonStateSpace.h"
#include "CodonSequenceAlignment.h"
#include "LocalParameters.h"
#include "EvolHistStatistics.h"


class SiteInterSubMatrix
{
public:

    // SiteInterSubMatrix
    double*** submatrixTreeSim;
    double*** mutmatrixTreeSim;
    double*** selmatrixTreeSim;
    double* TotalSubRate;
    double* TotalMutRate;
    double* TotalSubRateNonSyn;
    double* TotalMutRateNonSyn;
    double* TotalSubRateSyn;
    double* TotalMutRateSyn;

    double* PartialSubRate;
    double* PartialMutRate;
    double* PartialSubRateNonSyn;
    double* PartialMutRateNonSyn;
    double* PartialSubRateSyn;
    double* PartialMutRateSyn;


    //parameters
    LocalParameters* lparam;

    //Constructor
    SiteInterSubMatrix(LocalParameters* lparam);
    SiteInterSubMatrix(LocalParameters* lparam, string seq);
    virtual ~SiteInterSubMatrix();


    //Getters


    double GetSubRate(int NodeIndex, int site_codon);
    double GetMutRate(int NodeIndex, int site_codon);

    double GetSubRateNonSyn(int NodeIndex, int site_codon,int**CurrentNodeNucSequence);
    double GetMutRateNonSyn(int NodeIndex, int site_codon,int**CurrentNodeNucSequence);

    double GetSubRateSyn(int NodeIndex, int site_codon,int**CurrentNodeNucSequence);
    double GetMutRateSyn(int NodeIndex, int site_codon,int**CurrentNodeNucSequence);

    double GetMutRateCpG(int NodeIndex, int**CurrentNodeNucSequence);
    double GetSubRateCpG(int NodeIndex, int**CurrentNodeNucSequence);


    //Setters
    void transfertTotalRate(int sourceNodeIndex, int sinkNodeIndex);
    void findCodonContext(int NodeIndex, int site_nuc,int nucFrom, int nucTo, int &pos1From, int &pos2From, int &pos3From, int &pos1To, int &pos2To, int &pos3To,int** CurrentNodeNucSequence);
    void ComputePartialRates(int NodeIndex, int site_codon,int** CurrentNodeNucSequence);
    void UpdateSubMatrixTreeSim(int NnodeIndex, int site_codon,int**CurrentNodeNucSequence);
    void UpdateSubMatrixSeq(int NnodeIndex,int**CurrentNodeNucSequence);
    int testGpTcontext(int NnodeIndex, int site, int nucFrom, int nucTo,int**CurrentNodeNucSequence);
    int testCpGcontext(int NnodeIndex, int site, int nucFrom, int nucTo,int**CurrentNodeNucSequence);
    int testTpAcontext(int NnodeIndex, int site, int nucFrom, int nucTo,int**CurrentNodeNucSequence);
    int testGCPref(int innucFrom, int innucTo);
    int testContextDinuc(int NodeIndex, int site_nuc, int* context, int nucTo, int**CurrentNodeNucSequence);
    void resetSubMatrix();
    void resetSubMatrixSeq();
    void transfertNodeMatrix(int sourceNodeIndex, int sinkNodeIndex, int site_nuc);
    
    //writers
    void WriteSubMatrix(ostream &mutation_os, ostream &selection_os, int nucTo);


    ////
    // Getters TotalRates
    ////

    double GetTotalMutRate(int NodeIndex)
    {

        return TotalMutRate[NodeIndex];
    }

    double GetTotalSubRate(int NodeIndex)
    {

        return TotalSubRate[NodeIndex];
    }


    double GetTotalMutRateNonSyn(int NodeIndex)
    {

        return TotalMutRateNonSyn[NodeIndex];
    }

    double GetTotalSubRateNonSyn(int NodeIndex)
    {

        return TotalSubRateNonSyn[NodeIndex];
    }

    double GetTotalMutRateSyn(int NodeIndex)
    {

        return TotalMutRateSyn[NodeIndex];
    }

    double GetTotalSubRateSyn(int NodeIndex)
    {

        return TotalSubRateSyn[NodeIndex];
    }

    ////
    // Getters PartialRates
    ////

    double GetPartialMutRate(int NodeIndex)
    {

        return  PartialMutRate[NodeIndex];
    }

    double GetPartialSubRate(int NodeIndex)
    {

        return  PartialSubRate[NodeIndex];
    }


    double GetPartialMutRateNonSyn(int NodeIndex)
    {

        return  PartialMutRateNonSyn[NodeIndex];
    }

    double GetPartialSubRateNonSyn(int NodeIndex)
    {

        return  PartialSubRateNonSyn[NodeIndex];
    }

    double GetPartialMutRateSyn(int NodeIndex)
    {

        return  PartialMutRateSyn[NodeIndex];
    }

    double GetPartialSubRateSyn(int NodeIndex)
    {

        return  PartialSubRateSyn[NodeIndex];
    }



    double GetSubRate(int NodeIndex, int site_nuc, int nucTo)
    {
        return submatrixTreeSim[NodeIndex][site_nuc][nucTo];

    }

    double GetMutRate(int NodeIndex, int site_nuc, int nucTo)
    {
        return mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];

    }

protected:
private:
};

#endif // SITEINTERSUBMATRIX_H
