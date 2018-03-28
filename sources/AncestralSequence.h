/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ANCESTRALSEQUENCE_H
#define ANCESTRALSEQUENCE_H


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



class AncestralSequence
{
public:

    double** CurrentStationaryCodonSequence;
    int*  CurrentCodonSequence;
    int*  CurrentNucSequence;
    int*  CurrentAncestralCodonSequence;
    int*  CurrentAncestralNucSequence;

    LocalParameters* lparam;






    AncestralSequence(LocalParameters* lparam);
    ~AncestralSequence();

    void WriteStationary();
    void GetNewStationaryCodonSequence();
    int GetCurrentAncestralCodonSequence(int site_codon);
    int GetCurrentAncestralNucSequence(int site_nuc);

protected:

private:
};

#endif // ANCESTRALSEQUENCE_H
