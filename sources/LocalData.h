#ifndef DATASTATS_H
#define DATASTATS_H

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
#include <sstream>
#include <map>

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

#include "Tree.h"
#include "GlobalParameters.h"
#include "SequenceAlignment.h"
#include "CodonStateSpace.h"
#include "CodonSequenceAlignment.h"
#include "StringStreamUtils.h"


class  LocalData
{
public:
    static const int Nnucp = 4;
    static const int Nnucrr = 6;
    static const int Ndinuc = 16;
    static const int Nstate_aa = 20;

    //to be set by global param
    double TOOSMALL;
    double TOOLARGE;
    double TOOLARGENEGATIVE;

    int NSummaries;
    int NParam;
    int NMapStats;

    string* listParam;
    string* listSummaries;
    string* listMapStats;

    std::map<string,int> mapUsedParam;
    std::map<string,int> mapUsedSummaries;
    std::map<string,int> mapUsedMapStats;
    std::map<string,int> mapUsedMapAncStats;

    int NusedMapStats;
    int NusedMapAncStats;
    int NusedParam;
    int NusedSummaries;
    int Ngenes;

    std::vector<string> listGenes;

    string localcontrolfile, output, model;
    std::vector<double> summariesRealData;


    void writeRealDataSummaries(ofstream&os,bool headers= true);
    void toFasta(ofstream &os, int** currentNodeleafCodonSequence);
    void toAli(ofstream &os, int** currentNodeleafCodonSequence);

    //GlobalParameters* gparam;
    string data;
    string tree;
    FileSequenceAlignment* dnadata;
    CodonStateSpace* codonstatespace;
    CodonSequenceAlignment* codondata;
    Tree* refTree;
    const TaxonSet*  taxonset;
    bool iscodon;
    string code;
    int Nsite_codon, Nsite_nuc, Ntaxa, Nnode, Nstate_codon;

    void readLocalData(int k);
    LocalData(GlobalParameters* gparam);
    virtual ~ LocalData();

protected:

private:
};

#endif // DATASTATS_H
