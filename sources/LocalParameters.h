#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

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


#include "Random.h"
#include "Tree.h"
#include "SequenceAlignment.h"
#include "CodonStateSpace.h"
#include "CodonSequenceAlignment.h"
#include "StringStreamUtils.h"
#include "GlobalParameters.h"
#include "BranchSpecificParameters.h"


class LocalParameters
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
    int NEvoStats;
    int NSiteSpecificEvoStats;

    string* listParam;
    string* listSummaries;
    string* listEvoStats;
    string* listSiteSpecificEvoStats;
    string* listSpecies;

    std::map<string,int> mapUsedParam;
    std::map<string,int> mapUsedSummaries;
    std::map<string,int> mapUsedAccessorySummaries;
    std::map<string,int> mapUsedEvoStats;
    std::map<string,int> mapUsedEvoAncStats;
    std::map<string,int> mapUsedSiteSpecificEvoStats;

    int NusedEvoStats;
    int NusedSiteSpecificEvoStats;
    int NusedEvoAncStats;
    int NusedParam;
    int NusedSummaries;
    int NusedAccessorySummaries;
    int Ngenes;

    string data, localcontrolfile, output, model;



    BranchSpecificParameters** bparam;


    FileSequenceAlignment* dnadata;
    CodonStateSpace* codonstatespace;
    CodonSequenceAlignment* codondata;
    Tree* refTree;
    const TaxonSet*  taxonset;
    Link* outgroupLink;
    Link* newlink;
    Link* newnext;
    Branch* branchToInGroup;
    Branch* branchToOutGroup;
    Node* newnode;
    Random* rnd;

    bool iscodon, isdata;
    int optCpG, opt, optTpA;
    std::map <int,int> gtrMap;
    int gtr1NodeIndex;
    int gtr2NodeIndex;
    double rootlength, percentFromOutGroup, branchLengthBetweenInAndOutGroup ;


    struct params
    {


    };

    int Nsite_codon, Nsite_nuc, Ntaxa, Nstate_codon, startPoint, endPoint, everyPoint, tofasta, Nrep;
    string controlfile, code, chain, posteriorfile, taxa_a, taxa_b,taxa_gtr1_a, taxa_gtr1_b, taxa_gtr1_c, taxa_gtr2_a, taxa_gtr2_b,taxa_gtr2_c, distance, transformation;

    bool getrate, getrate1, getrate2;

    double omega;
    double *nucp, *nucrr, *nucp1, *nucrr1, *nucp2, *nucrr2; // 6 param
    double **nucrrnr, **nucrrnr1, **nucrrnr2; //12 param
    double **gtnr, **gtnr1, **gtnr2;
    double **ssaaprofiles;
    double *codonprofile;
    int *alloc;
    double lambda_TBL, lambda_omega, lambda_CpG, lambdaTG, lambdaCA, lambda_TpA, MutationNormFactor, MutationNormFactor1, MutationNormFactor2;
    double* muBranch;
    int fixNsite, fixomega, fixlambda_omega, fixlambda_TBL, fixlambda_CpG, fixlambda_TpA, fixgtr, fixgtr1, fixgtr2, fixgtnr, fixstat, fixts, fixtr, fixrr, fixkappa, fixhky, randomseed, verbose, rooted, fixroot, fixss;
    int MCMCpointID;

    std::vector<double> summariesRealData;
    std::vector<double> accessorysummariesRealData;
    std::vector<double> summariesSimulatedData;
    std::vector<double> accessorysummariesSimulatedData;
    std::vector<double> evostats;
    std::vector<double> ancevostats;
    std::vector<double> sitespecificevostats;
    std::vector<double> weights;

    //Constructor
    LocalParameters(GlobalParameters* gparam);
    void initContainers();
    virtual ~LocalParameters();

    // Writers
    void writeRealDataSummaries(ofstream&os,bool headers= true);
    void writeParam(ofstream& os);



    void tobstats(ofstream& os);
    void tobstats(ofstream& os,const Link* from);
    void toSsstats(ofstream& os);
    void toFasta(ofstream &os, int** currentNodeleafCodonSequence);
    void toAli(ofstream &os, int** currentNodeleafCodonSequence);

    // Readers
    void readChainCodonMutSelSBDP(int pt_i);
    void readChainCodonMutSelFinite(int pt_i);
    void readChainCodonMutSelSBDP();
    void readChainCodonMutSelFinite();
    void readFMutSelCodeML();
    void readLocalInstructions();


    //Setters
    void SetCurrentParametersFromPosterior(std::vector<std::vector<double>>posterior,int it);
    void SetTree();
    void SetBranchesLengthsBetweenInAndOutGroup();
    void SetRootBetweenInAndOutGroup();
    void SetRootLCA();
    void SetTreeStuff();

    void SetTreeStuffRecursively(Link* from, int notNodeIndex, int gtrIndex);

    //Getters
    int GetPointID();
    std::vector<double> GetCurrentParameters();
    std::vector<double> GetCurrentSummaries();
    std::vector<double> GetCurrentAccessorySummaries();
    std::vector<double> GetCurrentEvoStats();
    std::vector<double> GetCurrentAncEvoStats();
    std::vector<double> GetCurrentSiteSpecificEvoStats();
    std::vector<double> GetCurrentDistances();
    std::vector<double> GetCurrentWeights();


    void GetGTR1();
    void GetGTR2();


    void incrementStartPoint()
    {
        this->startPoint++;

    }

    int GetNsiteCodon()
    {
        return Nsite_codon;

    }

    double GetGTNR(int i, int j)
    {
        double norm = GetRateGTNR();
        return ((this->nucp[j] * this->nucrrnr[i][j])/norm);
    }

    double GetGTR(int i, int j )
    {

        return ((this->nucp[j] * this->nucrrnr[i][j]) / GetRate());

    }

    double GetGTRCodeML(int i, int j )
    {

        return (this->nucp[j] * this->nucrrnr[i][j]);

    }

    int    GetNucRRIndex(int i, int j)
    {
        return (i<j) ? (2 * 4 - i - 1) * i / 2 + j - i - 1 : (2 * 4 - j - 1) * j / 2 + i - j - 1 ;
    }

    double GetRateGTNR()
    {
        double norm = 0.0;
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
            {
                if (i!=j)
                {
                    norm += this->nucp[i] * this->nucp[j] * this->nucrrnr[i][j];
                }
            }
        }
        return norm * 3;
    }

    double GetRate()
    {
        if (getrate)
        {
            return MutationNormFactor;
        }
        double norm = 0.0;
        for (int i=0; i<4-1; i++)
        {
            for (int j=i+1; j<4; j++)
            {
                //norm += nucp[i] * nucp[j] * nucrr[GetNucRRIndex(i,j)];
                norm += this->nucp[i] * this->nucp[j] * this->nucrrnr[i][j];
                //cerr << this->nucp[i] << " " << this->nucp[j] << " " << this->nucrrnr[i][j] << " " << norm << "\n";
            }
        }
        // 2 for the symetry of the matrix??, and 3 for the number of codon positons
        getrate = true;
        MutationNormFactor = 2 * (norm * 3);
        return   MutationNormFactor;
    }



protected:

private:
};

#endif // PARAMETERS_H
