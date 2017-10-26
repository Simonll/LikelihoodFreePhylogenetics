#ifndef SUMMARYSTATISTICS_H
#define SUMMARYSTATISTICS_H


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



#include "Tree.h"
#include "CodonStateSpace.h"
#include "CodonSequenceAlignment.h"
#include "SequenceAlignment.h"
#include "Random.h"
#include "StringStreamUtils.h"
#include "LocalParameters.h"
#include "LocalData.h"

class SummaryStatistics
{
public:
    //Parameters
    LocalParameters* lparam;
    LocalData* ldata;
    typedef double (SummaryStatistics::*funcpt)(CodonSequenceAlignment* codondata);
    std::map<string,funcpt> GetSummariesMap;



    //Summaries
    double** dinuc_usage;
    double** dinuc12_usage;
    double** dinuc23_usage;
    double** dinuc31_usage;
    double** dicodon_usage;
    double** diaa_usage;
    double* codon_usage;
    double* fcodon_usage;
    double* codon_usage_wonR;
    double* aa_usage;
    double* aa_usage_wonR;
    double* nuc_usage;
    double* nuc1_usage;
    double* nuc2_usage;
    double* nuc3_usage;
    double* nuc_meandiff;
    double* nuc1_meandiff;
    double* nuc2_meandiff;
    double* nuc3_meandiff;
    double* codon_meandiff;
    double* aa_meandiff;
    double* aa_wonR_meandiff;
    double* CGNAGR;
    int* nuc_pairwise;
    int* nuc1_pairwise;
    int* nuc2_pairwise;
    int* nuc3_pairwise;
    int* aa_pairwise;
    int* dinucCpG_pairwise;

    double nuc_site_comphet;
    double nuc1_site_comphet;
    double nuc2_site_comphet;
    double nuc3_site_comphet;
    double nuc_taxa_comphet;
    double nuc1_taxa_comphet;
    double nuc2_taxa_comphet;
    double nuc3_taxa_comphet;
    double codon_site_comphet;
    double codon_taxa_comphet;
    double aa_site_comphet;
    double aa_taxa_comphet;
    double RSCUentropy;


    bool codon_bool;
    bool codon_wonR_bool;
    bool dinuc_bool;
    bool dinuc12_bool;
    bool dinuc23_bool;
    bool dinuc31_bool;
    bool aa_bool;
    bool aa_wonR_bool;
    bool dicodon_bool;
    bool diaa_bool;
    bool nuc_bool;
    bool nuc1_bool;
    bool nuc2_bool;
    bool nuc3_bool;
    bool nuc_meandiff_bool;
    bool nuc1_meandiff_bool;
    bool nuc2_meandiff_bool;
    bool nuc3_meandiff_bool;
    bool codon_meandiff_bool;
    bool aa_meandiff_bool;
    bool aa_wonR_meandiff_bool;
    bool CGNAGR_bool;
    bool nuc_pairwise_bool;
    bool nuc1_pairwise_bool;
    bool nuc2_pairwise_bool;
    bool nuc3_pairwise_bool;
    bool aa_pairwise_bool;
    bool dinucCpG_pairwise_bool;
    bool nuc_site_comphet_bool;
    bool nuc1_site_comphet_bool;
    bool nuc2_site_comphet_bool;
    bool nuc3_site_comphet_bool;
    bool nuc_taxa_comphet_bool;
    bool nuc1_taxa_comphet_bool;
    bool nuc2_taxa_comphet_bool;
    bool nuc3_taxa_comphet_bool;
    bool codon_site_comphet_bool;
    bool codon_taxa_comphet_bool;
    bool aa_site_comphet_bool;
    bool aa_taxa_comphet_bool;
    bool RSCUentropy_bool;


    //Constructors
    SummaryStatistics(LocalParameters *lparam);
    SummaryStatistics(LocalData *ldata);


    //Destructors
    virtual ~SummaryStatistics();



    void MapFunctions();
    void computeSummariesFromData();
    void computeSummaries();
    void computeSummaries(int** CurrentNodeLeafCodonSequence);
    void computeSummariesAncestralSequence(int** CurrentAncestralCodonSequence);
    double SummaryStatistics::transformSummaryStatistics(double s);
    void GetRealStat();
//        std::vector<double> ReadRealStat(string inrealstat) {
//            std::vector<double> cur_real_statistics;
//            ostringstream buffer;
//            buffer << inrealstat;
//            ifstream is(buffer.str());
//            if (!is)       {
//                cerr << "error: did not find " << buffer.str() << "\n";
//                exit(1);
//            }
//
//            ReadLine(is); // Read headers
//            double cur_dl;
//            while(is.eof()) {
//                is >> cur_dl;
//                cur_real_statistics.push_back(cur_dl);
//
//            }
//
//            is.close();
//            return cur_real_statistics;
//
//        }





protected:

private:
    /////////////////
    // compositional heterogenity
    /////////////////
    double Getnuc_site_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc_site_comphet_bool)
        {
            nuc_site_comphet = codondata->nuc_site_comphet();
            nuc_site_comphet_bool = true;
        }
        return nuc_site_comphet;
    }
    double Getnuc1_site_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_site_comphet_bool)
        {
            nuc1_site_comphet = codondata->nuc1_site_comphet();
            nuc1_site_comphet_bool = true;
        }
        return nuc1_site_comphet;
    }
    double Getnuc2_site_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_site_comphet_bool)
        {
            nuc2_site_comphet = codondata->nuc2_site_comphet();
            nuc2_site_comphet_bool = true;
        }
        return nuc2_site_comphet;
    }
    double Getnuc3_site_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_site_comphet_bool)
        {
            nuc3_site_comphet = codondata->nuc3_site_comphet();
            nuc3_site_comphet_bool = true;
        }
        return nuc3_site_comphet;
    }
    double Getnuc_taxa_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc_taxa_comphet_bool)
        {
            nuc_taxa_comphet = codondata->nuc_taxa_comphet();
            nuc_taxa_comphet_bool = true;
        }
        return nuc_taxa_comphet;
    }
    double Getnuc1_taxa_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_taxa_comphet_bool)
        {
            nuc1_taxa_comphet = codondata->nuc1_taxa_comphet();
            nuc1_taxa_comphet_bool = true;
        }
        return nuc1_taxa_comphet;
    }
    double Getnuc2_taxa_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_taxa_comphet_bool)
        {
            nuc2_taxa_comphet = codondata->nuc2_taxa_comphet();
            nuc2_taxa_comphet_bool = true;
        }
        return nuc2_taxa_comphet;
    }
    double Getnuc3_taxa_comphet(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_taxa_comphet_bool)
        {
            nuc3_taxa_comphet = codondata->nuc3_taxa_comphet();
            nuc3_taxa_comphet_bool = true;
        }
        return nuc3_taxa_comphet;
    }
    double Getcodon_site_comphet(CodonSequenceAlignment* codondata)
    {
        if(!codon_site_comphet_bool)
        {
            codon_site_comphet = codondata->codon_site_comphet();
            codon_site_comphet_bool = true;
        }
        return codon_site_comphet;
    }
    double Getcodon_taxa_comphet(CodonSequenceAlignment* codondata)
    {
        if(!codon_taxa_comphet_bool)
        {
            codon_taxa_comphet = codondata->codon_taxa_comphet();
            codon_taxa_comphet_bool = true;
        }
        return codon_taxa_comphet;
    }
    double Getaa_site_comphet(CodonSequenceAlignment* codondata)
    {
        if(!aa_site_comphet_bool)
        {
            aa_site_comphet = codondata->aa_site_comphet();
            aa_site_comphet_bool = true;
        }
        return aa_site_comphet;
    }
    double Getaa_taxa_comphet(CodonSequenceAlignment* codondata)
    {
        if(!aa_taxa_comphet_bool)
        {
            aa_taxa_comphet = codondata->aa_taxa_comphet();
            aa_taxa_comphet_bool = true;
        }
        return aa_taxa_comphet;
    }


    /////////////////
    // codon_usage
    /////////////////

    double GetGGG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGG")];
    }
    double GetGGA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGA")];
    }
    double GetGGC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGC")];
    }
    double GetGGT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGT")];
    }
    double GetGAG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAG")];
    }
    double GetGAA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAA")];
    }
    double GetGAC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAC")];
    }
    double GetGAT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAT")];
    }
    double GetGCG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCG")];
    }
    double GetGCA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCA")];
    }
    double GetGCC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCC")];
    }
    double GetGCT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCT")];
    }
    double GetGTG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTG")];
    }
    double GetGTA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTA")];
    }
    double GetGTC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTC")];
    }
    double GetGTT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTT")];
    }
    double GetAGG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGG")];
    }
    double GetAGA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGA")];
    }
    double GetAGC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGC")];
    }
    double GetAGT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGT")];
    }
    double GetAAG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAG")];
    }
    double GetAAA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAA")];
    }
    double GetAAC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAC")];
    }
    double GetAAT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAT")];
    }
    double GetACG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACG")];
    }
    double GetACA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACA")];
    }
    double GetACC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACC")];
    }
    double GetACT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACT")];
    }
    double GetATG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATG")];
    }
    double GetATA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATA")];
    }
    double GetATC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATC")];
    }
    double GetATT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATT")];
    }
    double GetCGG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGG")];
    }
    double GetCGA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGA")];
    }
    double GetCGC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGC")];
    }
    double GetCGT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGT")];
    }
    double GetCAG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAG")];
    }
    double GetCAA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAA")];
    }
    double GetCAC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAC")];
    }
    double GetCAT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAT")];
    }
    double GetCCG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCG")];
    }
    double GetCCA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCA")];
    }
    double GetCCC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCC")];
    }
    double GetCCT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCT")];
    }
    double GetCTG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTG")];
    }
    double GetCTA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTA")];
    }
    double GetCTC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTC")];
    }
    double GetCTT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTT")];
    }
    double GetTGG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGG")];
    }
    double GetTGA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGA")];
    }
    double GetTGC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGC")];
    }
    double GetTGT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGT")];
    }
    double GetTAG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAG")];
    }

    double GetTAA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAA")];
    }

    double GetTAC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAC")];
    }

    double GetTAT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAT")];
    }
    double GetTCG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCG")];
    }

    double GetTCA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCA")];
    }
    double GetTCC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCC")];
    }
    double GetTCT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCT")];
    }

    double GetTTG(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTG")];
    }

    double GetTTA(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTA")];
    }
    double GetTTC(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTC")];
    }
    double GetTTT(CodonSequenceAlignment* codondata)
    {
        if(!codon_bool)
        {
            codondata->CodonUsagePerAA(codon_usage);
            codon_bool = true;
        }
        return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTT")];
    }



/////////////////
// codon_usage_wonR
/////////////////

    double GetGGGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GGG")];
    }
    double GetGGAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GGA")];
    }
    double GetGGCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GGC")];
    }
    double GetGGTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GGT")];
    }
    double GetGAGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GAG")];
    }
    double GetGAAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GAA")];
    }
    double GetGACwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GAC")];
    }
    double GetGATwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GAT")];
    }
    double GetGCGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GCG")];
    }
    double GetGCAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GCA")];
    }
    double GetGCCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GCC")];
    }
    double GetGCTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GCT")];
    }
    double GetGTGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GTG")];
    }
    double GetGTAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GTA")];
    }
    double GetGTCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GTC")];
    }
    double GetGTTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("GTT")];
    }
    double GetAGGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AGG")];
    }
    double GetAGAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AGA")];
    }
    double GetAGCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AGC")];
    }
    double GetAGTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AGT")];
    }
    double GetAAGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AAG")];
    }
    double GetAAAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AAA")];
    }
    double GetAACwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AAC")];
    }
    double GetAATwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("AAT")];
    }
    double GetACGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ACG")];
    }
    double GetACAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ACA")];
    }
    double GetACCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ACC")];
    }
    double GetACTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ACT")];
    }
    double GetATGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ATG")];
    }
    double GetATAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ATA")];
    }
    double GetATCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ATC")];
    }
    double GetATTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("ATT")];
    }
    double GetCGGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CGG")];
    }
    double GetCGAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CGA")];
    }
    double GetCGCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CGC")];
    }
    double GetCGTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CGT")];
    }
    double GetCAGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CAG")];
    }
    double GetCAAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CAA")];
    }
    double GetCACwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CAC")];
    }
    double GetCATwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CAT")];
    }
    double GetCCGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CCG")];
    }
    double GetCCAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CCA")];
    }
    double GetCCCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CCC")];
    }
    double GetCCTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CCT")];
    }
    double GetCTGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CTG")];
    }
    double GetCTAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CTA")];
    }
    double GetCTCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CTC")];
    }
    double GetCTTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("CTT")];
    }
    double GetTGGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TGG")];
    }
    double GetTGAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TGA")];
    }
    double GetTGCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TGC")];
    }
    double GetTGTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TGT")];
    }
    double GetTAGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TAG")];
    }

    double GetTAAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TAA")];
    }

    double GetTACwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TAC")];
    }

    double GetTATwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TAT")];
    }
    double GetTCGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TCG")];
    }

    double GetTCAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TCA")];
    }
    double GetTCCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TCC")];
    }
    double GetTCTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TCT")];
    }

    double GetTTGwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TTG")];
    }

    double GetTTAwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TTA")];
    }
    double GetTTCwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TTC")];
    }
    double GetTTTwonR(CodonSequenceAlignment* codondata)
    {
        if(!codon_wonR_bool)
        {
            codondata->CodonUsagePerAAwonR(codon_usage_wonR);
            codon_wonR_bool = true;
        }
        return (double) codon_usage_wonR[codondata->GetCodonStateSpace()->GetState("TTT")];
    }




    /////////////////
    // aa_usage
    /////////////////

    double GetA(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[0];
    }
    double GetC(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[1];
    }
    double GetD(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[2];
    }
    double GetE(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[3];
    }
    double GetF(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[4];
    }
    double GetG(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[5];
    }
    double GetH(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[6];
    }
    double GetI(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[7];
    }
    double GetK(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[8];
    }
    double GetL(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[9];
    }
    double GetM(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[10];
    }
    double GetN(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[11];
    }
    double GetP(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[12];
    }
    double GetQ(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[13];
    }
    double GetR(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[14];
    }
    double GetS(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[15];
    }
    double GetT(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[16];
    }
    double GetV(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[17];
    }
    double GetW(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[18];
    }
    double GetY(CodonSequenceAlignment* codondata)
    {
        if(!aa_bool)
        {
            codondata->aa_usage(aa_usage);
            aa_bool = true;
        }
        return (double) aa_usage[19];
    }



//        /////////////////
//        // aa_usagewo_nr
//        /////////////////

    double GetAwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[0];
    }
    double GetCwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[1];
    }
    double GetDwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[2];
    }
    double GetEwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[3];
    }
    double GetFwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[4];
    }
    double GetGwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[5];
    }
    double GetHwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[6];
    }
    double GetIwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[7];
    }
    double GetKwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[8];
    }
    double GetLwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[9];
    }
    double GetMwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[10];
    }
    double GetNwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[11];
    }
    double GetPwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[12];
    }
    double GetQwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[13];
    }
    double GetRwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[14];
    }
    double GetSwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[15];
    }
    double GetTwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[16];
    }
    double GetVwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[17];
    }
    double GetWwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[18];
    }
    double GetYwonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_wonR_bool)
        {
            codondata->aa_usage_wonr(aa_usage_wonR);
            aa_wonR_bool = true;
        }
        return (double) aa_usage_wonR[19];
    }



    /////////////////
    // dinuc31
    /////////////////
    double GetDinuc31AA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[0][0];
    }

    double GetDinuc31AC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[0][1];
    }

    double GetDinuc31AG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[0][2];

    }
    double GetDinuc31AT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[0][3];
    }
    double GetDinuc31CA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[1][0];
    }

    double GetDinuc31CC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[1][1];
    }
    double GetDinuc31CG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[1][2];

    }
    double GetDinuc31CT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[1][3];
    }
    double GetDinuc31GA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[2][0];
    }

    double GetDinuc31GC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[2][1];
    }
    double GetDinuc31GG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[2][2];

    }
    double GetDinuc31GT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[2][3];
    }
    double GetDinuc31TA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[3][0];
    }

    double GetDinuc31TC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[3][1];
    }
    double GetDinuc31TG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[3][2];

    }
    double GetDinuc31TT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc31_bool)
        {
            codondata->dinuc31_usage(dinuc31_usage);
            dinuc31_bool = true;
        }
        return (double) dinuc31_usage[3][3];
    }
    /////////////////
    // dinuc23
    /////////////////
    double GetDinuc23AA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[0][0];
    }

    double GetDinuc23AC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[0][1];
    }
    double GetDinuc23AG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[0][2];

    }
    double GetDinuc23AT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[0][3];
    }
    double GetDinuc23CA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[1][0];
    }

    double GetDinuc23CC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[1][1];
    }
    double GetDinuc23CG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[1][2];

    }
    double GetDinuc23CT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[1][3];
    }
    double GetDinuc23GA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[2][0];
    }

    double GetDinuc23GC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[2][1];
    }
    double GetDinuc23GG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[2][2];

    }
    double GetDinuc23GT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[2][3];
    }
    double GetDinuc23TA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[3][0];
    }

    double GetDinuc23TC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[3][1];
    }
    double GetDinuc23TG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[3][2];

    }
    double GetDinuc23TT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc23_bool)
        {
            codondata->dinuc23_usage(dinuc23_usage);
            dinuc23_bool = true;
        }
        return (double) dinuc23_usage[3][3];
    }

    /////////////////
    // dinuc12
    /////////////////
    double GetDinuc12AA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[0][0];
    }

    double GetDinuc12AC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[0][1];
    }
    double GetDinuc12AG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[0][2];

    }
    double GetDinuc12AT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[0][3];
    }
    double GetDinuc12CA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[1][0];
    }

    double GetDinuc12CC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[1][1];
    }
    double GetDinuc12CG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[1][2];

    }
    double GetDinuc12CT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[1][3];
    }
    double GetDinuc12GA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[2][0];
    }

    double GetDinuc12GC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[2][1];
    }
    double GetDinuc12GG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[2][2];

    }
    double GetDinuc12GT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[2][3];
    }
    double GetDinuc12TA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[3][0];
    }

    double GetDinuc12TC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[3][1];
    }
    double GetDinuc12TG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[3][2];

    }
    double GetDinuc12TT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc12_bool)
        {
            codondata->dinuc12_usage(dinuc12_usage);
            dinuc12_bool = true;
        }
        return (double) dinuc12_usage[3][3];
    }
    /////////////////
    // dinuc
    /////////////////
    double GetDinucAA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[0][0];
    }

    double GetDinucAC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[0][1];
    }
    double GetDinucAG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[0][2];

    }
    double GetDinucAT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[0][3];
    }
    double GetDinucCA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[1][0];
    }

    double GetDinucCC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[1][1];
    }
    double GetDinucCG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[1][2];

    }
    double GetDinucCT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[1][3];
    }
    double GetDinucGA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[2][0];
    }

    double GetDinucGC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[2][1];
    }
    double GetDinucGG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[2][2];

    }
    double GetDinucGT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[2][3];
    }
    double GetDinucTA(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[3][0];
    }

    double GetDinucTC(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[3][1];
    }
    double GetDinucTG(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[3][2];

    }
    double GetDinucTT(CodonSequenceAlignment* codondata)
    {
        if(!dinuc_bool)
        {
            codondata->dinuc_usage(dinuc_usage);
            dinuc_bool = true;
        }
        return (double) dinuc_usage[3][3];
    }

    /////////////////
    // nuc1
    /////////////////
    double GetNuc1A(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_bool)
        {
            codondata->nuc1_usage(nuc1_usage);
            nuc1_bool = true;
        }
        return (double) nuc1_usage[0];
    }

    double GetNuc1C(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_bool)
        {
            codondata->nuc1_usage(nuc1_usage);
            nuc1_bool = true;
        }
        return (double) nuc1_usage[1];
    }

    double GetNuc1G(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_bool)
        {
            codondata->nuc1_usage(nuc1_usage);
            nuc1_bool = true;
        }
        return (double) nuc1_usage[2];
    }

    double GetNuc1T(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_bool)
        {
            codondata->nuc1_usage(nuc1_usage);
            nuc1_bool = true;
        }
        return (double) nuc1_usage[3];
    }
    /////////////////
    // nuc2
    /////////////////
    double GetNuc2A(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_bool)
        {
            codondata->nuc2_usage(nuc2_usage);
            nuc2_bool = true;
        }
        return (double) nuc2_usage[0];
    }

    double GetNuc2C(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_bool)
        {
            codondata->nuc2_usage(nuc2_usage);
            nuc2_bool = true;
        }
        return (double) nuc2_usage[1];
    }

    double GetNuc2G(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_bool)
        {
            codondata->nuc2_usage(nuc2_usage);
            nuc2_bool = true;
        }
        return (double) nuc2_usage[2];
    }

    double GetNuc2T(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_bool)
        {
            codondata->nuc2_usage(nuc2_usage);
            nuc2_bool = true;
        }
        return (double) nuc2_usage[3];
    }
    /////////////////
    // nuc3
    /////////////////
    double GetNuc3A(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_bool)
        {
            codondata->nuc3_usage(nuc3_usage);
            nuc3_bool = true;
        }
        return (double) nuc3_usage[0];
    }

    double GetNuc3C(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_bool)
        {
            codondata->nuc3_usage(nuc3_usage);
            nuc3_bool = true;
        }
        return (double) nuc3_usage[1];
    }

    double GetNuc3G(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_bool)
        {
            codondata->nuc3_usage(nuc3_usage);
            nuc3_bool = true;
        }
        return (double) nuc3_usage[2];
    }

    double GetNuc3T(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_bool)
        {
            codondata->nuc3_usage(nuc3_usage);
            nuc3_bool = true;
        }
        return (double) nuc3_usage[3];
    }
    /////////////////
    // nuc
    /////////////////
    double GetNucA(CodonSequenceAlignment* codondata)
    {
        if(!nuc_bool)
        {
            codondata->nuc_usage(nuc_usage);
            nuc_bool = true;
        }
        return (double) nuc_usage[0];
    }

    double GetNucC(CodonSequenceAlignment* codondata)
    {
        if(!nuc_bool)
        {
            codondata->nuc_usage(nuc_usage);
            nuc_bool = true;
        }
        return (double) nuc_usage[1];
    }

    double GetNucG(CodonSequenceAlignment* codondata)
    {
        if(!nuc_bool)
        {
            codondata->nuc_usage(nuc_usage);
            nuc_bool = true;
        }
        return (double) nuc_usage[2];
    }

    double GetNucT(CodonSequenceAlignment* codondata)
    {
        if(!nuc_bool)
        {
            codondata->nuc_usage(nuc_usage);
            nuc_bool = true;
        }
        return (double) nuc_usage[3];
    }


    /////////////////
    // nuc3_meandiff
    /////////////////
    double GetNuc3mean(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_meandiff_bool)
        {
            codondata->nuc3_meandiff(nuc3_meandiff);
            nuc3_meandiff_bool = true;
        }
        return (double) nuc3_meandiff[0];
    }

    double GetNuc3var(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_meandiff_bool)
        {
            codondata->nuc3_meandiff(nuc3_meandiff);
            nuc3_meandiff_bool = true;
        }
        return (double) nuc3_meandiff[3];
    }
    /////////////////
    // nuc2_meandiff
    /////////////////
    double GetNuc2mean(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_meandiff_bool)
        {
            codondata->nuc2_meandiff(nuc2_meandiff);
            nuc2_meandiff_bool = true;
        }
        return (double) nuc2_meandiff[0];
    }

    double GetNuc2var(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_meandiff_bool)
        {
            codondata->nuc2_meandiff(nuc2_meandiff);
            nuc2_meandiff_bool = true;
        }
        return (double) nuc2_meandiff[2];
    }
    /////////////////
    // nuc1_meandiff
    /////////////////
    double GetNuc1mean(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_meandiff_bool)
        {
            codondata->nuc1_meandiff(nuc1_meandiff);
            nuc1_meandiff_bool = true;
        }
        return (double) nuc1_meandiff[0];
    }

    double GetNuc1var(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_meandiff_bool)
        {
            codondata->nuc1_meandiff(nuc1_meandiff);
            nuc1_meandiff_bool = true;
        }
        return (double) nuc1_meandiff[1];
    }
    /////////////////
    // nuc_meandiff
    /////////////////
    double GetNucmean(CodonSequenceAlignment* codondata)
    {
        if(!nuc_meandiff_bool)
        {
            codondata->nuc_meandiff(nuc_meandiff);
            nuc_meandiff_bool = true;
        }
        return (double) nuc_meandiff[0];
    }

    double GetNucvar(CodonSequenceAlignment* codondata)
    {
        if(!nuc_meandiff_bool)
        {
            codondata->nuc_meandiff(nuc_meandiff);
            nuc_meandiff_bool = true;
        }
        return (double) nuc_meandiff[1];
    }

    /////////////////
    // codon_meandiff
    /////////////////
    double GetCodonmean(CodonSequenceAlignment* codondata)
    {
        if(!codon_meandiff_bool)
        {
            codondata->codon_meandiff(codon_meandiff);
            codon_meandiff_bool = true;
        }
        return (double) codon_meandiff[0];
    }

    double GetCodonvar(CodonSequenceAlignment* codondata)
    {
        if(!codon_meandiff_bool)
        {
            codondata->codon_meandiff(codon_meandiff);
            codon_meandiff_bool = true;
        }
        return (double) codon_meandiff[1];
    }
    /////////////////
    // aa_meandiff
    /////////////////
    double GetAAmean(CodonSequenceAlignment* codondata)
    {
        if(!aa_meandiff_bool)
        {
            codondata->aa_meandiff(aa_meandiff);
            aa_meandiff_bool = true;
        }
        return (double) aa_meandiff[0];
    }

    double GetAAvar(CodonSequenceAlignment* codondata)
    {
        if(!aa_meandiff_bool)
        {
            codondata->aa_meandiff(aa_meandiff);
            aa_meandiff_bool = true;
        }
        return (double) aa_meandiff[1];
    }


    double GetAAmean_wonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_meandiff_bool)
        {
            codondata->aa_meandiff_wonr(aa_wonR_meandiff);
            aa_wonR_meandiff_bool = true;
        }
        return (double) aa_wonR_meandiff[0];
    }

    double GetAAvar_wonR(CodonSequenceAlignment* codondata)
    {
        if(!aa_meandiff_bool)
        {
            codondata->aa_meandiff_wonr(aa_wonR_meandiff);
            aa_wonR_meandiff_bool = true;
        }
        return (double) aa_wonR_meandiff[1];
    }

    /////////////////
    // dinucCpG_pairwise
    /////////////////
    double GetdinucCpG_TpG(CodonSequenceAlignment* codondata)
    {
        if(!dinucCpG_pairwise_bool)
        {
            codondata->dinucCpG_pairwise(dinucCpG_pairwise);
            dinucCpG_pairwise_bool = true;
        }
        return (double) dinucCpG_pairwise[0];
    }

    double GetdinucCpG_CpA(CodonSequenceAlignment* codondata)
    {
        if(!dinucCpG_pairwise_bool)
        {
            codondata->dinucCpG_pairwise(dinucCpG_pairwise);
            dinucCpG_pairwise_bool = true;
        }
        return (double) dinucCpG_pairwise[1];
    }

    double GetdinucApG_TpG(CodonSequenceAlignment* codondata)
    {
        if(!dinucCpG_pairwise_bool)
        {
            codondata->dinucCpG_pairwise(dinucCpG_pairwise);
            dinucCpG_pairwise_bool = true;
        }
        return (double) dinucCpG_pairwise[2];
    }


    /////////////////
    // aa_pairwise
    /////////////////
    double GetpwAA(CodonSequenceAlignment* codondata)
    {
        if(!aa_pairwise_bool)
        {
            codondata->aa_pairwise(aa_pairwise);
            aa_pairwise_bool = true;
        }
        double sum = 0;
        for (int i = 0 ; i < 190; i++)
        {
            sum += aa_pairwise[i];
        }
        return sum;
    }




    /////////////////
    // nuc3_pairwise
    /////////////////

    double Getpw3GT(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc3_pairwise);
            nuc3_pairwise_bool = true;
        }
        return (double) nuc3_pairwise[5];
    }

    double Getpw3CT(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc3_pairwise);
            nuc3_pairwise_bool = true;
        }
        return (double) nuc3_pairwise[4];
    }

    double Getpw3CG(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc3_pairwise);
            nuc3_pairwise_bool = true;
        }
        return (double) nuc3_pairwise[3];
    }


    double Getpw3AT(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc3_pairwise);
            nuc3_pairwise_bool = true;
        }
        return (double) nuc3_pairwise[2];
    }

    double Getpw3AG(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc3_pairwise);
            nuc3_pairwise_bool = true;
        }
        return (double) nuc3_pairwise[1];
    }


    double Getpw3AC(CodonSequenceAlignment* codondata)
    {
        if(!nuc3_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc3_pairwise);
            nuc3_pairwise_bool = true;
        }
        return (double) nuc3_pairwise[0];
    }

    /////////////////
    // nuc2_pairwise
    /////////////////

    double Getpw2GT(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc2_pairwise);
            nuc2_pairwise_bool = true;
        }
        return (double) nuc2_pairwise[5];
    }

    double Getpw2CT(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc2_pairwise);
            nuc2_pairwise_bool = true;
        }
        return (double) nuc2_pairwise[4];
    }

    double Getpw2CG(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc2_pairwise);
            nuc2_pairwise_bool = true;
        }
        return (double) nuc2_pairwise[3];
    }


    double Getpw2AT(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc2_pairwise);
            nuc2_pairwise_bool = true;
        }
        return (double) nuc2_pairwise[2];
    }

    double Getpw2AG(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc2_pairwise);
            nuc2_pairwise_bool = true;
        }
        return (double) nuc2_pairwise[1];
    }


    double Getpw2AC(CodonSequenceAlignment* codondata)
    {
        if(!nuc2_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc2_pairwise);
            nuc2_pairwise_bool = true;
        }
        return (double) nuc2_pairwise[0];
    }

    /////////////////
    // nuc1_pairwise
    /////////////////

    double Getpw1GT(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc1_pairwise);
            nuc1_pairwise_bool = true;
        }
        return (double) nuc1_pairwise[5];
    }

    double Getpw1CT(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc1_pairwise);
            nuc1_pairwise_bool = true;
        }
        return (double) nuc1_pairwise[4];
    }

    double Getpw1CG(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc1_pairwise);
            nuc1_pairwise_bool = true;
        }
        return (double) nuc1_pairwise[3];
    }


    double Getpw1AT(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc1_pairwise);
            nuc1_pairwise_bool = true;
        }
        return (double) nuc1_pairwise[2];
    }

    double Getpw1AG(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc1_pairwise);
            nuc1_pairwise_bool = true;
        }
        return (double) nuc1_pairwise[1];
    }


    double Getpw1AC(CodonSequenceAlignment* codondata)
    {
        if(!nuc1_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc1_pairwise);
            nuc1_pairwise_bool = true;
        }
        return (double) nuc1_pairwise[0];
    }


    /////////////////
    // nuc_pairwise
    /////////////////


    double GetpwGT(CodonSequenceAlignment* codondata)
    {
        if(!nuc_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc_pairwise);
            nuc_pairwise_bool = true;
        }
        return (double) nuc_pairwise[5];
    }

    double GetpwCT(CodonSequenceAlignment* codondata)
    {
        if(!nuc_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc_pairwise);
            nuc_pairwise_bool = true;
        }
        return (double) nuc_pairwise[4];
    }

    double GetpwCG(CodonSequenceAlignment* codondata)
    {
        if(!nuc_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc_pairwise);
            nuc_pairwise_bool = true;
        }
        return (double) nuc_pairwise[3];
    }


    double GetpwAT(CodonSequenceAlignment* codondata)
    {
        if(!nuc_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc_pairwise);
            nuc_pairwise_bool = true;
        }
        return (double) nuc_pairwise[2];
    }

    double GetpwAG(CodonSequenceAlignment* codondata)
    {
        if(!nuc_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc_pairwise);
            nuc_pairwise_bool = true;
        }
        return (double) nuc_pairwise[1];
    }


    double GetpwAC(CodonSequenceAlignment* codondata)
    {
        if(!nuc_pairwise_bool)
        {
            codondata->nuc_pairwise(nuc_pairwise);
            nuc_pairwise_bool = true;
        }
        return (double) nuc_pairwise[0];
    }


    double GetCGNAGRcgnconst(CodonSequenceAlignment* codondata)
    {
        if(!CGNAGR_bool)
        {
            codondata->CGNAGR(CGNAGR);
            CGNAGR_bool = true;
        }
        return (double) CGNAGR[0];

    }

    double GetCGNAGRagrconst(CodonSequenceAlignment* codondata)
    {
        if(!CGNAGR_bool)
        {
            codondata->CGNAGR(CGNAGR);
            CGNAGR_bool = true;
        }
        return (double) CGNAGR[1];

    }


    double GetCGNAGRcgnvar(CodonSequenceAlignment* codondata)
    {
        if(!CGNAGR_bool)
        {
            codondata->CGNAGR(CGNAGR);
            CGNAGR_bool = true;
        }
        return (double) CGNAGR[2];

    }

    double GetCGNAGRagrvar(CodonSequenceAlignment* codondata)
    {
        if(!CGNAGR_bool)
        {
            codondata->CGNAGR(CGNAGR);
            CGNAGR_bool = true;
        }
        return (double) CGNAGR[3];

    }


    double GetDIAA_AA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][0];
    }
    double GetDIAA_AC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][1];
    }
    double GetDIAA_AD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][2];
    }
    double GetDIAA_AE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][3];
    }
    double GetDIAA_AF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][4];
    }
    double GetDIAA_AG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][5];
    }
    double GetDIAA_AH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][6];
    }
    double GetDIAA_AI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][7];
    }
    double GetDIAA_AK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][8];
    }
    double GetDIAA_AL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][9];
    }
    double GetDIAA_AM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][10];
    }
    double GetDIAA_AN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][11];
    }
    double GetDIAA_AP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][12];
    }
    double GetDIAA_AQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][13];
    }
    double GetDIAA_AR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][14];
    }
    double GetDIAA_AS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][15];
    }
    double GetDIAA_AT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][16];
    }
    double GetDIAA_AV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][17];
    }
    double GetDIAA_AW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][18];
    }
    double GetDIAA_AY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[0][19];
    }
    double GetDIAA_CA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][0];
    }
    double GetDIAA_CC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][1];
    }
    double GetDIAA_CD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][2];
    }
    double GetDIAA_CE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][3];
    }
    double GetDIAA_CF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][4];
    }
    double GetDIAA_CG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][5];
    }
    double GetDIAA_CH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][6];
    }
    double GetDIAA_CI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][7];
    }
    double GetDIAA_CK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][8];
    }
    double GetDIAA_CL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][9];
    }
    double GetDIAA_CM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][10];
    }
    double GetDIAA_CN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][11];
    }
    double GetDIAA_CP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][12];
    }
    double GetDIAA_CQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][13];
    }
    double GetDIAA_CR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][14];
    }
    double GetDIAA_CS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][15];
    }
    double GetDIAA_CT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][16];
    }
    double GetDIAA_CV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][17];
    }
    double GetDIAA_CW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][18];
    }
    double GetDIAA_CY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[1][19];
    }
    double GetDIAA_DA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][0];
    }
    double GetDIAA_DC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][1];
    }
    double GetDIAA_DD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][2];
    }
    double GetDIAA_DE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][3];
    }
    double GetDIAA_DF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][4];
    }
    double GetDIAA_DG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][5];
    }
    double GetDIAA_DH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][6];
    }
    double GetDIAA_DI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][7];
    }
    double GetDIAA_DK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][8];
    }
    double GetDIAA_DL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][9];
    }
    double GetDIAA_DM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][10];
    }
    double GetDIAA_DN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][11];
    }
    double GetDIAA_DP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][12];
    }
    double GetDIAA_DQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][13];
    }
    double GetDIAA_DR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][14];
    }
    double GetDIAA_DS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][15];
    }
    double GetDIAA_DT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][16];
    }
    double GetDIAA_DV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][17];
    }
    double GetDIAA_DW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][18];
    }
    double GetDIAA_DY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[2][19];
    }
    double GetDIAA_EA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][0];
    }
    double GetDIAA_EC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][1];
    }
    double GetDIAA_ED(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][2];
    }
    double GetDIAA_EE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][3];
    }
    double GetDIAA_EF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][4];
    }
    double GetDIAA_EG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][5];
    }
    double GetDIAA_EH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][6];
    }
    double GetDIAA_EI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][7];
    }
    double GetDIAA_EK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][8];
    }
    double GetDIAA_EL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][9];
    }
    double GetDIAA_EM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][10];
    }
    double GetDIAA_EN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][11];
    }
    double GetDIAA_EP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][12];
    }
    double GetDIAA_EQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][13];
    }
    double GetDIAA_ER(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][14];
    }
    double GetDIAA_ES(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][15];
    }
    double GetDIAA_ET(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][16];
    }
    double GetDIAA_EV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][17];
    }
    double GetDIAA_EW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][18];
    }
    double GetDIAA_EY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[3][19];
    }
    double GetDIAA_FA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][0];
    }
    double GetDIAA_FC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][1];
    }
    double GetDIAA_FD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][2];
    }
    double GetDIAA_FE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][3];
    }
    double GetDIAA_FF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][4];
    }
    double GetDIAA_FG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][5];
    }
    double GetDIAA_FH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][6];
    }
    double GetDIAA_FI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][7];
    }
    double GetDIAA_FK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][8];
    }
    double GetDIAA_FL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][9];
    }
    double GetDIAA_FM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][10];
    }
    double GetDIAA_FN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][11];
    }
    double GetDIAA_FP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][12];
    }
    double GetDIAA_FQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][13];
    }
    double GetDIAA_FR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][14];
    }
    double GetDIAA_FS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][15];
    }
    double GetDIAA_FT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][16];
    }
    double GetDIAA_FV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][17];
    }
    double GetDIAA_FW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][18];
    }
    double GetDIAA_FY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[4][19];
    }
    double GetDIAA_GA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][0];
    }
    double GetDIAA_GC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][1];
    }
    double GetDIAA_GD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][2];
    }
    double GetDIAA_GE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][3];
    }
    double GetDIAA_GF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][4];
    }
    double GetDIAA_GG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][5];
    }
    double GetDIAA_GH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][6];
    }
    double GetDIAA_GI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][7];
    }
    double GetDIAA_GK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][8];
    }
    double GetDIAA_GL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][9];
    }
    double GetDIAA_GM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][10];
    }
    double GetDIAA_GN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][11];
    }
    double GetDIAA_GP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][12];
    }
    double GetDIAA_GQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][13];
    }
    double GetDIAA_GR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][14];
    }
    double GetDIAA_GS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][15];
    }
    double GetDIAA_GT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][16];
    }
    double GetDIAA_GV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][17];
    }
    double GetDIAA_GW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][18];
    }
    double GetDIAA_GY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[5][19];
    }
    double GetDIAA_HA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][0];
    }
    double GetDIAA_HC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][1];
    }
    double GetDIAA_HD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][2];
    }
    double GetDIAA_HE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][3];
    }
    double GetDIAA_HF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][4];
    }
    double GetDIAA_HG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][5];
    }
    double GetDIAA_HH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][6];
    }
    double GetDIAA_HI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][7];
    }
    double GetDIAA_HK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][8];
    }
    double GetDIAA_HL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][9];
    }
    double GetDIAA_HM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][10];
    }
    double GetDIAA_HN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][11];
    }
    double GetDIAA_HP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][12];
    }
    double GetDIAA_HQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][13];
    }
    double GetDIAA_HR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][14];
    }
    double GetDIAA_HS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][15];
    }
    double GetDIAA_HT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][16];
    }
    double GetDIAA_HV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][17];
    }
    double GetDIAA_HW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][18];
    }
    double GetDIAA_HY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[6][19];
    }
    double GetDIAA_IA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][0];
    }
    double GetDIAA_IC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][1];
    }
    double GetDIAA_ID(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][2];
    }
    double GetDIAA_IE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][3];
    }
    double GetDIAA_IF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][4];
    }
    double GetDIAA_IG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][5];
    }
    double GetDIAA_IH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][6];
    }
    double GetDIAA_II(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][7];
    }
    double GetDIAA_IK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][8];
    }
    double GetDIAA_IL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][9];
    }
    double GetDIAA_IM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][10];
    }
    double GetDIAA_IN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][11];
    }
    double GetDIAA_IP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][12];
    }
    double GetDIAA_IQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][13];
    }
    double GetDIAA_IR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][14];
    }
    double GetDIAA_IS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][15];
    }
    double GetDIAA_IT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][16];
    }
    double GetDIAA_IV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][17];
    }
    double GetDIAA_IW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][18];
    }
    double GetDIAA_IY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[7][19];
    }
    double GetDIAA_KA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][0];
    }
    double GetDIAA_KC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][1];
    }
    double GetDIAA_KD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][2];
    }
    double GetDIAA_KE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][3];
    }
    double GetDIAA_KF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][4];
    }
    double GetDIAA_KG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][5];
    }
    double GetDIAA_KH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][6];
    }
    double GetDIAA_KI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][7];
    }
    double GetDIAA_KK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][8];
    }
    double GetDIAA_KL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][9];
    }
    double GetDIAA_KM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][10];
    }
    double GetDIAA_KN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][11];
    }
    double GetDIAA_KP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][12];
    }
    double GetDIAA_KQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][13];
    }
    double GetDIAA_KR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][14];
    }
    double GetDIAA_KS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][15];
    }
    double GetDIAA_KT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][16];
    }
    double GetDIAA_KV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][17];
    }
    double GetDIAA_KW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][18];
    }
    double GetDIAA_KY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[8][19];
    }
    double GetDIAA_LA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][0];
    }
    double GetDIAA_LC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][1];
    }
    double GetDIAA_LD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][2];
    }
    double GetDIAA_LE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][3];
    }
    double GetDIAA_LF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][4];
    }
    double GetDIAA_LG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][5];
    }
    double GetDIAA_LH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][6];
    }
    double GetDIAA_LI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][7];
    }
    double GetDIAA_LK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][8];
    }
    double GetDIAA_LL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][9];
    }
    double GetDIAA_LM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][10];
    }
    double GetDIAA_LN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][11];
    }
    double GetDIAA_LP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][12];
    }
    double GetDIAA_LQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][13];
    }
    double GetDIAA_LR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][14];
    }
    double GetDIAA_LS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][15];
    }
    double GetDIAA_LT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][16];
    }
    double GetDIAA_LV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][17];
    }
    double GetDIAA_LW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][18];
    }
    double GetDIAA_LY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[9][19];
    }
    double GetDIAA_MA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][0];
    }
    double GetDIAA_MC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][1];
    }
    double GetDIAA_MD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][2];
    }
    double GetDIAA_ME(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][3];
    }
    double GetDIAA_MF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][4];
    }
    double GetDIAA_MG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][5];
    }
    double GetDIAA_MH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][6];
    }
    double GetDIAA_MI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][7];
    }
    double GetDIAA_MK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][8];
    }
    double GetDIAA_ML(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][9];
    }
    double GetDIAA_MM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][10];
    }
    double GetDIAA_MN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][11];
    }
    double GetDIAA_MP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][12];
    }
    double GetDIAA_MQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][13];
    }
    double GetDIAA_MR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][14];
    }
    double GetDIAA_MS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][15];
    }
    double GetDIAA_MT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][16];
    }
    double GetDIAA_MV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][17];
    }
    double GetDIAA_MW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][18];
    }
    double GetDIAA_MY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[10][19];
    }
    double GetDIAA_NA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][0];
    }
    double GetDIAA_NC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][1];
    }
    double GetDIAA_ND(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][2];
    }
    double GetDIAA_NE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][3];
    }
    double GetDIAA_NF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][4];
    }
    double GetDIAA_NG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][5];
    }
    double GetDIAA_NH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][6];
    }
    double GetDIAA_NI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][7];
    }
    double GetDIAA_NK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][8];
    }
    double GetDIAA_NL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][9];
    }
    double GetDIAA_NM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][10];
    }
    double GetDIAA_NN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][11];
    }
    double GetDIAA_NP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][12];
    }
    double GetDIAA_NQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][13];
    }
    double GetDIAA_NR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][14];
    }
    double GetDIAA_NS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][15];
    }
    double GetDIAA_NT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][16];
    }
    double GetDIAA_NV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][17];
    }
    double GetDIAA_NW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][18];
    }
    double GetDIAA_NY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[11][19];
    }
    double GetDIAA_PA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][0];
    }
    double GetDIAA_PC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][1];
    }
    double GetDIAA_PD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][2];
    }
    double GetDIAA_PE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][3];
    }
    double GetDIAA_PF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][4];
    }
    double GetDIAA_PG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][5];
    }
    double GetDIAA_PH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][6];
    }
    double GetDIAA_PI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][7];
    }
    double GetDIAA_PK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][8];
    }
    double GetDIAA_PL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][9];
    }
    double GetDIAA_PM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][10];
    }
    double GetDIAA_PN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][11];
    }
    double GetDIAA_PP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][12];
    }
    double GetDIAA_PQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][13];
    }
    double GetDIAA_PR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][14];
    }
    double GetDIAA_PS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][15];
    }
    double GetDIAA_PT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][16];
    }
    double GetDIAA_PV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][17];
    }
    double GetDIAA_PW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][18];
    }
    double GetDIAA_PY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[12][19];
    }
    double GetDIAA_QA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][0];
    }
    double GetDIAA_QC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][1];
    }
    double GetDIAA_QD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][2];
    }
    double GetDIAA_QE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][3];
    }
    double GetDIAA_QF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][4];
    }
    double GetDIAA_QG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][5];
    }
    double GetDIAA_QH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][6];
    }
    double GetDIAA_QI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][7];
    }
    double GetDIAA_QK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][8];
    }
    double GetDIAA_QL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][9];
    }
    double GetDIAA_QM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][10];
    }
    double GetDIAA_QN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][11];
    }
    double GetDIAA_QP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][12];
    }
    double GetDIAA_QQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][13];
    }
    double GetDIAA_QR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][14];
    }
    double GetDIAA_QS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][15];
    }
    double GetDIAA_QT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][16];
    }
    double GetDIAA_QV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][17];
    }
    double GetDIAA_QW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][18];
    }
    double GetDIAA_QY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[13][19];
    }
    double GetDIAA_RA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][0];
    }
    double GetDIAA_RC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][1];
    }
    double GetDIAA_RD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][2];
    }
    double GetDIAA_RE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][3];
    }
    double GetDIAA_RF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][4];
    }
    double GetDIAA_RG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][5];
    }
    double GetDIAA_RH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][6];
    }
    double GetDIAA_RI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][7];
    }
    double GetDIAA_RK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][8];
    }
    double GetDIAA_RL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][9];
    }
    double GetDIAA_RM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][10];
    }
    double GetDIAA_RN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][11];
    }
    double GetDIAA_RP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][12];
    }
    double GetDIAA_RQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][13];
    }
    double GetDIAA_RR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][14];
    }
    double GetDIAA_RS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][15];
    }
    double GetDIAA_RT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][16];
    }
    double GetDIAA_RV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][17];
    }
    double GetDIAA_RW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][18];
    }
    double GetDIAA_RY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[14][19];
    }
    double GetDIAA_SA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][0];
    }
    double GetDIAA_SC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][1];
    }
    double GetDIAA_SD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][2];
    }
    double GetDIAA_SE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][3];
    }
    double GetDIAA_SF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][4];
    }
    double GetDIAA_SG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][5];
    }
    double GetDIAA_SH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][6];
    }
    double GetDIAA_SI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][7];
    }
    double GetDIAA_SK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][8];
    }
    double GetDIAA_SL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][9];
    }
    double GetDIAA_SM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][10];
    }
    double GetDIAA_SN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][11];
    }
    double GetDIAA_SP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][12];
    }
    double GetDIAA_SQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][13];
    }
    double GetDIAA_SR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][14];
    }
    double GetDIAA_SS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][15];
    }
    double GetDIAA_ST(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][16];
    }
    double GetDIAA_SV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][17];
    }
    double GetDIAA_SW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][18];
    }
    double GetDIAA_SY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[15][19];
    }
    double GetDIAA_TA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][0];
    }
    double GetDIAA_TC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][1];
    }
    double GetDIAA_TD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][2];
    }
    double GetDIAA_TE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][3];
    }
    double GetDIAA_TF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][4];
    }
    double GetDIAA_TG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][5];
    }
    double GetDIAA_TH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][6];
    }
    double GetDIAA_TI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][7];
    }
    double GetDIAA_TK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][8];
    }
    double GetDIAA_TL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][9];
    }
    double GetDIAA_TM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][10];
    }
    double GetDIAA_TN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][11];
    }
    double GetDIAA_TP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][12];
    }
    double GetDIAA_TQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][13];
    }
    double GetDIAA_TR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][14];
    }
    double GetDIAA_TS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][15];
    }
    double GetDIAA_TT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][16];
    }
    double GetDIAA_TV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][17];
    }
    double GetDIAA_TW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][18];
    }
    double GetDIAA_TY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[16][19];
    }
    double GetDIAA_VA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][0];
    }
    double GetDIAA_VC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][1];
    }
    double GetDIAA_VD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][2];
    }
    double GetDIAA_VE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][3];
    }
    double GetDIAA_VF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][4];
    }
    double GetDIAA_VG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][5];
    }
    double GetDIAA_VH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][6];
    }
    double GetDIAA_VI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][7];
    }
    double GetDIAA_VK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][8];
    }
    double GetDIAA_VL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][9];
    }
    double GetDIAA_VM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][10];
    }
    double GetDIAA_VN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][11];
    }
    double GetDIAA_VP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][12];
    }
    double GetDIAA_VQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][13];
    }
    double GetDIAA_VR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][14];
    }
    double GetDIAA_VS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][15];
    }
    double GetDIAA_VT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][16];
    }
    double GetDIAA_VV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][17];
    }
    double GetDIAA_VW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][18];
    }
    double GetDIAA_VY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[17][19];
    }
    double GetDIAA_WA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][0];
    }
    double GetDIAA_WC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][1];
    }
    double GetDIAA_WD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][2];
    }
    double GetDIAA_WE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][3];
    }
    double GetDIAA_WF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][4];
    }
    double GetDIAA_WG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][5];
    }
    double GetDIAA_WH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][6];
    }
    double GetDIAA_WI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][7];
    }
    double GetDIAA_WK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][8];
    }
    double GetDIAA_WL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][9];
    }
    double GetDIAA_WM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][10];
    }
    double GetDIAA_WN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][11];
    }
    double GetDIAA_WP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][12];
    }
    double GetDIAA_WQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][13];
    }
    double GetDIAA_WR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][14];
    }
    double GetDIAA_WS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][15];
    }
    double GetDIAA_WT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][16];
    }
    double GetDIAA_WV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][17];
    }
    double GetDIAA_WW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][18];
    }
    double GetDIAA_WY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[18][19];
    }
    double GetDIAA_YA(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][0];
    }
    double GetDIAA_YC(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][1];
    }
    double GetDIAA_YD(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][2];
    }
    double GetDIAA_YE(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][3];
    }
    double GetDIAA_YF(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][4];
    }
    double GetDIAA_YG(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][5];
    }
    double GetDIAA_YH(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][6];
    }
    double GetDIAA_YI(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][7];
    }
    double GetDIAA_YK(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][8];
    }
    double GetDIAA_YL(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][9];
    }
    double GetDIAA_YM(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][10];
    }
    double GetDIAA_YN(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][11];
    }
    double GetDIAA_YP(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][12];
    }
    double GetDIAA_YQ(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][13];
    }
    double GetDIAA_YR(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][14];
    }
    double GetDIAA_YS(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][15];
    }
    double GetDIAA_YT(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][16];
    }
    double GetDIAA_YV(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][17];
    }
    double GetDIAA_YW(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][18];
    }
    double GetDIAA_YY(CodonSequenceAlignment* codondata)
    {
        if(!diaa_bool)
        {
            codondata->diaa_usage(diaa_usage);
            diaa_bool = true;
        }
        return diaa_usage[19][19];
    }

    double GetRSCUentropy(CodonSequenceAlignment* codondata)
    {
        if(!RSCUentropy_bool)
        {
            codondata->RSCUEntropy(RSCUentropy);
            RSCUentropy_bool = true;
        }

    }



};

#endif // SUMMARYSTATISTICS_H
