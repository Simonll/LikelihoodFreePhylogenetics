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
        double* aa_usage;
        double* aa_usagewo_nr;
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


        bool codon_bool;
        bool dinuc_bool;
        bool dinuc12_bool;
        bool dinuc23_bool;
        bool dinuc31_bool;
        bool aa_bool;
        bool aa_wo_nr;
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


        //Constructors
        SummaryStatistics(LocalParameters *lparam);
        SummaryStatistics(LocalData *ldata);


        //Destructors
        virtual ~SummaryStatistics();



        void MapFunctions();
        void computeSummariesFromData();
        void computeSummaries();
        void computeSummaries(int** curent_nodeleaf_sequence_codon);
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
        double Getnuc_site_comphet(CodonSequenceAlignment* codondata){
            if(!nuc_site_comphet_bool){
                nuc_site_comphet = codondata->nuc_site_comphet();
                nuc_site_comphet_bool = true;
            }
            return nuc_site_comphet;
        }
        double Getnuc1_site_comphet(CodonSequenceAlignment* codondata){
            if(!nuc1_site_comphet_bool){
                nuc1_site_comphet = codondata->nuc1_site_comphet();
                nuc1_site_comphet_bool = true;
            }
            return nuc1_site_comphet;
        }
        double Getnuc2_site_comphet(CodonSequenceAlignment* codondata){
            if(!nuc2_site_comphet_bool){
                nuc2_site_comphet = codondata->nuc2_site_comphet();
                nuc2_site_comphet_bool = true;
            }
            return nuc2_site_comphet;
        }
        double Getnuc3_site_comphet(CodonSequenceAlignment* codondata){
            if(!nuc3_site_comphet_bool){
                nuc3_site_comphet = codondata->nuc3_site_comphet();
                nuc3_site_comphet_bool = true;
            }
            return nuc3_site_comphet;
        }
        double Getnuc_taxa_comphet(CodonSequenceAlignment* codondata){
            if(!nuc_taxa_comphet_bool){
                nuc_taxa_comphet = codondata->nuc_taxa_comphet();
                nuc_taxa_comphet_bool = true;
            }
            return nuc_taxa_comphet;
        }
        double Getnuc1_taxa_comphet(CodonSequenceAlignment* codondata){
            if(!nuc1_taxa_comphet_bool){
                nuc1_taxa_comphet = codondata->nuc1_taxa_comphet();
                nuc1_taxa_comphet_bool = true;
            }
            return nuc1_taxa_comphet;
        }
        double Getnuc2_taxa_comphet(CodonSequenceAlignment* codondata){
            if(!nuc2_taxa_comphet_bool){
                nuc2_taxa_comphet = codondata->nuc2_taxa_comphet();
                nuc2_taxa_comphet_bool = true;
            }
            return nuc2_taxa_comphet;
        }
        double Getnuc3_taxa_comphet(CodonSequenceAlignment* codondata){
            if(!nuc3_taxa_comphet_bool){
                nuc3_taxa_comphet = codondata->nuc3_taxa_comphet();
                nuc3_taxa_comphet_bool = true;
            }
            return nuc3_taxa_comphet;
        }
        double Getcodon_site_comphet(CodonSequenceAlignment* codondata){
            if(!codon_site_comphet_bool){
                codon_site_comphet = codondata->codon_site_comphet();
                codon_site_comphet_bool = true;
            }
            return codon_site_comphet;
        }
        double Getcodon_taxa_comphet(CodonSequenceAlignment* codondata){
            if(!codon_taxa_comphet_bool){
                codon_taxa_comphet = codondata->codon_taxa_comphet();
                codon_taxa_comphet_bool = true;
            }
            return codon_taxa_comphet;
        }
        double Getaa_site_comphet(CodonSequenceAlignment* codondata){
            if(!aa_site_comphet_bool){
                aa_site_comphet = codondata->aa_site_comphet();
                aa_site_comphet_bool = true;
            }
            return aa_site_comphet;
        }
        double Getaa_taxa_comphet(CodonSequenceAlignment* codondata){
            if(!aa_taxa_comphet_bool){
                aa_taxa_comphet = codondata->aa_taxa_comphet();
                aa_taxa_comphet_bool = true;
            }
            return aa_taxa_comphet;
        }


        /////////////////
        // codon_usage
        /////////////////

         double GetGGG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGG")];
         }
          double GetGGA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGA")];
         }
         double GetGGC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGC")];
         }
         double GetGGT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GGT")];
         }
         double GetGAG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAG")];
         }
         double GetGAA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAA")];
         }
          double GetGAC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAC")];
         }
         double GetGAT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GAT")];
         }
         double GetGCG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCG")];
         }
         double GetGCA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCA")];
         }
         double GetGCC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCC")];
         }
         double GetGCT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GCT")];
         }
         double GetGTG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTG")];
         }
         double GetGTA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTA")];
         }
         double GetGTC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTC")];
         }
         double GetGTT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("GTT")];
         }
         double GetAGG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGG")];
         }
         double GetAGA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGA")];
         }
         double GetAGC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGC")];
         }
         double GetAGT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AGT")];
         }
         double GetAAG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAG")];
         }
         double GetAAA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAA")];
         }
         double GetAAC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAC")];
         }
         double GetAAT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("AAT")];
         }
         double GetACG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACG")];
         }
         double GetACA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACA")];
         }
         double GetACC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACC")];
         }
         double GetACT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ACT")];
         }
         double GetATG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATG")];
         }
         double GetATA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATA")];
         }
         double GetATC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATC")];
         }
         double GetATT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("ATT")];
         }
         double GetCGG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGG")];
         }
         double GetCGA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGA")];
         }
         double GetCGC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGC")];
         }
         double GetCGT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CGT")];
         }
         double GetCAG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAG")];
         }
         double GetCAA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAA")];
         }
         double GetCAC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAC")];
         }
        double GetCAT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CAT")];
         }
        double GetCCG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCG")];
         }
        double GetCCA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCA")];
         }
        double GetCCC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCC")];
         }
        double GetCCT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CCT")];
         }
        double GetCTG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTG")];
         }
        double GetCTA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTA")];
         }
        double GetCTC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTC")];
         }
        double GetCTT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("CTT")];
         }
         double GetTGG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGG")];
         }
         double GetTGA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGA")];
         }
         double GetTGC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGC")];
         }
         double GetTGT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TGT")];
         }
         double GetTAG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAG")];
         }

         double GetTAA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAA")];
         }

         double GetTAC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAC")];
         }

        double GetTAT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TAT")];
        }
        double GetTCG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCG")];
        }

        double GetTCA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCA")];
        }
        double GetTCC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCC")];
        }
        double GetTCT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TCT")];
        }

        double GetTTG(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTG")];
        }

        double GetTTA(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTA")];
        }
        double GetTTC(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTC")];
        }
        double GetTTT(CodonSequenceAlignment* codondata){
            if(!codon_bool) {
                codondata->CodonUsagePerAA(codon_usage);
                codon_bool = true;
            }
            return (double) codon_usage[codondata->GetCodonStateSpace()->GetState("TTT")];
        }





        /////////////////
        // aa_usage
        /////////////////

        double GetA(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[0];
        }
        double GetC(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[1];
        }
        double GetD(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[2];
        }
        double GetE(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[3];
        }
        double GetF(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[4];
        }
        double GetG(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[5];
        }
        double GetH(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[6];
        }
        double GetI(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[7];
        }
        double GetK(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[8];
        }
        double GetL(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[9];
        }
        double GetM(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[10];
        }
        double GetN(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[11];
        }
        double GetP(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[12];
        }
        double GetQ(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[13];
        }
        double GetR(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[14];
        }
        double GetS(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[15];
        }
        double GetT(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[16];
        }
        double GetV(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[17];
        }
        double GetW(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[18];
        }
        double GetY(CodonSequenceAlignment* codondata){
            if(!aa_bool) {
                codondata->aa_usage(aa_usage);
                aa_bool = true;
            }
            return (double) aa_usage[19];
        }



//        /////////////////
//        // aa_usagewo_nr
//        /////////////////
//        double GetAwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[0];
//        }
//        double GetCwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[1];
//        }
//        double GetDwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[2];
//        }
//        double GetEwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[3];
//        }
//        double GetFwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[4];
//        }
//        double GetGwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[5];
//        }
//        double GetHwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[6];
//        }
//        double GetIwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[7];
//        }
//        double GetKwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[8];
//        }
//        double GetLwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[9];
//        }
//        double GetMwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[10];
//        }
//        double GetNwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[11];
//        }
//        double GetPwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[12];
//        }
//        double GetQwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[13];
//        }
//        double GetRwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[14];
//        }
//        double GetSwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[15];
//        }
//        double GetTwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[16];
//        }
//        double GetVwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[17];
//        }
//        double GetWwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[18];
//        }
//        double GetYwonR(CodonSequenceAlignment* codondata){
//            if(!aa_bool) {
//                codondata->aa_usagewo_nr(aa_usagewo_nr);
//                aa_bool = true;
//            }
//            return (double) aa_usagewo_nr[19];
//        }




        /////////////////
        // dinuc31
        /////////////////
        double GetDinuc31AA(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[0][0];
        }

        double GetDinuc31AC(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[0][1];
        }

        double GetDinuc31AG(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[0][2];

        }
        double GetDinuc31AT(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[0][3];
        }
        double GetDinuc31CA(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[1][0];
        }

        double GetDinuc31CC(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[1][1];
        }
        double GetDinuc31CG(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[1][2];

        }
        double GetDinuc31CT(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[1][3];
        }
        double GetDinuc31GA(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[2][0];
        }

        double GetDinuc31GC(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[2][1];
        }
        double GetDinuc31GG(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[2][2];

        }
        double GetDinuc31GT(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[2][3];
        }
        double GetDinuc31TA(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[3][0];
        }

        double GetDinuc31TC(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[3][1];
        }
        double GetDinuc31TG(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[3][2];

        }
        double GetDinuc31TT(CodonSequenceAlignment* codondata){
            if(!dinuc31_bool) {
                codondata->dinuc31_usage(dinuc31_usage);
                dinuc31_bool = true;
            }
            return (double) dinuc31_usage[3][3];
        }
        /////////////////
        // dinuc23
        /////////////////
        double GetDinuc23AA(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[0][0];
        }

        double GetDinuc23AC(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[0][1];
        }
        double GetDinuc23AG(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[0][2];

        }
        double GetDinuc23AT(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[0][3];
        }
        double GetDinuc23CA(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[1][0];
        }

        double GetDinuc23CC(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[1][1];
        }
        double GetDinuc23CG(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[1][2];

        }
        double GetDinuc23CT(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[1][3];
        }
        double GetDinuc23GA(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[2][0];
        }

        double GetDinuc23GC(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[2][1];
        }
        double GetDinuc23GG(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[2][2];

        }
        double GetDinuc23GT(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[2][3];
        }
        double GetDinuc23TA(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[3][0];
        }

        double GetDinuc23TC(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[3][1];
        }
        double GetDinuc23TG(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[3][2];

        }
        double GetDinuc23TT(CodonSequenceAlignment* codondata){
            if(!dinuc23_bool) {
                codondata->dinuc23_usage(dinuc23_usage);
                dinuc23_bool = true;
            }
            return (double) dinuc23_usage[3][3];
        }

        /////////////////
        // dinuc12
        /////////////////
        double GetDinuc12AA(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[0][0];
        }

        double GetDinuc12AC(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[0][1];
        }
        double GetDinuc12AG(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[0][2];

        }
        double GetDinuc12AT(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[0][3];
        }
        double GetDinuc12CA(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[1][0];
        }

        double GetDinuc12CC(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[1][1];
        }
        double GetDinuc12CG(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[1][2];

        }
        double GetDinuc12CT(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[1][3];
        }
        double GetDinuc12GA(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[2][0];
        }

        double GetDinuc12GC(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[2][1];
        }
        double GetDinuc12GG(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[2][2];

        }
        double GetDinuc12GT(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[2][3];
        }
        double GetDinuc12TA(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[3][0];
        }

        double GetDinuc12TC(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[3][1];
        }
        double GetDinuc12TG(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[3][2];

        }
        double GetDinuc12TT(CodonSequenceAlignment* codondata){
            if(!dinuc12_bool) {
                codondata->dinuc12_usage(dinuc12_usage);
                dinuc12_bool = true;
            }
            return (double) dinuc12_usage[3][3];
        }
        /////////////////
        // dinuc
        /////////////////
        double GetDinucAA(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[0][0];
        }

        double GetDinucAC(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[0][1];
        }
        double GetDinucAG(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[0][2];

        }
        double GetDinucAT(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[0][3];
        }
        double GetDinucCA(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[1][0];
        }

        double GetDinucCC(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[1][1];
        }
        double GetDinucCG(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[1][2];

        }
        double GetDinucCT(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[1][3];
        }
        double GetDinucGA(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[2][0];
        }

        double GetDinucGC(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[2][1];
        }
        double GetDinucGG(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[2][2];

        }
        double GetDinucGT(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[2][3];
        }
        double GetDinucTA(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[3][0];
        }

        double GetDinucTC(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[3][1];
        }
        double GetDinucTG(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[3][2];

        }
        double GetDinucTT(CodonSequenceAlignment* codondata){
            if(!dinuc_bool) {
                codondata->dinuc_usage(dinuc_usage);
                dinuc_bool = true;
            }
            return (double) dinuc_usage[3][3];
        }

        /////////////////
        // nuc1
        /////////////////
        double GetNuc1A(CodonSequenceAlignment* codondata){
            if(!nuc1_bool) {
                codondata->nuc1_usage(nuc1_usage);
                nuc1_bool = true;
            }
            return (double) nuc1_usage[0];
        }

         double GetNuc1C(CodonSequenceAlignment* codondata){
            if(!nuc1_bool) {
                codondata->nuc1_usage(nuc1_usage);
                nuc1_bool = true;
            }
            return (double) nuc1_usage[1];
        }

         double GetNuc1G(CodonSequenceAlignment* codondata){
            if(!nuc1_bool) {
                codondata->nuc1_usage(nuc1_usage);
                nuc1_bool = true;
            }
            return (double) nuc1_usage[2];
        }

         double GetNuc1T(CodonSequenceAlignment* codondata){
            if(!nuc1_bool) {
                codondata->nuc1_usage(nuc1_usage);
                nuc1_bool = true;
            }
            return (double) nuc1_usage[3];
        }
        /////////////////
        // nuc2
        /////////////////
        double GetNuc2A(CodonSequenceAlignment* codondata){
            if(!nuc2_bool) {
                codondata->nuc2_usage(nuc2_usage);
                nuc2_bool = true;
            }
            return (double) nuc2_usage[0];
        }

         double GetNuc2C(CodonSequenceAlignment* codondata){
            if(!nuc2_bool) {
                codondata->nuc2_usage(nuc2_usage);
                nuc2_bool = true;
            }
            return (double) nuc2_usage[1];
        }

         double GetNuc2G(CodonSequenceAlignment* codondata){
            if(!nuc2_bool) {
                codondata->nuc2_usage(nuc2_usage);
                nuc2_bool = true;
            }
            return (double) nuc2_usage[2];
        }

         double GetNuc2T(CodonSequenceAlignment* codondata){
            if(!nuc2_bool) {
                codondata->nuc2_usage(nuc2_usage);
                nuc2_bool = true;
            }
            return (double) nuc2_usage[3];
        }
        /////////////////
        // nuc3
        /////////////////
        double GetNuc3A(CodonSequenceAlignment* codondata){
            if(!nuc3_bool) {
                codondata->nuc3_usage(nuc3_usage);
                nuc3_bool = true;
            }
            return (double) nuc3_usage[0];
        }

         double GetNuc3C(CodonSequenceAlignment* codondata){
            if(!nuc3_bool) {
                codondata->nuc3_usage(nuc3_usage);
                nuc3_bool = true;
            }
            return (double) nuc3_usage[1];
        }

         double GetNuc3G(CodonSequenceAlignment* codondata){
            if(!nuc3_bool) {
                codondata->nuc3_usage(nuc3_usage);
                nuc3_bool = true;
            }
            return (double) nuc3_usage[2];
        }

         double GetNuc3T(CodonSequenceAlignment* codondata){
            if(!nuc3_bool) {
                codondata->nuc3_usage(nuc3_usage);
                nuc3_bool = true;
            }
            return (double) nuc3_usage[3];
        }
        /////////////////
        // nuc
        /////////////////
        double GetNucA(CodonSequenceAlignment* codondata){
            if(!nuc_bool) {
                codondata->nuc_usage(nuc_usage);
                nuc_bool = true;
            }
            return (double) nuc_usage[0];
        }

         double GetNucC(CodonSequenceAlignment* codondata){
            if(!nuc_bool) {
                codondata->nuc_usage(nuc_usage);
                nuc_bool = true;
            }
            return (double) nuc_usage[1];
        }

         double GetNucG(CodonSequenceAlignment* codondata){
            if(!nuc_bool) {
                codondata->nuc_usage(nuc_usage);
                nuc_bool = true;
            }
            return (double) nuc_usage[2];
        }

         double GetNucT(CodonSequenceAlignment* codondata){
            if(!nuc_bool) {
                codondata->nuc_usage(nuc_usage);
                nuc_bool = true;
            }
            return (double) nuc_usage[3];
        }


        /////////////////
        // nuc3_meandiff
        /////////////////
        double GetNuc3mean(CodonSequenceAlignment* codondata){
            if(!nuc3_meandiff_bool) {
                codondata->nuc3_meandiff(nuc3_meandiff);
                nuc3_meandiff_bool = true;
            }
            return (double) nuc3_meandiff[0];
        }

        double GetNuc3var(CodonSequenceAlignment* codondata){
            if(!nuc3_meandiff_bool) {
                codondata->nuc3_meandiff(nuc3_meandiff);
                nuc3_meandiff_bool = true;
            }
            return (double) nuc3_meandiff[3];
        }
        /////////////////
        // nuc2_meandiff
        /////////////////
        double GetNuc2mean(CodonSequenceAlignment* codondata){
            if(!nuc2_meandiff_bool) {
                codondata->nuc2_meandiff(nuc2_meandiff);
                nuc2_meandiff_bool = true;
            }
            return (double) nuc2_meandiff[0];
        }

        double GetNuc2var(CodonSequenceAlignment* codondata){
            if(!nuc2_meandiff_bool) {
                codondata->nuc2_meandiff(nuc2_meandiff);
                nuc2_meandiff_bool = true;
            }
            return (double) nuc2_meandiff[2];
        }
        /////////////////
        // nuc1_meandiff
        /////////////////
        double GetNuc1mean(CodonSequenceAlignment* codondata){
            if(!nuc1_meandiff_bool) {
                codondata->nuc1_meandiff(nuc1_meandiff);
                nuc1_meandiff_bool = true;
            }
            return (double) nuc1_meandiff[0];
        }

        double GetNuc1var(CodonSequenceAlignment* codondata){
            if(!nuc1_meandiff_bool) {
                codondata->nuc1_meandiff(nuc1_meandiff);
                nuc1_meandiff_bool = true;
            }
            return (double) nuc1_meandiff[1];
        }
        /////////////////
       // nuc_meandiff
       /////////////////
       double GetNucmean(CodonSequenceAlignment* codondata){
           if(!nuc_meandiff_bool) {
               codondata->nuc_meandiff(nuc_meandiff);
               nuc_meandiff_bool = true;
           }
           return (double) nuc_meandiff[0];
       }

       double GetNucvar(CodonSequenceAlignment* codondata){
           if(!nuc_meandiff_bool) {
               codondata->nuc_meandiff(nuc_meandiff);
               nuc_meandiff_bool = true;
           }
           return (double) nuc_meandiff[1];
       }

        /////////////////
       // codon_meandiff
       /////////////////
       double GetCodonmean(CodonSequenceAlignment* codondata){
           if(!codon_meandiff_bool) {
               codondata->codon_meandiff(codon_meandiff);
               codon_meandiff_bool = true;
           }
           return (double) codon_meandiff[0];
       }

       double GetCodonvar(CodonSequenceAlignment* codondata){
           if(!codon_meandiff_bool) {
               codondata->codon_meandiff(codon_meandiff);
               codon_meandiff_bool = true;
           }
           return (double) codon_meandiff[1];
       }
       /////////////////
       // aa_meandiff
       /////////////////
       double GetAAmean(CodonSequenceAlignment* codondata){
           if(!aa_meandiff_bool) {
               codondata->aa_meandiff(aa_meandiff);
               aa_meandiff_bool = true;
           }
           return (double) aa_meandiff[0];
       }

       double GetAAvar(CodonSequenceAlignment* codondata){
           if(!aa_meandiff_bool) {
               codondata->aa_meandiff(aa_meandiff);
               aa_meandiff_bool = true;
           }
           return (double) aa_meandiff[1];
       }


       double GetAAmean_wonR(CodonSequenceAlignment* codondata){
           if(!aa_meandiff_bool) {
               codondata->aa_meandiff_wonr(aa_wonR_meandiff);
               aa_wonR_meandiff_bool = true;
           }
           return (double) aa_wonR_meandiff[0];
       }

       double GetAAvar_wonR(CodonSequenceAlignment* codondata){
           if(!aa_meandiff_bool) {
               codondata->aa_meandiff_wonr(aa_wonR_meandiff);
               aa_wonR_meandiff_bool = true;
           }
           return (double) aa_wonR_meandiff[1];
       }

       /////////////////
       // dinucCpG_pairwise
       /////////////////
       double GetdinucCpG_TpG(CodonSequenceAlignment* codondata){
           if(!dinucCpG_pairwise_bool) {
               codondata->dinucCpG_pairwise(dinucCpG_pairwise);
               dinucCpG_pairwise_bool = true;
           }
           return (double) dinucCpG_pairwise[0];
       }

       double GetdinucCpG_CpA(CodonSequenceAlignment* codondata){
           if(!dinucCpG_pairwise_bool) {
               codondata->dinucCpG_pairwise(dinucCpG_pairwise);
               dinucCpG_pairwise_bool = true;
           }
           return (double) dinucCpG_pairwise[1];
       }

       double GetdinucApG_TpG(CodonSequenceAlignment* codondata){
           if(!dinucCpG_pairwise_bool) {
               codondata->dinucCpG_pairwise(dinucCpG_pairwise);
               dinucCpG_pairwise_bool = true;
           }
           return (double) dinucCpG_pairwise[2];
       }


       /////////////////
       // aa_pairwise
       /////////////////
       double GetpwAA(CodonSequenceAlignment* codondata){
            if(!aa_pairwise_bool) {
                codondata->aa_pairwise(aa_pairwise);
                aa_pairwise_bool = true;
            }
            double sum = 0;
            for (int i = 0 ; i < 190; i++) {
                sum += aa_pairwise[i];
            }
            return sum;
        }




       /////////////////
       // nuc3_pairwise
       /////////////////

       double Getpw3GT(CodonSequenceAlignment* codondata){
           if(!nuc3_pairwise_bool) {
               codondata->nuc_pairwise(nuc3_pairwise);
               nuc3_pairwise_bool = true;
           }
           return (double) nuc3_pairwise[5];
       }

       double Getpw3CT(CodonSequenceAlignment* codondata){
           if(!nuc3_pairwise_bool) {
               codondata->nuc_pairwise(nuc3_pairwise);
               nuc3_pairwise_bool = true;
           }
           return (double) nuc3_pairwise[4];
       }

       double Getpw3CG(CodonSequenceAlignment* codondata){
           if(!nuc3_pairwise_bool) {
               codondata->nuc_pairwise(nuc3_pairwise);
               nuc3_pairwise_bool = true;
           }
           return (double) nuc3_pairwise[3];
       }


       double Getpw3AT(CodonSequenceAlignment* codondata){
           if(!nuc3_pairwise_bool) {
               codondata->nuc_pairwise(nuc3_pairwise);
               nuc3_pairwise_bool = true;
           }
           return (double) nuc3_pairwise[2];
       }

       double Getpw3AG(CodonSequenceAlignment* codondata){
           if(!nuc3_pairwise_bool) {
               codondata->nuc_pairwise(nuc3_pairwise);
               nuc3_pairwise_bool = true;
           }
           return (double) nuc3_pairwise[1];
       }


       double Getpw3AC(CodonSequenceAlignment* codondata){
           if(!nuc3_pairwise_bool) {
               codondata->nuc_pairwise(nuc3_pairwise);
               nuc3_pairwise_bool = true;
           }
           return (double) nuc3_pairwise[0];
       }

       /////////////////
       // nuc2_pairwise
       /////////////////

     double Getpw2GT(CodonSequenceAlignment* codondata){
          if(!nuc2_pairwise_bool) {
              codondata->nuc_pairwise(nuc2_pairwise);
              nuc2_pairwise_bool = true;
          }
          return (double) nuc2_pairwise[5];
      }

      double Getpw2CT(CodonSequenceAlignment* codondata){
          if(!nuc2_pairwise_bool) {
              codondata->nuc_pairwise(nuc2_pairwise);
              nuc2_pairwise_bool = true;
          }
          return (double) nuc2_pairwise[4];
      }

      double Getpw2CG(CodonSequenceAlignment* codondata){
          if(!nuc2_pairwise_bool) {
              codondata->nuc_pairwise(nuc2_pairwise);
              nuc2_pairwise_bool = true;
          }
          return (double) nuc2_pairwise[3];
      }


      double Getpw2AT(CodonSequenceAlignment* codondata){
          if(!nuc2_pairwise_bool) {
              codondata->nuc_pairwise(nuc2_pairwise);
              nuc2_pairwise_bool = true;
          }
          return (double) nuc2_pairwise[2];
      }

      double Getpw2AG(CodonSequenceAlignment* codondata){
          if(!nuc2_pairwise_bool) {
              codondata->nuc_pairwise(nuc2_pairwise);
              nuc2_pairwise_bool = true;
          }
          return (double) nuc2_pairwise[1];
      }


      double Getpw2AC(CodonSequenceAlignment* codondata){
          if(!nuc2_pairwise_bool) {
              codondata->nuc_pairwise(nuc2_pairwise);
              nuc2_pairwise_bool = true;
          }
          return (double) nuc2_pairwise[0];
      }

        /////////////////
        // nuc1_pairwise
        /////////////////

     double Getpw1GT(CodonSequenceAlignment* codondata){
          if(!nuc1_pairwise_bool) {
              codondata->nuc_pairwise(nuc1_pairwise);
              nuc1_pairwise_bool = true;
          }
          return (double) nuc1_pairwise[5];
      }

      double Getpw1CT(CodonSequenceAlignment* codondata){
          if(!nuc1_pairwise_bool) {
              codondata->nuc_pairwise(nuc1_pairwise);
              nuc1_pairwise_bool = true;
          }
          return (double) nuc1_pairwise[4];
      }

      double Getpw1CG(CodonSequenceAlignment* codondata){
          if(!nuc1_pairwise_bool) {
              codondata->nuc_pairwise(nuc1_pairwise);
              nuc1_pairwise_bool = true;
          }
          return (double) nuc1_pairwise[3];
      }


      double Getpw1AT(CodonSequenceAlignment* codondata){
          if(!nuc1_pairwise_bool) {
              codondata->nuc_pairwise(nuc1_pairwise);
              nuc1_pairwise_bool = true;
          }
          return (double) nuc1_pairwise[2];
      }

      double Getpw1AG(CodonSequenceAlignment* codondata){
          if(!nuc1_pairwise_bool) {
              codondata->nuc_pairwise(nuc1_pairwise);
              nuc1_pairwise_bool = true;
          }
          return (double) nuc1_pairwise[1];
      }


      double Getpw1AC(CodonSequenceAlignment* codondata){
          if(!nuc1_pairwise_bool) {
              codondata->nuc_pairwise(nuc1_pairwise);
              nuc1_pairwise_bool = true;
          }
          return (double) nuc1_pairwise[0];
      }


        /////////////////
        // nuc_pairwise
        /////////////////


        double GetpwGT(CodonSequenceAlignment* codondata){
            if(!nuc_pairwise_bool) {
                codondata->nuc_pairwise(nuc_pairwise);
                nuc_pairwise_bool = true;
            }
            return (double) nuc_pairwise[5];
        }

        double GetpwCT(CodonSequenceAlignment* codondata){
            if(!nuc_pairwise_bool) {
                codondata->nuc_pairwise(nuc_pairwise);
                nuc_pairwise_bool = true;
            }
            return (double) nuc_pairwise[4];
        }

        double GetpwCG(CodonSequenceAlignment* codondata){
            if(!nuc_pairwise_bool) {
                codondata->nuc_pairwise(nuc_pairwise);
                nuc_pairwise_bool = true;
            }
            return (double) nuc_pairwise[3];
        }


        double GetpwAT(CodonSequenceAlignment* codondata){
            if(!nuc_pairwise_bool) {
                codondata->nuc_pairwise(nuc_pairwise);
                nuc_pairwise_bool = true;
            }
            return (double) nuc_pairwise[2];
        }

        double GetpwAG(CodonSequenceAlignment* codondata){
            if(!nuc_pairwise_bool) {
                codondata->nuc_pairwise(nuc_pairwise);
                nuc_pairwise_bool = true;
            }
            return (double) nuc_pairwise[1];
        }


        double GetpwAC(CodonSequenceAlignment* codondata){
            if(!nuc_pairwise_bool) {
                codondata->nuc_pairwise(nuc_pairwise);
                nuc_pairwise_bool = true;
            }
            return (double) nuc_pairwise[0];
        }


        double GetCGNAGRcgnconst(CodonSequenceAlignment* codondata) {
                if(!CGNAGR_bool) {
                    codondata->CGNAGR(CGNAGR);
                    CGNAGR_bool = true;
                }
            return (double) CGNAGR[0];

        }

        double GetCGNAGRagrconst(CodonSequenceAlignment* codondata) {
                if(!CGNAGR_bool) {
                    codondata->CGNAGR(CGNAGR);
                    CGNAGR_bool = true;
                }
            return (double) CGNAGR[1];

        }


        double GetCGNAGRcgnvar(CodonSequenceAlignment* codondata) {
                if(!CGNAGR_bool) {
                    codondata->CGNAGR(CGNAGR);
                    CGNAGR_bool = true;
                }
            return (double) CGNAGR[2];

        }

        double GetCGNAGRagrvar(CodonSequenceAlignment* codondata) {
                if(!CGNAGR_bool) {
                    codondata->CGNAGR(CGNAGR);
                    CGNAGR_bool = true;
                }
            return (double) CGNAGR[3];

        }

};

#endif // SUMMARYSTATISTICS_H
