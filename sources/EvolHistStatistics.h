#ifndef EVOLHISTSTATISTICS_H
#define EVOLHISTSTATISTICS_H



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

#include "LocalParameters.h"


class EvolHistStatistics
{
public:


    LocalParameters* lparam;


    //Mappings Statistics
    double Nsub;
    double Nsynsub;
    int* ssNsub;
    int* ssNnsynsub;
    int* ssNsynsub;

    int* ssNsubbin;
    int* ssNsynsubbin;


    double** branch_stat;
    double** MutRate;
    double** SubRate;
    double*** gtnr_stat;
    double*** dinuc_stat;
    double*** codon_stat;
    //Synonymous
    double*** gtnrSyn_stat;
    double*** dinucSyn_stat;
    double*** codonSyn_stat;
    //Nonsynonymous
    double*** gtnrNSyn_stat;
    double*** dinucNSyn_stat;
    double*** codonNSyn_stat;

    int**** ssgtnr_stat;
    int**** ssdinuc_stat;
    int**** sscodon_stat;
    //Synonymous
    int**** ssgtnrSyn_stat;
    int**** ssdinucSyn_stat;
    int**** sscodonSyn_stat;
    //Nonsynonymous
    int**** ssgtnrNSyn_stat;
    int**** ssdinucNSyn_stat;
    int**** sscodonNSyn_stat;



    typedef double (EvolHistStatistics::*funcpt)();
    typedef std::vector<double> (EvolHistStatistics::*funcptwtvector)();

    std::map<string,funcpt> GetEvoStatMap;
    std::map<string,funcptwtvector> GetSiteSpecificEvoStatsMap;


    EvolHistStatistics(LocalParameters* inparam);

    virtual ~EvolHistStatistics();

    //Setters
    void resetEvoStats();


    //Getters
    void GetEvoStats();
    void GetSiteSpecificEvoStats();
    void GetEvoAncStats();
    int  GetDinucContext(int nuc_a, int nuc_b)
    {
        int k = 0;
        for (int i = 0; i < 4; i++)
        {
            for (int j =0 ; j <4; j++)
            {
                if(nuc_a == i && nuc_b == j)
                {
                    return k;
                }
                k++;
            }
        }

    }


//        void tobstats(ostream & os,Tree* refTree) {
//            tobstats(os, refTree->GetRoot(),refTree);
//        }
//
//        void tobstats(ostream & os, const Link* from,Tree* refTree) {
//          if (!from)  {
//              from = refTree->GetRoot();
//          }
//          if (from->isLeaf()) {
//              os << from->GetNode()->GetName();
//          }
//          else    {
//              os << '(';
//              for (const Link* link = from->Next(); link!=from; link=link->Next())    {
//                  tobstats(os,link->Out(),refTree);
//                  if (link->Next() != from)   {
//                      os << ',';
//                  }
//              }
//              os << ')';
//          }
//          if (from->isRoot()) {
//              os << ";\n";
//          }
//          else {
//              os <<":"  << branch_stat[from->GetNode()->GetIndex()][0] << "_"  << branch_stat[from->GetNode()->GetIndex()][1] << "_"  <<  branch_stat[from->GetNode()->GetIndex()][2] << "_"  << branch_stat[from->GetNode()->GetIndex()][3] << "_"  << branch_stat[from->GetNode()->GetIndex()][4] << "_"  << branch_stat[from->GetNode()->GetIndex()][5];
//          }
//
//        }
//
//        void oSsstats(ostream & os){
//            os << "Nsub\t" << "Nnsynsub\t" << "Nsynsub\n";
//            for (int site_codon = 0 ; site_codon < param->Nsite_codon; site_codon++) {
//                os << ssNsub[site_codon] <<  "\t" << (ssNsub[site_codon] - ssNsynsub[site_codon]) << "\t" << ssNsynsub[site_codon] << "\n";
//            }
//        }
protected:

private:


    std::vector<double> GetssNsub()
    {
        std::vector<double> vec;


        int* bin  = new int[100];
        for (int i = 0 ; i < 100; i++)
        {
            bin[i] = 0;

        }

        for (int i = 0; i < lparam->Nsite_codon; i++)
        {
            int c  = ssNsub[i];

            if (ssNsub[i] >= 100 )
            {
                c = 100-1;

            }
            bin[c]++;

        }

        for (int i = 0; i < 100; i++)
        {

            vec.push_back(static_cast<double> (bin[i]));

        }


        delete [] bin;

        return vec;

    }

    std::vector<double> GetssNsynsub()
    {
        std::vector<double> vec;


        int* bin  = new int[100];
        for (int i = 0 ; i < 100; i++)
        {
            bin[i] = 0;

        }

        for (int i = 0; i<lparam->Nsite_codon; i++)
        {
            int c  = ssNsynsub[i];

            if (ssNsynsub[i] >= 100 )
            {
                c = 100-1;

            }
            bin[c]++;

        }

        for (int i = 0; i<100; i++)
        {

            vec.push_back(static_cast<double> (bin[i]));

        }


        delete [] bin;

        return vec;

    }




    double GetGTNR_AA()
    {
        return gtnr_stat[3][0][0];
    }
    double GetGTNR_AC()
    {
        return gtnr_stat[3][0][1];
    }
    double GetGTNR_AG()
    {
        return gtnr_stat[3][0][2];
    }
    double GetGTNR_AT()
    {
        return gtnr_stat[3][0][3];
    }
    double GetGTNR_CA()
    {
        return gtnr_stat[3][1][0];
    }
    double GetGTNR_CC()
    {
        return gtnr_stat[3][1][1];
    }
    double GetGTNR_CG()
    {
        return gtnr_stat[3][1][2];
    }
    double GetGTNR_CT()
    {
        return gtnr_stat[3][1][3];
    }
    double GetGTNR_GA()
    {
        return gtnr_stat[3][2][0];
    }
    double GetGTNR_GC()
    {
        return gtnr_stat[3][2][1];
    }
    double GetGTNR_GG()
    {
        return gtnr_stat[3][2][2];
    }
    double GetGTNR_GT()
    {
        return gtnr_stat[3][2][3];
    }
    double GetGTNR_TA()
    {
        return gtnr_stat[3][3][0];
    }
    double GetGTNR_TC()
    {
        return gtnr_stat[3][3][1];
    }
    double GetGTNR_TG()
    {
        return gtnr_stat[3][3][2];
    }
    double GetGTNR_TT()
    {
        return gtnr_stat[3][3][3];
    }
    double GetGTNR1_AA()
    {
        return gtnr_stat[0][0][0];
    }
    double GetGTNR1_AC()
    {
        return gtnr_stat[0][0][1];
    }
    double GetGTNR1_AG()
    {
        return gtnr_stat[0][0][2];
    }
    double GetGTNR1_AT()
    {
        return gtnr_stat[0][0][3];
    }
    double GetGTNR1_CA()
    {
        return gtnr_stat[0][1][0];
    }
    double GetGTNR1_CC()
    {
        return gtnr_stat[0][1][1];
    }
    double GetGTNR1_CG()
    {
        return gtnr_stat[0][1][2];
    }
    double GetGTNR1_CT()
    {
        return gtnr_stat[0][1][3];
    }
    double GetGTNR1_GA()
    {
        return gtnr_stat[0][2][0];
    }
    double GetGTNR1_GC()
    {
        return gtnr_stat[0][2][1];
    }
    double GetGTNR1_GG()
    {
        return gtnr_stat[0][2][2];
    }
    double GetGTNR1_GT()
    {
        return gtnr_stat[0][2][3];
    }
    double GetGTNR1_TA()
    {
        return gtnr_stat[0][3][0];
    }
    double GetGTNR1_TC()
    {
        return gtnr_stat[0][3][1];
    }
    double GetGTNR1_TG()
    {
        return gtnr_stat[0][3][2];
    }
    double GetGTNR1_TT()
    {
        return gtnr_stat[0][3][3];
    }



    double GetGTNR2_AA()
    {
        return gtnr_stat[1][0][0];
    }
    double GetGTNR2_AC()
    {
        return gtnr_stat[1][0][1];
    }
    double GetGTNR2_AG()
    {
        return gtnr_stat[1][0][2];
    }
    double GetGTNR2_AT()
    {
        return gtnr_stat[1][0][3];
    }
    double GetGTNR2_CA()
    {
        return gtnr_stat[1][1][0];
    }
    double GetGTNR2_CC()
    {
        return gtnr_stat[1][1][1];
    }
    double GetGTNR2_CG()
    {
        return gtnr_stat[1][1][2];
    }
    double GetGTNR2_CT()
    {
        return gtnr_stat[1][1][3];
    }
    double GetGTNR2_GA()
    {
        return gtnr_stat[1][2][0];
    }
    double GetGTNR2_GC()
    {
        return gtnr_stat[1][2][1];
    }
    double GetGTNR2_GG()
    {
        return gtnr_stat[1][2][2];
    }
    double GetGTNR2_GT()
    {
        return gtnr_stat[1][2][3];
    }
    double GetGTNR2_TA()
    {
        return gtnr_stat[1][3][0];
    }
    double GetGTNR2_TC()
    {
        return gtnr_stat[1][3][1];
    }
    double GetGTNR2_TG()
    {
        return gtnr_stat[1][3][2];
    }
    double GetGTNR2_TT()
    {
        return gtnr_stat[1][3][3];
    }
    double GetGTNR3_AA()
    {
        return gtnr_stat[2][0][0];
    }
    double GetGTNR3_AC()
    {
        return gtnr_stat[2][0][1];
    }
    double GetGTNR3_AG()
    {
        return gtnr_stat[2][0][2];
    }
    double GetGTNR3_AT()
    {
        return gtnr_stat[2][0][3];
    }
    double GetGTNR3_CA()
    {
        return gtnr_stat[2][1][0];
    }
    double GetGTNR3_CC()
    {
        return gtnr_stat[2][1][1];
    }
    double GetGTNR3_CG()
    {
        return gtnr_stat[2][1][2];
    }
    double GetGTNR3_CT()
    {
        return gtnr_stat[2][1][3];
    }
    double GetGTNR3_GA()
    {
        return gtnr_stat[2][2][0];
    }
    double GetGTNR3_GC()
    {
        return gtnr_stat[2][2][1];
    }
    double GetGTNR3_GG()
    {
        return gtnr_stat[2][2][2];
    }
    double GetGTNR3_GT()
    {
        return gtnr_stat[2][2][3];
    }
    double GetGTNR3_TA()
    {
        return gtnr_stat[2][3][0];
    }
    double GetGTNR3_TC()
    {
        return gtnr_stat[2][3][1];
    }
    double GetGTNR3_TG()
    {
        return gtnr_stat[2][3][2];
    }
    double GetGTNR3_TT()
    {
        return gtnr_stat[2][3][3];
    }

    double GetGTNRSyn_AA()
    {
        return gtnrSyn_stat[3][0][0];
    }
    double GetGTNRSyn_AC()
    {
        return gtnrSyn_stat[3][0][1];
    }
    double GetGTNRSyn_AG()
    {
        return gtnrSyn_stat[3][0][2];
    }
    double GetGTNRSyn_AT()
    {
        return gtnrSyn_stat[3][0][3];
    }
    double GetGTNRSyn_CA()
    {
        return gtnrSyn_stat[3][1][0];
    }
    double GetGTNRSyn_CC()
    {
        return gtnrSyn_stat[3][1][1];
    }
    double GetGTNRSyn_CG()
    {
        return gtnrSyn_stat[3][1][2];
    }
    double GetGTNRSyn_CT()
    {
        return gtnrSyn_stat[3][1][3];
    }
    double GetGTNRSyn_GA()
    {
        return gtnrSyn_stat[3][2][0];
    }
    double GetGTNRSyn_GC()
    {
        return gtnrSyn_stat[3][2][1];
    }
    double GetGTNRSyn_GG()
    {
        return gtnrSyn_stat[3][2][2];
    }
    double GetGTNRSyn_GT()
    {
        return gtnrSyn_stat[3][2][3];
    }
    double GetGTNRSyn_TA()
    {
        return gtnrSyn_stat[3][3][0];
    }
    double GetGTNRSyn_TC()
    {
        return gtnrSyn_stat[3][3][1];
    }
    double GetGTNRSyn_TG()
    {
        return gtnrSyn_stat[3][3][2];
    }
    double GetGTNRSyn_TT()
    {
        return gtnrSyn_stat[3][3][3];
    }
    double GetGTNRSyn1_AA()
    {
        return gtnrSyn_stat[0][0][0];
    }
    double GetGTNRSyn1_AC()
    {
        return gtnrSyn_stat[0][0][1];
    }
    double GetGTNRSyn1_AG()
    {
        return gtnrSyn_stat[0][0][2];
    }
    double GetGTNRSyn1_AT()
    {
        return gtnrSyn_stat[0][0][3];
    }
    double GetGTNRSyn1_CA()
    {
        return gtnrSyn_stat[0][1][0];
    }
    double GetGTNRSyn1_CC()
    {
        return gtnrSyn_stat[0][1][1];
    }
    double GetGTNRSyn1_CG()
    {
        return gtnrSyn_stat[0][1][2];
    }
    double GetGTNRSyn1_CT()
    {
        return gtnrSyn_stat[0][1][3];
    }
    double GetGTNRSyn1_GA()
    {
        return gtnrSyn_stat[0][2][0];
    }
    double GetGTNRSyn1_GC()
    {
        return gtnrSyn_stat[0][2][1];
    }
    double GetGTNRSyn1_GG()
    {
        return gtnrSyn_stat[0][2][2];
    }
    double GetGTNRSyn1_GT()
    {
        return gtnrSyn_stat[0][2][3];
    }
    double GetGTNRSyn1_TA()
    {
        return gtnrSyn_stat[0][3][0];
    }
    double GetGTNRSyn1_TC()
    {
        return gtnrSyn_stat[0][3][1];
    }
    double GetGTNRSyn1_TG()
    {
        return gtnrSyn_stat[0][3][2];
    }
    double GetGTNRSyn1_TT()
    {
        return gtnrSyn_stat[0][3][3];
    }



    double GetGTNRSyn2_AA()
    {
        return gtnrSyn_stat[1][0][0];
    }
    double GetGTNRSyn2_AC()
    {
        return gtnrSyn_stat[1][0][1];
    }
    double GetGTNRSyn2_AG()
    {
        return gtnrSyn_stat[1][0][2];
    }
    double GetGTNRSyn2_AT()
    {
        return gtnrSyn_stat[1][0][3];
    }
    double GetGTNRSyn2_CA()
    {
        return gtnrSyn_stat[1][1][0];
    }
    double GetGTNRSyn2_CC()
    {
        return gtnrSyn_stat[1][1][1];
    }
    double GetGTNRSyn2_CG()
    {
        return gtnrSyn_stat[1][1][2];
    }
    double GetGTNRSyn2_CT()
    {
        return gtnrSyn_stat[1][1][3];
    }
    double GetGTNRSyn2_GA()
    {
        return gtnrSyn_stat[1][2][0];
    }
    double GetGTNRSyn2_GC()
    {
        return gtnrSyn_stat[1][2][1];
    }
    double GetGTNRSyn2_GG()
    {
        return gtnrSyn_stat[1][2][2];
    }
    double GetGTNRSyn2_GT()
    {
        return gtnrSyn_stat[1][2][3];
    }
    double GetGTNRSyn2_TA()
    {
        return gtnrSyn_stat[1][3][0];
    }
    double GetGTNRSyn2_TC()
    {
        return gtnrSyn_stat[1][3][1];
    }
    double GetGTNRSyn2_TG()
    {
        return gtnrSyn_stat[1][3][2];
    }
    double GetGTNRSyn2_TT()
    {
        return gtnrSyn_stat[1][3][3];
    }
    double GetGTNRSyn3_AA()
    {
        return gtnrSyn_stat[2][0][0];
    }
    double GetGTNRSyn3_AC()
    {
        return gtnrSyn_stat[2][0][1];
    }
    double GetGTNRSyn3_AG()
    {
        return gtnrSyn_stat[2][0][2];
    }
    double GetGTNRSyn3_AT()
    {
        return gtnrSyn_stat[2][0][3];
    }
    double GetGTNRSyn3_CA()
    {
        return gtnrSyn_stat[2][1][0];
    }
    double GetGTNRSyn3_CC()
    {
        return gtnrSyn_stat[2][1][1];
    }
    double GetGTNRSyn3_CG()
    {
        return gtnrSyn_stat[2][1][2];
    }
    double GetGTNRSyn3_CT()
    {
        return gtnrSyn_stat[2][1][3];
    }
    double GetGTNRSyn3_GA()
    {
        return gtnrSyn_stat[2][2][0];
    }
    double GetGTNRSyn3_GC()
    {
        return gtnrSyn_stat[2][2][1];
    }
    double GetGTNRSyn3_GG()
    {
        return gtnrSyn_stat[2][2][2];
    }
    double GetGTNRSyn3_GT()
    {
        return gtnrSyn_stat[2][2][3];
    }
    double GetGTNRSyn3_TA()
    {
        return gtnrSyn_stat[2][3][0];
    }
    double GetGTNRSyn3_TC()
    {
        return gtnrSyn_stat[2][3][1];
    }
    double GetGTNRSyn3_TG()
    {
        return gtnrSyn_stat[2][3][2];
    }
    double GetGTNRSyn3_TT()
    {
        return gtnrSyn_stat[2][3][3];
    }
    double GetGTNRNSyn_AA()
    {
        return gtnrNSyn_stat[3][0][0];
    }
    double GetGTNRNSyn_AC()
    {
        return gtnrNSyn_stat[3][0][1];
    }
    double GetGTNRNSyn_AG()
    {
        return gtnrNSyn_stat[3][0][2];
    }
    double GetGTNRNSyn_AT()
    {
        return gtnrNSyn_stat[3][0][3];
    }
    double GetGTNRNSyn_CA()
    {
        return gtnrNSyn_stat[3][1][0];
    }
    double GetGTNRNSyn_CC()
    {
        return gtnrNSyn_stat[3][1][1];
    }
    double GetGTNRNSyn_CG()
    {
        return gtnrNSyn_stat[3][1][2];
    }
    double GetGTNRNSyn_CT()
    {
        return gtnrNSyn_stat[3][1][3];
    }
    double GetGTNRNSyn_GA()
    {
        return gtnrNSyn_stat[3][2][0];
    }
    double GetGTNRNSyn_GC()
    {
        return gtnrNSyn_stat[3][2][1];
    }
    double GetGTNRNSyn_GG()
    {
        return gtnrNSyn_stat[3][2][2];
    }
    double GetGTNRNSyn_GT()
    {
        return gtnrNSyn_stat[3][2][3];
    }
    double GetGTNRNSyn_TA()
    {
        return gtnrNSyn_stat[3][3][0];
    }
    double GetGTNRNSyn_TC()
    {
        return gtnrNSyn_stat[3][3][1];
    }
    double GetGTNRNSyn_TG()
    {
        return gtnrNSyn_stat[3][3][2];
    }
    double GetGTNRNSyn_TT()
    {
        return gtnrNSyn_stat[3][3][3];
    }
    double GetGTNRNSyn1_AA()
    {
        return gtnrNSyn_stat[0][0][0];
    }
    double GetGTNRNSyn1_AC()
    {
        return gtnrNSyn_stat[0][0][1];
    }
    double GetGTNRNSyn1_AG()
    {
        return gtnrNSyn_stat[0][0][2];
    }
    double GetGTNRNSyn1_AT()
    {
        return gtnrNSyn_stat[0][0][3];
    }
    double GetGTNRNSyn1_CA()
    {
        return gtnrNSyn_stat[0][1][0];
    }
    double GetGTNRNSyn1_CC()
    {
        return gtnrNSyn_stat[0][1][1];
    }
    double GetGTNRNSyn1_CG()
    {
        return gtnrNSyn_stat[0][1][2];
    }
    double GetGTNRNSyn1_CT()
    {
        return gtnrNSyn_stat[0][1][3];
    }
    double GetGTNRNSyn1_GA()
    {
        return gtnrNSyn_stat[0][2][0];
    }
    double GetGTNRNSyn1_GC()
    {
        return gtnrNSyn_stat[0][2][1];
    }
    double GetGTNRNSyn1_GG()
    {
        return gtnrNSyn_stat[0][2][2];
    }
    double GetGTNRNSyn1_GT()
    {
        return gtnrNSyn_stat[0][2][3];
    }
    double GetGTNRNSyn1_TA()
    {
        return gtnrNSyn_stat[0][3][0];
    }
    double GetGTNRNSyn1_TC()
    {
        return gtnrNSyn_stat[0][3][1];
    }
    double GetGTNRNSyn1_TG()
    {
        return gtnrNSyn_stat[0][3][2];
    }
    double GetGTNRNSyn1_TT()
    {
        return gtnrNSyn_stat[0][3][3];
    }



    double GetGTNRNSyn2_AA()
    {
        return gtnrNSyn_stat[1][0][0];
    }
    double GetGTNRNSyn2_AC()
    {
        return gtnrNSyn_stat[1][0][1];
    }
    double GetGTNRNSyn2_AG()
    {
        return gtnrNSyn_stat[1][0][2];
    }
    double GetGTNRNSyn2_AT()
    {
        return gtnrNSyn_stat[1][0][3];
    }
    double GetGTNRNSyn2_CA()
    {
        return gtnrNSyn_stat[1][1][0];
    }
    double GetGTNRNSyn2_CC()
    {
        return gtnrNSyn_stat[1][1][1];
    }
    double GetGTNRNSyn2_CG()
    {
        return gtnrNSyn_stat[1][1][2];
    }
    double GetGTNRNSyn2_CT()
    {
        return gtnrNSyn_stat[1][1][3];
    }
    double GetGTNRNSyn2_GA()
    {
        return gtnrNSyn_stat[1][2][0];
    }
    double GetGTNRNSyn2_GC()
    {
        return gtnrNSyn_stat[1][2][1];
    }
    double GetGTNRNSyn2_GG()
    {
        return gtnrNSyn_stat[1][2][2];
    }
    double GetGTNRNSyn2_GT()
    {
        return gtnrNSyn_stat[1][2][3];
    }
    double GetGTNRNSyn2_TA()
    {
        return gtnrNSyn_stat[1][3][0];
    }
    double GetGTNRNSyn2_TC()
    {
        return gtnrNSyn_stat[1][3][1];
    }
    double GetGTNRNSyn2_TG()
    {
        return gtnrNSyn_stat[1][3][2];
    }
    double GetGTNRNSyn2_TT()
    {
        return gtnrNSyn_stat[1][3][3];
    }
    double GetGTNRNSyn3_AA()
    {
        return gtnrNSyn_stat[2][0][0];
    }
    double GetGTNRNSyn3_AC()
    {
        return gtnrNSyn_stat[2][0][1];
    }
    double GetGTNRNSyn3_AG()
    {
        return gtnrNSyn_stat[2][0][2];
    }
    double GetGTNRNSyn3_AT()
    {
        return gtnrNSyn_stat[2][0][3];
    }
    double GetGTNRNSyn3_CA()
    {
        return gtnrNSyn_stat[2][1][0];
    }
    double GetGTNRNSyn3_CC()
    {
        return gtnrNSyn_stat[2][1][1];
    }
    double GetGTNRNSyn3_CG()
    {
        return gtnrNSyn_stat[2][1][2];
    }
    double GetGTNRNSyn3_CT()
    {
        return gtnrNSyn_stat[2][1][3];
    }
    double GetGTNRNSyn3_GA()
    {
        return gtnrNSyn_stat[2][2][0];
    }
    double GetGTNRNSyn3_GC()
    {
        return gtnrNSyn_stat[2][2][1];
    }
    double GetGTNRNSyn3_GG()
    {
        return gtnrNSyn_stat[2][2][2];
    }
    double GetGTNRNSyn3_GT()
    {
        return gtnrNSyn_stat[2][2][3];
    }
    double GetGTNRNSyn3_TA()
    {
        return gtnrNSyn_stat[2][3][0];
    }
    double GetGTNRNSyn3_TC()
    {
        return gtnrNSyn_stat[2][3][1];
    }
    double GetGTNRNSyn3_TG()
    {
        return gtnrNSyn_stat[2][3][2];
    }
    double GetGTNRNSyn3_TT()
    {
        return gtnrNSyn_stat[2][3][3];
    }

    double GetDinuc_AAAA()
    {
        return dinuc_stat[3][0][0];
    }
    double GetDinuc_AAAC()
    {
        return dinuc_stat[3][0][1];
    }
    double GetDinuc_AAAG()
    {
        return dinuc_stat[3][0][2];
    }
    double GetDinuc_AAAT()
    {
        return dinuc_stat[3][0][3];
    }
    double GetDinuc_AACA()
    {
        return dinuc_stat[3][0][4];
    }
    double GetDinuc_AACC()
    {
        return dinuc_stat[3][0][5];
    }
    double GetDinuc_AACG()
    {
        return dinuc_stat[3][0][6];
    }
    double GetDinuc_AACT()
    {
        return dinuc_stat[3][0][7];
    }
    double GetDinuc_AAGA()
    {
        return dinuc_stat[3][0][8];
    }
    double GetDinuc_AAGC()
    {
        return dinuc_stat[3][0][9];
    }
    double GetDinuc_AAGG()
    {
        return dinuc_stat[3][0][10];
    }
    double GetDinuc_AAGT()
    {
        return dinuc_stat[3][0][11];
    }
    double GetDinuc_AATA()
    {
        return dinuc_stat[3][0][12];
    }
    double GetDinuc_AATC()
    {
        return dinuc_stat[3][0][13];
    }
    double GetDinuc_AATG()
    {
        return dinuc_stat[3][0][14];
    }
    double GetDinuc_AATT()
    {
        return dinuc_stat[3][0][15];
    }
    double GetDinuc_ACAA()
    {
        return dinuc_stat[3][1][0];
    }
    double GetDinuc_ACAC()
    {
        return dinuc_stat[3][1][1];
    }
    double GetDinuc_ACAG()
    {
        return dinuc_stat[3][1][2];
    }
    double GetDinuc_ACAT()
    {
        return dinuc_stat[3][1][3];
    }
    double GetDinuc_ACCA()
    {
        return dinuc_stat[3][1][4];
    }
    double GetDinuc_ACCC()
    {
        return dinuc_stat[3][1][5];
    }
    double GetDinuc_ACCG()
    {
        return dinuc_stat[3][1][6];
    }
    double GetDinuc_ACCT()
    {
        return dinuc_stat[3][1][7];
    }
    double GetDinuc_ACGA()
    {
        return dinuc_stat[3][1][8];
    }
    double GetDinuc_ACGC()
    {
        return dinuc_stat[3][1][9];
    }
    double GetDinuc_ACGG()
    {
        return dinuc_stat[3][1][10];
    }
    double GetDinuc_ACGT()
    {
        return dinuc_stat[3][1][11];
    }
    double GetDinuc_ACTA()
    {
        return dinuc_stat[3][1][12];
    }
    double GetDinuc_ACTC()
    {
        return dinuc_stat[3][1][13];
    }
    double GetDinuc_ACTG()
    {
        return dinuc_stat[3][1][14];
    }
    double GetDinuc_ACTT()
    {
        return dinuc_stat[3][1][15];
    }
    double GetDinuc_AGAA()
    {
        return dinuc_stat[3][2][0];
    }
    double GetDinuc_AGAC()
    {
        return dinuc_stat[3][2][1];
    }
    double GetDinuc_AGAG()
    {
        return dinuc_stat[3][2][2];
    }
    double GetDinuc_AGAT()
    {
        return dinuc_stat[3][2][3];
    }
    double GetDinuc_AGCA()
    {
        return dinuc_stat[3][2][4];
    }
    double GetDinuc_AGCC()
    {
        return dinuc_stat[3][2][5];
    }
    double GetDinuc_AGCG()
    {
        return dinuc_stat[3][2][6];
    }
    double GetDinuc_AGCT()
    {
        return dinuc_stat[3][2][7];
    }
    double GetDinuc_AGGA()
    {
        return dinuc_stat[3][2][8];
    }
    double GetDinuc_AGGC()
    {
        return dinuc_stat[3][2][9];
    }
    double GetDinuc_AGGG()
    {
        return dinuc_stat[3][2][10];
    }
    double GetDinuc_AGGT()
    {
        return dinuc_stat[3][2][11];
    }
    double GetDinuc_AGTA()
    {
        return dinuc_stat[3][2][12];
    }
    double GetDinuc_AGTC()
    {
        return dinuc_stat[3][2][13];
    }
    double GetDinuc_AGTG()
    {
        return dinuc_stat[3][2][14];
    }
    double GetDinuc_AGTT()
    {
        return dinuc_stat[3][2][15];
    }
    double GetDinuc_ATAA()
    {
        return dinuc_stat[3][3][0];
    }
    double GetDinuc_ATAC()
    {
        return dinuc_stat[3][3][1];
    }
    double GetDinuc_ATAG()
    {
        return dinuc_stat[3][3][2];
    }
    double GetDinuc_ATAT()
    {
        return dinuc_stat[3][3][3];
    }
    double GetDinuc_ATCA()
    {
        return dinuc_stat[3][3][4];
    }
    double GetDinuc_ATCC()
    {
        return dinuc_stat[3][3][5];
    }
    double GetDinuc_ATCG()
    {
        return dinuc_stat[3][3][6];
    }
    double GetDinuc_ATCT()
    {
        return dinuc_stat[3][3][7];
    }
    double GetDinuc_ATGA()
    {
        return dinuc_stat[3][3][8];
    }
    double GetDinuc_ATGC()
    {
        return dinuc_stat[3][3][9];
    }
    double GetDinuc_ATGG()
    {
        return dinuc_stat[3][3][10];
    }
    double GetDinuc_ATGT()
    {
        return dinuc_stat[3][3][11];
    }
    double GetDinuc_ATTA()
    {
        return dinuc_stat[3][3][12];
    }
    double GetDinuc_ATTC()
    {
        return dinuc_stat[3][3][13];
    }
    double GetDinuc_ATTG()
    {
        return dinuc_stat[3][3][14];
    }
    double GetDinuc_ATTT()
    {
        return dinuc_stat[3][3][15];
    }
    double GetDinuc_CAAA()
    {
        return dinuc_stat[3][4][0];
    }
    double GetDinuc_CAAC()
    {
        return dinuc_stat[3][4][1];
    }
    double GetDinuc_CAAG()
    {
        return dinuc_stat[3][4][2];
    }
    double GetDinuc_CAAT()
    {
        return dinuc_stat[3][4][3];
    }
    double GetDinuc_CACA()
    {
        return dinuc_stat[3][4][4];
    }
    double GetDinuc_CACC()
    {
        return dinuc_stat[3][4][5];
    }
    double GetDinuc_CACG()
    {
        return dinuc_stat[3][4][6];
    }
    double GetDinuc_CACT()
    {
        return dinuc_stat[3][4][7];
    }
    double GetDinuc_CAGA()
    {
        return dinuc_stat[3][4][8];
    }
    double GetDinuc_CAGC()
    {
        return dinuc_stat[3][4][9];
    }
    double GetDinuc_CAGG()
    {
        return dinuc_stat[3][4][10];
    }
    double GetDinuc_CAGT()
    {
        return dinuc_stat[3][4][11];
    }
    double GetDinuc_CATA()
    {
        return dinuc_stat[3][4][12];
    }
    double GetDinuc_CATC()
    {
        return dinuc_stat[3][4][13];
    }
    double GetDinuc_CATG()
    {
        return dinuc_stat[3][4][14];
    }
    double GetDinuc_CATT()
    {
        return dinuc_stat[3][4][15];
    }
    double GetDinuc_CCAA()
    {
        return dinuc_stat[3][5][0];
    }
    double GetDinuc_CCAC()
    {
        return dinuc_stat[3][5][1];
    }
    double GetDinuc_CCAG()
    {
        return dinuc_stat[3][5][2];
    }
    double GetDinuc_CCAT()
    {
        return dinuc_stat[3][5][3];
    }
    double GetDinuc_CCCA()
    {
        return dinuc_stat[3][5][4];
    }
    double GetDinuc_CCCC()
    {
        return dinuc_stat[3][5][5];
    }
    double GetDinuc_CCCG()
    {
        return dinuc_stat[3][5][6];
    }
    double GetDinuc_CCCT()
    {
        return dinuc_stat[3][5][7];
    }
    double GetDinuc_CCGA()
    {
        return dinuc_stat[3][5][8];
    }
    double GetDinuc_CCGC()
    {
        return dinuc_stat[3][5][9];
    }
    double GetDinuc_CCGG()
    {
        return dinuc_stat[3][5][10];
    }
    double GetDinuc_CCGT()
    {
        return dinuc_stat[3][5][11];
    }
    double GetDinuc_CCTA()
    {
        return dinuc_stat[3][5][12];
    }
    double GetDinuc_CCTC()
    {
        return dinuc_stat[3][5][13];
    }
    double GetDinuc_CCTG()
    {
        return dinuc_stat[3][5][14];
    }
    double GetDinuc_CCTT()
    {
        return dinuc_stat[3][5][15];
    }
    double GetDinuc_CGAA()
    {
        return dinuc_stat[3][6][0];
    }
    double GetDinuc_CGAC()
    {
        return dinuc_stat[3][6][1];
    }
    double GetDinuc_CGAG()
    {
        return dinuc_stat[3][6][2];
    }
    double GetDinuc_CGAT()
    {
        return dinuc_stat[3][6][3];
    }
    double GetDinuc_CGCA()
    {
        return dinuc_stat[3][6][4];
    }
    double GetDinuc_CGCC()
    {
        return dinuc_stat[3][6][5];
    }
    double GetDinuc_CGCG()
    {
        return dinuc_stat[3][6][6];
    }
    double GetDinuc_CGCT()
    {
        return dinuc_stat[3][6][7];
    }
    double GetDinuc_CGGA()
    {
        return dinuc_stat[3][6][8];
    }
    double GetDinuc_CGGC()
    {
        return dinuc_stat[3][6][9];
    }
    double GetDinuc_CGGG()
    {
        return dinuc_stat[3][6][10];
    }
    double GetDinuc_CGGT()
    {
        return dinuc_stat[3][6][11];
    }
    double GetDinuc_CGTA()
    {
        return dinuc_stat[3][6][12];
    }
    double GetDinuc_CGTC()
    {
        return dinuc_stat[3][6][13];
    }
    double GetDinuc_CGTG()
    {
        return dinuc_stat[3][6][14];
    }
    double GetDinuc_CGTT()
    {
        return dinuc_stat[3][6][15];
    }
    double GetDinuc_CTAA()
    {
        return dinuc_stat[3][7][0];
    }
    double GetDinuc_CTAC()
    {
        return dinuc_stat[3][7][1];
    }
    double GetDinuc_CTAG()
    {
        return dinuc_stat[3][7][2];
    }
    double GetDinuc_CTAT()
    {
        return dinuc_stat[3][7][3];
    }
    double GetDinuc_CTCA()
    {
        return dinuc_stat[3][7][4];
    }
    double GetDinuc_CTCC()
    {
        return dinuc_stat[3][7][5];
    }
    double GetDinuc_CTCG()
    {
        return dinuc_stat[3][7][6];
    }
    double GetDinuc_CTCT()
    {
        return dinuc_stat[3][7][7];
    }
    double GetDinuc_CTGA()
    {
        return dinuc_stat[3][7][8];
    }
    double GetDinuc_CTGC()
    {
        return dinuc_stat[3][7][9];
    }
    double GetDinuc_CTGG()
    {
        return dinuc_stat[3][7][10];
    }
    double GetDinuc_CTGT()
    {
        return dinuc_stat[3][7][11];
    }
    double GetDinuc_CTTA()
    {
        return dinuc_stat[3][7][12];
    }
    double GetDinuc_CTTC()
    {
        return dinuc_stat[3][7][13];
    }
    double GetDinuc_CTTG()
    {
        return dinuc_stat[3][7][14];
    }
    double GetDinuc_CTTT()
    {
        return dinuc_stat[3][7][15];
    }
    double GetDinuc_GAAA()
    {
        return dinuc_stat[3][8][0];
    }
    double GetDinuc_GAAC()
    {
        return dinuc_stat[3][8][1];
    }
    double GetDinuc_GAAG()
    {
        return dinuc_stat[3][8][2];
    }
    double GetDinuc_GAAT()
    {
        return dinuc_stat[3][8][3];
    }
    double GetDinuc_GACA()
    {
        return dinuc_stat[3][8][4];
    }
    double GetDinuc_GACC()
    {
        return dinuc_stat[3][8][5];
    }
    double GetDinuc_GACG()
    {
        return dinuc_stat[3][8][6];
    }
    double GetDinuc_GACT()
    {
        return dinuc_stat[3][8][7];
    }
    double GetDinuc_GAGA()
    {
        return dinuc_stat[3][8][8];
    }
    double GetDinuc_GAGC()
    {
        return dinuc_stat[3][8][9];
    }
    double GetDinuc_GAGG()
    {
        return dinuc_stat[3][8][10];
    }
    double GetDinuc_GAGT()
    {
        return dinuc_stat[3][8][11];
    }
    double GetDinuc_GATA()
    {
        return dinuc_stat[3][8][12];
    }
    double GetDinuc_GATC()
    {
        return dinuc_stat[3][8][13];
    }
    double GetDinuc_GATG()
    {
        return dinuc_stat[3][8][14];
    }
    double GetDinuc_GATT()
    {
        return dinuc_stat[3][8][15];
    }
    double GetDinuc_GCAA()
    {
        return dinuc_stat[3][9][0];
    }
    double GetDinuc_GCAC()
    {
        return dinuc_stat[3][9][1];
    }
    double GetDinuc_GCAG()
    {
        return dinuc_stat[3][9][2];
    }
    double GetDinuc_GCAT()
    {
        return dinuc_stat[3][9][3];
    }
    double GetDinuc_GCCA()
    {
        return dinuc_stat[3][9][4];
    }
    double GetDinuc_GCCC()
    {
        return dinuc_stat[3][9][5];
    }
    double GetDinuc_GCCG()
    {
        return dinuc_stat[3][9][6];
    }
    double GetDinuc_GCCT()
    {
        return dinuc_stat[3][9][7];
    }
    double GetDinuc_GCGA()
    {
        return dinuc_stat[3][9][8];
    }
    double GetDinuc_GCGC()
    {
        return dinuc_stat[3][9][9];
    }
    double GetDinuc_GCGG()
    {
        return dinuc_stat[3][9][10];
    }
    double GetDinuc_GCGT()
    {
        return dinuc_stat[3][9][11];
    }
    double GetDinuc_GCTA()
    {
        return dinuc_stat[3][9][12];
    }
    double GetDinuc_GCTC()
    {
        return dinuc_stat[3][9][13];
    }
    double GetDinuc_GCTG()
    {
        return dinuc_stat[3][9][14];
    }
    double GetDinuc_GCTT()
    {
        return dinuc_stat[3][9][15];
    }
    double GetDinuc_GGAA()
    {
        return dinuc_stat[3][10][0];
    }
    double GetDinuc_GGAC()
    {
        return dinuc_stat[3][10][1];
    }
    double GetDinuc_GGAG()
    {
        return dinuc_stat[3][10][2];
    }
    double GetDinuc_GGAT()
    {
        return dinuc_stat[3][10][3];
    }
    double GetDinuc_GGCA()
    {
        return dinuc_stat[3][10][4];
    }
    double GetDinuc_GGCC()
    {
        return dinuc_stat[3][10][5];
    }
    double GetDinuc_GGCG()
    {
        return dinuc_stat[3][10][6];
    }
    double GetDinuc_GGCT()
    {
        return dinuc_stat[3][10][7];
    }
    double GetDinuc_GGGA()
    {
        return dinuc_stat[3][10][8];
    }
    double GetDinuc_GGGC()
    {
        return dinuc_stat[3][10][9];
    }
    double GetDinuc_GGGG()
    {
        return dinuc_stat[3][10][10];
    }
    double GetDinuc_GGGT()
    {
        return dinuc_stat[3][10][11];
    }
    double GetDinuc_GGTA()
    {
        return dinuc_stat[3][10][12];
    }
    double GetDinuc_GGTC()
    {
        return dinuc_stat[3][10][13];
    }
    double GetDinuc_GGTG()
    {
        return dinuc_stat[3][10][14];
    }
    double GetDinuc_GGTT()
    {
        return dinuc_stat[3][10][15];
    }
    double GetDinuc_GTAA()
    {
        return dinuc_stat[3][11][0];
    }
    double GetDinuc_GTAC()
    {
        return dinuc_stat[3][11][1];
    }
    double GetDinuc_GTAG()
    {
        return dinuc_stat[3][11][2];
    }
    double GetDinuc_GTAT()
    {
        return dinuc_stat[3][11][3];
    }
    double GetDinuc_GTCA()
    {
        return dinuc_stat[3][11][4];
    }
    double GetDinuc_GTCC()
    {
        return dinuc_stat[3][11][5];
    }
    double GetDinuc_GTCG()
    {
        return dinuc_stat[3][11][6];
    }
    double GetDinuc_GTCT()
    {
        return dinuc_stat[3][11][7];
    }
    double GetDinuc_GTGA()
    {
        return dinuc_stat[3][11][8];
    }
    double GetDinuc_GTGC()
    {
        return dinuc_stat[3][11][9];
    }
    double GetDinuc_GTGG()
    {
        return dinuc_stat[3][11][10];
    }
    double GetDinuc_GTGT()
    {
        return dinuc_stat[3][11][11];
    }
    double GetDinuc_GTTA()
    {
        return dinuc_stat[3][11][12];
    }
    double GetDinuc_GTTC()
    {
        return dinuc_stat[3][11][13];
    }
    double GetDinuc_GTTG()
    {
        return dinuc_stat[3][11][14];
    }
    double GetDinuc_GTTT()
    {
        return dinuc_stat[3][11][15];
    }
    double GetDinuc_TAAA()
    {
        return dinuc_stat[3][12][0];
    }
    double GetDinuc_TAAC()
    {
        return dinuc_stat[3][12][1];
    }
    double GetDinuc_TAAG()
    {
        return dinuc_stat[3][12][2];
    }
    double GetDinuc_TAAT()
    {
        return dinuc_stat[3][12][3];
    }
    double GetDinuc_TACA()
    {
        return dinuc_stat[3][12][4];
    }
    double GetDinuc_TACC()
    {
        return dinuc_stat[3][12][5];
    }
    double GetDinuc_TACG()
    {
        return dinuc_stat[3][12][6];
    }
    double GetDinuc_TACT()
    {
        return dinuc_stat[3][12][7];
    }
    double GetDinuc_TAGA()
    {
        return dinuc_stat[3][12][8];
    }
    double GetDinuc_TAGC()
    {
        return dinuc_stat[3][12][9];
    }
    double GetDinuc_TAGG()
    {
        return dinuc_stat[3][12][10];
    }
    double GetDinuc_TAGT()
    {
        return dinuc_stat[3][12][11];
    }
    double GetDinuc_TATA()
    {
        return dinuc_stat[3][12][12];
    }
    double GetDinuc_TATC()
    {
        return dinuc_stat[3][12][13];
    }
    double GetDinuc_TATG()
    {
        return dinuc_stat[3][12][14];
    }
    double GetDinuc_TATT()
    {
        return dinuc_stat[3][12][15];
    }
    double GetDinuc_TCAA()
    {
        return dinuc_stat[3][13][0];
    }
    double GetDinuc_TCAC()
    {
        return dinuc_stat[3][13][1];
    }
    double GetDinuc_TCAG()
    {
        return dinuc_stat[3][13][2];
    }
    double GetDinuc_TCAT()
    {
        return dinuc_stat[3][13][3];
    }
    double GetDinuc_TCCA()
    {
        return dinuc_stat[3][13][4];
    }
    double GetDinuc_TCCC()
    {
        return dinuc_stat[3][13][5];
    }
    double GetDinuc_TCCG()
    {
        return dinuc_stat[3][13][6];
    }
    double GetDinuc_TCCT()
    {
        return dinuc_stat[3][13][7];
    }
    double GetDinuc_TCGA()
    {
        return dinuc_stat[3][13][8];
    }
    double GetDinuc_TCGC()
    {
        return dinuc_stat[3][13][9];
    }
    double GetDinuc_TCGG()
    {
        return dinuc_stat[3][13][10];
    }
    double GetDinuc_TCGT()
    {
        return dinuc_stat[3][13][11];
    }
    double GetDinuc_TCTA()
    {
        return dinuc_stat[3][13][12];
    }
    double GetDinuc_TCTC()
    {
        return dinuc_stat[3][13][13];
    }
    double GetDinuc_TCTG()
    {
        return dinuc_stat[3][13][14];
    }
    double GetDinuc_TCTT()
    {
        return dinuc_stat[3][13][15];
    }
    double GetDinuc_TGAA()
    {
        return dinuc_stat[3][14][0];
    }
    double GetDinuc_TGAC()
    {
        return dinuc_stat[3][14][1];
    }
    double GetDinuc_TGAG()
    {
        return dinuc_stat[3][14][2];
    }
    double GetDinuc_TGAT()
    {
        return dinuc_stat[3][14][3];
    }
    double GetDinuc_TGCA()
    {
        return dinuc_stat[3][14][4];
    }
    double GetDinuc_TGCC()
    {
        return dinuc_stat[3][14][5];
    }
    double GetDinuc_TGCG()
    {
        return dinuc_stat[3][14][6];
    }
    double GetDinuc_TGCT()
    {
        return dinuc_stat[3][14][7];
    }
    double GetDinuc_TGGA()
    {
        return dinuc_stat[3][14][8];
    }
    double GetDinuc_TGGC()
    {
        return dinuc_stat[3][14][9];
    }
    double GetDinuc_TGGG()
    {
        return dinuc_stat[3][14][10];
    }
    double GetDinuc_TGGT()
    {
        return dinuc_stat[3][14][11];
    }
    double GetDinuc_TGTA()
    {
        return dinuc_stat[3][14][12];
    }
    double GetDinuc_TGTC()
    {
        return dinuc_stat[3][14][13];
    }
    double GetDinuc_TGTG()
    {
        return dinuc_stat[3][14][14];
    }
    double GetDinuc_TGTT()
    {
        return dinuc_stat[3][14][15];
    }
    double GetDinuc_TTAA()
    {
        return dinuc_stat[3][15][0];
    }
    double GetDinuc_TTAC()
    {
        return dinuc_stat[3][15][1];
    }
    double GetDinuc_TTAG()
    {
        return dinuc_stat[3][15][2];
    }
    double GetDinuc_TTAT()
    {
        return dinuc_stat[3][15][3];
    }
    double GetDinuc_TTCA()
    {
        return dinuc_stat[3][15][4];
    }
    double GetDinuc_TTCC()
    {
        return dinuc_stat[3][15][5];
    }
    double GetDinuc_TTCG()
    {
        return dinuc_stat[3][15][6];
    }
    double GetDinuc_TTCT()
    {
        return dinuc_stat[3][15][7];
    }
    double GetDinuc_TTGA()
    {
        return dinuc_stat[3][15][8];
    }
    double GetDinuc_TTGC()
    {
        return dinuc_stat[3][15][9];
    }
    double GetDinuc_TTGG()
    {
        return dinuc_stat[3][15][10];
    }
    double GetDinuc_TTGT()
    {
        return dinuc_stat[3][15][11];
    }
    double GetDinuc_TTTA()
    {
        return dinuc_stat[3][15][12];
    }
    double GetDinuc_TTTC()
    {
        return dinuc_stat[3][15][13];
    }
    double GetDinuc_TTTG()
    {
        return dinuc_stat[3][15][14];
    }
    double GetDinuc_TTTT()
    {
        return dinuc_stat[3][15][15];
    }
    double GetDinuc12_AAAA()
    {
        return dinuc_stat[0][0][0];
    }
    double GetDinuc12_AAAC()
    {
        return dinuc_stat[0][0][1];
    }
    double GetDinuc12_AAAG()
    {
        return dinuc_stat[0][0][2];
    }
    double GetDinuc12_AAAT()
    {
        return dinuc_stat[0][0][3];
    }
    double GetDinuc12_AACA()
    {
        return dinuc_stat[0][0][4];
    }
    double GetDinuc12_AACC()
    {
        return dinuc_stat[0][0][5];
    }
    double GetDinuc12_AACG()
    {
        return dinuc_stat[0][0][6];
    }
    double GetDinuc12_AACT()
    {
        return dinuc_stat[0][0][7];
    }
    double GetDinuc12_AAGA()
    {
        return dinuc_stat[0][0][8];
    }
    double GetDinuc12_AAGC()
    {
        return dinuc_stat[0][0][9];
    }
    double GetDinuc12_AAGG()
    {
        return dinuc_stat[0][0][10];
    }
    double GetDinuc12_AAGT()
    {
        return dinuc_stat[0][0][11];
    }
    double GetDinuc12_AATA()
    {
        return dinuc_stat[0][0][12];
    }
    double GetDinuc12_AATC()
    {
        return dinuc_stat[0][0][13];
    }
    double GetDinuc12_AATG()
    {
        return dinuc_stat[0][0][14];
    }
    double GetDinuc12_AATT()
    {
        return dinuc_stat[0][0][15];
    }
    double GetDinuc12_ACAA()
    {
        return dinuc_stat[0][1][0];
    }
    double GetDinuc12_ACAC()
    {
        return dinuc_stat[0][1][1];
    }
    double GetDinuc12_ACAG()
    {
        return dinuc_stat[0][1][2];
    }
    double GetDinuc12_ACAT()
    {
        return dinuc_stat[0][1][3];
    }
    double GetDinuc12_ACCA()
    {
        return dinuc_stat[0][1][4];
    }
    double GetDinuc12_ACCC()
    {
        return dinuc_stat[0][1][5];
    }
    double GetDinuc12_ACCG()
    {
        return dinuc_stat[0][1][6];
    }
    double GetDinuc12_ACCT()
    {
        return dinuc_stat[0][1][7];
    }
    double GetDinuc12_ACGA()
    {
        return dinuc_stat[0][1][8];
    }
    double GetDinuc12_ACGC()
    {
        return dinuc_stat[0][1][9];
    }
    double GetDinuc12_ACGG()
    {
        return dinuc_stat[0][1][10];
    }
    double GetDinuc12_ACGT()
    {
        return dinuc_stat[0][1][11];
    }
    double GetDinuc12_ACTA()
    {
        return dinuc_stat[0][1][12];
    }
    double GetDinuc12_ACTC()
    {
        return dinuc_stat[0][1][13];
    }
    double GetDinuc12_ACTG()
    {
        return dinuc_stat[0][1][14];
    }
    double GetDinuc12_ACTT()
    {
        return dinuc_stat[0][1][15];
    }
    double GetDinuc12_AGAA()
    {
        return dinuc_stat[0][2][0];
    }
    double GetDinuc12_AGAC()
    {
        return dinuc_stat[0][2][1];
    }
    double GetDinuc12_AGAG()
    {
        return dinuc_stat[0][2][2];
    }
    double GetDinuc12_AGAT()
    {
        return dinuc_stat[0][2][3];
    }
    double GetDinuc12_AGCA()
    {
        return dinuc_stat[0][2][4];
    }
    double GetDinuc12_AGCC()
    {
        return dinuc_stat[0][2][5];
    }
    double GetDinuc12_AGCG()
    {
        return dinuc_stat[0][2][6];
    }
    double GetDinuc12_AGCT()
    {
        return dinuc_stat[0][2][7];
    }
    double GetDinuc12_AGGA()
    {
        return dinuc_stat[0][2][8];
    }
    double GetDinuc12_AGGC()
    {
        return dinuc_stat[0][2][9];
    }
    double GetDinuc12_AGGG()
    {
        return dinuc_stat[0][2][10];
    }
    double GetDinuc12_AGGT()
    {
        return dinuc_stat[0][2][11];
    }
    double GetDinuc12_AGTA()
    {
        return dinuc_stat[0][2][12];
    }
    double GetDinuc12_AGTC()
    {
        return dinuc_stat[0][2][13];
    }
    double GetDinuc12_AGTG()
    {
        return dinuc_stat[0][2][14];
    }
    double GetDinuc12_AGTT()
    {
        return dinuc_stat[0][2][15];
    }
    double GetDinuc12_ATAA()
    {
        return dinuc_stat[0][3][0];
    }
    double GetDinuc12_ATAC()
    {
        return dinuc_stat[0][3][1];
    }
    double GetDinuc12_ATAG()
    {
        return dinuc_stat[0][3][2];
    }
    double GetDinuc12_ATAT()
    {
        return dinuc_stat[0][3][3];
    }
    double GetDinuc12_ATCA()
    {
        return dinuc_stat[0][3][4];
    }
    double GetDinuc12_ATCC()
    {
        return dinuc_stat[0][3][5];
    }
    double GetDinuc12_ATCG()
    {
        return dinuc_stat[0][3][6];
    }
    double GetDinuc12_ATCT()
    {
        return dinuc_stat[0][3][7];
    }
    double GetDinuc12_ATGA()
    {
        return dinuc_stat[0][3][8];
    }
    double GetDinuc12_ATGC()
    {
        return dinuc_stat[0][3][9];
    }
    double GetDinuc12_ATGG()
    {
        return dinuc_stat[0][3][10];
    }
    double GetDinuc12_ATGT()
    {
        return dinuc_stat[0][3][11];
    }
    double GetDinuc12_ATTA()
    {
        return dinuc_stat[0][3][12];
    }
    double GetDinuc12_ATTC()
    {
        return dinuc_stat[0][3][13];
    }
    double GetDinuc12_ATTG()
    {
        return dinuc_stat[0][3][14];
    }
    double GetDinuc12_ATTT()
    {
        return dinuc_stat[0][3][15];
    }
    double GetDinuc12_CAAA()
    {
        return dinuc_stat[0][4][0];
    }
    double GetDinuc12_CAAC()
    {
        return dinuc_stat[0][4][1];
    }
    double GetDinuc12_CAAG()
    {
        return dinuc_stat[0][4][2];
    }
    double GetDinuc12_CAAT()
    {
        return dinuc_stat[0][4][3];
    }
    double GetDinuc12_CACA()
    {
        return dinuc_stat[0][4][4];
    }
    double GetDinuc12_CACC()
    {
        return dinuc_stat[0][4][5];
    }
    double GetDinuc12_CACG()
    {
        return dinuc_stat[0][4][6];
    }
    double GetDinuc12_CACT()
    {
        return dinuc_stat[0][4][7];
    }
    double GetDinuc12_CAGA()
    {
        return dinuc_stat[0][4][8];
    }
    double GetDinuc12_CAGC()
    {
        return dinuc_stat[0][4][9];
    }
    double GetDinuc12_CAGG()
    {
        return dinuc_stat[0][4][10];
    }
    double GetDinuc12_CAGT()
    {
        return dinuc_stat[0][4][11];
    }
    double GetDinuc12_CATA()
    {
        return dinuc_stat[0][4][12];
    }
    double GetDinuc12_CATC()
    {
        return dinuc_stat[0][4][13];
    }
    double GetDinuc12_CATG()
    {
        return dinuc_stat[0][4][14];
    }
    double GetDinuc12_CATT()
    {
        return dinuc_stat[0][4][15];
    }
    double GetDinuc12_CCAA()
    {
        return dinuc_stat[0][5][0];
    }
    double GetDinuc12_CCAC()
    {
        return dinuc_stat[0][5][1];
    }
    double GetDinuc12_CCAG()
    {
        return dinuc_stat[0][5][2];
    }
    double GetDinuc12_CCAT()
    {
        return dinuc_stat[0][5][3];
    }
    double GetDinuc12_CCCA()
    {
        return dinuc_stat[0][5][4];
    }
    double GetDinuc12_CCCC()
    {
        return dinuc_stat[0][5][5];
    }
    double GetDinuc12_CCCG()
    {
        return dinuc_stat[0][5][6];
    }
    double GetDinuc12_CCCT()
    {
        return dinuc_stat[0][5][7];
    }
    double GetDinuc12_CCGA()
    {
        return dinuc_stat[0][5][8];
    }
    double GetDinuc12_CCGC()
    {
        return dinuc_stat[0][5][9];
    }
    double GetDinuc12_CCGG()
    {
        return dinuc_stat[0][5][10];
    }
    double GetDinuc12_CCGT()
    {
        return dinuc_stat[0][5][11];
    }
    double GetDinuc12_CCTA()
    {
        return dinuc_stat[0][5][12];
    }
    double GetDinuc12_CCTC()
    {
        return dinuc_stat[0][5][13];
    }
    double GetDinuc12_CCTG()
    {
        return dinuc_stat[0][5][14];
    }
    double GetDinuc12_CCTT()
    {
        return dinuc_stat[0][5][15];
    }
    double GetDinuc12_CGAA()
    {
        return dinuc_stat[0][6][0];
    }
    double GetDinuc12_CGAC()
    {
        return dinuc_stat[0][6][1];
    }
    double GetDinuc12_CGAG()
    {
        return dinuc_stat[0][6][2];
    }
    double GetDinuc12_CGAT()
    {
        return dinuc_stat[0][6][3];
    }
    double GetDinuc12_CGCA()
    {
        return dinuc_stat[0][6][4];
    }
    double GetDinuc12_CGCC()
    {
        return dinuc_stat[0][6][5];
    }
    double GetDinuc12_CGCG()
    {
        return dinuc_stat[0][6][6];
    }
    double GetDinuc12_CGCT()
    {
        return dinuc_stat[0][6][7];
    }
    double GetDinuc12_CGGA()
    {
        return dinuc_stat[0][6][8];
    }
    double GetDinuc12_CGGC()
    {
        return dinuc_stat[0][6][9];
    }
    double GetDinuc12_CGGG()
    {
        return dinuc_stat[0][6][10];
    }
    double GetDinuc12_CGGT()
    {
        return dinuc_stat[0][6][11];
    }
    double GetDinuc12_CGTA()
    {
        return dinuc_stat[0][6][12];
    }
    double GetDinuc12_CGTC()
    {
        return dinuc_stat[0][6][13];
    }
    double GetDinuc12_CGTG()
    {
        return dinuc_stat[0][6][14];
    }
    double GetDinuc12_CGTT()
    {
        return dinuc_stat[0][6][15];
    }
    double GetDinuc12_CTAA()
    {
        return dinuc_stat[0][7][0];
    }
    double GetDinuc12_CTAC()
    {
        return dinuc_stat[0][7][1];
    }
    double GetDinuc12_CTAG()
    {
        return dinuc_stat[0][7][2];
    }
    double GetDinuc12_CTAT()
    {
        return dinuc_stat[0][7][3];
    }
    double GetDinuc12_CTCA()
    {
        return dinuc_stat[0][7][4];
    }
    double GetDinuc12_CTCC()
    {
        return dinuc_stat[0][7][5];
    }
    double GetDinuc12_CTCG()
    {
        return dinuc_stat[0][7][6];
    }
    double GetDinuc12_CTCT()
    {
        return dinuc_stat[0][7][7];
    }
    double GetDinuc12_CTGA()
    {
        return dinuc_stat[0][7][8];
    }
    double GetDinuc12_CTGC()
    {
        return dinuc_stat[0][7][9];
    }
    double GetDinuc12_CTGG()
    {
        return dinuc_stat[0][7][10];
    }
    double GetDinuc12_CTGT()
    {
        return dinuc_stat[0][7][11];
    }
    double GetDinuc12_CTTA()
    {
        return dinuc_stat[0][7][12];
    }
    double GetDinuc12_CTTC()
    {
        return dinuc_stat[0][7][13];
    }
    double GetDinuc12_CTTG()
    {
        return dinuc_stat[0][7][14];
    }
    double GetDinuc12_CTTT()
    {
        return dinuc_stat[0][7][15];
    }
    double GetDinuc12_GAAA()
    {
        return dinuc_stat[0][8][0];
    }
    double GetDinuc12_GAAC()
    {
        return dinuc_stat[0][8][1];
    }
    double GetDinuc12_GAAG()
    {
        return dinuc_stat[0][8][2];
    }
    double GetDinuc12_GAAT()
    {
        return dinuc_stat[0][8][3];
    }
    double GetDinuc12_GACA()
    {
        return dinuc_stat[0][8][4];
    }
    double GetDinuc12_GACC()
    {
        return dinuc_stat[0][8][5];
    }
    double GetDinuc12_GACG()
    {
        return dinuc_stat[0][8][6];
    }
    double GetDinuc12_GACT()
    {
        return dinuc_stat[0][8][7];
    }
    double GetDinuc12_GAGA()
    {
        return dinuc_stat[0][8][8];
    }
    double GetDinuc12_GAGC()
    {
        return dinuc_stat[0][8][9];
    }
    double GetDinuc12_GAGG()
    {
        return dinuc_stat[0][8][10];
    }
    double GetDinuc12_GAGT()
    {
        return dinuc_stat[0][8][11];
    }
    double GetDinuc12_GATA()
    {
        return dinuc_stat[0][8][12];
    }
    double GetDinuc12_GATC()
    {
        return dinuc_stat[0][8][13];
    }
    double GetDinuc12_GATG()
    {
        return dinuc_stat[0][8][14];
    }
    double GetDinuc12_GATT()
    {
        return dinuc_stat[0][8][15];
    }
    double GetDinuc12_GCAA()
    {
        return dinuc_stat[0][9][0];
    }
    double GetDinuc12_GCAC()
    {
        return dinuc_stat[0][9][1];
    }
    double GetDinuc12_GCAG()
    {
        return dinuc_stat[0][9][2];
    }
    double GetDinuc12_GCAT()
    {
        return dinuc_stat[0][9][3];
    }
    double GetDinuc12_GCCA()
    {
        return dinuc_stat[0][9][4];
    }
    double GetDinuc12_GCCC()
    {
        return dinuc_stat[0][9][5];
    }
    double GetDinuc12_GCCG()
    {
        return dinuc_stat[0][9][6];
    }
    double GetDinuc12_GCCT()
    {
        return dinuc_stat[0][9][7];
    }
    double GetDinuc12_GCGA()
    {
        return dinuc_stat[0][9][8];
    }
    double GetDinuc12_GCGC()
    {
        return dinuc_stat[0][9][9];
    }
    double GetDinuc12_GCGG()
    {
        return dinuc_stat[0][9][10];
    }
    double GetDinuc12_GCGT()
    {
        return dinuc_stat[0][9][11];
    }
    double GetDinuc12_GCTA()
    {
        return dinuc_stat[0][9][12];
    }
    double GetDinuc12_GCTC()
    {
        return dinuc_stat[0][9][13];
    }
    double GetDinuc12_GCTG()
    {
        return dinuc_stat[0][9][14];
    }
    double GetDinuc12_GCTT()
    {
        return dinuc_stat[0][9][15];
    }
    double GetDinuc12_GGAA()
    {
        return dinuc_stat[0][10][0];
    }
    double GetDinuc12_GGAC()
    {
        return dinuc_stat[0][10][1];
    }
    double GetDinuc12_GGAG()
    {
        return dinuc_stat[0][10][2];
    }
    double GetDinuc12_GGAT()
    {
        return dinuc_stat[0][10][3];
    }
    double GetDinuc12_GGCA()
    {
        return dinuc_stat[0][10][4];
    }
    double GetDinuc12_GGCC()
    {
        return dinuc_stat[0][10][5];
    }
    double GetDinuc12_GGCG()
    {
        return dinuc_stat[0][10][6];
    }
    double GetDinuc12_GGCT()
    {
        return dinuc_stat[0][10][7];
    }
    double GetDinuc12_GGGA()
    {
        return dinuc_stat[0][10][8];
    }
    double GetDinuc12_GGGC()
    {
        return dinuc_stat[0][10][9];
    }
    double GetDinuc12_GGGG()
    {
        return dinuc_stat[0][10][10];
    }
    double GetDinuc12_GGGT()
    {
        return dinuc_stat[0][10][11];
    }
    double GetDinuc12_GGTA()
    {
        return dinuc_stat[0][10][12];
    }
    double GetDinuc12_GGTC()
    {
        return dinuc_stat[0][10][13];
    }
    double GetDinuc12_GGTG()
    {
        return dinuc_stat[0][10][14];
    }
    double GetDinuc12_GGTT()
    {
        return dinuc_stat[0][10][15];
    }
    double GetDinuc12_GTAA()
    {
        return dinuc_stat[0][11][0];
    }
    double GetDinuc12_GTAC()
    {
        return dinuc_stat[0][11][1];
    }
    double GetDinuc12_GTAG()
    {
        return dinuc_stat[0][11][2];
    }
    double GetDinuc12_GTAT()
    {
        return dinuc_stat[0][11][3];
    }
    double GetDinuc12_GTCA()
    {
        return dinuc_stat[0][11][4];
    }
    double GetDinuc12_GTCC()
    {
        return dinuc_stat[0][11][5];
    }
    double GetDinuc12_GTCG()
    {
        return dinuc_stat[0][11][6];
    }
    double GetDinuc12_GTCT()
    {
        return dinuc_stat[0][11][7];
    }
    double GetDinuc12_GTGA()
    {
        return dinuc_stat[0][11][8];
    }
    double GetDinuc12_GTGC()
    {
        return dinuc_stat[0][11][9];
    }
    double GetDinuc12_GTGG()
    {
        return dinuc_stat[0][11][10];
    }
    double GetDinuc12_GTGT()
    {
        return dinuc_stat[0][11][11];
    }
    double GetDinuc12_GTTA()
    {
        return dinuc_stat[0][11][12];
    }
    double GetDinuc12_GTTC()
    {
        return dinuc_stat[0][11][13];
    }
    double GetDinuc12_GTTG()
    {
        return dinuc_stat[0][11][14];
    }
    double GetDinuc12_GTTT()
    {
        return dinuc_stat[0][11][15];
    }
    double GetDinuc12_TAAA()
    {
        return dinuc_stat[0][12][0];
    }
    double GetDinuc12_TAAC()
    {
        return dinuc_stat[0][12][1];
    }
    double GetDinuc12_TAAG()
    {
        return dinuc_stat[0][12][2];
    }
    double GetDinuc12_TAAT()
    {
        return dinuc_stat[0][12][3];
    }
    double GetDinuc12_TACA()
    {
        return dinuc_stat[0][12][4];
    }
    double GetDinuc12_TACC()
    {
        return dinuc_stat[0][12][5];
    }
    double GetDinuc12_TACG()
    {
        return dinuc_stat[0][12][6];
    }
    double GetDinuc12_TACT()
    {
        return dinuc_stat[0][12][7];
    }
    double GetDinuc12_TAGA()
    {
        return dinuc_stat[0][12][8];
    }
    double GetDinuc12_TAGC()
    {
        return dinuc_stat[0][12][9];
    }
    double GetDinuc12_TAGG()
    {
        return dinuc_stat[0][12][10];
    }
    double GetDinuc12_TAGT()
    {
        return dinuc_stat[0][12][11];
    }
    double GetDinuc12_TATA()
    {
        return dinuc_stat[0][12][12];
    }
    double GetDinuc12_TATC()
    {
        return dinuc_stat[0][12][13];
    }
    double GetDinuc12_TATG()
    {
        return dinuc_stat[0][12][14];
    }
    double GetDinuc12_TATT()
    {
        return dinuc_stat[0][12][15];
    }
    double GetDinuc12_TCAA()
    {
        return dinuc_stat[0][13][0];
    }
    double GetDinuc12_TCAC()
    {
        return dinuc_stat[0][13][1];
    }
    double GetDinuc12_TCAG()
    {
        return dinuc_stat[0][13][2];
    }
    double GetDinuc12_TCAT()
    {
        return dinuc_stat[0][13][3];
    }
    double GetDinuc12_TCCA()
    {
        return dinuc_stat[0][13][4];
    }
    double GetDinuc12_TCCC()
    {
        return dinuc_stat[0][13][5];
    }
    double GetDinuc12_TCCG()
    {
        return dinuc_stat[0][13][6];
    }
    double GetDinuc12_TCCT()
    {
        return dinuc_stat[0][13][7];
    }
    double GetDinuc12_TCGA()
    {
        return dinuc_stat[0][13][8];
    }
    double GetDinuc12_TCGC()
    {
        return dinuc_stat[0][13][9];
    }
    double GetDinuc12_TCGG()
    {
        return dinuc_stat[0][13][10];
    }
    double GetDinuc12_TCGT()
    {
        return dinuc_stat[0][13][11];
    }
    double GetDinuc12_TCTA()
    {
        return dinuc_stat[0][13][12];
    }
    double GetDinuc12_TCTC()
    {
        return dinuc_stat[0][13][13];
    }
    double GetDinuc12_TCTG()
    {
        return dinuc_stat[0][13][14];
    }
    double GetDinuc12_TCTT()
    {
        return dinuc_stat[0][13][15];
    }
    double GetDinuc12_TGAA()
    {
        return dinuc_stat[0][14][0];
    }
    double GetDinuc12_TGAC()
    {
        return dinuc_stat[0][14][1];
    }
    double GetDinuc12_TGAG()
    {
        return dinuc_stat[0][14][2];
    }
    double GetDinuc12_TGAT()
    {
        return dinuc_stat[0][14][3];
    }
    double GetDinuc12_TGCA()
    {
        return dinuc_stat[0][14][4];
    }
    double GetDinuc12_TGCC()
    {
        return dinuc_stat[0][14][5];
    }
    double GetDinuc12_TGCG()
    {
        return dinuc_stat[0][14][6];
    }
    double GetDinuc12_TGCT()
    {
        return dinuc_stat[0][14][7];
    }
    double GetDinuc12_TGGA()
    {
        return dinuc_stat[0][14][8];
    }
    double GetDinuc12_TGGC()
    {
        return dinuc_stat[0][14][9];
    }
    double GetDinuc12_TGGG()
    {
        return dinuc_stat[0][14][10];
    }
    double GetDinuc12_TGGT()
    {
        return dinuc_stat[0][14][11];
    }
    double GetDinuc12_TGTA()
    {
        return dinuc_stat[0][14][12];
    }
    double GetDinuc12_TGTC()
    {
        return dinuc_stat[0][14][13];
    }
    double GetDinuc12_TGTG()
    {
        return dinuc_stat[0][14][14];
    }
    double GetDinuc12_TGTT()
    {
        return dinuc_stat[0][14][15];
    }
    double GetDinuc12_TTAA()
    {
        return dinuc_stat[0][15][0];
    }
    double GetDinuc12_TTAC()
    {
        return dinuc_stat[0][15][1];
    }
    double GetDinuc12_TTAG()
    {
        return dinuc_stat[0][15][2];
    }
    double GetDinuc12_TTAT()
    {
        return dinuc_stat[0][15][3];
    }
    double GetDinuc12_TTCA()
    {
        return dinuc_stat[0][15][4];
    }
    double GetDinuc12_TTCC()
    {
        return dinuc_stat[0][15][5];
    }
    double GetDinuc12_TTCG()
    {
        return dinuc_stat[0][15][6];
    }
    double GetDinuc12_TTCT()
    {
        return dinuc_stat[0][15][7];
    }
    double GetDinuc12_TTGA()
    {
        return dinuc_stat[0][15][8];
    }
    double GetDinuc12_TTGC()
    {
        return dinuc_stat[0][15][9];
    }
    double GetDinuc12_TTGG()
    {
        return dinuc_stat[0][15][10];
    }
    double GetDinuc12_TTGT()
    {
        return dinuc_stat[0][15][11];
    }
    double GetDinuc12_TTTA()
    {
        return dinuc_stat[0][15][12];
    }
    double GetDinuc12_TTTC()
    {
        return dinuc_stat[0][15][13];
    }
    double GetDinuc12_TTTG()
    {
        return dinuc_stat[0][15][14];
    }
    double GetDinuc12_TTTT()
    {
        return dinuc_stat[0][15][15];
    }

    double GetDinuc23_AAAA()
    {
        return dinuc_stat[1][0][0];
    }
    double GetDinuc23_AAAC()
    {
        return dinuc_stat[1][0][1];
    }
    double GetDinuc23_AAAG()
    {
        return dinuc_stat[1][0][2];
    }
    double GetDinuc23_AAAT()
    {
        return dinuc_stat[1][0][3];
    }
    double GetDinuc23_AACA()
    {
        return dinuc_stat[1][0][4];
    }
    double GetDinuc23_AACC()
    {
        return dinuc_stat[1][0][5];
    }
    double GetDinuc23_AACG()
    {
        return dinuc_stat[1][0][6];
    }
    double GetDinuc23_AACT()
    {
        return dinuc_stat[1][0][7];
    }
    double GetDinuc23_AAGA()
    {
        return dinuc_stat[1][0][8];
    }
    double GetDinuc23_AAGC()
    {
        return dinuc_stat[1][0][9];
    }
    double GetDinuc23_AAGG()
    {
        return dinuc_stat[1][0][10];
    }
    double GetDinuc23_AAGT()
    {
        return dinuc_stat[1][0][11];
    }
    double GetDinuc23_AATA()
    {
        return dinuc_stat[1][0][12];
    }
    double GetDinuc23_AATC()
    {
        return dinuc_stat[1][0][13];
    }
    double GetDinuc23_AATG()
    {
        return dinuc_stat[1][0][14];
    }
    double GetDinuc23_AATT()
    {
        return dinuc_stat[1][0][15];
    }
    double GetDinuc23_ACAA()
    {
        return dinuc_stat[1][1][0];
    }
    double GetDinuc23_ACAC()
    {
        return dinuc_stat[1][1][1];
    }
    double GetDinuc23_ACAG()
    {
        return dinuc_stat[1][1][2];
    }
    double GetDinuc23_ACAT()
    {
        return dinuc_stat[1][1][3];
    }
    double GetDinuc23_ACCA()
    {
        return dinuc_stat[1][1][4];
    }
    double GetDinuc23_ACCC()
    {
        return dinuc_stat[1][1][5];
    }
    double GetDinuc23_ACCG()
    {
        return dinuc_stat[1][1][6];
    }
    double GetDinuc23_ACCT()
    {
        return dinuc_stat[1][1][7];
    }
    double GetDinuc23_ACGA()
    {
        return dinuc_stat[1][1][8];
    }
    double GetDinuc23_ACGC()
    {
        return dinuc_stat[1][1][9];
    }
    double GetDinuc23_ACGG()
    {
        return dinuc_stat[1][1][10];
    }
    double GetDinuc23_ACGT()
    {
        return dinuc_stat[1][1][11];
    }
    double GetDinuc23_ACTA()
    {
        return dinuc_stat[1][1][12];
    }
    double GetDinuc23_ACTC()
    {
        return dinuc_stat[1][1][13];
    }
    double GetDinuc23_ACTG()
    {
        return dinuc_stat[1][1][14];
    }
    double GetDinuc23_ACTT()
    {
        return dinuc_stat[1][1][15];
    }
    double GetDinuc23_AGAA()
    {
        return dinuc_stat[1][2][0];
    }
    double GetDinuc23_AGAC()
    {
        return dinuc_stat[1][2][1];
    }
    double GetDinuc23_AGAG()
    {
        return dinuc_stat[1][2][2];
    }
    double GetDinuc23_AGAT()
    {
        return dinuc_stat[1][2][3];
    }
    double GetDinuc23_AGCA()
    {
        return dinuc_stat[1][2][4];
    }
    double GetDinuc23_AGCC()
    {
        return dinuc_stat[1][2][5];
    }
    double GetDinuc23_AGCG()
    {
        return dinuc_stat[1][2][6];
    }
    double GetDinuc23_AGCT()
    {
        return dinuc_stat[1][2][7];
    }
    double GetDinuc23_AGGA()
    {
        return dinuc_stat[1][2][8];
    }
    double GetDinuc23_AGGC()
    {
        return dinuc_stat[1][2][9];
    }
    double GetDinuc23_AGGG()
    {
        return dinuc_stat[1][2][10];
    }
    double GetDinuc23_AGGT()
    {
        return dinuc_stat[1][2][11];
    }
    double GetDinuc23_AGTA()
    {
        return dinuc_stat[1][2][12];
    }
    double GetDinuc23_AGTC()
    {
        return dinuc_stat[1][2][13];
    }
    double GetDinuc23_AGTG()
    {
        return dinuc_stat[1][2][14];
    }
    double GetDinuc23_AGTT()
    {
        return dinuc_stat[1][2][15];
    }
    double GetDinuc23_ATAA()
    {
        return dinuc_stat[1][3][0];
    }
    double GetDinuc23_ATAC()
    {
        return dinuc_stat[1][3][1];
    }
    double GetDinuc23_ATAG()
    {
        return dinuc_stat[1][3][2];
    }
    double GetDinuc23_ATAT()
    {
        return dinuc_stat[1][3][3];
    }
    double GetDinuc23_ATCA()
    {
        return dinuc_stat[1][3][4];
    }
    double GetDinuc23_ATCC()
    {
        return dinuc_stat[1][3][5];
    }
    double GetDinuc23_ATCG()
    {
        return dinuc_stat[1][3][6];
    }
    double GetDinuc23_ATCT()
    {
        return dinuc_stat[1][3][7];
    }
    double GetDinuc23_ATGA()
    {
        return dinuc_stat[1][3][8];
    }
    double GetDinuc23_ATGC()
    {
        return dinuc_stat[1][3][9];
    }
    double GetDinuc23_ATGG()
    {
        return dinuc_stat[1][3][10];
    }
    double GetDinuc23_ATGT()
    {
        return dinuc_stat[1][3][11];
    }
    double GetDinuc23_ATTA()
    {
        return dinuc_stat[1][3][12];
    }
    double GetDinuc23_ATTC()
    {
        return dinuc_stat[1][3][13];
    }
    double GetDinuc23_ATTG()
    {
        return dinuc_stat[1][3][14];
    }
    double GetDinuc23_ATTT()
    {
        return dinuc_stat[1][3][15];
    }
    double GetDinuc23_CAAA()
    {
        return dinuc_stat[1][4][0];
    }
    double GetDinuc23_CAAC()
    {
        return dinuc_stat[1][4][1];
    }
    double GetDinuc23_CAAG()
    {
        return dinuc_stat[1][4][2];
    }
    double GetDinuc23_CAAT()
    {
        return dinuc_stat[1][4][3];
    }
    double GetDinuc23_CACA()
    {
        return dinuc_stat[1][4][4];
    }
    double GetDinuc23_CACC()
    {
        return dinuc_stat[1][4][5];
    }
    double GetDinuc23_CACG()
    {
        return dinuc_stat[1][4][6];
    }
    double GetDinuc23_CACT()
    {
        return dinuc_stat[1][4][7];
    }
    double GetDinuc23_CAGA()
    {
        return dinuc_stat[1][4][8];
    }
    double GetDinuc23_CAGC()
    {
        return dinuc_stat[1][4][9];
    }
    double GetDinuc23_CAGG()
    {
        return dinuc_stat[1][4][10];
    }
    double GetDinuc23_CAGT()
    {
        return dinuc_stat[1][4][11];
    }
    double GetDinuc23_CATA()
    {
        return dinuc_stat[1][4][12];
    }
    double GetDinuc23_CATC()
    {
        return dinuc_stat[1][4][13];
    }
    double GetDinuc23_CATG()
    {
        return dinuc_stat[1][4][14];
    }
    double GetDinuc23_CATT()
    {
        return dinuc_stat[1][4][15];
    }
    double GetDinuc23_CCAA()
    {
        return dinuc_stat[1][5][0];
    }
    double GetDinuc23_CCAC()
    {
        return dinuc_stat[1][5][1];
    }
    double GetDinuc23_CCAG()
    {
        return dinuc_stat[1][5][2];
    }
    double GetDinuc23_CCAT()
    {
        return dinuc_stat[1][5][3];
    }
    double GetDinuc23_CCCA()
    {
        return dinuc_stat[1][5][4];
    }
    double GetDinuc23_CCCC()
    {
        return dinuc_stat[1][5][5];
    }
    double GetDinuc23_CCCG()
    {
        return dinuc_stat[1][5][6];
    }
    double GetDinuc23_CCCT()
    {
        return dinuc_stat[1][5][7];
    }
    double GetDinuc23_CCGA()
    {
        return dinuc_stat[1][5][8];
    }
    double GetDinuc23_CCGC()
    {
        return dinuc_stat[1][5][9];
    }
    double GetDinuc23_CCGG()
    {
        return dinuc_stat[1][5][10];
    }
    double GetDinuc23_CCGT()
    {
        return dinuc_stat[1][5][11];
    }
    double GetDinuc23_CCTA()
    {
        return dinuc_stat[1][5][12];
    }
    double GetDinuc23_CCTC()
    {
        return dinuc_stat[1][5][13];
    }
    double GetDinuc23_CCTG()
    {
        return dinuc_stat[1][5][14];
    }
    double GetDinuc23_CCTT()
    {
        return dinuc_stat[1][5][15];
    }
    double GetDinuc23_CGAA()
    {
        return dinuc_stat[1][6][0];
    }
    double GetDinuc23_CGAC()
    {
        return dinuc_stat[1][6][1];
    }
    double GetDinuc23_CGAG()
    {
        return dinuc_stat[1][6][2];
    }
    double GetDinuc23_CGAT()
    {
        return dinuc_stat[1][6][3];
    }
    double GetDinuc23_CGCA()
    {
        return dinuc_stat[1][6][4];
    }
    double GetDinuc23_CGCC()
    {
        return dinuc_stat[1][6][5];
    }
    double GetDinuc23_CGCG()
    {
        return dinuc_stat[1][6][6];
    }
    double GetDinuc23_CGCT()
    {
        return dinuc_stat[1][6][7];
    }
    double GetDinuc23_CGGA()
    {
        return dinuc_stat[1][6][8];
    }
    double GetDinuc23_CGGC()
    {
        return dinuc_stat[1][6][9];
    }
    double GetDinuc23_CGGG()
    {
        return dinuc_stat[1][6][10];
    }
    double GetDinuc23_CGGT()
    {
        return dinuc_stat[1][6][11];
    }
    double GetDinuc23_CGTA()
    {
        return dinuc_stat[1][6][12];
    }
    double GetDinuc23_CGTC()
    {
        return dinuc_stat[1][6][13];
    }
    double GetDinuc23_CGTG()
    {
        return dinuc_stat[1][6][14];
    }
    double GetDinuc23_CGTT()
    {
        return dinuc_stat[1][6][15];
    }
    double GetDinuc23_CTAA()
    {
        return dinuc_stat[1][7][0];
    }
    double GetDinuc23_CTAC()
    {
        return dinuc_stat[1][7][1];
    }
    double GetDinuc23_CTAG()
    {
        return dinuc_stat[1][7][2];
    }
    double GetDinuc23_CTAT()
    {
        return dinuc_stat[1][7][3];
    }
    double GetDinuc23_CTCA()
    {
        return dinuc_stat[1][7][4];
    }
    double GetDinuc23_CTCC()
    {
        return dinuc_stat[1][7][5];
    }
    double GetDinuc23_CTCG()
    {
        return dinuc_stat[1][7][6];
    }
    double GetDinuc23_CTCT()
    {
        return dinuc_stat[1][7][7];
    }
    double GetDinuc23_CTGA()
    {
        return dinuc_stat[1][7][8];
    }
    double GetDinuc23_CTGC()
    {
        return dinuc_stat[1][7][9];
    }
    double GetDinuc23_CTGG()
    {
        return dinuc_stat[1][7][10];
    }
    double GetDinuc23_CTGT()
    {
        return dinuc_stat[1][7][11];
    }
    double GetDinuc23_CTTA()
    {
        return dinuc_stat[1][7][12];
    }
    double GetDinuc23_CTTC()
    {
        return dinuc_stat[1][7][13];
    }
    double GetDinuc23_CTTG()
    {
        return dinuc_stat[1][7][14];
    }
    double GetDinuc23_CTTT()
    {
        return dinuc_stat[1][7][15];
    }
    double GetDinuc23_GAAA()
    {
        return dinuc_stat[1][8][0];
    }
    double GetDinuc23_GAAC()
    {
        return dinuc_stat[1][8][1];
    }
    double GetDinuc23_GAAG()
    {
        return dinuc_stat[1][8][2];
    }
    double GetDinuc23_GAAT()
    {
        return dinuc_stat[1][8][3];
    }
    double GetDinuc23_GACA()
    {
        return dinuc_stat[1][8][4];
    }
    double GetDinuc23_GACC()
    {
        return dinuc_stat[1][8][5];
    }
    double GetDinuc23_GACG()
    {
        return dinuc_stat[1][8][6];
    }
    double GetDinuc23_GACT()
    {
        return dinuc_stat[1][8][7];
    }
    double GetDinuc23_GAGA()
    {
        return dinuc_stat[1][8][8];
    }
    double GetDinuc23_GAGC()
    {
        return dinuc_stat[1][8][9];
    }
    double GetDinuc23_GAGG()
    {
        return dinuc_stat[1][8][10];
    }
    double GetDinuc23_GAGT()
    {
        return dinuc_stat[1][8][11];
    }
    double GetDinuc23_GATA()
    {
        return dinuc_stat[1][8][12];
    }
    double GetDinuc23_GATC()
    {
        return dinuc_stat[1][8][13];
    }
    double GetDinuc23_GATG()
    {
        return dinuc_stat[1][8][14];
    }
    double GetDinuc23_GATT()
    {
        return dinuc_stat[1][8][15];
    }
    double GetDinuc23_GCAA()
    {
        return dinuc_stat[1][9][0];
    }
    double GetDinuc23_GCAC()
    {
        return dinuc_stat[1][9][1];
    }
    double GetDinuc23_GCAG()
    {
        return dinuc_stat[1][9][2];
    }
    double GetDinuc23_GCAT()
    {
        return dinuc_stat[1][9][3];
    }
    double GetDinuc23_GCCA()
    {
        return dinuc_stat[1][9][4];
    }
    double GetDinuc23_GCCC()
    {
        return dinuc_stat[1][9][5];
    }
    double GetDinuc23_GCCG()
    {
        return dinuc_stat[1][9][6];
    }
    double GetDinuc23_GCCT()
    {
        return dinuc_stat[1][9][7];
    }
    double GetDinuc23_GCGA()
    {
        return dinuc_stat[1][9][8];
    }
    double GetDinuc23_GCGC()
    {
        return dinuc_stat[1][9][9];
    }
    double GetDinuc23_GCGG()
    {
        return dinuc_stat[1][9][10];
    }
    double GetDinuc23_GCGT()
    {
        return dinuc_stat[1][9][11];
    }
    double GetDinuc23_GCTA()
    {
        return dinuc_stat[1][9][12];
    }
    double GetDinuc23_GCTC()
    {
        return dinuc_stat[1][9][13];
    }
    double GetDinuc23_GCTG()
    {
        return dinuc_stat[1][9][14];
    }
    double GetDinuc23_GCTT()
    {
        return dinuc_stat[1][9][15];
    }
    double GetDinuc23_GGAA()
    {
        return dinuc_stat[1][10][0];
    }
    double GetDinuc23_GGAC()
    {
        return dinuc_stat[1][10][1];
    }
    double GetDinuc23_GGAG()
    {
        return dinuc_stat[1][10][2];
    }
    double GetDinuc23_GGAT()
    {
        return dinuc_stat[1][10][3];
    }
    double GetDinuc23_GGCA()
    {
        return dinuc_stat[1][10][4];
    }
    double GetDinuc23_GGCC()
    {
        return dinuc_stat[1][10][5];
    }
    double GetDinuc23_GGCG()
    {
        return dinuc_stat[1][10][6];
    }
    double GetDinuc23_GGCT()
    {
        return dinuc_stat[1][10][7];
    }
    double GetDinuc23_GGGA()
    {
        return dinuc_stat[1][10][8];
    }
    double GetDinuc23_GGGC()
    {
        return dinuc_stat[1][10][9];
    }
    double GetDinuc23_GGGG()
    {
        return dinuc_stat[1][10][10];
    }
    double GetDinuc23_GGGT()
    {
        return dinuc_stat[1][10][11];
    }
    double GetDinuc23_GGTA()
    {
        return dinuc_stat[1][10][12];
    }
    double GetDinuc23_GGTC()
    {
        return dinuc_stat[1][10][13];
    }
    double GetDinuc23_GGTG()
    {
        return dinuc_stat[1][10][14];
    }
    double GetDinuc23_GGTT()
    {
        return dinuc_stat[1][10][15];
    }
    double GetDinuc23_GTAA()
    {
        return dinuc_stat[1][11][0];
    }
    double GetDinuc23_GTAC()
    {
        return dinuc_stat[1][11][1];
    }
    double GetDinuc23_GTAG()
    {
        return dinuc_stat[1][11][2];
    }
    double GetDinuc23_GTAT()
    {
        return dinuc_stat[1][11][3];
    }
    double GetDinuc23_GTCA()
    {
        return dinuc_stat[1][11][4];
    }
    double GetDinuc23_GTCC()
    {
        return dinuc_stat[1][11][5];
    }
    double GetDinuc23_GTCG()
    {
        return dinuc_stat[1][11][6];
    }
    double GetDinuc23_GTCT()
    {
        return dinuc_stat[1][11][7];
    }
    double GetDinuc23_GTGA()
    {
        return dinuc_stat[1][11][8];
    }
    double GetDinuc23_GTGC()
    {
        return dinuc_stat[1][11][9];
    }
    double GetDinuc23_GTGG()
    {
        return dinuc_stat[1][11][10];
    }
    double GetDinuc23_GTGT()
    {
        return dinuc_stat[1][11][11];
    }
    double GetDinuc23_GTTA()
    {
        return dinuc_stat[1][11][12];
    }
    double GetDinuc23_GTTC()
    {
        return dinuc_stat[1][11][13];
    }
    double GetDinuc23_GTTG()
    {
        return dinuc_stat[1][11][14];
    }
    double GetDinuc23_GTTT()
    {
        return dinuc_stat[1][11][15];
    }
    double GetDinuc23_TAAA()
    {
        return dinuc_stat[1][12][0];
    }
    double GetDinuc23_TAAC()
    {
        return dinuc_stat[1][12][1];
    }
    double GetDinuc23_TAAG()
    {
        return dinuc_stat[1][12][2];
    }
    double GetDinuc23_TAAT()
    {
        return dinuc_stat[1][12][3];
    }
    double GetDinuc23_TACA()
    {
        return dinuc_stat[1][12][4];
    }
    double GetDinuc23_TACC()
    {
        return dinuc_stat[1][12][5];
    }
    double GetDinuc23_TACG()
    {
        return dinuc_stat[1][12][6];
    }
    double GetDinuc23_TACT()
    {
        return dinuc_stat[1][12][7];
    }
    double GetDinuc23_TAGA()
    {
        return dinuc_stat[1][12][8];
    }
    double GetDinuc23_TAGC()
    {
        return dinuc_stat[1][12][9];
    }
    double GetDinuc23_TAGG()
    {
        return dinuc_stat[1][12][10];
    }
    double GetDinuc23_TAGT()
    {
        return dinuc_stat[1][12][11];
    }
    double GetDinuc23_TATA()
    {
        return dinuc_stat[1][12][12];
    }
    double GetDinuc23_TATC()
    {
        return dinuc_stat[1][12][13];
    }
    double GetDinuc23_TATG()
    {
        return dinuc_stat[1][12][14];
    }
    double GetDinuc23_TATT()
    {
        return dinuc_stat[1][12][15];
    }
    double GetDinuc23_TCAA()
    {
        return dinuc_stat[1][13][0];
    }
    double GetDinuc23_TCAC()
    {
        return dinuc_stat[1][13][1];
    }
    double GetDinuc23_TCAG()
    {
        return dinuc_stat[1][13][2];
    }
    double GetDinuc23_TCAT()
    {
        return dinuc_stat[1][13][3];
    }
    double GetDinuc23_TCCA()
    {
        return dinuc_stat[1][13][4];
    }
    double GetDinuc23_TCCC()
    {
        return dinuc_stat[1][13][5];
    }
    double GetDinuc23_TCCG()
    {
        return dinuc_stat[1][13][6];
    }
    double GetDinuc23_TCCT()
    {
        return dinuc_stat[1][13][7];
    }
    double GetDinuc23_TCGA()
    {
        return dinuc_stat[1][13][8];
    }
    double GetDinuc23_TCGC()
    {
        return dinuc_stat[1][13][9];
    }
    double GetDinuc23_TCGG()
    {
        return dinuc_stat[1][13][10];
    }
    double GetDinuc23_TCGT()
    {
        return dinuc_stat[1][13][11];
    }
    double GetDinuc23_TCTA()
    {
        return dinuc_stat[1][13][12];
    }
    double GetDinuc23_TCTC()
    {
        return dinuc_stat[1][13][13];
    }
    double GetDinuc23_TCTG()
    {
        return dinuc_stat[1][13][14];
    }
    double GetDinuc23_TCTT()
    {
        return dinuc_stat[1][13][15];
    }
    double GetDinuc23_TGAA()
    {
        return dinuc_stat[1][14][0];
    }
    double GetDinuc23_TGAC()
    {
        return dinuc_stat[1][14][1];
    }
    double GetDinuc23_TGAG()
    {
        return dinuc_stat[1][14][2];
    }
    double GetDinuc23_TGAT()
    {
        return dinuc_stat[1][14][3];
    }
    double GetDinuc23_TGCA()
    {
        return dinuc_stat[1][14][4];
    }
    double GetDinuc23_TGCC()
    {
        return dinuc_stat[1][14][5];
    }
    double GetDinuc23_TGCG()
    {
        return dinuc_stat[1][14][6];
    }
    double GetDinuc23_TGCT()
    {
        return dinuc_stat[1][14][7];
    }
    double GetDinuc23_TGGA()
    {
        return dinuc_stat[1][14][8];
    }
    double GetDinuc23_TGGC()
    {
        return dinuc_stat[1][14][9];
    }
    double GetDinuc23_TGGG()
    {
        return dinuc_stat[1][14][10];
    }
    double GetDinuc23_TGGT()
    {
        return dinuc_stat[1][14][11];
    }
    double GetDinuc23_TGTA()
    {
        return dinuc_stat[1][14][12];
    }
    double GetDinuc23_TGTC()
    {
        return dinuc_stat[1][14][13];
    }
    double GetDinuc23_TGTG()
    {
        return dinuc_stat[1][14][14];
    }
    double GetDinuc23_TGTT()
    {
        return dinuc_stat[1][14][15];
    }
    double GetDinuc23_TTAA()
    {
        return dinuc_stat[1][15][0];
    }
    double GetDinuc23_TTAC()
    {
        return dinuc_stat[1][15][1];
    }
    double GetDinuc23_TTAG()
    {
        return dinuc_stat[1][15][2];
    }
    double GetDinuc23_TTAT()
    {
        return dinuc_stat[1][15][3];
    }
    double GetDinuc23_TTCA()
    {
        return dinuc_stat[1][15][4];
    }
    double GetDinuc23_TTCC()
    {
        return dinuc_stat[1][15][5];
    }
    double GetDinuc23_TTCG()
    {
        return dinuc_stat[1][15][6];
    }
    double GetDinuc23_TTCT()
    {
        return dinuc_stat[1][15][7];
    }
    double GetDinuc23_TTGA()
    {
        return dinuc_stat[1][15][8];
    }
    double GetDinuc23_TTGC()
    {
        return dinuc_stat[1][15][9];
    }
    double GetDinuc23_TTGG()
    {
        return dinuc_stat[1][15][10];
    }
    double GetDinuc23_TTGT()
    {
        return dinuc_stat[1][15][11];
    }
    double GetDinuc23_TTTA()
    {
        return dinuc_stat[1][15][12];
    }
    double GetDinuc23_TTTC()
    {
        return dinuc_stat[1][15][13];
    }
    double GetDinuc23_TTTG()
    {
        return dinuc_stat[1][15][14];
    }
    double GetDinuc23_TTTT()
    {
        return dinuc_stat[1][15][15];
    }

    double GetDinuc31_AAAA()
    {
        return dinuc_stat[2][0][0];
    }
    double GetDinuc31_AAAC()
    {
        return dinuc_stat[2][0][1];
    }
    double GetDinuc31_AAAG()
    {
        return dinuc_stat[2][0][2];
    }
    double GetDinuc31_AAAT()
    {
        return dinuc_stat[2][0][3];
    }
    double GetDinuc31_AACA()
    {
        return dinuc_stat[2][0][4];
    }
    double GetDinuc31_AACC()
    {
        return dinuc_stat[2][0][5];
    }
    double GetDinuc31_AACG()
    {
        return dinuc_stat[2][0][6];
    }
    double GetDinuc31_AACT()
    {
        return dinuc_stat[2][0][7];
    }
    double GetDinuc31_AAGA()
    {
        return dinuc_stat[2][0][8];
    }
    double GetDinuc31_AAGC()
    {
        return dinuc_stat[2][0][9];
    }
    double GetDinuc31_AAGG()
    {
        return dinuc_stat[2][0][10];
    }
    double GetDinuc31_AAGT()
    {
        return dinuc_stat[2][0][11];
    }
    double GetDinuc31_AATA()
    {
        return dinuc_stat[2][0][12];
    }
    double GetDinuc31_AATC()
    {
        return dinuc_stat[2][0][13];
    }
    double GetDinuc31_AATG()
    {
        return dinuc_stat[2][0][14];
    }
    double GetDinuc31_AATT()
    {
        return dinuc_stat[2][0][15];
    }
    double GetDinuc31_ACAA()
    {
        return dinuc_stat[2][1][0];
    }
    double GetDinuc31_ACAC()
    {
        return dinuc_stat[2][1][1];
    }
    double GetDinuc31_ACAG()
    {
        return dinuc_stat[2][1][2];
    }
    double GetDinuc31_ACAT()
    {
        return dinuc_stat[2][1][3];
    }
    double GetDinuc31_ACCA()
    {
        return dinuc_stat[2][1][4];
    }
    double GetDinuc31_ACCC()
    {
        return dinuc_stat[2][1][5];
    }
    double GetDinuc31_ACCG()
    {
        return dinuc_stat[2][1][6];
    }
    double GetDinuc31_ACCT()
    {
        return dinuc_stat[2][1][7];
    }
    double GetDinuc31_ACGA()
    {
        return dinuc_stat[2][1][8];
    }
    double GetDinuc31_ACGC()
    {
        return dinuc_stat[2][1][9];
    }
    double GetDinuc31_ACGG()
    {
        return dinuc_stat[2][1][10];
    }
    double GetDinuc31_ACGT()
    {
        return dinuc_stat[2][1][11];
    }
    double GetDinuc31_ACTA()
    {
        return dinuc_stat[2][1][12];
    }
    double GetDinuc31_ACTC()
    {
        return dinuc_stat[2][1][13];
    }
    double GetDinuc31_ACTG()
    {
        return dinuc_stat[2][1][14];
    }
    double GetDinuc31_ACTT()
    {
        return dinuc_stat[2][1][15];
    }
    double GetDinuc31_AGAA()
    {
        return dinuc_stat[2][2][0];
    }
    double GetDinuc31_AGAC()
    {
        return dinuc_stat[2][2][1];
    }
    double GetDinuc31_AGAG()
    {
        return dinuc_stat[2][2][2];
    }
    double GetDinuc31_AGAT()
    {
        return dinuc_stat[2][2][3];
    }
    double GetDinuc31_AGCA()
    {
        return dinuc_stat[2][2][4];
    }
    double GetDinuc31_AGCC()
    {
        return dinuc_stat[2][2][5];
    }
    double GetDinuc31_AGCG()
    {
        return dinuc_stat[2][2][6];
    }
    double GetDinuc31_AGCT()
    {
        return dinuc_stat[2][2][7];
    }
    double GetDinuc31_AGGA()
    {
        return dinuc_stat[2][2][8];
    }
    double GetDinuc31_AGGC()
    {
        return dinuc_stat[2][2][9];
    }
    double GetDinuc31_AGGG()
    {
        return dinuc_stat[2][2][10];
    }
    double GetDinuc31_AGGT()
    {
        return dinuc_stat[2][2][11];
    }
    double GetDinuc31_AGTA()
    {
        return dinuc_stat[2][2][12];
    }
    double GetDinuc31_AGTC()
    {
        return dinuc_stat[2][2][13];
    }
    double GetDinuc31_AGTG()
    {
        return dinuc_stat[2][2][14];
    }
    double GetDinuc31_AGTT()
    {
        return dinuc_stat[2][2][15];
    }
    double GetDinuc31_ATAA()
    {
        return dinuc_stat[2][3][0];
    }
    double GetDinuc31_ATAC()
    {
        return dinuc_stat[2][3][1];
    }
    double GetDinuc31_ATAG()
    {
        return dinuc_stat[2][3][2];
    }
    double GetDinuc31_ATAT()
    {
        return dinuc_stat[2][3][3];
    }
    double GetDinuc31_ATCA()
    {
        return dinuc_stat[2][3][4];
    }
    double GetDinuc31_ATCC()
    {
        return dinuc_stat[2][3][5];
    }
    double GetDinuc31_ATCG()
    {
        return dinuc_stat[2][3][6];
    }
    double GetDinuc31_ATCT()
    {
        return dinuc_stat[2][3][7];
    }
    double GetDinuc31_ATGA()
    {
        return dinuc_stat[2][3][8];
    }
    double GetDinuc31_ATGC()
    {
        return dinuc_stat[2][3][9];
    }
    double GetDinuc31_ATGG()
    {
        return dinuc_stat[2][3][10];
    }
    double GetDinuc31_ATGT()
    {
        return dinuc_stat[2][3][11];
    }
    double GetDinuc31_ATTA()
    {
        return dinuc_stat[2][3][12];
    }
    double GetDinuc31_ATTC()
    {
        return dinuc_stat[2][3][13];
    }
    double GetDinuc31_ATTG()
    {
        return dinuc_stat[2][3][14];
    }
    double GetDinuc31_ATTT()
    {
        return dinuc_stat[2][3][15];
    }
    double GetDinuc31_CAAA()
    {
        return dinuc_stat[2][4][0];
    }
    double GetDinuc31_CAAC()
    {
        return dinuc_stat[2][4][1];
    }
    double GetDinuc31_CAAG()
    {
        return dinuc_stat[2][4][2];
    }
    double GetDinuc31_CAAT()
    {
        return dinuc_stat[2][4][3];
    }
    double GetDinuc31_CACA()
    {
        return dinuc_stat[2][4][4];
    }
    double GetDinuc31_CACC()
    {
        return dinuc_stat[2][4][5];
    }
    double GetDinuc31_CACG()
    {
        return dinuc_stat[2][4][6];
    }
    double GetDinuc31_CACT()
    {
        return dinuc_stat[2][4][7];
    }
    double GetDinuc31_CAGA()
    {
        return dinuc_stat[2][4][8];
    }
    double GetDinuc31_CAGC()
    {
        return dinuc_stat[2][4][9];
    }
    double GetDinuc31_CAGG()
    {
        return dinuc_stat[2][4][10];
    }
    double GetDinuc31_CAGT()
    {
        return dinuc_stat[2][4][11];
    }
    double GetDinuc31_CATA()
    {
        return dinuc_stat[2][4][12];
    }
    double GetDinuc31_CATC()
    {
        return dinuc_stat[2][4][13];
    }
    double GetDinuc31_CATG()
    {
        return dinuc_stat[2][4][14];
    }
    double GetDinuc31_CATT()
    {
        return dinuc_stat[2][4][15];
    }
    double GetDinuc31_CCAA()
    {
        return dinuc_stat[2][5][0];
    }
    double GetDinuc31_CCAC()
    {
        return dinuc_stat[2][5][1];
    }
    double GetDinuc31_CCAG()
    {
        return dinuc_stat[2][5][2];
    }
    double GetDinuc31_CCAT()
    {
        return dinuc_stat[2][5][3];
    }
    double GetDinuc31_CCCA()
    {
        return dinuc_stat[2][5][4];
    }
    double GetDinuc31_CCCC()
    {
        return dinuc_stat[2][5][5];
    }
    double GetDinuc31_CCCG()
    {
        return dinuc_stat[2][5][6];
    }
    double GetDinuc31_CCCT()
    {
        return dinuc_stat[2][5][7];
    }
    double GetDinuc31_CCGA()
    {
        return dinuc_stat[2][5][8];
    }
    double GetDinuc31_CCGC()
    {
        return dinuc_stat[2][5][9];
    }
    double GetDinuc31_CCGG()
    {
        return dinuc_stat[2][5][10];
    }
    double GetDinuc31_CCGT()
    {
        return dinuc_stat[2][5][11];
    }
    double GetDinuc31_CCTA()
    {
        return dinuc_stat[2][5][12];
    }
    double GetDinuc31_CCTC()
    {
        return dinuc_stat[2][5][13];
    }
    double GetDinuc31_CCTG()
    {
        return dinuc_stat[2][5][14];
    }
    double GetDinuc31_CCTT()
    {
        return dinuc_stat[2][5][15];
    }
    double GetDinuc31_CGAA()
    {
        return dinuc_stat[2][6][0];
    }
    double GetDinuc31_CGAC()
    {
        return dinuc_stat[2][6][1];
    }
    double GetDinuc31_CGAG()
    {
        return dinuc_stat[2][6][2];
    }
    double GetDinuc31_CGAT()
    {
        return dinuc_stat[2][6][3];
    }
    double GetDinuc31_CGCA()
    {
        return dinuc_stat[2][6][4];
    }
    double GetDinuc31_CGCC()
    {
        return dinuc_stat[2][6][5];
    }
    double GetDinuc31_CGCG()
    {
        return dinuc_stat[2][6][6];
    }
    double GetDinuc31_CGCT()
    {
        return dinuc_stat[2][6][7];
    }
    double GetDinuc31_CGGA()
    {
        return dinuc_stat[2][6][8];
    }
    double GetDinuc31_CGGC()
    {
        return dinuc_stat[2][6][9];
    }
    double GetDinuc31_CGGG()
    {
        return dinuc_stat[2][6][10];
    }
    double GetDinuc31_CGGT()
    {
        return dinuc_stat[2][6][11];
    }
    double GetDinuc31_CGTA()
    {
        return dinuc_stat[2][6][12];
    }
    double GetDinuc31_CGTC()
    {
        return dinuc_stat[2][6][13];
    }
    double GetDinuc31_CGTG()
    {
        return dinuc_stat[2][6][14];
    }
    double GetDinuc31_CGTT()
    {
        return dinuc_stat[2][6][15];
    }
    double GetDinuc31_CTAA()
    {
        return dinuc_stat[2][7][0];
    }
    double GetDinuc31_CTAC()
    {
        return dinuc_stat[2][7][1];
    }
    double GetDinuc31_CTAG()
    {
        return dinuc_stat[2][7][2];
    }
    double GetDinuc31_CTAT()
    {
        return dinuc_stat[2][7][3];
    }
    double GetDinuc31_CTCA()
    {
        return dinuc_stat[2][7][4];
    }
    double GetDinuc31_CTCC()
    {
        return dinuc_stat[2][7][5];
    }
    double GetDinuc31_CTCG()
    {
        return dinuc_stat[2][7][6];
    }
    double GetDinuc31_CTCT()
    {
        return dinuc_stat[2][7][7];
    }
    double GetDinuc31_CTGA()
    {
        return dinuc_stat[2][7][8];
    }
    double GetDinuc31_CTGC()
    {
        return dinuc_stat[2][7][9];
    }
    double GetDinuc31_CTGG()
    {
        return dinuc_stat[2][7][10];
    }
    double GetDinuc31_CTGT()
    {
        return dinuc_stat[2][7][11];
    }
    double GetDinuc31_CTTA()
    {
        return dinuc_stat[2][7][12];
    }
    double GetDinuc31_CTTC()
    {
        return dinuc_stat[2][7][13];
    }
    double GetDinuc31_CTTG()
    {
        return dinuc_stat[2][7][14];
    }
    double GetDinuc31_CTTT()
    {
        return dinuc_stat[2][7][15];
    }
    double GetDinuc31_GAAA()
    {
        return dinuc_stat[2][8][0];
    }
    double GetDinuc31_GAAC()
    {
        return dinuc_stat[2][8][1];
    }
    double GetDinuc31_GAAG()
    {
        return dinuc_stat[2][8][2];
    }
    double GetDinuc31_GAAT()
    {
        return dinuc_stat[2][8][3];
    }
    double GetDinuc31_GACA()
    {
        return dinuc_stat[2][8][4];
    }
    double GetDinuc31_GACC()
    {
        return dinuc_stat[2][8][5];
    }
    double GetDinuc31_GACG()
    {
        return dinuc_stat[2][8][6];
    }
    double GetDinuc31_GACT()
    {
        return dinuc_stat[2][8][7];
    }
    double GetDinuc31_GAGA()
    {
        return dinuc_stat[2][8][8];
    }
    double GetDinuc31_GAGC()
    {
        return dinuc_stat[2][8][9];
    }
    double GetDinuc31_GAGG()
    {
        return dinuc_stat[2][8][10];
    }
    double GetDinuc31_GAGT()
    {
        return dinuc_stat[2][8][11];
    }
    double GetDinuc31_GATA()
    {
        return dinuc_stat[2][8][12];
    }
    double GetDinuc31_GATC()
    {
        return dinuc_stat[2][8][13];
    }
    double GetDinuc31_GATG()
    {
        return dinuc_stat[2][8][14];
    }
    double GetDinuc31_GATT()
    {
        return dinuc_stat[2][8][15];
    }
    double GetDinuc31_GCAA()
    {
        return dinuc_stat[2][9][0];
    }
    double GetDinuc31_GCAC()
    {
        return dinuc_stat[2][9][1];
    }
    double GetDinuc31_GCAG()
    {
        return dinuc_stat[2][9][2];
    }
    double GetDinuc31_GCAT()
    {
        return dinuc_stat[2][9][3];
    }
    double GetDinuc31_GCCA()
    {
        return dinuc_stat[2][9][4];
    }
    double GetDinuc31_GCCC()
    {
        return dinuc_stat[2][9][5];
    }
    double GetDinuc31_GCCG()
    {
        return dinuc_stat[2][9][6];
    }
    double GetDinuc31_GCCT()
    {
        return dinuc_stat[2][9][7];
    }
    double GetDinuc31_GCGA()
    {
        return dinuc_stat[2][9][8];
    }
    double GetDinuc31_GCGC()
    {
        return dinuc_stat[2][9][9];
    }
    double GetDinuc31_GCGG()
    {
        return dinuc_stat[2][9][10];
    }
    double GetDinuc31_GCGT()
    {
        return dinuc_stat[2][9][11];
    }
    double GetDinuc31_GCTA()
    {
        return dinuc_stat[2][9][12];
    }
    double GetDinuc31_GCTC()
    {
        return dinuc_stat[2][9][13];
    }
    double GetDinuc31_GCTG()
    {
        return dinuc_stat[2][9][14];
    }
    double GetDinuc31_GCTT()
    {
        return dinuc_stat[2][9][15];
    }
    double GetDinuc31_GGAA()
    {
        return dinuc_stat[2][10][0];
    }
    double GetDinuc31_GGAC()
    {
        return dinuc_stat[2][10][1];
    }
    double GetDinuc31_GGAG()
    {
        return dinuc_stat[2][10][2];
    }
    double GetDinuc31_GGAT()
    {
        return dinuc_stat[2][10][3];
    }
    double GetDinuc31_GGCA()
    {
        return dinuc_stat[2][10][4];
    }
    double GetDinuc31_GGCC()
    {
        return dinuc_stat[2][10][5];
    }
    double GetDinuc31_GGCG()
    {
        return dinuc_stat[2][10][6];
    }
    double GetDinuc31_GGCT()
    {
        return dinuc_stat[2][10][7];
    }
    double GetDinuc31_GGGA()
    {
        return dinuc_stat[2][10][8];
    }
    double GetDinuc31_GGGC()
    {
        return dinuc_stat[2][10][9];
    }
    double GetDinuc31_GGGG()
    {
        return dinuc_stat[2][10][10];
    }
    double GetDinuc31_GGGT()
    {
        return dinuc_stat[2][10][11];
    }
    double GetDinuc31_GGTA()
    {
        return dinuc_stat[2][10][12];
    }
    double GetDinuc31_GGTC()
    {
        return dinuc_stat[2][10][13];
    }
    double GetDinuc31_GGTG()
    {
        return dinuc_stat[2][10][14];
    }
    double GetDinuc31_GGTT()
    {
        return dinuc_stat[2][10][15];
    }
    double GetDinuc31_GTAA()
    {
        return dinuc_stat[2][11][0];
    }
    double GetDinuc31_GTAC()
    {
        return dinuc_stat[2][11][1];
    }
    double GetDinuc31_GTAG()
    {
        return dinuc_stat[2][11][2];
    }
    double GetDinuc31_GTAT()
    {
        return dinuc_stat[2][11][3];
    }
    double GetDinuc31_GTCA()
    {
        return dinuc_stat[2][11][4];
    }
    double GetDinuc31_GTCC()
    {
        return dinuc_stat[2][11][5];
    }
    double GetDinuc31_GTCG()
    {
        return dinuc_stat[2][11][6];
    }
    double GetDinuc31_GTCT()
    {
        return dinuc_stat[2][11][7];
    }
    double GetDinuc31_GTGA()
    {
        return dinuc_stat[2][11][8];
    }
    double GetDinuc31_GTGC()
    {
        return dinuc_stat[2][11][9];
    }
    double GetDinuc31_GTGG()
    {
        return dinuc_stat[2][11][10];
    }
    double GetDinuc31_GTGT()
    {
        return dinuc_stat[2][11][11];
    }
    double GetDinuc31_GTTA()
    {
        return dinuc_stat[2][11][12];
    }
    double GetDinuc31_GTTC()
    {
        return dinuc_stat[2][11][13];
    }
    double GetDinuc31_GTTG()
    {
        return dinuc_stat[2][11][14];
    }
    double GetDinuc31_GTTT()
    {
        return dinuc_stat[2][11][15];
    }
    double GetDinuc31_TAAA()
    {
        return dinuc_stat[2][12][0];
    }
    double GetDinuc31_TAAC()
    {
        return dinuc_stat[2][12][1];
    }
    double GetDinuc31_TAAG()
    {
        return dinuc_stat[2][12][2];
    }
    double GetDinuc31_TAAT()
    {
        return dinuc_stat[2][12][3];
    }
    double GetDinuc31_TACA()
    {
        return dinuc_stat[2][12][4];
    }
    double GetDinuc31_TACC()
    {
        return dinuc_stat[2][12][5];
    }
    double GetDinuc31_TACG()
    {
        return dinuc_stat[2][12][6];
    }
    double GetDinuc31_TACT()
    {
        return dinuc_stat[2][12][7];
    }
    double GetDinuc31_TAGA()
    {
        return dinuc_stat[2][12][8];
    }
    double GetDinuc31_TAGC()
    {
        return dinuc_stat[2][12][9];
    }
    double GetDinuc31_TAGG()
    {
        return dinuc_stat[2][12][10];
    }
    double GetDinuc31_TAGT()
    {
        return dinuc_stat[2][12][11];
    }
    double GetDinuc31_TATA()
    {
        return dinuc_stat[2][12][12];
    }
    double GetDinuc31_TATC()
    {
        return dinuc_stat[2][12][13];
    }
    double GetDinuc31_TATG()
    {
        return dinuc_stat[2][12][14];
    }
    double GetDinuc31_TATT()
    {
        return dinuc_stat[2][12][15];
    }
    double GetDinuc31_TCAA()
    {
        return dinuc_stat[2][13][0];
    }
    double GetDinuc31_TCAC()
    {
        return dinuc_stat[2][13][1];
    }
    double GetDinuc31_TCAG()
    {
        return dinuc_stat[2][13][2];
    }
    double GetDinuc31_TCAT()
    {
        return dinuc_stat[2][13][3];
    }
    double GetDinuc31_TCCA()
    {
        return dinuc_stat[2][13][4];
    }
    double GetDinuc31_TCCC()
    {
        return dinuc_stat[2][13][5];
    }
    double GetDinuc31_TCCG()
    {
        return dinuc_stat[2][13][6];
    }
    double GetDinuc31_TCCT()
    {
        return dinuc_stat[2][13][7];
    }
    double GetDinuc31_TCGA()
    {
        return dinuc_stat[2][13][8];
    }
    double GetDinuc31_TCGC()
    {
        return dinuc_stat[2][13][9];
    }
    double GetDinuc31_TCGG()
    {
        return dinuc_stat[2][13][10];
    }
    double GetDinuc31_TCGT()
    {
        return dinuc_stat[2][13][11];
    }
    double GetDinuc31_TCTA()
    {
        return dinuc_stat[2][13][12];
    }
    double GetDinuc31_TCTC()
    {
        return dinuc_stat[2][13][13];
    }
    double GetDinuc31_TCTG()
    {
        return dinuc_stat[2][13][14];
    }
    double GetDinuc31_TCTT()
    {
        return dinuc_stat[2][13][15];
    }
    double GetDinuc31_TGAA()
    {
        return dinuc_stat[2][14][0];
    }
    double GetDinuc31_TGAC()
    {
        return dinuc_stat[2][14][1];
    }
    double GetDinuc31_TGAG()
    {
        return dinuc_stat[2][14][2];
    }
    double GetDinuc31_TGAT()
    {
        return dinuc_stat[2][14][3];
    }
    double GetDinuc31_TGCA()
    {
        return dinuc_stat[2][14][4];
    }
    double GetDinuc31_TGCC()
    {
        return dinuc_stat[2][14][5];
    }
    double GetDinuc31_TGCG()
    {
        return dinuc_stat[2][14][6];
    }
    double GetDinuc31_TGCT()
    {
        return dinuc_stat[2][14][7];
    }
    double GetDinuc31_TGGA()
    {
        return dinuc_stat[2][14][8];
    }
    double GetDinuc31_TGGC()
    {
        return dinuc_stat[2][14][9];
    }
    double GetDinuc31_TGGG()
    {
        return dinuc_stat[2][14][10];
    }
    double GetDinuc31_TGGT()
    {
        return dinuc_stat[2][14][11];
    }
    double GetDinuc31_TGTA()
    {
        return dinuc_stat[2][14][12];
    }
    double GetDinuc31_TGTC()
    {
        return dinuc_stat[2][14][13];
    }
    double GetDinuc31_TGTG()
    {
        return dinuc_stat[2][14][14];
    }
    double GetDinuc31_TGTT()
    {
        return dinuc_stat[2][14][15];
    }
    double GetDinuc31_TTAA()
    {
        return dinuc_stat[2][15][0];
    }
    double GetDinuc31_TTAC()
    {
        return dinuc_stat[2][15][1];
    }
    double GetDinuc31_TTAG()
    {
        return dinuc_stat[2][15][2];
    }
    double GetDinuc31_TTAT()
    {
        return dinuc_stat[2][15][3];
    }
    double GetDinuc31_TTCA()
    {
        return dinuc_stat[2][15][4];
    }
    double GetDinuc31_TTCC()
    {
        return dinuc_stat[2][15][5];
    }
    double GetDinuc31_TTCG()
    {
        return dinuc_stat[2][15][6];
    }
    double GetDinuc31_TTCT()
    {
        return dinuc_stat[2][15][7];
    }
    double GetDinuc31_TTGA()
    {
        return dinuc_stat[2][15][8];
    }
    double GetDinuc31_TTGC()
    {
        return dinuc_stat[2][15][9];
    }
    double GetDinuc31_TTGG()
    {
        return dinuc_stat[2][15][10];
    }
    double GetDinuc31_TTGT()
    {
        return dinuc_stat[2][15][11];
    }
    double GetDinuc31_TTTA()
    {
        return dinuc_stat[2][15][12];
    }
    double GetDinuc31_TTTC()
    {
        return dinuc_stat[2][15][13];
    }
    double GetDinuc31_TTTG()
    {
        return dinuc_stat[2][15][14];
    }
    double GetDinuc31_TTTT()
    {
        return dinuc_stat[2][15][15];
    }

    double GetDinucSyn_AAAA()
    {
        return dinucSyn_stat[3][0][0];
    }
    double GetDinucSyn_AAAC()
    {
        return dinucSyn_stat[3][0][1];
    }
    double GetDinucSyn_AAAG()
    {
        return dinucSyn_stat[3][0][2];
    }
    double GetDinucSyn_AAAT()
    {
        return dinucSyn_stat[3][0][3];
    }
    double GetDinucSyn_AACA()
    {
        return dinucSyn_stat[3][0][4];
    }
    double GetDinucSyn_AACC()
    {
        return dinucSyn_stat[3][0][5];
    }
    double GetDinucSyn_AACG()
    {
        return dinucSyn_stat[3][0][6];
    }
    double GetDinucSyn_AACT()
    {
        return dinucSyn_stat[3][0][7];
    }
    double GetDinucSyn_AAGA()
    {
        return dinucSyn_stat[3][0][8];
    }
    double GetDinucSyn_AAGC()
    {
        return dinucSyn_stat[3][0][9];
    }
    double GetDinucSyn_AAGG()
    {
        return dinucSyn_stat[3][0][10];
    }
    double GetDinucSyn_AAGT()
    {
        return dinucSyn_stat[3][0][11];
    }
    double GetDinucSyn_AATA()
    {
        return dinucSyn_stat[3][0][12];
    }
    double GetDinucSyn_AATC()
    {
        return dinucSyn_stat[3][0][13];
    }
    double GetDinucSyn_AATG()
    {
        return dinucSyn_stat[3][0][14];
    }
    double GetDinucSyn_AATT()
    {
        return dinucSyn_stat[3][0][15];
    }
    double GetDinucSyn_ACAA()
    {
        return dinucSyn_stat[3][1][0];
    }
    double GetDinucSyn_ACAC()
    {
        return dinucSyn_stat[3][1][1];
    }
    double GetDinucSyn_ACAG()
    {
        return dinucSyn_stat[3][1][2];
    }
    double GetDinucSyn_ACAT()
    {
        return dinucSyn_stat[3][1][3];
    }
    double GetDinucSyn_ACCA()
    {
        return dinucSyn_stat[3][1][4];
    }
    double GetDinucSyn_ACCC()
    {
        return dinucSyn_stat[3][1][5];
    }
    double GetDinucSyn_ACCG()
    {
        return dinucSyn_stat[3][1][6];
    }
    double GetDinucSyn_ACCT()
    {
        return dinucSyn_stat[3][1][7];
    }
    double GetDinucSyn_ACGA()
    {
        return dinucSyn_stat[3][1][8];
    }
    double GetDinucSyn_ACGC()
    {
        return dinucSyn_stat[3][1][9];
    }
    double GetDinucSyn_ACGG()
    {
        return dinucSyn_stat[3][1][10];
    }
    double GetDinucSyn_ACGT()
    {
        return dinucSyn_stat[3][1][11];
    }
    double GetDinucSyn_ACTA()
    {
        return dinucSyn_stat[3][1][12];
    }
    double GetDinucSyn_ACTC()
    {
        return dinucSyn_stat[3][1][13];
    }
    double GetDinucSyn_ACTG()
    {
        return dinucSyn_stat[3][1][14];
    }
    double GetDinucSyn_ACTT()
    {
        return dinucSyn_stat[3][1][15];
    }
    double GetDinucSyn_AGAA()
    {
        return dinucSyn_stat[3][2][0];
    }
    double GetDinucSyn_AGAC()
    {
        return dinucSyn_stat[3][2][1];
    }
    double GetDinucSyn_AGAG()
    {
        return dinucSyn_stat[3][2][2];
    }
    double GetDinucSyn_AGAT()
    {
        return dinucSyn_stat[3][2][3];
    }
    double GetDinucSyn_AGCA()
    {
        return dinucSyn_stat[3][2][4];
    }
    double GetDinucSyn_AGCC()
    {
        return dinucSyn_stat[3][2][5];
    }
    double GetDinucSyn_AGCG()
    {
        return dinucSyn_stat[3][2][6];
    }
    double GetDinucSyn_AGCT()
    {
        return dinucSyn_stat[3][2][7];
    }
    double GetDinucSyn_AGGA()
    {
        return dinucSyn_stat[3][2][8];
    }
    double GetDinucSyn_AGGC()
    {
        return dinucSyn_stat[3][2][9];
    }
    double GetDinucSyn_AGGG()
    {
        return dinucSyn_stat[3][2][10];
    }
    double GetDinucSyn_AGGT()
    {
        return dinucSyn_stat[3][2][11];
    }
    double GetDinucSyn_AGTA()
    {
        return dinucSyn_stat[3][2][12];
    }
    double GetDinucSyn_AGTC()
    {
        return dinucSyn_stat[3][2][13];
    }
    double GetDinucSyn_AGTG()
    {
        return dinucSyn_stat[3][2][14];
    }
    double GetDinucSyn_AGTT()
    {
        return dinucSyn_stat[3][2][15];
    }
    double GetDinucSyn_ATAA()
    {
        return dinucSyn_stat[3][3][0];
    }
    double GetDinucSyn_ATAC()
    {
        return dinucSyn_stat[3][3][1];
    }
    double GetDinucSyn_ATAG()
    {
        return dinucSyn_stat[3][3][2];
    }
    double GetDinucSyn_ATAT()
    {
        return dinucSyn_stat[3][3][3];
    }
    double GetDinucSyn_ATCA()
    {
        return dinucSyn_stat[3][3][4];
    }
    double GetDinucSyn_ATCC()
    {
        return dinucSyn_stat[3][3][5];
    }
    double GetDinucSyn_ATCG()
    {
        return dinucSyn_stat[3][3][6];
    }
    double GetDinucSyn_ATCT()
    {
        return dinucSyn_stat[3][3][7];
    }
    double GetDinucSyn_ATGA()
    {
        return dinucSyn_stat[3][3][8];
    }
    double GetDinucSyn_ATGC()
    {
        return dinucSyn_stat[3][3][9];
    }
    double GetDinucSyn_ATGG()
    {
        return dinucSyn_stat[3][3][10];
    }
    double GetDinucSyn_ATGT()
    {
        return dinucSyn_stat[3][3][11];
    }
    double GetDinucSyn_ATTA()
    {
        return dinucSyn_stat[3][3][12];
    }
    double GetDinucSyn_ATTC()
    {
        return dinucSyn_stat[3][3][13];
    }
    double GetDinucSyn_ATTG()
    {
        return dinucSyn_stat[3][3][14];
    }
    double GetDinucSyn_ATTT()
    {
        return dinucSyn_stat[3][3][15];
    }
    double GetDinucSyn_CAAA()
    {
        return dinucSyn_stat[3][4][0];
    }
    double GetDinucSyn_CAAC()
    {
        return dinucSyn_stat[3][4][1];
    }
    double GetDinucSyn_CAAG()
    {
        return dinucSyn_stat[3][4][2];
    }
    double GetDinucSyn_CAAT()
    {
        return dinucSyn_stat[3][4][3];
    }
    double GetDinucSyn_CACA()
    {
        return dinucSyn_stat[3][4][4];
    }
    double GetDinucSyn_CACC()
    {
        return dinucSyn_stat[3][4][5];
    }
    double GetDinucSyn_CACG()
    {
        return dinucSyn_stat[3][4][6];
    }
    double GetDinucSyn_CACT()
    {
        return dinucSyn_stat[3][4][7];
    }
    double GetDinucSyn_CAGA()
    {
        return dinucSyn_stat[3][4][8];
    }
    double GetDinucSyn_CAGC()
    {
        return dinucSyn_stat[3][4][9];
    }
    double GetDinucSyn_CAGG()
    {
        return dinucSyn_stat[3][4][10];
    }
    double GetDinucSyn_CAGT()
    {
        return dinucSyn_stat[3][4][11];
    }
    double GetDinucSyn_CATA()
    {
        return dinucSyn_stat[3][4][12];
    }
    double GetDinucSyn_CATC()
    {
        return dinucSyn_stat[3][4][13];
    }
    double GetDinucSyn_CATG()
    {
        return dinucSyn_stat[3][4][14];
    }
    double GetDinucSyn_CATT()
    {
        return dinucSyn_stat[3][4][15];
    }
    double GetDinucSyn_CCAA()
    {
        return dinucSyn_stat[3][5][0];
    }
    double GetDinucSyn_CCAC()
    {
        return dinucSyn_stat[3][5][1];
    }
    double GetDinucSyn_CCAG()
    {
        return dinucSyn_stat[3][5][2];
    }
    double GetDinucSyn_CCAT()
    {
        return dinucSyn_stat[3][5][3];
    }
    double GetDinucSyn_CCCA()
    {
        return dinucSyn_stat[3][5][4];
    }
    double GetDinucSyn_CCCC()
    {
        return dinucSyn_stat[3][5][5];
    }
    double GetDinucSyn_CCCG()
    {
        return dinucSyn_stat[3][5][6];
    }
    double GetDinucSyn_CCCT()
    {
        return dinucSyn_stat[3][5][7];
    }
    double GetDinucSyn_CCGA()
    {
        return dinucSyn_stat[3][5][8];
    }
    double GetDinucSyn_CCGC()
    {
        return dinucSyn_stat[3][5][9];
    }
    double GetDinucSyn_CCGG()
    {
        return dinucSyn_stat[3][5][10];
    }
    double GetDinucSyn_CCGT()
    {
        return dinucSyn_stat[3][5][11];
    }
    double GetDinucSyn_CCTA()
    {
        return dinucSyn_stat[3][5][12];
    }
    double GetDinucSyn_CCTC()
    {
        return dinucSyn_stat[3][5][13];
    }
    double GetDinucSyn_CCTG()
    {
        return dinucSyn_stat[3][5][14];
    }
    double GetDinucSyn_CCTT()
    {
        return dinucSyn_stat[3][5][15];
    }
    double GetDinucSyn_CGAA()
    {
        return dinucSyn_stat[3][6][0];
    }
    double GetDinucSyn_CGAC()
    {
        return dinucSyn_stat[3][6][1];
    }
    double GetDinucSyn_CGAG()
    {
        return dinucSyn_stat[3][6][2];
    }
    double GetDinucSyn_CGAT()
    {
        return dinucSyn_stat[3][6][3];
    }
    double GetDinucSyn_CGCA()
    {
        return dinucSyn_stat[3][6][4];
    }
    double GetDinucSyn_CGCC()
    {
        return dinucSyn_stat[3][6][5];
    }
    double GetDinucSyn_CGCG()
    {
        return dinucSyn_stat[3][6][6];
    }
    double GetDinucSyn_CGCT()
    {
        return dinucSyn_stat[3][6][7];
    }
    double GetDinucSyn_CGGA()
    {
        return dinucSyn_stat[3][6][8];
    }
    double GetDinucSyn_CGGC()
    {
        return dinucSyn_stat[3][6][9];
    }
    double GetDinucSyn_CGGG()
    {
        return dinucSyn_stat[3][6][10];
    }
    double GetDinucSyn_CGGT()
    {
        return dinucSyn_stat[3][6][11];
    }
    double GetDinucSyn_CGTA()
    {
        return dinucSyn_stat[3][6][12];
    }
    double GetDinucSyn_CGTC()
    {
        return dinucSyn_stat[3][6][13];
    }
    double GetDinucSyn_CGTG()
    {
        return dinucSyn_stat[3][6][14];
    }
    double GetDinucSyn_CGTT()
    {
        return dinucSyn_stat[3][6][15];
    }
    double GetDinucSyn_CTAA()
    {
        return dinucSyn_stat[3][7][0];
    }
    double GetDinucSyn_CTAC()
    {
        return dinucSyn_stat[3][7][1];
    }
    double GetDinucSyn_CTAG()
    {
        return dinucSyn_stat[3][7][2];
    }
    double GetDinucSyn_CTAT()
    {
        return dinucSyn_stat[3][7][3];
    }
    double GetDinucSyn_CTCA()
    {
        return dinucSyn_stat[3][7][4];
    }
    double GetDinucSyn_CTCC()
    {
        return dinucSyn_stat[3][7][5];
    }
    double GetDinucSyn_CTCG()
    {
        return dinucSyn_stat[3][7][6];
    }
    double GetDinucSyn_CTCT()
    {
        return dinucSyn_stat[3][7][7];
    }
    double GetDinucSyn_CTGA()
    {
        return dinucSyn_stat[3][7][8];
    }
    double GetDinucSyn_CTGC()
    {
        return dinucSyn_stat[3][7][9];
    }
    double GetDinucSyn_CTGG()
    {
        return dinucSyn_stat[3][7][10];
    }
    double GetDinucSyn_CTGT()
    {
        return dinucSyn_stat[3][7][11];
    }
    double GetDinucSyn_CTTA()
    {
        return dinucSyn_stat[3][7][12];
    }
    double GetDinucSyn_CTTC()
    {
        return dinucSyn_stat[3][7][13];
    }
    double GetDinucSyn_CTTG()
    {
        return dinucSyn_stat[3][7][14];
    }
    double GetDinucSyn_CTTT()
    {
        return dinucSyn_stat[3][7][15];
    }
    double GetDinucSyn_GAAA()
    {
        return dinucSyn_stat[3][8][0];
    }
    double GetDinucSyn_GAAC()
    {
        return dinucSyn_stat[3][8][1];
    }
    double GetDinucSyn_GAAG()
    {
        return dinucSyn_stat[3][8][2];
    }
    double GetDinucSyn_GAAT()
    {
        return dinucSyn_stat[3][8][3];
    }
    double GetDinucSyn_GACA()
    {
        return dinucSyn_stat[3][8][4];
    }
    double GetDinucSyn_GACC()
    {
        return dinucSyn_stat[3][8][5];
    }
    double GetDinucSyn_GACG()
    {
        return dinucSyn_stat[3][8][6];
    }
    double GetDinucSyn_GACT()
    {
        return dinucSyn_stat[3][8][7];
    }
    double GetDinucSyn_GAGA()
    {
        return dinucSyn_stat[3][8][8];
    }
    double GetDinucSyn_GAGC()
    {
        return dinucSyn_stat[3][8][9];
    }
    double GetDinucSyn_GAGG()
    {
        return dinucSyn_stat[3][8][10];
    }
    double GetDinucSyn_GAGT()
    {
        return dinucSyn_stat[3][8][11];
    }
    double GetDinucSyn_GATA()
    {
        return dinucSyn_stat[3][8][12];
    }
    double GetDinucSyn_GATC()
    {
        return dinucSyn_stat[3][8][13];
    }
    double GetDinucSyn_GATG()
    {
        return dinucSyn_stat[3][8][14];
    }
    double GetDinucSyn_GATT()
    {
        return dinucSyn_stat[3][8][15];
    }
    double GetDinucSyn_GCAA()
    {
        return dinucSyn_stat[3][9][0];
    }
    double GetDinucSyn_GCAC()
    {
        return dinucSyn_stat[3][9][1];
    }
    double GetDinucSyn_GCAG()
    {
        return dinucSyn_stat[3][9][2];
    }
    double GetDinucSyn_GCAT()
    {
        return dinucSyn_stat[3][9][3];
    }
    double GetDinucSyn_GCCA()
    {
        return dinucSyn_stat[3][9][4];
    }
    double GetDinucSyn_GCCC()
    {
        return dinucSyn_stat[3][9][5];
    }
    double GetDinucSyn_GCCG()
    {
        return dinucSyn_stat[3][9][6];
    }
    double GetDinucSyn_GCCT()
    {
        return dinucSyn_stat[3][9][7];
    }
    double GetDinucSyn_GCGA()
    {
        return dinucSyn_stat[3][9][8];
    }
    double GetDinucSyn_GCGC()
    {
        return dinucSyn_stat[3][9][9];
    }
    double GetDinucSyn_GCGG()
    {
        return dinucSyn_stat[3][9][10];
    }
    double GetDinucSyn_GCGT()
    {
        return dinucSyn_stat[3][9][11];
    }
    double GetDinucSyn_GCTA()
    {
        return dinucSyn_stat[3][9][12];
    }
    double GetDinucSyn_GCTC()
    {
        return dinucSyn_stat[3][9][13];
    }
    double GetDinucSyn_GCTG()
    {
        return dinucSyn_stat[3][9][14];
    }
    double GetDinucSyn_GCTT()
    {
        return dinucSyn_stat[3][9][15];
    }
    double GetDinucSyn_GGAA()
    {
        return dinucSyn_stat[3][10][0];
    }
    double GetDinucSyn_GGAC()
    {
        return dinucSyn_stat[3][10][1];
    }
    double GetDinucSyn_GGAG()
    {
        return dinucSyn_stat[3][10][2];
    }
    double GetDinucSyn_GGAT()
    {
        return dinucSyn_stat[3][10][3];
    }
    double GetDinucSyn_GGCA()
    {
        return dinucSyn_stat[3][10][4];
    }
    double GetDinucSyn_GGCC()
    {
        return dinucSyn_stat[3][10][5];
    }
    double GetDinucSyn_GGCG()
    {
        return dinucSyn_stat[3][10][6];
    }
    double GetDinucSyn_GGCT()
    {
        return dinucSyn_stat[3][10][7];
    }
    double GetDinucSyn_GGGA()
    {
        return dinucSyn_stat[3][10][8];
    }
    double GetDinucSyn_GGGC()
    {
        return dinucSyn_stat[3][10][9];
    }
    double GetDinucSyn_GGGG()
    {
        return dinucSyn_stat[3][10][10];
    }
    double GetDinucSyn_GGGT()
    {
        return dinucSyn_stat[3][10][11];
    }
    double GetDinucSyn_GGTA()
    {
        return dinucSyn_stat[3][10][12];
    }
    double GetDinucSyn_GGTC()
    {
        return dinucSyn_stat[3][10][13];
    }
    double GetDinucSyn_GGTG()
    {
        return dinucSyn_stat[3][10][14];
    }
    double GetDinucSyn_GGTT()
    {
        return dinucSyn_stat[3][10][15];
    }
    double GetDinucSyn_GTAA()
    {
        return dinucSyn_stat[3][11][0];
    }
    double GetDinucSyn_GTAC()
    {
        return dinucSyn_stat[3][11][1];
    }
    double GetDinucSyn_GTAG()
    {
        return dinucSyn_stat[3][11][2];
    }
    double GetDinucSyn_GTAT()
    {
        return dinucSyn_stat[3][11][3];
    }
    double GetDinucSyn_GTCA()
    {
        return dinucSyn_stat[3][11][4];
    }
    double GetDinucSyn_GTCC()
    {
        return dinucSyn_stat[3][11][5];
    }
    double GetDinucSyn_GTCG()
    {
        return dinucSyn_stat[3][11][6];
    }
    double GetDinucSyn_GTCT()
    {
        return dinucSyn_stat[3][11][7];
    }
    double GetDinucSyn_GTGA()
    {
        return dinucSyn_stat[3][11][8];
    }
    double GetDinucSyn_GTGC()
    {
        return dinucSyn_stat[3][11][9];
    }
    double GetDinucSyn_GTGG()
    {
        return dinucSyn_stat[3][11][10];
    }
    double GetDinucSyn_GTGT()
    {
        return dinucSyn_stat[3][11][11];
    }
    double GetDinucSyn_GTTA()
    {
        return dinucSyn_stat[3][11][12];
    }
    double GetDinucSyn_GTTC()
    {
        return dinucSyn_stat[3][11][13];
    }
    double GetDinucSyn_GTTG()
    {
        return dinucSyn_stat[3][11][14];
    }
    double GetDinucSyn_GTTT()
    {
        return dinucSyn_stat[3][11][15];
    }
    double GetDinucSyn_TAAA()
    {
        return dinucSyn_stat[3][12][0];
    }
    double GetDinucSyn_TAAC()
    {
        return dinucSyn_stat[3][12][1];
    }
    double GetDinucSyn_TAAG()
    {
        return dinucSyn_stat[3][12][2];
    }
    double GetDinucSyn_TAAT()
    {
        return dinucSyn_stat[3][12][3];
    }
    double GetDinucSyn_TACA()
    {
        return dinucSyn_stat[3][12][4];
    }
    double GetDinucSyn_TACC()
    {
        return dinucSyn_stat[3][12][5];
    }
    double GetDinucSyn_TACG()
    {
        return dinucSyn_stat[3][12][6];
    }
    double GetDinucSyn_TACT()
    {
        return dinucSyn_stat[3][12][7];
    }
    double GetDinucSyn_TAGA()
    {
        return dinucSyn_stat[3][12][8];
    }
    double GetDinucSyn_TAGC()
    {
        return dinucSyn_stat[3][12][9];
    }
    double GetDinucSyn_TAGG()
    {
        return dinucSyn_stat[3][12][10];
    }
    double GetDinucSyn_TAGT()
    {
        return dinucSyn_stat[3][12][11];
    }
    double GetDinucSyn_TATA()
    {
        return dinucSyn_stat[3][12][12];
    }
    double GetDinucSyn_TATC()
    {
        return dinucSyn_stat[3][12][13];
    }
    double GetDinucSyn_TATG()
    {
        return dinucSyn_stat[3][12][14];
    }
    double GetDinucSyn_TATT()
    {
        return dinucSyn_stat[3][12][15];
    }
    double GetDinucSyn_TCAA()
    {
        return dinucSyn_stat[3][13][0];
    }
    double GetDinucSyn_TCAC()
    {
        return dinucSyn_stat[3][13][1];
    }
    double GetDinucSyn_TCAG()
    {
        return dinucSyn_stat[3][13][2];
    }
    double GetDinucSyn_TCAT()
    {
        return dinucSyn_stat[3][13][3];
    }
    double GetDinucSyn_TCCA()
    {
        return dinucSyn_stat[3][13][4];
    }
    double GetDinucSyn_TCCC()
    {
        return dinucSyn_stat[3][13][5];
    }
    double GetDinucSyn_TCCG()
    {
        return dinucSyn_stat[3][13][6];
    }
    double GetDinucSyn_TCCT()
    {
        return dinucSyn_stat[3][13][7];
    }
    double GetDinucSyn_TCGA()
    {
        return dinucSyn_stat[3][13][8];
    }
    double GetDinucSyn_TCGC()
    {
        return dinucSyn_stat[3][13][9];
    }
    double GetDinucSyn_TCGG()
    {
        return dinucSyn_stat[3][13][10];
    }
    double GetDinucSyn_TCGT()
    {
        return dinucSyn_stat[3][13][11];
    }
    double GetDinucSyn_TCTA()
    {
        return dinucSyn_stat[3][13][12];
    }
    double GetDinucSyn_TCTC()
    {
        return dinucSyn_stat[3][13][13];
    }
    double GetDinucSyn_TCTG()
    {
        return dinucSyn_stat[3][13][14];
    }
    double GetDinucSyn_TCTT()
    {
        return dinucSyn_stat[3][13][15];
    }
    double GetDinucSyn_TGAA()
    {
        return dinucSyn_stat[3][14][0];
    }
    double GetDinucSyn_TGAC()
    {
        return dinucSyn_stat[3][14][1];
    }
    double GetDinucSyn_TGAG()
    {
        return dinucSyn_stat[3][14][2];
    }
    double GetDinucSyn_TGAT()
    {
        return dinucSyn_stat[3][14][3];
    }
    double GetDinucSyn_TGCA()
    {
        return dinucSyn_stat[3][14][4];
    }
    double GetDinucSyn_TGCC()
    {
        return dinucSyn_stat[3][14][5];
    }
    double GetDinucSyn_TGCG()
    {
        return dinucSyn_stat[3][14][6];
    }
    double GetDinucSyn_TGCT()
    {
        return dinucSyn_stat[3][14][7];
    }
    double GetDinucSyn_TGGA()
    {
        return dinucSyn_stat[3][14][8];
    }
    double GetDinucSyn_TGGC()
    {
        return dinucSyn_stat[3][14][9];
    }
    double GetDinucSyn_TGGG()
    {
        return dinucSyn_stat[3][14][10];
    }
    double GetDinucSyn_TGGT()
    {
        return dinucSyn_stat[3][14][11];
    }
    double GetDinucSyn_TGTA()
    {
        return dinucSyn_stat[3][14][12];
    }
    double GetDinucSyn_TGTC()
    {
        return dinucSyn_stat[3][14][13];
    }
    double GetDinucSyn_TGTG()
    {
        return dinucSyn_stat[3][14][14];
    }
    double GetDinucSyn_TGTT()
    {
        return dinucSyn_stat[3][14][15];
    }
    double GetDinucSyn_TTAA()
    {
        return dinucSyn_stat[3][15][0];
    }
    double GetDinucSyn_TTAC()
    {
        return dinucSyn_stat[3][15][1];
    }
    double GetDinucSyn_TTAG()
    {
        return dinucSyn_stat[3][15][2];
    }
    double GetDinucSyn_TTAT()
    {
        return dinucSyn_stat[3][15][3];
    }
    double GetDinucSyn_TTCA()
    {
        return dinucSyn_stat[3][15][4];
    }
    double GetDinucSyn_TTCC()
    {
        return dinucSyn_stat[3][15][5];
    }
    double GetDinucSyn_TTCG()
    {
        return dinucSyn_stat[3][15][6];
    }
    double GetDinucSyn_TTCT()
    {
        return dinucSyn_stat[3][15][7];
    }
    double GetDinucSyn_TTGA()
    {
        return dinucSyn_stat[3][15][8];
    }
    double GetDinucSyn_TTGC()
    {
        return dinucSyn_stat[3][15][9];
    }
    double GetDinucSyn_TTGG()
    {
        return dinucSyn_stat[3][15][10];
    }
    double GetDinucSyn_TTGT()
    {
        return dinucSyn_stat[3][15][11];
    }
    double GetDinucSyn_TTTA()
    {
        return dinucSyn_stat[3][15][12];
    }
    double GetDinucSyn_TTTC()
    {
        return dinucSyn_stat[3][15][13];
    }
    double GetDinucSyn_TTTG()
    {
        return dinucSyn_stat[3][15][14];
    }
    double GetDinucSyn_TTTT()
    {
        return dinucSyn_stat[3][15][15];
    }
    double GetDinucSyn12_AAAA()
    {
        return dinucSyn_stat[0][0][0];
    }
    double GetDinucSyn12_AAAC()
    {
        return dinucSyn_stat[0][0][1];
    }
    double GetDinucSyn12_AAAG()
    {
        return dinucSyn_stat[0][0][2];
    }
    double GetDinucSyn12_AAAT()
    {
        return dinucSyn_stat[0][0][3];
    }
    double GetDinucSyn12_AACA()
    {
        return dinucSyn_stat[0][0][4];
    }
    double GetDinucSyn12_AACC()
    {
        return dinucSyn_stat[0][0][5];
    }
    double GetDinucSyn12_AACG()
    {
        return dinucSyn_stat[0][0][6];
    }
    double GetDinucSyn12_AACT()
    {
        return dinucSyn_stat[0][0][7];
    }
    double GetDinucSyn12_AAGA()
    {
        return dinucSyn_stat[0][0][8];
    }
    double GetDinucSyn12_AAGC()
    {
        return dinucSyn_stat[0][0][9];
    }
    double GetDinucSyn12_AAGG()
    {
        return dinucSyn_stat[0][0][10];
    }
    double GetDinucSyn12_AAGT()
    {
        return dinucSyn_stat[0][0][11];
    }
    double GetDinucSyn12_AATA()
    {
        return dinucSyn_stat[0][0][12];
    }
    double GetDinucSyn12_AATC()
    {
        return dinucSyn_stat[0][0][13];
    }
    double GetDinucSyn12_AATG()
    {
        return dinucSyn_stat[0][0][14];
    }
    double GetDinucSyn12_AATT()
    {
        return dinucSyn_stat[0][0][15];
    }
    double GetDinucSyn12_ACAA()
    {
        return dinucSyn_stat[0][1][0];
    }
    double GetDinucSyn12_ACAC()
    {
        return dinucSyn_stat[0][1][1];
    }
    double GetDinucSyn12_ACAG()
    {
        return dinucSyn_stat[0][1][2];
    }
    double GetDinucSyn12_ACAT()
    {
        return dinucSyn_stat[0][1][3];
    }
    double GetDinucSyn12_ACCA()
    {
        return dinucSyn_stat[0][1][4];
    }
    double GetDinucSyn12_ACCC()
    {
        return dinucSyn_stat[0][1][5];
    }
    double GetDinucSyn12_ACCG()
    {
        return dinucSyn_stat[0][1][6];
    }
    double GetDinucSyn12_ACCT()
    {
        return dinucSyn_stat[0][1][7];
    }
    double GetDinucSyn12_ACGA()
    {
        return dinucSyn_stat[0][1][8];
    }
    double GetDinucSyn12_ACGC()
    {
        return dinucSyn_stat[0][1][9];
    }
    double GetDinucSyn12_ACGG()
    {
        return dinucSyn_stat[0][1][10];
    }
    double GetDinucSyn12_ACGT()
    {
        return dinucSyn_stat[0][1][11];
    }
    double GetDinucSyn12_ACTA()
    {
        return dinucSyn_stat[0][1][12];
    }
    double GetDinucSyn12_ACTC()
    {
        return dinucSyn_stat[0][1][13];
    }
    double GetDinucSyn12_ACTG()
    {
        return dinucSyn_stat[0][1][14];
    }
    double GetDinucSyn12_ACTT()
    {
        return dinucSyn_stat[0][1][15];
    }
    double GetDinucSyn12_AGAA()
    {
        return dinucSyn_stat[0][2][0];
    }
    double GetDinucSyn12_AGAC()
    {
        return dinucSyn_stat[0][2][1];
    }
    double GetDinucSyn12_AGAG()
    {
        return dinucSyn_stat[0][2][2];
    }
    double GetDinucSyn12_AGAT()
    {
        return dinucSyn_stat[0][2][3];
    }
    double GetDinucSyn12_AGCA()
    {
        return dinucSyn_stat[0][2][4];
    }
    double GetDinucSyn12_AGCC()
    {
        return dinucSyn_stat[0][2][5];
    }
    double GetDinucSyn12_AGCG()
    {
        return dinucSyn_stat[0][2][6];
    }
    double GetDinucSyn12_AGCT()
    {
        return dinucSyn_stat[0][2][7];
    }
    double GetDinucSyn12_AGGA()
    {
        return dinucSyn_stat[0][2][8];
    }
    double GetDinucSyn12_AGGC()
    {
        return dinucSyn_stat[0][2][9];
    }
    double GetDinucSyn12_AGGG()
    {
        return dinucSyn_stat[0][2][10];
    }
    double GetDinucSyn12_AGGT()
    {
        return dinucSyn_stat[0][2][11];
    }
    double GetDinucSyn12_AGTA()
    {
        return dinucSyn_stat[0][2][12];
    }
    double GetDinucSyn12_AGTC()
    {
        return dinucSyn_stat[0][2][13];
    }
    double GetDinucSyn12_AGTG()
    {
        return dinucSyn_stat[0][2][14];
    }
    double GetDinucSyn12_AGTT()
    {
        return dinucSyn_stat[0][2][15];
    }
    double GetDinucSyn12_ATAA()
    {
        return dinucSyn_stat[0][3][0];
    }
    double GetDinucSyn12_ATAC()
    {
        return dinucSyn_stat[0][3][1];
    }
    double GetDinucSyn12_ATAG()
    {
        return dinucSyn_stat[0][3][2];
    }
    double GetDinucSyn12_ATAT()
    {
        return dinucSyn_stat[0][3][3];
    }
    double GetDinucSyn12_ATCA()
    {
        return dinucSyn_stat[0][3][4];
    }
    double GetDinucSyn12_ATCC()
    {
        return dinucSyn_stat[0][3][5];
    }
    double GetDinucSyn12_ATCG()
    {
        return dinucSyn_stat[0][3][6];
    }
    double GetDinucSyn12_ATCT()
    {
        return dinucSyn_stat[0][3][7];
    }
    double GetDinucSyn12_ATGA()
    {
        return dinucSyn_stat[0][3][8];
    }
    double GetDinucSyn12_ATGC()
    {
        return dinucSyn_stat[0][3][9];
    }
    double GetDinucSyn12_ATGG()
    {
        return dinucSyn_stat[0][3][10];
    }
    double GetDinucSyn12_ATGT()
    {
        return dinucSyn_stat[0][3][11];
    }
    double GetDinucSyn12_ATTA()
    {
        return dinucSyn_stat[0][3][12];
    }
    double GetDinucSyn12_ATTC()
    {
        return dinucSyn_stat[0][3][13];
    }
    double GetDinucSyn12_ATTG()
    {
        return dinucSyn_stat[0][3][14];
    }
    double GetDinucSyn12_ATTT()
    {
        return dinucSyn_stat[0][3][15];
    }
    double GetDinucSyn12_CAAA()
    {
        return dinucSyn_stat[0][4][0];
    }
    double GetDinucSyn12_CAAC()
    {
        return dinucSyn_stat[0][4][1];
    }
    double GetDinucSyn12_CAAG()
    {
        return dinucSyn_stat[0][4][2];
    }
    double GetDinucSyn12_CAAT()
    {
        return dinucSyn_stat[0][4][3];
    }
    double GetDinucSyn12_CACA()
    {
        return dinucSyn_stat[0][4][4];
    }
    double GetDinucSyn12_CACC()
    {
        return dinucSyn_stat[0][4][5];
    }
    double GetDinucSyn12_CACG()
    {
        return dinucSyn_stat[0][4][6];
    }
    double GetDinucSyn12_CACT()
    {
        return dinucSyn_stat[0][4][7];
    }
    double GetDinucSyn12_CAGA()
    {
        return dinucSyn_stat[0][4][8];
    }
    double GetDinucSyn12_CAGC()
    {
        return dinucSyn_stat[0][4][9];
    }
    double GetDinucSyn12_CAGG()
    {
        return dinucSyn_stat[0][4][10];
    }
    double GetDinucSyn12_CAGT()
    {
        return dinucSyn_stat[0][4][11];
    }
    double GetDinucSyn12_CATA()
    {
        return dinucSyn_stat[0][4][12];
    }
    double GetDinucSyn12_CATC()
    {
        return dinucSyn_stat[0][4][13];
    }
    double GetDinucSyn12_CATG()
    {
        return dinucSyn_stat[0][4][14];
    }
    double GetDinucSyn12_CATT()
    {
        return dinucSyn_stat[0][4][15];
    }
    double GetDinucSyn12_CCAA()
    {
        return dinucSyn_stat[0][5][0];
    }
    double GetDinucSyn12_CCAC()
    {
        return dinucSyn_stat[0][5][1];
    }
    double GetDinucSyn12_CCAG()
    {
        return dinucSyn_stat[0][5][2];
    }
    double GetDinucSyn12_CCAT()
    {
        return dinucSyn_stat[0][5][3];
    }
    double GetDinucSyn12_CCCA()
    {
        return dinucSyn_stat[0][5][4];
    }
    double GetDinucSyn12_CCCC()
    {
        return dinucSyn_stat[0][5][5];
    }
    double GetDinucSyn12_CCCG()
    {
        return dinucSyn_stat[0][5][6];
    }
    double GetDinucSyn12_CCCT()
    {
        return dinucSyn_stat[0][5][7];
    }
    double GetDinucSyn12_CCGA()
    {
        return dinucSyn_stat[0][5][8];
    }
    double GetDinucSyn12_CCGC()
    {
        return dinucSyn_stat[0][5][9];
    }
    double GetDinucSyn12_CCGG()
    {
        return dinucSyn_stat[0][5][10];
    }
    double GetDinucSyn12_CCGT()
    {
        return dinucSyn_stat[0][5][11];
    }
    double GetDinucSyn12_CCTA()
    {
        return dinucSyn_stat[0][5][12];
    }
    double GetDinucSyn12_CCTC()
    {
        return dinucSyn_stat[0][5][13];
    }
    double GetDinucSyn12_CCTG()
    {
        return dinucSyn_stat[0][5][14];
    }
    double GetDinucSyn12_CCTT()
    {
        return dinucSyn_stat[0][5][15];
    }
    double GetDinucSyn12_CGAA()
    {
        return dinucSyn_stat[0][6][0];
    }
    double GetDinucSyn12_CGAC()
    {
        return dinucSyn_stat[0][6][1];
    }
    double GetDinucSyn12_CGAG()
    {
        return dinucSyn_stat[0][6][2];
    }
    double GetDinucSyn12_CGAT()
    {
        return dinucSyn_stat[0][6][3];
    }
    double GetDinucSyn12_CGCA()
    {
        return dinucSyn_stat[0][6][4];
    }
    double GetDinucSyn12_CGCC()
    {
        return dinucSyn_stat[0][6][5];
    }
    double GetDinucSyn12_CGCG()
    {
        return dinucSyn_stat[0][6][6];
    }
    double GetDinucSyn12_CGCT()
    {
        return dinucSyn_stat[0][6][7];
    }
    double GetDinucSyn12_CGGA()
    {
        return dinucSyn_stat[0][6][8];
    }
    double GetDinucSyn12_CGGC()
    {
        return dinucSyn_stat[0][6][9];
    }
    double GetDinucSyn12_CGGG()
    {
        return dinucSyn_stat[0][6][10];
    }
    double GetDinucSyn12_CGGT()
    {
        return dinucSyn_stat[0][6][11];
    }
    double GetDinucSyn12_CGTA()
    {
        return dinucSyn_stat[0][6][12];
    }
    double GetDinucSyn12_CGTC()
    {
        return dinucSyn_stat[0][6][13];
    }
    double GetDinucSyn12_CGTG()
    {
        return dinucSyn_stat[0][6][14];
    }
    double GetDinucSyn12_CGTT()
    {
        return dinucSyn_stat[0][6][15];
    }
    double GetDinucSyn12_CTAA()
    {
        return dinucSyn_stat[0][7][0];
    }
    double GetDinucSyn12_CTAC()
    {
        return dinucSyn_stat[0][7][1];
    }
    double GetDinucSyn12_CTAG()
    {
        return dinucSyn_stat[0][7][2];
    }
    double GetDinucSyn12_CTAT()
    {
        return dinucSyn_stat[0][7][3];
    }
    double GetDinucSyn12_CTCA()
    {
        return dinucSyn_stat[0][7][4];
    }
    double GetDinucSyn12_CTCC()
    {
        return dinucSyn_stat[0][7][5];
    }
    double GetDinucSyn12_CTCG()
    {
        return dinucSyn_stat[0][7][6];
    }
    double GetDinucSyn12_CTCT()
    {
        return dinucSyn_stat[0][7][7];
    }
    double GetDinucSyn12_CTGA()
    {
        return dinucSyn_stat[0][7][8];
    }
    double GetDinucSyn12_CTGC()
    {
        return dinucSyn_stat[0][7][9];
    }
    double GetDinucSyn12_CTGG()
    {
        return dinucSyn_stat[0][7][10];
    }
    double GetDinucSyn12_CTGT()
    {
        return dinucSyn_stat[0][7][11];
    }
    double GetDinucSyn12_CTTA()
    {
        return dinucSyn_stat[0][7][12];
    }
    double GetDinucSyn12_CTTC()
    {
        return dinucSyn_stat[0][7][13];
    }
    double GetDinucSyn12_CTTG()
    {
        return dinucSyn_stat[0][7][14];
    }
    double GetDinucSyn12_CTTT()
    {
        return dinucSyn_stat[0][7][15];
    }
    double GetDinucSyn12_GAAA()
    {
        return dinucSyn_stat[0][8][0];
    }
    double GetDinucSyn12_GAAC()
    {
        return dinucSyn_stat[0][8][1];
    }
    double GetDinucSyn12_GAAG()
    {
        return dinucSyn_stat[0][8][2];
    }
    double GetDinucSyn12_GAAT()
    {
        return dinucSyn_stat[0][8][3];
    }
    double GetDinucSyn12_GACA()
    {
        return dinucSyn_stat[0][8][4];
    }
    double GetDinucSyn12_GACC()
    {
        return dinucSyn_stat[0][8][5];
    }
    double GetDinucSyn12_GACG()
    {
        return dinucSyn_stat[0][8][6];
    }
    double GetDinucSyn12_GACT()
    {
        return dinucSyn_stat[0][8][7];
    }
    double GetDinucSyn12_GAGA()
    {
        return dinucSyn_stat[0][8][8];
    }
    double GetDinucSyn12_GAGC()
    {
        return dinucSyn_stat[0][8][9];
    }
    double GetDinucSyn12_GAGG()
    {
        return dinucSyn_stat[0][8][10];
    }
    double GetDinucSyn12_GAGT()
    {
        return dinucSyn_stat[0][8][11];
    }
    double GetDinucSyn12_GATA()
    {
        return dinucSyn_stat[0][8][12];
    }
    double GetDinucSyn12_GATC()
    {
        return dinucSyn_stat[0][8][13];
    }
    double GetDinucSyn12_GATG()
    {
        return dinucSyn_stat[0][8][14];
    }
    double GetDinucSyn12_GATT()
    {
        return dinucSyn_stat[0][8][15];
    }
    double GetDinucSyn12_GCAA()
    {
        return dinucSyn_stat[0][9][0];
    }
    double GetDinucSyn12_GCAC()
    {
        return dinucSyn_stat[0][9][1];
    }
    double GetDinucSyn12_GCAG()
    {
        return dinucSyn_stat[0][9][2];
    }
    double GetDinucSyn12_GCAT()
    {
        return dinucSyn_stat[0][9][3];
    }
    double GetDinucSyn12_GCCA()
    {
        return dinucSyn_stat[0][9][4];
    }
    double GetDinucSyn12_GCCC()
    {
        return dinucSyn_stat[0][9][5];
    }
    double GetDinucSyn12_GCCG()
    {
        return dinucSyn_stat[0][9][6];
    }
    double GetDinucSyn12_GCCT()
    {
        return dinucSyn_stat[0][9][7];
    }
    double GetDinucSyn12_GCGA()
    {
        return dinucSyn_stat[0][9][8];
    }
    double GetDinucSyn12_GCGC()
    {
        return dinucSyn_stat[0][9][9];
    }
    double GetDinucSyn12_GCGG()
    {
        return dinucSyn_stat[0][9][10];
    }
    double GetDinucSyn12_GCGT()
    {
        return dinucSyn_stat[0][9][11];
    }
    double GetDinucSyn12_GCTA()
    {
        return dinucSyn_stat[0][9][12];
    }
    double GetDinucSyn12_GCTC()
    {
        return dinucSyn_stat[0][9][13];
    }
    double GetDinucSyn12_GCTG()
    {
        return dinucSyn_stat[0][9][14];
    }
    double GetDinucSyn12_GCTT()
    {
        return dinucSyn_stat[0][9][15];
    }
    double GetDinucSyn12_GGAA()
    {
        return dinucSyn_stat[0][10][0];
    }
    double GetDinucSyn12_GGAC()
    {
        return dinucSyn_stat[0][10][1];
    }
    double GetDinucSyn12_GGAG()
    {
        return dinucSyn_stat[0][10][2];
    }
    double GetDinucSyn12_GGAT()
    {
        return dinucSyn_stat[0][10][3];
    }
    double GetDinucSyn12_GGCA()
    {
        return dinucSyn_stat[0][10][4];
    }
    double GetDinucSyn12_GGCC()
    {
        return dinucSyn_stat[0][10][5];
    }
    double GetDinucSyn12_GGCG()
    {
        return dinucSyn_stat[0][10][6];
    }
    double GetDinucSyn12_GGCT()
    {
        return dinucSyn_stat[0][10][7];
    }
    double GetDinucSyn12_GGGA()
    {
        return dinucSyn_stat[0][10][8];
    }
    double GetDinucSyn12_GGGC()
    {
        return dinucSyn_stat[0][10][9];
    }
    double GetDinucSyn12_GGGG()
    {
        return dinucSyn_stat[0][10][10];
    }
    double GetDinucSyn12_GGGT()
    {
        return dinucSyn_stat[0][10][11];
    }
    double GetDinucSyn12_GGTA()
    {
        return dinucSyn_stat[0][10][12];
    }
    double GetDinucSyn12_GGTC()
    {
        return dinucSyn_stat[0][10][13];
    }
    double GetDinucSyn12_GGTG()
    {
        return dinucSyn_stat[0][10][14];
    }
    double GetDinucSyn12_GGTT()
    {
        return dinucSyn_stat[0][10][15];
    }
    double GetDinucSyn12_GTAA()
    {
        return dinucSyn_stat[0][11][0];
    }
    double GetDinucSyn12_GTAC()
    {
        return dinucSyn_stat[0][11][1];
    }
    double GetDinucSyn12_GTAG()
    {
        return dinucSyn_stat[0][11][2];
    }
    double GetDinucSyn12_GTAT()
    {
        return dinucSyn_stat[0][11][3];
    }
    double GetDinucSyn12_GTCA()
    {
        return dinucSyn_stat[0][11][4];
    }
    double GetDinucSyn12_GTCC()
    {
        return dinucSyn_stat[0][11][5];
    }
    double GetDinucSyn12_GTCG()
    {
        return dinucSyn_stat[0][11][6];
    }
    double GetDinucSyn12_GTCT()
    {
        return dinucSyn_stat[0][11][7];
    }
    double GetDinucSyn12_GTGA()
    {
        return dinucSyn_stat[0][11][8];
    }
    double GetDinucSyn12_GTGC()
    {
        return dinucSyn_stat[0][11][9];
    }
    double GetDinucSyn12_GTGG()
    {
        return dinucSyn_stat[0][11][10];
    }
    double GetDinucSyn12_GTGT()
    {
        return dinucSyn_stat[0][11][11];
    }
    double GetDinucSyn12_GTTA()
    {
        return dinucSyn_stat[0][11][12];
    }
    double GetDinucSyn12_GTTC()
    {
        return dinucSyn_stat[0][11][13];
    }
    double GetDinucSyn12_GTTG()
    {
        return dinucSyn_stat[0][11][14];
    }
    double GetDinucSyn12_GTTT()
    {
        return dinucSyn_stat[0][11][15];
    }
    double GetDinucSyn12_TAAA()
    {
        return dinucSyn_stat[0][12][0];
    }
    double GetDinucSyn12_TAAC()
    {
        return dinucSyn_stat[0][12][1];
    }
    double GetDinucSyn12_TAAG()
    {
        return dinucSyn_stat[0][12][2];
    }
    double GetDinucSyn12_TAAT()
    {
        return dinucSyn_stat[0][12][3];
    }
    double GetDinucSyn12_TACA()
    {
        return dinucSyn_stat[0][12][4];
    }
    double GetDinucSyn12_TACC()
    {
        return dinucSyn_stat[0][12][5];
    }
    double GetDinucSyn12_TACG()
    {
        return dinucSyn_stat[0][12][6];
    }
    double GetDinucSyn12_TACT()
    {
        return dinucSyn_stat[0][12][7];
    }
    double GetDinucSyn12_TAGA()
    {
        return dinucSyn_stat[0][12][8];
    }
    double GetDinucSyn12_TAGC()
    {
        return dinucSyn_stat[0][12][9];
    }
    double GetDinucSyn12_TAGG()
    {
        return dinucSyn_stat[0][12][10];
    }
    double GetDinucSyn12_TAGT()
    {
        return dinucSyn_stat[0][12][11];
    }
    double GetDinucSyn12_TATA()
    {
        return dinucSyn_stat[0][12][12];
    }
    double GetDinucSyn12_TATC()
    {
        return dinucSyn_stat[0][12][13];
    }
    double GetDinucSyn12_TATG()
    {
        return dinucSyn_stat[0][12][14];
    }
    double GetDinucSyn12_TATT()
    {
        return dinucSyn_stat[0][12][15];
    }
    double GetDinucSyn12_TCAA()
    {
        return dinucSyn_stat[0][13][0];
    }
    double GetDinucSyn12_TCAC()
    {
        return dinucSyn_stat[0][13][1];
    }
    double GetDinucSyn12_TCAG()
    {
        return dinucSyn_stat[0][13][2];
    }
    double GetDinucSyn12_TCAT()
    {
        return dinucSyn_stat[0][13][3];
    }
    double GetDinucSyn12_TCCA()
    {
        return dinucSyn_stat[0][13][4];
    }
    double GetDinucSyn12_TCCC()
    {
        return dinucSyn_stat[0][13][5];
    }
    double GetDinucSyn12_TCCG()
    {
        return dinucSyn_stat[0][13][6];
    }
    double GetDinucSyn12_TCCT()
    {
        return dinucSyn_stat[0][13][7];
    }
    double GetDinucSyn12_TCGA()
    {
        return dinucSyn_stat[0][13][8];
    }
    double GetDinucSyn12_TCGC()
    {
        return dinucSyn_stat[0][13][9];
    }
    double GetDinucSyn12_TCGG()
    {
        return dinucSyn_stat[0][13][10];
    }
    double GetDinucSyn12_TCGT()
    {
        return dinucSyn_stat[0][13][11];
    }
    double GetDinucSyn12_TCTA()
    {
        return dinucSyn_stat[0][13][12];
    }
    double GetDinucSyn12_TCTC()
    {
        return dinucSyn_stat[0][13][13];
    }
    double GetDinucSyn12_TCTG()
    {
        return dinucSyn_stat[0][13][14];
    }
    double GetDinucSyn12_TCTT()
    {
        return dinucSyn_stat[0][13][15];
    }
    double GetDinucSyn12_TGAA()
    {
        return dinucSyn_stat[0][14][0];
    }
    double GetDinucSyn12_TGAC()
    {
        return dinucSyn_stat[0][14][1];
    }
    double GetDinucSyn12_TGAG()
    {
        return dinucSyn_stat[0][14][2];
    }
    double GetDinucSyn12_TGAT()
    {
        return dinucSyn_stat[0][14][3];
    }
    double GetDinucSyn12_TGCA()
    {
        return dinucSyn_stat[0][14][4];
    }
    double GetDinucSyn12_TGCC()
    {
        return dinucSyn_stat[0][14][5];
    }
    double GetDinucSyn12_TGCG()
    {
        return dinucSyn_stat[0][14][6];
    }
    double GetDinucSyn12_TGCT()
    {
        return dinucSyn_stat[0][14][7];
    }
    double GetDinucSyn12_TGGA()
    {
        return dinucSyn_stat[0][14][8];
    }
    double GetDinucSyn12_TGGC()
    {
        return dinucSyn_stat[0][14][9];
    }
    double GetDinucSyn12_TGGG()
    {
        return dinucSyn_stat[0][14][10];
    }
    double GetDinucSyn12_TGGT()
    {
        return dinucSyn_stat[0][14][11];
    }
    double GetDinucSyn12_TGTA()
    {
        return dinucSyn_stat[0][14][12];
    }
    double GetDinucSyn12_TGTC()
    {
        return dinucSyn_stat[0][14][13];
    }
    double GetDinucSyn12_TGTG()
    {
        return dinucSyn_stat[0][14][14];
    }
    double GetDinucSyn12_TGTT()
    {
        return dinucSyn_stat[0][14][15];
    }
    double GetDinucSyn12_TTAA()
    {
        return dinucSyn_stat[0][15][0];
    }
    double GetDinucSyn12_TTAC()
    {
        return dinucSyn_stat[0][15][1];
    }
    double GetDinucSyn12_TTAG()
    {
        return dinucSyn_stat[0][15][2];
    }
    double GetDinucSyn12_TTAT()
    {
        return dinucSyn_stat[0][15][3];
    }
    double GetDinucSyn12_TTCA()
    {
        return dinucSyn_stat[0][15][4];
    }
    double GetDinucSyn12_TTCC()
    {
        return dinucSyn_stat[0][15][5];
    }
    double GetDinucSyn12_TTCG()
    {
        return dinucSyn_stat[0][15][6];
    }
    double GetDinucSyn12_TTCT()
    {
        return dinucSyn_stat[0][15][7];
    }
    double GetDinucSyn12_TTGA()
    {
        return dinucSyn_stat[0][15][8];
    }
    double GetDinucSyn12_TTGC()
    {
        return dinucSyn_stat[0][15][9];
    }
    double GetDinucSyn12_TTGG()
    {
        return dinucSyn_stat[0][15][10];
    }
    double GetDinucSyn12_TTGT()
    {
        return dinucSyn_stat[0][15][11];
    }
    double GetDinucSyn12_TTTA()
    {
        return dinucSyn_stat[0][15][12];
    }
    double GetDinucSyn12_TTTC()
    {
        return dinucSyn_stat[0][15][13];
    }
    double GetDinucSyn12_TTTG()
    {
        return dinucSyn_stat[0][15][14];
    }
    double GetDinucSyn12_TTTT()
    {
        return dinucSyn_stat[0][15][15];
    }

    double GetDinucSyn23_AAAA()
    {
        return dinucSyn_stat[1][0][0];
    }
    double GetDinucSyn23_AAAC()
    {
        return dinucSyn_stat[1][0][1];
    }
    double GetDinucSyn23_AAAG()
    {
        return dinucSyn_stat[1][0][2];
    }
    double GetDinucSyn23_AAAT()
    {
        return dinucSyn_stat[1][0][3];
    }
    double GetDinucSyn23_AACA()
    {
        return dinucSyn_stat[1][0][4];
    }
    double GetDinucSyn23_AACC()
    {
        return dinucSyn_stat[1][0][5];
    }
    double GetDinucSyn23_AACG()
    {
        return dinucSyn_stat[1][0][6];
    }
    double GetDinucSyn23_AACT()
    {
        return dinucSyn_stat[1][0][7];
    }
    double GetDinucSyn23_AAGA()
    {
        return dinucSyn_stat[1][0][8];
    }
    double GetDinucSyn23_AAGC()
    {
        return dinucSyn_stat[1][0][9];
    }
    double GetDinucSyn23_AAGG()
    {
        return dinucSyn_stat[1][0][10];
    }
    double GetDinucSyn23_AAGT()
    {
        return dinucSyn_stat[1][0][11];
    }
    double GetDinucSyn23_AATA()
    {
        return dinucSyn_stat[1][0][12];
    }
    double GetDinucSyn23_AATC()
    {
        return dinucSyn_stat[1][0][13];
    }
    double GetDinucSyn23_AATG()
    {
        return dinucSyn_stat[1][0][14];
    }
    double GetDinucSyn23_AATT()
    {
        return dinucSyn_stat[1][0][15];
    }
    double GetDinucSyn23_ACAA()
    {
        return dinucSyn_stat[1][1][0];
    }
    double GetDinucSyn23_ACAC()
    {
        return dinucSyn_stat[1][1][1];
    }
    double GetDinucSyn23_ACAG()
    {
        return dinucSyn_stat[1][1][2];
    }
    double GetDinucSyn23_ACAT()
    {
        return dinucSyn_stat[1][1][3];
    }
    double GetDinucSyn23_ACCA()
    {
        return dinucSyn_stat[1][1][4];
    }
    double GetDinucSyn23_ACCC()
    {
        return dinucSyn_stat[1][1][5];
    }
    double GetDinucSyn23_ACCG()
    {
        return dinucSyn_stat[1][1][6];
    }
    double GetDinucSyn23_ACCT()
    {
        return dinucSyn_stat[1][1][7];
    }
    double GetDinucSyn23_ACGA()
    {
        return dinucSyn_stat[1][1][8];
    }
    double GetDinucSyn23_ACGC()
    {
        return dinucSyn_stat[1][1][9];
    }
    double GetDinucSyn23_ACGG()
    {
        return dinucSyn_stat[1][1][10];
    }
    double GetDinucSyn23_ACGT()
    {
        return dinucSyn_stat[1][1][11];
    }
    double GetDinucSyn23_ACTA()
    {
        return dinucSyn_stat[1][1][12];
    }
    double GetDinucSyn23_ACTC()
    {
        return dinucSyn_stat[1][1][13];
    }
    double GetDinucSyn23_ACTG()
    {
        return dinucSyn_stat[1][1][14];
    }
    double GetDinucSyn23_ACTT()
    {
        return dinucSyn_stat[1][1][15];
    }
    double GetDinucSyn23_AGAA()
    {
        return dinucSyn_stat[1][2][0];
    }
    double GetDinucSyn23_AGAC()
    {
        return dinucSyn_stat[1][2][1];
    }
    double GetDinucSyn23_AGAG()
    {
        return dinucSyn_stat[1][2][2];
    }
    double GetDinucSyn23_AGAT()
    {
        return dinucSyn_stat[1][2][3];
    }
    double GetDinucSyn23_AGCA()
    {
        return dinucSyn_stat[1][2][4];
    }
    double GetDinucSyn23_AGCC()
    {
        return dinucSyn_stat[1][2][5];
    }
    double GetDinucSyn23_AGCG()
    {
        return dinucSyn_stat[1][2][6];
    }
    double GetDinucSyn23_AGCT()
    {
        return dinucSyn_stat[1][2][7];
    }
    double GetDinucSyn23_AGGA()
    {
        return dinucSyn_stat[1][2][8];
    }
    double GetDinucSyn23_AGGC()
    {
        return dinucSyn_stat[1][2][9];
    }
    double GetDinucSyn23_AGGG()
    {
        return dinucSyn_stat[1][2][10];
    }
    double GetDinucSyn23_AGGT()
    {
        return dinucSyn_stat[1][2][11];
    }
    double GetDinucSyn23_AGTA()
    {
        return dinucSyn_stat[1][2][12];
    }
    double GetDinucSyn23_AGTC()
    {
        return dinucSyn_stat[1][2][13];
    }
    double GetDinucSyn23_AGTG()
    {
        return dinucSyn_stat[1][2][14];
    }
    double GetDinucSyn23_AGTT()
    {
        return dinucSyn_stat[1][2][15];
    }
    double GetDinucSyn23_ATAA()
    {
        return dinucSyn_stat[1][3][0];
    }
    double GetDinucSyn23_ATAC()
    {
        return dinucSyn_stat[1][3][1];
    }
    double GetDinucSyn23_ATAG()
    {
        return dinucSyn_stat[1][3][2];
    }
    double GetDinucSyn23_ATAT()
    {
        return dinucSyn_stat[1][3][3];
    }
    double GetDinucSyn23_ATCA()
    {
        return dinucSyn_stat[1][3][4];
    }
    double GetDinucSyn23_ATCC()
    {
        return dinucSyn_stat[1][3][5];
    }
    double GetDinucSyn23_ATCG()
    {
        return dinucSyn_stat[1][3][6];
    }
    double GetDinucSyn23_ATCT()
    {
        return dinucSyn_stat[1][3][7];
    }
    double GetDinucSyn23_ATGA()
    {
        return dinucSyn_stat[1][3][8];
    }
    double GetDinucSyn23_ATGC()
    {
        return dinucSyn_stat[1][3][9];
    }
    double GetDinucSyn23_ATGG()
    {
        return dinucSyn_stat[1][3][10];
    }
    double GetDinucSyn23_ATGT()
    {
        return dinucSyn_stat[1][3][11];
    }
    double GetDinucSyn23_ATTA()
    {
        return dinucSyn_stat[1][3][12];
    }
    double GetDinucSyn23_ATTC()
    {
        return dinucSyn_stat[1][3][13];
    }
    double GetDinucSyn23_ATTG()
    {
        return dinucSyn_stat[1][3][14];
    }
    double GetDinucSyn23_ATTT()
    {
        return dinucSyn_stat[1][3][15];
    }
    double GetDinucSyn23_CAAA()
    {
        return dinucSyn_stat[1][4][0];
    }
    double GetDinucSyn23_CAAC()
    {
        return dinucSyn_stat[1][4][1];
    }
    double GetDinucSyn23_CAAG()
    {
        return dinucSyn_stat[1][4][2];
    }
    double GetDinucSyn23_CAAT()
    {
        return dinucSyn_stat[1][4][3];
    }
    double GetDinucSyn23_CACA()
    {
        return dinucSyn_stat[1][4][4];
    }
    double GetDinucSyn23_CACC()
    {
        return dinucSyn_stat[1][4][5];
    }
    double GetDinucSyn23_CACG()
    {
        return dinucSyn_stat[1][4][6];
    }
    double GetDinucSyn23_CACT()
    {
        return dinucSyn_stat[1][4][7];
    }
    double GetDinucSyn23_CAGA()
    {
        return dinucSyn_stat[1][4][8];
    }
    double GetDinucSyn23_CAGC()
    {
        return dinucSyn_stat[1][4][9];
    }
    double GetDinucSyn23_CAGG()
    {
        return dinucSyn_stat[1][4][10];
    }
    double GetDinucSyn23_CAGT()
    {
        return dinucSyn_stat[1][4][11];
    }
    double GetDinucSyn23_CATA()
    {
        return dinucSyn_stat[1][4][12];
    }
    double GetDinucSyn23_CATC()
    {
        return dinucSyn_stat[1][4][13];
    }
    double GetDinucSyn23_CATG()
    {
        return dinucSyn_stat[1][4][14];
    }
    double GetDinucSyn23_CATT()
    {
        return dinucSyn_stat[1][4][15];
    }
    double GetDinucSyn23_CCAA()
    {
        return dinucSyn_stat[1][5][0];
    }
    double GetDinucSyn23_CCAC()
    {
        return dinucSyn_stat[1][5][1];
    }
    double GetDinucSyn23_CCAG()
    {
        return dinucSyn_stat[1][5][2];
    }
    double GetDinucSyn23_CCAT()
    {
        return dinucSyn_stat[1][5][3];
    }
    double GetDinucSyn23_CCCA()
    {
        return dinucSyn_stat[1][5][4];
    }
    double GetDinucSyn23_CCCC()
    {
        return dinucSyn_stat[1][5][5];
    }
    double GetDinucSyn23_CCCG()
    {
        return dinucSyn_stat[1][5][6];
    }
    double GetDinucSyn23_CCCT()
    {
        return dinucSyn_stat[1][5][7];
    }
    double GetDinucSyn23_CCGA()
    {
        return dinucSyn_stat[1][5][8];
    }
    double GetDinucSyn23_CCGC()
    {
        return dinucSyn_stat[1][5][9];
    }
    double GetDinucSyn23_CCGG()
    {
        return dinucSyn_stat[1][5][10];
    }
    double GetDinucSyn23_CCGT()
    {
        return dinucSyn_stat[1][5][11];
    }
    double GetDinucSyn23_CCTA()
    {
        return dinucSyn_stat[1][5][12];
    }
    double GetDinucSyn23_CCTC()
    {
        return dinucSyn_stat[1][5][13];
    }
    double GetDinucSyn23_CCTG()
    {
        return dinucSyn_stat[1][5][14];
    }
    double GetDinucSyn23_CCTT()
    {
        return dinucSyn_stat[1][5][15];
    }
    double GetDinucSyn23_CGAA()
    {
        return dinucSyn_stat[1][6][0];
    }
    double GetDinucSyn23_CGAC()
    {
        return dinucSyn_stat[1][6][1];
    }
    double GetDinucSyn23_CGAG()
    {
        return dinucSyn_stat[1][6][2];
    }
    double GetDinucSyn23_CGAT()
    {
        return dinucSyn_stat[1][6][3];
    }
    double GetDinucSyn23_CGCA()
    {
        return dinucSyn_stat[1][6][4];
    }
    double GetDinucSyn23_CGCC()
    {
        return dinucSyn_stat[1][6][5];
    }
    double GetDinucSyn23_CGCG()
    {
        return dinucSyn_stat[1][6][6];
    }
    double GetDinucSyn23_CGCT()
    {
        return dinucSyn_stat[1][6][7];
    }
    double GetDinucSyn23_CGGA()
    {
        return dinucSyn_stat[1][6][8];
    }
    double GetDinucSyn23_CGGC()
    {
        return dinucSyn_stat[1][6][9];
    }
    double GetDinucSyn23_CGGG()
    {
        return dinucSyn_stat[1][6][10];
    }
    double GetDinucSyn23_CGGT()
    {
        return dinucSyn_stat[1][6][11];
    }
    double GetDinucSyn23_CGTA()
    {
        return dinucSyn_stat[1][6][12];
    }
    double GetDinucSyn23_CGTC()
    {
        return dinucSyn_stat[1][6][13];
    }
    double GetDinucSyn23_CGTG()
    {
        return dinucSyn_stat[1][6][14];
    }
    double GetDinucSyn23_CGTT()
    {
        return dinucSyn_stat[1][6][15];
    }
    double GetDinucSyn23_CTAA()
    {
        return dinucSyn_stat[1][7][0];
    }
    double GetDinucSyn23_CTAC()
    {
        return dinucSyn_stat[1][7][1];
    }
    double GetDinucSyn23_CTAG()
    {
        return dinucSyn_stat[1][7][2];
    }
    double GetDinucSyn23_CTAT()
    {
        return dinucSyn_stat[1][7][3];
    }
    double GetDinucSyn23_CTCA()
    {
        return dinucSyn_stat[1][7][4];
    }
    double GetDinucSyn23_CTCC()
    {
        return dinucSyn_stat[1][7][5];
    }
    double GetDinucSyn23_CTCG()
    {
        return dinucSyn_stat[1][7][6];
    }
    double GetDinucSyn23_CTCT()
    {
        return dinucSyn_stat[1][7][7];
    }
    double GetDinucSyn23_CTGA()
    {
        return dinucSyn_stat[1][7][8];
    }
    double GetDinucSyn23_CTGC()
    {
        return dinucSyn_stat[1][7][9];
    }
    double GetDinucSyn23_CTGG()
    {
        return dinucSyn_stat[1][7][10];
    }
    double GetDinucSyn23_CTGT()
    {
        return dinucSyn_stat[1][7][11];
    }
    double GetDinucSyn23_CTTA()
    {
        return dinucSyn_stat[1][7][12];
    }
    double GetDinucSyn23_CTTC()
    {
        return dinucSyn_stat[1][7][13];
    }
    double GetDinucSyn23_CTTG()
    {
        return dinucSyn_stat[1][7][14];
    }
    double GetDinucSyn23_CTTT()
    {
        return dinucSyn_stat[1][7][15];
    }
    double GetDinucSyn23_GAAA()
    {
        return dinucSyn_stat[1][8][0];
    }
    double GetDinucSyn23_GAAC()
    {
        return dinucSyn_stat[1][8][1];
    }
    double GetDinucSyn23_GAAG()
    {
        return dinucSyn_stat[1][8][2];
    }
    double GetDinucSyn23_GAAT()
    {
        return dinucSyn_stat[1][8][3];
    }
    double GetDinucSyn23_GACA()
    {
        return dinucSyn_stat[1][8][4];
    }
    double GetDinucSyn23_GACC()
    {
        return dinucSyn_stat[1][8][5];
    }
    double GetDinucSyn23_GACG()
    {
        return dinucSyn_stat[1][8][6];
    }
    double GetDinucSyn23_GACT()
    {
        return dinucSyn_stat[1][8][7];
    }
    double GetDinucSyn23_GAGA()
    {
        return dinucSyn_stat[1][8][8];
    }
    double GetDinucSyn23_GAGC()
    {
        return dinucSyn_stat[1][8][9];
    }
    double GetDinucSyn23_GAGG()
    {
        return dinucSyn_stat[1][8][10];
    }
    double GetDinucSyn23_GAGT()
    {
        return dinucSyn_stat[1][8][11];
    }
    double GetDinucSyn23_GATA()
    {
        return dinucSyn_stat[1][8][12];
    }
    double GetDinucSyn23_GATC()
    {
        return dinucSyn_stat[1][8][13];
    }
    double GetDinucSyn23_GATG()
    {
        return dinucSyn_stat[1][8][14];
    }
    double GetDinucSyn23_GATT()
    {
        return dinucSyn_stat[1][8][15];
    }
    double GetDinucSyn23_GCAA()
    {
        return dinucSyn_stat[1][9][0];
    }
    double GetDinucSyn23_GCAC()
    {
        return dinucSyn_stat[1][9][1];
    }
    double GetDinucSyn23_GCAG()
    {
        return dinucSyn_stat[1][9][2];
    }
    double GetDinucSyn23_GCAT()
    {
        return dinucSyn_stat[1][9][3];
    }
    double GetDinucSyn23_GCCA()
    {
        return dinucSyn_stat[1][9][4];
    }
    double GetDinucSyn23_GCCC()
    {
        return dinucSyn_stat[1][9][5];
    }
    double GetDinucSyn23_GCCG()
    {
        return dinucSyn_stat[1][9][6];
    }
    double GetDinucSyn23_GCCT()
    {
        return dinucSyn_stat[1][9][7];
    }
    double GetDinucSyn23_GCGA()
    {
        return dinucSyn_stat[1][9][8];
    }
    double GetDinucSyn23_GCGC()
    {
        return dinucSyn_stat[1][9][9];
    }
    double GetDinucSyn23_GCGG()
    {
        return dinucSyn_stat[1][9][10];
    }
    double GetDinucSyn23_GCGT()
    {
        return dinucSyn_stat[1][9][11];
    }
    double GetDinucSyn23_GCTA()
    {
        return dinucSyn_stat[1][9][12];
    }
    double GetDinucSyn23_GCTC()
    {
        return dinucSyn_stat[1][9][13];
    }
    double GetDinucSyn23_GCTG()
    {
        return dinucSyn_stat[1][9][14];
    }
    double GetDinucSyn23_GCTT()
    {
        return dinucSyn_stat[1][9][15];
    }
    double GetDinucSyn23_GGAA()
    {
        return dinucSyn_stat[1][10][0];
    }
    double GetDinucSyn23_GGAC()
    {
        return dinucSyn_stat[1][10][1];
    }
    double GetDinucSyn23_GGAG()
    {
        return dinucSyn_stat[1][10][2];
    }
    double GetDinucSyn23_GGAT()
    {
        return dinucSyn_stat[1][10][3];
    }
    double GetDinucSyn23_GGCA()
    {
        return dinucSyn_stat[1][10][4];
    }
    double GetDinucSyn23_GGCC()
    {
        return dinucSyn_stat[1][10][5];
    }
    double GetDinucSyn23_GGCG()
    {
        return dinucSyn_stat[1][10][6];
    }
    double GetDinucSyn23_GGCT()
    {
        return dinucSyn_stat[1][10][7];
    }
    double GetDinucSyn23_GGGA()
    {
        return dinucSyn_stat[1][10][8];
    }
    double GetDinucSyn23_GGGC()
    {
        return dinucSyn_stat[1][10][9];
    }
    double GetDinucSyn23_GGGG()
    {
        return dinucSyn_stat[1][10][10];
    }
    double GetDinucSyn23_GGGT()
    {
        return dinucSyn_stat[1][10][11];
    }
    double GetDinucSyn23_GGTA()
    {
        return dinucSyn_stat[1][10][12];
    }
    double GetDinucSyn23_GGTC()
    {
        return dinucSyn_stat[1][10][13];
    }
    double GetDinucSyn23_GGTG()
    {
        return dinucSyn_stat[1][10][14];
    }
    double GetDinucSyn23_GGTT()
    {
        return dinucSyn_stat[1][10][15];
    }
    double GetDinucSyn23_GTAA()
    {
        return dinucSyn_stat[1][11][0];
    }
    double GetDinucSyn23_GTAC()
    {
        return dinucSyn_stat[1][11][1];
    }
    double GetDinucSyn23_GTAG()
    {
        return dinucSyn_stat[1][11][2];
    }
    double GetDinucSyn23_GTAT()
    {
        return dinucSyn_stat[1][11][3];
    }
    double GetDinucSyn23_GTCA()
    {
        return dinucSyn_stat[1][11][4];
    }
    double GetDinucSyn23_GTCC()
    {
        return dinucSyn_stat[1][11][5];
    }
    double GetDinucSyn23_GTCG()
    {
        return dinucSyn_stat[1][11][6];
    }
    double GetDinucSyn23_GTCT()
    {
        return dinucSyn_stat[1][11][7];
    }
    double GetDinucSyn23_GTGA()
    {
        return dinucSyn_stat[1][11][8];
    }
    double GetDinucSyn23_GTGC()
    {
        return dinucSyn_stat[1][11][9];
    }
    double GetDinucSyn23_GTGG()
    {
        return dinucSyn_stat[1][11][10];
    }
    double GetDinucSyn23_GTGT()
    {
        return dinucSyn_stat[1][11][11];
    }
    double GetDinucSyn23_GTTA()
    {
        return dinucSyn_stat[1][11][12];
    }
    double GetDinucSyn23_GTTC()
    {
        return dinucSyn_stat[1][11][13];
    }
    double GetDinucSyn23_GTTG()
    {
        return dinucSyn_stat[1][11][14];
    }
    double GetDinucSyn23_GTTT()
    {
        return dinucSyn_stat[1][11][15];
    }
    double GetDinucSyn23_TAAA()
    {
        return dinucSyn_stat[1][12][0];
    }
    double GetDinucSyn23_TAAC()
    {
        return dinucSyn_stat[1][12][1];
    }
    double GetDinucSyn23_TAAG()
    {
        return dinucSyn_stat[1][12][2];
    }
    double GetDinucSyn23_TAAT()
    {
        return dinucSyn_stat[1][12][3];
    }
    double GetDinucSyn23_TACA()
    {
        return dinucSyn_stat[1][12][4];
    }
    double GetDinucSyn23_TACC()
    {
        return dinucSyn_stat[1][12][5];
    }
    double GetDinucSyn23_TACG()
    {
        return dinucSyn_stat[1][12][6];
    }
    double GetDinucSyn23_TACT()
    {
        return dinucSyn_stat[1][12][7];
    }
    double GetDinucSyn23_TAGA()
    {
        return dinucSyn_stat[1][12][8];
    }
    double GetDinucSyn23_TAGC()
    {
        return dinucSyn_stat[1][12][9];
    }
    double GetDinucSyn23_TAGG()
    {
        return dinucSyn_stat[1][12][10];
    }
    double GetDinucSyn23_TAGT()
    {
        return dinucSyn_stat[1][12][11];
    }
    double GetDinucSyn23_TATA()
    {
        return dinucSyn_stat[1][12][12];
    }
    double GetDinucSyn23_TATC()
    {
        return dinucSyn_stat[1][12][13];
    }
    double GetDinucSyn23_TATG()
    {
        return dinucSyn_stat[1][12][14];
    }
    double GetDinucSyn23_TATT()
    {
        return dinucSyn_stat[1][12][15];
    }
    double GetDinucSyn23_TCAA()
    {
        return dinucSyn_stat[1][13][0];
    }
    double GetDinucSyn23_TCAC()
    {
        return dinucSyn_stat[1][13][1];
    }
    double GetDinucSyn23_TCAG()
    {
        return dinucSyn_stat[1][13][2];
    }
    double GetDinucSyn23_TCAT()
    {
        return dinucSyn_stat[1][13][3];
    }
    double GetDinucSyn23_TCCA()
    {
        return dinucSyn_stat[1][13][4];
    }
    double GetDinucSyn23_TCCC()
    {
        return dinucSyn_stat[1][13][5];
    }
    double GetDinucSyn23_TCCG()
    {
        return dinucSyn_stat[1][13][6];
    }
    double GetDinucSyn23_TCCT()
    {
        return dinucSyn_stat[1][13][7];
    }
    double GetDinucSyn23_TCGA()
    {
        return dinucSyn_stat[1][13][8];
    }
    double GetDinucSyn23_TCGC()
    {
        return dinucSyn_stat[1][13][9];
    }
    double GetDinucSyn23_TCGG()
    {
        return dinucSyn_stat[1][13][10];
    }
    double GetDinucSyn23_TCGT()
    {
        return dinucSyn_stat[1][13][11];
    }
    double GetDinucSyn23_TCTA()
    {
        return dinucSyn_stat[1][13][12];
    }
    double GetDinucSyn23_TCTC()
    {
        return dinucSyn_stat[1][13][13];
    }
    double GetDinucSyn23_TCTG()
    {
        return dinucSyn_stat[1][13][14];
    }
    double GetDinucSyn23_TCTT()
    {
        return dinucSyn_stat[1][13][15];
    }
    double GetDinucSyn23_TGAA()
    {
        return dinucSyn_stat[1][14][0];
    }
    double GetDinucSyn23_TGAC()
    {
        return dinucSyn_stat[1][14][1];
    }
    double GetDinucSyn23_TGAG()
    {
        return dinucSyn_stat[1][14][2];
    }
    double GetDinucSyn23_TGAT()
    {
        return dinucSyn_stat[1][14][3];
    }
    double GetDinucSyn23_TGCA()
    {
        return dinucSyn_stat[1][14][4];
    }
    double GetDinucSyn23_TGCC()
    {
        return dinucSyn_stat[1][14][5];
    }
    double GetDinucSyn23_TGCG()
    {
        return dinucSyn_stat[1][14][6];
    }
    double GetDinucSyn23_TGCT()
    {
        return dinucSyn_stat[1][14][7];
    }
    double GetDinucSyn23_TGGA()
    {
        return dinucSyn_stat[1][14][8];
    }
    double GetDinucSyn23_TGGC()
    {
        return dinucSyn_stat[1][14][9];
    }
    double GetDinucSyn23_TGGG()
    {
        return dinucSyn_stat[1][14][10];
    }
    double GetDinucSyn23_TGGT()
    {
        return dinucSyn_stat[1][14][11];
    }
    double GetDinucSyn23_TGTA()
    {
        return dinucSyn_stat[1][14][12];
    }
    double GetDinucSyn23_TGTC()
    {
        return dinucSyn_stat[1][14][13];
    }
    double GetDinucSyn23_TGTG()
    {
        return dinucSyn_stat[1][14][14];
    }
    double GetDinucSyn23_TGTT()
    {
        return dinucSyn_stat[1][14][15];
    }
    double GetDinucSyn23_TTAA()
    {
        return dinucSyn_stat[1][15][0];
    }
    double GetDinucSyn23_TTAC()
    {
        return dinucSyn_stat[1][15][1];
    }
    double GetDinucSyn23_TTAG()
    {
        return dinucSyn_stat[1][15][2];
    }
    double GetDinucSyn23_TTAT()
    {
        return dinucSyn_stat[1][15][3];
    }
    double GetDinucSyn23_TTCA()
    {
        return dinucSyn_stat[1][15][4];
    }
    double GetDinucSyn23_TTCC()
    {
        return dinucSyn_stat[1][15][5];
    }
    double GetDinucSyn23_TTCG()
    {
        return dinucSyn_stat[1][15][6];
    }
    double GetDinucSyn23_TTCT()
    {
        return dinucSyn_stat[1][15][7];
    }
    double GetDinucSyn23_TTGA()
    {
        return dinucSyn_stat[1][15][8];
    }
    double GetDinucSyn23_TTGC()
    {
        return dinucSyn_stat[1][15][9];
    }
    double GetDinucSyn23_TTGG()
    {
        return dinucSyn_stat[1][15][10];
    }
    double GetDinucSyn23_TTGT()
    {
        return dinucSyn_stat[1][15][11];
    }
    double GetDinucSyn23_TTTA()
    {
        return dinucSyn_stat[1][15][12];
    }
    double GetDinucSyn23_TTTC()
    {
        return dinucSyn_stat[1][15][13];
    }
    double GetDinucSyn23_TTTG()
    {
        return dinucSyn_stat[1][15][14];
    }
    double GetDinucSyn23_TTTT()
    {
        return dinucSyn_stat[1][15][15];
    }

    double GetDinucSyn31_AAAA()
    {
        return dinucSyn_stat[2][0][0];
    }
    double GetDinucSyn31_AAAC()
    {
        return dinucSyn_stat[2][0][1];
    }
    double GetDinucSyn31_AAAG()
    {
        return dinucSyn_stat[2][0][2];
    }
    double GetDinucSyn31_AAAT()
    {
        return dinucSyn_stat[2][0][3];
    }
    double GetDinucSyn31_AACA()
    {
        return dinucSyn_stat[2][0][4];
    }
    double GetDinucSyn31_AACC()
    {
        return dinucSyn_stat[2][0][5];
    }
    double GetDinucSyn31_AACG()
    {
        return dinucSyn_stat[2][0][6];
    }
    double GetDinucSyn31_AACT()
    {
        return dinucSyn_stat[2][0][7];
    }
    double GetDinucSyn31_AAGA()
    {
        return dinucSyn_stat[2][0][8];
    }
    double GetDinucSyn31_AAGC()
    {
        return dinucSyn_stat[2][0][9];
    }
    double GetDinucSyn31_AAGG()
    {
        return dinucSyn_stat[2][0][10];
    }
    double GetDinucSyn31_AAGT()
    {
        return dinucSyn_stat[2][0][11];
    }
    double GetDinucSyn31_AATA()
    {
        return dinucSyn_stat[2][0][12];
    }
    double GetDinucSyn31_AATC()
    {
        return dinucSyn_stat[2][0][13];
    }
    double GetDinucSyn31_AATG()
    {
        return dinucSyn_stat[2][0][14];
    }
    double GetDinucSyn31_AATT()
    {
        return dinucSyn_stat[2][0][15];
    }
    double GetDinucSyn31_ACAA()
    {
        return dinucSyn_stat[2][1][0];
    }
    double GetDinucSyn31_ACAC()
    {
        return dinucSyn_stat[2][1][1];
    }
    double GetDinucSyn31_ACAG()
    {
        return dinucSyn_stat[2][1][2];
    }
    double GetDinucSyn31_ACAT()
    {
        return dinucSyn_stat[2][1][3];
    }
    double GetDinucSyn31_ACCA()
    {
        return dinucSyn_stat[2][1][4];
    }
    double GetDinucSyn31_ACCC()
    {
        return dinucSyn_stat[2][1][5];
    }
    double GetDinucSyn31_ACCG()
    {
        return dinucSyn_stat[2][1][6];
    }
    double GetDinucSyn31_ACCT()
    {
        return dinucSyn_stat[2][1][7];
    }
    double GetDinucSyn31_ACGA()
    {
        return dinucSyn_stat[2][1][8];
    }
    double GetDinucSyn31_ACGC()
    {
        return dinucSyn_stat[2][1][9];
    }
    double GetDinucSyn31_ACGG()
    {
        return dinucSyn_stat[2][1][10];
    }
    double GetDinucSyn31_ACGT()
    {
        return dinucSyn_stat[2][1][11];
    }
    double GetDinucSyn31_ACTA()
    {
        return dinucSyn_stat[2][1][12];
    }
    double GetDinucSyn31_ACTC()
    {
        return dinucSyn_stat[2][1][13];
    }
    double GetDinucSyn31_ACTG()
    {
        return dinucSyn_stat[2][1][14];
    }
    double GetDinucSyn31_ACTT()
    {
        return dinucSyn_stat[2][1][15];
    }
    double GetDinucSyn31_AGAA()
    {
        return dinucSyn_stat[2][2][0];
    }
    double GetDinucSyn31_AGAC()
    {
        return dinucSyn_stat[2][2][1];
    }
    double GetDinucSyn31_AGAG()
    {
        return dinucSyn_stat[2][2][2];
    }
    double GetDinucSyn31_AGAT()
    {
        return dinucSyn_stat[2][2][3];
    }
    double GetDinucSyn31_AGCA()
    {
        return dinucSyn_stat[2][2][4];
    }
    double GetDinucSyn31_AGCC()
    {
        return dinucSyn_stat[2][2][5];
    }
    double GetDinucSyn31_AGCG()
    {
        return dinucSyn_stat[2][2][6];
    }
    double GetDinucSyn31_AGCT()
    {
        return dinucSyn_stat[2][2][7];
    }
    double GetDinucSyn31_AGGA()
    {
        return dinucSyn_stat[2][2][8];
    }
    double GetDinucSyn31_AGGC()
    {
        return dinucSyn_stat[2][2][9];
    }
    double GetDinucSyn31_AGGG()
    {
        return dinucSyn_stat[2][2][10];
    }
    double GetDinucSyn31_AGGT()
    {
        return dinucSyn_stat[2][2][11];
    }
    double GetDinucSyn31_AGTA()
    {
        return dinucSyn_stat[2][2][12];
    }
    double GetDinucSyn31_AGTC()
    {
        return dinucSyn_stat[2][2][13];
    }
    double GetDinucSyn31_AGTG()
    {
        return dinucSyn_stat[2][2][14];
    }
    double GetDinucSyn31_AGTT()
    {
        return dinucSyn_stat[2][2][15];
    }
    double GetDinucSyn31_ATAA()
    {
        return dinucSyn_stat[2][3][0];
    }
    double GetDinucSyn31_ATAC()
    {
        return dinucSyn_stat[2][3][1];
    }
    double GetDinucSyn31_ATAG()
    {
        return dinucSyn_stat[2][3][2];
    }
    double GetDinucSyn31_ATAT()
    {
        return dinucSyn_stat[2][3][3];
    }
    double GetDinucSyn31_ATCA()
    {
        return dinucSyn_stat[2][3][4];
    }
    double GetDinucSyn31_ATCC()
    {
        return dinucSyn_stat[2][3][5];
    }
    double GetDinucSyn31_ATCG()
    {
        return dinucSyn_stat[2][3][6];
    }
    double GetDinucSyn31_ATCT()
    {
        return dinucSyn_stat[2][3][7];
    }
    double GetDinucSyn31_ATGA()
    {
        return dinucSyn_stat[2][3][8];
    }
    double GetDinucSyn31_ATGC()
    {
        return dinucSyn_stat[2][3][9];
    }
    double GetDinucSyn31_ATGG()
    {
        return dinucSyn_stat[2][3][10];
    }
    double GetDinucSyn31_ATGT()
    {
        return dinucSyn_stat[2][3][11];
    }
    double GetDinucSyn31_ATTA()
    {
        return dinucSyn_stat[2][3][12];
    }
    double GetDinucSyn31_ATTC()
    {
        return dinucSyn_stat[2][3][13];
    }
    double GetDinucSyn31_ATTG()
    {
        return dinucSyn_stat[2][3][14];
    }
    double GetDinucSyn31_ATTT()
    {
        return dinucSyn_stat[2][3][15];
    }
    double GetDinucSyn31_CAAA()
    {
        return dinucSyn_stat[2][4][0];
    }
    double GetDinucSyn31_CAAC()
    {
        return dinucSyn_stat[2][4][1];
    }
    double GetDinucSyn31_CAAG()
    {
        return dinucSyn_stat[2][4][2];
    }
    double GetDinucSyn31_CAAT()
    {
        return dinucSyn_stat[2][4][3];
    }
    double GetDinucSyn31_CACA()
    {
        return dinucSyn_stat[2][4][4];
    }
    double GetDinucSyn31_CACC()
    {
        return dinucSyn_stat[2][4][5];
    }
    double GetDinucSyn31_CACG()
    {
        return dinucSyn_stat[2][4][6];
    }
    double GetDinucSyn31_CACT()
    {
        return dinucSyn_stat[2][4][7];
    }
    double GetDinucSyn31_CAGA()
    {
        return dinucSyn_stat[2][4][8];
    }
    double GetDinucSyn31_CAGC()
    {
        return dinucSyn_stat[2][4][9];
    }
    double GetDinucSyn31_CAGG()
    {
        return dinucSyn_stat[2][4][10];
    }
    double GetDinucSyn31_CAGT()
    {
        return dinucSyn_stat[2][4][11];
    }
    double GetDinucSyn31_CATA()
    {
        return dinucSyn_stat[2][4][12];
    }
    double GetDinucSyn31_CATC()
    {
        return dinucSyn_stat[2][4][13];
    }
    double GetDinucSyn31_CATG()
    {
        return dinucSyn_stat[2][4][14];
    }
    double GetDinucSyn31_CATT()
    {
        return dinucSyn_stat[2][4][15];
    }
    double GetDinucSyn31_CCAA()
    {
        return dinucSyn_stat[2][5][0];
    }
    double GetDinucSyn31_CCAC()
    {
        return dinucSyn_stat[2][5][1];
    }
    double GetDinucSyn31_CCAG()
    {
        return dinucSyn_stat[2][5][2];
    }
    double GetDinucSyn31_CCAT()
    {
        return dinucSyn_stat[2][5][3];
    }
    double GetDinucSyn31_CCCA()
    {
        return dinucSyn_stat[2][5][4];
    }
    double GetDinucSyn31_CCCC()
    {
        return dinucSyn_stat[2][5][5];
    }
    double GetDinucSyn31_CCCG()
    {
        return dinucSyn_stat[2][5][6];
    }
    double GetDinucSyn31_CCCT()
    {
        return dinucSyn_stat[2][5][7];
    }
    double GetDinucSyn31_CCGA()
    {
        return dinucSyn_stat[2][5][8];
    }
    double GetDinucSyn31_CCGC()
    {
        return dinucSyn_stat[2][5][9];
    }
    double GetDinucSyn31_CCGG()
    {
        return dinucSyn_stat[2][5][10];
    }
    double GetDinucSyn31_CCGT()
    {
        return dinucSyn_stat[2][5][11];
    }
    double GetDinucSyn31_CCTA()
    {
        return dinucSyn_stat[2][5][12];
    }
    double GetDinucSyn31_CCTC()
    {
        return dinucSyn_stat[2][5][13];
    }
    double GetDinucSyn31_CCTG()
    {
        return dinucSyn_stat[2][5][14];
    }
    double GetDinucSyn31_CCTT()
    {
        return dinucSyn_stat[2][5][15];
    }
    double GetDinucSyn31_CGAA()
    {
        return dinucSyn_stat[2][6][0];
    }
    double GetDinucSyn31_CGAC()
    {
        return dinucSyn_stat[2][6][1];
    }
    double GetDinucSyn31_CGAG()
    {
        return dinucSyn_stat[2][6][2];
    }
    double GetDinucSyn31_CGAT()
    {
        return dinucSyn_stat[2][6][3];
    }
    double GetDinucSyn31_CGCA()
    {
        return dinucSyn_stat[2][6][4];
    }
    double GetDinucSyn31_CGCC()
    {
        return dinucSyn_stat[2][6][5];
    }
    double GetDinucSyn31_CGCG()
    {
        return dinucSyn_stat[2][6][6];
    }
    double GetDinucSyn31_CGCT()
    {
        return dinucSyn_stat[2][6][7];
    }
    double GetDinucSyn31_CGGA()
    {
        return dinucSyn_stat[2][6][8];
    }
    double GetDinucSyn31_CGGC()
    {
        return dinucSyn_stat[2][6][9];
    }
    double GetDinucSyn31_CGGG()
    {
        return dinucSyn_stat[2][6][10];
    }
    double GetDinucSyn31_CGGT()
    {
        return dinucSyn_stat[2][6][11];
    }
    double GetDinucSyn31_CGTA()
    {
        return dinucSyn_stat[2][6][12];
    }
    double GetDinucSyn31_CGTC()
    {
        return dinucSyn_stat[2][6][13];
    }
    double GetDinucSyn31_CGTG()
    {
        return dinucSyn_stat[2][6][14];
    }
    double GetDinucSyn31_CGTT()
    {
        return dinucSyn_stat[2][6][15];
    }
    double GetDinucSyn31_CTAA()
    {
        return dinucSyn_stat[2][7][0];
    }
    double GetDinucSyn31_CTAC()
    {
        return dinucSyn_stat[2][7][1];
    }
    double GetDinucSyn31_CTAG()
    {
        return dinucSyn_stat[2][7][2];
    }
    double GetDinucSyn31_CTAT()
    {
        return dinucSyn_stat[2][7][3];
    }
    double GetDinucSyn31_CTCA()
    {
        return dinucSyn_stat[2][7][4];
    }
    double GetDinucSyn31_CTCC()
    {
        return dinucSyn_stat[2][7][5];
    }
    double GetDinucSyn31_CTCG()
    {
        return dinucSyn_stat[2][7][6];
    }
    double GetDinucSyn31_CTCT()
    {
        return dinucSyn_stat[2][7][7];
    }
    double GetDinucSyn31_CTGA()
    {
        return dinucSyn_stat[2][7][8];
    }
    double GetDinucSyn31_CTGC()
    {
        return dinucSyn_stat[2][7][9];
    }
    double GetDinucSyn31_CTGG()
    {
        return dinucSyn_stat[2][7][10];
    }
    double GetDinucSyn31_CTGT()
    {
        return dinucSyn_stat[2][7][11];
    }
    double GetDinucSyn31_CTTA()
    {
        return dinucSyn_stat[2][7][12];
    }
    double GetDinucSyn31_CTTC()
    {
        return dinucSyn_stat[2][7][13];
    }
    double GetDinucSyn31_CTTG()
    {
        return dinucSyn_stat[2][7][14];
    }
    double GetDinucSyn31_CTTT()
    {
        return dinucSyn_stat[2][7][15];
    }
    double GetDinucSyn31_GAAA()
    {
        return dinucSyn_stat[2][8][0];
    }
    double GetDinucSyn31_GAAC()
    {
        return dinucSyn_stat[2][8][1];
    }
    double GetDinucSyn31_GAAG()
    {
        return dinucSyn_stat[2][8][2];
    }
    double GetDinucSyn31_GAAT()
    {
        return dinucSyn_stat[2][8][3];
    }
    double GetDinucSyn31_GACA()
    {
        return dinucSyn_stat[2][8][4];
    }
    double GetDinucSyn31_GACC()
    {
        return dinucSyn_stat[2][8][5];
    }
    double GetDinucSyn31_GACG()
    {
        return dinucSyn_stat[2][8][6];
    }
    double GetDinucSyn31_GACT()
    {
        return dinucSyn_stat[2][8][7];
    }
    double GetDinucSyn31_GAGA()
    {
        return dinucSyn_stat[2][8][8];
    }
    double GetDinucSyn31_GAGC()
    {
        return dinucSyn_stat[2][8][9];
    }
    double GetDinucSyn31_GAGG()
    {
        return dinucSyn_stat[2][8][10];
    }
    double GetDinucSyn31_GAGT()
    {
        return dinucSyn_stat[2][8][11];
    }
    double GetDinucSyn31_GATA()
    {
        return dinucSyn_stat[2][8][12];
    }
    double GetDinucSyn31_GATC()
    {
        return dinucSyn_stat[2][8][13];
    }
    double GetDinucSyn31_GATG()
    {
        return dinucSyn_stat[2][8][14];
    }
    double GetDinucSyn31_GATT()
    {
        return dinucSyn_stat[2][8][15];
    }
    double GetDinucSyn31_GCAA()
    {
        return dinucSyn_stat[2][9][0];
    }
    double GetDinucSyn31_GCAC()
    {
        return dinucSyn_stat[2][9][1];
    }
    double GetDinucSyn31_GCAG()
    {
        return dinucSyn_stat[2][9][2];
    }
    double GetDinucSyn31_GCAT()
    {
        return dinucSyn_stat[2][9][3];
    }
    double GetDinucSyn31_GCCA()
    {
        return dinucSyn_stat[2][9][4];
    }
    double GetDinucSyn31_GCCC()
    {
        return dinucSyn_stat[2][9][5];
    }
    double GetDinucSyn31_GCCG()
    {
        return dinucSyn_stat[2][9][6];
    }
    double GetDinucSyn31_GCCT()
    {
        return dinucSyn_stat[2][9][7];
    }
    double GetDinucSyn31_GCGA()
    {
        return dinucSyn_stat[2][9][8];
    }
    double GetDinucSyn31_GCGC()
    {
        return dinucSyn_stat[2][9][9];
    }
    double GetDinucSyn31_GCGG()
    {
        return dinucSyn_stat[2][9][10];
    }
    double GetDinucSyn31_GCGT()
    {
        return dinucSyn_stat[2][9][11];
    }
    double GetDinucSyn31_GCTA()
    {
        return dinucSyn_stat[2][9][12];
    }
    double GetDinucSyn31_GCTC()
    {
        return dinucSyn_stat[2][9][13];
    }
    double GetDinucSyn31_GCTG()
    {
        return dinucSyn_stat[2][9][14];
    }
    double GetDinucSyn31_GCTT()
    {
        return dinucSyn_stat[2][9][15];
    }
    double GetDinucSyn31_GGAA()
    {
        return dinucSyn_stat[2][10][0];
    }
    double GetDinucSyn31_GGAC()
    {
        return dinucSyn_stat[2][10][1];
    }
    double GetDinucSyn31_GGAG()
    {
        return dinucSyn_stat[2][10][2];
    }
    double GetDinucSyn31_GGAT()
    {
        return dinucSyn_stat[2][10][3];
    }
    double GetDinucSyn31_GGCA()
    {
        return dinucSyn_stat[2][10][4];
    }
    double GetDinucSyn31_GGCC()
    {
        return dinucSyn_stat[2][10][5];
    }
    double GetDinucSyn31_GGCG()
    {
        return dinucSyn_stat[2][10][6];
    }
    double GetDinucSyn31_GGCT()
    {
        return dinucSyn_stat[2][10][7];
    }
    double GetDinucSyn31_GGGA()
    {
        return dinucSyn_stat[2][10][8];
    }
    double GetDinucSyn31_GGGC()
    {
        return dinucSyn_stat[2][10][9];
    }
    double GetDinucSyn31_GGGG()
    {
        return dinucSyn_stat[2][10][10];
    }
    double GetDinucSyn31_GGGT()
    {
        return dinucSyn_stat[2][10][11];
    }
    double GetDinucSyn31_GGTA()
    {
        return dinucSyn_stat[2][10][12];
    }
    double GetDinucSyn31_GGTC()
    {
        return dinucSyn_stat[2][10][13];
    }
    double GetDinucSyn31_GGTG()
    {
        return dinucSyn_stat[2][10][14];
    }
    double GetDinucSyn31_GGTT()
    {
        return dinucSyn_stat[2][10][15];
    }
    double GetDinucSyn31_GTAA()
    {
        return dinucSyn_stat[2][11][0];
    }
    double GetDinucSyn31_GTAC()
    {
        return dinucSyn_stat[2][11][1];
    }
    double GetDinucSyn31_GTAG()
    {
        return dinucSyn_stat[2][11][2];
    }
    double GetDinucSyn31_GTAT()
    {
        return dinucSyn_stat[2][11][3];
    }
    double GetDinucSyn31_GTCA()
    {
        return dinucSyn_stat[2][11][4];
    }
    double GetDinucSyn31_GTCC()
    {
        return dinucSyn_stat[2][11][5];
    }
    double GetDinucSyn31_GTCG()
    {
        return dinucSyn_stat[2][11][6];
    }
    double GetDinucSyn31_GTCT()
    {
        return dinucSyn_stat[2][11][7];
    }
    double GetDinucSyn31_GTGA()
    {
        return dinucSyn_stat[2][11][8];
    }
    double GetDinucSyn31_GTGC()
    {
        return dinucSyn_stat[2][11][9];
    }
    double GetDinucSyn31_GTGG()
    {
        return dinucSyn_stat[2][11][10];
    }
    double GetDinucSyn31_GTGT()
    {
        return dinucSyn_stat[2][11][11];
    }
    double GetDinucSyn31_GTTA()
    {
        return dinucSyn_stat[2][11][12];
    }
    double GetDinucSyn31_GTTC()
    {
        return dinucSyn_stat[2][11][13];
    }
    double GetDinucSyn31_GTTG()
    {
        return dinucSyn_stat[2][11][14];
    }
    double GetDinucSyn31_GTTT()
    {
        return dinucSyn_stat[2][11][15];
    }
    double GetDinucSyn31_TAAA()
    {
        return dinucSyn_stat[2][12][0];
    }
    double GetDinucSyn31_TAAC()
    {
        return dinucSyn_stat[2][12][1];
    }
    double GetDinucSyn31_TAAG()
    {
        return dinucSyn_stat[2][12][2];
    }
    double GetDinucSyn31_TAAT()
    {
        return dinucSyn_stat[2][12][3];
    }
    double GetDinucSyn31_TACA()
    {
        return dinucSyn_stat[2][12][4];
    }
    double GetDinucSyn31_TACC()
    {
        return dinucSyn_stat[2][12][5];
    }
    double GetDinucSyn31_TACG()
    {
        return dinucSyn_stat[2][12][6];
    }
    double GetDinucSyn31_TACT()
    {
        return dinucSyn_stat[2][12][7];
    }
    double GetDinucSyn31_TAGA()
    {
        return dinucSyn_stat[2][12][8];
    }
    double GetDinucSyn31_TAGC()
    {
        return dinucSyn_stat[2][12][9];
    }
    double GetDinucSyn31_TAGG()
    {
        return dinucSyn_stat[2][12][10];
    }
    double GetDinucSyn31_TAGT()
    {
        return dinucSyn_stat[2][12][11];
    }
    double GetDinucSyn31_TATA()
    {
        return dinucSyn_stat[2][12][12];
    }
    double GetDinucSyn31_TATC()
    {
        return dinucSyn_stat[2][12][13];
    }
    double GetDinucSyn31_TATG()
    {
        return dinucSyn_stat[2][12][14];
    }
    double GetDinucSyn31_TATT()
    {
        return dinucSyn_stat[2][12][15];
    }
    double GetDinucSyn31_TCAA()
    {
        return dinucSyn_stat[2][13][0];
    }
    double GetDinucSyn31_TCAC()
    {
        return dinucSyn_stat[2][13][1];
    }
    double GetDinucSyn31_TCAG()
    {
        return dinucSyn_stat[2][13][2];
    }
    double GetDinucSyn31_TCAT()
    {
        return dinucSyn_stat[2][13][3];
    }
    double GetDinucSyn31_TCCA()
    {
        return dinucSyn_stat[2][13][4];
    }
    double GetDinucSyn31_TCCC()
    {
        return dinucSyn_stat[2][13][5];
    }
    double GetDinucSyn31_TCCG()
    {
        return dinucSyn_stat[2][13][6];
    }
    double GetDinucSyn31_TCCT()
    {
        return dinucSyn_stat[2][13][7];
    }
    double GetDinucSyn31_TCGA()
    {
        return dinucSyn_stat[2][13][8];
    }
    double GetDinucSyn31_TCGC()
    {
        return dinucSyn_stat[2][13][9];
    }
    double GetDinucSyn31_TCGG()
    {
        return dinucSyn_stat[2][13][10];
    }
    double GetDinucSyn31_TCGT()
    {
        return dinucSyn_stat[2][13][11];
    }
    double GetDinucSyn31_TCTA()
    {
        return dinucSyn_stat[2][13][12];
    }
    double GetDinucSyn31_TCTC()
    {
        return dinucSyn_stat[2][13][13];
    }
    double GetDinucSyn31_TCTG()
    {
        return dinucSyn_stat[2][13][14];
    }
    double GetDinucSyn31_TCTT()
    {
        return dinucSyn_stat[2][13][15];
    }
    double GetDinucSyn31_TGAA()
    {
        return dinucSyn_stat[2][14][0];
    }
    double GetDinucSyn31_TGAC()
    {
        return dinucSyn_stat[2][14][1];
    }
    double GetDinucSyn31_TGAG()
    {
        return dinucSyn_stat[2][14][2];
    }
    double GetDinucSyn31_TGAT()
    {
        return dinucSyn_stat[2][14][3];
    }
    double GetDinucSyn31_TGCA()
    {
        return dinucSyn_stat[2][14][4];
    }
    double GetDinucSyn31_TGCC()
    {
        return dinucSyn_stat[2][14][5];
    }
    double GetDinucSyn31_TGCG()
    {
        return dinucSyn_stat[2][14][6];
    }
    double GetDinucSyn31_TGCT()
    {
        return dinucSyn_stat[2][14][7];
    }
    double GetDinucSyn31_TGGA()
    {
        return dinucSyn_stat[2][14][8];
    }
    double GetDinucSyn31_TGGC()
    {
        return dinucSyn_stat[2][14][9];
    }
    double GetDinucSyn31_TGGG()
    {
        return dinucSyn_stat[2][14][10];
    }
    double GetDinucSyn31_TGGT()
    {
        return dinucSyn_stat[2][14][11];
    }
    double GetDinucSyn31_TGTA()
    {
        return dinucSyn_stat[2][14][12];
    }
    double GetDinucSyn31_TGTC()
    {
        return dinucSyn_stat[2][14][13];
    }
    double GetDinucSyn31_TGTG()
    {
        return dinucSyn_stat[2][14][14];
    }
    double GetDinucSyn31_TGTT()
    {
        return dinucSyn_stat[2][14][15];
    }
    double GetDinucSyn31_TTAA()
    {
        return dinucSyn_stat[2][15][0];
    }
    double GetDinucSyn31_TTAC()
    {
        return dinucSyn_stat[2][15][1];
    }
    double GetDinucSyn31_TTAG()
    {
        return dinucSyn_stat[2][15][2];
    }
    double GetDinucSyn31_TTAT()
    {
        return dinucSyn_stat[2][15][3];
    }
    double GetDinucSyn31_TTCA()
    {
        return dinucSyn_stat[2][15][4];
    }
    double GetDinucSyn31_TTCC()
    {
        return dinucSyn_stat[2][15][5];
    }
    double GetDinucSyn31_TTCG()
    {
        return dinucSyn_stat[2][15][6];
    }
    double GetDinucSyn31_TTCT()
    {
        return dinucSyn_stat[2][15][7];
    }
    double GetDinucSyn31_TTGA()
    {
        return dinucSyn_stat[2][15][8];
    }
    double GetDinucSyn31_TTGC()
    {
        return dinucSyn_stat[2][15][9];
    }
    double GetDinucSyn31_TTGG()
    {
        return dinucSyn_stat[2][15][10];
    }
    double GetDinucSyn31_TTGT()
    {
        return dinucSyn_stat[2][15][11];
    }
    double GetDinucSyn31_TTTA()
    {
        return dinucSyn_stat[2][15][12];
    }
    double GetDinucSyn31_TTTC()
    {
        return dinucSyn_stat[2][15][13];
    }
    double GetDinucSyn31_TTTG()
    {
        return dinucSyn_stat[2][15][14];
    }
    double GetDinucSyn31_TTTT()
    {
        return dinucSyn_stat[2][15][15];
    }

    double GetDinucNSyn_AAAA()
    {
        return dinucNSyn_stat[3][0][0];
    }
    double GetDinucNSyn_AAAC()
    {
        return dinucNSyn_stat[3][0][1];
    }
    double GetDinucNSyn_AAAG()
    {
        return dinucNSyn_stat[3][0][2];
    }
    double GetDinucNSyn_AAAT()
    {
        return dinucNSyn_stat[3][0][3];
    }
    double GetDinucNSyn_AACA()
    {
        return dinucNSyn_stat[3][0][4];
    }
    double GetDinucNSyn_AACC()
    {
        return dinucNSyn_stat[3][0][5];
    }
    double GetDinucNSyn_AACG()
    {
        return dinucNSyn_stat[3][0][6];
    }
    double GetDinucNSyn_AACT()
    {
        return dinucNSyn_stat[3][0][7];
    }
    double GetDinucNSyn_AAGA()
    {
        return dinucNSyn_stat[3][0][8];
    }
    double GetDinucNSyn_AAGC()
    {
        return dinucNSyn_stat[3][0][9];
    }
    double GetDinucNSyn_AAGG()
    {
        return dinucNSyn_stat[3][0][10];
    }
    double GetDinucNSyn_AAGT()
    {
        return dinucNSyn_stat[3][0][11];
    }
    double GetDinucNSyn_AATA()
    {
        return dinucNSyn_stat[3][0][12];
    }
    double GetDinucNSyn_AATC()
    {
        return dinucNSyn_stat[3][0][13];
    }
    double GetDinucNSyn_AATG()
    {
        return dinucNSyn_stat[3][0][14];
    }
    double GetDinucNSyn_AATT()
    {
        return dinucNSyn_stat[3][0][15];
    }
    double GetDinucNSyn_ACAA()
    {
        return dinucNSyn_stat[3][1][0];
    }
    double GetDinucNSyn_ACAC()
    {
        return dinucNSyn_stat[3][1][1];
    }
    double GetDinucNSyn_ACAG()
    {
        return dinucNSyn_stat[3][1][2];
    }
    double GetDinucNSyn_ACAT()
    {
        return dinucNSyn_stat[3][1][3];
    }
    double GetDinucNSyn_ACCA()
    {
        return dinucNSyn_stat[3][1][4];
    }
    double GetDinucNSyn_ACCC()
    {
        return dinucNSyn_stat[3][1][5];
    }
    double GetDinucNSyn_ACCG()
    {
        return dinucNSyn_stat[3][1][6];
    }
    double GetDinucNSyn_ACCT()
    {
        return dinucNSyn_stat[3][1][7];
    }
    double GetDinucNSyn_ACGA()
    {
        return dinucNSyn_stat[3][1][8];
    }
    double GetDinucNSyn_ACGC()
    {
        return dinucNSyn_stat[3][1][9];
    }
    double GetDinucNSyn_ACGG()
    {
        return dinucNSyn_stat[3][1][10];
    }
    double GetDinucNSyn_ACGT()
    {
        return dinucNSyn_stat[3][1][11];
    }
    double GetDinucNSyn_ACTA()
    {
        return dinucNSyn_stat[3][1][12];
    }
    double GetDinucNSyn_ACTC()
    {
        return dinucNSyn_stat[3][1][13];
    }
    double GetDinucNSyn_ACTG()
    {
        return dinucNSyn_stat[3][1][14];
    }
    double GetDinucNSyn_ACTT()
    {
        return dinucNSyn_stat[3][1][15];
    }
    double GetDinucNSyn_AGAA()
    {
        return dinucNSyn_stat[3][2][0];
    }
    double GetDinucNSyn_AGAC()
    {
        return dinucNSyn_stat[3][2][1];
    }
    double GetDinucNSyn_AGAG()
    {
        return dinucNSyn_stat[3][2][2];
    }
    double GetDinucNSyn_AGAT()
    {
        return dinucNSyn_stat[3][2][3];
    }
    double GetDinucNSyn_AGCA()
    {
        return dinucNSyn_stat[3][2][4];
    }
    double GetDinucNSyn_AGCC()
    {
        return dinucNSyn_stat[3][2][5];
    }
    double GetDinucNSyn_AGCG()
    {
        return dinucNSyn_stat[3][2][6];
    }
    double GetDinucNSyn_AGCT()
    {
        return dinucNSyn_stat[3][2][7];
    }
    double GetDinucNSyn_AGGA()
    {
        return dinucNSyn_stat[3][2][8];
    }
    double GetDinucNSyn_AGGC()
    {
        return dinucNSyn_stat[3][2][9];
    }
    double GetDinucNSyn_AGGG()
    {
        return dinucNSyn_stat[3][2][10];
    }
    double GetDinucNSyn_AGGT()
    {
        return dinucNSyn_stat[3][2][11];
    }
    double GetDinucNSyn_AGTA()
    {
        return dinucNSyn_stat[3][2][12];
    }
    double GetDinucNSyn_AGTC()
    {
        return dinucNSyn_stat[3][2][13];
    }
    double GetDinucNSyn_AGTG()
    {
        return dinucNSyn_stat[3][2][14];
    }
    double GetDinucNSyn_AGTT()
    {
        return dinucNSyn_stat[3][2][15];
    }
    double GetDinucNSyn_ATAA()
    {
        return dinucNSyn_stat[3][3][0];
    }
    double GetDinucNSyn_ATAC()
    {
        return dinucNSyn_stat[3][3][1];
    }
    double GetDinucNSyn_ATAG()
    {
        return dinucNSyn_stat[3][3][2];
    }
    double GetDinucNSyn_ATAT()
    {
        return dinucNSyn_stat[3][3][3];
    }
    double GetDinucNSyn_ATCA()
    {
        return dinucNSyn_stat[3][3][4];
    }
    double GetDinucNSyn_ATCC()
    {
        return dinucNSyn_stat[3][3][5];
    }
    double GetDinucNSyn_ATCG()
    {
        return dinucNSyn_stat[3][3][6];
    }
    double GetDinucNSyn_ATCT()
    {
        return dinucNSyn_stat[3][3][7];
    }
    double GetDinucNSyn_ATGA()
    {
        return dinucNSyn_stat[3][3][8];
    }
    double GetDinucNSyn_ATGC()
    {
        return dinucNSyn_stat[3][3][9];
    }
    double GetDinucNSyn_ATGG()
    {
        return dinucNSyn_stat[3][3][10];
    }
    double GetDinucNSyn_ATGT()
    {
        return dinucNSyn_stat[3][3][11];
    }
    double GetDinucNSyn_ATTA()
    {
        return dinucNSyn_stat[3][3][12];
    }
    double GetDinucNSyn_ATTC()
    {
        return dinucNSyn_stat[3][3][13];
    }
    double GetDinucNSyn_ATTG()
    {
        return dinucNSyn_stat[3][3][14];
    }
    double GetDinucNSyn_ATTT()
    {
        return dinucNSyn_stat[3][3][15];
    }
    double GetDinucNSyn_CAAA()
    {
        return dinucNSyn_stat[3][4][0];
    }
    double GetDinucNSyn_CAAC()
    {
        return dinucNSyn_stat[3][4][1];
    }
    double GetDinucNSyn_CAAG()
    {
        return dinucNSyn_stat[3][4][2];
    }
    double GetDinucNSyn_CAAT()
    {
        return dinucNSyn_stat[3][4][3];
    }
    double GetDinucNSyn_CACA()
    {
        return dinucNSyn_stat[3][4][4];
    }
    double GetDinucNSyn_CACC()
    {
        return dinucNSyn_stat[3][4][5];
    }
    double GetDinucNSyn_CACG()
    {
        return dinucNSyn_stat[3][4][6];
    }
    double GetDinucNSyn_CACT()
    {
        return dinucNSyn_stat[3][4][7];
    }
    double GetDinucNSyn_CAGA()
    {
        return dinucNSyn_stat[3][4][8];
    }
    double GetDinucNSyn_CAGC()
    {
        return dinucNSyn_stat[3][4][9];
    }
    double GetDinucNSyn_CAGG()
    {
        return dinucNSyn_stat[3][4][10];
    }
    double GetDinucNSyn_CAGT()
    {
        return dinucNSyn_stat[3][4][11];
    }
    double GetDinucNSyn_CATA()
    {
        return dinucNSyn_stat[3][4][12];
    }
    double GetDinucNSyn_CATC()
    {
        return dinucNSyn_stat[3][4][13];
    }
    double GetDinucNSyn_CATG()
    {
        return dinucNSyn_stat[3][4][14];
    }
    double GetDinucNSyn_CATT()
    {
        return dinucNSyn_stat[3][4][15];
    }
    double GetDinucNSyn_CCAA()
    {
        return dinucNSyn_stat[3][5][0];
    }
    double GetDinucNSyn_CCAC()
    {
        return dinucNSyn_stat[3][5][1];
    }
    double GetDinucNSyn_CCAG()
    {
        return dinucNSyn_stat[3][5][2];
    }
    double GetDinucNSyn_CCAT()
    {
        return dinucNSyn_stat[3][5][3];
    }
    double GetDinucNSyn_CCCA()
    {
        return dinucNSyn_stat[3][5][4];
    }
    double GetDinucNSyn_CCCC()
    {
        return dinucNSyn_stat[3][5][5];
    }
    double GetDinucNSyn_CCCG()
    {
        return dinucNSyn_stat[3][5][6];
    }
    double GetDinucNSyn_CCCT()
    {
        return dinucNSyn_stat[3][5][7];
    }
    double GetDinucNSyn_CCGA()
    {
        return dinucNSyn_stat[3][5][8];
    }
    double GetDinucNSyn_CCGC()
    {
        return dinucNSyn_stat[3][5][9];
    }
    double GetDinucNSyn_CCGG()
    {
        return dinucNSyn_stat[3][5][10];
    }
    double GetDinucNSyn_CCGT()
    {
        return dinucNSyn_stat[3][5][11];
    }
    double GetDinucNSyn_CCTA()
    {
        return dinucNSyn_stat[3][5][12];
    }
    double GetDinucNSyn_CCTC()
    {
        return dinucNSyn_stat[3][5][13];
    }
    double GetDinucNSyn_CCTG()
    {
        return dinucNSyn_stat[3][5][14];
    }
    double GetDinucNSyn_CCTT()
    {
        return dinucNSyn_stat[3][5][15];
    }
    double GetDinucNSyn_CGAA()
    {
        return dinucNSyn_stat[3][6][0];
    }
    double GetDinucNSyn_CGAC()
    {
        return dinucNSyn_stat[3][6][1];
    }
    double GetDinucNSyn_CGAG()
    {
        return dinucNSyn_stat[3][6][2];
    }
    double GetDinucNSyn_CGAT()
    {
        return dinucNSyn_stat[3][6][3];
    }
    double GetDinucNSyn_CGCA()
    {
        return dinucNSyn_stat[3][6][4];
    }
    double GetDinucNSyn_CGCC()
    {
        return dinucNSyn_stat[3][6][5];
    }
    double GetDinucNSyn_CGCG()
    {
        return dinucNSyn_stat[3][6][6];
    }
    double GetDinucNSyn_CGCT()
    {
        return dinucNSyn_stat[3][6][7];
    }
    double GetDinucNSyn_CGGA()
    {
        return dinucNSyn_stat[3][6][8];
    }
    double GetDinucNSyn_CGGC()
    {
        return dinucNSyn_stat[3][6][9];
    }
    double GetDinucNSyn_CGGG()
    {
        return dinucNSyn_stat[3][6][10];
    }
    double GetDinucNSyn_CGGT()
    {
        return dinucNSyn_stat[3][6][11];
    }
    double GetDinucNSyn_CGTA()
    {
        return dinucNSyn_stat[3][6][12];
    }
    double GetDinucNSyn_CGTC()
    {
        return dinucNSyn_stat[3][6][13];
    }
    double GetDinucNSyn_CGTG()
    {
        return dinucNSyn_stat[3][6][14];
    }
    double GetDinucNSyn_CGTT()
    {
        return dinucNSyn_stat[3][6][15];
    }
    double GetDinucNSyn_CTAA()
    {
        return dinucNSyn_stat[3][7][0];
    }
    double GetDinucNSyn_CTAC()
    {
        return dinucNSyn_stat[3][7][1];
    }
    double GetDinucNSyn_CTAG()
    {
        return dinucNSyn_stat[3][7][2];
    }
    double GetDinucNSyn_CTAT()
    {
        return dinucNSyn_stat[3][7][3];
    }
    double GetDinucNSyn_CTCA()
    {
        return dinucNSyn_stat[3][7][4];
    }
    double GetDinucNSyn_CTCC()
    {
        return dinucNSyn_stat[3][7][5];
    }
    double GetDinucNSyn_CTCG()
    {
        return dinucNSyn_stat[3][7][6];
    }
    double GetDinucNSyn_CTCT()
    {
        return dinucNSyn_stat[3][7][7];
    }
    double GetDinucNSyn_CTGA()
    {
        return dinucNSyn_stat[3][7][8];
    }
    double GetDinucNSyn_CTGC()
    {
        return dinucNSyn_stat[3][7][9];
    }
    double GetDinucNSyn_CTGG()
    {
        return dinucNSyn_stat[3][7][10];
    }
    double GetDinucNSyn_CTGT()
    {
        return dinucNSyn_stat[3][7][11];
    }
    double GetDinucNSyn_CTTA()
    {
        return dinucNSyn_stat[3][7][12];
    }
    double GetDinucNSyn_CTTC()
    {
        return dinucNSyn_stat[3][7][13];
    }
    double GetDinucNSyn_CTTG()
    {
        return dinucNSyn_stat[3][7][14];
    }
    double GetDinucNSyn_CTTT()
    {
        return dinucNSyn_stat[3][7][15];
    }
    double GetDinucNSyn_GAAA()
    {
        return dinucNSyn_stat[3][8][0];
    }
    double GetDinucNSyn_GAAC()
    {
        return dinucNSyn_stat[3][8][1];
    }
    double GetDinucNSyn_GAAG()
    {
        return dinucNSyn_stat[3][8][2];
    }
    double GetDinucNSyn_GAAT()
    {
        return dinucNSyn_stat[3][8][3];
    }
    double GetDinucNSyn_GACA()
    {
        return dinucNSyn_stat[3][8][4];
    }
    double GetDinucNSyn_GACC()
    {
        return dinucNSyn_stat[3][8][5];
    }
    double GetDinucNSyn_GACG()
    {
        return dinucNSyn_stat[3][8][6];
    }
    double GetDinucNSyn_GACT()
    {
        return dinucNSyn_stat[3][8][7];
    }
    double GetDinucNSyn_GAGA()
    {
        return dinucNSyn_stat[3][8][8];
    }
    double GetDinucNSyn_GAGC()
    {
        return dinucNSyn_stat[3][8][9];
    }
    double GetDinucNSyn_GAGG()
    {
        return dinucNSyn_stat[3][8][10];
    }
    double GetDinucNSyn_GAGT()
    {
        return dinucNSyn_stat[3][8][11];
    }
    double GetDinucNSyn_GATA()
    {
        return dinucNSyn_stat[3][8][12];
    }
    double GetDinucNSyn_GATC()
    {
        return dinucNSyn_stat[3][8][13];
    }
    double GetDinucNSyn_GATG()
    {
        return dinucNSyn_stat[3][8][14];
    }
    double GetDinucNSyn_GATT()
    {
        return dinucNSyn_stat[3][8][15];
    }
    double GetDinucNSyn_GCAA()
    {
        return dinucNSyn_stat[3][9][0];
    }
    double GetDinucNSyn_GCAC()
    {
        return dinucNSyn_stat[3][9][1];
    }
    double GetDinucNSyn_GCAG()
    {
        return dinucNSyn_stat[3][9][2];
    }
    double GetDinucNSyn_GCAT()
    {
        return dinucNSyn_stat[3][9][3];
    }
    double GetDinucNSyn_GCCA()
    {
        return dinucNSyn_stat[3][9][4];
    }
    double GetDinucNSyn_GCCC()
    {
        return dinucNSyn_stat[3][9][5];
    }
    double GetDinucNSyn_GCCG()
    {
        return dinucNSyn_stat[3][9][6];
    }
    double GetDinucNSyn_GCCT()
    {
        return dinucNSyn_stat[3][9][7];
    }
    double GetDinucNSyn_GCGA()
    {
        return dinucNSyn_stat[3][9][8];
    }
    double GetDinucNSyn_GCGC()
    {
        return dinucNSyn_stat[3][9][9];
    }
    double GetDinucNSyn_GCGG()
    {
        return dinucNSyn_stat[3][9][10];
    }
    double GetDinucNSyn_GCGT()
    {
        return dinucNSyn_stat[3][9][11];
    }
    double GetDinucNSyn_GCTA()
    {
        return dinucNSyn_stat[3][9][12];
    }
    double GetDinucNSyn_GCTC()
    {
        return dinucNSyn_stat[3][9][13];
    }
    double GetDinucNSyn_GCTG()
    {
        return dinucNSyn_stat[3][9][14];
    }
    double GetDinucNSyn_GCTT()
    {
        return dinucNSyn_stat[3][9][15];
    }
    double GetDinucNSyn_GGAA()
    {
        return dinucNSyn_stat[3][10][0];
    }
    double GetDinucNSyn_GGAC()
    {
        return dinucNSyn_stat[3][10][1];
    }
    double GetDinucNSyn_GGAG()
    {
        return dinucNSyn_stat[3][10][2];
    }
    double GetDinucNSyn_GGAT()
    {
        return dinucNSyn_stat[3][10][3];
    }
    double GetDinucNSyn_GGCA()
    {
        return dinucNSyn_stat[3][10][4];
    }
    double GetDinucNSyn_GGCC()
    {
        return dinucNSyn_stat[3][10][5];
    }
    double GetDinucNSyn_GGCG()
    {
        return dinucNSyn_stat[3][10][6];
    }
    double GetDinucNSyn_GGCT()
    {
        return dinucNSyn_stat[3][10][7];
    }
    double GetDinucNSyn_GGGA()
    {
        return dinucNSyn_stat[3][10][8];
    }
    double GetDinucNSyn_GGGC()
    {
        return dinucNSyn_stat[3][10][9];
    }
    double GetDinucNSyn_GGGG()
    {
        return dinucNSyn_stat[3][10][10];
    }
    double GetDinucNSyn_GGGT()
    {
        return dinucNSyn_stat[3][10][11];
    }
    double GetDinucNSyn_GGTA()
    {
        return dinucNSyn_stat[3][10][12];
    }
    double GetDinucNSyn_GGTC()
    {
        return dinucNSyn_stat[3][10][13];
    }
    double GetDinucNSyn_GGTG()
    {
        return dinucNSyn_stat[3][10][14];
    }
    double GetDinucNSyn_GGTT()
    {
        return dinucNSyn_stat[3][10][15];
    }
    double GetDinucNSyn_GTAA()
    {
        return dinucNSyn_stat[3][11][0];
    }
    double GetDinucNSyn_GTAC()
    {
        return dinucNSyn_stat[3][11][1];
    }
    double GetDinucNSyn_GTAG()
    {
        return dinucNSyn_stat[3][11][2];
    }
    double GetDinucNSyn_GTAT()
    {
        return dinucNSyn_stat[3][11][3];
    }
    double GetDinucNSyn_GTCA()
    {
        return dinucNSyn_stat[3][11][4];
    }
    double GetDinucNSyn_GTCC()
    {
        return dinucNSyn_stat[3][11][5];
    }
    double GetDinucNSyn_GTCG()
    {
        return dinucNSyn_stat[3][11][6];
    }
    double GetDinucNSyn_GTCT()
    {
        return dinucNSyn_stat[3][11][7];
    }
    double GetDinucNSyn_GTGA()
    {
        return dinucNSyn_stat[3][11][8];
    }
    double GetDinucNSyn_GTGC()
    {
        return dinucNSyn_stat[3][11][9];
    }
    double GetDinucNSyn_GTGG()
    {
        return dinucNSyn_stat[3][11][10];
    }
    double GetDinucNSyn_GTGT()
    {
        return dinucNSyn_stat[3][11][11];
    }
    double GetDinucNSyn_GTTA()
    {
        return dinucNSyn_stat[3][11][12];
    }
    double GetDinucNSyn_GTTC()
    {
        return dinucNSyn_stat[3][11][13];
    }
    double GetDinucNSyn_GTTG()
    {
        return dinucNSyn_stat[3][11][14];
    }
    double GetDinucNSyn_GTTT()
    {
        return dinucNSyn_stat[3][11][15];
    }
    double GetDinucNSyn_TAAA()
    {
        return dinucNSyn_stat[3][12][0];
    }
    double GetDinucNSyn_TAAC()
    {
        return dinucNSyn_stat[3][12][1];
    }
    double GetDinucNSyn_TAAG()
    {
        return dinucNSyn_stat[3][12][2];
    }
    double GetDinucNSyn_TAAT()
    {
        return dinucNSyn_stat[3][12][3];
    }
    double GetDinucNSyn_TACA()
    {
        return dinucNSyn_stat[3][12][4];
    }
    double GetDinucNSyn_TACC()
    {
        return dinucNSyn_stat[3][12][5];
    }
    double GetDinucNSyn_TACG()
    {
        return dinucNSyn_stat[3][12][6];
    }
    double GetDinucNSyn_TACT()
    {
        return dinucNSyn_stat[3][12][7];
    }
    double GetDinucNSyn_TAGA()
    {
        return dinucNSyn_stat[3][12][8];
    }
    double GetDinucNSyn_TAGC()
    {
        return dinucNSyn_stat[3][12][9];
    }
    double GetDinucNSyn_TAGG()
    {
        return dinucNSyn_stat[3][12][10];
    }
    double GetDinucNSyn_TAGT()
    {
        return dinucNSyn_stat[3][12][11];
    }
    double GetDinucNSyn_TATA()
    {
        return dinucNSyn_stat[3][12][12];
    }
    double GetDinucNSyn_TATC()
    {
        return dinucNSyn_stat[3][12][13];
    }
    double GetDinucNSyn_TATG()
    {
        return dinucNSyn_stat[3][12][14];
    }
    double GetDinucNSyn_TATT()
    {
        return dinucNSyn_stat[3][12][15];
    }
    double GetDinucNSyn_TCAA()
    {
        return dinucNSyn_stat[3][13][0];
    }
    double GetDinucNSyn_TCAC()
    {
        return dinucNSyn_stat[3][13][1];
    }
    double GetDinucNSyn_TCAG()
    {
        return dinucNSyn_stat[3][13][2];
    }
    double GetDinucNSyn_TCAT()
    {
        return dinucNSyn_stat[3][13][3];
    }
    double GetDinucNSyn_TCCA()
    {
        return dinucNSyn_stat[3][13][4];
    }
    double GetDinucNSyn_TCCC()
    {
        return dinucNSyn_stat[3][13][5];
    }
    double GetDinucNSyn_TCCG()
    {
        return dinucNSyn_stat[3][13][6];
    }
    double GetDinucNSyn_TCCT()
    {
        return dinucNSyn_stat[3][13][7];
    }
    double GetDinucNSyn_TCGA()
    {
        return dinucNSyn_stat[3][13][8];
    }
    double GetDinucNSyn_TCGC()
    {
        return dinucNSyn_stat[3][13][9];
    }
    double GetDinucNSyn_TCGG()
    {
        return dinucNSyn_stat[3][13][10];
    }
    double GetDinucNSyn_TCGT()
    {
        return dinucNSyn_stat[3][13][11];
    }
    double GetDinucNSyn_TCTA()
    {
        return dinucNSyn_stat[3][13][12];
    }
    double GetDinucNSyn_TCTC()
    {
        return dinucNSyn_stat[3][13][13];
    }
    double GetDinucNSyn_TCTG()
    {
        return dinucNSyn_stat[3][13][14];
    }
    double GetDinucNSyn_TCTT()
    {
        return dinucNSyn_stat[3][13][15];
    }
    double GetDinucNSyn_TGAA()
    {
        return dinucNSyn_stat[3][14][0];
    }
    double GetDinucNSyn_TGAC()
    {
        return dinucNSyn_stat[3][14][1];
    }
    double GetDinucNSyn_TGAG()
    {
        return dinucNSyn_stat[3][14][2];
    }
    double GetDinucNSyn_TGAT()
    {
        return dinucNSyn_stat[3][14][3];
    }
    double GetDinucNSyn_TGCA()
    {
        return dinucNSyn_stat[3][14][4];
    }
    double GetDinucNSyn_TGCC()
    {
        return dinucNSyn_stat[3][14][5];
    }
    double GetDinucNSyn_TGCG()
    {
        return dinucNSyn_stat[3][14][6];
    }
    double GetDinucNSyn_TGCT()
    {
        return dinucNSyn_stat[3][14][7];
    }
    double GetDinucNSyn_TGGA()
    {
        return dinucNSyn_stat[3][14][8];
    }
    double GetDinucNSyn_TGGC()
    {
        return dinucNSyn_stat[3][14][9];
    }
    double GetDinucNSyn_TGGG()
    {
        return dinucNSyn_stat[3][14][10];
    }
    double GetDinucNSyn_TGGT()
    {
        return dinucNSyn_stat[3][14][11];
    }
    double GetDinucNSyn_TGTA()
    {
        return dinucNSyn_stat[3][14][12];
    }
    double GetDinucNSyn_TGTC()
    {
        return dinucNSyn_stat[3][14][13];
    }
    double GetDinucNSyn_TGTG()
    {
        return dinucNSyn_stat[3][14][14];
    }
    double GetDinucNSyn_TGTT()
    {
        return dinucNSyn_stat[3][14][15];
    }
    double GetDinucNSyn_TTAA()
    {
        return dinucNSyn_stat[3][15][0];
    }
    double GetDinucNSyn_TTAC()
    {
        return dinucNSyn_stat[3][15][1];
    }
    double GetDinucNSyn_TTAG()
    {
        return dinucNSyn_stat[3][15][2];
    }
    double GetDinucNSyn_TTAT()
    {
        return dinucNSyn_stat[3][15][3];
    }
    double GetDinucNSyn_TTCA()
    {
        return dinucNSyn_stat[3][15][4];
    }
    double GetDinucNSyn_TTCC()
    {
        return dinucNSyn_stat[3][15][5];
    }
    double GetDinucNSyn_TTCG()
    {
        return dinucNSyn_stat[3][15][6];
    }
    double GetDinucNSyn_TTCT()
    {
        return dinucNSyn_stat[3][15][7];
    }
    double GetDinucNSyn_TTGA()
    {
        return dinucNSyn_stat[3][15][8];
    }
    double GetDinucNSyn_TTGC()
    {
        return dinucNSyn_stat[3][15][9];
    }
    double GetDinucNSyn_TTGG()
    {
        return dinucNSyn_stat[3][15][10];
    }
    double GetDinucNSyn_TTGT()
    {
        return dinucNSyn_stat[3][15][11];
    }
    double GetDinucNSyn_TTTA()
    {
        return dinucNSyn_stat[3][15][12];
    }
    double GetDinucNSyn_TTTC()
    {
        return dinucNSyn_stat[3][15][13];
    }
    double GetDinucNSyn_TTTG()
    {
        return dinucNSyn_stat[3][15][14];
    }
    double GetDinucNSyn_TTTT()
    {
        return dinucNSyn_stat[3][15][15];
    }
    double GetDinucNSyn12_AAAA()
    {
        return dinucNSyn_stat[0][0][0];
    }
    double GetDinucNSyn12_AAAC()
    {
        return dinucNSyn_stat[0][0][1];
    }
    double GetDinucNSyn12_AAAG()
    {
        return dinucNSyn_stat[0][0][2];
    }
    double GetDinucNSyn12_AAAT()
    {
        return dinucNSyn_stat[0][0][3];
    }
    double GetDinucNSyn12_AACA()
    {
        return dinucNSyn_stat[0][0][4];
    }
    double GetDinucNSyn12_AACC()
    {
        return dinucNSyn_stat[0][0][5];
    }
    double GetDinucNSyn12_AACG()
    {
        return dinucNSyn_stat[0][0][6];
    }
    double GetDinucNSyn12_AACT()
    {
        return dinucNSyn_stat[0][0][7];
    }
    double GetDinucNSyn12_AAGA()
    {
        return dinucNSyn_stat[0][0][8];
    }
    double GetDinucNSyn12_AAGC()
    {
        return dinucNSyn_stat[0][0][9];
    }
    double GetDinucNSyn12_AAGG()
    {
        return dinucNSyn_stat[0][0][10];
    }
    double GetDinucNSyn12_AAGT()
    {
        return dinucNSyn_stat[0][0][11];
    }
    double GetDinucNSyn12_AATA()
    {
        return dinucNSyn_stat[0][0][12];
    }
    double GetDinucNSyn12_AATC()
    {
        return dinucNSyn_stat[0][0][13];
    }
    double GetDinucNSyn12_AATG()
    {
        return dinucNSyn_stat[0][0][14];
    }
    double GetDinucNSyn12_AATT()
    {
        return dinucNSyn_stat[0][0][15];
    }
    double GetDinucNSyn12_ACAA()
    {
        return dinucNSyn_stat[0][1][0];
    }
    double GetDinucNSyn12_ACAC()
    {
        return dinucNSyn_stat[0][1][1];
    }
    double GetDinucNSyn12_ACAG()
    {
        return dinucNSyn_stat[0][1][2];
    }
    double GetDinucNSyn12_ACAT()
    {
        return dinucNSyn_stat[0][1][3];
    }
    double GetDinucNSyn12_ACCA()
    {
        return dinucNSyn_stat[0][1][4];
    }
    double GetDinucNSyn12_ACCC()
    {
        return dinucNSyn_stat[0][1][5];
    }
    double GetDinucNSyn12_ACCG()
    {
        return dinucNSyn_stat[0][1][6];
    }
    double GetDinucNSyn12_ACCT()
    {
        return dinucNSyn_stat[0][1][7];
    }
    double GetDinucNSyn12_ACGA()
    {
        return dinucNSyn_stat[0][1][8];
    }
    double GetDinucNSyn12_ACGC()
    {
        return dinucNSyn_stat[0][1][9];
    }
    double GetDinucNSyn12_ACGG()
    {
        return dinucNSyn_stat[0][1][10];
    }
    double GetDinucNSyn12_ACGT()
    {
        return dinucNSyn_stat[0][1][11];
    }
    double GetDinucNSyn12_ACTA()
    {
        return dinucNSyn_stat[0][1][12];
    }
    double GetDinucNSyn12_ACTC()
    {
        return dinucNSyn_stat[0][1][13];
    }
    double GetDinucNSyn12_ACTG()
    {
        return dinucNSyn_stat[0][1][14];
    }
    double GetDinucNSyn12_ACTT()
    {
        return dinucNSyn_stat[0][1][15];
    }
    double GetDinucNSyn12_AGAA()
    {
        return dinucNSyn_stat[0][2][0];
    }
    double GetDinucNSyn12_AGAC()
    {
        return dinucNSyn_stat[0][2][1];
    }
    double GetDinucNSyn12_AGAG()
    {
        return dinucNSyn_stat[0][2][2];
    }
    double GetDinucNSyn12_AGAT()
    {
        return dinucNSyn_stat[0][2][3];
    }
    double GetDinucNSyn12_AGCA()
    {
        return dinucNSyn_stat[0][2][4];
    }
    double GetDinucNSyn12_AGCC()
    {
        return dinucNSyn_stat[0][2][5];
    }
    double GetDinucNSyn12_AGCG()
    {
        return dinucNSyn_stat[0][2][6];
    }
    double GetDinucNSyn12_AGCT()
    {
        return dinucNSyn_stat[0][2][7];
    }
    double GetDinucNSyn12_AGGA()
    {
        return dinucNSyn_stat[0][2][8];
    }
    double GetDinucNSyn12_AGGC()
    {
        return dinucNSyn_stat[0][2][9];
    }
    double GetDinucNSyn12_AGGG()
    {
        return dinucNSyn_stat[0][2][10];
    }
    double GetDinucNSyn12_AGGT()
    {
        return dinucNSyn_stat[0][2][11];
    }
    double GetDinucNSyn12_AGTA()
    {
        return dinucNSyn_stat[0][2][12];
    }
    double GetDinucNSyn12_AGTC()
    {
        return dinucNSyn_stat[0][2][13];
    }
    double GetDinucNSyn12_AGTG()
    {
        return dinucNSyn_stat[0][2][14];
    }
    double GetDinucNSyn12_AGTT()
    {
        return dinucNSyn_stat[0][2][15];
    }
    double GetDinucNSyn12_ATAA()
    {
        return dinucNSyn_stat[0][3][0];
    }
    double GetDinucNSyn12_ATAC()
    {
        return dinucNSyn_stat[0][3][1];
    }
    double GetDinucNSyn12_ATAG()
    {
        return dinucNSyn_stat[0][3][2];
    }
    double GetDinucNSyn12_ATAT()
    {
        return dinucNSyn_stat[0][3][3];
    }
    double GetDinucNSyn12_ATCA()
    {
        return dinucNSyn_stat[0][3][4];
    }
    double GetDinucNSyn12_ATCC()
    {
        return dinucNSyn_stat[0][3][5];
    }
    double GetDinucNSyn12_ATCG()
    {
        return dinucNSyn_stat[0][3][6];
    }
    double GetDinucNSyn12_ATCT()
    {
        return dinucNSyn_stat[0][3][7];
    }
    double GetDinucNSyn12_ATGA()
    {
        return dinucNSyn_stat[0][3][8];
    }
    double GetDinucNSyn12_ATGC()
    {
        return dinucNSyn_stat[0][3][9];
    }
    double GetDinucNSyn12_ATGG()
    {
        return dinucNSyn_stat[0][3][10];
    }
    double GetDinucNSyn12_ATGT()
    {
        return dinucNSyn_stat[0][3][11];
    }
    double GetDinucNSyn12_ATTA()
    {
        return dinucNSyn_stat[0][3][12];
    }
    double GetDinucNSyn12_ATTC()
    {
        return dinucNSyn_stat[0][3][13];
    }
    double GetDinucNSyn12_ATTG()
    {
        return dinucNSyn_stat[0][3][14];
    }
    double GetDinucNSyn12_ATTT()
    {
        return dinucNSyn_stat[0][3][15];
    }
    double GetDinucNSyn12_CAAA()
    {
        return dinucNSyn_stat[0][4][0];
    }
    double GetDinucNSyn12_CAAC()
    {
        return dinucNSyn_stat[0][4][1];
    }
    double GetDinucNSyn12_CAAG()
    {
        return dinucNSyn_stat[0][4][2];
    }
    double GetDinucNSyn12_CAAT()
    {
        return dinucNSyn_stat[0][4][3];
    }
    double GetDinucNSyn12_CACA()
    {
        return dinucNSyn_stat[0][4][4];
    }
    double GetDinucNSyn12_CACC()
    {
        return dinucNSyn_stat[0][4][5];
    }
    double GetDinucNSyn12_CACG()
    {
        return dinucNSyn_stat[0][4][6];
    }
    double GetDinucNSyn12_CACT()
    {
        return dinucNSyn_stat[0][4][7];
    }
    double GetDinucNSyn12_CAGA()
    {
        return dinucNSyn_stat[0][4][8];
    }
    double GetDinucNSyn12_CAGC()
    {
        return dinucNSyn_stat[0][4][9];
    }
    double GetDinucNSyn12_CAGG()
    {
        return dinucNSyn_stat[0][4][10];
    }
    double GetDinucNSyn12_CAGT()
    {
        return dinucNSyn_stat[0][4][11];
    }
    double GetDinucNSyn12_CATA()
    {
        return dinucNSyn_stat[0][4][12];
    }
    double GetDinucNSyn12_CATC()
    {
        return dinucNSyn_stat[0][4][13];
    }
    double GetDinucNSyn12_CATG()
    {
        return dinucNSyn_stat[0][4][14];
    }
    double GetDinucNSyn12_CATT()
    {
        return dinucNSyn_stat[0][4][15];
    }
    double GetDinucNSyn12_CCAA()
    {
        return dinucNSyn_stat[0][5][0];
    }
    double GetDinucNSyn12_CCAC()
    {
        return dinucNSyn_stat[0][5][1];
    }
    double GetDinucNSyn12_CCAG()
    {
        return dinucNSyn_stat[0][5][2];
    }
    double GetDinucNSyn12_CCAT()
    {
        return dinucNSyn_stat[0][5][3];
    }
    double GetDinucNSyn12_CCCA()
    {
        return dinucNSyn_stat[0][5][4];
    }
    double GetDinucNSyn12_CCCC()
    {
        return dinucNSyn_stat[0][5][5];
    }
    double GetDinucNSyn12_CCCG()
    {
        return dinucNSyn_stat[0][5][6];
    }
    double GetDinucNSyn12_CCCT()
    {
        return dinucNSyn_stat[0][5][7];
    }
    double GetDinucNSyn12_CCGA()
    {
        return dinucNSyn_stat[0][5][8];
    }
    double GetDinucNSyn12_CCGC()
    {
        return dinucNSyn_stat[0][5][9];
    }
    double GetDinucNSyn12_CCGG()
    {
        return dinucNSyn_stat[0][5][10];
    }
    double GetDinucNSyn12_CCGT()
    {
        return dinucNSyn_stat[0][5][11];
    }
    double GetDinucNSyn12_CCTA()
    {
        return dinucNSyn_stat[0][5][12];
    }
    double GetDinucNSyn12_CCTC()
    {
        return dinucNSyn_stat[0][5][13];
    }
    double GetDinucNSyn12_CCTG()
    {
        return dinucNSyn_stat[0][5][14];
    }
    double GetDinucNSyn12_CCTT()
    {
        return dinucNSyn_stat[0][5][15];
    }
    double GetDinucNSyn12_CGAA()
    {
        return dinucNSyn_stat[0][6][0];
    }
    double GetDinucNSyn12_CGAC()
    {
        return dinucNSyn_stat[0][6][1];
    }
    double GetDinucNSyn12_CGAG()
    {
        return dinucNSyn_stat[0][6][2];
    }
    double GetDinucNSyn12_CGAT()
    {
        return dinucNSyn_stat[0][6][3];
    }
    double GetDinucNSyn12_CGCA()
    {
        return dinucNSyn_stat[0][6][4];
    }
    double GetDinucNSyn12_CGCC()
    {
        return dinucNSyn_stat[0][6][5];
    }
    double GetDinucNSyn12_CGCG()
    {
        return dinucNSyn_stat[0][6][6];
    }
    double GetDinucNSyn12_CGCT()
    {
        return dinucNSyn_stat[0][6][7];
    }
    double GetDinucNSyn12_CGGA()
    {
        return dinucNSyn_stat[0][6][8];
    }
    double GetDinucNSyn12_CGGC()
    {
        return dinucNSyn_stat[0][6][9];
    }
    double GetDinucNSyn12_CGGG()
    {
        return dinucNSyn_stat[0][6][10];
    }
    double GetDinucNSyn12_CGGT()
    {
        return dinucNSyn_stat[0][6][11];
    }
    double GetDinucNSyn12_CGTA()
    {
        return dinucNSyn_stat[0][6][12];
    }
    double GetDinucNSyn12_CGTC()
    {
        return dinucNSyn_stat[0][6][13];
    }
    double GetDinucNSyn12_CGTG()
    {
        return dinucNSyn_stat[0][6][14];
    }
    double GetDinucNSyn12_CGTT()
    {
        return dinucNSyn_stat[0][6][15];
    }
    double GetDinucNSyn12_CTAA()
    {
        return dinucNSyn_stat[0][7][0];
    }
    double GetDinucNSyn12_CTAC()
    {
        return dinucNSyn_stat[0][7][1];
    }
    double GetDinucNSyn12_CTAG()
    {
        return dinucNSyn_stat[0][7][2];
    }
    double GetDinucNSyn12_CTAT()
    {
        return dinucNSyn_stat[0][7][3];
    }
    double GetDinucNSyn12_CTCA()
    {
        return dinucNSyn_stat[0][7][4];
    }
    double GetDinucNSyn12_CTCC()
    {
        return dinucNSyn_stat[0][7][5];
    }
    double GetDinucNSyn12_CTCG()
    {
        return dinucNSyn_stat[0][7][6];
    }
    double GetDinucNSyn12_CTCT()
    {
        return dinucNSyn_stat[0][7][7];
    }
    double GetDinucNSyn12_CTGA()
    {
        return dinucNSyn_stat[0][7][8];
    }
    double GetDinucNSyn12_CTGC()
    {
        return dinucNSyn_stat[0][7][9];
    }
    double GetDinucNSyn12_CTGG()
    {
        return dinucNSyn_stat[0][7][10];
    }
    double GetDinucNSyn12_CTGT()
    {
        return dinucNSyn_stat[0][7][11];
    }
    double GetDinucNSyn12_CTTA()
    {
        return dinucNSyn_stat[0][7][12];
    }
    double GetDinucNSyn12_CTTC()
    {
        return dinucNSyn_stat[0][7][13];
    }
    double GetDinucNSyn12_CTTG()
    {
        return dinucNSyn_stat[0][7][14];
    }
    double GetDinucNSyn12_CTTT()
    {
        return dinucNSyn_stat[0][7][15];
    }
    double GetDinucNSyn12_GAAA()
    {
        return dinucNSyn_stat[0][8][0];
    }
    double GetDinucNSyn12_GAAC()
    {
        return dinucNSyn_stat[0][8][1];
    }
    double GetDinucNSyn12_GAAG()
    {
        return dinucNSyn_stat[0][8][2];
    }
    double GetDinucNSyn12_GAAT()
    {
        return dinucNSyn_stat[0][8][3];
    }
    double GetDinucNSyn12_GACA()
    {
        return dinucNSyn_stat[0][8][4];
    }
    double GetDinucNSyn12_GACC()
    {
        return dinucNSyn_stat[0][8][5];
    }
    double GetDinucNSyn12_GACG()
    {
        return dinucNSyn_stat[0][8][6];
    }
    double GetDinucNSyn12_GACT()
    {
        return dinucNSyn_stat[0][8][7];
    }
    double GetDinucNSyn12_GAGA()
    {
        return dinucNSyn_stat[0][8][8];
    }
    double GetDinucNSyn12_GAGC()
    {
        return dinucNSyn_stat[0][8][9];
    }
    double GetDinucNSyn12_GAGG()
    {
        return dinucNSyn_stat[0][8][10];
    }
    double GetDinucNSyn12_GAGT()
    {
        return dinucNSyn_stat[0][8][11];
    }
    double GetDinucNSyn12_GATA()
    {
        return dinucNSyn_stat[0][8][12];
    }
    double GetDinucNSyn12_GATC()
    {
        return dinucNSyn_stat[0][8][13];
    }
    double GetDinucNSyn12_GATG()
    {
        return dinucNSyn_stat[0][8][14];
    }
    double GetDinucNSyn12_GATT()
    {
        return dinucNSyn_stat[0][8][15];
    }
    double GetDinucNSyn12_GCAA()
    {
        return dinucNSyn_stat[0][9][0];
    }
    double GetDinucNSyn12_GCAC()
    {
        return dinucNSyn_stat[0][9][1];
    }
    double GetDinucNSyn12_GCAG()
    {
        return dinucNSyn_stat[0][9][2];
    }
    double GetDinucNSyn12_GCAT()
    {
        return dinucNSyn_stat[0][9][3];
    }
    double GetDinucNSyn12_GCCA()
    {
        return dinucNSyn_stat[0][9][4];
    }
    double GetDinucNSyn12_GCCC()
    {
        return dinucNSyn_stat[0][9][5];
    }
    double GetDinucNSyn12_GCCG()
    {
        return dinucNSyn_stat[0][9][6];
    }
    double GetDinucNSyn12_GCCT()
    {
        return dinucNSyn_stat[0][9][7];
    }
    double GetDinucNSyn12_GCGA()
    {
        return dinucNSyn_stat[0][9][8];
    }
    double GetDinucNSyn12_GCGC()
    {
        return dinucNSyn_stat[0][9][9];
    }
    double GetDinucNSyn12_GCGG()
    {
        return dinucNSyn_stat[0][9][10];
    }
    double GetDinucNSyn12_GCGT()
    {
        return dinucNSyn_stat[0][9][11];
    }
    double GetDinucNSyn12_GCTA()
    {
        return dinucNSyn_stat[0][9][12];
    }
    double GetDinucNSyn12_GCTC()
    {
        return dinucNSyn_stat[0][9][13];
    }
    double GetDinucNSyn12_GCTG()
    {
        return dinucNSyn_stat[0][9][14];
    }
    double GetDinucNSyn12_GCTT()
    {
        return dinucNSyn_stat[0][9][15];
    }
    double GetDinucNSyn12_GGAA()
    {
        return dinucNSyn_stat[0][10][0];
    }
    double GetDinucNSyn12_GGAC()
    {
        return dinucNSyn_stat[0][10][1];
    }
    double GetDinucNSyn12_GGAG()
    {
        return dinucNSyn_stat[0][10][2];
    }
    double GetDinucNSyn12_GGAT()
    {
        return dinucNSyn_stat[0][10][3];
    }
    double GetDinucNSyn12_GGCA()
    {
        return dinucNSyn_stat[0][10][4];
    }
    double GetDinucNSyn12_GGCC()
    {
        return dinucNSyn_stat[0][10][5];
    }
    double GetDinucNSyn12_GGCG()
    {
        return dinucNSyn_stat[0][10][6];
    }
    double GetDinucNSyn12_GGCT()
    {
        return dinucNSyn_stat[0][10][7];
    }
    double GetDinucNSyn12_GGGA()
    {
        return dinucNSyn_stat[0][10][8];
    }
    double GetDinucNSyn12_GGGC()
    {
        return dinucNSyn_stat[0][10][9];
    }
    double GetDinucNSyn12_GGGG()
    {
        return dinucNSyn_stat[0][10][10];
    }
    double GetDinucNSyn12_GGGT()
    {
        return dinucNSyn_stat[0][10][11];
    }
    double GetDinucNSyn12_GGTA()
    {
        return dinucNSyn_stat[0][10][12];
    }
    double GetDinucNSyn12_GGTC()
    {
        return dinucNSyn_stat[0][10][13];
    }
    double GetDinucNSyn12_GGTG()
    {
        return dinucNSyn_stat[0][10][14];
    }
    double GetDinucNSyn12_GGTT()
    {
        return dinucNSyn_stat[0][10][15];
    }
    double GetDinucNSyn12_GTAA()
    {
        return dinucNSyn_stat[0][11][0];
    }
    double GetDinucNSyn12_GTAC()
    {
        return dinucNSyn_stat[0][11][1];
    }
    double GetDinucNSyn12_GTAG()
    {
        return dinucNSyn_stat[0][11][2];
    }
    double GetDinucNSyn12_GTAT()
    {
        return dinucNSyn_stat[0][11][3];
    }
    double GetDinucNSyn12_GTCA()
    {
        return dinucNSyn_stat[0][11][4];
    }
    double GetDinucNSyn12_GTCC()
    {
        return dinucNSyn_stat[0][11][5];
    }
    double GetDinucNSyn12_GTCG()
    {
        return dinucNSyn_stat[0][11][6];
    }
    double GetDinucNSyn12_GTCT()
    {
        return dinucNSyn_stat[0][11][7];
    }
    double GetDinucNSyn12_GTGA()
    {
        return dinucNSyn_stat[0][11][8];
    }
    double GetDinucNSyn12_GTGC()
    {
        return dinucNSyn_stat[0][11][9];
    }
    double GetDinucNSyn12_GTGG()
    {
        return dinucNSyn_stat[0][11][10];
    }
    double GetDinucNSyn12_GTGT()
    {
        return dinucNSyn_stat[0][11][11];
    }
    double GetDinucNSyn12_GTTA()
    {
        return dinucNSyn_stat[0][11][12];
    }
    double GetDinucNSyn12_GTTC()
    {
        return dinucNSyn_stat[0][11][13];
    }
    double GetDinucNSyn12_GTTG()
    {
        return dinucNSyn_stat[0][11][14];
    }
    double GetDinucNSyn12_GTTT()
    {
        return dinucNSyn_stat[0][11][15];
    }
    double GetDinucNSyn12_TAAA()
    {
        return dinucNSyn_stat[0][12][0];
    }
    double GetDinucNSyn12_TAAC()
    {
        return dinucNSyn_stat[0][12][1];
    }
    double GetDinucNSyn12_TAAG()
    {
        return dinucNSyn_stat[0][12][2];
    }
    double GetDinucNSyn12_TAAT()
    {
        return dinucNSyn_stat[0][12][3];
    }
    double GetDinucNSyn12_TACA()
    {
        return dinucNSyn_stat[0][12][4];
    }
    double GetDinucNSyn12_TACC()
    {
        return dinucNSyn_stat[0][12][5];
    }
    double GetDinucNSyn12_TACG()
    {
        return dinucNSyn_stat[0][12][6];
    }
    double GetDinucNSyn12_TACT()
    {
        return dinucNSyn_stat[0][12][7];
    }
    double GetDinucNSyn12_TAGA()
    {
        return dinucNSyn_stat[0][12][8];
    }
    double GetDinucNSyn12_TAGC()
    {
        return dinucNSyn_stat[0][12][9];
    }
    double GetDinucNSyn12_TAGG()
    {
        return dinucNSyn_stat[0][12][10];
    }
    double GetDinucNSyn12_TAGT()
    {
        return dinucNSyn_stat[0][12][11];
    }
    double GetDinucNSyn12_TATA()
    {
        return dinucNSyn_stat[0][12][12];
    }
    double GetDinucNSyn12_TATC()
    {
        return dinucNSyn_stat[0][12][13];
    }
    double GetDinucNSyn12_TATG()
    {
        return dinucNSyn_stat[0][12][14];
    }
    double GetDinucNSyn12_TATT()
    {
        return dinucNSyn_stat[0][12][15];
    }
    double GetDinucNSyn12_TCAA()
    {
        return dinucNSyn_stat[0][13][0];
    }
    double GetDinucNSyn12_TCAC()
    {
        return dinucNSyn_stat[0][13][1];
    }
    double GetDinucNSyn12_TCAG()
    {
        return dinucNSyn_stat[0][13][2];
    }
    double GetDinucNSyn12_TCAT()
    {
        return dinucNSyn_stat[0][13][3];
    }
    double GetDinucNSyn12_TCCA()
    {
        return dinucNSyn_stat[0][13][4];
    }
    double GetDinucNSyn12_TCCC()
    {
        return dinucNSyn_stat[0][13][5];
    }
    double GetDinucNSyn12_TCCG()
    {
        return dinucNSyn_stat[0][13][6];
    }
    double GetDinucNSyn12_TCCT()
    {
        return dinucNSyn_stat[0][13][7];
    }
    double GetDinucNSyn12_TCGA()
    {
        return dinucNSyn_stat[0][13][8];
    }
    double GetDinucNSyn12_TCGC()
    {
        return dinucNSyn_stat[0][13][9];
    }
    double GetDinucNSyn12_TCGG()
    {
        return dinucNSyn_stat[0][13][10];
    }
    double GetDinucNSyn12_TCGT()
    {
        return dinucNSyn_stat[0][13][11];
    }
    double GetDinucNSyn12_TCTA()
    {
        return dinucNSyn_stat[0][13][12];
    }
    double GetDinucNSyn12_TCTC()
    {
        return dinucNSyn_stat[0][13][13];
    }
    double GetDinucNSyn12_TCTG()
    {
        return dinucNSyn_stat[0][13][14];
    }
    double GetDinucNSyn12_TCTT()
    {
        return dinucNSyn_stat[0][13][15];
    }
    double GetDinucNSyn12_TGAA()
    {
        return dinucNSyn_stat[0][14][0];
    }
    double GetDinucNSyn12_TGAC()
    {
        return dinucNSyn_stat[0][14][1];
    }
    double GetDinucNSyn12_TGAG()
    {
        return dinucNSyn_stat[0][14][2];
    }
    double GetDinucNSyn12_TGAT()
    {
        return dinucNSyn_stat[0][14][3];
    }
    double GetDinucNSyn12_TGCA()
    {
        return dinucNSyn_stat[0][14][4];
    }
    double GetDinucNSyn12_TGCC()
    {
        return dinucNSyn_stat[0][14][5];
    }
    double GetDinucNSyn12_TGCG()
    {
        return dinucNSyn_stat[0][14][6];
    }
    double GetDinucNSyn12_TGCT()
    {
        return dinucNSyn_stat[0][14][7];
    }
    double GetDinucNSyn12_TGGA()
    {
        return dinucNSyn_stat[0][14][8];
    }
    double GetDinucNSyn12_TGGC()
    {
        return dinucNSyn_stat[0][14][9];
    }
    double GetDinucNSyn12_TGGG()
    {
        return dinucNSyn_stat[0][14][10];
    }
    double GetDinucNSyn12_TGGT()
    {
        return dinucNSyn_stat[0][14][11];
    }
    double GetDinucNSyn12_TGTA()
    {
        return dinucNSyn_stat[0][14][12];
    }
    double GetDinucNSyn12_TGTC()
    {
        return dinucNSyn_stat[0][14][13];
    }
    double GetDinucNSyn12_TGTG()
    {
        return dinucNSyn_stat[0][14][14];
    }
    double GetDinucNSyn12_TGTT()
    {
        return dinucNSyn_stat[0][14][15];
    }
    double GetDinucNSyn12_TTAA()
    {
        return dinucNSyn_stat[0][15][0];
    }
    double GetDinucNSyn12_TTAC()
    {
        return dinucNSyn_stat[0][15][1];
    }
    double GetDinucNSyn12_TTAG()
    {
        return dinucNSyn_stat[0][15][2];
    }
    double GetDinucNSyn12_TTAT()
    {
        return dinucNSyn_stat[0][15][3];
    }
    double GetDinucNSyn12_TTCA()
    {
        return dinucNSyn_stat[0][15][4];
    }
    double GetDinucNSyn12_TTCC()
    {
        return dinucNSyn_stat[0][15][5];
    }
    double GetDinucNSyn12_TTCG()
    {
        return dinucNSyn_stat[0][15][6];
    }
    double GetDinucNSyn12_TTCT()
    {
        return dinucNSyn_stat[0][15][7];
    }
    double GetDinucNSyn12_TTGA()
    {
        return dinucNSyn_stat[0][15][8];
    }
    double GetDinucNSyn12_TTGC()
    {
        return dinucNSyn_stat[0][15][9];
    }
    double GetDinucNSyn12_TTGG()
    {
        return dinucNSyn_stat[0][15][10];
    }
    double GetDinucNSyn12_TTGT()
    {
        return dinucNSyn_stat[0][15][11];
    }
    double GetDinucNSyn12_TTTA()
    {
        return dinucNSyn_stat[0][15][12];
    }
    double GetDinucNSyn12_TTTC()
    {
        return dinucNSyn_stat[0][15][13];
    }
    double GetDinucNSyn12_TTTG()
    {
        return dinucNSyn_stat[0][15][14];
    }
    double GetDinucNSyn12_TTTT()
    {
        return dinucNSyn_stat[0][15][15];
    }

    double GetDinucNSyn23_AAAA()
    {
        return dinucNSyn_stat[1][0][0];
    }
    double GetDinucNSyn23_AAAC()
    {
        return dinucNSyn_stat[1][0][1];
    }
    double GetDinucNSyn23_AAAG()
    {
        return dinucNSyn_stat[1][0][2];
    }
    double GetDinucNSyn23_AAAT()
    {
        return dinucNSyn_stat[1][0][3];
    }
    double GetDinucNSyn23_AACA()
    {
        return dinucNSyn_stat[1][0][4];
    }
    double GetDinucNSyn23_AACC()
    {
        return dinucNSyn_stat[1][0][5];
    }
    double GetDinucNSyn23_AACG()
    {
        return dinucNSyn_stat[1][0][6];
    }
    double GetDinucNSyn23_AACT()
    {
        return dinucNSyn_stat[1][0][7];
    }
    double GetDinucNSyn23_AAGA()
    {
        return dinucNSyn_stat[1][0][8];
    }
    double GetDinucNSyn23_AAGC()
    {
        return dinucNSyn_stat[1][0][9];
    }
    double GetDinucNSyn23_AAGG()
    {
        return dinucNSyn_stat[1][0][10];
    }
    double GetDinucNSyn23_AAGT()
    {
        return dinucNSyn_stat[1][0][11];
    }
    double GetDinucNSyn23_AATA()
    {
        return dinucNSyn_stat[1][0][12];
    }
    double GetDinucNSyn23_AATC()
    {
        return dinucNSyn_stat[1][0][13];
    }
    double GetDinucNSyn23_AATG()
    {
        return dinucNSyn_stat[1][0][14];
    }
    double GetDinucNSyn23_AATT()
    {
        return dinucNSyn_stat[1][0][15];
    }
    double GetDinucNSyn23_ACAA()
    {
        return dinucNSyn_stat[1][1][0];
    }
    double GetDinucNSyn23_ACAC()
    {
        return dinucNSyn_stat[1][1][1];
    }
    double GetDinucNSyn23_ACAG()
    {
        return dinucNSyn_stat[1][1][2];
    }
    double GetDinucNSyn23_ACAT()
    {
        return dinucNSyn_stat[1][1][3];
    }
    double GetDinucNSyn23_ACCA()
    {
        return dinucNSyn_stat[1][1][4];
    }
    double GetDinucNSyn23_ACCC()
    {
        return dinucNSyn_stat[1][1][5];
    }
    double GetDinucNSyn23_ACCG()
    {
        return dinucNSyn_stat[1][1][6];
    }
    double GetDinucNSyn23_ACCT()
    {
        return dinucNSyn_stat[1][1][7];
    }
    double GetDinucNSyn23_ACGA()
    {
        return dinucNSyn_stat[1][1][8];
    }
    double GetDinucNSyn23_ACGC()
    {
        return dinucNSyn_stat[1][1][9];
    }
    double GetDinucNSyn23_ACGG()
    {
        return dinucNSyn_stat[1][1][10];
    }
    double GetDinucNSyn23_ACGT()
    {
        return dinucNSyn_stat[1][1][11];
    }
    double GetDinucNSyn23_ACTA()
    {
        return dinucNSyn_stat[1][1][12];
    }
    double GetDinucNSyn23_ACTC()
    {
        return dinucNSyn_stat[1][1][13];
    }
    double GetDinucNSyn23_ACTG()
    {
        return dinucNSyn_stat[1][1][14];
    }
    double GetDinucNSyn23_ACTT()
    {
        return dinucNSyn_stat[1][1][15];
    }
    double GetDinucNSyn23_AGAA()
    {
        return dinucNSyn_stat[1][2][0];
    }
    double GetDinucNSyn23_AGAC()
    {
        return dinucNSyn_stat[1][2][1];
    }
    double GetDinucNSyn23_AGAG()
    {
        return dinucNSyn_stat[1][2][2];
    }
    double GetDinucNSyn23_AGAT()
    {
        return dinucNSyn_stat[1][2][3];
    }
    double GetDinucNSyn23_AGCA()
    {
        return dinucNSyn_stat[1][2][4];
    }
    double GetDinucNSyn23_AGCC()
    {
        return dinucNSyn_stat[1][2][5];
    }
    double GetDinucNSyn23_AGCG()
    {
        return dinucNSyn_stat[1][2][6];
    }
    double GetDinucNSyn23_AGCT()
    {
        return dinucNSyn_stat[1][2][7];
    }
    double GetDinucNSyn23_AGGA()
    {
        return dinucNSyn_stat[1][2][8];
    }
    double GetDinucNSyn23_AGGC()
    {
        return dinucNSyn_stat[1][2][9];
    }
    double GetDinucNSyn23_AGGG()
    {
        return dinucNSyn_stat[1][2][10];
    }
    double GetDinucNSyn23_AGGT()
    {
        return dinucNSyn_stat[1][2][11];
    }
    double GetDinucNSyn23_AGTA()
    {
        return dinucNSyn_stat[1][2][12];
    }
    double GetDinucNSyn23_AGTC()
    {
        return dinucNSyn_stat[1][2][13];
    }
    double GetDinucNSyn23_AGTG()
    {
        return dinucNSyn_stat[1][2][14];
    }
    double GetDinucNSyn23_AGTT()
    {
        return dinucNSyn_stat[1][2][15];
    }
    double GetDinucNSyn23_ATAA()
    {
        return dinucNSyn_stat[1][3][0];
    }
    double GetDinucNSyn23_ATAC()
    {
        return dinucNSyn_stat[1][3][1];
    }
    double GetDinucNSyn23_ATAG()
    {
        return dinucNSyn_stat[1][3][2];
    }
    double GetDinucNSyn23_ATAT()
    {
        return dinucNSyn_stat[1][3][3];
    }
    double GetDinucNSyn23_ATCA()
    {
        return dinucNSyn_stat[1][3][4];
    }
    double GetDinucNSyn23_ATCC()
    {
        return dinucNSyn_stat[1][3][5];
    }
    double GetDinucNSyn23_ATCG()
    {
        return dinucNSyn_stat[1][3][6];
    }
    double GetDinucNSyn23_ATCT()
    {
        return dinucNSyn_stat[1][3][7];
    }
    double GetDinucNSyn23_ATGA()
    {
        return dinucNSyn_stat[1][3][8];
    }
    double GetDinucNSyn23_ATGC()
    {
        return dinucNSyn_stat[1][3][9];
    }
    double GetDinucNSyn23_ATGG()
    {
        return dinucNSyn_stat[1][3][10];
    }
    double GetDinucNSyn23_ATGT()
    {
        return dinucNSyn_stat[1][3][11];
    }
    double GetDinucNSyn23_ATTA()
    {
        return dinucNSyn_stat[1][3][12];
    }
    double GetDinucNSyn23_ATTC()
    {
        return dinucNSyn_stat[1][3][13];
    }
    double GetDinucNSyn23_ATTG()
    {
        return dinucNSyn_stat[1][3][14];
    }
    double GetDinucNSyn23_ATTT()
    {
        return dinucNSyn_stat[1][3][15];
    }
    double GetDinucNSyn23_CAAA()
    {
        return dinucNSyn_stat[1][4][0];
    }
    double GetDinucNSyn23_CAAC()
    {
        return dinucNSyn_stat[1][4][1];
    }
    double GetDinucNSyn23_CAAG()
    {
        return dinucNSyn_stat[1][4][2];
    }
    double GetDinucNSyn23_CAAT()
    {
        return dinucNSyn_stat[1][4][3];
    }
    double GetDinucNSyn23_CACA()
    {
        return dinucNSyn_stat[1][4][4];
    }
    double GetDinucNSyn23_CACC()
    {
        return dinucNSyn_stat[1][4][5];
    }
    double GetDinucNSyn23_CACG()
    {
        return dinucNSyn_stat[1][4][6];
    }
    double GetDinucNSyn23_CACT()
    {
        return dinucNSyn_stat[1][4][7];
    }
    double GetDinucNSyn23_CAGA()
    {
        return dinucNSyn_stat[1][4][8];
    }
    double GetDinucNSyn23_CAGC()
    {
        return dinucNSyn_stat[1][4][9];
    }
    double GetDinucNSyn23_CAGG()
    {
        return dinucNSyn_stat[1][4][10];
    }
    double GetDinucNSyn23_CAGT()
    {
        return dinucNSyn_stat[1][4][11];
    }
    double GetDinucNSyn23_CATA()
    {
        return dinucNSyn_stat[1][4][12];
    }
    double GetDinucNSyn23_CATC()
    {
        return dinucNSyn_stat[1][4][13];
    }
    double GetDinucNSyn23_CATG()
    {
        return dinucNSyn_stat[1][4][14];
    }
    double GetDinucNSyn23_CATT()
    {
        return dinucNSyn_stat[1][4][15];
    }
    double GetDinucNSyn23_CCAA()
    {
        return dinucNSyn_stat[1][5][0];
    }
    double GetDinucNSyn23_CCAC()
    {
        return dinucNSyn_stat[1][5][1];
    }
    double GetDinucNSyn23_CCAG()
    {
        return dinucNSyn_stat[1][5][2];
    }
    double GetDinucNSyn23_CCAT()
    {
        return dinucNSyn_stat[1][5][3];
    }
    double GetDinucNSyn23_CCCA()
    {
        return dinucNSyn_stat[1][5][4];
    }
    double GetDinucNSyn23_CCCC()
    {
        return dinucNSyn_stat[1][5][5];
    }
    double GetDinucNSyn23_CCCG()
    {
        return dinucNSyn_stat[1][5][6];
    }
    double GetDinucNSyn23_CCCT()
    {
        return dinucNSyn_stat[1][5][7];
    }
    double GetDinucNSyn23_CCGA()
    {
        return dinucNSyn_stat[1][5][8];
    }
    double GetDinucNSyn23_CCGC()
    {
        return dinucNSyn_stat[1][5][9];
    }
    double GetDinucNSyn23_CCGG()
    {
        return dinucNSyn_stat[1][5][10];
    }
    double GetDinucNSyn23_CCGT()
    {
        return dinucNSyn_stat[1][5][11];
    }
    double GetDinucNSyn23_CCTA()
    {
        return dinucNSyn_stat[1][5][12];
    }
    double GetDinucNSyn23_CCTC()
    {
        return dinucNSyn_stat[1][5][13];
    }
    double GetDinucNSyn23_CCTG()
    {
        return dinucNSyn_stat[1][5][14];
    }
    double GetDinucNSyn23_CCTT()
    {
        return dinucNSyn_stat[1][5][15];
    }
    double GetDinucNSyn23_CGAA()
    {
        return dinucNSyn_stat[1][6][0];
    }
    double GetDinucNSyn23_CGAC()
    {
        return dinucNSyn_stat[1][6][1];
    }
    double GetDinucNSyn23_CGAG()
    {
        return dinucNSyn_stat[1][6][2];
    }
    double GetDinucNSyn23_CGAT()
    {
        return dinucNSyn_stat[1][6][3];
    }
    double GetDinucNSyn23_CGCA()
    {
        return dinucNSyn_stat[1][6][4];
    }
    double GetDinucNSyn23_CGCC()
    {
        return dinucNSyn_stat[1][6][5];
    }
    double GetDinucNSyn23_CGCG()
    {
        return dinucNSyn_stat[1][6][6];
    }
    double GetDinucNSyn23_CGCT()
    {
        return dinucNSyn_stat[1][6][7];
    }
    double GetDinucNSyn23_CGGA()
    {
        return dinucNSyn_stat[1][6][8];
    }
    double GetDinucNSyn23_CGGC()
    {
        return dinucNSyn_stat[1][6][9];
    }
    double GetDinucNSyn23_CGGG()
    {
        return dinucNSyn_stat[1][6][10];
    }
    double GetDinucNSyn23_CGGT()
    {
        return dinucNSyn_stat[1][6][11];
    }
    double GetDinucNSyn23_CGTA()
    {
        return dinucNSyn_stat[1][6][12];
    }
    double GetDinucNSyn23_CGTC()
    {
        return dinucNSyn_stat[1][6][13];
    }
    double GetDinucNSyn23_CGTG()
    {
        return dinucNSyn_stat[1][6][14];
    }
    double GetDinucNSyn23_CGTT()
    {
        return dinucNSyn_stat[1][6][15];
    }
    double GetDinucNSyn23_CTAA()
    {
        return dinucNSyn_stat[1][7][0];
    }
    double GetDinucNSyn23_CTAC()
    {
        return dinucNSyn_stat[1][7][1];
    }
    double GetDinucNSyn23_CTAG()
    {
        return dinucNSyn_stat[1][7][2];
    }
    double GetDinucNSyn23_CTAT()
    {
        return dinucNSyn_stat[1][7][3];
    }
    double GetDinucNSyn23_CTCA()
    {
        return dinucNSyn_stat[1][7][4];
    }
    double GetDinucNSyn23_CTCC()
    {
        return dinucNSyn_stat[1][7][5];
    }
    double GetDinucNSyn23_CTCG()
    {
        return dinucNSyn_stat[1][7][6];
    }
    double GetDinucNSyn23_CTCT()
    {
        return dinucNSyn_stat[1][7][7];
    }
    double GetDinucNSyn23_CTGA()
    {
        return dinucNSyn_stat[1][7][8];
    }
    double GetDinucNSyn23_CTGC()
    {
        return dinucNSyn_stat[1][7][9];
    }
    double GetDinucNSyn23_CTGG()
    {
        return dinucNSyn_stat[1][7][10];
    }
    double GetDinucNSyn23_CTGT()
    {
        return dinucNSyn_stat[1][7][11];
    }
    double GetDinucNSyn23_CTTA()
    {
        return dinucNSyn_stat[1][7][12];
    }
    double GetDinucNSyn23_CTTC()
    {
        return dinucNSyn_stat[1][7][13];
    }
    double GetDinucNSyn23_CTTG()
    {
        return dinucNSyn_stat[1][7][14];
    }
    double GetDinucNSyn23_CTTT()
    {
        return dinucNSyn_stat[1][7][15];
    }
    double GetDinucNSyn23_GAAA()
    {
        return dinucNSyn_stat[1][8][0];
    }
    double GetDinucNSyn23_GAAC()
    {
        return dinucNSyn_stat[1][8][1];
    }
    double GetDinucNSyn23_GAAG()
    {
        return dinucNSyn_stat[1][8][2];
    }
    double GetDinucNSyn23_GAAT()
    {
        return dinucNSyn_stat[1][8][3];
    }
    double GetDinucNSyn23_GACA()
    {
        return dinucNSyn_stat[1][8][4];
    }
    double GetDinucNSyn23_GACC()
    {
        return dinucNSyn_stat[1][8][5];
    }
    double GetDinucNSyn23_GACG()
    {
        return dinucNSyn_stat[1][8][6];
    }
    double GetDinucNSyn23_GACT()
    {
        return dinucNSyn_stat[1][8][7];
    }
    double GetDinucNSyn23_GAGA()
    {
        return dinucNSyn_stat[1][8][8];
    }
    double GetDinucNSyn23_GAGC()
    {
        return dinucNSyn_stat[1][8][9];
    }
    double GetDinucNSyn23_GAGG()
    {
        return dinucNSyn_stat[1][8][10];
    }
    double GetDinucNSyn23_GAGT()
    {
        return dinucNSyn_stat[1][8][11];
    }
    double GetDinucNSyn23_GATA()
    {
        return dinucNSyn_stat[1][8][12];
    }
    double GetDinucNSyn23_GATC()
    {
        return dinucNSyn_stat[1][8][13];
    }
    double GetDinucNSyn23_GATG()
    {
        return dinucNSyn_stat[1][8][14];
    }
    double GetDinucNSyn23_GATT()
    {
        return dinucNSyn_stat[1][8][15];
    }
    double GetDinucNSyn23_GCAA()
    {
        return dinucNSyn_stat[1][9][0];
    }
    double GetDinucNSyn23_GCAC()
    {
        return dinucNSyn_stat[1][9][1];
    }
    double GetDinucNSyn23_GCAG()
    {
        return dinucNSyn_stat[1][9][2];
    }
    double GetDinucNSyn23_GCAT()
    {
        return dinucNSyn_stat[1][9][3];
    }
    double GetDinucNSyn23_GCCA()
    {
        return dinucNSyn_stat[1][9][4];
    }
    double GetDinucNSyn23_GCCC()
    {
        return dinucNSyn_stat[1][9][5];
    }
    double GetDinucNSyn23_GCCG()
    {
        return dinucNSyn_stat[1][9][6];
    }
    double GetDinucNSyn23_GCCT()
    {
        return dinucNSyn_stat[1][9][7];
    }
    double GetDinucNSyn23_GCGA()
    {
        return dinucNSyn_stat[1][9][8];
    }
    double GetDinucNSyn23_GCGC()
    {
        return dinucNSyn_stat[1][9][9];
    }
    double GetDinucNSyn23_GCGG()
    {
        return dinucNSyn_stat[1][9][10];
    }
    double GetDinucNSyn23_GCGT()
    {
        return dinucNSyn_stat[1][9][11];
    }
    double GetDinucNSyn23_GCTA()
    {
        return dinucNSyn_stat[1][9][12];
    }
    double GetDinucNSyn23_GCTC()
    {
        return dinucNSyn_stat[1][9][13];
    }
    double GetDinucNSyn23_GCTG()
    {
        return dinucNSyn_stat[1][9][14];
    }
    double GetDinucNSyn23_GCTT()
    {
        return dinucNSyn_stat[1][9][15];
    }
    double GetDinucNSyn23_GGAA()
    {
        return dinucNSyn_stat[1][10][0];
    }
    double GetDinucNSyn23_GGAC()
    {
        return dinucNSyn_stat[1][10][1];
    }
    double GetDinucNSyn23_GGAG()
    {
        return dinucNSyn_stat[1][10][2];
    }
    double GetDinucNSyn23_GGAT()
    {
        return dinucNSyn_stat[1][10][3];
    }
    double GetDinucNSyn23_GGCA()
    {
        return dinucNSyn_stat[1][10][4];
    }
    double GetDinucNSyn23_GGCC()
    {
        return dinucNSyn_stat[1][10][5];
    }
    double GetDinucNSyn23_GGCG()
    {
        return dinucNSyn_stat[1][10][6];
    }
    double GetDinucNSyn23_GGCT()
    {
        return dinucNSyn_stat[1][10][7];
    }
    double GetDinucNSyn23_GGGA()
    {
        return dinucNSyn_stat[1][10][8];
    }
    double GetDinucNSyn23_GGGC()
    {
        return dinucNSyn_stat[1][10][9];
    }
    double GetDinucNSyn23_GGGG()
    {
        return dinucNSyn_stat[1][10][10];
    }
    double GetDinucNSyn23_GGGT()
    {
        return dinucNSyn_stat[1][10][11];
    }
    double GetDinucNSyn23_GGTA()
    {
        return dinucNSyn_stat[1][10][12];
    }
    double GetDinucNSyn23_GGTC()
    {
        return dinucNSyn_stat[1][10][13];
    }
    double GetDinucNSyn23_GGTG()
    {
        return dinucNSyn_stat[1][10][14];
    }
    double GetDinucNSyn23_GGTT()
    {
        return dinucNSyn_stat[1][10][15];
    }
    double GetDinucNSyn23_GTAA()
    {
        return dinucNSyn_stat[1][11][0];
    }
    double GetDinucNSyn23_GTAC()
    {
        return dinucNSyn_stat[1][11][1];
    }
    double GetDinucNSyn23_GTAG()
    {
        return dinucNSyn_stat[1][11][2];
    }
    double GetDinucNSyn23_GTAT()
    {
        return dinucNSyn_stat[1][11][3];
    }
    double GetDinucNSyn23_GTCA()
    {
        return dinucNSyn_stat[1][11][4];
    }
    double GetDinucNSyn23_GTCC()
    {
        return dinucNSyn_stat[1][11][5];
    }
    double GetDinucNSyn23_GTCG()
    {
        return dinucNSyn_stat[1][11][6];
    }
    double GetDinucNSyn23_GTCT()
    {
        return dinucNSyn_stat[1][11][7];
    }
    double GetDinucNSyn23_GTGA()
    {
        return dinucNSyn_stat[1][11][8];
    }
    double GetDinucNSyn23_GTGC()
    {
        return dinucNSyn_stat[1][11][9];
    }
    double GetDinucNSyn23_GTGG()
    {
        return dinucNSyn_stat[1][11][10];
    }
    double GetDinucNSyn23_GTGT()
    {
        return dinucNSyn_stat[1][11][11];
    }
    double GetDinucNSyn23_GTTA()
    {
        return dinucNSyn_stat[1][11][12];
    }
    double GetDinucNSyn23_GTTC()
    {
        return dinucNSyn_stat[1][11][13];
    }
    double GetDinucNSyn23_GTTG()
    {
        return dinucNSyn_stat[1][11][14];
    }
    double GetDinucNSyn23_GTTT()
    {
        return dinucNSyn_stat[1][11][15];
    }
    double GetDinucNSyn23_TAAA()
    {
        return dinucNSyn_stat[1][12][0];
    }
    double GetDinucNSyn23_TAAC()
    {
        return dinucNSyn_stat[1][12][1];
    }
    double GetDinucNSyn23_TAAG()
    {
        return dinucNSyn_stat[1][12][2];
    }
    double GetDinucNSyn23_TAAT()
    {
        return dinucNSyn_stat[1][12][3];
    }
    double GetDinucNSyn23_TACA()
    {
        return dinucNSyn_stat[1][12][4];
    }
    double GetDinucNSyn23_TACC()
    {
        return dinucNSyn_stat[1][12][5];
    }
    double GetDinucNSyn23_TACG()
    {
        return dinucNSyn_stat[1][12][6];
    }
    double GetDinucNSyn23_TACT()
    {
        return dinucNSyn_stat[1][12][7];
    }
    double GetDinucNSyn23_TAGA()
    {
        return dinucNSyn_stat[1][12][8];
    }
    double GetDinucNSyn23_TAGC()
    {
        return dinucNSyn_stat[1][12][9];
    }
    double GetDinucNSyn23_TAGG()
    {
        return dinucNSyn_stat[1][12][10];
    }
    double GetDinucNSyn23_TAGT()
    {
        return dinucNSyn_stat[1][12][11];
    }
    double GetDinucNSyn23_TATA()
    {
        return dinucNSyn_stat[1][12][12];
    }
    double GetDinucNSyn23_TATC()
    {
        return dinucNSyn_stat[1][12][13];
    }
    double GetDinucNSyn23_TATG()
    {
        return dinucNSyn_stat[1][12][14];
    }
    double GetDinucNSyn23_TATT()
    {
        return dinucNSyn_stat[1][12][15];
    }
    double GetDinucNSyn23_TCAA()
    {
        return dinucNSyn_stat[1][13][0];
    }
    double GetDinucNSyn23_TCAC()
    {
        return dinucNSyn_stat[1][13][1];
    }
    double GetDinucNSyn23_TCAG()
    {
        return dinucNSyn_stat[1][13][2];
    }
    double GetDinucNSyn23_TCAT()
    {
        return dinucNSyn_stat[1][13][3];
    }
    double GetDinucNSyn23_TCCA()
    {
        return dinucNSyn_stat[1][13][4];
    }
    double GetDinucNSyn23_TCCC()
    {
        return dinucNSyn_stat[1][13][5];
    }
    double GetDinucNSyn23_TCCG()
    {
        return dinucNSyn_stat[1][13][6];
    }
    double GetDinucNSyn23_TCCT()
    {
        return dinucNSyn_stat[1][13][7];
    }
    double GetDinucNSyn23_TCGA()
    {
        return dinucNSyn_stat[1][13][8];
    }
    double GetDinucNSyn23_TCGC()
    {
        return dinucNSyn_stat[1][13][9];
    }
    double GetDinucNSyn23_TCGG()
    {
        return dinucNSyn_stat[1][13][10];
    }
    double GetDinucNSyn23_TCGT()
    {
        return dinucNSyn_stat[1][13][11];
    }
    double GetDinucNSyn23_TCTA()
    {
        return dinucNSyn_stat[1][13][12];
    }
    double GetDinucNSyn23_TCTC()
    {
        return dinucNSyn_stat[1][13][13];
    }
    double GetDinucNSyn23_TCTG()
    {
        return dinucNSyn_stat[1][13][14];
    }
    double GetDinucNSyn23_TCTT()
    {
        return dinucNSyn_stat[1][13][15];
    }
    double GetDinucNSyn23_TGAA()
    {
        return dinucNSyn_stat[1][14][0];
    }
    double GetDinucNSyn23_TGAC()
    {
        return dinucNSyn_stat[1][14][1];
    }
    double GetDinucNSyn23_TGAG()
    {
        return dinucNSyn_stat[1][14][2];
    }
    double GetDinucNSyn23_TGAT()
    {
        return dinucNSyn_stat[1][14][3];
    }
    double GetDinucNSyn23_TGCA()
    {
        return dinucNSyn_stat[1][14][4];
    }
    double GetDinucNSyn23_TGCC()
    {
        return dinucNSyn_stat[1][14][5];
    }
    double GetDinucNSyn23_TGCG()
    {
        return dinucNSyn_stat[1][14][6];
    }
    double GetDinucNSyn23_TGCT()
    {
        return dinucNSyn_stat[1][14][7];
    }
    double GetDinucNSyn23_TGGA()
    {
        return dinucNSyn_stat[1][14][8];
    }
    double GetDinucNSyn23_TGGC()
    {
        return dinucNSyn_stat[1][14][9];
    }
    double GetDinucNSyn23_TGGG()
    {
        return dinucNSyn_stat[1][14][10];
    }
    double GetDinucNSyn23_TGGT()
    {
        return dinucNSyn_stat[1][14][11];
    }
    double GetDinucNSyn23_TGTA()
    {
        return dinucNSyn_stat[1][14][12];
    }
    double GetDinucNSyn23_TGTC()
    {
        return dinucNSyn_stat[1][14][13];
    }
    double GetDinucNSyn23_TGTG()
    {
        return dinucNSyn_stat[1][14][14];
    }
    double GetDinucNSyn23_TGTT()
    {
        return dinucNSyn_stat[1][14][15];
    }
    double GetDinucNSyn23_TTAA()
    {
        return dinucNSyn_stat[1][15][0];
    }
    double GetDinucNSyn23_TTAC()
    {
        return dinucNSyn_stat[1][15][1];
    }
    double GetDinucNSyn23_TTAG()
    {
        return dinucNSyn_stat[1][15][2];
    }
    double GetDinucNSyn23_TTAT()
    {
        return dinucNSyn_stat[1][15][3];
    }
    double GetDinucNSyn23_TTCA()
    {
        return dinucNSyn_stat[1][15][4];
    }
    double GetDinucNSyn23_TTCC()
    {
        return dinucNSyn_stat[1][15][5];
    }
    double GetDinucNSyn23_TTCG()
    {
        return dinucNSyn_stat[1][15][6];
    }
    double GetDinucNSyn23_TTCT()
    {
        return dinucNSyn_stat[1][15][7];
    }
    double GetDinucNSyn23_TTGA()
    {
        return dinucNSyn_stat[1][15][8];
    }
    double GetDinucNSyn23_TTGC()
    {
        return dinucNSyn_stat[1][15][9];
    }
    double GetDinucNSyn23_TTGG()
    {
        return dinucNSyn_stat[1][15][10];
    }
    double GetDinucNSyn23_TTGT()
    {
        return dinucNSyn_stat[1][15][11];
    }
    double GetDinucNSyn23_TTTA()
    {
        return dinucNSyn_stat[1][15][12];
    }
    double GetDinucNSyn23_TTTC()
    {
        return dinucNSyn_stat[1][15][13];
    }
    double GetDinucNSyn23_TTTG()
    {
        return dinucNSyn_stat[1][15][14];
    }
    double GetDinucNSyn23_TTTT()
    {
        return dinucNSyn_stat[1][15][15];
    }

    double GetDinucNSyn31_AAAA()
    {
        return dinucNSyn_stat[2][0][0];
    }
    double GetDinucNSyn31_AAAC()
    {
        return dinucNSyn_stat[2][0][1];
    }
    double GetDinucNSyn31_AAAG()
    {
        return dinucNSyn_stat[2][0][2];
    }
    double GetDinucNSyn31_AAAT()
    {
        return dinucNSyn_stat[2][0][3];
    }
    double GetDinucNSyn31_AACA()
    {
        return dinucNSyn_stat[2][0][4];
    }
    double GetDinucNSyn31_AACC()
    {
        return dinucNSyn_stat[2][0][5];
    }
    double GetDinucNSyn31_AACG()
    {
        return dinucNSyn_stat[2][0][6];
    }
    double GetDinucNSyn31_AACT()
    {
        return dinucNSyn_stat[2][0][7];
    }
    double GetDinucNSyn31_AAGA()
    {
        return dinucNSyn_stat[2][0][8];
    }
    double GetDinucNSyn31_AAGC()
    {
        return dinucNSyn_stat[2][0][9];
    }
    double GetDinucNSyn31_AAGG()
    {
        return dinucNSyn_stat[2][0][10];
    }
    double GetDinucNSyn31_AAGT()
    {
        return dinucNSyn_stat[2][0][11];
    }
    double GetDinucNSyn31_AATA()
    {
        return dinucNSyn_stat[2][0][12];
    }
    double GetDinucNSyn31_AATC()
    {
        return dinucNSyn_stat[2][0][13];
    }
    double GetDinucNSyn31_AATG()
    {
        return dinucNSyn_stat[2][0][14];
    }
    double GetDinucNSyn31_AATT()
    {
        return dinucNSyn_stat[2][0][15];
    }
    double GetDinucNSyn31_ACAA()
    {
        return dinucNSyn_stat[2][1][0];
    }
    double GetDinucNSyn31_ACAC()
    {
        return dinucNSyn_stat[2][1][1];
    }
    double GetDinucNSyn31_ACAG()
    {
        return dinucNSyn_stat[2][1][2];
    }
    double GetDinucNSyn31_ACAT()
    {
        return dinucNSyn_stat[2][1][3];
    }
    double GetDinucNSyn31_ACCA()
    {
        return dinucNSyn_stat[2][1][4];
    }
    double GetDinucNSyn31_ACCC()
    {
        return dinucNSyn_stat[2][1][5];
    }
    double GetDinucNSyn31_ACCG()
    {
        return dinucNSyn_stat[2][1][6];
    }
    double GetDinucNSyn31_ACCT()
    {
        return dinucNSyn_stat[2][1][7];
    }
    double GetDinucNSyn31_ACGA()
    {
        return dinucNSyn_stat[2][1][8];
    }
    double GetDinucNSyn31_ACGC()
    {
        return dinucNSyn_stat[2][1][9];
    }
    double GetDinucNSyn31_ACGG()
    {
        return dinucNSyn_stat[2][1][10];
    }
    double GetDinucNSyn31_ACGT()
    {
        return dinucNSyn_stat[2][1][11];
    }
    double GetDinucNSyn31_ACTA()
    {
        return dinucNSyn_stat[2][1][12];
    }
    double GetDinucNSyn31_ACTC()
    {
        return dinucNSyn_stat[2][1][13];
    }
    double GetDinucNSyn31_ACTG()
    {
        return dinucNSyn_stat[2][1][14];
    }
    double GetDinucNSyn31_ACTT()
    {
        return dinucNSyn_stat[2][1][15];
    }
    double GetDinucNSyn31_AGAA()
    {
        return dinucNSyn_stat[2][2][0];
    }
    double GetDinucNSyn31_AGAC()
    {
        return dinucNSyn_stat[2][2][1];
    }
    double GetDinucNSyn31_AGAG()
    {
        return dinucNSyn_stat[2][2][2];
    }
    double GetDinucNSyn31_AGAT()
    {
        return dinucNSyn_stat[2][2][3];
    }
    double GetDinucNSyn31_AGCA()
    {
        return dinucNSyn_stat[2][2][4];
    }
    double GetDinucNSyn31_AGCC()
    {
        return dinucNSyn_stat[2][2][5];
    }
    double GetDinucNSyn31_AGCG()
    {
        return dinucNSyn_stat[2][2][6];
    }
    double GetDinucNSyn31_AGCT()
    {
        return dinucNSyn_stat[2][2][7];
    }
    double GetDinucNSyn31_AGGA()
    {
        return dinucNSyn_stat[2][2][8];
    }
    double GetDinucNSyn31_AGGC()
    {
        return dinucNSyn_stat[2][2][9];
    }
    double GetDinucNSyn31_AGGG()
    {
        return dinucNSyn_stat[2][2][10];
    }
    double GetDinucNSyn31_AGGT()
    {
        return dinucNSyn_stat[2][2][11];
    }
    double GetDinucNSyn31_AGTA()
    {
        return dinucNSyn_stat[2][2][12];
    }
    double GetDinucNSyn31_AGTC()
    {
        return dinucNSyn_stat[2][2][13];
    }
    double GetDinucNSyn31_AGTG()
    {
        return dinucNSyn_stat[2][2][14];
    }
    double GetDinucNSyn31_AGTT()
    {
        return dinucNSyn_stat[2][2][15];
    }
    double GetDinucNSyn31_ATAA()
    {
        return dinucNSyn_stat[2][3][0];
    }
    double GetDinucNSyn31_ATAC()
    {
        return dinucNSyn_stat[2][3][1];
    }
    double GetDinucNSyn31_ATAG()
    {
        return dinucNSyn_stat[2][3][2];
    }
    double GetDinucNSyn31_ATAT()
    {
        return dinucNSyn_stat[2][3][3];
    }
    double GetDinucNSyn31_ATCA()
    {
        return dinucNSyn_stat[2][3][4];
    }
    double GetDinucNSyn31_ATCC()
    {
        return dinucNSyn_stat[2][3][5];
    }
    double GetDinucNSyn31_ATCG()
    {
        return dinucNSyn_stat[2][3][6];
    }
    double GetDinucNSyn31_ATCT()
    {
        return dinucNSyn_stat[2][3][7];
    }
    double GetDinucNSyn31_ATGA()
    {
        return dinucNSyn_stat[2][3][8];
    }
    double GetDinucNSyn31_ATGC()
    {
        return dinucNSyn_stat[2][3][9];
    }
    double GetDinucNSyn31_ATGG()
    {
        return dinucNSyn_stat[2][3][10];
    }
    double GetDinucNSyn31_ATGT()
    {
        return dinucNSyn_stat[2][3][11];
    }
    double GetDinucNSyn31_ATTA()
    {
        return dinucNSyn_stat[2][3][12];
    }
    double GetDinucNSyn31_ATTC()
    {
        return dinucNSyn_stat[2][3][13];
    }
    double GetDinucNSyn31_ATTG()
    {
        return dinucNSyn_stat[2][3][14];
    }
    double GetDinucNSyn31_ATTT()
    {
        return dinucNSyn_stat[2][3][15];
    }
    double GetDinucNSyn31_CAAA()
    {
        return dinucNSyn_stat[2][4][0];
    }
    double GetDinucNSyn31_CAAC()
    {
        return dinucNSyn_stat[2][4][1];
    }
    double GetDinucNSyn31_CAAG()
    {
        return dinucNSyn_stat[2][4][2];
    }
    double GetDinucNSyn31_CAAT()
    {
        return dinucNSyn_stat[2][4][3];
    }
    double GetDinucNSyn31_CACA()
    {
        return dinucNSyn_stat[2][4][4];
    }
    double GetDinucNSyn31_CACC()
    {
        return dinucNSyn_stat[2][4][5];
    }
    double GetDinucNSyn31_CACG()
    {
        return dinucNSyn_stat[2][4][6];
    }
    double GetDinucNSyn31_CACT()
    {
        return dinucNSyn_stat[2][4][7];
    }
    double GetDinucNSyn31_CAGA()
    {
        return dinucNSyn_stat[2][4][8];
    }
    double GetDinucNSyn31_CAGC()
    {
        return dinucNSyn_stat[2][4][9];
    }
    double GetDinucNSyn31_CAGG()
    {
        return dinucNSyn_stat[2][4][10];
    }
    double GetDinucNSyn31_CAGT()
    {
        return dinucNSyn_stat[2][4][11];
    }
    double GetDinucNSyn31_CATA()
    {
        return dinucNSyn_stat[2][4][12];
    }
    double GetDinucNSyn31_CATC()
    {
        return dinucNSyn_stat[2][4][13];
    }
    double GetDinucNSyn31_CATG()
    {
        return dinucNSyn_stat[2][4][14];
    }
    double GetDinucNSyn31_CATT()
    {
        return dinucNSyn_stat[2][4][15];
    }
    double GetDinucNSyn31_CCAA()
    {
        return dinucNSyn_stat[2][5][0];
    }
    double GetDinucNSyn31_CCAC()
    {
        return dinucNSyn_stat[2][5][1];
    }
    double GetDinucNSyn31_CCAG()
    {
        return dinucNSyn_stat[2][5][2];
    }
    double GetDinucNSyn31_CCAT()
    {
        return dinucNSyn_stat[2][5][3];
    }
    double GetDinucNSyn31_CCCA()
    {
        return dinucNSyn_stat[2][5][4];
    }
    double GetDinucNSyn31_CCCC()
    {
        return dinucNSyn_stat[2][5][5];
    }
    double GetDinucNSyn31_CCCG()
    {
        return dinucNSyn_stat[2][5][6];
    }
    double GetDinucNSyn31_CCCT()
    {
        return dinucNSyn_stat[2][5][7];
    }
    double GetDinucNSyn31_CCGA()
    {
        return dinucNSyn_stat[2][5][8];
    }
    double GetDinucNSyn31_CCGC()
    {
        return dinucNSyn_stat[2][5][9];
    }
    double GetDinucNSyn31_CCGG()
    {
        return dinucNSyn_stat[2][5][10];
    }
    double GetDinucNSyn31_CCGT()
    {
        return dinucNSyn_stat[2][5][11];
    }
    double GetDinucNSyn31_CCTA()
    {
        return dinucNSyn_stat[2][5][12];
    }
    double GetDinucNSyn31_CCTC()
    {
        return dinucNSyn_stat[2][5][13];
    }
    double GetDinucNSyn31_CCTG()
    {
        return dinucNSyn_stat[2][5][14];
    }
    double GetDinucNSyn31_CCTT()
    {
        return dinucNSyn_stat[2][5][15];
    }
    double GetDinucNSyn31_CGAA()
    {
        return dinucNSyn_stat[2][6][0];
    }
    double GetDinucNSyn31_CGAC()
    {
        return dinucNSyn_stat[2][6][1];
    }
    double GetDinucNSyn31_CGAG()
    {
        return dinucNSyn_stat[2][6][2];
    }
    double GetDinucNSyn31_CGAT()
    {
        return dinucNSyn_stat[2][6][3];
    }
    double GetDinucNSyn31_CGCA()
    {
        return dinucNSyn_stat[2][6][4];
    }
    double GetDinucNSyn31_CGCC()
    {
        return dinucNSyn_stat[2][6][5];
    }
    double GetDinucNSyn31_CGCG()
    {
        return dinucNSyn_stat[2][6][6];
    }
    double GetDinucNSyn31_CGCT()
    {
        return dinucNSyn_stat[2][6][7];
    }
    double GetDinucNSyn31_CGGA()
    {
        return dinucNSyn_stat[2][6][8];
    }
    double GetDinucNSyn31_CGGC()
    {
        return dinucNSyn_stat[2][6][9];
    }
    double GetDinucNSyn31_CGGG()
    {
        return dinucNSyn_stat[2][6][10];
    }
    double GetDinucNSyn31_CGGT()
    {
        return dinucNSyn_stat[2][6][11];
    }
    double GetDinucNSyn31_CGTA()
    {
        return dinucNSyn_stat[2][6][12];
    }
    double GetDinucNSyn31_CGTC()
    {
        return dinucNSyn_stat[2][6][13];
    }
    double GetDinucNSyn31_CGTG()
    {
        return dinucNSyn_stat[2][6][14];
    }
    double GetDinucNSyn31_CGTT()
    {
        return dinucNSyn_stat[2][6][15];
    }
    double GetDinucNSyn31_CTAA()
    {
        return dinucNSyn_stat[2][7][0];
    }
    double GetDinucNSyn31_CTAC()
    {
        return dinucNSyn_stat[2][7][1];
    }
    double GetDinucNSyn31_CTAG()
    {
        return dinucNSyn_stat[2][7][2];
    }
    double GetDinucNSyn31_CTAT()
    {
        return dinucNSyn_stat[2][7][3];
    }
    double GetDinucNSyn31_CTCA()
    {
        return dinucNSyn_stat[2][7][4];
    }
    double GetDinucNSyn31_CTCC()
    {
        return dinucNSyn_stat[2][7][5];
    }
    double GetDinucNSyn31_CTCG()
    {
        return dinucNSyn_stat[2][7][6];
    }
    double GetDinucNSyn31_CTCT()
    {
        return dinucNSyn_stat[2][7][7];
    }
    double GetDinucNSyn31_CTGA()
    {
        return dinucNSyn_stat[2][7][8];
    }
    double GetDinucNSyn31_CTGC()
    {
        return dinucNSyn_stat[2][7][9];
    }
    double GetDinucNSyn31_CTGG()
    {
        return dinucNSyn_stat[2][7][10];
    }
    double GetDinucNSyn31_CTGT()
    {
        return dinucNSyn_stat[2][7][11];
    }
    double GetDinucNSyn31_CTTA()
    {
        return dinucNSyn_stat[2][7][12];
    }
    double GetDinucNSyn31_CTTC()
    {
        return dinucNSyn_stat[2][7][13];
    }
    double GetDinucNSyn31_CTTG()
    {
        return dinucNSyn_stat[2][7][14];
    }
    double GetDinucNSyn31_CTTT()
    {
        return dinucNSyn_stat[2][7][15];
    }
    double GetDinucNSyn31_GAAA()
    {
        return dinucNSyn_stat[2][8][0];
    }
    double GetDinucNSyn31_GAAC()
    {
        return dinucNSyn_stat[2][8][1];
    }
    double GetDinucNSyn31_GAAG()
    {
        return dinucNSyn_stat[2][8][2];
    }
    double GetDinucNSyn31_GAAT()
    {
        return dinucNSyn_stat[2][8][3];
    }
    double GetDinucNSyn31_GACA()
    {
        return dinucNSyn_stat[2][8][4];
    }
    double GetDinucNSyn31_GACC()
    {
        return dinucNSyn_stat[2][8][5];
    }
    double GetDinucNSyn31_GACG()
    {
        return dinucNSyn_stat[2][8][6];
    }
    double GetDinucNSyn31_GACT()
    {
        return dinucNSyn_stat[2][8][7];
    }
    double GetDinucNSyn31_GAGA()
    {
        return dinucNSyn_stat[2][8][8];
    }
    double GetDinucNSyn31_GAGC()
    {
        return dinucNSyn_stat[2][8][9];
    }
    double GetDinucNSyn31_GAGG()
    {
        return dinucNSyn_stat[2][8][10];
    }
    double GetDinucNSyn31_GAGT()
    {
        return dinucNSyn_stat[2][8][11];
    }
    double GetDinucNSyn31_GATA()
    {
        return dinucNSyn_stat[2][8][12];
    }
    double GetDinucNSyn31_GATC()
    {
        return dinucNSyn_stat[2][8][13];
    }
    double GetDinucNSyn31_GATG()
    {
        return dinucNSyn_stat[2][8][14];
    }
    double GetDinucNSyn31_GATT()
    {
        return dinucNSyn_stat[2][8][15];
    }
    double GetDinucNSyn31_GCAA()
    {
        return dinucNSyn_stat[2][9][0];
    }
    double GetDinucNSyn31_GCAC()
    {
        return dinucNSyn_stat[2][9][1];
    }
    double GetDinucNSyn31_GCAG()
    {
        return dinucNSyn_stat[2][9][2];
    }
    double GetDinucNSyn31_GCAT()
    {
        return dinucNSyn_stat[2][9][3];
    }
    double GetDinucNSyn31_GCCA()
    {
        return dinucNSyn_stat[2][9][4];
    }
    double GetDinucNSyn31_GCCC()
    {
        return dinucNSyn_stat[2][9][5];
    }
    double GetDinucNSyn31_GCCG()
    {
        return dinucNSyn_stat[2][9][6];
    }
    double GetDinucNSyn31_GCCT()
    {
        return dinucNSyn_stat[2][9][7];
    }
    double GetDinucNSyn31_GCGA()
    {
        return dinucNSyn_stat[2][9][8];
    }
    double GetDinucNSyn31_GCGC()
    {
        return dinucNSyn_stat[2][9][9];
    }
    double GetDinucNSyn31_GCGG()
    {
        return dinucNSyn_stat[2][9][10];
    }
    double GetDinucNSyn31_GCGT()
    {
        return dinucNSyn_stat[2][9][11];
    }
    double GetDinucNSyn31_GCTA()
    {
        return dinucNSyn_stat[2][9][12];
    }
    double GetDinucNSyn31_GCTC()
    {
        return dinucNSyn_stat[2][9][13];
    }
    double GetDinucNSyn31_GCTG()
    {
        return dinucNSyn_stat[2][9][14];
    }
    double GetDinucNSyn31_GCTT()
    {
        return dinucNSyn_stat[2][9][15];
    }
    double GetDinucNSyn31_GGAA()
    {
        return dinucNSyn_stat[2][10][0];
    }
    double GetDinucNSyn31_GGAC()
    {
        return dinucNSyn_stat[2][10][1];
    }
    double GetDinucNSyn31_GGAG()
    {
        return dinucNSyn_stat[2][10][2];
    }
    double GetDinucNSyn31_GGAT()
    {
        return dinucNSyn_stat[2][10][3];
    }
    double GetDinucNSyn31_GGCA()
    {
        return dinucNSyn_stat[2][10][4];
    }
    double GetDinucNSyn31_GGCC()
    {
        return dinucNSyn_stat[2][10][5];
    }
    double GetDinucNSyn31_GGCG()
    {
        return dinucNSyn_stat[2][10][6];
    }
    double GetDinucNSyn31_GGCT()
    {
        return dinucNSyn_stat[2][10][7];
    }
    double GetDinucNSyn31_GGGA()
    {
        return dinucNSyn_stat[2][10][8];
    }
    double GetDinucNSyn31_GGGC()
    {
        return dinucNSyn_stat[2][10][9];
    }
    double GetDinucNSyn31_GGGG()
    {
        return dinucNSyn_stat[2][10][10];
    }
    double GetDinucNSyn31_GGGT()
    {
        return dinucNSyn_stat[2][10][11];
    }
    double GetDinucNSyn31_GGTA()
    {
        return dinucNSyn_stat[2][10][12];
    }
    double GetDinucNSyn31_GGTC()
    {
        return dinucNSyn_stat[2][10][13];
    }
    double GetDinucNSyn31_GGTG()
    {
        return dinucNSyn_stat[2][10][14];
    }
    double GetDinucNSyn31_GGTT()
    {
        return dinucNSyn_stat[2][10][15];
    }
    double GetDinucNSyn31_GTAA()
    {
        return dinucNSyn_stat[2][11][0];
    }
    double GetDinucNSyn31_GTAC()
    {
        return dinucNSyn_stat[2][11][1];
    }
    double GetDinucNSyn31_GTAG()
    {
        return dinucNSyn_stat[2][11][2];
    }
    double GetDinucNSyn31_GTAT()
    {
        return dinucNSyn_stat[2][11][3];
    }
    double GetDinucNSyn31_GTCA()
    {
        return dinucNSyn_stat[2][11][4];
    }
    double GetDinucNSyn31_GTCC()
    {
        return dinucNSyn_stat[2][11][5];
    }
    double GetDinucNSyn31_GTCG()
    {
        return dinucNSyn_stat[2][11][6];
    }
    double GetDinucNSyn31_GTCT()
    {
        return dinucNSyn_stat[2][11][7];
    }
    double GetDinucNSyn31_GTGA()
    {
        return dinucNSyn_stat[2][11][8];
    }
    double GetDinucNSyn31_GTGC()
    {
        return dinucNSyn_stat[2][11][9];
    }
    double GetDinucNSyn31_GTGG()
    {
        return dinucNSyn_stat[2][11][10];
    }
    double GetDinucNSyn31_GTGT()
    {
        return dinucNSyn_stat[2][11][11];
    }
    double GetDinucNSyn31_GTTA()
    {
        return dinucNSyn_stat[2][11][12];
    }
    double GetDinucNSyn31_GTTC()
    {
        return dinucNSyn_stat[2][11][13];
    }
    double GetDinucNSyn31_GTTG()
    {
        return dinucNSyn_stat[2][11][14];
    }
    double GetDinucNSyn31_GTTT()
    {
        return dinucNSyn_stat[2][11][15];
    }
    double GetDinucNSyn31_TAAA()
    {
        return dinucNSyn_stat[2][12][0];
    }
    double GetDinucNSyn31_TAAC()
    {
        return dinucNSyn_stat[2][12][1];
    }
    double GetDinucNSyn31_TAAG()
    {
        return dinucNSyn_stat[2][12][2];
    }
    double GetDinucNSyn31_TAAT()
    {
        return dinucNSyn_stat[2][12][3];
    }
    double GetDinucNSyn31_TACA()
    {
        return dinucNSyn_stat[2][12][4];
    }
    double GetDinucNSyn31_TACC()
    {
        return dinucNSyn_stat[2][12][5];
    }
    double GetDinucNSyn31_TACG()
    {
        return dinucNSyn_stat[2][12][6];
    }
    double GetDinucNSyn31_TACT()
    {
        return dinucNSyn_stat[2][12][7];
    }
    double GetDinucNSyn31_TAGA()
    {
        return dinucNSyn_stat[2][12][8];
    }
    double GetDinucNSyn31_TAGC()
    {
        return dinucNSyn_stat[2][12][9];
    }
    double GetDinucNSyn31_TAGG()
    {
        return dinucNSyn_stat[2][12][10];
    }
    double GetDinucNSyn31_TAGT()
    {
        return dinucNSyn_stat[2][12][11];
    }
    double GetDinucNSyn31_TATA()
    {
        return dinucNSyn_stat[2][12][12];
    }
    double GetDinucNSyn31_TATC()
    {
        return dinucNSyn_stat[2][12][13];
    }
    double GetDinucNSyn31_TATG()
    {
        return dinucNSyn_stat[2][12][14];
    }
    double GetDinucNSyn31_TATT()
    {
        return dinucNSyn_stat[2][12][15];
    }
    double GetDinucNSyn31_TCAA()
    {
        return dinucNSyn_stat[2][13][0];
    }
    double GetDinucNSyn31_TCAC()
    {
        return dinucNSyn_stat[2][13][1];
    }
    double GetDinucNSyn31_TCAG()
    {
        return dinucNSyn_stat[2][13][2];
    }
    double GetDinucNSyn31_TCAT()
    {
        return dinucNSyn_stat[2][13][3];
    }
    double GetDinucNSyn31_TCCA()
    {
        return dinucNSyn_stat[2][13][4];
    }
    double GetDinucNSyn31_TCCC()
    {
        return dinucNSyn_stat[2][13][5];
    }
    double GetDinucNSyn31_TCCG()
    {
        return dinucNSyn_stat[2][13][6];
    }
    double GetDinucNSyn31_TCCT()
    {
        return dinucNSyn_stat[2][13][7];
    }
    double GetDinucNSyn31_TCGA()
    {
        return dinucNSyn_stat[2][13][8];
    }
    double GetDinucNSyn31_TCGC()
    {
        return dinucNSyn_stat[2][13][9];
    }
    double GetDinucNSyn31_TCGG()
    {
        return dinucNSyn_stat[2][13][10];
    }
    double GetDinucNSyn31_TCGT()
    {
        return dinucNSyn_stat[2][13][11];
    }
    double GetDinucNSyn31_TCTA()
    {
        return dinucNSyn_stat[2][13][12];
    }
    double GetDinucNSyn31_TCTC()
    {
        return dinucNSyn_stat[2][13][13];
    }
    double GetDinucNSyn31_TCTG()
    {
        return dinucNSyn_stat[2][13][14];
    }
    double GetDinucNSyn31_TCTT()
    {
        return dinucNSyn_stat[2][13][15];
    }
    double GetDinucNSyn31_TGAA()
    {
        return dinucNSyn_stat[2][14][0];
    }
    double GetDinucNSyn31_TGAC()
    {
        return dinucNSyn_stat[2][14][1];
    }
    double GetDinucNSyn31_TGAG()
    {
        return dinucNSyn_stat[2][14][2];
    }
    double GetDinucNSyn31_TGAT()
    {
        return dinucNSyn_stat[2][14][3];
    }
    double GetDinucNSyn31_TGCA()
    {
        return dinucNSyn_stat[2][14][4];
    }
    double GetDinucNSyn31_TGCC()
    {
        return dinucNSyn_stat[2][14][5];
    }
    double GetDinucNSyn31_TGCG()
    {
        return dinucNSyn_stat[2][14][6];
    }
    double GetDinucNSyn31_TGCT()
    {
        return dinucNSyn_stat[2][14][7];
    }
    double GetDinucNSyn31_TGGA()
    {
        return dinucNSyn_stat[2][14][8];
    }
    double GetDinucNSyn31_TGGC()
    {
        return dinucNSyn_stat[2][14][9];
    }
    double GetDinucNSyn31_TGGG()
    {
        return dinucNSyn_stat[2][14][10];
    }
    double GetDinucNSyn31_TGGT()
    {
        return dinucNSyn_stat[2][14][11];
    }
    double GetDinucNSyn31_TGTA()
    {
        return dinucNSyn_stat[2][14][12];
    }
    double GetDinucNSyn31_TGTC()
    {
        return dinucNSyn_stat[2][14][13];
    }
    double GetDinucNSyn31_TGTG()
    {
        return dinucNSyn_stat[2][14][14];
    }
    double GetDinucNSyn31_TGTT()
    {
        return dinucNSyn_stat[2][14][15];
    }
    double GetDinucNSyn31_TTAA()
    {
        return dinucNSyn_stat[2][15][0];
    }
    double GetDinucNSyn31_TTAC()
    {
        return dinucNSyn_stat[2][15][1];
    }
    double GetDinucNSyn31_TTAG()
    {
        return dinucNSyn_stat[2][15][2];
    }
    double GetDinucNSyn31_TTAT()
    {
        return dinucNSyn_stat[2][15][3];
    }
    double GetDinucNSyn31_TTCA()
    {
        return dinucNSyn_stat[2][15][4];
    }
    double GetDinucNSyn31_TTCC()
    {
        return dinucNSyn_stat[2][15][5];
    }
    double GetDinucNSyn31_TTCG()
    {
        return dinucNSyn_stat[2][15][6];
    }
    double GetDinucNSyn31_TTCT()
    {
        return dinucNSyn_stat[2][15][7];
    }
    double GetDinucNSyn31_TTGA()
    {
        return dinucNSyn_stat[2][15][8];
    }
    double GetDinucNSyn31_TTGC()
    {
        return dinucNSyn_stat[2][15][9];
    }
    double GetDinucNSyn31_TTGG()
    {
        return dinucNSyn_stat[2][15][10];
    }
    double GetDinucNSyn31_TTGT()
    {
        return dinucNSyn_stat[2][15][11];
    }
    double GetDinucNSyn31_TTTA()
    {
        return dinucNSyn_stat[2][15][12];
    }
    double GetDinucNSyn31_TTTC()
    {
        return dinucNSyn_stat[2][15][13];
    }
    double GetDinucNSyn31_TTTG()
    {
        return dinucNSyn_stat[2][15][14];
    }
    double GetDinucNSyn31_TTTT()
    {
        return dinucNSyn_stat[2][15][15];
    }

    double GetNsub()
    {
        return (Nsub/lparam->Nsite_codon);
    }
    double GetNSynsub()
    {
        return (Nsynsub/lparam->Nsite_codon);
    }

    double GetMutRateStart()
    {
        return MutRate[0][0];
    }

    double GetSubRateStart()
    {
        return SubRate[0][0];
    }

    double GetMutRateNonSynStart()
    {
        return MutRate[0][1];
    }

    double GetSubRateNonSynStart()
    {
        return SubRate[0][1];
    }

    double GetMutRateSynStart()
    {
        return MutRate[0][2];
    }

    double GetSubRateSynStart()
    {
        return SubRate[0][2];
    }

    double GetMutRateEnd()
    {
        return MutRate[1][0];
    }

    double GetSubRateEnd()
    {
        return SubRate[1][0];
    }

    double GetMutRateNonSynEnd()
    {
        return MutRate[1][1];
    }

    double GetSubRateNonSynEnd()
    {
        return SubRate[1][1];
    }
    double GetMutRateSynEnd()
    {
        return MutRate[1][2];
    }

    double GetSubRateSynEnd()
    {
        return SubRate[1][2];
    }


};

#endif // EVOLHISTSTATISTICS_H
