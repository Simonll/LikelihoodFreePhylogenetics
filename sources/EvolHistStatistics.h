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
        double* ssNnsynsub;
        double* ssNsynsub;

        double** branch_stat;
        double*** gtnr_stat;
        double*** dinuc_stat;
        //Synonymous
        double*** gtnrSyn_stat;
        double*** dinucSyn_stat;
        //Nonsynonymous
        double*** gtnrNSyn_stat;
        double*** dinucNSyn_stat;

        typedef double (EvolHistStatistics::*funcpt)();
        std::map<string,funcpt> GetEvoStatMap;



        EvolHistStatistics(LocalParameters* inparam);

        virtual ~EvolHistStatistics();

        void resetMappingStat();
        void GetMapStats();
        void GetMapAncStats();
        int GetDinucContext(int nuc_a, int nuc_b) {
            int k = 0;
            for (int i = 0; i < 4; i++){
                for (int j =0 ; j <4; j++){
                    if(nuc_a == i && nuc_b == j){
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

   double GetGTNR_AA(){
       return gtnr_stat[3][0][0];
   }
   double GetGTNR_AC(){
       return gtnr_stat[3][0][1];
   }
   double GetGTNR_AG(){
       return gtnr_stat[3][0][2];
   }
   double GetGTNR_AT(){
       return gtnr_stat[3][0][3];
   }
   double GetGTNR_CA(){
       return gtnr_stat[3][1][0];
   }
   double GetGTNR_CC(){
       return gtnr_stat[3][1][1];
   }
   double GetGTNR_CG(){
       return gtnr_stat[3][1][2];
   }
   double GetGTNR_CT(){
       return gtnr_stat[3][1][3];
   }
   double GetGTNR_GA(){
       return gtnr_stat[3][2][0];
   }
   double GetGTNR_GC(){
       return gtnr_stat[3][2][1];
   }
   double GetGTNR_GG(){
       return gtnr_stat[3][2][2];
   }
   double GetGTNR_GT(){
       return gtnr_stat[3][2][3];
   }
   double GetGTNR_TA(){
       return gtnr_stat[3][3][0];
   }
   double GetGTNR_TC(){
       return gtnr_stat[3][3][1];
   }
   double GetGTNR_TG(){
       return gtnr_stat[3][3][2];
   }
   double GetGTNR_TT(){
       return gtnr_stat[3][3][3];
   }
   double GetGTNR1_AA(){
   return gtnr_stat[0][0][0];
   }
   double GetGTNR1_AC(){
       return gtnr_stat[0][0][1];
   }
   double GetGTNR1_AG(){
       return gtnr_stat[0][0][2];
   }
   double GetGTNR1_AT(){
       return gtnr_stat[0][0][3];
   }
   double GetGTNR1_CA(){
       return gtnr_stat[0][1][0];
   }
   double GetGTNR1_CC(){
       return gtnr_stat[0][1][1];
   }
   double GetGTNR1_CG(){
       return gtnr_stat[0][1][2];
   }
   double GetGTNR1_CT(){
       return gtnr_stat[0][1][3];
   }
   double GetGTNR1_GA(){
       return gtnr_stat[0][2][0];
   }
   double GetGTNR1_GC(){
       return gtnr_stat[0][2][1];
   }
   double GetGTNR1_GG(){
       return gtnr_stat[0][2][2];
   }
   double GetGTNR1_GT(){
       return gtnr_stat[0][2][3];
   }
   double GetGTNR1_TA(){
       return gtnr_stat[0][3][0];
   }
   double GetGTNR1_TC(){
       return gtnr_stat[0][3][1];
   }
   double GetGTNR1_TG(){
       return gtnr_stat[0][3][2];
   }
   double GetGTNR1_TT(){
       return gtnr_stat[0][3][3];
   }



   double GetGTNR2_AA(){
       return gtnr_stat[1][0][0];
   }
   double GetGTNR2_AC(){
       return gtnr_stat[1][0][1];
   }
   double GetGTNR2_AG(){
       return gtnr_stat[1][0][2];
   }
   double GetGTNR2_AT(){
       return gtnr_stat[1][0][3];
   }
   double GetGTNR2_CA(){
       return gtnr_stat[1][1][0];
   }
   double GetGTNR2_CC(){
       return gtnr_stat[1][1][1];
   }
   double GetGTNR2_CG(){
       return gtnr_stat[1][1][2];
   }
   double GetGTNR2_CT(){
       return gtnr_stat[1][1][3];
   }
   double GetGTNR2_GA(){
       return gtnr_stat[1][2][0];
   }
   double GetGTNR2_GC(){
       return gtnr_stat[1][2][1];
   }
   double GetGTNR2_GG(){
       return gtnr_stat[1][2][2];
   }
   double GetGTNR2_GT(){
       return gtnr_stat[1][2][3];
   }
   double GetGTNR2_TA(){
       return gtnr_stat[1][3][0];
   }
   double GetGTNR2_TC(){
       return gtnr_stat[1][3][1];
   }
   double GetGTNR2_TG(){
       return gtnr_stat[1][3][2];
   }
   double GetGTNR2_TT(){
       return gtnr_stat[1][3][3];
   }
   double GetGTNR3_AA(){
   return gtnr_stat[2][0][0];
}
   double GetGTNR3_AC(){
       return gtnr_stat[2][0][1];
   }
   double GetGTNR3_AG(){
       return gtnr_stat[2][0][2];
   }
   double GetGTNR3_AT(){
       return gtnr_stat[2][0][3];
   }
   double GetGTNR3_CA(){
       return gtnr_stat[2][1][0];
   }
   double GetGTNR3_CC(){
       return gtnr_stat[2][1][1];
   }
   double GetGTNR3_CG(){
       return gtnr_stat[2][1][2];
   }
   double GetGTNR3_CT(){
       return gtnr_stat[2][1][3];
   }
   double GetGTNR3_GA(){
       return gtnr_stat[2][2][0];
   }
   double GetGTNR3_GC(){
       return gtnr_stat[2][2][1];
   }
   double GetGTNR3_GG(){
       return gtnr_stat[2][2][2];
   }
   double GetGTNR3_GT(){
       return gtnr_stat[2][2][3];
   }
   double GetGTNR3_TA(){
       return gtnr_stat[2][3][0];
   }
   double GetGTNR3_TC(){
       return gtnr_stat[2][3][1];
   }
   double GetGTNR3_TG(){
       return gtnr_stat[2][3][2];
   }
   double GetGTNR3_TT(){
       return gtnr_stat[2][3][3];
   }

   double GetGTNRSyn_AA(){
       return gtnrSyn_stat[3][0][0];
   }
   double GetGTNRSyn_AC(){
       return gtnrSyn_stat[3][0][1];
   }
   double GetGTNRSyn_AG(){
       return gtnrSyn_stat[3][0][2];
   }
   double GetGTNRSyn_AT(){
       return gtnrSyn_stat[3][0][3];
   }
   double GetGTNRSyn_CA(){
       return gtnrSyn_stat[3][1][0];
   }
   double GetGTNRSyn_CC(){
       return gtnrSyn_stat[3][1][1];
   }
   double GetGTNRSyn_CG(){
       return gtnrSyn_stat[3][1][2];
   }
   double GetGTNRSyn_CT(){
       return gtnrSyn_stat[3][1][3];
   }
   double GetGTNRSyn_GA(){
       return gtnrSyn_stat[3][2][0];
   }
   double GetGTNRSyn_GC(){
       return gtnrSyn_stat[3][2][1];
   }
   double GetGTNRSyn_GG(){
       return gtnrSyn_stat[3][2][2];
   }
   double GetGTNRSyn_GT(){
       return gtnrSyn_stat[3][2][3];
   }
   double GetGTNRSyn_TA(){
       return gtnrSyn_stat[3][3][0];
   }
   double GetGTNRSyn_TC(){
       return gtnrSyn_stat[3][3][1];
   }
   double GetGTNRSyn_TG(){
       return gtnrSyn_stat[3][3][2];
   }
   double GetGTNRSyn_TT(){
       return gtnrSyn_stat[3][3][3];
   }
   double GetGTNRSyn1_AA(){
   return gtnrSyn_stat[0][0][0];
   }
   double GetGTNRSyn1_AC(){
       return gtnrSyn_stat[0][0][1];
   }
   double GetGTNRSyn1_AG(){
       return gtnrSyn_stat[0][0][2];
   }
   double GetGTNRSyn1_AT(){
       return gtnrSyn_stat[0][0][3];
   }
   double GetGTNRSyn1_CA(){
       return gtnrSyn_stat[0][1][0];
   }
   double GetGTNRSyn1_CC(){
       return gtnrSyn_stat[0][1][1];
   }
   double GetGTNRSyn1_CG(){
       return gtnrSyn_stat[0][1][2];
   }
   double GetGTNRSyn1_CT(){
       return gtnrSyn_stat[0][1][3];
   }
   double GetGTNRSyn1_GA(){
       return gtnrSyn_stat[0][2][0];
   }
   double GetGTNRSyn1_GC(){
       return gtnrSyn_stat[0][2][1];
   }
   double GetGTNRSyn1_GG(){
       return gtnrSyn_stat[0][2][2];
   }
   double GetGTNRSyn1_GT(){
       return gtnrSyn_stat[0][2][3];
   }
   double GetGTNRSyn1_TA(){
       return gtnrSyn_stat[0][3][0];
   }
   double GetGTNRSyn1_TC(){
       return gtnrSyn_stat[0][3][1];
   }
   double GetGTNRSyn1_TG(){
       return gtnrSyn_stat[0][3][2];
   }
   double GetGTNRSyn1_TT(){
       return gtnrSyn_stat[0][3][3];
   }



   double GetGTNRSyn2_AA(){
       return gtnrSyn_stat[1][0][0];
   }
   double GetGTNRSyn2_AC(){
       return gtnrSyn_stat[1][0][1];
   }
   double GetGTNRSyn2_AG(){
       return gtnrSyn_stat[1][0][2];
   }
   double GetGTNRSyn2_AT(){
       return gtnrSyn_stat[1][0][3];
   }
   double GetGTNRSyn2_CA(){
       return gtnrSyn_stat[1][1][0];
   }
   double GetGTNRSyn2_CC(){
       return gtnrSyn_stat[1][1][1];
   }
   double GetGTNRSyn2_CG(){
       return gtnrSyn_stat[1][1][2];
   }
   double GetGTNRSyn2_CT(){
       return gtnrSyn_stat[1][1][3];
   }
   double GetGTNRSyn2_GA(){
       return gtnrSyn_stat[1][2][0];
   }
   double GetGTNRSyn2_GC(){
       return gtnrSyn_stat[1][2][1];
   }
   double GetGTNRSyn2_GG(){
       return gtnrSyn_stat[1][2][2];
   }
   double GetGTNRSyn2_GT(){
       return gtnrSyn_stat[1][2][3];
   }
   double GetGTNRSyn2_TA(){
       return gtnrSyn_stat[1][3][0];
   }
   double GetGTNRSyn2_TC(){
       return gtnrSyn_stat[1][3][1];
   }
   double GetGTNRSyn2_TG(){
       return gtnrSyn_stat[1][3][2];
   }
   double GetGTNRSyn2_TT(){
       return gtnrSyn_stat[1][3][3];
   }
   double GetGTNRSyn3_AA(){
   return gtnrSyn_stat[2][0][0];
}
   double GetGTNRSyn3_AC(){
       return gtnrSyn_stat[2][0][1];
   }
   double GetGTNRSyn3_AG(){
       return gtnrSyn_stat[2][0][2];
   }
   double GetGTNRSyn3_AT(){
       return gtnrSyn_stat[2][0][3];
   }
   double GetGTNRSyn3_CA(){
       return gtnrSyn_stat[2][1][0];
   }
   double GetGTNRSyn3_CC(){
       return gtnrSyn_stat[2][1][1];
   }
   double GetGTNRSyn3_CG(){
       return gtnrSyn_stat[2][1][2];
   }
   double GetGTNRSyn3_CT(){
       return gtnrSyn_stat[2][1][3];
   }
   double GetGTNRSyn3_GA(){
       return gtnrSyn_stat[2][2][0];
   }
   double GetGTNRSyn3_GC(){
       return gtnrSyn_stat[2][2][1];
   }
   double GetGTNRSyn3_GG(){
       return gtnrSyn_stat[2][2][2];
   }
   double GetGTNRSyn3_GT(){
       return gtnrSyn_stat[2][2][3];
   }
   double GetGTNRSyn3_TA(){
       return gtnrSyn_stat[2][3][0];
   }
   double GetGTNRSyn3_TC(){
       return gtnrSyn_stat[2][3][1];
   }
   double GetGTNRSyn3_TG(){
       return gtnrSyn_stat[2][3][2];
   }
   double GetGTNRSyn3_TT(){
       return gtnrSyn_stat[2][3][3];
   }
   double GetGTNRNSyn_AA(){
       return gtnrNSyn_stat[3][0][0];
   }
   double GetGTNRNSyn_AC(){
       return gtnrNSyn_stat[3][0][1];
   }
   double GetGTNRNSyn_AG(){
       return gtnrNSyn_stat[3][0][2];
   }
   double GetGTNRNSyn_AT(){
       return gtnrNSyn_stat[3][0][3];
   }
   double GetGTNRNSyn_CA(){
       return gtnrNSyn_stat[3][1][0];
   }
   double GetGTNRNSyn_CC(){
       return gtnrNSyn_stat[3][1][1];
   }
   double GetGTNRNSyn_CG(){
       return gtnrNSyn_stat[3][1][2];
   }
   double GetGTNRNSyn_CT(){
       return gtnrNSyn_stat[3][1][3];
   }
   double GetGTNRNSyn_GA(){
       return gtnrNSyn_stat[3][2][0];
   }
   double GetGTNRNSyn_GC(){
       return gtnrNSyn_stat[3][2][1];
   }
   double GetGTNRNSyn_GG(){
       return gtnrNSyn_stat[3][2][2];
   }
   double GetGTNRNSyn_GT(){
       return gtnrNSyn_stat[3][2][3];
   }
   double GetGTNRNSyn_TA(){
       return gtnrNSyn_stat[3][3][0];
   }
   double GetGTNRNSyn_TC(){
       return gtnrNSyn_stat[3][3][1];
   }
   double GetGTNRNSyn_TG(){
       return gtnrNSyn_stat[3][3][2];
   }
   double GetGTNRNSyn_TT(){
       return gtnrNSyn_stat[3][3][3];
   }
   double GetGTNRNSyn1_AA(){
   return gtnrNSyn_stat[0][0][0];
   }
   double GetGTNRNSyn1_AC(){
       return gtnrNSyn_stat[0][0][1];
   }
   double GetGTNRNSyn1_AG(){
       return gtnrNSyn_stat[0][0][2];
   }
   double GetGTNRNSyn1_AT(){
       return gtnrNSyn_stat[0][0][3];
   }
   double GetGTNRNSyn1_CA(){
       return gtnrNSyn_stat[0][1][0];
   }
   double GetGTNRNSyn1_CC(){
       return gtnrNSyn_stat[0][1][1];
   }
   double GetGTNRNSyn1_CG(){
       return gtnrNSyn_stat[0][1][2];
   }
   double GetGTNRNSyn1_CT(){
       return gtnrNSyn_stat[0][1][3];
   }
   double GetGTNRNSyn1_GA(){
       return gtnrNSyn_stat[0][2][0];
   }
   double GetGTNRNSyn1_GC(){
       return gtnrNSyn_stat[0][2][1];
   }
   double GetGTNRNSyn1_GG(){
       return gtnrNSyn_stat[0][2][2];
   }
   double GetGTNRNSyn1_GT(){
       return gtnrNSyn_stat[0][2][3];
   }
   double GetGTNRNSyn1_TA(){
       return gtnrNSyn_stat[0][3][0];
   }
   double GetGTNRNSyn1_TC(){
       return gtnrNSyn_stat[0][3][1];
   }
   double GetGTNRNSyn1_TG(){
       return gtnrNSyn_stat[0][3][2];
   }
   double GetGTNRNSyn1_TT(){
       return gtnrNSyn_stat[0][3][3];
   }



   double GetGTNRNSyn2_AA(){
       return gtnrNSyn_stat[1][0][0];
   }
   double GetGTNRNSyn2_AC(){
       return gtnrNSyn_stat[1][0][1];
   }
   double GetGTNRNSyn2_AG(){
       return gtnrNSyn_stat[1][0][2];
   }
   double GetGTNRNSyn2_AT(){
       return gtnrNSyn_stat[1][0][3];
   }
   double GetGTNRNSyn2_CA(){
       return gtnrNSyn_stat[1][1][0];
   }
   double GetGTNRNSyn2_CC(){
       return gtnrNSyn_stat[1][1][1];
   }
   double GetGTNRNSyn2_CG(){
       return gtnrNSyn_stat[1][1][2];
   }
   double GetGTNRNSyn2_CT(){
       return gtnrNSyn_stat[1][1][3];
   }
   double GetGTNRNSyn2_GA(){
       return gtnrNSyn_stat[1][2][0];
   }
   double GetGTNRNSyn2_GC(){
       return gtnrNSyn_stat[1][2][1];
   }
   double GetGTNRNSyn2_GG(){
       return gtnrNSyn_stat[1][2][2];
   }
   double GetGTNRNSyn2_GT(){
       return gtnrNSyn_stat[1][2][3];
   }
   double GetGTNRNSyn2_TA(){
       return gtnrNSyn_stat[1][3][0];
   }
   double GetGTNRNSyn2_TC(){
       return gtnrNSyn_stat[1][3][1];
   }
   double GetGTNRNSyn2_TG(){
       return gtnrNSyn_stat[1][3][2];
   }
   double GetGTNRNSyn2_TT(){
       return gtnrNSyn_stat[1][3][3];
   }
   double GetGTNRNSyn3_AA(){
   return gtnrNSyn_stat[2][0][0];
}
   double GetGTNRNSyn3_AC(){
       return gtnrNSyn_stat[2][0][1];
   }
   double GetGTNRNSyn3_AG(){
       return gtnrNSyn_stat[2][0][2];
   }
   double GetGTNRNSyn3_AT(){
       return gtnrNSyn_stat[2][0][3];
   }
   double GetGTNRNSyn3_CA(){
       return gtnrNSyn_stat[2][1][0];
   }
   double GetGTNRNSyn3_CC(){
       return gtnrNSyn_stat[2][1][1];
   }
   double GetGTNRNSyn3_CG(){
       return gtnrNSyn_stat[2][1][2];
   }
   double GetGTNRNSyn3_CT(){
       return gtnrNSyn_stat[2][1][3];
   }
   double GetGTNRNSyn3_GA(){
       return gtnrNSyn_stat[2][2][0];
   }
   double GetGTNRNSyn3_GC(){
       return gtnrNSyn_stat[2][2][1];
   }
   double GetGTNRNSyn3_GG(){
       return gtnrNSyn_stat[2][2][2];
   }
   double GetGTNRNSyn3_GT(){
       return gtnrNSyn_stat[2][2][3];
   }
   double GetGTNRNSyn3_TA(){
       return gtnrNSyn_stat[2][3][0];
   }
   double GetGTNRNSyn3_TC(){
       return gtnrNSyn_stat[2][3][1];
   }
   double GetGTNRNSyn3_TG(){
       return gtnrNSyn_stat[2][3][2];
   }
   double GetGTNRNSyn3_TT(){
       return gtnrNSyn_stat[2][3][3];
   }

 double GetDinuc_AAAA(){
     return dinuc_stat[3][0][0];
 }
 double GetDinuc_AAAC(){
     return dinuc_stat[3][0][1];
 }
 double GetDinuc_AAAG(){
     return dinuc_stat[3][0][2];
 }
 double GetDinuc_AAAT(){
     return dinuc_stat[3][0][3];
 }
 double GetDinuc_AACA(){
     return dinuc_stat[3][0][4];
 }
 double GetDinuc_AACC(){
     return dinuc_stat[3][0][5];
 }
 double GetDinuc_AACG(){
     return dinuc_stat[3][0][6];
 }
 double GetDinuc_AACT(){
     return dinuc_stat[3][0][7];
 }
 double GetDinuc_AAGA(){
     return dinuc_stat[3][0][8];
 }
 double GetDinuc_AAGC(){
     return dinuc_stat[3][0][9];
 }
 double GetDinuc_AAGG(){
     return dinuc_stat[3][0][10];
 }
 double GetDinuc_AAGT(){
     return dinuc_stat[3][0][11];
 }
 double GetDinuc_AATA(){
     return dinuc_stat[3][0][12];
 }
 double GetDinuc_AATC(){
     return dinuc_stat[3][0][13];
 }
 double GetDinuc_AATG(){
     return dinuc_stat[3][0][14];
 }
 double GetDinuc_AATT(){
     return dinuc_stat[3][0][15];
 }
 double GetDinuc_ACAA(){
     return dinuc_stat[3][1][0];
 }
 double GetDinuc_ACAC(){
     return dinuc_stat[3][1][1];
 }
 double GetDinuc_ACAG(){
     return dinuc_stat[3][1][2];
 }
 double GetDinuc_ACAT(){
     return dinuc_stat[3][1][3];
 }
 double GetDinuc_ACCA(){
     return dinuc_stat[3][1][4];
 }
 double GetDinuc_ACCC(){
     return dinuc_stat[3][1][5];
 }
 double GetDinuc_ACCG(){
     return dinuc_stat[3][1][6];
 }
 double GetDinuc_ACCT(){
     return dinuc_stat[3][1][7];
 }
 double GetDinuc_ACGA(){
     return dinuc_stat[3][1][8];
 }
 double GetDinuc_ACGC(){
     return dinuc_stat[3][1][9];
 }
 double GetDinuc_ACGG(){
     return dinuc_stat[3][1][10];
 }
 double GetDinuc_ACGT(){
     return dinuc_stat[3][1][11];
 }
 double GetDinuc_ACTA(){
     return dinuc_stat[3][1][12];
 }
 double GetDinuc_ACTC(){
     return dinuc_stat[3][1][13];
 }
 double GetDinuc_ACTG(){
     return dinuc_stat[3][1][14];
 }
 double GetDinuc_ACTT(){
     return dinuc_stat[3][1][15];
 }
 double GetDinuc_AGAA(){
     return dinuc_stat[3][2][0];
 }
 double GetDinuc_AGAC(){
     return dinuc_stat[3][2][1];
 }
 double GetDinuc_AGAG(){
     return dinuc_stat[3][2][2];
 }
 double GetDinuc_AGAT(){
     return dinuc_stat[3][2][3];
 }
 double GetDinuc_AGCA(){
     return dinuc_stat[3][2][4];
 }
 double GetDinuc_AGCC(){
     return dinuc_stat[3][2][5];
 }
 double GetDinuc_AGCG(){
     return dinuc_stat[3][2][6];
 }
 double GetDinuc_AGCT(){
     return dinuc_stat[3][2][7];
 }
 double GetDinuc_AGGA(){
     return dinuc_stat[3][2][8];
 }
 double GetDinuc_AGGC(){
     return dinuc_stat[3][2][9];
 }
 double GetDinuc_AGGG(){
     return dinuc_stat[3][2][10];
 }
 double GetDinuc_AGGT(){
     return dinuc_stat[3][2][11];
 }
 double GetDinuc_AGTA(){
     return dinuc_stat[3][2][12];
 }
 double GetDinuc_AGTC(){
     return dinuc_stat[3][2][13];
 }
 double GetDinuc_AGTG(){
     return dinuc_stat[3][2][14];
 }
 double GetDinuc_AGTT(){
     return dinuc_stat[3][2][15];
 }
 double GetDinuc_ATAA(){
     return dinuc_stat[3][3][0];
 }
 double GetDinuc_ATAC(){
     return dinuc_stat[3][3][1];
 }
 double GetDinuc_ATAG(){
     return dinuc_stat[3][3][2];
 }
 double GetDinuc_ATAT(){
     return dinuc_stat[3][3][3];
 }
 double GetDinuc_ATCA(){
     return dinuc_stat[3][3][4];
 }
 double GetDinuc_ATCC(){
     return dinuc_stat[3][3][5];
 }
 double GetDinuc_ATCG(){
     return dinuc_stat[3][3][6];
 }
 double GetDinuc_ATCT(){
     return dinuc_stat[3][3][7];
 }
 double GetDinuc_ATGA(){
     return dinuc_stat[3][3][8];
 }
 double GetDinuc_ATGC(){
     return dinuc_stat[3][3][9];
 }
 double GetDinuc_ATGG(){
     return dinuc_stat[3][3][10];
 }
 double GetDinuc_ATGT(){
     return dinuc_stat[3][3][11];
 }
 double GetDinuc_ATTA(){
     return dinuc_stat[3][3][12];
 }
 double GetDinuc_ATTC(){
     return dinuc_stat[3][3][13];
 }
 double GetDinuc_ATTG(){
     return dinuc_stat[3][3][14];
 }
 double GetDinuc_ATTT(){
     return dinuc_stat[3][3][15];
 }
 double GetDinuc_CAAA(){
     return dinuc_stat[3][4][0];
 }
 double GetDinuc_CAAC(){
     return dinuc_stat[3][4][1];
 }
 double GetDinuc_CAAG(){
     return dinuc_stat[3][4][2];
 }
 double GetDinuc_CAAT(){
     return dinuc_stat[3][4][3];
 }
 double GetDinuc_CACA(){
     return dinuc_stat[3][4][4];
 }
 double GetDinuc_CACC(){
     return dinuc_stat[3][4][5];
 }
 double GetDinuc_CACG(){
     return dinuc_stat[3][4][6];
 }
 double GetDinuc_CACT(){
     return dinuc_stat[3][4][7];
 }
 double GetDinuc_CAGA(){
     return dinuc_stat[3][4][8];
 }
 double GetDinuc_CAGC(){
     return dinuc_stat[3][4][9];
 }
 double GetDinuc_CAGG(){
     return dinuc_stat[3][4][10];
 }
 double GetDinuc_CAGT(){
     return dinuc_stat[3][4][11];
 }
 double GetDinuc_CATA(){
     return dinuc_stat[3][4][12];
 }
 double GetDinuc_CATC(){
     return dinuc_stat[3][4][13];
 }
 double GetDinuc_CATG(){
     return dinuc_stat[3][4][14];
 }
 double GetDinuc_CATT(){
     return dinuc_stat[3][4][15];
 }
 double GetDinuc_CCAA(){
     return dinuc_stat[3][5][0];
 }
 double GetDinuc_CCAC(){
     return dinuc_stat[3][5][1];
 }
 double GetDinuc_CCAG(){
     return dinuc_stat[3][5][2];
 }
 double GetDinuc_CCAT(){
     return dinuc_stat[3][5][3];
 }
 double GetDinuc_CCCA(){
     return dinuc_stat[3][5][4];
 }
 double GetDinuc_CCCC(){
     return dinuc_stat[3][5][5];
 }
 double GetDinuc_CCCG(){
     return dinuc_stat[3][5][6];
 }
 double GetDinuc_CCCT(){
     return dinuc_stat[3][5][7];
 }
 double GetDinuc_CCGA(){
     return dinuc_stat[3][5][8];
 }
 double GetDinuc_CCGC(){
     return dinuc_stat[3][5][9];
 }
 double GetDinuc_CCGG(){
     return dinuc_stat[3][5][10];
 }
 double GetDinuc_CCGT(){
     return dinuc_stat[3][5][11];
 }
 double GetDinuc_CCTA(){
     return dinuc_stat[3][5][12];
 }
 double GetDinuc_CCTC(){
     return dinuc_stat[3][5][13];
 }
 double GetDinuc_CCTG(){
     return dinuc_stat[3][5][14];
 }
 double GetDinuc_CCTT(){
     return dinuc_stat[3][5][15];
 }
 double GetDinuc_CGAA(){
     return dinuc_stat[3][6][0];
 }
 double GetDinuc_CGAC(){
     return dinuc_stat[3][6][1];
 }
 double GetDinuc_CGAG(){
     return dinuc_stat[3][6][2];
 }
 double GetDinuc_CGAT(){
     return dinuc_stat[3][6][3];
 }
 double GetDinuc_CGCA(){
     return dinuc_stat[3][6][4];
 }
 double GetDinuc_CGCC(){
     return dinuc_stat[3][6][5];
 }
 double GetDinuc_CGCG(){
     return dinuc_stat[3][6][6];
 }
 double GetDinuc_CGCT(){
     return dinuc_stat[3][6][7];
 }
 double GetDinuc_CGGA(){
     return dinuc_stat[3][6][8];
 }
 double GetDinuc_CGGC(){
     return dinuc_stat[3][6][9];
 }
 double GetDinuc_CGGG(){
     return dinuc_stat[3][6][10];
 }
 double GetDinuc_CGGT(){
     return dinuc_stat[3][6][11];
 }
 double GetDinuc_CGTA(){
     return dinuc_stat[3][6][12];
 }
 double GetDinuc_CGTC(){
     return dinuc_stat[3][6][13];
 }
 double GetDinuc_CGTG(){
     return dinuc_stat[3][6][14];
 }
 double GetDinuc_CGTT(){
     return dinuc_stat[3][6][15];
 }
 double GetDinuc_CTAA(){
     return dinuc_stat[3][7][0];
 }
 double GetDinuc_CTAC(){
     return dinuc_stat[3][7][1];
 }
 double GetDinuc_CTAG(){
     return dinuc_stat[3][7][2];
 }
 double GetDinuc_CTAT(){
     return dinuc_stat[3][7][3];
 }
 double GetDinuc_CTCA(){
     return dinuc_stat[3][7][4];
 }
 double GetDinuc_CTCC(){
     return dinuc_stat[3][7][5];
 }
 double GetDinuc_CTCG(){
     return dinuc_stat[3][7][6];
 }
 double GetDinuc_CTCT(){
     return dinuc_stat[3][7][7];
 }
 double GetDinuc_CTGA(){
     return dinuc_stat[3][7][8];
 }
 double GetDinuc_CTGC(){
     return dinuc_stat[3][7][9];
 }
 double GetDinuc_CTGG(){
     return dinuc_stat[3][7][10];
 }
 double GetDinuc_CTGT(){
     return dinuc_stat[3][7][11];
 }
 double GetDinuc_CTTA(){
     return dinuc_stat[3][7][12];
 }
 double GetDinuc_CTTC(){
     return dinuc_stat[3][7][13];
 }
 double GetDinuc_CTTG(){
     return dinuc_stat[3][7][14];
 }
 double GetDinuc_CTTT(){
     return dinuc_stat[3][7][15];
 }
 double GetDinuc_GAAA(){
     return dinuc_stat[3][8][0];
 }
 double GetDinuc_GAAC(){
     return dinuc_stat[3][8][1];
 }
 double GetDinuc_GAAG(){
     return dinuc_stat[3][8][2];
 }
 double GetDinuc_GAAT(){
     return dinuc_stat[3][8][3];
 }
 double GetDinuc_GACA(){
     return dinuc_stat[3][8][4];
 }
 double GetDinuc_GACC(){
     return dinuc_stat[3][8][5];
 }
 double GetDinuc_GACG(){
     return dinuc_stat[3][8][6];
 }
 double GetDinuc_GACT(){
     return dinuc_stat[3][8][7];
 }
 double GetDinuc_GAGA(){
     return dinuc_stat[3][8][8];
 }
 double GetDinuc_GAGC(){
     return dinuc_stat[3][8][9];
 }
 double GetDinuc_GAGG(){
     return dinuc_stat[3][8][10];
 }
 double GetDinuc_GAGT(){
     return dinuc_stat[3][8][11];
 }
 double GetDinuc_GATA(){
     return dinuc_stat[3][8][12];
 }
 double GetDinuc_GATC(){
     return dinuc_stat[3][8][13];
 }
 double GetDinuc_GATG(){
     return dinuc_stat[3][8][14];
 }
 double GetDinuc_GATT(){
     return dinuc_stat[3][8][15];
 }
 double GetDinuc_GCAA(){
     return dinuc_stat[3][9][0];
 }
 double GetDinuc_GCAC(){
     return dinuc_stat[3][9][1];
 }
 double GetDinuc_GCAG(){
     return dinuc_stat[3][9][2];
 }
 double GetDinuc_GCAT(){
     return dinuc_stat[3][9][3];
 }
 double GetDinuc_GCCA(){
     return dinuc_stat[3][9][4];
 }
 double GetDinuc_GCCC(){
     return dinuc_stat[3][9][5];
 }
 double GetDinuc_GCCG(){
     return dinuc_stat[3][9][6];
 }
 double GetDinuc_GCCT(){
     return dinuc_stat[3][9][7];
 }
 double GetDinuc_GCGA(){
     return dinuc_stat[3][9][8];
 }
 double GetDinuc_GCGC(){
     return dinuc_stat[3][9][9];
 }
 double GetDinuc_GCGG(){
     return dinuc_stat[3][9][10];
 }
 double GetDinuc_GCGT(){
     return dinuc_stat[3][9][11];
 }
 double GetDinuc_GCTA(){
     return dinuc_stat[3][9][12];
 }
 double GetDinuc_GCTC(){
     return dinuc_stat[3][9][13];
 }
 double GetDinuc_GCTG(){
     return dinuc_stat[3][9][14];
 }
 double GetDinuc_GCTT(){
     return dinuc_stat[3][9][15];
 }
 double GetDinuc_GGAA(){
     return dinuc_stat[3][10][0];
 }
 double GetDinuc_GGAC(){
     return dinuc_stat[3][10][1];
 }
 double GetDinuc_GGAG(){
     return dinuc_stat[3][10][2];
 }
 double GetDinuc_GGAT(){
     return dinuc_stat[3][10][3];
 }
 double GetDinuc_GGCA(){
     return dinuc_stat[3][10][4];
 }
 double GetDinuc_GGCC(){
     return dinuc_stat[3][10][5];
 }
 double GetDinuc_GGCG(){
     return dinuc_stat[3][10][6];
 }
 double GetDinuc_GGCT(){
     return dinuc_stat[3][10][7];
 }
 double GetDinuc_GGGA(){
     return dinuc_stat[3][10][8];
 }
 double GetDinuc_GGGC(){
     return dinuc_stat[3][10][9];
 }
 double GetDinuc_GGGG(){
     return dinuc_stat[3][10][10];
 }
 double GetDinuc_GGGT(){
     return dinuc_stat[3][10][11];
 }
 double GetDinuc_GGTA(){
     return dinuc_stat[3][10][12];
 }
 double GetDinuc_GGTC(){
     return dinuc_stat[3][10][13];
 }
 double GetDinuc_GGTG(){
     return dinuc_stat[3][10][14];
 }
 double GetDinuc_GGTT(){
     return dinuc_stat[3][10][15];
 }
 double GetDinuc_GTAA(){
     return dinuc_stat[3][11][0];
 }
 double GetDinuc_GTAC(){
     return dinuc_stat[3][11][1];
 }
 double GetDinuc_GTAG(){
     return dinuc_stat[3][11][2];
 }
 double GetDinuc_GTAT(){
     return dinuc_stat[3][11][3];
 }
 double GetDinuc_GTCA(){
     return dinuc_stat[3][11][4];
 }
 double GetDinuc_GTCC(){
     return dinuc_stat[3][11][5];
 }
 double GetDinuc_GTCG(){
     return dinuc_stat[3][11][6];
 }
 double GetDinuc_GTCT(){
     return dinuc_stat[3][11][7];
 }
 double GetDinuc_GTGA(){
     return dinuc_stat[3][11][8];
 }
 double GetDinuc_GTGC(){
     return dinuc_stat[3][11][9];
 }
 double GetDinuc_GTGG(){
     return dinuc_stat[3][11][10];
 }
 double GetDinuc_GTGT(){
     return dinuc_stat[3][11][11];
 }
 double GetDinuc_GTTA(){
     return dinuc_stat[3][11][12];
 }
 double GetDinuc_GTTC(){
     return dinuc_stat[3][11][13];
 }
 double GetDinuc_GTTG(){
     return dinuc_stat[3][11][14];
 }
 double GetDinuc_GTTT(){
     return dinuc_stat[3][11][15];
 }
 double GetDinuc_TAAA(){
     return dinuc_stat[3][12][0];
 }
 double GetDinuc_TAAC(){
     return dinuc_stat[3][12][1];
 }
 double GetDinuc_TAAG(){
     return dinuc_stat[3][12][2];
 }
 double GetDinuc_TAAT(){
     return dinuc_stat[3][12][3];
 }
 double GetDinuc_TACA(){
     return dinuc_stat[3][12][4];
 }
 double GetDinuc_TACC(){
     return dinuc_stat[3][12][5];
 }
 double GetDinuc_TACG(){
     return dinuc_stat[3][12][6];
 }
 double GetDinuc_TACT(){
     return dinuc_stat[3][12][7];
 }
 double GetDinuc_TAGA(){
     return dinuc_stat[3][12][8];
 }
 double GetDinuc_TAGC(){
     return dinuc_stat[3][12][9];
 }
 double GetDinuc_TAGG(){
     return dinuc_stat[3][12][10];
 }
 double GetDinuc_TAGT(){
     return dinuc_stat[3][12][11];
 }
 double GetDinuc_TATA(){
     return dinuc_stat[3][12][12];
 }
 double GetDinuc_TATC(){
     return dinuc_stat[3][12][13];
 }
 double GetDinuc_TATG(){
     return dinuc_stat[3][12][14];
 }
 double GetDinuc_TATT(){
     return dinuc_stat[3][12][15];
 }
 double GetDinuc_TCAA(){
     return dinuc_stat[3][13][0];
 }
 double GetDinuc_TCAC(){
     return dinuc_stat[3][13][1];
 }
 double GetDinuc_TCAG(){
     return dinuc_stat[3][13][2];
 }
 double GetDinuc_TCAT(){
     return dinuc_stat[3][13][3];
 }
 double GetDinuc_TCCA(){
     return dinuc_stat[3][13][4];
 }
 double GetDinuc_TCCC(){
     return dinuc_stat[3][13][5];
 }
 double GetDinuc_TCCG(){
     return dinuc_stat[3][13][6];
 }
 double GetDinuc_TCCT(){
     return dinuc_stat[3][13][7];
 }
 double GetDinuc_TCGA(){
     return dinuc_stat[3][13][8];
 }
 double GetDinuc_TCGC(){
     return dinuc_stat[3][13][9];
 }
 double GetDinuc_TCGG(){
     return dinuc_stat[3][13][10];
 }
 double GetDinuc_TCGT(){
     return dinuc_stat[3][13][11];
 }
 double GetDinuc_TCTA(){
     return dinuc_stat[3][13][12];
 }
 double GetDinuc_TCTC(){
     return dinuc_stat[3][13][13];
 }
 double GetDinuc_TCTG(){
     return dinuc_stat[3][13][14];
 }
 double GetDinuc_TCTT(){
     return dinuc_stat[3][13][15];
 }
 double GetDinuc_TGAA(){
     return dinuc_stat[3][14][0];
 }
 double GetDinuc_TGAC(){
     return dinuc_stat[3][14][1];
 }
 double GetDinuc_TGAG(){
     return dinuc_stat[3][14][2];
 }
 double GetDinuc_TGAT(){
     return dinuc_stat[3][14][3];
 }
 double GetDinuc_TGCA(){
     return dinuc_stat[3][14][4];
 }
 double GetDinuc_TGCC(){
     return dinuc_stat[3][14][5];
 }
 double GetDinuc_TGCG(){
     return dinuc_stat[3][14][6];
 }
 double GetDinuc_TGCT(){
     return dinuc_stat[3][14][7];
 }
 double GetDinuc_TGGA(){
     return dinuc_stat[3][14][8];
 }
 double GetDinuc_TGGC(){
     return dinuc_stat[3][14][9];
 }
 double GetDinuc_TGGG(){
     return dinuc_stat[3][14][10];
 }
 double GetDinuc_TGGT(){
     return dinuc_stat[3][14][11];
 }
 double GetDinuc_TGTA(){
     return dinuc_stat[3][14][12];
 }
 double GetDinuc_TGTC(){
     return dinuc_stat[3][14][13];
 }
 double GetDinuc_TGTG(){
     return dinuc_stat[3][14][14];
 }
 double GetDinuc_TGTT(){
     return dinuc_stat[3][14][15];
 }
 double GetDinuc_TTAA(){
     return dinuc_stat[3][15][0];
 }
 double GetDinuc_TTAC(){
     return dinuc_stat[3][15][1];
 }
 double GetDinuc_TTAG(){
     return dinuc_stat[3][15][2];
 }
 double GetDinuc_TTAT(){
     return dinuc_stat[3][15][3];
 }
 double GetDinuc_TTCA(){
     return dinuc_stat[3][15][4];
 }
 double GetDinuc_TTCC(){
     return dinuc_stat[3][15][5];
 }
 double GetDinuc_TTCG(){
     return dinuc_stat[3][15][6];
 }
 double GetDinuc_TTCT(){
     return dinuc_stat[3][15][7];
 }
 double GetDinuc_TTGA(){
     return dinuc_stat[3][15][8];
 }
 double GetDinuc_TTGC(){
     return dinuc_stat[3][15][9];
 }
 double GetDinuc_TTGG(){
     return dinuc_stat[3][15][10];
 }
 double GetDinuc_TTGT(){
     return dinuc_stat[3][15][11];
 }
 double GetDinuc_TTTA(){
     return dinuc_stat[3][15][12];
 }
 double GetDinuc_TTTC(){
     return dinuc_stat[3][15][13];
 }
 double GetDinuc_TTTG(){
     return dinuc_stat[3][15][14];
 }
 double GetDinuc_TTTT(){
     return dinuc_stat[3][15][15];
 }
 double GetDinuc12_AAAA(){
     return dinuc_stat[0][0][0];
 }
 double GetDinuc12_AAAC(){
     return dinuc_stat[0][0][1];
 }
 double GetDinuc12_AAAG(){
     return dinuc_stat[0][0][2];
 }
 double GetDinuc12_AAAT(){
     return dinuc_stat[0][0][3];
 }
 double GetDinuc12_AACA(){
     return dinuc_stat[0][0][4];
 }
 double GetDinuc12_AACC(){
     return dinuc_stat[0][0][5];
 }
 double GetDinuc12_AACG(){
     return dinuc_stat[0][0][6];
 }
 double GetDinuc12_AACT(){
     return dinuc_stat[0][0][7];
 }
 double GetDinuc12_AAGA(){
     return dinuc_stat[0][0][8];
 }
 double GetDinuc12_AAGC(){
     return dinuc_stat[0][0][9];
 }
 double GetDinuc12_AAGG(){
     return dinuc_stat[0][0][10];
 }
 double GetDinuc12_AAGT(){
     return dinuc_stat[0][0][11];
 }
 double GetDinuc12_AATA(){
     return dinuc_stat[0][0][12];
 }
 double GetDinuc12_AATC(){
     return dinuc_stat[0][0][13];
 }
 double GetDinuc12_AATG(){
     return dinuc_stat[0][0][14];
 }
 double GetDinuc12_AATT(){
     return dinuc_stat[0][0][15];
 }
 double GetDinuc12_ACAA(){
     return dinuc_stat[0][1][0];
 }
 double GetDinuc12_ACAC(){
     return dinuc_stat[0][1][1];
 }
 double GetDinuc12_ACAG(){
     return dinuc_stat[0][1][2];
 }
 double GetDinuc12_ACAT(){
     return dinuc_stat[0][1][3];
 }
 double GetDinuc12_ACCA(){
     return dinuc_stat[0][1][4];
 }
 double GetDinuc12_ACCC(){
     return dinuc_stat[0][1][5];
 }
 double GetDinuc12_ACCG(){
     return dinuc_stat[0][1][6];
 }
 double GetDinuc12_ACCT(){
     return dinuc_stat[0][1][7];
 }
 double GetDinuc12_ACGA(){
     return dinuc_stat[0][1][8];
 }
 double GetDinuc12_ACGC(){
     return dinuc_stat[0][1][9];
 }
 double GetDinuc12_ACGG(){
     return dinuc_stat[0][1][10];
 }
 double GetDinuc12_ACGT(){
     return dinuc_stat[0][1][11];
 }
 double GetDinuc12_ACTA(){
     return dinuc_stat[0][1][12];
 }
 double GetDinuc12_ACTC(){
     return dinuc_stat[0][1][13];
 }
 double GetDinuc12_ACTG(){
     return dinuc_stat[0][1][14];
 }
 double GetDinuc12_ACTT(){
     return dinuc_stat[0][1][15];
 }
 double GetDinuc12_AGAA(){
     return dinuc_stat[0][2][0];
 }
 double GetDinuc12_AGAC(){
     return dinuc_stat[0][2][1];
 }
 double GetDinuc12_AGAG(){
     return dinuc_stat[0][2][2];
 }
 double GetDinuc12_AGAT(){
     return dinuc_stat[0][2][3];
 }
 double GetDinuc12_AGCA(){
     return dinuc_stat[0][2][4];
 }
 double GetDinuc12_AGCC(){
     return dinuc_stat[0][2][5];
 }
 double GetDinuc12_AGCG(){
     return dinuc_stat[0][2][6];
 }
 double GetDinuc12_AGCT(){
     return dinuc_stat[0][2][7];
 }
 double GetDinuc12_AGGA(){
     return dinuc_stat[0][2][8];
 }
 double GetDinuc12_AGGC(){
     return dinuc_stat[0][2][9];
 }
 double GetDinuc12_AGGG(){
     return dinuc_stat[0][2][10];
 }
 double GetDinuc12_AGGT(){
     return dinuc_stat[0][2][11];
 }
 double GetDinuc12_AGTA(){
     return dinuc_stat[0][2][12];
 }
 double GetDinuc12_AGTC(){
     return dinuc_stat[0][2][13];
 }
 double GetDinuc12_AGTG(){
     return dinuc_stat[0][2][14];
 }
 double GetDinuc12_AGTT(){
     return dinuc_stat[0][2][15];
 }
 double GetDinuc12_ATAA(){
     return dinuc_stat[0][3][0];
 }
 double GetDinuc12_ATAC(){
     return dinuc_stat[0][3][1];
 }
 double GetDinuc12_ATAG(){
     return dinuc_stat[0][3][2];
 }
 double GetDinuc12_ATAT(){
     return dinuc_stat[0][3][3];
 }
 double GetDinuc12_ATCA(){
     return dinuc_stat[0][3][4];
 }
 double GetDinuc12_ATCC(){
     return dinuc_stat[0][3][5];
 }
 double GetDinuc12_ATCG(){
     return dinuc_stat[0][3][6];
 }
 double GetDinuc12_ATCT(){
     return dinuc_stat[0][3][7];
 }
 double GetDinuc12_ATGA(){
     return dinuc_stat[0][3][8];
 }
 double GetDinuc12_ATGC(){
     return dinuc_stat[0][3][9];
 }
 double GetDinuc12_ATGG(){
     return dinuc_stat[0][3][10];
 }
 double GetDinuc12_ATGT(){
     return dinuc_stat[0][3][11];
 }
 double GetDinuc12_ATTA(){
     return dinuc_stat[0][3][12];
 }
 double GetDinuc12_ATTC(){
     return dinuc_stat[0][3][13];
 }
 double GetDinuc12_ATTG(){
     return dinuc_stat[0][3][14];
 }
 double GetDinuc12_ATTT(){
     return dinuc_stat[0][3][15];
 }
 double GetDinuc12_CAAA(){
     return dinuc_stat[0][4][0];
 }
 double GetDinuc12_CAAC(){
     return dinuc_stat[0][4][1];
 }
 double GetDinuc12_CAAG(){
     return dinuc_stat[0][4][2];
 }
 double GetDinuc12_CAAT(){
     return dinuc_stat[0][4][3];
 }
 double GetDinuc12_CACA(){
     return dinuc_stat[0][4][4];
 }
 double GetDinuc12_CACC(){
     return dinuc_stat[0][4][5];
 }
 double GetDinuc12_CACG(){
     return dinuc_stat[0][4][6];
 }
 double GetDinuc12_CACT(){
     return dinuc_stat[0][4][7];
 }
 double GetDinuc12_CAGA(){
     return dinuc_stat[0][4][8];
 }
 double GetDinuc12_CAGC(){
     return dinuc_stat[0][4][9];
 }
 double GetDinuc12_CAGG(){
     return dinuc_stat[0][4][10];
 }
 double GetDinuc12_CAGT(){
     return dinuc_stat[0][4][11];
 }
 double GetDinuc12_CATA(){
     return dinuc_stat[0][4][12];
 }
 double GetDinuc12_CATC(){
     return dinuc_stat[0][4][13];
 }
 double GetDinuc12_CATG(){
     return dinuc_stat[0][4][14];
 }
 double GetDinuc12_CATT(){
     return dinuc_stat[0][4][15];
 }
 double GetDinuc12_CCAA(){
     return dinuc_stat[0][5][0];
 }
 double GetDinuc12_CCAC(){
     return dinuc_stat[0][5][1];
 }
 double GetDinuc12_CCAG(){
     return dinuc_stat[0][5][2];
 }
 double GetDinuc12_CCAT(){
     return dinuc_stat[0][5][3];
 }
 double GetDinuc12_CCCA(){
     return dinuc_stat[0][5][4];
 }
 double GetDinuc12_CCCC(){
     return dinuc_stat[0][5][5];
 }
 double GetDinuc12_CCCG(){
     return dinuc_stat[0][5][6];
 }
 double GetDinuc12_CCCT(){
     return dinuc_stat[0][5][7];
 }
 double GetDinuc12_CCGA(){
     return dinuc_stat[0][5][8];
 }
 double GetDinuc12_CCGC(){
     return dinuc_stat[0][5][9];
 }
 double GetDinuc12_CCGG(){
     return dinuc_stat[0][5][10];
 }
 double GetDinuc12_CCGT(){
     return dinuc_stat[0][5][11];
 }
 double GetDinuc12_CCTA(){
     return dinuc_stat[0][5][12];
 }
 double GetDinuc12_CCTC(){
     return dinuc_stat[0][5][13];
 }
 double GetDinuc12_CCTG(){
     return dinuc_stat[0][5][14];
 }
 double GetDinuc12_CCTT(){
     return dinuc_stat[0][5][15];
 }
 double GetDinuc12_CGAA(){
     return dinuc_stat[0][6][0];
 }
 double GetDinuc12_CGAC(){
     return dinuc_stat[0][6][1];
 }
 double GetDinuc12_CGAG(){
     return dinuc_stat[0][6][2];
 }
 double GetDinuc12_CGAT(){
     return dinuc_stat[0][6][3];
 }
 double GetDinuc12_CGCA(){
     return dinuc_stat[0][6][4];
 }
 double GetDinuc12_CGCC(){
     return dinuc_stat[0][6][5];
 }
 double GetDinuc12_CGCG(){
     return dinuc_stat[0][6][6];
 }
 double GetDinuc12_CGCT(){
     return dinuc_stat[0][6][7];
 }
 double GetDinuc12_CGGA(){
     return dinuc_stat[0][6][8];
 }
 double GetDinuc12_CGGC(){
     return dinuc_stat[0][6][9];
 }
 double GetDinuc12_CGGG(){
     return dinuc_stat[0][6][10];
 }
 double GetDinuc12_CGGT(){
     return dinuc_stat[0][6][11];
 }
 double GetDinuc12_CGTA(){
     return dinuc_stat[0][6][12];
 }
 double GetDinuc12_CGTC(){
     return dinuc_stat[0][6][13];
 }
 double GetDinuc12_CGTG(){
     return dinuc_stat[0][6][14];
 }
 double GetDinuc12_CGTT(){
     return dinuc_stat[0][6][15];
 }
 double GetDinuc12_CTAA(){
     return dinuc_stat[0][7][0];
 }
 double GetDinuc12_CTAC(){
     return dinuc_stat[0][7][1];
 }
 double GetDinuc12_CTAG(){
     return dinuc_stat[0][7][2];
 }
 double GetDinuc12_CTAT(){
     return dinuc_stat[0][7][3];
 }
 double GetDinuc12_CTCA(){
     return dinuc_stat[0][7][4];
 }
 double GetDinuc12_CTCC(){
     return dinuc_stat[0][7][5];
 }
 double GetDinuc12_CTCG(){
     return dinuc_stat[0][7][6];
 }
 double GetDinuc12_CTCT(){
     return dinuc_stat[0][7][7];
 }
 double GetDinuc12_CTGA(){
     return dinuc_stat[0][7][8];
 }
 double GetDinuc12_CTGC(){
     return dinuc_stat[0][7][9];
 }
 double GetDinuc12_CTGG(){
     return dinuc_stat[0][7][10];
 }
 double GetDinuc12_CTGT(){
     return dinuc_stat[0][7][11];
 }
 double GetDinuc12_CTTA(){
     return dinuc_stat[0][7][12];
 }
 double GetDinuc12_CTTC(){
     return dinuc_stat[0][7][13];
 }
 double GetDinuc12_CTTG(){
     return dinuc_stat[0][7][14];
 }
 double GetDinuc12_CTTT(){
     return dinuc_stat[0][7][15];
 }
 double GetDinuc12_GAAA(){
     return dinuc_stat[0][8][0];
 }
 double GetDinuc12_GAAC(){
     return dinuc_stat[0][8][1];
 }
 double GetDinuc12_GAAG(){
     return dinuc_stat[0][8][2];
 }
 double GetDinuc12_GAAT(){
     return dinuc_stat[0][8][3];
 }
 double GetDinuc12_GACA(){
     return dinuc_stat[0][8][4];
 }
 double GetDinuc12_GACC(){
     return dinuc_stat[0][8][5];
 }
 double GetDinuc12_GACG(){
     return dinuc_stat[0][8][6];
 }
 double GetDinuc12_GACT(){
     return dinuc_stat[0][8][7];
 }
 double GetDinuc12_GAGA(){
     return dinuc_stat[0][8][8];
 }
 double GetDinuc12_GAGC(){
     return dinuc_stat[0][8][9];
 }
 double GetDinuc12_GAGG(){
     return dinuc_stat[0][8][10];
 }
 double GetDinuc12_GAGT(){
     return dinuc_stat[0][8][11];
 }
 double GetDinuc12_GATA(){
     return dinuc_stat[0][8][12];
 }
 double GetDinuc12_GATC(){
     return dinuc_stat[0][8][13];
 }
 double GetDinuc12_GATG(){
     return dinuc_stat[0][8][14];
 }
 double GetDinuc12_GATT(){
     return dinuc_stat[0][8][15];
 }
 double GetDinuc12_GCAA(){
     return dinuc_stat[0][9][0];
 }
 double GetDinuc12_GCAC(){
     return dinuc_stat[0][9][1];
 }
 double GetDinuc12_GCAG(){
     return dinuc_stat[0][9][2];
 }
 double GetDinuc12_GCAT(){
     return dinuc_stat[0][9][3];
 }
 double GetDinuc12_GCCA(){
     return dinuc_stat[0][9][4];
 }
 double GetDinuc12_GCCC(){
     return dinuc_stat[0][9][5];
 }
 double GetDinuc12_GCCG(){
     return dinuc_stat[0][9][6];
 }
 double GetDinuc12_GCCT(){
     return dinuc_stat[0][9][7];
 }
 double GetDinuc12_GCGA(){
     return dinuc_stat[0][9][8];
 }
 double GetDinuc12_GCGC(){
     return dinuc_stat[0][9][9];
 }
 double GetDinuc12_GCGG(){
     return dinuc_stat[0][9][10];
 }
 double GetDinuc12_GCGT(){
     return dinuc_stat[0][9][11];
 }
 double GetDinuc12_GCTA(){
     return dinuc_stat[0][9][12];
 }
 double GetDinuc12_GCTC(){
     return dinuc_stat[0][9][13];
 }
 double GetDinuc12_GCTG(){
     return dinuc_stat[0][9][14];
 }
 double GetDinuc12_GCTT(){
     return dinuc_stat[0][9][15];
 }
 double GetDinuc12_GGAA(){
     return dinuc_stat[0][10][0];
 }
 double GetDinuc12_GGAC(){
     return dinuc_stat[0][10][1];
 }
 double GetDinuc12_GGAG(){
     return dinuc_stat[0][10][2];
 }
 double GetDinuc12_GGAT(){
     return dinuc_stat[0][10][3];
 }
 double GetDinuc12_GGCA(){
     return dinuc_stat[0][10][4];
 }
 double GetDinuc12_GGCC(){
     return dinuc_stat[0][10][5];
 }
 double GetDinuc12_GGCG(){
     return dinuc_stat[0][10][6];
 }
 double GetDinuc12_GGCT(){
     return dinuc_stat[0][10][7];
 }
 double GetDinuc12_GGGA(){
     return dinuc_stat[0][10][8];
 }
 double GetDinuc12_GGGC(){
     return dinuc_stat[0][10][9];
 }
 double GetDinuc12_GGGG(){
     return dinuc_stat[0][10][10];
 }
 double GetDinuc12_GGGT(){
     return dinuc_stat[0][10][11];
 }
 double GetDinuc12_GGTA(){
     return dinuc_stat[0][10][12];
 }
 double GetDinuc12_GGTC(){
     return dinuc_stat[0][10][13];
 }
 double GetDinuc12_GGTG(){
     return dinuc_stat[0][10][14];
 }
 double GetDinuc12_GGTT(){
     return dinuc_stat[0][10][15];
 }
 double GetDinuc12_GTAA(){
     return dinuc_stat[0][11][0];
 }
 double GetDinuc12_GTAC(){
     return dinuc_stat[0][11][1];
 }
 double GetDinuc12_GTAG(){
     return dinuc_stat[0][11][2];
 }
 double GetDinuc12_GTAT(){
     return dinuc_stat[0][11][3];
 }
 double GetDinuc12_GTCA(){
     return dinuc_stat[0][11][4];
 }
 double GetDinuc12_GTCC(){
     return dinuc_stat[0][11][5];
 }
 double GetDinuc12_GTCG(){
     return dinuc_stat[0][11][6];
 }
 double GetDinuc12_GTCT(){
     return dinuc_stat[0][11][7];
 }
 double GetDinuc12_GTGA(){
     return dinuc_stat[0][11][8];
 }
 double GetDinuc12_GTGC(){
     return dinuc_stat[0][11][9];
 }
 double GetDinuc12_GTGG(){
     return dinuc_stat[0][11][10];
 }
 double GetDinuc12_GTGT(){
     return dinuc_stat[0][11][11];
 }
 double GetDinuc12_GTTA(){
     return dinuc_stat[0][11][12];
 }
 double GetDinuc12_GTTC(){
     return dinuc_stat[0][11][13];
 }
 double GetDinuc12_GTTG(){
     return dinuc_stat[0][11][14];
 }
 double GetDinuc12_GTTT(){
     return dinuc_stat[0][11][15];
 }
 double GetDinuc12_TAAA(){
     return dinuc_stat[0][12][0];
 }
 double GetDinuc12_TAAC(){
     return dinuc_stat[0][12][1];
 }
 double GetDinuc12_TAAG(){
     return dinuc_stat[0][12][2];
 }
 double GetDinuc12_TAAT(){
     return dinuc_stat[0][12][3];
 }
 double GetDinuc12_TACA(){
     return dinuc_stat[0][12][4];
 }
 double GetDinuc12_TACC(){
     return dinuc_stat[0][12][5];
 }
 double GetDinuc12_TACG(){
     return dinuc_stat[0][12][6];
 }
 double GetDinuc12_TACT(){
     return dinuc_stat[0][12][7];
 }
 double GetDinuc12_TAGA(){
     return dinuc_stat[0][12][8];
 }
 double GetDinuc12_TAGC(){
     return dinuc_stat[0][12][9];
 }
 double GetDinuc12_TAGG(){
     return dinuc_stat[0][12][10];
 }
 double GetDinuc12_TAGT(){
     return dinuc_stat[0][12][11];
 }
 double GetDinuc12_TATA(){
     return dinuc_stat[0][12][12];
 }
 double GetDinuc12_TATC(){
     return dinuc_stat[0][12][13];
 }
 double GetDinuc12_TATG(){
     return dinuc_stat[0][12][14];
 }
 double GetDinuc12_TATT(){
     return dinuc_stat[0][12][15];
 }
 double GetDinuc12_TCAA(){
     return dinuc_stat[0][13][0];
 }
 double GetDinuc12_TCAC(){
     return dinuc_stat[0][13][1];
 }
 double GetDinuc12_TCAG(){
     return dinuc_stat[0][13][2];
 }
 double GetDinuc12_TCAT(){
     return dinuc_stat[0][13][3];
 }
 double GetDinuc12_TCCA(){
     return dinuc_stat[0][13][4];
 }
 double GetDinuc12_TCCC(){
     return dinuc_stat[0][13][5];
 }
 double GetDinuc12_TCCG(){
     return dinuc_stat[0][13][6];
 }
 double GetDinuc12_TCCT(){
     return dinuc_stat[0][13][7];
 }
 double GetDinuc12_TCGA(){
     return dinuc_stat[0][13][8];
 }
 double GetDinuc12_TCGC(){
     return dinuc_stat[0][13][9];
 }
 double GetDinuc12_TCGG(){
     return dinuc_stat[0][13][10];
 }
 double GetDinuc12_TCGT(){
     return dinuc_stat[0][13][11];
 }
 double GetDinuc12_TCTA(){
     return dinuc_stat[0][13][12];
 }
 double GetDinuc12_TCTC(){
     return dinuc_stat[0][13][13];
 }
 double GetDinuc12_TCTG(){
     return dinuc_stat[0][13][14];
 }
 double GetDinuc12_TCTT(){
     return dinuc_stat[0][13][15];
 }
 double GetDinuc12_TGAA(){
     return dinuc_stat[0][14][0];
 }
 double GetDinuc12_TGAC(){
     return dinuc_stat[0][14][1];
 }
 double GetDinuc12_TGAG(){
     return dinuc_stat[0][14][2];
 }
 double GetDinuc12_TGAT(){
     return dinuc_stat[0][14][3];
 }
 double GetDinuc12_TGCA(){
     return dinuc_stat[0][14][4];
 }
 double GetDinuc12_TGCC(){
     return dinuc_stat[0][14][5];
 }
 double GetDinuc12_TGCG(){
     return dinuc_stat[0][14][6];
 }
 double GetDinuc12_TGCT(){
     return dinuc_stat[0][14][7];
 }
 double GetDinuc12_TGGA(){
     return dinuc_stat[0][14][8];
 }
 double GetDinuc12_TGGC(){
     return dinuc_stat[0][14][9];
 }
 double GetDinuc12_TGGG(){
     return dinuc_stat[0][14][10];
 }
 double GetDinuc12_TGGT(){
     return dinuc_stat[0][14][11];
 }
 double GetDinuc12_TGTA(){
     return dinuc_stat[0][14][12];
 }
 double GetDinuc12_TGTC(){
     return dinuc_stat[0][14][13];
 }
 double GetDinuc12_TGTG(){
     return dinuc_stat[0][14][14];
 }
 double GetDinuc12_TGTT(){
     return dinuc_stat[0][14][15];
 }
 double GetDinuc12_TTAA(){
     return dinuc_stat[0][15][0];
 }
 double GetDinuc12_TTAC(){
     return dinuc_stat[0][15][1];
 }
 double GetDinuc12_TTAG(){
     return dinuc_stat[0][15][2];
 }
 double GetDinuc12_TTAT(){
     return dinuc_stat[0][15][3];
 }
 double GetDinuc12_TTCA(){
     return dinuc_stat[0][15][4];
 }
 double GetDinuc12_TTCC(){
     return dinuc_stat[0][15][5];
 }
 double GetDinuc12_TTCG(){
     return dinuc_stat[0][15][6];
 }
 double GetDinuc12_TTCT(){
     return dinuc_stat[0][15][7];
 }
 double GetDinuc12_TTGA(){
     return dinuc_stat[0][15][8];
 }
 double GetDinuc12_TTGC(){
     return dinuc_stat[0][15][9];
 }
 double GetDinuc12_TTGG(){
     return dinuc_stat[0][15][10];
 }
 double GetDinuc12_TTGT(){
     return dinuc_stat[0][15][11];
 }
 double GetDinuc12_TTTA(){
     return dinuc_stat[0][15][12];
 }
 double GetDinuc12_TTTC(){
     return dinuc_stat[0][15][13];
 }
 double GetDinuc12_TTTG(){
     return dinuc_stat[0][15][14];
 }
 double GetDinuc12_TTTT(){
     return dinuc_stat[0][15][15];
 }

 double GetDinuc23_AAAA(){
     return dinuc_stat[1][0][0];
 }
 double GetDinuc23_AAAC(){
     return dinuc_stat[1][0][1];
 }
 double GetDinuc23_AAAG(){
     return dinuc_stat[1][0][2];
 }
 double GetDinuc23_AAAT(){
     return dinuc_stat[1][0][3];
 }
 double GetDinuc23_AACA(){
     return dinuc_stat[1][0][4];
 }
 double GetDinuc23_AACC(){
     return dinuc_stat[1][0][5];
 }
 double GetDinuc23_AACG(){
     return dinuc_stat[1][0][6];
 }
 double GetDinuc23_AACT(){
     return dinuc_stat[1][0][7];
 }
 double GetDinuc23_AAGA(){
     return dinuc_stat[1][0][8];
 }
 double GetDinuc23_AAGC(){
     return dinuc_stat[1][0][9];
 }
 double GetDinuc23_AAGG(){
     return dinuc_stat[1][0][10];
 }
 double GetDinuc23_AAGT(){
     return dinuc_stat[1][0][11];
 }
 double GetDinuc23_AATA(){
     return dinuc_stat[1][0][12];
 }
 double GetDinuc23_AATC(){
     return dinuc_stat[1][0][13];
 }
 double GetDinuc23_AATG(){
     return dinuc_stat[1][0][14];
 }
 double GetDinuc23_AATT(){
     return dinuc_stat[1][0][15];
 }
 double GetDinuc23_ACAA(){
     return dinuc_stat[1][1][0];
 }
 double GetDinuc23_ACAC(){
     return dinuc_stat[1][1][1];
 }
 double GetDinuc23_ACAG(){
     return dinuc_stat[1][1][2];
 }
 double GetDinuc23_ACAT(){
     return dinuc_stat[1][1][3];
 }
 double GetDinuc23_ACCA(){
     return dinuc_stat[1][1][4];
 }
 double GetDinuc23_ACCC(){
     return dinuc_stat[1][1][5];
 }
 double GetDinuc23_ACCG(){
     return dinuc_stat[1][1][6];
 }
 double GetDinuc23_ACCT(){
     return dinuc_stat[1][1][7];
 }
 double GetDinuc23_ACGA(){
     return dinuc_stat[1][1][8];
 }
 double GetDinuc23_ACGC(){
     return dinuc_stat[1][1][9];
 }
 double GetDinuc23_ACGG(){
     return dinuc_stat[1][1][10];
 }
 double GetDinuc23_ACGT(){
     return dinuc_stat[1][1][11];
 }
 double GetDinuc23_ACTA(){
     return dinuc_stat[1][1][12];
 }
 double GetDinuc23_ACTC(){
     return dinuc_stat[1][1][13];
 }
 double GetDinuc23_ACTG(){
     return dinuc_stat[1][1][14];
 }
 double GetDinuc23_ACTT(){
     return dinuc_stat[1][1][15];
 }
 double GetDinuc23_AGAA(){
     return dinuc_stat[1][2][0];
 }
 double GetDinuc23_AGAC(){
     return dinuc_stat[1][2][1];
 }
 double GetDinuc23_AGAG(){
     return dinuc_stat[1][2][2];
 }
 double GetDinuc23_AGAT(){
     return dinuc_stat[1][2][3];
 }
 double GetDinuc23_AGCA(){
     return dinuc_stat[1][2][4];
 }
 double GetDinuc23_AGCC(){
     return dinuc_stat[1][2][5];
 }
 double GetDinuc23_AGCG(){
     return dinuc_stat[1][2][6];
 }
 double GetDinuc23_AGCT(){
     return dinuc_stat[1][2][7];
 }
 double GetDinuc23_AGGA(){
     return dinuc_stat[1][2][8];
 }
 double GetDinuc23_AGGC(){
     return dinuc_stat[1][2][9];
 }
 double GetDinuc23_AGGG(){
     return dinuc_stat[1][2][10];
 }
 double GetDinuc23_AGGT(){
     return dinuc_stat[1][2][11];
 }
 double GetDinuc23_AGTA(){
     return dinuc_stat[1][2][12];
 }
 double GetDinuc23_AGTC(){
     return dinuc_stat[1][2][13];
 }
 double GetDinuc23_AGTG(){
     return dinuc_stat[1][2][14];
 }
 double GetDinuc23_AGTT(){
     return dinuc_stat[1][2][15];
 }
 double GetDinuc23_ATAA(){
     return dinuc_stat[1][3][0];
 }
 double GetDinuc23_ATAC(){
     return dinuc_stat[1][3][1];
 }
 double GetDinuc23_ATAG(){
     return dinuc_stat[1][3][2];
 }
 double GetDinuc23_ATAT(){
     return dinuc_stat[1][3][3];
 }
 double GetDinuc23_ATCA(){
     return dinuc_stat[1][3][4];
 }
 double GetDinuc23_ATCC(){
     return dinuc_stat[1][3][5];
 }
 double GetDinuc23_ATCG(){
     return dinuc_stat[1][3][6];
 }
 double GetDinuc23_ATCT(){
     return dinuc_stat[1][3][7];
 }
 double GetDinuc23_ATGA(){
     return dinuc_stat[1][3][8];
 }
 double GetDinuc23_ATGC(){
     return dinuc_stat[1][3][9];
 }
 double GetDinuc23_ATGG(){
     return dinuc_stat[1][3][10];
 }
 double GetDinuc23_ATGT(){
     return dinuc_stat[1][3][11];
 }
 double GetDinuc23_ATTA(){
     return dinuc_stat[1][3][12];
 }
 double GetDinuc23_ATTC(){
     return dinuc_stat[1][3][13];
 }
 double GetDinuc23_ATTG(){
     return dinuc_stat[1][3][14];
 }
 double GetDinuc23_ATTT(){
     return dinuc_stat[1][3][15];
 }
 double GetDinuc23_CAAA(){
     return dinuc_stat[1][4][0];
 }
 double GetDinuc23_CAAC(){
     return dinuc_stat[1][4][1];
 }
 double GetDinuc23_CAAG(){
     return dinuc_stat[1][4][2];
 }
 double GetDinuc23_CAAT(){
     return dinuc_stat[1][4][3];
 }
 double GetDinuc23_CACA(){
     return dinuc_stat[1][4][4];
 }
 double GetDinuc23_CACC(){
     return dinuc_stat[1][4][5];
 }
 double GetDinuc23_CACG(){
     return dinuc_stat[1][4][6];
 }
 double GetDinuc23_CACT(){
     return dinuc_stat[1][4][7];
 }
 double GetDinuc23_CAGA(){
     return dinuc_stat[1][4][8];
 }
 double GetDinuc23_CAGC(){
     return dinuc_stat[1][4][9];
 }
 double GetDinuc23_CAGG(){
     return dinuc_stat[1][4][10];
 }
 double GetDinuc23_CAGT(){
     return dinuc_stat[1][4][11];
 }
 double GetDinuc23_CATA(){
     return dinuc_stat[1][4][12];
 }
 double GetDinuc23_CATC(){
     return dinuc_stat[1][4][13];
 }
 double GetDinuc23_CATG(){
     return dinuc_stat[1][4][14];
 }
 double GetDinuc23_CATT(){
     return dinuc_stat[1][4][15];
 }
 double GetDinuc23_CCAA(){
     return dinuc_stat[1][5][0];
 }
 double GetDinuc23_CCAC(){
     return dinuc_stat[1][5][1];
 }
 double GetDinuc23_CCAG(){
     return dinuc_stat[1][5][2];
 }
 double GetDinuc23_CCAT(){
     return dinuc_stat[1][5][3];
 }
 double GetDinuc23_CCCA(){
     return dinuc_stat[1][5][4];
 }
 double GetDinuc23_CCCC(){
     return dinuc_stat[1][5][5];
 }
 double GetDinuc23_CCCG(){
     return dinuc_stat[1][5][6];
 }
 double GetDinuc23_CCCT(){
     return dinuc_stat[1][5][7];
 }
 double GetDinuc23_CCGA(){
     return dinuc_stat[1][5][8];
 }
 double GetDinuc23_CCGC(){
     return dinuc_stat[1][5][9];
 }
 double GetDinuc23_CCGG(){
     return dinuc_stat[1][5][10];
 }
 double GetDinuc23_CCGT(){
     return dinuc_stat[1][5][11];
 }
 double GetDinuc23_CCTA(){
     return dinuc_stat[1][5][12];
 }
 double GetDinuc23_CCTC(){
     return dinuc_stat[1][5][13];
 }
 double GetDinuc23_CCTG(){
     return dinuc_stat[1][5][14];
 }
 double GetDinuc23_CCTT(){
     return dinuc_stat[1][5][15];
 }
 double GetDinuc23_CGAA(){
     return dinuc_stat[1][6][0];
 }
 double GetDinuc23_CGAC(){
     return dinuc_stat[1][6][1];
 }
 double GetDinuc23_CGAG(){
     return dinuc_stat[1][6][2];
 }
 double GetDinuc23_CGAT(){
     return dinuc_stat[1][6][3];
 }
 double GetDinuc23_CGCA(){
     return dinuc_stat[1][6][4];
 }
 double GetDinuc23_CGCC(){
     return dinuc_stat[1][6][5];
 }
 double GetDinuc23_CGCG(){
     return dinuc_stat[1][6][6];
 }
 double GetDinuc23_CGCT(){
     return dinuc_stat[1][6][7];
 }
 double GetDinuc23_CGGA(){
     return dinuc_stat[1][6][8];
 }
 double GetDinuc23_CGGC(){
     return dinuc_stat[1][6][9];
 }
 double GetDinuc23_CGGG(){
     return dinuc_stat[1][6][10];
 }
 double GetDinuc23_CGGT(){
     return dinuc_stat[1][6][11];
 }
 double GetDinuc23_CGTA(){
     return dinuc_stat[1][6][12];
 }
 double GetDinuc23_CGTC(){
     return dinuc_stat[1][6][13];
 }
 double GetDinuc23_CGTG(){
     return dinuc_stat[1][6][14];
 }
 double GetDinuc23_CGTT(){
     return dinuc_stat[1][6][15];
 }
 double GetDinuc23_CTAA(){
     return dinuc_stat[1][7][0];
 }
 double GetDinuc23_CTAC(){
     return dinuc_stat[1][7][1];
 }
 double GetDinuc23_CTAG(){
     return dinuc_stat[1][7][2];
 }
 double GetDinuc23_CTAT(){
     return dinuc_stat[1][7][3];
 }
 double GetDinuc23_CTCA(){
     return dinuc_stat[1][7][4];
 }
 double GetDinuc23_CTCC(){
     return dinuc_stat[1][7][5];
 }
 double GetDinuc23_CTCG(){
     return dinuc_stat[1][7][6];
 }
 double GetDinuc23_CTCT(){
     return dinuc_stat[1][7][7];
 }
 double GetDinuc23_CTGA(){
     return dinuc_stat[1][7][8];
 }
 double GetDinuc23_CTGC(){
     return dinuc_stat[1][7][9];
 }
 double GetDinuc23_CTGG(){
     return dinuc_stat[1][7][10];
 }
 double GetDinuc23_CTGT(){
     return dinuc_stat[1][7][11];
 }
 double GetDinuc23_CTTA(){
     return dinuc_stat[1][7][12];
 }
 double GetDinuc23_CTTC(){
     return dinuc_stat[1][7][13];
 }
 double GetDinuc23_CTTG(){
     return dinuc_stat[1][7][14];
 }
 double GetDinuc23_CTTT(){
     return dinuc_stat[1][7][15];
 }
 double GetDinuc23_GAAA(){
     return dinuc_stat[1][8][0];
 }
 double GetDinuc23_GAAC(){
     return dinuc_stat[1][8][1];
 }
 double GetDinuc23_GAAG(){
     return dinuc_stat[1][8][2];
 }
 double GetDinuc23_GAAT(){
     return dinuc_stat[1][8][3];
 }
 double GetDinuc23_GACA(){
     return dinuc_stat[1][8][4];
 }
 double GetDinuc23_GACC(){
     return dinuc_stat[1][8][5];
 }
 double GetDinuc23_GACG(){
     return dinuc_stat[1][8][6];
 }
 double GetDinuc23_GACT(){
     return dinuc_stat[1][8][7];
 }
 double GetDinuc23_GAGA(){
     return dinuc_stat[1][8][8];
 }
 double GetDinuc23_GAGC(){
     return dinuc_stat[1][8][9];
 }
 double GetDinuc23_GAGG(){
     return dinuc_stat[1][8][10];
 }
 double GetDinuc23_GAGT(){
     return dinuc_stat[1][8][11];
 }
 double GetDinuc23_GATA(){
     return dinuc_stat[1][8][12];
 }
 double GetDinuc23_GATC(){
     return dinuc_stat[1][8][13];
 }
 double GetDinuc23_GATG(){
     return dinuc_stat[1][8][14];
 }
 double GetDinuc23_GATT(){
     return dinuc_stat[1][8][15];
 }
 double GetDinuc23_GCAA(){
     return dinuc_stat[1][9][0];
 }
 double GetDinuc23_GCAC(){
     return dinuc_stat[1][9][1];
 }
 double GetDinuc23_GCAG(){
     return dinuc_stat[1][9][2];
 }
 double GetDinuc23_GCAT(){
     return dinuc_stat[1][9][3];
 }
 double GetDinuc23_GCCA(){
     return dinuc_stat[1][9][4];
 }
 double GetDinuc23_GCCC(){
     return dinuc_stat[1][9][5];
 }
 double GetDinuc23_GCCG(){
     return dinuc_stat[1][9][6];
 }
 double GetDinuc23_GCCT(){
     return dinuc_stat[1][9][7];
 }
 double GetDinuc23_GCGA(){
     return dinuc_stat[1][9][8];
 }
 double GetDinuc23_GCGC(){
     return dinuc_stat[1][9][9];
 }
 double GetDinuc23_GCGG(){
     return dinuc_stat[1][9][10];
 }
 double GetDinuc23_GCGT(){
     return dinuc_stat[1][9][11];
 }
 double GetDinuc23_GCTA(){
     return dinuc_stat[1][9][12];
 }
 double GetDinuc23_GCTC(){
     return dinuc_stat[1][9][13];
 }
 double GetDinuc23_GCTG(){
     return dinuc_stat[1][9][14];
 }
 double GetDinuc23_GCTT(){
     return dinuc_stat[1][9][15];
 }
 double GetDinuc23_GGAA(){
     return dinuc_stat[1][10][0];
 }
 double GetDinuc23_GGAC(){
     return dinuc_stat[1][10][1];
 }
 double GetDinuc23_GGAG(){
     return dinuc_stat[1][10][2];
 }
 double GetDinuc23_GGAT(){
     return dinuc_stat[1][10][3];
 }
 double GetDinuc23_GGCA(){
     return dinuc_stat[1][10][4];
 }
 double GetDinuc23_GGCC(){
     return dinuc_stat[1][10][5];
 }
 double GetDinuc23_GGCG(){
     return dinuc_stat[1][10][6];
 }
 double GetDinuc23_GGCT(){
     return dinuc_stat[1][10][7];
 }
 double GetDinuc23_GGGA(){
     return dinuc_stat[1][10][8];
 }
 double GetDinuc23_GGGC(){
     return dinuc_stat[1][10][9];
 }
 double GetDinuc23_GGGG(){
     return dinuc_stat[1][10][10];
 }
 double GetDinuc23_GGGT(){
     return dinuc_stat[1][10][11];
 }
 double GetDinuc23_GGTA(){
     return dinuc_stat[1][10][12];
 }
 double GetDinuc23_GGTC(){
     return dinuc_stat[1][10][13];
 }
 double GetDinuc23_GGTG(){
     return dinuc_stat[1][10][14];
 }
 double GetDinuc23_GGTT(){
     return dinuc_stat[1][10][15];
 }
 double GetDinuc23_GTAA(){
     return dinuc_stat[1][11][0];
 }
 double GetDinuc23_GTAC(){
     return dinuc_stat[1][11][1];
 }
 double GetDinuc23_GTAG(){
     return dinuc_stat[1][11][2];
 }
 double GetDinuc23_GTAT(){
     return dinuc_stat[1][11][3];
 }
 double GetDinuc23_GTCA(){
     return dinuc_stat[1][11][4];
 }
 double GetDinuc23_GTCC(){
     return dinuc_stat[1][11][5];
 }
 double GetDinuc23_GTCG(){
     return dinuc_stat[1][11][6];
 }
 double GetDinuc23_GTCT(){
     return dinuc_stat[1][11][7];
 }
 double GetDinuc23_GTGA(){
     return dinuc_stat[1][11][8];
 }
 double GetDinuc23_GTGC(){
     return dinuc_stat[1][11][9];
 }
 double GetDinuc23_GTGG(){
     return dinuc_stat[1][11][10];
 }
 double GetDinuc23_GTGT(){
     return dinuc_stat[1][11][11];
 }
 double GetDinuc23_GTTA(){
     return dinuc_stat[1][11][12];
 }
 double GetDinuc23_GTTC(){
     return dinuc_stat[1][11][13];
 }
 double GetDinuc23_GTTG(){
     return dinuc_stat[1][11][14];
 }
 double GetDinuc23_GTTT(){
     return dinuc_stat[1][11][15];
 }
 double GetDinuc23_TAAA(){
     return dinuc_stat[1][12][0];
 }
 double GetDinuc23_TAAC(){
     return dinuc_stat[1][12][1];
 }
 double GetDinuc23_TAAG(){
     return dinuc_stat[1][12][2];
 }
 double GetDinuc23_TAAT(){
     return dinuc_stat[1][12][3];
 }
 double GetDinuc23_TACA(){
     return dinuc_stat[1][12][4];
 }
 double GetDinuc23_TACC(){
     return dinuc_stat[1][12][5];
 }
 double GetDinuc23_TACG(){
     return dinuc_stat[1][12][6];
 }
 double GetDinuc23_TACT(){
     return dinuc_stat[1][12][7];
 }
 double GetDinuc23_TAGA(){
     return dinuc_stat[1][12][8];
 }
 double GetDinuc23_TAGC(){
     return dinuc_stat[1][12][9];
 }
 double GetDinuc23_TAGG(){
     return dinuc_stat[1][12][10];
 }
 double GetDinuc23_TAGT(){
     return dinuc_stat[1][12][11];
 }
 double GetDinuc23_TATA(){
     return dinuc_stat[1][12][12];
 }
 double GetDinuc23_TATC(){
     return dinuc_stat[1][12][13];
 }
 double GetDinuc23_TATG(){
     return dinuc_stat[1][12][14];
 }
 double GetDinuc23_TATT(){
     return dinuc_stat[1][12][15];
 }
 double GetDinuc23_TCAA(){
     return dinuc_stat[1][13][0];
 }
 double GetDinuc23_TCAC(){
     return dinuc_stat[1][13][1];
 }
 double GetDinuc23_TCAG(){
     return dinuc_stat[1][13][2];
 }
 double GetDinuc23_TCAT(){
     return dinuc_stat[1][13][3];
 }
 double GetDinuc23_TCCA(){
     return dinuc_stat[1][13][4];
 }
 double GetDinuc23_TCCC(){
     return dinuc_stat[1][13][5];
 }
 double GetDinuc23_TCCG(){
     return dinuc_stat[1][13][6];
 }
 double GetDinuc23_TCCT(){
     return dinuc_stat[1][13][7];
 }
 double GetDinuc23_TCGA(){
     return dinuc_stat[1][13][8];
 }
 double GetDinuc23_TCGC(){
     return dinuc_stat[1][13][9];
 }
 double GetDinuc23_TCGG(){
     return dinuc_stat[1][13][10];
 }
 double GetDinuc23_TCGT(){
     return dinuc_stat[1][13][11];
 }
 double GetDinuc23_TCTA(){
     return dinuc_stat[1][13][12];
 }
 double GetDinuc23_TCTC(){
     return dinuc_stat[1][13][13];
 }
 double GetDinuc23_TCTG(){
     return dinuc_stat[1][13][14];
 }
 double GetDinuc23_TCTT(){
     return dinuc_stat[1][13][15];
 }
 double GetDinuc23_TGAA(){
     return dinuc_stat[1][14][0];
 }
 double GetDinuc23_TGAC(){
     return dinuc_stat[1][14][1];
 }
 double GetDinuc23_TGAG(){
     return dinuc_stat[1][14][2];
 }
 double GetDinuc23_TGAT(){
     return dinuc_stat[1][14][3];
 }
 double GetDinuc23_TGCA(){
     return dinuc_stat[1][14][4];
 }
 double GetDinuc23_TGCC(){
     return dinuc_stat[1][14][5];
 }
 double GetDinuc23_TGCG(){
     return dinuc_stat[1][14][6];
 }
 double GetDinuc23_TGCT(){
     return dinuc_stat[1][14][7];
 }
 double GetDinuc23_TGGA(){
     return dinuc_stat[1][14][8];
 }
 double GetDinuc23_TGGC(){
     return dinuc_stat[1][14][9];
 }
 double GetDinuc23_TGGG(){
     return dinuc_stat[1][14][10];
 }
 double GetDinuc23_TGGT(){
     return dinuc_stat[1][14][11];
 }
 double GetDinuc23_TGTA(){
     return dinuc_stat[1][14][12];
 }
 double GetDinuc23_TGTC(){
     return dinuc_stat[1][14][13];
 }
 double GetDinuc23_TGTG(){
     return dinuc_stat[1][14][14];
 }
 double GetDinuc23_TGTT(){
     return dinuc_stat[1][14][15];
 }
 double GetDinuc23_TTAA(){
     return dinuc_stat[1][15][0];
 }
 double GetDinuc23_TTAC(){
     return dinuc_stat[1][15][1];
 }
 double GetDinuc23_TTAG(){
     return dinuc_stat[1][15][2];
 }
 double GetDinuc23_TTAT(){
     return dinuc_stat[1][15][3];
 }
 double GetDinuc23_TTCA(){
     return dinuc_stat[1][15][4];
 }
 double GetDinuc23_TTCC(){
     return dinuc_stat[1][15][5];
 }
 double GetDinuc23_TTCG(){
     return dinuc_stat[1][15][6];
 }
 double GetDinuc23_TTCT(){
     return dinuc_stat[1][15][7];
 }
 double GetDinuc23_TTGA(){
     return dinuc_stat[1][15][8];
 }
 double GetDinuc23_TTGC(){
     return dinuc_stat[1][15][9];
 }
 double GetDinuc23_TTGG(){
     return dinuc_stat[1][15][10];
 }
 double GetDinuc23_TTGT(){
     return dinuc_stat[1][15][11];
 }
 double GetDinuc23_TTTA(){
     return dinuc_stat[1][15][12];
 }
 double GetDinuc23_TTTC(){
     return dinuc_stat[1][15][13];
 }
 double GetDinuc23_TTTG(){
     return dinuc_stat[1][15][14];
 }
 double GetDinuc23_TTTT(){
     return dinuc_stat[1][15][15];
 }

 double GetDinuc31_AAAA(){
     return dinuc_stat[2][0][0];
  }
  double GetDinuc31_AAAC(){
     return dinuc_stat[2][0][1];
  }
  double GetDinuc31_AAAG(){
     return dinuc_stat[2][0][2];
  }
  double GetDinuc31_AAAT(){
     return dinuc_stat[2][0][3];
  }
  double GetDinuc31_AACA(){
     return dinuc_stat[2][0][4];
  }
  double GetDinuc31_AACC(){
     return dinuc_stat[2][0][5];
  }
  double GetDinuc31_AACG(){
     return dinuc_stat[2][0][6];
  }
  double GetDinuc31_AACT(){
     return dinuc_stat[2][0][7];
  }
  double GetDinuc31_AAGA(){
     return dinuc_stat[2][0][8];
  }
  double GetDinuc31_AAGC(){
     return dinuc_stat[2][0][9];
  }
  double GetDinuc31_AAGG(){
     return dinuc_stat[2][0][10];
  }
  double GetDinuc31_AAGT(){
     return dinuc_stat[2][0][11];
  }
  double GetDinuc31_AATA(){
     return dinuc_stat[2][0][12];
  }
  double GetDinuc31_AATC(){
     return dinuc_stat[2][0][13];
  }
  double GetDinuc31_AATG(){
     return dinuc_stat[2][0][14];
  }
  double GetDinuc31_AATT(){
     return dinuc_stat[2][0][15];
  }
  double GetDinuc31_ACAA(){
     return dinuc_stat[2][1][0];
  }
  double GetDinuc31_ACAC(){
     return dinuc_stat[2][1][1];
  }
  double GetDinuc31_ACAG(){
     return dinuc_stat[2][1][2];
  }
  double GetDinuc31_ACAT(){
     return dinuc_stat[2][1][3];
  }
  double GetDinuc31_ACCA(){
     return dinuc_stat[2][1][4];
  }
  double GetDinuc31_ACCC(){
     return dinuc_stat[2][1][5];
  }
  double GetDinuc31_ACCG(){
     return dinuc_stat[2][1][6];
  }
  double GetDinuc31_ACCT(){
     return dinuc_stat[2][1][7];
  }
  double GetDinuc31_ACGA(){
     return dinuc_stat[2][1][8];
  }
  double GetDinuc31_ACGC(){
     return dinuc_stat[2][1][9];
  }
  double GetDinuc31_ACGG(){
     return dinuc_stat[2][1][10];
  }
  double GetDinuc31_ACGT(){
     return dinuc_stat[2][1][11];
  }
  double GetDinuc31_ACTA(){
     return dinuc_stat[2][1][12];
  }
  double GetDinuc31_ACTC(){
     return dinuc_stat[2][1][13];
  }
  double GetDinuc31_ACTG(){
     return dinuc_stat[2][1][14];
  }
  double GetDinuc31_ACTT(){
     return dinuc_stat[2][1][15];
  }
  double GetDinuc31_AGAA(){
     return dinuc_stat[2][2][0];
  }
  double GetDinuc31_AGAC(){
     return dinuc_stat[2][2][1];
  }
  double GetDinuc31_AGAG(){
     return dinuc_stat[2][2][2];
  }
  double GetDinuc31_AGAT(){
     return dinuc_stat[2][2][3];
  }
  double GetDinuc31_AGCA(){
     return dinuc_stat[2][2][4];
  }
  double GetDinuc31_AGCC(){
     return dinuc_stat[2][2][5];
  }
  double GetDinuc31_AGCG(){
     return dinuc_stat[2][2][6];
  }
  double GetDinuc31_AGCT(){
     return dinuc_stat[2][2][7];
  }
  double GetDinuc31_AGGA(){
     return dinuc_stat[2][2][8];
  }
  double GetDinuc31_AGGC(){
     return dinuc_stat[2][2][9];
  }
  double GetDinuc31_AGGG(){
     return dinuc_stat[2][2][10];
  }
  double GetDinuc31_AGGT(){
     return dinuc_stat[2][2][11];
  }
  double GetDinuc31_AGTA(){
     return dinuc_stat[2][2][12];
  }
  double GetDinuc31_AGTC(){
     return dinuc_stat[2][2][13];
  }
  double GetDinuc31_AGTG(){
     return dinuc_stat[2][2][14];
  }
  double GetDinuc31_AGTT(){
     return dinuc_stat[2][2][15];
  }
  double GetDinuc31_ATAA(){
     return dinuc_stat[2][3][0];
  }
  double GetDinuc31_ATAC(){
     return dinuc_stat[2][3][1];
  }
  double GetDinuc31_ATAG(){
     return dinuc_stat[2][3][2];
  }
  double GetDinuc31_ATAT(){
     return dinuc_stat[2][3][3];
  }
  double GetDinuc31_ATCA(){
     return dinuc_stat[2][3][4];
  }
  double GetDinuc31_ATCC(){
     return dinuc_stat[2][3][5];
  }
  double GetDinuc31_ATCG(){
     return dinuc_stat[2][3][6];
  }
  double GetDinuc31_ATCT(){
     return dinuc_stat[2][3][7];
  }
  double GetDinuc31_ATGA(){
     return dinuc_stat[2][3][8];
  }
  double GetDinuc31_ATGC(){
     return dinuc_stat[2][3][9];
  }
  double GetDinuc31_ATGG(){
     return dinuc_stat[2][3][10];
  }
  double GetDinuc31_ATGT(){
     return dinuc_stat[2][3][11];
  }
  double GetDinuc31_ATTA(){
     return dinuc_stat[2][3][12];
  }
  double GetDinuc31_ATTC(){
     return dinuc_stat[2][3][13];
  }
  double GetDinuc31_ATTG(){
     return dinuc_stat[2][3][14];
  }
  double GetDinuc31_ATTT(){
     return dinuc_stat[2][3][15];
  }
  double GetDinuc31_CAAA(){
     return dinuc_stat[2][4][0];
  }
  double GetDinuc31_CAAC(){
     return dinuc_stat[2][4][1];
  }
  double GetDinuc31_CAAG(){
     return dinuc_stat[2][4][2];
  }
  double GetDinuc31_CAAT(){
     return dinuc_stat[2][4][3];
  }
  double GetDinuc31_CACA(){
     return dinuc_stat[2][4][4];
  }
  double GetDinuc31_CACC(){
     return dinuc_stat[2][4][5];
  }
  double GetDinuc31_CACG(){
     return dinuc_stat[2][4][6];
  }
  double GetDinuc31_CACT(){
     return dinuc_stat[2][4][7];
  }
  double GetDinuc31_CAGA(){
     return dinuc_stat[2][4][8];
  }
  double GetDinuc31_CAGC(){
     return dinuc_stat[2][4][9];
  }
  double GetDinuc31_CAGG(){
     return dinuc_stat[2][4][10];
  }
  double GetDinuc31_CAGT(){
     return dinuc_stat[2][4][11];
  }
  double GetDinuc31_CATA(){
     return dinuc_stat[2][4][12];
  }
  double GetDinuc31_CATC(){
     return dinuc_stat[2][4][13];
  }
  double GetDinuc31_CATG(){
     return dinuc_stat[2][4][14];
  }
  double GetDinuc31_CATT(){
     return dinuc_stat[2][4][15];
  }
  double GetDinuc31_CCAA(){
     return dinuc_stat[2][5][0];
  }
  double GetDinuc31_CCAC(){
     return dinuc_stat[2][5][1];
  }
  double GetDinuc31_CCAG(){
     return dinuc_stat[2][5][2];
  }
  double GetDinuc31_CCAT(){
     return dinuc_stat[2][5][3];
  }
  double GetDinuc31_CCCA(){
     return dinuc_stat[2][5][4];
  }
  double GetDinuc31_CCCC(){
     return dinuc_stat[2][5][5];
  }
  double GetDinuc31_CCCG(){
     return dinuc_stat[2][5][6];
  }
  double GetDinuc31_CCCT(){
     return dinuc_stat[2][5][7];
  }
  double GetDinuc31_CCGA(){
     return dinuc_stat[2][5][8];
  }
  double GetDinuc31_CCGC(){
     return dinuc_stat[2][5][9];
  }
  double GetDinuc31_CCGG(){
     return dinuc_stat[2][5][10];
  }
  double GetDinuc31_CCGT(){
     return dinuc_stat[2][5][11];
  }
  double GetDinuc31_CCTA(){
     return dinuc_stat[2][5][12];
  }
  double GetDinuc31_CCTC(){
     return dinuc_stat[2][5][13];
  }
  double GetDinuc31_CCTG(){
     return dinuc_stat[2][5][14];
  }
  double GetDinuc31_CCTT(){
     return dinuc_stat[2][5][15];
  }
  double GetDinuc31_CGAA(){
     return dinuc_stat[2][6][0];
  }
  double GetDinuc31_CGAC(){
     return dinuc_stat[2][6][1];
  }
  double GetDinuc31_CGAG(){
     return dinuc_stat[2][6][2];
  }
  double GetDinuc31_CGAT(){
     return dinuc_stat[2][6][3];
  }
  double GetDinuc31_CGCA(){
     return dinuc_stat[2][6][4];
  }
  double GetDinuc31_CGCC(){
     return dinuc_stat[2][6][5];
  }
  double GetDinuc31_CGCG(){
     return dinuc_stat[2][6][6];
  }
  double GetDinuc31_CGCT(){
     return dinuc_stat[2][6][7];
  }
  double GetDinuc31_CGGA(){
     return dinuc_stat[2][6][8];
  }
  double GetDinuc31_CGGC(){
     return dinuc_stat[2][6][9];
  }
  double GetDinuc31_CGGG(){
     return dinuc_stat[2][6][10];
  }
  double GetDinuc31_CGGT(){
     return dinuc_stat[2][6][11];
  }
  double GetDinuc31_CGTA(){
     return dinuc_stat[2][6][12];
  }
  double GetDinuc31_CGTC(){
     return dinuc_stat[2][6][13];
  }
  double GetDinuc31_CGTG(){
     return dinuc_stat[2][6][14];
  }
  double GetDinuc31_CGTT(){
     return dinuc_stat[2][6][15];
  }
  double GetDinuc31_CTAA(){
     return dinuc_stat[2][7][0];
  }
  double GetDinuc31_CTAC(){
     return dinuc_stat[2][7][1];
  }
  double GetDinuc31_CTAG(){
     return dinuc_stat[2][7][2];
  }
  double GetDinuc31_CTAT(){
     return dinuc_stat[2][7][3];
  }
  double GetDinuc31_CTCA(){
     return dinuc_stat[2][7][4];
  }
  double GetDinuc31_CTCC(){
     return dinuc_stat[2][7][5];
  }
  double GetDinuc31_CTCG(){
     return dinuc_stat[2][7][6];
  }
  double GetDinuc31_CTCT(){
     return dinuc_stat[2][7][7];
  }
  double GetDinuc31_CTGA(){
     return dinuc_stat[2][7][8];
  }
  double GetDinuc31_CTGC(){
     return dinuc_stat[2][7][9];
  }
  double GetDinuc31_CTGG(){
     return dinuc_stat[2][7][10];
  }
  double GetDinuc31_CTGT(){
     return dinuc_stat[2][7][11];
  }
  double GetDinuc31_CTTA(){
     return dinuc_stat[2][7][12];
  }
  double GetDinuc31_CTTC(){
     return dinuc_stat[2][7][13];
  }
  double GetDinuc31_CTTG(){
     return dinuc_stat[2][7][14];
  }
  double GetDinuc31_CTTT(){
     return dinuc_stat[2][7][15];
  }
  double GetDinuc31_GAAA(){
     return dinuc_stat[2][8][0];
  }
  double GetDinuc31_GAAC(){
     return dinuc_stat[2][8][1];
  }
  double GetDinuc31_GAAG(){
     return dinuc_stat[2][8][2];
  }
  double GetDinuc31_GAAT(){
     return dinuc_stat[2][8][3];
  }
  double GetDinuc31_GACA(){
     return dinuc_stat[2][8][4];
  }
  double GetDinuc31_GACC(){
     return dinuc_stat[2][8][5];
  }
  double GetDinuc31_GACG(){
     return dinuc_stat[2][8][6];
  }
  double GetDinuc31_GACT(){
     return dinuc_stat[2][8][7];
  }
  double GetDinuc31_GAGA(){
     return dinuc_stat[2][8][8];
  }
  double GetDinuc31_GAGC(){
     return dinuc_stat[2][8][9];
  }
  double GetDinuc31_GAGG(){
     return dinuc_stat[2][8][10];
  }
  double GetDinuc31_GAGT(){
     return dinuc_stat[2][8][11];
  }
  double GetDinuc31_GATA(){
     return dinuc_stat[2][8][12];
  }
  double GetDinuc31_GATC(){
     return dinuc_stat[2][8][13];
  }
  double GetDinuc31_GATG(){
     return dinuc_stat[2][8][14];
  }
  double GetDinuc31_GATT(){
     return dinuc_stat[2][8][15];
  }
  double GetDinuc31_GCAA(){
     return dinuc_stat[2][9][0];
  }
  double GetDinuc31_GCAC(){
     return dinuc_stat[2][9][1];
  }
  double GetDinuc31_GCAG(){
     return dinuc_stat[2][9][2];
  }
  double GetDinuc31_GCAT(){
     return dinuc_stat[2][9][3];
  }
  double GetDinuc31_GCCA(){
     return dinuc_stat[2][9][4];
  }
  double GetDinuc31_GCCC(){
     return dinuc_stat[2][9][5];
  }
  double GetDinuc31_GCCG(){
     return dinuc_stat[2][9][6];
  }
  double GetDinuc31_GCCT(){
     return dinuc_stat[2][9][7];
  }
  double GetDinuc31_GCGA(){
     return dinuc_stat[2][9][8];
  }
  double GetDinuc31_GCGC(){
     return dinuc_stat[2][9][9];
  }
  double GetDinuc31_GCGG(){
     return dinuc_stat[2][9][10];
  }
  double GetDinuc31_GCGT(){
     return dinuc_stat[2][9][11];
  }
  double GetDinuc31_GCTA(){
     return dinuc_stat[2][9][12];
  }
  double GetDinuc31_GCTC(){
     return dinuc_stat[2][9][13];
  }
  double GetDinuc31_GCTG(){
     return dinuc_stat[2][9][14];
  }
  double GetDinuc31_GCTT(){
     return dinuc_stat[2][9][15];
  }
  double GetDinuc31_GGAA(){
     return dinuc_stat[2][10][0];
  }
  double GetDinuc31_GGAC(){
     return dinuc_stat[2][10][1];
  }
  double GetDinuc31_GGAG(){
     return dinuc_stat[2][10][2];
  }
  double GetDinuc31_GGAT(){
     return dinuc_stat[2][10][3];
  }
  double GetDinuc31_GGCA(){
     return dinuc_stat[2][10][4];
  }
  double GetDinuc31_GGCC(){
     return dinuc_stat[2][10][5];
  }
  double GetDinuc31_GGCG(){
     return dinuc_stat[2][10][6];
  }
  double GetDinuc31_GGCT(){
     return dinuc_stat[2][10][7];
  }
  double GetDinuc31_GGGA(){
     return dinuc_stat[2][10][8];
  }
  double GetDinuc31_GGGC(){
     return dinuc_stat[2][10][9];
  }
  double GetDinuc31_GGGG(){
     return dinuc_stat[2][10][10];
  }
  double GetDinuc31_GGGT(){
     return dinuc_stat[2][10][11];
  }
  double GetDinuc31_GGTA(){
     return dinuc_stat[2][10][12];
  }
  double GetDinuc31_GGTC(){
     return dinuc_stat[2][10][13];
  }
  double GetDinuc31_GGTG(){
     return dinuc_stat[2][10][14];
  }
  double GetDinuc31_GGTT(){
     return dinuc_stat[2][10][15];
  }
  double GetDinuc31_GTAA(){
     return dinuc_stat[2][11][0];
  }
  double GetDinuc31_GTAC(){
     return dinuc_stat[2][11][1];
  }
  double GetDinuc31_GTAG(){
     return dinuc_stat[2][11][2];
  }
  double GetDinuc31_GTAT(){
     return dinuc_stat[2][11][3];
  }
  double GetDinuc31_GTCA(){
     return dinuc_stat[2][11][4];
  }
  double GetDinuc31_GTCC(){
     return dinuc_stat[2][11][5];
  }
  double GetDinuc31_GTCG(){
     return dinuc_stat[2][11][6];
  }
  double GetDinuc31_GTCT(){
     return dinuc_stat[2][11][7];
  }
  double GetDinuc31_GTGA(){
     return dinuc_stat[2][11][8];
  }
  double GetDinuc31_GTGC(){
     return dinuc_stat[2][11][9];
  }
  double GetDinuc31_GTGG(){
     return dinuc_stat[2][11][10];
  }
  double GetDinuc31_GTGT(){
     return dinuc_stat[2][11][11];
  }
  double GetDinuc31_GTTA(){
     return dinuc_stat[2][11][12];
  }
  double GetDinuc31_GTTC(){
     return dinuc_stat[2][11][13];
  }
  double GetDinuc31_GTTG(){
     return dinuc_stat[2][11][14];
  }
  double GetDinuc31_GTTT(){
     return dinuc_stat[2][11][15];
  }
  double GetDinuc31_TAAA(){
     return dinuc_stat[2][12][0];
  }
  double GetDinuc31_TAAC(){
     return dinuc_stat[2][12][1];
  }
  double GetDinuc31_TAAG(){
     return dinuc_stat[2][12][2];
  }
  double GetDinuc31_TAAT(){
     return dinuc_stat[2][12][3];
  }
  double GetDinuc31_TACA(){
     return dinuc_stat[2][12][4];
  }
  double GetDinuc31_TACC(){
     return dinuc_stat[2][12][5];
  }
  double GetDinuc31_TACG(){
     return dinuc_stat[2][12][6];
  }
  double GetDinuc31_TACT(){
     return dinuc_stat[2][12][7];
  }
  double GetDinuc31_TAGA(){
     return dinuc_stat[2][12][8];
  }
  double GetDinuc31_TAGC(){
     return dinuc_stat[2][12][9];
  }
  double GetDinuc31_TAGG(){
     return dinuc_stat[2][12][10];
  }
  double GetDinuc31_TAGT(){
     return dinuc_stat[2][12][11];
  }
  double GetDinuc31_TATA(){
     return dinuc_stat[2][12][12];
  }
  double GetDinuc31_TATC(){
     return dinuc_stat[2][12][13];
  }
  double GetDinuc31_TATG(){
     return dinuc_stat[2][12][14];
  }
  double GetDinuc31_TATT(){
     return dinuc_stat[2][12][15];
  }
  double GetDinuc31_TCAA(){
     return dinuc_stat[2][13][0];
  }
  double GetDinuc31_TCAC(){
     return dinuc_stat[2][13][1];
  }
  double GetDinuc31_TCAG(){
     return dinuc_stat[2][13][2];
  }
  double GetDinuc31_TCAT(){
     return dinuc_stat[2][13][3];
  }
  double GetDinuc31_TCCA(){
     return dinuc_stat[2][13][4];
  }
  double GetDinuc31_TCCC(){
     return dinuc_stat[2][13][5];
  }
  double GetDinuc31_TCCG(){
     return dinuc_stat[2][13][6];
  }
  double GetDinuc31_TCCT(){
     return dinuc_stat[2][13][7];
  }
  double GetDinuc31_TCGA(){
     return dinuc_stat[2][13][8];
  }
  double GetDinuc31_TCGC(){
     return dinuc_stat[2][13][9];
  }
  double GetDinuc31_TCGG(){
     return dinuc_stat[2][13][10];
  }
  double GetDinuc31_TCGT(){
     return dinuc_stat[2][13][11];
  }
  double GetDinuc31_TCTA(){
     return dinuc_stat[2][13][12];
  }
  double GetDinuc31_TCTC(){
     return dinuc_stat[2][13][13];
  }
  double GetDinuc31_TCTG(){
     return dinuc_stat[2][13][14];
  }
  double GetDinuc31_TCTT(){
     return dinuc_stat[2][13][15];
  }
  double GetDinuc31_TGAA(){
     return dinuc_stat[2][14][0];
  }
  double GetDinuc31_TGAC(){
     return dinuc_stat[2][14][1];
  }
  double GetDinuc31_TGAG(){
     return dinuc_stat[2][14][2];
  }
  double GetDinuc31_TGAT(){
     return dinuc_stat[2][14][3];
  }
  double GetDinuc31_TGCA(){
     return dinuc_stat[2][14][4];
  }
  double GetDinuc31_TGCC(){
     return dinuc_stat[2][14][5];
  }
  double GetDinuc31_TGCG(){
     return dinuc_stat[2][14][6];
  }
  double GetDinuc31_TGCT(){
     return dinuc_stat[2][14][7];
  }
  double GetDinuc31_TGGA(){
     return dinuc_stat[2][14][8];
  }
  double GetDinuc31_TGGC(){
     return dinuc_stat[2][14][9];
  }
  double GetDinuc31_TGGG(){
     return dinuc_stat[2][14][10];
  }
  double GetDinuc31_TGGT(){
     return dinuc_stat[2][14][11];
  }
  double GetDinuc31_TGTA(){
     return dinuc_stat[2][14][12];
  }
  double GetDinuc31_TGTC(){
     return dinuc_stat[2][14][13];
  }
  double GetDinuc31_TGTG(){
     return dinuc_stat[2][14][14];
  }
  double GetDinuc31_TGTT(){
     return dinuc_stat[2][14][15];
  }
  double GetDinuc31_TTAA(){
     return dinuc_stat[2][15][0];
  }
  double GetDinuc31_TTAC(){
     return dinuc_stat[2][15][1];
  }
  double GetDinuc31_TTAG(){
     return dinuc_stat[2][15][2];
  }
  double GetDinuc31_TTAT(){
     return dinuc_stat[2][15][3];
  }
  double GetDinuc31_TTCA(){
     return dinuc_stat[2][15][4];
  }
  double GetDinuc31_TTCC(){
     return dinuc_stat[2][15][5];
  }
  double GetDinuc31_TTCG(){
     return dinuc_stat[2][15][6];
  }
  double GetDinuc31_TTCT(){
     return dinuc_stat[2][15][7];
  }
  double GetDinuc31_TTGA(){
     return dinuc_stat[2][15][8];
  }
  double GetDinuc31_TTGC(){
     return dinuc_stat[2][15][9];
  }
  double GetDinuc31_TTGG(){
     return dinuc_stat[2][15][10];
  }
  double GetDinuc31_TTGT(){
     return dinuc_stat[2][15][11];
  }
  double GetDinuc31_TTTA(){
     return dinuc_stat[2][15][12];
  }
  double GetDinuc31_TTTC(){
     return dinuc_stat[2][15][13];
  }
  double GetDinuc31_TTTG(){
     return dinuc_stat[2][15][14];
  }
  double GetDinuc31_TTTT(){
     return dinuc_stat[2][15][15];
  }

double GetDinucSyn_AAAA(){
     return dinucSyn_stat[3][0][0];
 }
 double GetDinucSyn_AAAC(){
     return dinucSyn_stat[3][0][1];
 }
 double GetDinucSyn_AAAG(){
     return dinucSyn_stat[3][0][2];
 }
 double GetDinucSyn_AAAT(){
     return dinucSyn_stat[3][0][3];
 }
 double GetDinucSyn_AACA(){
     return dinucSyn_stat[3][0][4];
 }
 double GetDinucSyn_AACC(){
     return dinucSyn_stat[3][0][5];
 }
 double GetDinucSyn_AACG(){
     return dinucSyn_stat[3][0][6];
 }
 double GetDinucSyn_AACT(){
     return dinucSyn_stat[3][0][7];
 }
 double GetDinucSyn_AAGA(){
     return dinucSyn_stat[3][0][8];
 }
 double GetDinucSyn_AAGC(){
     return dinucSyn_stat[3][0][9];
 }
 double GetDinucSyn_AAGG(){
     return dinucSyn_stat[3][0][10];
 }
 double GetDinucSyn_AAGT(){
     return dinucSyn_stat[3][0][11];
 }
 double GetDinucSyn_AATA(){
     return dinucSyn_stat[3][0][12];
 }
 double GetDinucSyn_AATC(){
     return dinucSyn_stat[3][0][13];
 }
 double GetDinucSyn_AATG(){
     return dinucSyn_stat[3][0][14];
 }
 double GetDinucSyn_AATT(){
     return dinucSyn_stat[3][0][15];
 }
 double GetDinucSyn_ACAA(){
     return dinucSyn_stat[3][1][0];
 }
 double GetDinucSyn_ACAC(){
     return dinucSyn_stat[3][1][1];
 }
 double GetDinucSyn_ACAG(){
     return dinucSyn_stat[3][1][2];
 }
 double GetDinucSyn_ACAT(){
     return dinucSyn_stat[3][1][3];
 }
 double GetDinucSyn_ACCA(){
     return dinucSyn_stat[3][1][4];
 }
 double GetDinucSyn_ACCC(){
     return dinucSyn_stat[3][1][5];
 }
 double GetDinucSyn_ACCG(){
     return dinucSyn_stat[3][1][6];
 }
 double GetDinucSyn_ACCT(){
     return dinucSyn_stat[3][1][7];
 }
 double GetDinucSyn_ACGA(){
     return dinucSyn_stat[3][1][8];
 }
 double GetDinucSyn_ACGC(){
     return dinucSyn_stat[3][1][9];
 }
 double GetDinucSyn_ACGG(){
     return dinucSyn_stat[3][1][10];
 }
 double GetDinucSyn_ACGT(){
     return dinucSyn_stat[3][1][11];
 }
 double GetDinucSyn_ACTA(){
     return dinucSyn_stat[3][1][12];
 }
 double GetDinucSyn_ACTC(){
     return dinucSyn_stat[3][1][13];
 }
 double GetDinucSyn_ACTG(){
     return dinucSyn_stat[3][1][14];
 }
 double GetDinucSyn_ACTT(){
     return dinucSyn_stat[3][1][15];
 }
 double GetDinucSyn_AGAA(){
     return dinucSyn_stat[3][2][0];
 }
 double GetDinucSyn_AGAC(){
     return dinucSyn_stat[3][2][1];
 }
 double GetDinucSyn_AGAG(){
     return dinucSyn_stat[3][2][2];
 }
 double GetDinucSyn_AGAT(){
     return dinucSyn_stat[3][2][3];
 }
 double GetDinucSyn_AGCA(){
     return dinucSyn_stat[3][2][4];
 }
 double GetDinucSyn_AGCC(){
     return dinucSyn_stat[3][2][5];
 }
 double GetDinucSyn_AGCG(){
     return dinucSyn_stat[3][2][6];
 }
 double GetDinucSyn_AGCT(){
     return dinucSyn_stat[3][2][7];
 }
 double GetDinucSyn_AGGA(){
     return dinucSyn_stat[3][2][8];
 }
 double GetDinucSyn_AGGC(){
     return dinucSyn_stat[3][2][9];
 }
 double GetDinucSyn_AGGG(){
     return dinucSyn_stat[3][2][10];
 }
 double GetDinucSyn_AGGT(){
     return dinucSyn_stat[3][2][11];
 }
 double GetDinucSyn_AGTA(){
     return dinucSyn_stat[3][2][12];
 }
 double GetDinucSyn_AGTC(){
     return dinucSyn_stat[3][2][13];
 }
 double GetDinucSyn_AGTG(){
     return dinucSyn_stat[3][2][14];
 }
 double GetDinucSyn_AGTT(){
     return dinucSyn_stat[3][2][15];
 }
 double GetDinucSyn_ATAA(){
     return dinucSyn_stat[3][3][0];
 }
 double GetDinucSyn_ATAC(){
     return dinucSyn_stat[3][3][1];
 }
 double GetDinucSyn_ATAG(){
     return dinucSyn_stat[3][3][2];
 }
 double GetDinucSyn_ATAT(){
     return dinucSyn_stat[3][3][3];
 }
 double GetDinucSyn_ATCA(){
     return dinucSyn_stat[3][3][4];
 }
 double GetDinucSyn_ATCC(){
     return dinucSyn_stat[3][3][5];
 }
 double GetDinucSyn_ATCG(){
     return dinucSyn_stat[3][3][6];
 }
 double GetDinucSyn_ATCT(){
     return dinucSyn_stat[3][3][7];
 }
 double GetDinucSyn_ATGA(){
     return dinucSyn_stat[3][3][8];
 }
 double GetDinucSyn_ATGC(){
     return dinucSyn_stat[3][3][9];
 }
 double GetDinucSyn_ATGG(){
     return dinucSyn_stat[3][3][10];
 }
 double GetDinucSyn_ATGT(){
     return dinucSyn_stat[3][3][11];
 }
 double GetDinucSyn_ATTA(){
     return dinucSyn_stat[3][3][12];
 }
 double GetDinucSyn_ATTC(){
     return dinucSyn_stat[3][3][13];
 }
 double GetDinucSyn_ATTG(){
     return dinucSyn_stat[3][3][14];
 }
 double GetDinucSyn_ATTT(){
     return dinucSyn_stat[3][3][15];
 }
 double GetDinucSyn_CAAA(){
     return dinucSyn_stat[3][4][0];
 }
 double GetDinucSyn_CAAC(){
     return dinucSyn_stat[3][4][1];
 }
 double GetDinucSyn_CAAG(){
     return dinucSyn_stat[3][4][2];
 }
 double GetDinucSyn_CAAT(){
     return dinucSyn_stat[3][4][3];
 }
 double GetDinucSyn_CACA(){
     return dinucSyn_stat[3][4][4];
 }
 double GetDinucSyn_CACC(){
     return dinucSyn_stat[3][4][5];
 }
 double GetDinucSyn_CACG(){
     return dinucSyn_stat[3][4][6];
 }
 double GetDinucSyn_CACT(){
     return dinucSyn_stat[3][4][7];
 }
 double GetDinucSyn_CAGA(){
     return dinucSyn_stat[3][4][8];
 }
 double GetDinucSyn_CAGC(){
     return dinucSyn_stat[3][4][9];
 }
 double GetDinucSyn_CAGG(){
     return dinucSyn_stat[3][4][10];
 }
 double GetDinucSyn_CAGT(){
     return dinucSyn_stat[3][4][11];
 }
 double GetDinucSyn_CATA(){
     return dinucSyn_stat[3][4][12];
 }
 double GetDinucSyn_CATC(){
     return dinucSyn_stat[3][4][13];
 }
 double GetDinucSyn_CATG(){
     return dinucSyn_stat[3][4][14];
 }
 double GetDinucSyn_CATT(){
     return dinucSyn_stat[3][4][15];
 }
 double GetDinucSyn_CCAA(){
     return dinucSyn_stat[3][5][0];
 }
 double GetDinucSyn_CCAC(){
     return dinucSyn_stat[3][5][1];
 }
 double GetDinucSyn_CCAG(){
     return dinucSyn_stat[3][5][2];
 }
 double GetDinucSyn_CCAT(){
     return dinucSyn_stat[3][5][3];
 }
 double GetDinucSyn_CCCA(){
     return dinucSyn_stat[3][5][4];
 }
 double GetDinucSyn_CCCC(){
     return dinucSyn_stat[3][5][5];
 }
 double GetDinucSyn_CCCG(){
     return dinucSyn_stat[3][5][6];
 }
 double GetDinucSyn_CCCT(){
     return dinucSyn_stat[3][5][7];
 }
 double GetDinucSyn_CCGA(){
     return dinucSyn_stat[3][5][8];
 }
 double GetDinucSyn_CCGC(){
     return dinucSyn_stat[3][5][9];
 }
 double GetDinucSyn_CCGG(){
     return dinucSyn_stat[3][5][10];
 }
 double GetDinucSyn_CCGT(){
     return dinucSyn_stat[3][5][11];
 }
 double GetDinucSyn_CCTA(){
     return dinucSyn_stat[3][5][12];
 }
 double GetDinucSyn_CCTC(){
     return dinucSyn_stat[3][5][13];
 }
 double GetDinucSyn_CCTG(){
     return dinucSyn_stat[3][5][14];
 }
 double GetDinucSyn_CCTT(){
     return dinucSyn_stat[3][5][15];
 }
 double GetDinucSyn_CGAA(){
     return dinucSyn_stat[3][6][0];
 }
 double GetDinucSyn_CGAC(){
     return dinucSyn_stat[3][6][1];
 }
 double GetDinucSyn_CGAG(){
     return dinucSyn_stat[3][6][2];
 }
 double GetDinucSyn_CGAT(){
     return dinucSyn_stat[3][6][3];
 }
 double GetDinucSyn_CGCA(){
     return dinucSyn_stat[3][6][4];
 }
 double GetDinucSyn_CGCC(){
     return dinucSyn_stat[3][6][5];
 }
 double GetDinucSyn_CGCG(){
     return dinucSyn_stat[3][6][6];
 }
 double GetDinucSyn_CGCT(){
     return dinucSyn_stat[3][6][7];
 }
 double GetDinucSyn_CGGA(){
     return dinucSyn_stat[3][6][8];
 }
 double GetDinucSyn_CGGC(){
     return dinucSyn_stat[3][6][9];
 }
 double GetDinucSyn_CGGG(){
     return dinucSyn_stat[3][6][10];
 }
 double GetDinucSyn_CGGT(){
     return dinucSyn_stat[3][6][11];
 }
 double GetDinucSyn_CGTA(){
     return dinucSyn_stat[3][6][12];
 }
 double GetDinucSyn_CGTC(){
     return dinucSyn_stat[3][6][13];
 }
 double GetDinucSyn_CGTG(){
     return dinucSyn_stat[3][6][14];
 }
 double GetDinucSyn_CGTT(){
     return dinucSyn_stat[3][6][15];
 }
 double GetDinucSyn_CTAA(){
     return dinucSyn_stat[3][7][0];
 }
 double GetDinucSyn_CTAC(){
     return dinucSyn_stat[3][7][1];
 }
 double GetDinucSyn_CTAG(){
     return dinucSyn_stat[3][7][2];
 }
 double GetDinucSyn_CTAT(){
     return dinucSyn_stat[3][7][3];
 }
 double GetDinucSyn_CTCA(){
     return dinucSyn_stat[3][7][4];
 }
 double GetDinucSyn_CTCC(){
     return dinucSyn_stat[3][7][5];
 }
 double GetDinucSyn_CTCG(){
     return dinucSyn_stat[3][7][6];
 }
 double GetDinucSyn_CTCT(){
     return dinucSyn_stat[3][7][7];
 }
 double GetDinucSyn_CTGA(){
     return dinucSyn_stat[3][7][8];
 }
 double GetDinucSyn_CTGC(){
     return dinucSyn_stat[3][7][9];
 }
 double GetDinucSyn_CTGG(){
     return dinucSyn_stat[3][7][10];
 }
 double GetDinucSyn_CTGT(){
     return dinucSyn_stat[3][7][11];
 }
 double GetDinucSyn_CTTA(){
     return dinucSyn_stat[3][7][12];
 }
 double GetDinucSyn_CTTC(){
     return dinucSyn_stat[3][7][13];
 }
 double GetDinucSyn_CTTG(){
     return dinucSyn_stat[3][7][14];
 }
 double GetDinucSyn_CTTT(){
     return dinucSyn_stat[3][7][15];
 }
 double GetDinucSyn_GAAA(){
     return dinucSyn_stat[3][8][0];
 }
 double GetDinucSyn_GAAC(){
     return dinucSyn_stat[3][8][1];
 }
 double GetDinucSyn_GAAG(){
     return dinucSyn_stat[3][8][2];
 }
 double GetDinucSyn_GAAT(){
     return dinucSyn_stat[3][8][3];
 }
 double GetDinucSyn_GACA(){
     return dinucSyn_stat[3][8][4];
 }
 double GetDinucSyn_GACC(){
     return dinucSyn_stat[3][8][5];
 }
 double GetDinucSyn_GACG(){
     return dinucSyn_stat[3][8][6];
 }
 double GetDinucSyn_GACT(){
     return dinucSyn_stat[3][8][7];
 }
 double GetDinucSyn_GAGA(){
     return dinucSyn_stat[3][8][8];
 }
 double GetDinucSyn_GAGC(){
     return dinucSyn_stat[3][8][9];
 }
 double GetDinucSyn_GAGG(){
     return dinucSyn_stat[3][8][10];
 }
 double GetDinucSyn_GAGT(){
     return dinucSyn_stat[3][8][11];
 }
 double GetDinucSyn_GATA(){
     return dinucSyn_stat[3][8][12];
 }
 double GetDinucSyn_GATC(){
     return dinucSyn_stat[3][8][13];
 }
 double GetDinucSyn_GATG(){
     return dinucSyn_stat[3][8][14];
 }
 double GetDinucSyn_GATT(){
     return dinucSyn_stat[3][8][15];
 }
 double GetDinucSyn_GCAA(){
     return dinucSyn_stat[3][9][0];
 }
 double GetDinucSyn_GCAC(){
     return dinucSyn_stat[3][9][1];
 }
 double GetDinucSyn_GCAG(){
     return dinucSyn_stat[3][9][2];
 }
 double GetDinucSyn_GCAT(){
     return dinucSyn_stat[3][9][3];
 }
 double GetDinucSyn_GCCA(){
     return dinucSyn_stat[3][9][4];
 }
 double GetDinucSyn_GCCC(){
     return dinucSyn_stat[3][9][5];
 }
 double GetDinucSyn_GCCG(){
     return dinucSyn_stat[3][9][6];
 }
 double GetDinucSyn_GCCT(){
     return dinucSyn_stat[3][9][7];
 }
 double GetDinucSyn_GCGA(){
     return dinucSyn_stat[3][9][8];
 }
 double GetDinucSyn_GCGC(){
     return dinucSyn_stat[3][9][9];
 }
 double GetDinucSyn_GCGG(){
     return dinucSyn_stat[3][9][10];
 }
 double GetDinucSyn_GCGT(){
     return dinucSyn_stat[3][9][11];
 }
 double GetDinucSyn_GCTA(){
     return dinucSyn_stat[3][9][12];
 }
 double GetDinucSyn_GCTC(){
     return dinucSyn_stat[3][9][13];
 }
 double GetDinucSyn_GCTG(){
     return dinucSyn_stat[3][9][14];
 }
 double GetDinucSyn_GCTT(){
     return dinucSyn_stat[3][9][15];
 }
 double GetDinucSyn_GGAA(){
     return dinucSyn_stat[3][10][0];
 }
 double GetDinucSyn_GGAC(){
     return dinucSyn_stat[3][10][1];
 }
 double GetDinucSyn_GGAG(){
     return dinucSyn_stat[3][10][2];
 }
 double GetDinucSyn_GGAT(){
     return dinucSyn_stat[3][10][3];
 }
 double GetDinucSyn_GGCA(){
     return dinucSyn_stat[3][10][4];
 }
 double GetDinucSyn_GGCC(){
     return dinucSyn_stat[3][10][5];
 }
 double GetDinucSyn_GGCG(){
     return dinucSyn_stat[3][10][6];
 }
 double GetDinucSyn_GGCT(){
     return dinucSyn_stat[3][10][7];
 }
 double GetDinucSyn_GGGA(){
     return dinucSyn_stat[3][10][8];
 }
 double GetDinucSyn_GGGC(){
     return dinucSyn_stat[3][10][9];
 }
 double GetDinucSyn_GGGG(){
     return dinucSyn_stat[3][10][10];
 }
 double GetDinucSyn_GGGT(){
     return dinucSyn_stat[3][10][11];
 }
 double GetDinucSyn_GGTA(){
     return dinucSyn_stat[3][10][12];
 }
 double GetDinucSyn_GGTC(){
     return dinucSyn_stat[3][10][13];
 }
 double GetDinucSyn_GGTG(){
     return dinucSyn_stat[3][10][14];
 }
 double GetDinucSyn_GGTT(){
     return dinucSyn_stat[3][10][15];
 }
 double GetDinucSyn_GTAA(){
     return dinucSyn_stat[3][11][0];
 }
 double GetDinucSyn_GTAC(){
     return dinucSyn_stat[3][11][1];
 }
 double GetDinucSyn_GTAG(){
     return dinucSyn_stat[3][11][2];
 }
 double GetDinucSyn_GTAT(){
     return dinucSyn_stat[3][11][3];
 }
 double GetDinucSyn_GTCA(){
     return dinucSyn_stat[3][11][4];
 }
 double GetDinucSyn_GTCC(){
     return dinucSyn_stat[3][11][5];
 }
 double GetDinucSyn_GTCG(){
     return dinucSyn_stat[3][11][6];
 }
 double GetDinucSyn_GTCT(){
     return dinucSyn_stat[3][11][7];
 }
 double GetDinucSyn_GTGA(){
     return dinucSyn_stat[3][11][8];
 }
 double GetDinucSyn_GTGC(){
     return dinucSyn_stat[3][11][9];
 }
 double GetDinucSyn_GTGG(){
     return dinucSyn_stat[3][11][10];
 }
 double GetDinucSyn_GTGT(){
     return dinucSyn_stat[3][11][11];
 }
 double GetDinucSyn_GTTA(){
     return dinucSyn_stat[3][11][12];
 }
 double GetDinucSyn_GTTC(){
     return dinucSyn_stat[3][11][13];
 }
 double GetDinucSyn_GTTG(){
     return dinucSyn_stat[3][11][14];
 }
 double GetDinucSyn_GTTT(){
     return dinucSyn_stat[3][11][15];
 }
 double GetDinucSyn_TAAA(){
     return dinucSyn_stat[3][12][0];
 }
 double GetDinucSyn_TAAC(){
     return dinucSyn_stat[3][12][1];
 }
 double GetDinucSyn_TAAG(){
     return dinucSyn_stat[3][12][2];
 }
 double GetDinucSyn_TAAT(){
     return dinucSyn_stat[3][12][3];
 }
 double GetDinucSyn_TACA(){
     return dinucSyn_stat[3][12][4];
 }
 double GetDinucSyn_TACC(){
     return dinucSyn_stat[3][12][5];
 }
 double GetDinucSyn_TACG(){
     return dinucSyn_stat[3][12][6];
 }
 double GetDinucSyn_TACT(){
     return dinucSyn_stat[3][12][7];
 }
 double GetDinucSyn_TAGA(){
     return dinucSyn_stat[3][12][8];
 }
 double GetDinucSyn_TAGC(){
     return dinucSyn_stat[3][12][9];
 }
 double GetDinucSyn_TAGG(){
     return dinucSyn_stat[3][12][10];
 }
 double GetDinucSyn_TAGT(){
     return dinucSyn_stat[3][12][11];
 }
 double GetDinucSyn_TATA(){
     return dinucSyn_stat[3][12][12];
 }
 double GetDinucSyn_TATC(){
     return dinucSyn_stat[3][12][13];
 }
 double GetDinucSyn_TATG(){
     return dinucSyn_stat[3][12][14];
 }
 double GetDinucSyn_TATT(){
     return dinucSyn_stat[3][12][15];
 }
 double GetDinucSyn_TCAA(){
     return dinucSyn_stat[3][13][0];
 }
 double GetDinucSyn_TCAC(){
     return dinucSyn_stat[3][13][1];
 }
 double GetDinucSyn_TCAG(){
     return dinucSyn_stat[3][13][2];
 }
 double GetDinucSyn_TCAT(){
     return dinucSyn_stat[3][13][3];
 }
 double GetDinucSyn_TCCA(){
     return dinucSyn_stat[3][13][4];
 }
 double GetDinucSyn_TCCC(){
     return dinucSyn_stat[3][13][5];
 }
 double GetDinucSyn_TCCG(){
     return dinucSyn_stat[3][13][6];
 }
 double GetDinucSyn_TCCT(){
     return dinucSyn_stat[3][13][7];
 }
 double GetDinucSyn_TCGA(){
     return dinucSyn_stat[3][13][8];
 }
 double GetDinucSyn_TCGC(){
     return dinucSyn_stat[3][13][9];
 }
 double GetDinucSyn_TCGG(){
     return dinucSyn_stat[3][13][10];
 }
 double GetDinucSyn_TCGT(){
     return dinucSyn_stat[3][13][11];
 }
 double GetDinucSyn_TCTA(){
     return dinucSyn_stat[3][13][12];
 }
 double GetDinucSyn_TCTC(){
     return dinucSyn_stat[3][13][13];
 }
 double GetDinucSyn_TCTG(){
     return dinucSyn_stat[3][13][14];
 }
 double GetDinucSyn_TCTT(){
     return dinucSyn_stat[3][13][15];
 }
 double GetDinucSyn_TGAA(){
     return dinucSyn_stat[3][14][0];
 }
 double GetDinucSyn_TGAC(){
     return dinucSyn_stat[3][14][1];
 }
 double GetDinucSyn_TGAG(){
     return dinucSyn_stat[3][14][2];
 }
 double GetDinucSyn_TGAT(){
     return dinucSyn_stat[3][14][3];
 }
 double GetDinucSyn_TGCA(){
     return dinucSyn_stat[3][14][4];
 }
 double GetDinucSyn_TGCC(){
     return dinucSyn_stat[3][14][5];
 }
 double GetDinucSyn_TGCG(){
     return dinucSyn_stat[3][14][6];
 }
 double GetDinucSyn_TGCT(){
     return dinucSyn_stat[3][14][7];
 }
 double GetDinucSyn_TGGA(){
     return dinucSyn_stat[3][14][8];
 }
 double GetDinucSyn_TGGC(){
     return dinucSyn_stat[3][14][9];
 }
 double GetDinucSyn_TGGG(){
     return dinucSyn_stat[3][14][10];
 }
 double GetDinucSyn_TGGT(){
     return dinucSyn_stat[3][14][11];
 }
 double GetDinucSyn_TGTA(){
     return dinucSyn_stat[3][14][12];
 }
 double GetDinucSyn_TGTC(){
     return dinucSyn_stat[3][14][13];
 }
 double GetDinucSyn_TGTG(){
     return dinucSyn_stat[3][14][14];
 }
 double GetDinucSyn_TGTT(){
     return dinucSyn_stat[3][14][15];
 }
 double GetDinucSyn_TTAA(){
     return dinucSyn_stat[3][15][0];
 }
 double GetDinucSyn_TTAC(){
     return dinucSyn_stat[3][15][1];
 }
 double GetDinucSyn_TTAG(){
     return dinucSyn_stat[3][15][2];
 }
 double GetDinucSyn_TTAT(){
     return dinucSyn_stat[3][15][3];
 }
 double GetDinucSyn_TTCA(){
     return dinucSyn_stat[3][15][4];
 }
 double GetDinucSyn_TTCC(){
     return dinucSyn_stat[3][15][5];
 }
 double GetDinucSyn_TTCG(){
     return dinucSyn_stat[3][15][6];
 }
 double GetDinucSyn_TTCT(){
     return dinucSyn_stat[3][15][7];
 }
 double GetDinucSyn_TTGA(){
     return dinucSyn_stat[3][15][8];
 }
 double GetDinucSyn_TTGC(){
     return dinucSyn_stat[3][15][9];
 }
 double GetDinucSyn_TTGG(){
     return dinucSyn_stat[3][15][10];
 }
 double GetDinucSyn_TTGT(){
     return dinucSyn_stat[3][15][11];
 }
 double GetDinucSyn_TTTA(){
     return dinucSyn_stat[3][15][12];
 }
 double GetDinucSyn_TTTC(){
     return dinucSyn_stat[3][15][13];
 }
 double GetDinucSyn_TTTG(){
     return dinucSyn_stat[3][15][14];
 }
 double GetDinucSyn_TTTT(){
     return dinucSyn_stat[3][15][15];
 }
 double GetDinucSyn12_AAAA(){
     return dinucSyn_stat[0][0][0];
 }
 double GetDinucSyn12_AAAC(){
     return dinucSyn_stat[0][0][1];
 }
 double GetDinucSyn12_AAAG(){
     return dinucSyn_stat[0][0][2];
 }
 double GetDinucSyn12_AAAT(){
     return dinucSyn_stat[0][0][3];
 }
 double GetDinucSyn12_AACA(){
     return dinucSyn_stat[0][0][4];
 }
 double GetDinucSyn12_AACC(){
     return dinucSyn_stat[0][0][5];
 }
 double GetDinucSyn12_AACG(){
     return dinucSyn_stat[0][0][6];
 }
 double GetDinucSyn12_AACT(){
     return dinucSyn_stat[0][0][7];
 }
 double GetDinucSyn12_AAGA(){
     return dinucSyn_stat[0][0][8];
 }
 double GetDinucSyn12_AAGC(){
     return dinucSyn_stat[0][0][9];
 }
 double GetDinucSyn12_AAGG(){
     return dinucSyn_stat[0][0][10];
 }
 double GetDinucSyn12_AAGT(){
     return dinucSyn_stat[0][0][11];
 }
 double GetDinucSyn12_AATA(){
     return dinucSyn_stat[0][0][12];
 }
 double GetDinucSyn12_AATC(){
     return dinucSyn_stat[0][0][13];
 }
 double GetDinucSyn12_AATG(){
     return dinucSyn_stat[0][0][14];
 }
 double GetDinucSyn12_AATT(){
     return dinucSyn_stat[0][0][15];
 }
 double GetDinucSyn12_ACAA(){
     return dinucSyn_stat[0][1][0];
 }
 double GetDinucSyn12_ACAC(){
     return dinucSyn_stat[0][1][1];
 }
 double GetDinucSyn12_ACAG(){
     return dinucSyn_stat[0][1][2];
 }
 double GetDinucSyn12_ACAT(){
     return dinucSyn_stat[0][1][3];
 }
 double GetDinucSyn12_ACCA(){
     return dinucSyn_stat[0][1][4];
 }
 double GetDinucSyn12_ACCC(){
     return dinucSyn_stat[0][1][5];
 }
 double GetDinucSyn12_ACCG(){
     return dinucSyn_stat[0][1][6];
 }
 double GetDinucSyn12_ACCT(){
     return dinucSyn_stat[0][1][7];
 }
 double GetDinucSyn12_ACGA(){
     return dinucSyn_stat[0][1][8];
 }
 double GetDinucSyn12_ACGC(){
     return dinucSyn_stat[0][1][9];
 }
 double GetDinucSyn12_ACGG(){
     return dinucSyn_stat[0][1][10];
 }
 double GetDinucSyn12_ACGT(){
     return dinucSyn_stat[0][1][11];
 }
 double GetDinucSyn12_ACTA(){
     return dinucSyn_stat[0][1][12];
 }
 double GetDinucSyn12_ACTC(){
     return dinucSyn_stat[0][1][13];
 }
 double GetDinucSyn12_ACTG(){
     return dinucSyn_stat[0][1][14];
 }
 double GetDinucSyn12_ACTT(){
     return dinucSyn_stat[0][1][15];
 }
 double GetDinucSyn12_AGAA(){
     return dinucSyn_stat[0][2][0];
 }
 double GetDinucSyn12_AGAC(){
     return dinucSyn_stat[0][2][1];
 }
 double GetDinucSyn12_AGAG(){
     return dinucSyn_stat[0][2][2];
 }
 double GetDinucSyn12_AGAT(){
     return dinucSyn_stat[0][2][3];
 }
 double GetDinucSyn12_AGCA(){
     return dinucSyn_stat[0][2][4];
 }
 double GetDinucSyn12_AGCC(){
     return dinucSyn_stat[0][2][5];
 }
 double GetDinucSyn12_AGCG(){
     return dinucSyn_stat[0][2][6];
 }
 double GetDinucSyn12_AGCT(){
     return dinucSyn_stat[0][2][7];
 }
 double GetDinucSyn12_AGGA(){
     return dinucSyn_stat[0][2][8];
 }
 double GetDinucSyn12_AGGC(){
     return dinucSyn_stat[0][2][9];
 }
 double GetDinucSyn12_AGGG(){
     return dinucSyn_stat[0][2][10];
 }
 double GetDinucSyn12_AGGT(){
     return dinucSyn_stat[0][2][11];
 }
 double GetDinucSyn12_AGTA(){
     return dinucSyn_stat[0][2][12];
 }
 double GetDinucSyn12_AGTC(){
     return dinucSyn_stat[0][2][13];
 }
 double GetDinucSyn12_AGTG(){
     return dinucSyn_stat[0][2][14];
 }
 double GetDinucSyn12_AGTT(){
     return dinucSyn_stat[0][2][15];
 }
 double GetDinucSyn12_ATAA(){
     return dinucSyn_stat[0][3][0];
 }
 double GetDinucSyn12_ATAC(){
     return dinucSyn_stat[0][3][1];
 }
 double GetDinucSyn12_ATAG(){
     return dinucSyn_stat[0][3][2];
 }
 double GetDinucSyn12_ATAT(){
     return dinucSyn_stat[0][3][3];
 }
 double GetDinucSyn12_ATCA(){
     return dinucSyn_stat[0][3][4];
 }
 double GetDinucSyn12_ATCC(){
     return dinucSyn_stat[0][3][5];
 }
 double GetDinucSyn12_ATCG(){
     return dinucSyn_stat[0][3][6];
 }
 double GetDinucSyn12_ATCT(){
     return dinucSyn_stat[0][3][7];
 }
 double GetDinucSyn12_ATGA(){
     return dinucSyn_stat[0][3][8];
 }
 double GetDinucSyn12_ATGC(){
     return dinucSyn_stat[0][3][9];
 }
 double GetDinucSyn12_ATGG(){
     return dinucSyn_stat[0][3][10];
 }
 double GetDinucSyn12_ATGT(){
     return dinucSyn_stat[0][3][11];
 }
 double GetDinucSyn12_ATTA(){
     return dinucSyn_stat[0][3][12];
 }
 double GetDinucSyn12_ATTC(){
     return dinucSyn_stat[0][3][13];
 }
 double GetDinucSyn12_ATTG(){
     return dinucSyn_stat[0][3][14];
 }
 double GetDinucSyn12_ATTT(){
     return dinucSyn_stat[0][3][15];
 }
 double GetDinucSyn12_CAAA(){
     return dinucSyn_stat[0][4][0];
 }
 double GetDinucSyn12_CAAC(){
     return dinucSyn_stat[0][4][1];
 }
 double GetDinucSyn12_CAAG(){
     return dinucSyn_stat[0][4][2];
 }
 double GetDinucSyn12_CAAT(){
     return dinucSyn_stat[0][4][3];
 }
 double GetDinucSyn12_CACA(){
     return dinucSyn_stat[0][4][4];
 }
 double GetDinucSyn12_CACC(){
     return dinucSyn_stat[0][4][5];
 }
 double GetDinucSyn12_CACG(){
     return dinucSyn_stat[0][4][6];
 }
 double GetDinucSyn12_CACT(){
     return dinucSyn_stat[0][4][7];
 }
 double GetDinucSyn12_CAGA(){
     return dinucSyn_stat[0][4][8];
 }
 double GetDinucSyn12_CAGC(){
     return dinucSyn_stat[0][4][9];
 }
 double GetDinucSyn12_CAGG(){
     return dinucSyn_stat[0][4][10];
 }
 double GetDinucSyn12_CAGT(){
     return dinucSyn_stat[0][4][11];
 }
 double GetDinucSyn12_CATA(){
     return dinucSyn_stat[0][4][12];
 }
 double GetDinucSyn12_CATC(){
     return dinucSyn_stat[0][4][13];
 }
 double GetDinucSyn12_CATG(){
     return dinucSyn_stat[0][4][14];
 }
 double GetDinucSyn12_CATT(){
     return dinucSyn_stat[0][4][15];
 }
 double GetDinucSyn12_CCAA(){
     return dinucSyn_stat[0][5][0];
 }
 double GetDinucSyn12_CCAC(){
     return dinucSyn_stat[0][5][1];
 }
 double GetDinucSyn12_CCAG(){
     return dinucSyn_stat[0][5][2];
 }
 double GetDinucSyn12_CCAT(){
     return dinucSyn_stat[0][5][3];
 }
 double GetDinucSyn12_CCCA(){
     return dinucSyn_stat[0][5][4];
 }
 double GetDinucSyn12_CCCC(){
     return dinucSyn_stat[0][5][5];
 }
 double GetDinucSyn12_CCCG(){
     return dinucSyn_stat[0][5][6];
 }
 double GetDinucSyn12_CCCT(){
     return dinucSyn_stat[0][5][7];
 }
 double GetDinucSyn12_CCGA(){
     return dinucSyn_stat[0][5][8];
 }
 double GetDinucSyn12_CCGC(){
     return dinucSyn_stat[0][5][9];
 }
 double GetDinucSyn12_CCGG(){
     return dinucSyn_stat[0][5][10];
 }
 double GetDinucSyn12_CCGT(){
     return dinucSyn_stat[0][5][11];
 }
 double GetDinucSyn12_CCTA(){
     return dinucSyn_stat[0][5][12];
 }
 double GetDinucSyn12_CCTC(){
     return dinucSyn_stat[0][5][13];
 }
 double GetDinucSyn12_CCTG(){
     return dinucSyn_stat[0][5][14];
 }
 double GetDinucSyn12_CCTT(){
     return dinucSyn_stat[0][5][15];
 }
 double GetDinucSyn12_CGAA(){
     return dinucSyn_stat[0][6][0];
 }
 double GetDinucSyn12_CGAC(){
     return dinucSyn_stat[0][6][1];
 }
 double GetDinucSyn12_CGAG(){
     return dinucSyn_stat[0][6][2];
 }
 double GetDinucSyn12_CGAT(){
     return dinucSyn_stat[0][6][3];
 }
 double GetDinucSyn12_CGCA(){
     return dinucSyn_stat[0][6][4];
 }
 double GetDinucSyn12_CGCC(){
     return dinucSyn_stat[0][6][5];
 }
 double GetDinucSyn12_CGCG(){
     return dinucSyn_stat[0][6][6];
 }
 double GetDinucSyn12_CGCT(){
     return dinucSyn_stat[0][6][7];
 }
 double GetDinucSyn12_CGGA(){
     return dinucSyn_stat[0][6][8];
 }
 double GetDinucSyn12_CGGC(){
     return dinucSyn_stat[0][6][9];
 }
 double GetDinucSyn12_CGGG(){
     return dinucSyn_stat[0][6][10];
 }
 double GetDinucSyn12_CGGT(){
     return dinucSyn_stat[0][6][11];
 }
 double GetDinucSyn12_CGTA(){
     return dinucSyn_stat[0][6][12];
 }
 double GetDinucSyn12_CGTC(){
     return dinucSyn_stat[0][6][13];
 }
 double GetDinucSyn12_CGTG(){
     return dinucSyn_stat[0][6][14];
 }
 double GetDinucSyn12_CGTT(){
     return dinucSyn_stat[0][6][15];
 }
 double GetDinucSyn12_CTAA(){
     return dinucSyn_stat[0][7][0];
 }
 double GetDinucSyn12_CTAC(){
     return dinucSyn_stat[0][7][1];
 }
 double GetDinucSyn12_CTAG(){
     return dinucSyn_stat[0][7][2];
 }
 double GetDinucSyn12_CTAT(){
     return dinucSyn_stat[0][7][3];
 }
 double GetDinucSyn12_CTCA(){
     return dinucSyn_stat[0][7][4];
 }
 double GetDinucSyn12_CTCC(){
     return dinucSyn_stat[0][7][5];
 }
 double GetDinucSyn12_CTCG(){
     return dinucSyn_stat[0][7][6];
 }
 double GetDinucSyn12_CTCT(){
     return dinucSyn_stat[0][7][7];
 }
 double GetDinucSyn12_CTGA(){
     return dinucSyn_stat[0][7][8];
 }
 double GetDinucSyn12_CTGC(){
     return dinucSyn_stat[0][7][9];
 }
 double GetDinucSyn12_CTGG(){
     return dinucSyn_stat[0][7][10];
 }
 double GetDinucSyn12_CTGT(){
     return dinucSyn_stat[0][7][11];
 }
 double GetDinucSyn12_CTTA(){
     return dinucSyn_stat[0][7][12];
 }
 double GetDinucSyn12_CTTC(){
     return dinucSyn_stat[0][7][13];
 }
 double GetDinucSyn12_CTTG(){
     return dinucSyn_stat[0][7][14];
 }
 double GetDinucSyn12_CTTT(){
     return dinucSyn_stat[0][7][15];
 }
 double GetDinucSyn12_GAAA(){
     return dinucSyn_stat[0][8][0];
 }
 double GetDinucSyn12_GAAC(){
     return dinucSyn_stat[0][8][1];
 }
 double GetDinucSyn12_GAAG(){
     return dinucSyn_stat[0][8][2];
 }
 double GetDinucSyn12_GAAT(){
     return dinucSyn_stat[0][8][3];
 }
 double GetDinucSyn12_GACA(){
     return dinucSyn_stat[0][8][4];
 }
 double GetDinucSyn12_GACC(){
     return dinucSyn_stat[0][8][5];
 }
 double GetDinucSyn12_GACG(){
     return dinucSyn_stat[0][8][6];
 }
 double GetDinucSyn12_GACT(){
     return dinucSyn_stat[0][8][7];
 }
 double GetDinucSyn12_GAGA(){
     return dinucSyn_stat[0][8][8];
 }
 double GetDinucSyn12_GAGC(){
     return dinucSyn_stat[0][8][9];
 }
 double GetDinucSyn12_GAGG(){
     return dinucSyn_stat[0][8][10];
 }
 double GetDinucSyn12_GAGT(){
     return dinucSyn_stat[0][8][11];
 }
 double GetDinucSyn12_GATA(){
     return dinucSyn_stat[0][8][12];
 }
 double GetDinucSyn12_GATC(){
     return dinucSyn_stat[0][8][13];
 }
 double GetDinucSyn12_GATG(){
     return dinucSyn_stat[0][8][14];
 }
 double GetDinucSyn12_GATT(){
     return dinucSyn_stat[0][8][15];
 }
 double GetDinucSyn12_GCAA(){
     return dinucSyn_stat[0][9][0];
 }
 double GetDinucSyn12_GCAC(){
     return dinucSyn_stat[0][9][1];
 }
 double GetDinucSyn12_GCAG(){
     return dinucSyn_stat[0][9][2];
 }
 double GetDinucSyn12_GCAT(){
     return dinucSyn_stat[0][9][3];
 }
 double GetDinucSyn12_GCCA(){
     return dinucSyn_stat[0][9][4];
 }
 double GetDinucSyn12_GCCC(){
     return dinucSyn_stat[0][9][5];
 }
 double GetDinucSyn12_GCCG(){
     return dinucSyn_stat[0][9][6];
 }
 double GetDinucSyn12_GCCT(){
     return dinucSyn_stat[0][9][7];
 }
 double GetDinucSyn12_GCGA(){
     return dinucSyn_stat[0][9][8];
 }
 double GetDinucSyn12_GCGC(){
     return dinucSyn_stat[0][9][9];
 }
 double GetDinucSyn12_GCGG(){
     return dinucSyn_stat[0][9][10];
 }
 double GetDinucSyn12_GCGT(){
     return dinucSyn_stat[0][9][11];
 }
 double GetDinucSyn12_GCTA(){
     return dinucSyn_stat[0][9][12];
 }
 double GetDinucSyn12_GCTC(){
     return dinucSyn_stat[0][9][13];
 }
 double GetDinucSyn12_GCTG(){
     return dinucSyn_stat[0][9][14];
 }
 double GetDinucSyn12_GCTT(){
     return dinucSyn_stat[0][9][15];
 }
 double GetDinucSyn12_GGAA(){
     return dinucSyn_stat[0][10][0];
 }
 double GetDinucSyn12_GGAC(){
     return dinucSyn_stat[0][10][1];
 }
 double GetDinucSyn12_GGAG(){
     return dinucSyn_stat[0][10][2];
 }
 double GetDinucSyn12_GGAT(){
     return dinucSyn_stat[0][10][3];
 }
 double GetDinucSyn12_GGCA(){
     return dinucSyn_stat[0][10][4];
 }
 double GetDinucSyn12_GGCC(){
     return dinucSyn_stat[0][10][5];
 }
 double GetDinucSyn12_GGCG(){
     return dinucSyn_stat[0][10][6];
 }
 double GetDinucSyn12_GGCT(){
     return dinucSyn_stat[0][10][7];
 }
 double GetDinucSyn12_GGGA(){
     return dinucSyn_stat[0][10][8];
 }
 double GetDinucSyn12_GGGC(){
     return dinucSyn_stat[0][10][9];
 }
 double GetDinucSyn12_GGGG(){
     return dinucSyn_stat[0][10][10];
 }
 double GetDinucSyn12_GGGT(){
     return dinucSyn_stat[0][10][11];
 }
 double GetDinucSyn12_GGTA(){
     return dinucSyn_stat[0][10][12];
 }
 double GetDinucSyn12_GGTC(){
     return dinucSyn_stat[0][10][13];
 }
 double GetDinucSyn12_GGTG(){
     return dinucSyn_stat[0][10][14];
 }
 double GetDinucSyn12_GGTT(){
     return dinucSyn_stat[0][10][15];
 }
 double GetDinucSyn12_GTAA(){
     return dinucSyn_stat[0][11][0];
 }
 double GetDinucSyn12_GTAC(){
     return dinucSyn_stat[0][11][1];
 }
 double GetDinucSyn12_GTAG(){
     return dinucSyn_stat[0][11][2];
 }
 double GetDinucSyn12_GTAT(){
     return dinucSyn_stat[0][11][3];
 }
 double GetDinucSyn12_GTCA(){
     return dinucSyn_stat[0][11][4];
 }
 double GetDinucSyn12_GTCC(){
     return dinucSyn_stat[0][11][5];
 }
 double GetDinucSyn12_GTCG(){
     return dinucSyn_stat[0][11][6];
 }
 double GetDinucSyn12_GTCT(){
     return dinucSyn_stat[0][11][7];
 }
 double GetDinucSyn12_GTGA(){
     return dinucSyn_stat[0][11][8];
 }
 double GetDinucSyn12_GTGC(){
     return dinucSyn_stat[0][11][9];
 }
 double GetDinucSyn12_GTGG(){
     return dinucSyn_stat[0][11][10];
 }
 double GetDinucSyn12_GTGT(){
     return dinucSyn_stat[0][11][11];
 }
 double GetDinucSyn12_GTTA(){
     return dinucSyn_stat[0][11][12];
 }
 double GetDinucSyn12_GTTC(){
     return dinucSyn_stat[0][11][13];
 }
 double GetDinucSyn12_GTTG(){
     return dinucSyn_stat[0][11][14];
 }
 double GetDinucSyn12_GTTT(){
     return dinucSyn_stat[0][11][15];
 }
 double GetDinucSyn12_TAAA(){
     return dinucSyn_stat[0][12][0];
 }
 double GetDinucSyn12_TAAC(){
     return dinucSyn_stat[0][12][1];
 }
 double GetDinucSyn12_TAAG(){
     return dinucSyn_stat[0][12][2];
 }
 double GetDinucSyn12_TAAT(){
     return dinucSyn_stat[0][12][3];
 }
 double GetDinucSyn12_TACA(){
     return dinucSyn_stat[0][12][4];
 }
 double GetDinucSyn12_TACC(){
     return dinucSyn_stat[0][12][5];
 }
 double GetDinucSyn12_TACG(){
     return dinucSyn_stat[0][12][6];
 }
 double GetDinucSyn12_TACT(){
     return dinucSyn_stat[0][12][7];
 }
 double GetDinucSyn12_TAGA(){
     return dinucSyn_stat[0][12][8];
 }
 double GetDinucSyn12_TAGC(){
     return dinucSyn_stat[0][12][9];
 }
 double GetDinucSyn12_TAGG(){
     return dinucSyn_stat[0][12][10];
 }
 double GetDinucSyn12_TAGT(){
     return dinucSyn_stat[0][12][11];
 }
 double GetDinucSyn12_TATA(){
     return dinucSyn_stat[0][12][12];
 }
 double GetDinucSyn12_TATC(){
     return dinucSyn_stat[0][12][13];
 }
 double GetDinucSyn12_TATG(){
     return dinucSyn_stat[0][12][14];
 }
 double GetDinucSyn12_TATT(){
     return dinucSyn_stat[0][12][15];
 }
 double GetDinucSyn12_TCAA(){
     return dinucSyn_stat[0][13][0];
 }
 double GetDinucSyn12_TCAC(){
     return dinucSyn_stat[0][13][1];
 }
 double GetDinucSyn12_TCAG(){
     return dinucSyn_stat[0][13][2];
 }
 double GetDinucSyn12_TCAT(){
     return dinucSyn_stat[0][13][3];
 }
 double GetDinucSyn12_TCCA(){
     return dinucSyn_stat[0][13][4];
 }
 double GetDinucSyn12_TCCC(){
     return dinucSyn_stat[0][13][5];
 }
 double GetDinucSyn12_TCCG(){
     return dinucSyn_stat[0][13][6];
 }
 double GetDinucSyn12_TCCT(){
     return dinucSyn_stat[0][13][7];
 }
 double GetDinucSyn12_TCGA(){
     return dinucSyn_stat[0][13][8];
 }
 double GetDinucSyn12_TCGC(){
     return dinucSyn_stat[0][13][9];
 }
 double GetDinucSyn12_TCGG(){
     return dinucSyn_stat[0][13][10];
 }
 double GetDinucSyn12_TCGT(){
     return dinucSyn_stat[0][13][11];
 }
 double GetDinucSyn12_TCTA(){
     return dinucSyn_stat[0][13][12];
 }
 double GetDinucSyn12_TCTC(){
     return dinucSyn_stat[0][13][13];
 }
 double GetDinucSyn12_TCTG(){
     return dinucSyn_stat[0][13][14];
 }
 double GetDinucSyn12_TCTT(){
     return dinucSyn_stat[0][13][15];
 }
 double GetDinucSyn12_TGAA(){
     return dinucSyn_stat[0][14][0];
 }
 double GetDinucSyn12_TGAC(){
     return dinucSyn_stat[0][14][1];
 }
 double GetDinucSyn12_TGAG(){
     return dinucSyn_stat[0][14][2];
 }
 double GetDinucSyn12_TGAT(){
     return dinucSyn_stat[0][14][3];
 }
 double GetDinucSyn12_TGCA(){
     return dinucSyn_stat[0][14][4];
 }
 double GetDinucSyn12_TGCC(){
     return dinucSyn_stat[0][14][5];
 }
 double GetDinucSyn12_TGCG(){
     return dinucSyn_stat[0][14][6];
 }
 double GetDinucSyn12_TGCT(){
     return dinucSyn_stat[0][14][7];
 }
 double GetDinucSyn12_TGGA(){
     return dinucSyn_stat[0][14][8];
 }
 double GetDinucSyn12_TGGC(){
     return dinucSyn_stat[0][14][9];
 }
 double GetDinucSyn12_TGGG(){
     return dinucSyn_stat[0][14][10];
 }
 double GetDinucSyn12_TGGT(){
     return dinucSyn_stat[0][14][11];
 }
 double GetDinucSyn12_TGTA(){
     return dinucSyn_stat[0][14][12];
 }
 double GetDinucSyn12_TGTC(){
     return dinucSyn_stat[0][14][13];
 }
 double GetDinucSyn12_TGTG(){
     return dinucSyn_stat[0][14][14];
 }
 double GetDinucSyn12_TGTT(){
     return dinucSyn_stat[0][14][15];
 }
 double GetDinucSyn12_TTAA(){
     return dinucSyn_stat[0][15][0];
 }
 double GetDinucSyn12_TTAC(){
     return dinucSyn_stat[0][15][1];
 }
 double GetDinucSyn12_TTAG(){
     return dinucSyn_stat[0][15][2];
 }
 double GetDinucSyn12_TTAT(){
     return dinucSyn_stat[0][15][3];
 }
 double GetDinucSyn12_TTCA(){
     return dinucSyn_stat[0][15][4];
 }
 double GetDinucSyn12_TTCC(){
     return dinucSyn_stat[0][15][5];
 }
 double GetDinucSyn12_TTCG(){
     return dinucSyn_stat[0][15][6];
 }
 double GetDinucSyn12_TTCT(){
     return dinucSyn_stat[0][15][7];
 }
 double GetDinucSyn12_TTGA(){
     return dinucSyn_stat[0][15][8];
 }
 double GetDinucSyn12_TTGC(){
     return dinucSyn_stat[0][15][9];
 }
 double GetDinucSyn12_TTGG(){
     return dinucSyn_stat[0][15][10];
 }
 double GetDinucSyn12_TTGT(){
     return dinucSyn_stat[0][15][11];
 }
 double GetDinucSyn12_TTTA(){
     return dinucSyn_stat[0][15][12];
 }
 double GetDinucSyn12_TTTC(){
     return dinucSyn_stat[0][15][13];
 }
 double GetDinucSyn12_TTTG(){
     return dinucSyn_stat[0][15][14];
 }
 double GetDinucSyn12_TTTT(){
     return dinucSyn_stat[0][15][15];
 }

 double GetDinucSyn23_AAAA(){
     return dinucSyn_stat[1][0][0];
 }
 double GetDinucSyn23_AAAC(){
     return dinucSyn_stat[1][0][1];
 }
 double GetDinucSyn23_AAAG(){
     return dinucSyn_stat[1][0][2];
 }
 double GetDinucSyn23_AAAT(){
     return dinucSyn_stat[1][0][3];
 }
 double GetDinucSyn23_AACA(){
     return dinucSyn_stat[1][0][4];
 }
 double GetDinucSyn23_AACC(){
     return dinucSyn_stat[1][0][5];
 }
 double GetDinucSyn23_AACG(){
     return dinucSyn_stat[1][0][6];
 }
 double GetDinucSyn23_AACT(){
     return dinucSyn_stat[1][0][7];
 }
 double GetDinucSyn23_AAGA(){
     return dinucSyn_stat[1][0][8];
 }
 double GetDinucSyn23_AAGC(){
     return dinucSyn_stat[1][0][9];
 }
 double GetDinucSyn23_AAGG(){
     return dinucSyn_stat[1][0][10];
 }
 double GetDinucSyn23_AAGT(){
     return dinucSyn_stat[1][0][11];
 }
 double GetDinucSyn23_AATA(){
     return dinucSyn_stat[1][0][12];
 }
 double GetDinucSyn23_AATC(){
     return dinucSyn_stat[1][0][13];
 }
 double GetDinucSyn23_AATG(){
     return dinucSyn_stat[1][0][14];
 }
 double GetDinucSyn23_AATT(){
     return dinucSyn_stat[1][0][15];
 }
 double GetDinucSyn23_ACAA(){
     return dinucSyn_stat[1][1][0];
 }
 double GetDinucSyn23_ACAC(){
     return dinucSyn_stat[1][1][1];
 }
 double GetDinucSyn23_ACAG(){
     return dinucSyn_stat[1][1][2];
 }
 double GetDinucSyn23_ACAT(){
     return dinucSyn_stat[1][1][3];
 }
 double GetDinucSyn23_ACCA(){
     return dinucSyn_stat[1][1][4];
 }
 double GetDinucSyn23_ACCC(){
     return dinucSyn_stat[1][1][5];
 }
 double GetDinucSyn23_ACCG(){
     return dinucSyn_stat[1][1][6];
 }
 double GetDinucSyn23_ACCT(){
     return dinucSyn_stat[1][1][7];
 }
 double GetDinucSyn23_ACGA(){
     return dinucSyn_stat[1][1][8];
 }
 double GetDinucSyn23_ACGC(){
     return dinucSyn_stat[1][1][9];
 }
 double GetDinucSyn23_ACGG(){
     return dinucSyn_stat[1][1][10];
 }
 double GetDinucSyn23_ACGT(){
     return dinucSyn_stat[1][1][11];
 }
 double GetDinucSyn23_ACTA(){
     return dinucSyn_stat[1][1][12];
 }
 double GetDinucSyn23_ACTC(){
     return dinucSyn_stat[1][1][13];
 }
 double GetDinucSyn23_ACTG(){
     return dinucSyn_stat[1][1][14];
 }
 double GetDinucSyn23_ACTT(){
     return dinucSyn_stat[1][1][15];
 }
 double GetDinucSyn23_AGAA(){
     return dinucSyn_stat[1][2][0];
 }
 double GetDinucSyn23_AGAC(){
     return dinucSyn_stat[1][2][1];
 }
 double GetDinucSyn23_AGAG(){
     return dinucSyn_stat[1][2][2];
 }
 double GetDinucSyn23_AGAT(){
     return dinucSyn_stat[1][2][3];
 }
 double GetDinucSyn23_AGCA(){
     return dinucSyn_stat[1][2][4];
 }
 double GetDinucSyn23_AGCC(){
     return dinucSyn_stat[1][2][5];
 }
 double GetDinucSyn23_AGCG(){
     return dinucSyn_stat[1][2][6];
 }
 double GetDinucSyn23_AGCT(){
     return dinucSyn_stat[1][2][7];
 }
 double GetDinucSyn23_AGGA(){
     return dinucSyn_stat[1][2][8];
 }
 double GetDinucSyn23_AGGC(){
     return dinucSyn_stat[1][2][9];
 }
 double GetDinucSyn23_AGGG(){
     return dinucSyn_stat[1][2][10];
 }
 double GetDinucSyn23_AGGT(){
     return dinucSyn_stat[1][2][11];
 }
 double GetDinucSyn23_AGTA(){
     return dinucSyn_stat[1][2][12];
 }
 double GetDinucSyn23_AGTC(){
     return dinucSyn_stat[1][2][13];
 }
 double GetDinucSyn23_AGTG(){
     return dinucSyn_stat[1][2][14];
 }
 double GetDinucSyn23_AGTT(){
     return dinucSyn_stat[1][2][15];
 }
 double GetDinucSyn23_ATAA(){
     return dinucSyn_stat[1][3][0];
 }
 double GetDinucSyn23_ATAC(){
     return dinucSyn_stat[1][3][1];
 }
 double GetDinucSyn23_ATAG(){
     return dinucSyn_stat[1][3][2];
 }
 double GetDinucSyn23_ATAT(){
     return dinucSyn_stat[1][3][3];
 }
 double GetDinucSyn23_ATCA(){
     return dinucSyn_stat[1][3][4];
 }
 double GetDinucSyn23_ATCC(){
     return dinucSyn_stat[1][3][5];
 }
 double GetDinucSyn23_ATCG(){
     return dinucSyn_stat[1][3][6];
 }
 double GetDinucSyn23_ATCT(){
     return dinucSyn_stat[1][3][7];
 }
 double GetDinucSyn23_ATGA(){
     return dinucSyn_stat[1][3][8];
 }
 double GetDinucSyn23_ATGC(){
     return dinucSyn_stat[1][3][9];
 }
 double GetDinucSyn23_ATGG(){
     return dinucSyn_stat[1][3][10];
 }
 double GetDinucSyn23_ATGT(){
     return dinucSyn_stat[1][3][11];
 }
 double GetDinucSyn23_ATTA(){
     return dinucSyn_stat[1][3][12];
 }
 double GetDinucSyn23_ATTC(){
     return dinucSyn_stat[1][3][13];
 }
 double GetDinucSyn23_ATTG(){
     return dinucSyn_stat[1][3][14];
 }
 double GetDinucSyn23_ATTT(){
     return dinucSyn_stat[1][3][15];
 }
 double GetDinucSyn23_CAAA(){
     return dinucSyn_stat[1][4][0];
 }
 double GetDinucSyn23_CAAC(){
     return dinucSyn_stat[1][4][1];
 }
 double GetDinucSyn23_CAAG(){
     return dinucSyn_stat[1][4][2];
 }
 double GetDinucSyn23_CAAT(){
     return dinucSyn_stat[1][4][3];
 }
 double GetDinucSyn23_CACA(){
     return dinucSyn_stat[1][4][4];
 }
 double GetDinucSyn23_CACC(){
     return dinucSyn_stat[1][4][5];
 }
 double GetDinucSyn23_CACG(){
     return dinucSyn_stat[1][4][6];
 }
 double GetDinucSyn23_CACT(){
     return dinucSyn_stat[1][4][7];
 }
 double GetDinucSyn23_CAGA(){
     return dinucSyn_stat[1][4][8];
 }
 double GetDinucSyn23_CAGC(){
     return dinucSyn_stat[1][4][9];
 }
 double GetDinucSyn23_CAGG(){
     return dinucSyn_stat[1][4][10];
 }
 double GetDinucSyn23_CAGT(){
     return dinucSyn_stat[1][4][11];
 }
 double GetDinucSyn23_CATA(){
     return dinucSyn_stat[1][4][12];
 }
 double GetDinucSyn23_CATC(){
     return dinucSyn_stat[1][4][13];
 }
 double GetDinucSyn23_CATG(){
     return dinucSyn_stat[1][4][14];
 }
 double GetDinucSyn23_CATT(){
     return dinucSyn_stat[1][4][15];
 }
 double GetDinucSyn23_CCAA(){
     return dinucSyn_stat[1][5][0];
 }
 double GetDinucSyn23_CCAC(){
     return dinucSyn_stat[1][5][1];
 }
 double GetDinucSyn23_CCAG(){
     return dinucSyn_stat[1][5][2];
 }
 double GetDinucSyn23_CCAT(){
     return dinucSyn_stat[1][5][3];
 }
 double GetDinucSyn23_CCCA(){
     return dinucSyn_stat[1][5][4];
 }
 double GetDinucSyn23_CCCC(){
     return dinucSyn_stat[1][5][5];
 }
 double GetDinucSyn23_CCCG(){
     return dinucSyn_stat[1][5][6];
 }
 double GetDinucSyn23_CCCT(){
     return dinucSyn_stat[1][5][7];
 }
 double GetDinucSyn23_CCGA(){
     return dinucSyn_stat[1][5][8];
 }
 double GetDinucSyn23_CCGC(){
     return dinucSyn_stat[1][5][9];
 }
 double GetDinucSyn23_CCGG(){
     return dinucSyn_stat[1][5][10];
 }
 double GetDinucSyn23_CCGT(){
     return dinucSyn_stat[1][5][11];
 }
 double GetDinucSyn23_CCTA(){
     return dinucSyn_stat[1][5][12];
 }
 double GetDinucSyn23_CCTC(){
     return dinucSyn_stat[1][5][13];
 }
 double GetDinucSyn23_CCTG(){
     return dinucSyn_stat[1][5][14];
 }
 double GetDinucSyn23_CCTT(){
     return dinucSyn_stat[1][5][15];
 }
 double GetDinucSyn23_CGAA(){
     return dinucSyn_stat[1][6][0];
 }
 double GetDinucSyn23_CGAC(){
     return dinucSyn_stat[1][6][1];
 }
 double GetDinucSyn23_CGAG(){
     return dinucSyn_stat[1][6][2];
 }
 double GetDinucSyn23_CGAT(){
     return dinucSyn_stat[1][6][3];
 }
 double GetDinucSyn23_CGCA(){
     return dinucSyn_stat[1][6][4];
 }
 double GetDinucSyn23_CGCC(){
     return dinucSyn_stat[1][6][5];
 }
 double GetDinucSyn23_CGCG(){
     return dinucSyn_stat[1][6][6];
 }
 double GetDinucSyn23_CGCT(){
     return dinucSyn_stat[1][6][7];
 }
 double GetDinucSyn23_CGGA(){
     return dinucSyn_stat[1][6][8];
 }
 double GetDinucSyn23_CGGC(){
     return dinucSyn_stat[1][6][9];
 }
 double GetDinucSyn23_CGGG(){
     return dinucSyn_stat[1][6][10];
 }
 double GetDinucSyn23_CGGT(){
     return dinucSyn_stat[1][6][11];
 }
 double GetDinucSyn23_CGTA(){
     return dinucSyn_stat[1][6][12];
 }
 double GetDinucSyn23_CGTC(){
     return dinucSyn_stat[1][6][13];
 }
 double GetDinucSyn23_CGTG(){
     return dinucSyn_stat[1][6][14];
 }
 double GetDinucSyn23_CGTT(){
     return dinucSyn_stat[1][6][15];
 }
 double GetDinucSyn23_CTAA(){
     return dinucSyn_stat[1][7][0];
 }
 double GetDinucSyn23_CTAC(){
     return dinucSyn_stat[1][7][1];
 }
 double GetDinucSyn23_CTAG(){
     return dinucSyn_stat[1][7][2];
 }
 double GetDinucSyn23_CTAT(){
     return dinucSyn_stat[1][7][3];
 }
 double GetDinucSyn23_CTCA(){
     return dinucSyn_stat[1][7][4];
 }
 double GetDinucSyn23_CTCC(){
     return dinucSyn_stat[1][7][5];
 }
 double GetDinucSyn23_CTCG(){
     return dinucSyn_stat[1][7][6];
 }
 double GetDinucSyn23_CTCT(){
     return dinucSyn_stat[1][7][7];
 }
 double GetDinucSyn23_CTGA(){
     return dinucSyn_stat[1][7][8];
 }
 double GetDinucSyn23_CTGC(){
     return dinucSyn_stat[1][7][9];
 }
 double GetDinucSyn23_CTGG(){
     return dinucSyn_stat[1][7][10];
 }
 double GetDinucSyn23_CTGT(){
     return dinucSyn_stat[1][7][11];
 }
 double GetDinucSyn23_CTTA(){
     return dinucSyn_stat[1][7][12];
 }
 double GetDinucSyn23_CTTC(){
     return dinucSyn_stat[1][7][13];
 }
 double GetDinucSyn23_CTTG(){
     return dinucSyn_stat[1][7][14];
 }
 double GetDinucSyn23_CTTT(){
     return dinucSyn_stat[1][7][15];
 }
 double GetDinucSyn23_GAAA(){
     return dinucSyn_stat[1][8][0];
 }
 double GetDinucSyn23_GAAC(){
     return dinucSyn_stat[1][8][1];
 }
 double GetDinucSyn23_GAAG(){
     return dinucSyn_stat[1][8][2];
 }
 double GetDinucSyn23_GAAT(){
     return dinucSyn_stat[1][8][3];
 }
 double GetDinucSyn23_GACA(){
     return dinucSyn_stat[1][8][4];
 }
 double GetDinucSyn23_GACC(){
     return dinucSyn_stat[1][8][5];
 }
 double GetDinucSyn23_GACG(){
     return dinucSyn_stat[1][8][6];
 }
 double GetDinucSyn23_GACT(){
     return dinucSyn_stat[1][8][7];
 }
 double GetDinucSyn23_GAGA(){
     return dinucSyn_stat[1][8][8];
 }
 double GetDinucSyn23_GAGC(){
     return dinucSyn_stat[1][8][9];
 }
 double GetDinucSyn23_GAGG(){
     return dinucSyn_stat[1][8][10];
 }
 double GetDinucSyn23_GAGT(){
     return dinucSyn_stat[1][8][11];
 }
 double GetDinucSyn23_GATA(){
     return dinucSyn_stat[1][8][12];
 }
 double GetDinucSyn23_GATC(){
     return dinucSyn_stat[1][8][13];
 }
 double GetDinucSyn23_GATG(){
     return dinucSyn_stat[1][8][14];
 }
 double GetDinucSyn23_GATT(){
     return dinucSyn_stat[1][8][15];
 }
 double GetDinucSyn23_GCAA(){
     return dinucSyn_stat[1][9][0];
 }
 double GetDinucSyn23_GCAC(){
     return dinucSyn_stat[1][9][1];
 }
 double GetDinucSyn23_GCAG(){
     return dinucSyn_stat[1][9][2];
 }
 double GetDinucSyn23_GCAT(){
     return dinucSyn_stat[1][9][3];
 }
 double GetDinucSyn23_GCCA(){
     return dinucSyn_stat[1][9][4];
 }
 double GetDinucSyn23_GCCC(){
     return dinucSyn_stat[1][9][5];
 }
 double GetDinucSyn23_GCCG(){
     return dinucSyn_stat[1][9][6];
 }
 double GetDinucSyn23_GCCT(){
     return dinucSyn_stat[1][9][7];
 }
 double GetDinucSyn23_GCGA(){
     return dinucSyn_stat[1][9][8];
 }
 double GetDinucSyn23_GCGC(){
     return dinucSyn_stat[1][9][9];
 }
 double GetDinucSyn23_GCGG(){
     return dinucSyn_stat[1][9][10];
 }
 double GetDinucSyn23_GCGT(){
     return dinucSyn_stat[1][9][11];
 }
 double GetDinucSyn23_GCTA(){
     return dinucSyn_stat[1][9][12];
 }
 double GetDinucSyn23_GCTC(){
     return dinucSyn_stat[1][9][13];
 }
 double GetDinucSyn23_GCTG(){
     return dinucSyn_stat[1][9][14];
 }
 double GetDinucSyn23_GCTT(){
     return dinucSyn_stat[1][9][15];
 }
 double GetDinucSyn23_GGAA(){
     return dinucSyn_stat[1][10][0];
 }
 double GetDinucSyn23_GGAC(){
     return dinucSyn_stat[1][10][1];
 }
 double GetDinucSyn23_GGAG(){
     return dinucSyn_stat[1][10][2];
 }
 double GetDinucSyn23_GGAT(){
     return dinucSyn_stat[1][10][3];
 }
 double GetDinucSyn23_GGCA(){
     return dinucSyn_stat[1][10][4];
 }
 double GetDinucSyn23_GGCC(){
     return dinucSyn_stat[1][10][5];
 }
 double GetDinucSyn23_GGCG(){
     return dinucSyn_stat[1][10][6];
 }
 double GetDinucSyn23_GGCT(){
     return dinucSyn_stat[1][10][7];
 }
 double GetDinucSyn23_GGGA(){
     return dinucSyn_stat[1][10][8];
 }
 double GetDinucSyn23_GGGC(){
     return dinucSyn_stat[1][10][9];
 }
 double GetDinucSyn23_GGGG(){
     return dinucSyn_stat[1][10][10];
 }
 double GetDinucSyn23_GGGT(){
     return dinucSyn_stat[1][10][11];
 }
 double GetDinucSyn23_GGTA(){
     return dinucSyn_stat[1][10][12];
 }
 double GetDinucSyn23_GGTC(){
     return dinucSyn_stat[1][10][13];
 }
 double GetDinucSyn23_GGTG(){
     return dinucSyn_stat[1][10][14];
 }
 double GetDinucSyn23_GGTT(){
     return dinucSyn_stat[1][10][15];
 }
 double GetDinucSyn23_GTAA(){
     return dinucSyn_stat[1][11][0];
 }
 double GetDinucSyn23_GTAC(){
     return dinucSyn_stat[1][11][1];
 }
 double GetDinucSyn23_GTAG(){
     return dinucSyn_stat[1][11][2];
 }
 double GetDinucSyn23_GTAT(){
     return dinucSyn_stat[1][11][3];
 }
 double GetDinucSyn23_GTCA(){
     return dinucSyn_stat[1][11][4];
 }
 double GetDinucSyn23_GTCC(){
     return dinucSyn_stat[1][11][5];
 }
 double GetDinucSyn23_GTCG(){
     return dinucSyn_stat[1][11][6];
 }
 double GetDinucSyn23_GTCT(){
     return dinucSyn_stat[1][11][7];
 }
 double GetDinucSyn23_GTGA(){
     return dinucSyn_stat[1][11][8];
 }
 double GetDinucSyn23_GTGC(){
     return dinucSyn_stat[1][11][9];
 }
 double GetDinucSyn23_GTGG(){
     return dinucSyn_stat[1][11][10];
 }
 double GetDinucSyn23_GTGT(){
     return dinucSyn_stat[1][11][11];
 }
 double GetDinucSyn23_GTTA(){
     return dinucSyn_stat[1][11][12];
 }
 double GetDinucSyn23_GTTC(){
     return dinucSyn_stat[1][11][13];
 }
 double GetDinucSyn23_GTTG(){
     return dinucSyn_stat[1][11][14];
 }
 double GetDinucSyn23_GTTT(){
     return dinucSyn_stat[1][11][15];
 }
 double GetDinucSyn23_TAAA(){
     return dinucSyn_stat[1][12][0];
 }
 double GetDinucSyn23_TAAC(){
     return dinucSyn_stat[1][12][1];
 }
 double GetDinucSyn23_TAAG(){
     return dinucSyn_stat[1][12][2];
 }
 double GetDinucSyn23_TAAT(){
     return dinucSyn_stat[1][12][3];
 }
 double GetDinucSyn23_TACA(){
     return dinucSyn_stat[1][12][4];
 }
 double GetDinucSyn23_TACC(){
     return dinucSyn_stat[1][12][5];
 }
 double GetDinucSyn23_TACG(){
     return dinucSyn_stat[1][12][6];
 }
 double GetDinucSyn23_TACT(){
     return dinucSyn_stat[1][12][7];
 }
 double GetDinucSyn23_TAGA(){
     return dinucSyn_stat[1][12][8];
 }
 double GetDinucSyn23_TAGC(){
     return dinucSyn_stat[1][12][9];
 }
 double GetDinucSyn23_TAGG(){
     return dinucSyn_stat[1][12][10];
 }
 double GetDinucSyn23_TAGT(){
     return dinucSyn_stat[1][12][11];
 }
 double GetDinucSyn23_TATA(){
     return dinucSyn_stat[1][12][12];
 }
 double GetDinucSyn23_TATC(){
     return dinucSyn_stat[1][12][13];
 }
 double GetDinucSyn23_TATG(){
     return dinucSyn_stat[1][12][14];
 }
 double GetDinucSyn23_TATT(){
     return dinucSyn_stat[1][12][15];
 }
 double GetDinucSyn23_TCAA(){
     return dinucSyn_stat[1][13][0];
 }
 double GetDinucSyn23_TCAC(){
     return dinucSyn_stat[1][13][1];
 }
 double GetDinucSyn23_TCAG(){
     return dinucSyn_stat[1][13][2];
 }
 double GetDinucSyn23_TCAT(){
     return dinucSyn_stat[1][13][3];
 }
 double GetDinucSyn23_TCCA(){
     return dinucSyn_stat[1][13][4];
 }
 double GetDinucSyn23_TCCC(){
     return dinucSyn_stat[1][13][5];
 }
 double GetDinucSyn23_TCCG(){
     return dinucSyn_stat[1][13][6];
 }
 double GetDinucSyn23_TCCT(){
     return dinucSyn_stat[1][13][7];
 }
 double GetDinucSyn23_TCGA(){
     return dinucSyn_stat[1][13][8];
 }
 double GetDinucSyn23_TCGC(){
     return dinucSyn_stat[1][13][9];
 }
 double GetDinucSyn23_TCGG(){
     return dinucSyn_stat[1][13][10];
 }
 double GetDinucSyn23_TCGT(){
     return dinucSyn_stat[1][13][11];
 }
 double GetDinucSyn23_TCTA(){
     return dinucSyn_stat[1][13][12];
 }
 double GetDinucSyn23_TCTC(){
     return dinucSyn_stat[1][13][13];
 }
 double GetDinucSyn23_TCTG(){
     return dinucSyn_stat[1][13][14];
 }
 double GetDinucSyn23_TCTT(){
     return dinucSyn_stat[1][13][15];
 }
 double GetDinucSyn23_TGAA(){
     return dinucSyn_stat[1][14][0];
 }
 double GetDinucSyn23_TGAC(){
     return dinucSyn_stat[1][14][1];
 }
 double GetDinucSyn23_TGAG(){
     return dinucSyn_stat[1][14][2];
 }
 double GetDinucSyn23_TGAT(){
     return dinucSyn_stat[1][14][3];
 }
 double GetDinucSyn23_TGCA(){
     return dinucSyn_stat[1][14][4];
 }
 double GetDinucSyn23_TGCC(){
     return dinucSyn_stat[1][14][5];
 }
 double GetDinucSyn23_TGCG(){
     return dinucSyn_stat[1][14][6];
 }
 double GetDinucSyn23_TGCT(){
     return dinucSyn_stat[1][14][7];
 }
 double GetDinucSyn23_TGGA(){
     return dinucSyn_stat[1][14][8];
 }
 double GetDinucSyn23_TGGC(){
     return dinucSyn_stat[1][14][9];
 }
 double GetDinucSyn23_TGGG(){
     return dinucSyn_stat[1][14][10];
 }
 double GetDinucSyn23_TGGT(){
     return dinucSyn_stat[1][14][11];
 }
 double GetDinucSyn23_TGTA(){
     return dinucSyn_stat[1][14][12];
 }
 double GetDinucSyn23_TGTC(){
     return dinucSyn_stat[1][14][13];
 }
 double GetDinucSyn23_TGTG(){
     return dinucSyn_stat[1][14][14];
 }
 double GetDinucSyn23_TGTT(){
     return dinucSyn_stat[1][14][15];
 }
 double GetDinucSyn23_TTAA(){
     return dinucSyn_stat[1][15][0];
 }
 double GetDinucSyn23_TTAC(){
     return dinucSyn_stat[1][15][1];
 }
 double GetDinucSyn23_TTAG(){
     return dinucSyn_stat[1][15][2];
 }
 double GetDinucSyn23_TTAT(){
     return dinucSyn_stat[1][15][3];
 }
 double GetDinucSyn23_TTCA(){
     return dinucSyn_stat[1][15][4];
 }
 double GetDinucSyn23_TTCC(){
     return dinucSyn_stat[1][15][5];
 }
 double GetDinucSyn23_TTCG(){
     return dinucSyn_stat[1][15][6];
 }
 double GetDinucSyn23_TTCT(){
     return dinucSyn_stat[1][15][7];
 }
 double GetDinucSyn23_TTGA(){
     return dinucSyn_stat[1][15][8];
 }
 double GetDinucSyn23_TTGC(){
     return dinucSyn_stat[1][15][9];
 }
 double GetDinucSyn23_TTGG(){
     return dinucSyn_stat[1][15][10];
 }
 double GetDinucSyn23_TTGT(){
     return dinucSyn_stat[1][15][11];
 }
 double GetDinucSyn23_TTTA(){
     return dinucSyn_stat[1][15][12];
 }
 double GetDinucSyn23_TTTC(){
     return dinucSyn_stat[1][15][13];
 }
 double GetDinucSyn23_TTTG(){
     return dinucSyn_stat[1][15][14];
 }
 double GetDinucSyn23_TTTT(){
     return dinucSyn_stat[1][15][15];
 }

 double GetDinucSyn31_AAAA(){
     return dinucSyn_stat[2][0][0];
  }
  double GetDinucSyn31_AAAC(){
     return dinucSyn_stat[2][0][1];
  }
  double GetDinucSyn31_AAAG(){
     return dinucSyn_stat[2][0][2];
  }
  double GetDinucSyn31_AAAT(){
     return dinucSyn_stat[2][0][3];
  }
  double GetDinucSyn31_AACA(){
     return dinucSyn_stat[2][0][4];
  }
  double GetDinucSyn31_AACC(){
     return dinucSyn_stat[2][0][5];
  }
  double GetDinucSyn31_AACG(){
     return dinucSyn_stat[2][0][6];
  }
  double GetDinucSyn31_AACT(){
     return dinucSyn_stat[2][0][7];
  }
  double GetDinucSyn31_AAGA(){
     return dinucSyn_stat[2][0][8];
  }
  double GetDinucSyn31_AAGC(){
     return dinucSyn_stat[2][0][9];
  }
  double GetDinucSyn31_AAGG(){
     return dinucSyn_stat[2][0][10];
  }
  double GetDinucSyn31_AAGT(){
     return dinucSyn_stat[2][0][11];
  }
  double GetDinucSyn31_AATA(){
     return dinucSyn_stat[2][0][12];
  }
  double GetDinucSyn31_AATC(){
     return dinucSyn_stat[2][0][13];
  }
  double GetDinucSyn31_AATG(){
     return dinucSyn_stat[2][0][14];
  }
  double GetDinucSyn31_AATT(){
     return dinucSyn_stat[2][0][15];
  }
  double GetDinucSyn31_ACAA(){
     return dinucSyn_stat[2][1][0];
  }
  double GetDinucSyn31_ACAC(){
     return dinucSyn_stat[2][1][1];
  }
  double GetDinucSyn31_ACAG(){
     return dinucSyn_stat[2][1][2];
  }
  double GetDinucSyn31_ACAT(){
     return dinucSyn_stat[2][1][3];
  }
  double GetDinucSyn31_ACCA(){
     return dinucSyn_stat[2][1][4];
  }
  double GetDinucSyn31_ACCC(){
     return dinucSyn_stat[2][1][5];
  }
  double GetDinucSyn31_ACCG(){
     return dinucSyn_stat[2][1][6];
  }
  double GetDinucSyn31_ACCT(){
     return dinucSyn_stat[2][1][7];
  }
  double GetDinucSyn31_ACGA(){
     return dinucSyn_stat[2][1][8];
  }
  double GetDinucSyn31_ACGC(){
     return dinucSyn_stat[2][1][9];
  }
  double GetDinucSyn31_ACGG(){
     return dinucSyn_stat[2][1][10];
  }
  double GetDinucSyn31_ACGT(){
     return dinucSyn_stat[2][1][11];
  }
  double GetDinucSyn31_ACTA(){
     return dinucSyn_stat[2][1][12];
  }
  double GetDinucSyn31_ACTC(){
     return dinucSyn_stat[2][1][13];
  }
  double GetDinucSyn31_ACTG(){
     return dinucSyn_stat[2][1][14];
  }
  double GetDinucSyn31_ACTT(){
     return dinucSyn_stat[2][1][15];
  }
  double GetDinucSyn31_AGAA(){
     return dinucSyn_stat[2][2][0];
  }
  double GetDinucSyn31_AGAC(){
     return dinucSyn_stat[2][2][1];
  }
  double GetDinucSyn31_AGAG(){
     return dinucSyn_stat[2][2][2];
  }
  double GetDinucSyn31_AGAT(){
     return dinucSyn_stat[2][2][3];
  }
  double GetDinucSyn31_AGCA(){
     return dinucSyn_stat[2][2][4];
  }
  double GetDinucSyn31_AGCC(){
     return dinucSyn_stat[2][2][5];
  }
  double GetDinucSyn31_AGCG(){
     return dinucSyn_stat[2][2][6];
  }
  double GetDinucSyn31_AGCT(){
     return dinucSyn_stat[2][2][7];
  }
  double GetDinucSyn31_AGGA(){
     return dinucSyn_stat[2][2][8];
  }
  double GetDinucSyn31_AGGC(){
     return dinucSyn_stat[2][2][9];
  }
  double GetDinucSyn31_AGGG(){
     return dinucSyn_stat[2][2][10];
  }
  double GetDinucSyn31_AGGT(){
     return dinucSyn_stat[2][2][11];
  }
  double GetDinucSyn31_AGTA(){
     return dinucSyn_stat[2][2][12];
  }
  double GetDinucSyn31_AGTC(){
     return dinucSyn_stat[2][2][13];
  }
  double GetDinucSyn31_AGTG(){
     return dinucSyn_stat[2][2][14];
  }
  double GetDinucSyn31_AGTT(){
     return dinucSyn_stat[2][2][15];
  }
  double GetDinucSyn31_ATAA(){
     return dinucSyn_stat[2][3][0];
  }
  double GetDinucSyn31_ATAC(){
     return dinucSyn_stat[2][3][1];
  }
  double GetDinucSyn31_ATAG(){
     return dinucSyn_stat[2][3][2];
  }
  double GetDinucSyn31_ATAT(){
     return dinucSyn_stat[2][3][3];
  }
  double GetDinucSyn31_ATCA(){
     return dinucSyn_stat[2][3][4];
  }
  double GetDinucSyn31_ATCC(){
     return dinucSyn_stat[2][3][5];
  }
  double GetDinucSyn31_ATCG(){
     return dinucSyn_stat[2][3][6];
  }
  double GetDinucSyn31_ATCT(){
     return dinucSyn_stat[2][3][7];
  }
  double GetDinucSyn31_ATGA(){
     return dinucSyn_stat[2][3][8];
  }
  double GetDinucSyn31_ATGC(){
     return dinucSyn_stat[2][3][9];
  }
  double GetDinucSyn31_ATGG(){
     return dinucSyn_stat[2][3][10];
  }
  double GetDinucSyn31_ATGT(){
     return dinucSyn_stat[2][3][11];
  }
  double GetDinucSyn31_ATTA(){
     return dinucSyn_stat[2][3][12];
  }
  double GetDinucSyn31_ATTC(){
     return dinucSyn_stat[2][3][13];
  }
  double GetDinucSyn31_ATTG(){
     return dinucSyn_stat[2][3][14];
  }
  double GetDinucSyn31_ATTT(){
     return dinucSyn_stat[2][3][15];
  }
  double GetDinucSyn31_CAAA(){
     return dinucSyn_stat[2][4][0];
  }
  double GetDinucSyn31_CAAC(){
     return dinucSyn_stat[2][4][1];
  }
  double GetDinucSyn31_CAAG(){
     return dinucSyn_stat[2][4][2];
  }
  double GetDinucSyn31_CAAT(){
     return dinucSyn_stat[2][4][3];
  }
  double GetDinucSyn31_CACA(){
     return dinucSyn_stat[2][4][4];
  }
  double GetDinucSyn31_CACC(){
     return dinucSyn_stat[2][4][5];
  }
  double GetDinucSyn31_CACG(){
     return dinucSyn_stat[2][4][6];
  }
  double GetDinucSyn31_CACT(){
     return dinucSyn_stat[2][4][7];
  }
  double GetDinucSyn31_CAGA(){
     return dinucSyn_stat[2][4][8];
  }
  double GetDinucSyn31_CAGC(){
     return dinucSyn_stat[2][4][9];
  }
  double GetDinucSyn31_CAGG(){
     return dinucSyn_stat[2][4][10];
  }
  double GetDinucSyn31_CAGT(){
     return dinucSyn_stat[2][4][11];
  }
  double GetDinucSyn31_CATA(){
     return dinucSyn_stat[2][4][12];
  }
  double GetDinucSyn31_CATC(){
     return dinucSyn_stat[2][4][13];
  }
  double GetDinucSyn31_CATG(){
     return dinucSyn_stat[2][4][14];
  }
  double GetDinucSyn31_CATT(){
     return dinucSyn_stat[2][4][15];
  }
  double GetDinucSyn31_CCAA(){
     return dinucSyn_stat[2][5][0];
  }
  double GetDinucSyn31_CCAC(){
     return dinucSyn_stat[2][5][1];
  }
  double GetDinucSyn31_CCAG(){
     return dinucSyn_stat[2][5][2];
  }
  double GetDinucSyn31_CCAT(){
     return dinucSyn_stat[2][5][3];
  }
  double GetDinucSyn31_CCCA(){
     return dinucSyn_stat[2][5][4];
  }
  double GetDinucSyn31_CCCC(){
     return dinucSyn_stat[2][5][5];
  }
  double GetDinucSyn31_CCCG(){
     return dinucSyn_stat[2][5][6];
  }
  double GetDinucSyn31_CCCT(){
     return dinucSyn_stat[2][5][7];
  }
  double GetDinucSyn31_CCGA(){
     return dinucSyn_stat[2][5][8];
  }
  double GetDinucSyn31_CCGC(){
     return dinucSyn_stat[2][5][9];
  }
  double GetDinucSyn31_CCGG(){
     return dinucSyn_stat[2][5][10];
  }
  double GetDinucSyn31_CCGT(){
     return dinucSyn_stat[2][5][11];
  }
  double GetDinucSyn31_CCTA(){
     return dinucSyn_stat[2][5][12];
  }
  double GetDinucSyn31_CCTC(){
     return dinucSyn_stat[2][5][13];
  }
  double GetDinucSyn31_CCTG(){
     return dinucSyn_stat[2][5][14];
  }
  double GetDinucSyn31_CCTT(){
     return dinucSyn_stat[2][5][15];
  }
  double GetDinucSyn31_CGAA(){
     return dinucSyn_stat[2][6][0];
  }
  double GetDinucSyn31_CGAC(){
     return dinucSyn_stat[2][6][1];
  }
  double GetDinucSyn31_CGAG(){
     return dinucSyn_stat[2][6][2];
  }
  double GetDinucSyn31_CGAT(){
     return dinucSyn_stat[2][6][3];
  }
  double GetDinucSyn31_CGCA(){
     return dinucSyn_stat[2][6][4];
  }
  double GetDinucSyn31_CGCC(){
     return dinucSyn_stat[2][6][5];
  }
  double GetDinucSyn31_CGCG(){
     return dinucSyn_stat[2][6][6];
  }
  double GetDinucSyn31_CGCT(){
     return dinucSyn_stat[2][6][7];
  }
  double GetDinucSyn31_CGGA(){
     return dinucSyn_stat[2][6][8];
  }
  double GetDinucSyn31_CGGC(){
     return dinucSyn_stat[2][6][9];
  }
  double GetDinucSyn31_CGGG(){
     return dinucSyn_stat[2][6][10];
  }
  double GetDinucSyn31_CGGT(){
     return dinucSyn_stat[2][6][11];
  }
  double GetDinucSyn31_CGTA(){
     return dinucSyn_stat[2][6][12];
  }
  double GetDinucSyn31_CGTC(){
     return dinucSyn_stat[2][6][13];
  }
  double GetDinucSyn31_CGTG(){
     return dinucSyn_stat[2][6][14];
  }
  double GetDinucSyn31_CGTT(){
     return dinucSyn_stat[2][6][15];
  }
  double GetDinucSyn31_CTAA(){
     return dinucSyn_stat[2][7][0];
  }
  double GetDinucSyn31_CTAC(){
     return dinucSyn_stat[2][7][1];
  }
  double GetDinucSyn31_CTAG(){
     return dinucSyn_stat[2][7][2];
  }
  double GetDinucSyn31_CTAT(){
     return dinucSyn_stat[2][7][3];
  }
  double GetDinucSyn31_CTCA(){
     return dinucSyn_stat[2][7][4];
  }
  double GetDinucSyn31_CTCC(){
     return dinucSyn_stat[2][7][5];
  }
  double GetDinucSyn31_CTCG(){
     return dinucSyn_stat[2][7][6];
  }
  double GetDinucSyn31_CTCT(){
     return dinucSyn_stat[2][7][7];
  }
  double GetDinucSyn31_CTGA(){
     return dinucSyn_stat[2][7][8];
  }
  double GetDinucSyn31_CTGC(){
     return dinucSyn_stat[2][7][9];
  }
  double GetDinucSyn31_CTGG(){
     return dinucSyn_stat[2][7][10];
  }
  double GetDinucSyn31_CTGT(){
     return dinucSyn_stat[2][7][11];
  }
  double GetDinucSyn31_CTTA(){
     return dinucSyn_stat[2][7][12];
  }
  double GetDinucSyn31_CTTC(){
     return dinucSyn_stat[2][7][13];
  }
  double GetDinucSyn31_CTTG(){
     return dinucSyn_stat[2][7][14];
  }
  double GetDinucSyn31_CTTT(){
     return dinucSyn_stat[2][7][15];
  }
  double GetDinucSyn31_GAAA(){
     return dinucSyn_stat[2][8][0];
  }
  double GetDinucSyn31_GAAC(){
     return dinucSyn_stat[2][8][1];
  }
  double GetDinucSyn31_GAAG(){
     return dinucSyn_stat[2][8][2];
  }
  double GetDinucSyn31_GAAT(){
     return dinucSyn_stat[2][8][3];
  }
  double GetDinucSyn31_GACA(){
     return dinucSyn_stat[2][8][4];
  }
  double GetDinucSyn31_GACC(){
     return dinucSyn_stat[2][8][5];
  }
  double GetDinucSyn31_GACG(){
     return dinucSyn_stat[2][8][6];
  }
  double GetDinucSyn31_GACT(){
     return dinucSyn_stat[2][8][7];
  }
  double GetDinucSyn31_GAGA(){
     return dinucSyn_stat[2][8][8];
  }
  double GetDinucSyn31_GAGC(){
     return dinucSyn_stat[2][8][9];
  }
  double GetDinucSyn31_GAGG(){
     return dinucSyn_stat[2][8][10];
  }
  double GetDinucSyn31_GAGT(){
     return dinucSyn_stat[2][8][11];
  }
  double GetDinucSyn31_GATA(){
     return dinucSyn_stat[2][8][12];
  }
  double GetDinucSyn31_GATC(){
     return dinucSyn_stat[2][8][13];
  }
  double GetDinucSyn31_GATG(){
     return dinucSyn_stat[2][8][14];
  }
  double GetDinucSyn31_GATT(){
     return dinucSyn_stat[2][8][15];
  }
  double GetDinucSyn31_GCAA(){
     return dinucSyn_stat[2][9][0];
  }
  double GetDinucSyn31_GCAC(){
     return dinucSyn_stat[2][9][1];
  }
  double GetDinucSyn31_GCAG(){
     return dinucSyn_stat[2][9][2];
  }
  double GetDinucSyn31_GCAT(){
     return dinucSyn_stat[2][9][3];
  }
  double GetDinucSyn31_GCCA(){
     return dinucSyn_stat[2][9][4];
  }
  double GetDinucSyn31_GCCC(){
     return dinucSyn_stat[2][9][5];
  }
  double GetDinucSyn31_GCCG(){
     return dinucSyn_stat[2][9][6];
  }
  double GetDinucSyn31_GCCT(){
     return dinucSyn_stat[2][9][7];
  }
  double GetDinucSyn31_GCGA(){
     return dinucSyn_stat[2][9][8];
  }
  double GetDinucSyn31_GCGC(){
     return dinucSyn_stat[2][9][9];
  }
  double GetDinucSyn31_GCGG(){
     return dinucSyn_stat[2][9][10];
  }
  double GetDinucSyn31_GCGT(){
     return dinucSyn_stat[2][9][11];
  }
  double GetDinucSyn31_GCTA(){
     return dinucSyn_stat[2][9][12];
  }
  double GetDinucSyn31_GCTC(){
     return dinucSyn_stat[2][9][13];
  }
  double GetDinucSyn31_GCTG(){
     return dinucSyn_stat[2][9][14];
  }
  double GetDinucSyn31_GCTT(){
     return dinucSyn_stat[2][9][15];
  }
  double GetDinucSyn31_GGAA(){
     return dinucSyn_stat[2][10][0];
  }
  double GetDinucSyn31_GGAC(){
     return dinucSyn_stat[2][10][1];
  }
  double GetDinucSyn31_GGAG(){
     return dinucSyn_stat[2][10][2];
  }
  double GetDinucSyn31_GGAT(){
     return dinucSyn_stat[2][10][3];
  }
  double GetDinucSyn31_GGCA(){
     return dinucSyn_stat[2][10][4];
  }
  double GetDinucSyn31_GGCC(){
     return dinucSyn_stat[2][10][5];
  }
  double GetDinucSyn31_GGCG(){
     return dinucSyn_stat[2][10][6];
  }
  double GetDinucSyn31_GGCT(){
     return dinucSyn_stat[2][10][7];
  }
  double GetDinucSyn31_GGGA(){
     return dinucSyn_stat[2][10][8];
  }
  double GetDinucSyn31_GGGC(){
     return dinucSyn_stat[2][10][9];
  }
  double GetDinucSyn31_GGGG(){
     return dinucSyn_stat[2][10][10];
  }
  double GetDinucSyn31_GGGT(){
     return dinucSyn_stat[2][10][11];
  }
  double GetDinucSyn31_GGTA(){
     return dinucSyn_stat[2][10][12];
  }
  double GetDinucSyn31_GGTC(){
     return dinucSyn_stat[2][10][13];
  }
  double GetDinucSyn31_GGTG(){
     return dinucSyn_stat[2][10][14];
  }
  double GetDinucSyn31_GGTT(){
     return dinucSyn_stat[2][10][15];
  }
  double GetDinucSyn31_GTAA(){
     return dinucSyn_stat[2][11][0];
  }
  double GetDinucSyn31_GTAC(){
     return dinucSyn_stat[2][11][1];
  }
  double GetDinucSyn31_GTAG(){
     return dinucSyn_stat[2][11][2];
  }
  double GetDinucSyn31_GTAT(){
     return dinucSyn_stat[2][11][3];
  }
  double GetDinucSyn31_GTCA(){
     return dinucSyn_stat[2][11][4];
  }
  double GetDinucSyn31_GTCC(){
     return dinucSyn_stat[2][11][5];
  }
  double GetDinucSyn31_GTCG(){
     return dinucSyn_stat[2][11][6];
  }
  double GetDinucSyn31_GTCT(){
     return dinucSyn_stat[2][11][7];
  }
  double GetDinucSyn31_GTGA(){
     return dinucSyn_stat[2][11][8];
  }
  double GetDinucSyn31_GTGC(){
     return dinucSyn_stat[2][11][9];
  }
  double GetDinucSyn31_GTGG(){
     return dinucSyn_stat[2][11][10];
  }
  double GetDinucSyn31_GTGT(){
     return dinucSyn_stat[2][11][11];
  }
  double GetDinucSyn31_GTTA(){
     return dinucSyn_stat[2][11][12];
  }
  double GetDinucSyn31_GTTC(){
     return dinucSyn_stat[2][11][13];
  }
  double GetDinucSyn31_GTTG(){
     return dinucSyn_stat[2][11][14];
  }
  double GetDinucSyn31_GTTT(){
     return dinucSyn_stat[2][11][15];
  }
  double GetDinucSyn31_TAAA(){
     return dinucSyn_stat[2][12][0];
  }
  double GetDinucSyn31_TAAC(){
     return dinucSyn_stat[2][12][1];
  }
  double GetDinucSyn31_TAAG(){
     return dinucSyn_stat[2][12][2];
  }
  double GetDinucSyn31_TAAT(){
     return dinucSyn_stat[2][12][3];
  }
  double GetDinucSyn31_TACA(){
     return dinucSyn_stat[2][12][4];
  }
  double GetDinucSyn31_TACC(){
     return dinucSyn_stat[2][12][5];
  }
  double GetDinucSyn31_TACG(){
     return dinucSyn_stat[2][12][6];
  }
  double GetDinucSyn31_TACT(){
     return dinucSyn_stat[2][12][7];
  }
  double GetDinucSyn31_TAGA(){
     return dinucSyn_stat[2][12][8];
  }
  double GetDinucSyn31_TAGC(){
     return dinucSyn_stat[2][12][9];
  }
  double GetDinucSyn31_TAGG(){
     return dinucSyn_stat[2][12][10];
  }
  double GetDinucSyn31_TAGT(){
     return dinucSyn_stat[2][12][11];
  }
  double GetDinucSyn31_TATA(){
     return dinucSyn_stat[2][12][12];
  }
  double GetDinucSyn31_TATC(){
     return dinucSyn_stat[2][12][13];
  }
  double GetDinucSyn31_TATG(){
     return dinucSyn_stat[2][12][14];
  }
  double GetDinucSyn31_TATT(){
     return dinucSyn_stat[2][12][15];
  }
  double GetDinucSyn31_TCAA(){
     return dinucSyn_stat[2][13][0];
  }
  double GetDinucSyn31_TCAC(){
     return dinucSyn_stat[2][13][1];
  }
  double GetDinucSyn31_TCAG(){
     return dinucSyn_stat[2][13][2];
  }
  double GetDinucSyn31_TCAT(){
     return dinucSyn_stat[2][13][3];
  }
  double GetDinucSyn31_TCCA(){
     return dinucSyn_stat[2][13][4];
  }
  double GetDinucSyn31_TCCC(){
     return dinucSyn_stat[2][13][5];
  }
  double GetDinucSyn31_TCCG(){
     return dinucSyn_stat[2][13][6];
  }
  double GetDinucSyn31_TCCT(){
     return dinucSyn_stat[2][13][7];
  }
  double GetDinucSyn31_TCGA(){
     return dinucSyn_stat[2][13][8];
  }
  double GetDinucSyn31_TCGC(){
     return dinucSyn_stat[2][13][9];
  }
  double GetDinucSyn31_TCGG(){
     return dinucSyn_stat[2][13][10];
  }
  double GetDinucSyn31_TCGT(){
     return dinucSyn_stat[2][13][11];
  }
  double GetDinucSyn31_TCTA(){
     return dinucSyn_stat[2][13][12];
  }
  double GetDinucSyn31_TCTC(){
     return dinucSyn_stat[2][13][13];
  }
  double GetDinucSyn31_TCTG(){
     return dinucSyn_stat[2][13][14];
  }
  double GetDinucSyn31_TCTT(){
     return dinucSyn_stat[2][13][15];
  }
  double GetDinucSyn31_TGAA(){
     return dinucSyn_stat[2][14][0];
  }
  double GetDinucSyn31_TGAC(){
     return dinucSyn_stat[2][14][1];
  }
  double GetDinucSyn31_TGAG(){
     return dinucSyn_stat[2][14][2];
  }
  double GetDinucSyn31_TGAT(){
     return dinucSyn_stat[2][14][3];
  }
  double GetDinucSyn31_TGCA(){
     return dinucSyn_stat[2][14][4];
  }
  double GetDinucSyn31_TGCC(){
     return dinucSyn_stat[2][14][5];
  }
  double GetDinucSyn31_TGCG(){
     return dinucSyn_stat[2][14][6];
  }
  double GetDinucSyn31_TGCT(){
     return dinucSyn_stat[2][14][7];
  }
  double GetDinucSyn31_TGGA(){
     return dinucSyn_stat[2][14][8];
  }
  double GetDinucSyn31_TGGC(){
     return dinucSyn_stat[2][14][9];
  }
  double GetDinucSyn31_TGGG(){
     return dinucSyn_stat[2][14][10];
  }
  double GetDinucSyn31_TGGT(){
     return dinucSyn_stat[2][14][11];
  }
  double GetDinucSyn31_TGTA(){
     return dinucSyn_stat[2][14][12];
  }
  double GetDinucSyn31_TGTC(){
     return dinucSyn_stat[2][14][13];
  }
  double GetDinucSyn31_TGTG(){
     return dinucSyn_stat[2][14][14];
  }
  double GetDinucSyn31_TGTT(){
     return dinucSyn_stat[2][14][15];
  }
  double GetDinucSyn31_TTAA(){
     return dinucSyn_stat[2][15][0];
  }
  double GetDinucSyn31_TTAC(){
     return dinucSyn_stat[2][15][1];
  }
  double GetDinucSyn31_TTAG(){
     return dinucSyn_stat[2][15][2];
  }
  double GetDinucSyn31_TTAT(){
     return dinucSyn_stat[2][15][3];
  }
  double GetDinucSyn31_TTCA(){
     return dinucSyn_stat[2][15][4];
  }
  double GetDinucSyn31_TTCC(){
     return dinucSyn_stat[2][15][5];
  }
  double GetDinucSyn31_TTCG(){
     return dinucSyn_stat[2][15][6];
  }
  double GetDinucSyn31_TTCT(){
     return dinucSyn_stat[2][15][7];
  }
  double GetDinucSyn31_TTGA(){
     return dinucSyn_stat[2][15][8];
  }
  double GetDinucSyn31_TTGC(){
     return dinucSyn_stat[2][15][9];
  }
  double GetDinucSyn31_TTGG(){
     return dinucSyn_stat[2][15][10];
  }
  double GetDinucSyn31_TTGT(){
     return dinucSyn_stat[2][15][11];
  }
  double GetDinucSyn31_TTTA(){
     return dinucSyn_stat[2][15][12];
  }
  double GetDinucSyn31_TTTC(){
     return dinucSyn_stat[2][15][13];
  }
  double GetDinucSyn31_TTTG(){
     return dinucSyn_stat[2][15][14];
  }
  double GetDinucSyn31_TTTT(){
     return dinucSyn_stat[2][15][15];
  }

  double GetDinucNSyn_AAAA(){
     return dinucNSyn_stat[3][0][0];
 }
 double GetDinucNSyn_AAAC(){
     return dinucNSyn_stat[3][0][1];
 }
 double GetDinucNSyn_AAAG(){
     return dinucNSyn_stat[3][0][2];
 }
 double GetDinucNSyn_AAAT(){
     return dinucNSyn_stat[3][0][3];
 }
 double GetDinucNSyn_AACA(){
     return dinucNSyn_stat[3][0][4];
 }
 double GetDinucNSyn_AACC(){
     return dinucNSyn_stat[3][0][5];
 }
 double GetDinucNSyn_AACG(){
     return dinucNSyn_stat[3][0][6];
 }
 double GetDinucNSyn_AACT(){
     return dinucNSyn_stat[3][0][7];
 }
 double GetDinucNSyn_AAGA(){
     return dinucNSyn_stat[3][0][8];
 }
 double GetDinucNSyn_AAGC(){
     return dinucNSyn_stat[3][0][9];
 }
 double GetDinucNSyn_AAGG(){
     return dinucNSyn_stat[3][0][10];
 }
 double GetDinucNSyn_AAGT(){
     return dinucNSyn_stat[3][0][11];
 }
 double GetDinucNSyn_AATA(){
     return dinucNSyn_stat[3][0][12];
 }
 double GetDinucNSyn_AATC(){
     return dinucNSyn_stat[3][0][13];
 }
 double GetDinucNSyn_AATG(){
     return dinucNSyn_stat[3][0][14];
 }
 double GetDinucNSyn_AATT(){
     return dinucNSyn_stat[3][0][15];
 }
 double GetDinucNSyn_ACAA(){
     return dinucNSyn_stat[3][1][0];
 }
 double GetDinucNSyn_ACAC(){
     return dinucNSyn_stat[3][1][1];
 }
 double GetDinucNSyn_ACAG(){
     return dinucNSyn_stat[3][1][2];
 }
 double GetDinucNSyn_ACAT(){
     return dinucNSyn_stat[3][1][3];
 }
 double GetDinucNSyn_ACCA(){
     return dinucNSyn_stat[3][1][4];
 }
 double GetDinucNSyn_ACCC(){
     return dinucNSyn_stat[3][1][5];
 }
 double GetDinucNSyn_ACCG(){
     return dinucNSyn_stat[3][1][6];
 }
 double GetDinucNSyn_ACCT(){
     return dinucNSyn_stat[3][1][7];
 }
 double GetDinucNSyn_ACGA(){
     return dinucNSyn_stat[3][1][8];
 }
 double GetDinucNSyn_ACGC(){
     return dinucNSyn_stat[3][1][9];
 }
 double GetDinucNSyn_ACGG(){
     return dinucNSyn_stat[3][1][10];
 }
 double GetDinucNSyn_ACGT(){
     return dinucNSyn_stat[3][1][11];
 }
 double GetDinucNSyn_ACTA(){
     return dinucNSyn_stat[3][1][12];
 }
 double GetDinucNSyn_ACTC(){
     return dinucNSyn_stat[3][1][13];
 }
 double GetDinucNSyn_ACTG(){
     return dinucNSyn_stat[3][1][14];
 }
 double GetDinucNSyn_ACTT(){
     return dinucNSyn_stat[3][1][15];
 }
 double GetDinucNSyn_AGAA(){
     return dinucNSyn_stat[3][2][0];
 }
 double GetDinucNSyn_AGAC(){
     return dinucNSyn_stat[3][2][1];
 }
 double GetDinucNSyn_AGAG(){
     return dinucNSyn_stat[3][2][2];
 }
 double GetDinucNSyn_AGAT(){
     return dinucNSyn_stat[3][2][3];
 }
 double GetDinucNSyn_AGCA(){
     return dinucNSyn_stat[3][2][4];
 }
 double GetDinucNSyn_AGCC(){
     return dinucNSyn_stat[3][2][5];
 }
 double GetDinucNSyn_AGCG(){
     return dinucNSyn_stat[3][2][6];
 }
 double GetDinucNSyn_AGCT(){
     return dinucNSyn_stat[3][2][7];
 }
 double GetDinucNSyn_AGGA(){
     return dinucNSyn_stat[3][2][8];
 }
 double GetDinucNSyn_AGGC(){
     return dinucNSyn_stat[3][2][9];
 }
 double GetDinucNSyn_AGGG(){
     return dinucNSyn_stat[3][2][10];
 }
 double GetDinucNSyn_AGGT(){
     return dinucNSyn_stat[3][2][11];
 }
 double GetDinucNSyn_AGTA(){
     return dinucNSyn_stat[3][2][12];
 }
 double GetDinucNSyn_AGTC(){
     return dinucNSyn_stat[3][2][13];
 }
 double GetDinucNSyn_AGTG(){
     return dinucNSyn_stat[3][2][14];
 }
 double GetDinucNSyn_AGTT(){
     return dinucNSyn_stat[3][2][15];
 }
 double GetDinucNSyn_ATAA(){
     return dinucNSyn_stat[3][3][0];
 }
 double GetDinucNSyn_ATAC(){
     return dinucNSyn_stat[3][3][1];
 }
 double GetDinucNSyn_ATAG(){
     return dinucNSyn_stat[3][3][2];
 }
 double GetDinucNSyn_ATAT(){
     return dinucNSyn_stat[3][3][3];
 }
 double GetDinucNSyn_ATCA(){
     return dinucNSyn_stat[3][3][4];
 }
 double GetDinucNSyn_ATCC(){
     return dinucNSyn_stat[3][3][5];
 }
 double GetDinucNSyn_ATCG(){
     return dinucNSyn_stat[3][3][6];
 }
 double GetDinucNSyn_ATCT(){
     return dinucNSyn_stat[3][3][7];
 }
 double GetDinucNSyn_ATGA(){
     return dinucNSyn_stat[3][3][8];
 }
 double GetDinucNSyn_ATGC(){
     return dinucNSyn_stat[3][3][9];
 }
 double GetDinucNSyn_ATGG(){
     return dinucNSyn_stat[3][3][10];
 }
 double GetDinucNSyn_ATGT(){
     return dinucNSyn_stat[3][3][11];
 }
 double GetDinucNSyn_ATTA(){
     return dinucNSyn_stat[3][3][12];
 }
 double GetDinucNSyn_ATTC(){
     return dinucNSyn_stat[3][3][13];
 }
 double GetDinucNSyn_ATTG(){
     return dinucNSyn_stat[3][3][14];
 }
 double GetDinucNSyn_ATTT(){
     return dinucNSyn_stat[3][3][15];
 }
 double GetDinucNSyn_CAAA(){
     return dinucNSyn_stat[3][4][0];
 }
 double GetDinucNSyn_CAAC(){
     return dinucNSyn_stat[3][4][1];
 }
 double GetDinucNSyn_CAAG(){
     return dinucNSyn_stat[3][4][2];
 }
 double GetDinucNSyn_CAAT(){
     return dinucNSyn_stat[3][4][3];
 }
 double GetDinucNSyn_CACA(){
     return dinucNSyn_stat[3][4][4];
 }
 double GetDinucNSyn_CACC(){
     return dinucNSyn_stat[3][4][5];
 }
 double GetDinucNSyn_CACG(){
     return dinucNSyn_stat[3][4][6];
 }
 double GetDinucNSyn_CACT(){
     return dinucNSyn_stat[3][4][7];
 }
 double GetDinucNSyn_CAGA(){
     return dinucNSyn_stat[3][4][8];
 }
 double GetDinucNSyn_CAGC(){
     return dinucNSyn_stat[3][4][9];
 }
 double GetDinucNSyn_CAGG(){
     return dinucNSyn_stat[3][4][10];
 }
 double GetDinucNSyn_CAGT(){
     return dinucNSyn_stat[3][4][11];
 }
 double GetDinucNSyn_CATA(){
     return dinucNSyn_stat[3][4][12];
 }
 double GetDinucNSyn_CATC(){
     return dinucNSyn_stat[3][4][13];
 }
 double GetDinucNSyn_CATG(){
     return dinucNSyn_stat[3][4][14];
 }
 double GetDinucNSyn_CATT(){
     return dinucNSyn_stat[3][4][15];
 }
 double GetDinucNSyn_CCAA(){
     return dinucNSyn_stat[3][5][0];
 }
 double GetDinucNSyn_CCAC(){
     return dinucNSyn_stat[3][5][1];
 }
 double GetDinucNSyn_CCAG(){
     return dinucNSyn_stat[3][5][2];
 }
 double GetDinucNSyn_CCAT(){
     return dinucNSyn_stat[3][5][3];
 }
 double GetDinucNSyn_CCCA(){
     return dinucNSyn_stat[3][5][4];
 }
 double GetDinucNSyn_CCCC(){
     return dinucNSyn_stat[3][5][5];
 }
 double GetDinucNSyn_CCCG(){
     return dinucNSyn_stat[3][5][6];
 }
 double GetDinucNSyn_CCCT(){
     return dinucNSyn_stat[3][5][7];
 }
 double GetDinucNSyn_CCGA(){
     return dinucNSyn_stat[3][5][8];
 }
 double GetDinucNSyn_CCGC(){
     return dinucNSyn_stat[3][5][9];
 }
 double GetDinucNSyn_CCGG(){
     return dinucNSyn_stat[3][5][10];
 }
 double GetDinucNSyn_CCGT(){
     return dinucNSyn_stat[3][5][11];
 }
 double GetDinucNSyn_CCTA(){
     return dinucNSyn_stat[3][5][12];
 }
 double GetDinucNSyn_CCTC(){
     return dinucNSyn_stat[3][5][13];
 }
 double GetDinucNSyn_CCTG(){
     return dinucNSyn_stat[3][5][14];
 }
 double GetDinucNSyn_CCTT(){
     return dinucNSyn_stat[3][5][15];
 }
 double GetDinucNSyn_CGAA(){
     return dinucNSyn_stat[3][6][0];
 }
 double GetDinucNSyn_CGAC(){
     return dinucNSyn_stat[3][6][1];
 }
 double GetDinucNSyn_CGAG(){
     return dinucNSyn_stat[3][6][2];
 }
 double GetDinucNSyn_CGAT(){
     return dinucNSyn_stat[3][6][3];
 }
 double GetDinucNSyn_CGCA(){
     return dinucNSyn_stat[3][6][4];
 }
 double GetDinucNSyn_CGCC(){
     return dinucNSyn_stat[3][6][5];
 }
 double GetDinucNSyn_CGCG(){
     return dinucNSyn_stat[3][6][6];
 }
 double GetDinucNSyn_CGCT(){
     return dinucNSyn_stat[3][6][7];
 }
 double GetDinucNSyn_CGGA(){
     return dinucNSyn_stat[3][6][8];
 }
 double GetDinucNSyn_CGGC(){
     return dinucNSyn_stat[3][6][9];
 }
 double GetDinucNSyn_CGGG(){
     return dinucNSyn_stat[3][6][10];
 }
 double GetDinucNSyn_CGGT(){
     return dinucNSyn_stat[3][6][11];
 }
 double GetDinucNSyn_CGTA(){
     return dinucNSyn_stat[3][6][12];
 }
 double GetDinucNSyn_CGTC(){
     return dinucNSyn_stat[3][6][13];
 }
 double GetDinucNSyn_CGTG(){
     return dinucNSyn_stat[3][6][14];
 }
 double GetDinucNSyn_CGTT(){
     return dinucNSyn_stat[3][6][15];
 }
 double GetDinucNSyn_CTAA(){
     return dinucNSyn_stat[3][7][0];
 }
 double GetDinucNSyn_CTAC(){
     return dinucNSyn_stat[3][7][1];
 }
 double GetDinucNSyn_CTAG(){
     return dinucNSyn_stat[3][7][2];
 }
 double GetDinucNSyn_CTAT(){
     return dinucNSyn_stat[3][7][3];
 }
 double GetDinucNSyn_CTCA(){
     return dinucNSyn_stat[3][7][4];
 }
 double GetDinucNSyn_CTCC(){
     return dinucNSyn_stat[3][7][5];
 }
 double GetDinucNSyn_CTCG(){
     return dinucNSyn_stat[3][7][6];
 }
 double GetDinucNSyn_CTCT(){
     return dinucNSyn_stat[3][7][7];
 }
 double GetDinucNSyn_CTGA(){
     return dinucNSyn_stat[3][7][8];
 }
 double GetDinucNSyn_CTGC(){
     return dinucNSyn_stat[3][7][9];
 }
 double GetDinucNSyn_CTGG(){
     return dinucNSyn_stat[3][7][10];
 }
 double GetDinucNSyn_CTGT(){
     return dinucNSyn_stat[3][7][11];
 }
 double GetDinucNSyn_CTTA(){
     return dinucNSyn_stat[3][7][12];
 }
 double GetDinucNSyn_CTTC(){
     return dinucNSyn_stat[3][7][13];
 }
 double GetDinucNSyn_CTTG(){
     return dinucNSyn_stat[3][7][14];
 }
 double GetDinucNSyn_CTTT(){
     return dinucNSyn_stat[3][7][15];
 }
 double GetDinucNSyn_GAAA(){
     return dinucNSyn_stat[3][8][0];
 }
 double GetDinucNSyn_GAAC(){
     return dinucNSyn_stat[3][8][1];
 }
 double GetDinucNSyn_GAAG(){
     return dinucNSyn_stat[3][8][2];
 }
 double GetDinucNSyn_GAAT(){
     return dinucNSyn_stat[3][8][3];
 }
 double GetDinucNSyn_GACA(){
     return dinucNSyn_stat[3][8][4];
 }
 double GetDinucNSyn_GACC(){
     return dinucNSyn_stat[3][8][5];
 }
 double GetDinucNSyn_GACG(){
     return dinucNSyn_stat[3][8][6];
 }
 double GetDinucNSyn_GACT(){
     return dinucNSyn_stat[3][8][7];
 }
 double GetDinucNSyn_GAGA(){
     return dinucNSyn_stat[3][8][8];
 }
 double GetDinucNSyn_GAGC(){
     return dinucNSyn_stat[3][8][9];
 }
 double GetDinucNSyn_GAGG(){
     return dinucNSyn_stat[3][8][10];
 }
 double GetDinucNSyn_GAGT(){
     return dinucNSyn_stat[3][8][11];
 }
 double GetDinucNSyn_GATA(){
     return dinucNSyn_stat[3][8][12];
 }
 double GetDinucNSyn_GATC(){
     return dinucNSyn_stat[3][8][13];
 }
 double GetDinucNSyn_GATG(){
     return dinucNSyn_stat[3][8][14];
 }
 double GetDinucNSyn_GATT(){
     return dinucNSyn_stat[3][8][15];
 }
 double GetDinucNSyn_GCAA(){
     return dinucNSyn_stat[3][9][0];
 }
 double GetDinucNSyn_GCAC(){
     return dinucNSyn_stat[3][9][1];
 }
 double GetDinucNSyn_GCAG(){
     return dinucNSyn_stat[3][9][2];
 }
 double GetDinucNSyn_GCAT(){
     return dinucNSyn_stat[3][9][3];
 }
 double GetDinucNSyn_GCCA(){
     return dinucNSyn_stat[3][9][4];
 }
 double GetDinucNSyn_GCCC(){
     return dinucNSyn_stat[3][9][5];
 }
 double GetDinucNSyn_GCCG(){
     return dinucNSyn_stat[3][9][6];
 }
 double GetDinucNSyn_GCCT(){
     return dinucNSyn_stat[3][9][7];
 }
 double GetDinucNSyn_GCGA(){
     return dinucNSyn_stat[3][9][8];
 }
 double GetDinucNSyn_GCGC(){
     return dinucNSyn_stat[3][9][9];
 }
 double GetDinucNSyn_GCGG(){
     return dinucNSyn_stat[3][9][10];
 }
 double GetDinucNSyn_GCGT(){
     return dinucNSyn_stat[3][9][11];
 }
 double GetDinucNSyn_GCTA(){
     return dinucNSyn_stat[3][9][12];
 }
 double GetDinucNSyn_GCTC(){
     return dinucNSyn_stat[3][9][13];
 }
 double GetDinucNSyn_GCTG(){
     return dinucNSyn_stat[3][9][14];
 }
 double GetDinucNSyn_GCTT(){
     return dinucNSyn_stat[3][9][15];
 }
 double GetDinucNSyn_GGAA(){
     return dinucNSyn_stat[3][10][0];
 }
 double GetDinucNSyn_GGAC(){
     return dinucNSyn_stat[3][10][1];
 }
 double GetDinucNSyn_GGAG(){
     return dinucNSyn_stat[3][10][2];
 }
 double GetDinucNSyn_GGAT(){
     return dinucNSyn_stat[3][10][3];
 }
 double GetDinucNSyn_GGCA(){
     return dinucNSyn_stat[3][10][4];
 }
 double GetDinucNSyn_GGCC(){
     return dinucNSyn_stat[3][10][5];
 }
 double GetDinucNSyn_GGCG(){
     return dinucNSyn_stat[3][10][6];
 }
 double GetDinucNSyn_GGCT(){
     return dinucNSyn_stat[3][10][7];
 }
 double GetDinucNSyn_GGGA(){
     return dinucNSyn_stat[3][10][8];
 }
 double GetDinucNSyn_GGGC(){
     return dinucNSyn_stat[3][10][9];
 }
 double GetDinucNSyn_GGGG(){
     return dinucNSyn_stat[3][10][10];
 }
 double GetDinucNSyn_GGGT(){
     return dinucNSyn_stat[3][10][11];
 }
 double GetDinucNSyn_GGTA(){
     return dinucNSyn_stat[3][10][12];
 }
 double GetDinucNSyn_GGTC(){
     return dinucNSyn_stat[3][10][13];
 }
 double GetDinucNSyn_GGTG(){
     return dinucNSyn_stat[3][10][14];
 }
 double GetDinucNSyn_GGTT(){
     return dinucNSyn_stat[3][10][15];
 }
 double GetDinucNSyn_GTAA(){
     return dinucNSyn_stat[3][11][0];
 }
 double GetDinucNSyn_GTAC(){
     return dinucNSyn_stat[3][11][1];
 }
 double GetDinucNSyn_GTAG(){
     return dinucNSyn_stat[3][11][2];
 }
 double GetDinucNSyn_GTAT(){
     return dinucNSyn_stat[3][11][3];
 }
 double GetDinucNSyn_GTCA(){
     return dinucNSyn_stat[3][11][4];
 }
 double GetDinucNSyn_GTCC(){
     return dinucNSyn_stat[3][11][5];
 }
 double GetDinucNSyn_GTCG(){
     return dinucNSyn_stat[3][11][6];
 }
 double GetDinucNSyn_GTCT(){
     return dinucNSyn_stat[3][11][7];
 }
 double GetDinucNSyn_GTGA(){
     return dinucNSyn_stat[3][11][8];
 }
 double GetDinucNSyn_GTGC(){
     return dinucNSyn_stat[3][11][9];
 }
 double GetDinucNSyn_GTGG(){
     return dinucNSyn_stat[3][11][10];
 }
 double GetDinucNSyn_GTGT(){
     return dinucNSyn_stat[3][11][11];
 }
 double GetDinucNSyn_GTTA(){
     return dinucNSyn_stat[3][11][12];
 }
 double GetDinucNSyn_GTTC(){
     return dinucNSyn_stat[3][11][13];
 }
 double GetDinucNSyn_GTTG(){
     return dinucNSyn_stat[3][11][14];
 }
 double GetDinucNSyn_GTTT(){
     return dinucNSyn_stat[3][11][15];
 }
 double GetDinucNSyn_TAAA(){
     return dinucNSyn_stat[3][12][0];
 }
 double GetDinucNSyn_TAAC(){
     return dinucNSyn_stat[3][12][1];
 }
 double GetDinucNSyn_TAAG(){
     return dinucNSyn_stat[3][12][2];
 }
 double GetDinucNSyn_TAAT(){
     return dinucNSyn_stat[3][12][3];
 }
 double GetDinucNSyn_TACA(){
     return dinucNSyn_stat[3][12][4];
 }
 double GetDinucNSyn_TACC(){
     return dinucNSyn_stat[3][12][5];
 }
 double GetDinucNSyn_TACG(){
     return dinucNSyn_stat[3][12][6];
 }
 double GetDinucNSyn_TACT(){
     return dinucNSyn_stat[3][12][7];
 }
 double GetDinucNSyn_TAGA(){
     return dinucNSyn_stat[3][12][8];
 }
 double GetDinucNSyn_TAGC(){
     return dinucNSyn_stat[3][12][9];
 }
 double GetDinucNSyn_TAGG(){
     return dinucNSyn_stat[3][12][10];
 }
 double GetDinucNSyn_TAGT(){
     return dinucNSyn_stat[3][12][11];
 }
 double GetDinucNSyn_TATA(){
     return dinucNSyn_stat[3][12][12];
 }
 double GetDinucNSyn_TATC(){
     return dinucNSyn_stat[3][12][13];
 }
 double GetDinucNSyn_TATG(){
     return dinucNSyn_stat[3][12][14];
 }
 double GetDinucNSyn_TATT(){
     return dinucNSyn_stat[3][12][15];
 }
 double GetDinucNSyn_TCAA(){
     return dinucNSyn_stat[3][13][0];
 }
 double GetDinucNSyn_TCAC(){
     return dinucNSyn_stat[3][13][1];
 }
 double GetDinucNSyn_TCAG(){
     return dinucNSyn_stat[3][13][2];
 }
 double GetDinucNSyn_TCAT(){
     return dinucNSyn_stat[3][13][3];
 }
 double GetDinucNSyn_TCCA(){
     return dinucNSyn_stat[3][13][4];
 }
 double GetDinucNSyn_TCCC(){
     return dinucNSyn_stat[3][13][5];
 }
 double GetDinucNSyn_TCCG(){
     return dinucNSyn_stat[3][13][6];
 }
 double GetDinucNSyn_TCCT(){
     return dinucNSyn_stat[3][13][7];
 }
 double GetDinucNSyn_TCGA(){
     return dinucNSyn_stat[3][13][8];
 }
 double GetDinucNSyn_TCGC(){
     return dinucNSyn_stat[3][13][9];
 }
 double GetDinucNSyn_TCGG(){
     return dinucNSyn_stat[3][13][10];
 }
 double GetDinucNSyn_TCGT(){
     return dinucNSyn_stat[3][13][11];
 }
 double GetDinucNSyn_TCTA(){
     return dinucNSyn_stat[3][13][12];
 }
 double GetDinucNSyn_TCTC(){
     return dinucNSyn_stat[3][13][13];
 }
 double GetDinucNSyn_TCTG(){
     return dinucNSyn_stat[3][13][14];
 }
 double GetDinucNSyn_TCTT(){
     return dinucNSyn_stat[3][13][15];
 }
 double GetDinucNSyn_TGAA(){
     return dinucNSyn_stat[3][14][0];
 }
 double GetDinucNSyn_TGAC(){
     return dinucNSyn_stat[3][14][1];
 }
 double GetDinucNSyn_TGAG(){
     return dinucNSyn_stat[3][14][2];
 }
 double GetDinucNSyn_TGAT(){
     return dinucNSyn_stat[3][14][3];
 }
 double GetDinucNSyn_TGCA(){
     return dinucNSyn_stat[3][14][4];
 }
 double GetDinucNSyn_TGCC(){
     return dinucNSyn_stat[3][14][5];
 }
 double GetDinucNSyn_TGCG(){
     return dinucNSyn_stat[3][14][6];
 }
 double GetDinucNSyn_TGCT(){
     return dinucNSyn_stat[3][14][7];
 }
 double GetDinucNSyn_TGGA(){
     return dinucNSyn_stat[3][14][8];
 }
 double GetDinucNSyn_TGGC(){
     return dinucNSyn_stat[3][14][9];
 }
 double GetDinucNSyn_TGGG(){
     return dinucNSyn_stat[3][14][10];
 }
 double GetDinucNSyn_TGGT(){
     return dinucNSyn_stat[3][14][11];
 }
 double GetDinucNSyn_TGTA(){
     return dinucNSyn_stat[3][14][12];
 }
 double GetDinucNSyn_TGTC(){
     return dinucNSyn_stat[3][14][13];
 }
 double GetDinucNSyn_TGTG(){
     return dinucNSyn_stat[3][14][14];
 }
 double GetDinucNSyn_TGTT(){
     return dinucNSyn_stat[3][14][15];
 }
 double GetDinucNSyn_TTAA(){
     return dinucNSyn_stat[3][15][0];
 }
 double GetDinucNSyn_TTAC(){
     return dinucNSyn_stat[3][15][1];
 }
 double GetDinucNSyn_TTAG(){
     return dinucNSyn_stat[3][15][2];
 }
 double GetDinucNSyn_TTAT(){
     return dinucNSyn_stat[3][15][3];
 }
 double GetDinucNSyn_TTCA(){
     return dinucNSyn_stat[3][15][4];
 }
 double GetDinucNSyn_TTCC(){
     return dinucNSyn_stat[3][15][5];
 }
 double GetDinucNSyn_TTCG(){
     return dinucNSyn_stat[3][15][6];
 }
 double GetDinucNSyn_TTCT(){
     return dinucNSyn_stat[3][15][7];
 }
 double GetDinucNSyn_TTGA(){
     return dinucNSyn_stat[3][15][8];
 }
 double GetDinucNSyn_TTGC(){
     return dinucNSyn_stat[3][15][9];
 }
 double GetDinucNSyn_TTGG(){
     return dinucNSyn_stat[3][15][10];
 }
 double GetDinucNSyn_TTGT(){
     return dinucNSyn_stat[3][15][11];
 }
 double GetDinucNSyn_TTTA(){
     return dinucNSyn_stat[3][15][12];
 }
 double GetDinucNSyn_TTTC(){
     return dinucNSyn_stat[3][15][13];
 }
 double GetDinucNSyn_TTTG(){
     return dinucNSyn_stat[3][15][14];
 }
 double GetDinucNSyn_TTTT(){
     return dinucNSyn_stat[3][15][15];
 }
 double GetDinucNSyn12_AAAA(){
     return dinucNSyn_stat[0][0][0];
 }
 double GetDinucNSyn12_AAAC(){
     return dinucNSyn_stat[0][0][1];
 }
 double GetDinucNSyn12_AAAG(){
     return dinucNSyn_stat[0][0][2];
 }
 double GetDinucNSyn12_AAAT(){
     return dinucNSyn_stat[0][0][3];
 }
 double GetDinucNSyn12_AACA(){
     return dinucNSyn_stat[0][0][4];
 }
 double GetDinucNSyn12_AACC(){
     return dinucNSyn_stat[0][0][5];
 }
 double GetDinucNSyn12_AACG(){
     return dinucNSyn_stat[0][0][6];
 }
 double GetDinucNSyn12_AACT(){
     return dinucNSyn_stat[0][0][7];
 }
 double GetDinucNSyn12_AAGA(){
     return dinucNSyn_stat[0][0][8];
 }
 double GetDinucNSyn12_AAGC(){
     return dinucNSyn_stat[0][0][9];
 }
 double GetDinucNSyn12_AAGG(){
     return dinucNSyn_stat[0][0][10];
 }
 double GetDinucNSyn12_AAGT(){
     return dinucNSyn_stat[0][0][11];
 }
 double GetDinucNSyn12_AATA(){
     return dinucNSyn_stat[0][0][12];
 }
 double GetDinucNSyn12_AATC(){
     return dinucNSyn_stat[0][0][13];
 }
 double GetDinucNSyn12_AATG(){
     return dinucNSyn_stat[0][0][14];
 }
 double GetDinucNSyn12_AATT(){
     return dinucNSyn_stat[0][0][15];
 }
 double GetDinucNSyn12_ACAA(){
     return dinucNSyn_stat[0][1][0];
 }
 double GetDinucNSyn12_ACAC(){
     return dinucNSyn_stat[0][1][1];
 }
 double GetDinucNSyn12_ACAG(){
     return dinucNSyn_stat[0][1][2];
 }
 double GetDinucNSyn12_ACAT(){
     return dinucNSyn_stat[0][1][3];
 }
 double GetDinucNSyn12_ACCA(){
     return dinucNSyn_stat[0][1][4];
 }
 double GetDinucNSyn12_ACCC(){
     return dinucNSyn_stat[0][1][5];
 }
 double GetDinucNSyn12_ACCG(){
     return dinucNSyn_stat[0][1][6];
 }
 double GetDinucNSyn12_ACCT(){
     return dinucNSyn_stat[0][1][7];
 }
 double GetDinucNSyn12_ACGA(){
     return dinucNSyn_stat[0][1][8];
 }
 double GetDinucNSyn12_ACGC(){
     return dinucNSyn_stat[0][1][9];
 }
 double GetDinucNSyn12_ACGG(){
     return dinucNSyn_stat[0][1][10];
 }
 double GetDinucNSyn12_ACGT(){
     return dinucNSyn_stat[0][1][11];
 }
 double GetDinucNSyn12_ACTA(){
     return dinucNSyn_stat[0][1][12];
 }
 double GetDinucNSyn12_ACTC(){
     return dinucNSyn_stat[0][1][13];
 }
 double GetDinucNSyn12_ACTG(){
     return dinucNSyn_stat[0][1][14];
 }
 double GetDinucNSyn12_ACTT(){
     return dinucNSyn_stat[0][1][15];
 }
 double GetDinucNSyn12_AGAA(){
     return dinucNSyn_stat[0][2][0];
 }
 double GetDinucNSyn12_AGAC(){
     return dinucNSyn_stat[0][2][1];
 }
 double GetDinucNSyn12_AGAG(){
     return dinucNSyn_stat[0][2][2];
 }
 double GetDinucNSyn12_AGAT(){
     return dinucNSyn_stat[0][2][3];
 }
 double GetDinucNSyn12_AGCA(){
     return dinucNSyn_stat[0][2][4];
 }
 double GetDinucNSyn12_AGCC(){
     return dinucNSyn_stat[0][2][5];
 }
 double GetDinucNSyn12_AGCG(){
     return dinucNSyn_stat[0][2][6];
 }
 double GetDinucNSyn12_AGCT(){
     return dinucNSyn_stat[0][2][7];
 }
 double GetDinucNSyn12_AGGA(){
     return dinucNSyn_stat[0][2][8];
 }
 double GetDinucNSyn12_AGGC(){
     return dinucNSyn_stat[0][2][9];
 }
 double GetDinucNSyn12_AGGG(){
     return dinucNSyn_stat[0][2][10];
 }
 double GetDinucNSyn12_AGGT(){
     return dinucNSyn_stat[0][2][11];
 }
 double GetDinucNSyn12_AGTA(){
     return dinucNSyn_stat[0][2][12];
 }
 double GetDinucNSyn12_AGTC(){
     return dinucNSyn_stat[0][2][13];
 }
 double GetDinucNSyn12_AGTG(){
     return dinucNSyn_stat[0][2][14];
 }
 double GetDinucNSyn12_AGTT(){
     return dinucNSyn_stat[0][2][15];
 }
 double GetDinucNSyn12_ATAA(){
     return dinucNSyn_stat[0][3][0];
 }
 double GetDinucNSyn12_ATAC(){
     return dinucNSyn_stat[0][3][1];
 }
 double GetDinucNSyn12_ATAG(){
     return dinucNSyn_stat[0][3][2];
 }
 double GetDinucNSyn12_ATAT(){
     return dinucNSyn_stat[0][3][3];
 }
 double GetDinucNSyn12_ATCA(){
     return dinucNSyn_stat[0][3][4];
 }
 double GetDinucNSyn12_ATCC(){
     return dinucNSyn_stat[0][3][5];
 }
 double GetDinucNSyn12_ATCG(){
     return dinucNSyn_stat[0][3][6];
 }
 double GetDinucNSyn12_ATCT(){
     return dinucNSyn_stat[0][3][7];
 }
 double GetDinucNSyn12_ATGA(){
     return dinucNSyn_stat[0][3][8];
 }
 double GetDinucNSyn12_ATGC(){
     return dinucNSyn_stat[0][3][9];
 }
 double GetDinucNSyn12_ATGG(){
     return dinucNSyn_stat[0][3][10];
 }
 double GetDinucNSyn12_ATGT(){
     return dinucNSyn_stat[0][3][11];
 }
 double GetDinucNSyn12_ATTA(){
     return dinucNSyn_stat[0][3][12];
 }
 double GetDinucNSyn12_ATTC(){
     return dinucNSyn_stat[0][3][13];
 }
 double GetDinucNSyn12_ATTG(){
     return dinucNSyn_stat[0][3][14];
 }
 double GetDinucNSyn12_ATTT(){
     return dinucNSyn_stat[0][3][15];
 }
 double GetDinucNSyn12_CAAA(){
     return dinucNSyn_stat[0][4][0];
 }
 double GetDinucNSyn12_CAAC(){
     return dinucNSyn_stat[0][4][1];
 }
 double GetDinucNSyn12_CAAG(){
     return dinucNSyn_stat[0][4][2];
 }
 double GetDinucNSyn12_CAAT(){
     return dinucNSyn_stat[0][4][3];
 }
 double GetDinucNSyn12_CACA(){
     return dinucNSyn_stat[0][4][4];
 }
 double GetDinucNSyn12_CACC(){
     return dinucNSyn_stat[0][4][5];
 }
 double GetDinucNSyn12_CACG(){
     return dinucNSyn_stat[0][4][6];
 }
 double GetDinucNSyn12_CACT(){
     return dinucNSyn_stat[0][4][7];
 }
 double GetDinucNSyn12_CAGA(){
     return dinucNSyn_stat[0][4][8];
 }
 double GetDinucNSyn12_CAGC(){
     return dinucNSyn_stat[0][4][9];
 }
 double GetDinucNSyn12_CAGG(){
     return dinucNSyn_stat[0][4][10];
 }
 double GetDinucNSyn12_CAGT(){
     return dinucNSyn_stat[0][4][11];
 }
 double GetDinucNSyn12_CATA(){
     return dinucNSyn_stat[0][4][12];
 }
 double GetDinucNSyn12_CATC(){
     return dinucNSyn_stat[0][4][13];
 }
 double GetDinucNSyn12_CATG(){
     return dinucNSyn_stat[0][4][14];
 }
 double GetDinucNSyn12_CATT(){
     return dinucNSyn_stat[0][4][15];
 }
 double GetDinucNSyn12_CCAA(){
     return dinucNSyn_stat[0][5][0];
 }
 double GetDinucNSyn12_CCAC(){
     return dinucNSyn_stat[0][5][1];
 }
 double GetDinucNSyn12_CCAG(){
     return dinucNSyn_stat[0][5][2];
 }
 double GetDinucNSyn12_CCAT(){
     return dinucNSyn_stat[0][5][3];
 }
 double GetDinucNSyn12_CCCA(){
     return dinucNSyn_stat[0][5][4];
 }
 double GetDinucNSyn12_CCCC(){
     return dinucNSyn_stat[0][5][5];
 }
 double GetDinucNSyn12_CCCG(){
     return dinucNSyn_stat[0][5][6];
 }
 double GetDinucNSyn12_CCCT(){
     return dinucNSyn_stat[0][5][7];
 }
 double GetDinucNSyn12_CCGA(){
     return dinucNSyn_stat[0][5][8];
 }
 double GetDinucNSyn12_CCGC(){
     return dinucNSyn_stat[0][5][9];
 }
 double GetDinucNSyn12_CCGG(){
     return dinucNSyn_stat[0][5][10];
 }
 double GetDinucNSyn12_CCGT(){
     return dinucNSyn_stat[0][5][11];
 }
 double GetDinucNSyn12_CCTA(){
     return dinucNSyn_stat[0][5][12];
 }
 double GetDinucNSyn12_CCTC(){
     return dinucNSyn_stat[0][5][13];
 }
 double GetDinucNSyn12_CCTG(){
     return dinucNSyn_stat[0][5][14];
 }
 double GetDinucNSyn12_CCTT(){
     return dinucNSyn_stat[0][5][15];
 }
 double GetDinucNSyn12_CGAA(){
     return dinucNSyn_stat[0][6][0];
 }
 double GetDinucNSyn12_CGAC(){
     return dinucNSyn_stat[0][6][1];
 }
 double GetDinucNSyn12_CGAG(){
     return dinucNSyn_stat[0][6][2];
 }
 double GetDinucNSyn12_CGAT(){
     return dinucNSyn_stat[0][6][3];
 }
 double GetDinucNSyn12_CGCA(){
     return dinucNSyn_stat[0][6][4];
 }
 double GetDinucNSyn12_CGCC(){
     return dinucNSyn_stat[0][6][5];
 }
 double GetDinucNSyn12_CGCG(){
     return dinucNSyn_stat[0][6][6];
 }
 double GetDinucNSyn12_CGCT(){
     return dinucNSyn_stat[0][6][7];
 }
 double GetDinucNSyn12_CGGA(){
     return dinucNSyn_stat[0][6][8];
 }
 double GetDinucNSyn12_CGGC(){
     return dinucNSyn_stat[0][6][9];
 }
 double GetDinucNSyn12_CGGG(){
     return dinucNSyn_stat[0][6][10];
 }
 double GetDinucNSyn12_CGGT(){
     return dinucNSyn_stat[0][6][11];
 }
 double GetDinucNSyn12_CGTA(){
     return dinucNSyn_stat[0][6][12];
 }
 double GetDinucNSyn12_CGTC(){
     return dinucNSyn_stat[0][6][13];
 }
 double GetDinucNSyn12_CGTG(){
     return dinucNSyn_stat[0][6][14];
 }
 double GetDinucNSyn12_CGTT(){
     return dinucNSyn_stat[0][6][15];
 }
 double GetDinucNSyn12_CTAA(){
     return dinucNSyn_stat[0][7][0];
 }
 double GetDinucNSyn12_CTAC(){
     return dinucNSyn_stat[0][7][1];
 }
 double GetDinucNSyn12_CTAG(){
     return dinucNSyn_stat[0][7][2];
 }
 double GetDinucNSyn12_CTAT(){
     return dinucNSyn_stat[0][7][3];
 }
 double GetDinucNSyn12_CTCA(){
     return dinucNSyn_stat[0][7][4];
 }
 double GetDinucNSyn12_CTCC(){
     return dinucNSyn_stat[0][7][5];
 }
 double GetDinucNSyn12_CTCG(){
     return dinucNSyn_stat[0][7][6];
 }
 double GetDinucNSyn12_CTCT(){
     return dinucNSyn_stat[0][7][7];
 }
 double GetDinucNSyn12_CTGA(){
     return dinucNSyn_stat[0][7][8];
 }
 double GetDinucNSyn12_CTGC(){
     return dinucNSyn_stat[0][7][9];
 }
 double GetDinucNSyn12_CTGG(){
     return dinucNSyn_stat[0][7][10];
 }
 double GetDinucNSyn12_CTGT(){
     return dinucNSyn_stat[0][7][11];
 }
 double GetDinucNSyn12_CTTA(){
     return dinucNSyn_stat[0][7][12];
 }
 double GetDinucNSyn12_CTTC(){
     return dinucNSyn_stat[0][7][13];
 }
 double GetDinucNSyn12_CTTG(){
     return dinucNSyn_stat[0][7][14];
 }
 double GetDinucNSyn12_CTTT(){
     return dinucNSyn_stat[0][7][15];
 }
 double GetDinucNSyn12_GAAA(){
     return dinucNSyn_stat[0][8][0];
 }
 double GetDinucNSyn12_GAAC(){
     return dinucNSyn_stat[0][8][1];
 }
 double GetDinucNSyn12_GAAG(){
     return dinucNSyn_stat[0][8][2];
 }
 double GetDinucNSyn12_GAAT(){
     return dinucNSyn_stat[0][8][3];
 }
 double GetDinucNSyn12_GACA(){
     return dinucNSyn_stat[0][8][4];
 }
 double GetDinucNSyn12_GACC(){
     return dinucNSyn_stat[0][8][5];
 }
 double GetDinucNSyn12_GACG(){
     return dinucNSyn_stat[0][8][6];
 }
 double GetDinucNSyn12_GACT(){
     return dinucNSyn_stat[0][8][7];
 }
 double GetDinucNSyn12_GAGA(){
     return dinucNSyn_stat[0][8][8];
 }
 double GetDinucNSyn12_GAGC(){
     return dinucNSyn_stat[0][8][9];
 }
 double GetDinucNSyn12_GAGG(){
     return dinucNSyn_stat[0][8][10];
 }
 double GetDinucNSyn12_GAGT(){
     return dinucNSyn_stat[0][8][11];
 }
 double GetDinucNSyn12_GATA(){
     return dinucNSyn_stat[0][8][12];
 }
 double GetDinucNSyn12_GATC(){
     return dinucNSyn_stat[0][8][13];
 }
 double GetDinucNSyn12_GATG(){
     return dinucNSyn_stat[0][8][14];
 }
 double GetDinucNSyn12_GATT(){
     return dinucNSyn_stat[0][8][15];
 }
 double GetDinucNSyn12_GCAA(){
     return dinucNSyn_stat[0][9][0];
 }
 double GetDinucNSyn12_GCAC(){
     return dinucNSyn_stat[0][9][1];
 }
 double GetDinucNSyn12_GCAG(){
     return dinucNSyn_stat[0][9][2];
 }
 double GetDinucNSyn12_GCAT(){
     return dinucNSyn_stat[0][9][3];
 }
 double GetDinucNSyn12_GCCA(){
     return dinucNSyn_stat[0][9][4];
 }
 double GetDinucNSyn12_GCCC(){
     return dinucNSyn_stat[0][9][5];
 }
 double GetDinucNSyn12_GCCG(){
     return dinucNSyn_stat[0][9][6];
 }
 double GetDinucNSyn12_GCCT(){
     return dinucNSyn_stat[0][9][7];
 }
 double GetDinucNSyn12_GCGA(){
     return dinucNSyn_stat[0][9][8];
 }
 double GetDinucNSyn12_GCGC(){
     return dinucNSyn_stat[0][9][9];
 }
 double GetDinucNSyn12_GCGG(){
     return dinucNSyn_stat[0][9][10];
 }
 double GetDinucNSyn12_GCGT(){
     return dinucNSyn_stat[0][9][11];
 }
 double GetDinucNSyn12_GCTA(){
     return dinucNSyn_stat[0][9][12];
 }
 double GetDinucNSyn12_GCTC(){
     return dinucNSyn_stat[0][9][13];
 }
 double GetDinucNSyn12_GCTG(){
     return dinucNSyn_stat[0][9][14];
 }
 double GetDinucNSyn12_GCTT(){
     return dinucNSyn_stat[0][9][15];
 }
 double GetDinucNSyn12_GGAA(){
     return dinucNSyn_stat[0][10][0];
 }
 double GetDinucNSyn12_GGAC(){
     return dinucNSyn_stat[0][10][1];
 }
 double GetDinucNSyn12_GGAG(){
     return dinucNSyn_stat[0][10][2];
 }
 double GetDinucNSyn12_GGAT(){
     return dinucNSyn_stat[0][10][3];
 }
 double GetDinucNSyn12_GGCA(){
     return dinucNSyn_stat[0][10][4];
 }
 double GetDinucNSyn12_GGCC(){
     return dinucNSyn_stat[0][10][5];
 }
 double GetDinucNSyn12_GGCG(){
     return dinucNSyn_stat[0][10][6];
 }
 double GetDinucNSyn12_GGCT(){
     return dinucNSyn_stat[0][10][7];
 }
 double GetDinucNSyn12_GGGA(){
     return dinucNSyn_stat[0][10][8];
 }
 double GetDinucNSyn12_GGGC(){
     return dinucNSyn_stat[0][10][9];
 }
 double GetDinucNSyn12_GGGG(){
     return dinucNSyn_stat[0][10][10];
 }
 double GetDinucNSyn12_GGGT(){
     return dinucNSyn_stat[0][10][11];
 }
 double GetDinucNSyn12_GGTA(){
     return dinucNSyn_stat[0][10][12];
 }
 double GetDinucNSyn12_GGTC(){
     return dinucNSyn_stat[0][10][13];
 }
 double GetDinucNSyn12_GGTG(){
     return dinucNSyn_stat[0][10][14];
 }
 double GetDinucNSyn12_GGTT(){
     return dinucNSyn_stat[0][10][15];
 }
 double GetDinucNSyn12_GTAA(){
     return dinucNSyn_stat[0][11][0];
 }
 double GetDinucNSyn12_GTAC(){
     return dinucNSyn_stat[0][11][1];
 }
 double GetDinucNSyn12_GTAG(){
     return dinucNSyn_stat[0][11][2];
 }
 double GetDinucNSyn12_GTAT(){
     return dinucNSyn_stat[0][11][3];
 }
 double GetDinucNSyn12_GTCA(){
     return dinucNSyn_stat[0][11][4];
 }
 double GetDinucNSyn12_GTCC(){
     return dinucNSyn_stat[0][11][5];
 }
 double GetDinucNSyn12_GTCG(){
     return dinucNSyn_stat[0][11][6];
 }
 double GetDinucNSyn12_GTCT(){
     return dinucNSyn_stat[0][11][7];
 }
 double GetDinucNSyn12_GTGA(){
     return dinucNSyn_stat[0][11][8];
 }
 double GetDinucNSyn12_GTGC(){
     return dinucNSyn_stat[0][11][9];
 }
 double GetDinucNSyn12_GTGG(){
     return dinucNSyn_stat[0][11][10];
 }
 double GetDinucNSyn12_GTGT(){
     return dinucNSyn_stat[0][11][11];
 }
 double GetDinucNSyn12_GTTA(){
     return dinucNSyn_stat[0][11][12];
 }
 double GetDinucNSyn12_GTTC(){
     return dinucNSyn_stat[0][11][13];
 }
 double GetDinucNSyn12_GTTG(){
     return dinucNSyn_stat[0][11][14];
 }
 double GetDinucNSyn12_GTTT(){
     return dinucNSyn_stat[0][11][15];
 }
 double GetDinucNSyn12_TAAA(){
     return dinucNSyn_stat[0][12][0];
 }
 double GetDinucNSyn12_TAAC(){
     return dinucNSyn_stat[0][12][1];
 }
 double GetDinucNSyn12_TAAG(){
     return dinucNSyn_stat[0][12][2];
 }
 double GetDinucNSyn12_TAAT(){
     return dinucNSyn_stat[0][12][3];
 }
 double GetDinucNSyn12_TACA(){
     return dinucNSyn_stat[0][12][4];
 }
 double GetDinucNSyn12_TACC(){
     return dinucNSyn_stat[0][12][5];
 }
 double GetDinucNSyn12_TACG(){
     return dinucNSyn_stat[0][12][6];
 }
 double GetDinucNSyn12_TACT(){
     return dinucNSyn_stat[0][12][7];
 }
 double GetDinucNSyn12_TAGA(){
     return dinucNSyn_stat[0][12][8];
 }
 double GetDinucNSyn12_TAGC(){
     return dinucNSyn_stat[0][12][9];
 }
 double GetDinucNSyn12_TAGG(){
     return dinucNSyn_stat[0][12][10];
 }
 double GetDinucNSyn12_TAGT(){
     return dinucNSyn_stat[0][12][11];
 }
 double GetDinucNSyn12_TATA(){
     return dinucNSyn_stat[0][12][12];
 }
 double GetDinucNSyn12_TATC(){
     return dinucNSyn_stat[0][12][13];
 }
 double GetDinucNSyn12_TATG(){
     return dinucNSyn_stat[0][12][14];
 }
 double GetDinucNSyn12_TATT(){
     return dinucNSyn_stat[0][12][15];
 }
 double GetDinucNSyn12_TCAA(){
     return dinucNSyn_stat[0][13][0];
 }
 double GetDinucNSyn12_TCAC(){
     return dinucNSyn_stat[0][13][1];
 }
 double GetDinucNSyn12_TCAG(){
     return dinucNSyn_stat[0][13][2];
 }
 double GetDinucNSyn12_TCAT(){
     return dinucNSyn_stat[0][13][3];
 }
 double GetDinucNSyn12_TCCA(){
     return dinucNSyn_stat[0][13][4];
 }
 double GetDinucNSyn12_TCCC(){
     return dinucNSyn_stat[0][13][5];
 }
 double GetDinucNSyn12_TCCG(){
     return dinucNSyn_stat[0][13][6];
 }
 double GetDinucNSyn12_TCCT(){
     return dinucNSyn_stat[0][13][7];
 }
 double GetDinucNSyn12_TCGA(){
     return dinucNSyn_stat[0][13][8];
 }
 double GetDinucNSyn12_TCGC(){
     return dinucNSyn_stat[0][13][9];
 }
 double GetDinucNSyn12_TCGG(){
     return dinucNSyn_stat[0][13][10];
 }
 double GetDinucNSyn12_TCGT(){
     return dinucNSyn_stat[0][13][11];
 }
 double GetDinucNSyn12_TCTA(){
     return dinucNSyn_stat[0][13][12];
 }
 double GetDinucNSyn12_TCTC(){
     return dinucNSyn_stat[0][13][13];
 }
 double GetDinucNSyn12_TCTG(){
     return dinucNSyn_stat[0][13][14];
 }
 double GetDinucNSyn12_TCTT(){
     return dinucNSyn_stat[0][13][15];
 }
 double GetDinucNSyn12_TGAA(){
     return dinucNSyn_stat[0][14][0];
 }
 double GetDinucNSyn12_TGAC(){
     return dinucNSyn_stat[0][14][1];
 }
 double GetDinucNSyn12_TGAG(){
     return dinucNSyn_stat[0][14][2];
 }
 double GetDinucNSyn12_TGAT(){
     return dinucNSyn_stat[0][14][3];
 }
 double GetDinucNSyn12_TGCA(){
     return dinucNSyn_stat[0][14][4];
 }
 double GetDinucNSyn12_TGCC(){
     return dinucNSyn_stat[0][14][5];
 }
 double GetDinucNSyn12_TGCG(){
     return dinucNSyn_stat[0][14][6];
 }
 double GetDinucNSyn12_TGCT(){
     return dinucNSyn_stat[0][14][7];
 }
 double GetDinucNSyn12_TGGA(){
     return dinucNSyn_stat[0][14][8];
 }
 double GetDinucNSyn12_TGGC(){
     return dinucNSyn_stat[0][14][9];
 }
 double GetDinucNSyn12_TGGG(){
     return dinucNSyn_stat[0][14][10];
 }
 double GetDinucNSyn12_TGGT(){
     return dinucNSyn_stat[0][14][11];
 }
 double GetDinucNSyn12_TGTA(){
     return dinucNSyn_stat[0][14][12];
 }
 double GetDinucNSyn12_TGTC(){
     return dinucNSyn_stat[0][14][13];
 }
 double GetDinucNSyn12_TGTG(){
     return dinucNSyn_stat[0][14][14];
 }
 double GetDinucNSyn12_TGTT(){
     return dinucNSyn_stat[0][14][15];
 }
 double GetDinucNSyn12_TTAA(){
     return dinucNSyn_stat[0][15][0];
 }
 double GetDinucNSyn12_TTAC(){
     return dinucNSyn_stat[0][15][1];
 }
 double GetDinucNSyn12_TTAG(){
     return dinucNSyn_stat[0][15][2];
 }
 double GetDinucNSyn12_TTAT(){
     return dinucNSyn_stat[0][15][3];
 }
 double GetDinucNSyn12_TTCA(){
     return dinucNSyn_stat[0][15][4];
 }
 double GetDinucNSyn12_TTCC(){
     return dinucNSyn_stat[0][15][5];
 }
 double GetDinucNSyn12_TTCG(){
     return dinucNSyn_stat[0][15][6];
 }
 double GetDinucNSyn12_TTCT(){
     return dinucNSyn_stat[0][15][7];
 }
 double GetDinucNSyn12_TTGA(){
     return dinucNSyn_stat[0][15][8];
 }
 double GetDinucNSyn12_TTGC(){
     return dinucNSyn_stat[0][15][9];
 }
 double GetDinucNSyn12_TTGG(){
     return dinucNSyn_stat[0][15][10];
 }
 double GetDinucNSyn12_TTGT(){
     return dinucNSyn_stat[0][15][11];
 }
 double GetDinucNSyn12_TTTA(){
     return dinucNSyn_stat[0][15][12];
 }
 double GetDinucNSyn12_TTTC(){
     return dinucNSyn_stat[0][15][13];
 }
 double GetDinucNSyn12_TTTG(){
     return dinucNSyn_stat[0][15][14];
 }
 double GetDinucNSyn12_TTTT(){
     return dinucNSyn_stat[0][15][15];
 }

 double GetDinucNSyn23_AAAA(){
     return dinucNSyn_stat[1][0][0];
 }
 double GetDinucNSyn23_AAAC(){
     return dinucNSyn_stat[1][0][1];
 }
 double GetDinucNSyn23_AAAG(){
     return dinucNSyn_stat[1][0][2];
 }
 double GetDinucNSyn23_AAAT(){
     return dinucNSyn_stat[1][0][3];
 }
 double GetDinucNSyn23_AACA(){
     return dinucNSyn_stat[1][0][4];
 }
 double GetDinucNSyn23_AACC(){
     return dinucNSyn_stat[1][0][5];
 }
 double GetDinucNSyn23_AACG(){
     return dinucNSyn_stat[1][0][6];
 }
 double GetDinucNSyn23_AACT(){
     return dinucNSyn_stat[1][0][7];
 }
 double GetDinucNSyn23_AAGA(){
     return dinucNSyn_stat[1][0][8];
 }
 double GetDinucNSyn23_AAGC(){
     return dinucNSyn_stat[1][0][9];
 }
 double GetDinucNSyn23_AAGG(){
     return dinucNSyn_stat[1][0][10];
 }
 double GetDinucNSyn23_AAGT(){
     return dinucNSyn_stat[1][0][11];
 }
 double GetDinucNSyn23_AATA(){
     return dinucNSyn_stat[1][0][12];
 }
 double GetDinucNSyn23_AATC(){
     return dinucNSyn_stat[1][0][13];
 }
 double GetDinucNSyn23_AATG(){
     return dinucNSyn_stat[1][0][14];
 }
 double GetDinucNSyn23_AATT(){
     return dinucNSyn_stat[1][0][15];
 }
 double GetDinucNSyn23_ACAA(){
     return dinucNSyn_stat[1][1][0];
 }
 double GetDinucNSyn23_ACAC(){
     return dinucNSyn_stat[1][1][1];
 }
 double GetDinucNSyn23_ACAG(){
     return dinucNSyn_stat[1][1][2];
 }
 double GetDinucNSyn23_ACAT(){
     return dinucNSyn_stat[1][1][3];
 }
 double GetDinucNSyn23_ACCA(){
     return dinucNSyn_stat[1][1][4];
 }
 double GetDinucNSyn23_ACCC(){
     return dinucNSyn_stat[1][1][5];
 }
 double GetDinucNSyn23_ACCG(){
     return dinucNSyn_stat[1][1][6];
 }
 double GetDinucNSyn23_ACCT(){
     return dinucNSyn_stat[1][1][7];
 }
 double GetDinucNSyn23_ACGA(){
     return dinucNSyn_stat[1][1][8];
 }
 double GetDinucNSyn23_ACGC(){
     return dinucNSyn_stat[1][1][9];
 }
 double GetDinucNSyn23_ACGG(){
     return dinucNSyn_stat[1][1][10];
 }
 double GetDinucNSyn23_ACGT(){
     return dinucNSyn_stat[1][1][11];
 }
 double GetDinucNSyn23_ACTA(){
     return dinucNSyn_stat[1][1][12];
 }
 double GetDinucNSyn23_ACTC(){
     return dinucNSyn_stat[1][1][13];
 }
 double GetDinucNSyn23_ACTG(){
     return dinucNSyn_stat[1][1][14];
 }
 double GetDinucNSyn23_ACTT(){
     return dinucNSyn_stat[1][1][15];
 }
 double GetDinucNSyn23_AGAA(){
     return dinucNSyn_stat[1][2][0];
 }
 double GetDinucNSyn23_AGAC(){
     return dinucNSyn_stat[1][2][1];
 }
 double GetDinucNSyn23_AGAG(){
     return dinucNSyn_stat[1][2][2];
 }
 double GetDinucNSyn23_AGAT(){
     return dinucNSyn_stat[1][2][3];
 }
 double GetDinucNSyn23_AGCA(){
     return dinucNSyn_stat[1][2][4];
 }
 double GetDinucNSyn23_AGCC(){
     return dinucNSyn_stat[1][2][5];
 }
 double GetDinucNSyn23_AGCG(){
     return dinucNSyn_stat[1][2][6];
 }
 double GetDinucNSyn23_AGCT(){
     return dinucNSyn_stat[1][2][7];
 }
 double GetDinucNSyn23_AGGA(){
     return dinucNSyn_stat[1][2][8];
 }
 double GetDinucNSyn23_AGGC(){
     return dinucNSyn_stat[1][2][9];
 }
 double GetDinucNSyn23_AGGG(){
     return dinucNSyn_stat[1][2][10];
 }
 double GetDinucNSyn23_AGGT(){
     return dinucNSyn_stat[1][2][11];
 }
 double GetDinucNSyn23_AGTA(){
     return dinucNSyn_stat[1][2][12];
 }
 double GetDinucNSyn23_AGTC(){
     return dinucNSyn_stat[1][2][13];
 }
 double GetDinucNSyn23_AGTG(){
     return dinucNSyn_stat[1][2][14];
 }
 double GetDinucNSyn23_AGTT(){
     return dinucNSyn_stat[1][2][15];
 }
 double GetDinucNSyn23_ATAA(){
     return dinucNSyn_stat[1][3][0];
 }
 double GetDinucNSyn23_ATAC(){
     return dinucNSyn_stat[1][3][1];
 }
 double GetDinucNSyn23_ATAG(){
     return dinucNSyn_stat[1][3][2];
 }
 double GetDinucNSyn23_ATAT(){
     return dinucNSyn_stat[1][3][3];
 }
 double GetDinucNSyn23_ATCA(){
     return dinucNSyn_stat[1][3][4];
 }
 double GetDinucNSyn23_ATCC(){
     return dinucNSyn_stat[1][3][5];
 }
 double GetDinucNSyn23_ATCG(){
     return dinucNSyn_stat[1][3][6];
 }
 double GetDinucNSyn23_ATCT(){
     return dinucNSyn_stat[1][3][7];
 }
 double GetDinucNSyn23_ATGA(){
     return dinucNSyn_stat[1][3][8];
 }
 double GetDinucNSyn23_ATGC(){
     return dinucNSyn_stat[1][3][9];
 }
 double GetDinucNSyn23_ATGG(){
     return dinucNSyn_stat[1][3][10];
 }
 double GetDinucNSyn23_ATGT(){
     return dinucNSyn_stat[1][3][11];
 }
 double GetDinucNSyn23_ATTA(){
     return dinucNSyn_stat[1][3][12];
 }
 double GetDinucNSyn23_ATTC(){
     return dinucNSyn_stat[1][3][13];
 }
 double GetDinucNSyn23_ATTG(){
     return dinucNSyn_stat[1][3][14];
 }
 double GetDinucNSyn23_ATTT(){
     return dinucNSyn_stat[1][3][15];
 }
 double GetDinucNSyn23_CAAA(){
     return dinucNSyn_stat[1][4][0];
 }
 double GetDinucNSyn23_CAAC(){
     return dinucNSyn_stat[1][4][1];
 }
 double GetDinucNSyn23_CAAG(){
     return dinucNSyn_stat[1][4][2];
 }
 double GetDinucNSyn23_CAAT(){
     return dinucNSyn_stat[1][4][3];
 }
 double GetDinucNSyn23_CACA(){
     return dinucNSyn_stat[1][4][4];
 }
 double GetDinucNSyn23_CACC(){
     return dinucNSyn_stat[1][4][5];
 }
 double GetDinucNSyn23_CACG(){
     return dinucNSyn_stat[1][4][6];
 }
 double GetDinucNSyn23_CACT(){
     return dinucNSyn_stat[1][4][7];
 }
 double GetDinucNSyn23_CAGA(){
     return dinucNSyn_stat[1][4][8];
 }
 double GetDinucNSyn23_CAGC(){
     return dinucNSyn_stat[1][4][9];
 }
 double GetDinucNSyn23_CAGG(){
     return dinucNSyn_stat[1][4][10];
 }
 double GetDinucNSyn23_CAGT(){
     return dinucNSyn_stat[1][4][11];
 }
 double GetDinucNSyn23_CATA(){
     return dinucNSyn_stat[1][4][12];
 }
 double GetDinucNSyn23_CATC(){
     return dinucNSyn_stat[1][4][13];
 }
 double GetDinucNSyn23_CATG(){
     return dinucNSyn_stat[1][4][14];
 }
 double GetDinucNSyn23_CATT(){
     return dinucNSyn_stat[1][4][15];
 }
 double GetDinucNSyn23_CCAA(){
     return dinucNSyn_stat[1][5][0];
 }
 double GetDinucNSyn23_CCAC(){
     return dinucNSyn_stat[1][5][1];
 }
 double GetDinucNSyn23_CCAG(){
     return dinucNSyn_stat[1][5][2];
 }
 double GetDinucNSyn23_CCAT(){
     return dinucNSyn_stat[1][5][3];
 }
 double GetDinucNSyn23_CCCA(){
     return dinucNSyn_stat[1][5][4];
 }
 double GetDinucNSyn23_CCCC(){
     return dinucNSyn_stat[1][5][5];
 }
 double GetDinucNSyn23_CCCG(){
     return dinucNSyn_stat[1][5][6];
 }
 double GetDinucNSyn23_CCCT(){
     return dinucNSyn_stat[1][5][7];
 }
 double GetDinucNSyn23_CCGA(){
     return dinucNSyn_stat[1][5][8];
 }
 double GetDinucNSyn23_CCGC(){
     return dinucNSyn_stat[1][5][9];
 }
 double GetDinucNSyn23_CCGG(){
     return dinucNSyn_stat[1][5][10];
 }
 double GetDinucNSyn23_CCGT(){
     return dinucNSyn_stat[1][5][11];
 }
 double GetDinucNSyn23_CCTA(){
     return dinucNSyn_stat[1][5][12];
 }
 double GetDinucNSyn23_CCTC(){
     return dinucNSyn_stat[1][5][13];
 }
 double GetDinucNSyn23_CCTG(){
     return dinucNSyn_stat[1][5][14];
 }
 double GetDinucNSyn23_CCTT(){
     return dinucNSyn_stat[1][5][15];
 }
 double GetDinucNSyn23_CGAA(){
     return dinucNSyn_stat[1][6][0];
 }
 double GetDinucNSyn23_CGAC(){
     return dinucNSyn_stat[1][6][1];
 }
 double GetDinucNSyn23_CGAG(){
     return dinucNSyn_stat[1][6][2];
 }
 double GetDinucNSyn23_CGAT(){
     return dinucNSyn_stat[1][6][3];
 }
 double GetDinucNSyn23_CGCA(){
     return dinucNSyn_stat[1][6][4];
 }
 double GetDinucNSyn23_CGCC(){
     return dinucNSyn_stat[1][6][5];
 }
 double GetDinucNSyn23_CGCG(){
     return dinucNSyn_stat[1][6][6];
 }
 double GetDinucNSyn23_CGCT(){
     return dinucNSyn_stat[1][6][7];
 }
 double GetDinucNSyn23_CGGA(){
     return dinucNSyn_stat[1][6][8];
 }
 double GetDinucNSyn23_CGGC(){
     return dinucNSyn_stat[1][6][9];
 }
 double GetDinucNSyn23_CGGG(){
     return dinucNSyn_stat[1][6][10];
 }
 double GetDinucNSyn23_CGGT(){
     return dinucNSyn_stat[1][6][11];
 }
 double GetDinucNSyn23_CGTA(){
     return dinucNSyn_stat[1][6][12];
 }
 double GetDinucNSyn23_CGTC(){
     return dinucNSyn_stat[1][6][13];
 }
 double GetDinucNSyn23_CGTG(){
     return dinucNSyn_stat[1][6][14];
 }
 double GetDinucNSyn23_CGTT(){
     return dinucNSyn_stat[1][6][15];
 }
 double GetDinucNSyn23_CTAA(){
     return dinucNSyn_stat[1][7][0];
 }
 double GetDinucNSyn23_CTAC(){
     return dinucNSyn_stat[1][7][1];
 }
 double GetDinucNSyn23_CTAG(){
     return dinucNSyn_stat[1][7][2];
 }
 double GetDinucNSyn23_CTAT(){
     return dinucNSyn_stat[1][7][3];
 }
 double GetDinucNSyn23_CTCA(){
     return dinucNSyn_stat[1][7][4];
 }
 double GetDinucNSyn23_CTCC(){
     return dinucNSyn_stat[1][7][5];
 }
 double GetDinucNSyn23_CTCG(){
     return dinucNSyn_stat[1][7][6];
 }
 double GetDinucNSyn23_CTCT(){
     return dinucNSyn_stat[1][7][7];
 }
 double GetDinucNSyn23_CTGA(){
     return dinucNSyn_stat[1][7][8];
 }
 double GetDinucNSyn23_CTGC(){
     return dinucNSyn_stat[1][7][9];
 }
 double GetDinucNSyn23_CTGG(){
     return dinucNSyn_stat[1][7][10];
 }
 double GetDinucNSyn23_CTGT(){
     return dinucNSyn_stat[1][7][11];
 }
 double GetDinucNSyn23_CTTA(){
     return dinucNSyn_stat[1][7][12];
 }
 double GetDinucNSyn23_CTTC(){
     return dinucNSyn_stat[1][7][13];
 }
 double GetDinucNSyn23_CTTG(){
     return dinucNSyn_stat[1][7][14];
 }
 double GetDinucNSyn23_CTTT(){
     return dinucNSyn_stat[1][7][15];
 }
 double GetDinucNSyn23_GAAA(){
     return dinucNSyn_stat[1][8][0];
 }
 double GetDinucNSyn23_GAAC(){
     return dinucNSyn_stat[1][8][1];
 }
 double GetDinucNSyn23_GAAG(){
     return dinucNSyn_stat[1][8][2];
 }
 double GetDinucNSyn23_GAAT(){
     return dinucNSyn_stat[1][8][3];
 }
 double GetDinucNSyn23_GACA(){
     return dinucNSyn_stat[1][8][4];
 }
 double GetDinucNSyn23_GACC(){
     return dinucNSyn_stat[1][8][5];
 }
 double GetDinucNSyn23_GACG(){
     return dinucNSyn_stat[1][8][6];
 }
 double GetDinucNSyn23_GACT(){
     return dinucNSyn_stat[1][8][7];
 }
 double GetDinucNSyn23_GAGA(){
     return dinucNSyn_stat[1][8][8];
 }
 double GetDinucNSyn23_GAGC(){
     return dinucNSyn_stat[1][8][9];
 }
 double GetDinucNSyn23_GAGG(){
     return dinucNSyn_stat[1][8][10];
 }
 double GetDinucNSyn23_GAGT(){
     return dinucNSyn_stat[1][8][11];
 }
 double GetDinucNSyn23_GATA(){
     return dinucNSyn_stat[1][8][12];
 }
 double GetDinucNSyn23_GATC(){
     return dinucNSyn_stat[1][8][13];
 }
 double GetDinucNSyn23_GATG(){
     return dinucNSyn_stat[1][8][14];
 }
 double GetDinucNSyn23_GATT(){
     return dinucNSyn_stat[1][8][15];
 }
 double GetDinucNSyn23_GCAA(){
     return dinucNSyn_stat[1][9][0];
 }
 double GetDinucNSyn23_GCAC(){
     return dinucNSyn_stat[1][9][1];
 }
 double GetDinucNSyn23_GCAG(){
     return dinucNSyn_stat[1][9][2];
 }
 double GetDinucNSyn23_GCAT(){
     return dinucNSyn_stat[1][9][3];
 }
 double GetDinucNSyn23_GCCA(){
     return dinucNSyn_stat[1][9][4];
 }
 double GetDinucNSyn23_GCCC(){
     return dinucNSyn_stat[1][9][5];
 }
 double GetDinucNSyn23_GCCG(){
     return dinucNSyn_stat[1][9][6];
 }
 double GetDinucNSyn23_GCCT(){
     return dinucNSyn_stat[1][9][7];
 }
 double GetDinucNSyn23_GCGA(){
     return dinucNSyn_stat[1][9][8];
 }
 double GetDinucNSyn23_GCGC(){
     return dinucNSyn_stat[1][9][9];
 }
 double GetDinucNSyn23_GCGG(){
     return dinucNSyn_stat[1][9][10];
 }
 double GetDinucNSyn23_GCGT(){
     return dinucNSyn_stat[1][9][11];
 }
 double GetDinucNSyn23_GCTA(){
     return dinucNSyn_stat[1][9][12];
 }
 double GetDinucNSyn23_GCTC(){
     return dinucNSyn_stat[1][9][13];
 }
 double GetDinucNSyn23_GCTG(){
     return dinucNSyn_stat[1][9][14];
 }
 double GetDinucNSyn23_GCTT(){
     return dinucNSyn_stat[1][9][15];
 }
 double GetDinucNSyn23_GGAA(){
     return dinucNSyn_stat[1][10][0];
 }
 double GetDinucNSyn23_GGAC(){
     return dinucNSyn_stat[1][10][1];
 }
 double GetDinucNSyn23_GGAG(){
     return dinucNSyn_stat[1][10][2];
 }
 double GetDinucNSyn23_GGAT(){
     return dinucNSyn_stat[1][10][3];
 }
 double GetDinucNSyn23_GGCA(){
     return dinucNSyn_stat[1][10][4];
 }
 double GetDinucNSyn23_GGCC(){
     return dinucNSyn_stat[1][10][5];
 }
 double GetDinucNSyn23_GGCG(){
     return dinucNSyn_stat[1][10][6];
 }
 double GetDinucNSyn23_GGCT(){
     return dinucNSyn_stat[1][10][7];
 }
 double GetDinucNSyn23_GGGA(){
     return dinucNSyn_stat[1][10][8];
 }
 double GetDinucNSyn23_GGGC(){
     return dinucNSyn_stat[1][10][9];
 }
 double GetDinucNSyn23_GGGG(){
     return dinucNSyn_stat[1][10][10];
 }
 double GetDinucNSyn23_GGGT(){
     return dinucNSyn_stat[1][10][11];
 }
 double GetDinucNSyn23_GGTA(){
     return dinucNSyn_stat[1][10][12];
 }
 double GetDinucNSyn23_GGTC(){
     return dinucNSyn_stat[1][10][13];
 }
 double GetDinucNSyn23_GGTG(){
     return dinucNSyn_stat[1][10][14];
 }
 double GetDinucNSyn23_GGTT(){
     return dinucNSyn_stat[1][10][15];
 }
 double GetDinucNSyn23_GTAA(){
     return dinucNSyn_stat[1][11][0];
 }
 double GetDinucNSyn23_GTAC(){
     return dinucNSyn_stat[1][11][1];
 }
 double GetDinucNSyn23_GTAG(){
     return dinucNSyn_stat[1][11][2];
 }
 double GetDinucNSyn23_GTAT(){
     return dinucNSyn_stat[1][11][3];
 }
 double GetDinucNSyn23_GTCA(){
     return dinucNSyn_stat[1][11][4];
 }
 double GetDinucNSyn23_GTCC(){
     return dinucNSyn_stat[1][11][5];
 }
 double GetDinucNSyn23_GTCG(){
     return dinucNSyn_stat[1][11][6];
 }
 double GetDinucNSyn23_GTCT(){
     return dinucNSyn_stat[1][11][7];
 }
 double GetDinucNSyn23_GTGA(){
     return dinucNSyn_stat[1][11][8];
 }
 double GetDinucNSyn23_GTGC(){
     return dinucNSyn_stat[1][11][9];
 }
 double GetDinucNSyn23_GTGG(){
     return dinucNSyn_stat[1][11][10];
 }
 double GetDinucNSyn23_GTGT(){
     return dinucNSyn_stat[1][11][11];
 }
 double GetDinucNSyn23_GTTA(){
     return dinucNSyn_stat[1][11][12];
 }
 double GetDinucNSyn23_GTTC(){
     return dinucNSyn_stat[1][11][13];
 }
 double GetDinucNSyn23_GTTG(){
     return dinucNSyn_stat[1][11][14];
 }
 double GetDinucNSyn23_GTTT(){
     return dinucNSyn_stat[1][11][15];
 }
 double GetDinucNSyn23_TAAA(){
     return dinucNSyn_stat[1][12][0];
 }
 double GetDinucNSyn23_TAAC(){
     return dinucNSyn_stat[1][12][1];
 }
 double GetDinucNSyn23_TAAG(){
     return dinucNSyn_stat[1][12][2];
 }
 double GetDinucNSyn23_TAAT(){
     return dinucNSyn_stat[1][12][3];
 }
 double GetDinucNSyn23_TACA(){
     return dinucNSyn_stat[1][12][4];
 }
 double GetDinucNSyn23_TACC(){
     return dinucNSyn_stat[1][12][5];
 }
 double GetDinucNSyn23_TACG(){
     return dinucNSyn_stat[1][12][6];
 }
 double GetDinucNSyn23_TACT(){
     return dinucNSyn_stat[1][12][7];
 }
 double GetDinucNSyn23_TAGA(){
     return dinucNSyn_stat[1][12][8];
 }
 double GetDinucNSyn23_TAGC(){
     return dinucNSyn_stat[1][12][9];
 }
 double GetDinucNSyn23_TAGG(){
     return dinucNSyn_stat[1][12][10];
 }
 double GetDinucNSyn23_TAGT(){
     return dinucNSyn_stat[1][12][11];
 }
 double GetDinucNSyn23_TATA(){
     return dinucNSyn_stat[1][12][12];
 }
 double GetDinucNSyn23_TATC(){
     return dinucNSyn_stat[1][12][13];
 }
 double GetDinucNSyn23_TATG(){
     return dinucNSyn_stat[1][12][14];
 }
 double GetDinucNSyn23_TATT(){
     return dinucNSyn_stat[1][12][15];
 }
 double GetDinucNSyn23_TCAA(){
     return dinucNSyn_stat[1][13][0];
 }
 double GetDinucNSyn23_TCAC(){
     return dinucNSyn_stat[1][13][1];
 }
 double GetDinucNSyn23_TCAG(){
     return dinucNSyn_stat[1][13][2];
 }
 double GetDinucNSyn23_TCAT(){
     return dinucNSyn_stat[1][13][3];
 }
 double GetDinucNSyn23_TCCA(){
     return dinucNSyn_stat[1][13][4];
 }
 double GetDinucNSyn23_TCCC(){
     return dinucNSyn_stat[1][13][5];
 }
 double GetDinucNSyn23_TCCG(){
     return dinucNSyn_stat[1][13][6];
 }
 double GetDinucNSyn23_TCCT(){
     return dinucNSyn_stat[1][13][7];
 }
 double GetDinucNSyn23_TCGA(){
     return dinucNSyn_stat[1][13][8];
 }
 double GetDinucNSyn23_TCGC(){
     return dinucNSyn_stat[1][13][9];
 }
 double GetDinucNSyn23_TCGG(){
     return dinucNSyn_stat[1][13][10];
 }
 double GetDinucNSyn23_TCGT(){
     return dinucNSyn_stat[1][13][11];
 }
 double GetDinucNSyn23_TCTA(){
     return dinucNSyn_stat[1][13][12];
 }
 double GetDinucNSyn23_TCTC(){
     return dinucNSyn_stat[1][13][13];
 }
 double GetDinucNSyn23_TCTG(){
     return dinucNSyn_stat[1][13][14];
 }
 double GetDinucNSyn23_TCTT(){
     return dinucNSyn_stat[1][13][15];
 }
 double GetDinucNSyn23_TGAA(){
     return dinucNSyn_stat[1][14][0];
 }
 double GetDinucNSyn23_TGAC(){
     return dinucNSyn_stat[1][14][1];
 }
 double GetDinucNSyn23_TGAG(){
     return dinucNSyn_stat[1][14][2];
 }
 double GetDinucNSyn23_TGAT(){
     return dinucNSyn_stat[1][14][3];
 }
 double GetDinucNSyn23_TGCA(){
     return dinucNSyn_stat[1][14][4];
 }
 double GetDinucNSyn23_TGCC(){
     return dinucNSyn_stat[1][14][5];
 }
 double GetDinucNSyn23_TGCG(){
     return dinucNSyn_stat[1][14][6];
 }
 double GetDinucNSyn23_TGCT(){
     return dinucNSyn_stat[1][14][7];
 }
 double GetDinucNSyn23_TGGA(){
     return dinucNSyn_stat[1][14][8];
 }
 double GetDinucNSyn23_TGGC(){
     return dinucNSyn_stat[1][14][9];
 }
 double GetDinucNSyn23_TGGG(){
     return dinucNSyn_stat[1][14][10];
 }
 double GetDinucNSyn23_TGGT(){
     return dinucNSyn_stat[1][14][11];
 }
 double GetDinucNSyn23_TGTA(){
     return dinucNSyn_stat[1][14][12];
 }
 double GetDinucNSyn23_TGTC(){
     return dinucNSyn_stat[1][14][13];
 }
 double GetDinucNSyn23_TGTG(){
     return dinucNSyn_stat[1][14][14];
 }
 double GetDinucNSyn23_TGTT(){
     return dinucNSyn_stat[1][14][15];
 }
 double GetDinucNSyn23_TTAA(){
     return dinucNSyn_stat[1][15][0];
 }
 double GetDinucNSyn23_TTAC(){
     return dinucNSyn_stat[1][15][1];
 }
 double GetDinucNSyn23_TTAG(){
     return dinucNSyn_stat[1][15][2];
 }
 double GetDinucNSyn23_TTAT(){
     return dinucNSyn_stat[1][15][3];
 }
 double GetDinucNSyn23_TTCA(){
     return dinucNSyn_stat[1][15][4];
 }
 double GetDinucNSyn23_TTCC(){
     return dinucNSyn_stat[1][15][5];
 }
 double GetDinucNSyn23_TTCG(){
     return dinucNSyn_stat[1][15][6];
 }
 double GetDinucNSyn23_TTCT(){
     return dinucNSyn_stat[1][15][7];
 }
 double GetDinucNSyn23_TTGA(){
     return dinucNSyn_stat[1][15][8];
 }
 double GetDinucNSyn23_TTGC(){
     return dinucNSyn_stat[1][15][9];
 }
 double GetDinucNSyn23_TTGG(){
     return dinucNSyn_stat[1][15][10];
 }
 double GetDinucNSyn23_TTGT(){
     return dinucNSyn_stat[1][15][11];
 }
 double GetDinucNSyn23_TTTA(){
     return dinucNSyn_stat[1][15][12];
 }
 double GetDinucNSyn23_TTTC(){
     return dinucNSyn_stat[1][15][13];
 }
 double GetDinucNSyn23_TTTG(){
     return dinucNSyn_stat[1][15][14];
 }
 double GetDinucNSyn23_TTTT(){
     return dinucNSyn_stat[1][15][15];
 }

 double GetDinucNSyn31_AAAA(){
     return dinucNSyn_stat[2][0][0];
  }
  double GetDinucNSyn31_AAAC(){
     return dinucNSyn_stat[2][0][1];
  }
  double GetDinucNSyn31_AAAG(){
     return dinucNSyn_stat[2][0][2];
  }
  double GetDinucNSyn31_AAAT(){
     return dinucNSyn_stat[2][0][3];
  }
  double GetDinucNSyn31_AACA(){
     return dinucNSyn_stat[2][0][4];
  }
  double GetDinucNSyn31_AACC(){
     return dinucNSyn_stat[2][0][5];
  }
  double GetDinucNSyn31_AACG(){
     return dinucNSyn_stat[2][0][6];
  }
  double GetDinucNSyn31_AACT(){
     return dinucNSyn_stat[2][0][7];
  }
  double GetDinucNSyn31_AAGA(){
     return dinucNSyn_stat[2][0][8];
  }
  double GetDinucNSyn31_AAGC(){
     return dinucNSyn_stat[2][0][9];
  }
  double GetDinucNSyn31_AAGG(){
     return dinucNSyn_stat[2][0][10];
  }
  double GetDinucNSyn31_AAGT(){
     return dinucNSyn_stat[2][0][11];
  }
  double GetDinucNSyn31_AATA(){
     return dinucNSyn_stat[2][0][12];
  }
  double GetDinucNSyn31_AATC(){
     return dinucNSyn_stat[2][0][13];
  }
  double GetDinucNSyn31_AATG(){
     return dinucNSyn_stat[2][0][14];
  }
  double GetDinucNSyn31_AATT(){
     return dinucNSyn_stat[2][0][15];
  }
  double GetDinucNSyn31_ACAA(){
     return dinucNSyn_stat[2][1][0];
  }
  double GetDinucNSyn31_ACAC(){
     return dinucNSyn_stat[2][1][1];
  }
  double GetDinucNSyn31_ACAG(){
     return dinucNSyn_stat[2][1][2];
  }
  double GetDinucNSyn31_ACAT(){
     return dinucNSyn_stat[2][1][3];
  }
  double GetDinucNSyn31_ACCA(){
     return dinucNSyn_stat[2][1][4];
  }
  double GetDinucNSyn31_ACCC(){
     return dinucNSyn_stat[2][1][5];
  }
  double GetDinucNSyn31_ACCG(){
     return dinucNSyn_stat[2][1][6];
  }
  double GetDinucNSyn31_ACCT(){
     return dinucNSyn_stat[2][1][7];
  }
  double GetDinucNSyn31_ACGA(){
     return dinucNSyn_stat[2][1][8];
  }
  double GetDinucNSyn31_ACGC(){
     return dinucNSyn_stat[2][1][9];
  }
  double GetDinucNSyn31_ACGG(){
     return dinucNSyn_stat[2][1][10];
  }
  double GetDinucNSyn31_ACGT(){
     return dinucNSyn_stat[2][1][11];
  }
  double GetDinucNSyn31_ACTA(){
     return dinucNSyn_stat[2][1][12];
  }
  double GetDinucNSyn31_ACTC(){
     return dinucNSyn_stat[2][1][13];
  }
  double GetDinucNSyn31_ACTG(){
     return dinucNSyn_stat[2][1][14];
  }
  double GetDinucNSyn31_ACTT(){
     return dinucNSyn_stat[2][1][15];
  }
  double GetDinucNSyn31_AGAA(){
     return dinucNSyn_stat[2][2][0];
  }
  double GetDinucNSyn31_AGAC(){
     return dinucNSyn_stat[2][2][1];
  }
  double GetDinucNSyn31_AGAG(){
     return dinucNSyn_stat[2][2][2];
  }
  double GetDinucNSyn31_AGAT(){
     return dinucNSyn_stat[2][2][3];
  }
  double GetDinucNSyn31_AGCA(){
     return dinucNSyn_stat[2][2][4];
  }
  double GetDinucNSyn31_AGCC(){
     return dinucNSyn_stat[2][2][5];
  }
  double GetDinucNSyn31_AGCG(){
     return dinucNSyn_stat[2][2][6];
  }
  double GetDinucNSyn31_AGCT(){
     return dinucNSyn_stat[2][2][7];
  }
  double GetDinucNSyn31_AGGA(){
     return dinucNSyn_stat[2][2][8];
  }
  double GetDinucNSyn31_AGGC(){
     return dinucNSyn_stat[2][2][9];
  }
  double GetDinucNSyn31_AGGG(){
     return dinucNSyn_stat[2][2][10];
  }
  double GetDinucNSyn31_AGGT(){
     return dinucNSyn_stat[2][2][11];
  }
  double GetDinucNSyn31_AGTA(){
     return dinucNSyn_stat[2][2][12];
  }
  double GetDinucNSyn31_AGTC(){
     return dinucNSyn_stat[2][2][13];
  }
  double GetDinucNSyn31_AGTG(){
     return dinucNSyn_stat[2][2][14];
  }
  double GetDinucNSyn31_AGTT(){
     return dinucNSyn_stat[2][2][15];
  }
  double GetDinucNSyn31_ATAA(){
     return dinucNSyn_stat[2][3][0];
  }
  double GetDinucNSyn31_ATAC(){
     return dinucNSyn_stat[2][3][1];
  }
  double GetDinucNSyn31_ATAG(){
     return dinucNSyn_stat[2][3][2];
  }
  double GetDinucNSyn31_ATAT(){
     return dinucNSyn_stat[2][3][3];
  }
  double GetDinucNSyn31_ATCA(){
     return dinucNSyn_stat[2][3][4];
  }
  double GetDinucNSyn31_ATCC(){
     return dinucNSyn_stat[2][3][5];
  }
  double GetDinucNSyn31_ATCG(){
     return dinucNSyn_stat[2][3][6];
  }
  double GetDinucNSyn31_ATCT(){
     return dinucNSyn_stat[2][3][7];
  }
  double GetDinucNSyn31_ATGA(){
     return dinucNSyn_stat[2][3][8];
  }
  double GetDinucNSyn31_ATGC(){
     return dinucNSyn_stat[2][3][9];
  }
  double GetDinucNSyn31_ATGG(){
     return dinucNSyn_stat[2][3][10];
  }
  double GetDinucNSyn31_ATGT(){
     return dinucNSyn_stat[2][3][11];
  }
  double GetDinucNSyn31_ATTA(){
     return dinucNSyn_stat[2][3][12];
  }
  double GetDinucNSyn31_ATTC(){
     return dinucNSyn_stat[2][3][13];
  }
  double GetDinucNSyn31_ATTG(){
     return dinucNSyn_stat[2][3][14];
  }
  double GetDinucNSyn31_ATTT(){
     return dinucNSyn_stat[2][3][15];
  }
  double GetDinucNSyn31_CAAA(){
     return dinucNSyn_stat[2][4][0];
  }
  double GetDinucNSyn31_CAAC(){
     return dinucNSyn_stat[2][4][1];
  }
  double GetDinucNSyn31_CAAG(){
     return dinucNSyn_stat[2][4][2];
  }
  double GetDinucNSyn31_CAAT(){
     return dinucNSyn_stat[2][4][3];
  }
  double GetDinucNSyn31_CACA(){
     return dinucNSyn_stat[2][4][4];
  }
  double GetDinucNSyn31_CACC(){
     return dinucNSyn_stat[2][4][5];
  }
  double GetDinucNSyn31_CACG(){
     return dinucNSyn_stat[2][4][6];
  }
  double GetDinucNSyn31_CACT(){
     return dinucNSyn_stat[2][4][7];
  }
  double GetDinucNSyn31_CAGA(){
     return dinucNSyn_stat[2][4][8];
  }
  double GetDinucNSyn31_CAGC(){
     return dinucNSyn_stat[2][4][9];
  }
  double GetDinucNSyn31_CAGG(){
     return dinucNSyn_stat[2][4][10];
  }
  double GetDinucNSyn31_CAGT(){
     return dinucNSyn_stat[2][4][11];
  }
  double GetDinucNSyn31_CATA(){
     return dinucNSyn_stat[2][4][12];
  }
  double GetDinucNSyn31_CATC(){
     return dinucNSyn_stat[2][4][13];
  }
  double GetDinucNSyn31_CATG(){
     return dinucNSyn_stat[2][4][14];
  }
  double GetDinucNSyn31_CATT(){
     return dinucNSyn_stat[2][4][15];
  }
  double GetDinucNSyn31_CCAA(){
     return dinucNSyn_stat[2][5][0];
  }
  double GetDinucNSyn31_CCAC(){
     return dinucNSyn_stat[2][5][1];
  }
  double GetDinucNSyn31_CCAG(){
     return dinucNSyn_stat[2][5][2];
  }
  double GetDinucNSyn31_CCAT(){
     return dinucNSyn_stat[2][5][3];
  }
  double GetDinucNSyn31_CCCA(){
     return dinucNSyn_stat[2][5][4];
  }
  double GetDinucNSyn31_CCCC(){
     return dinucNSyn_stat[2][5][5];
  }
  double GetDinucNSyn31_CCCG(){
     return dinucNSyn_stat[2][5][6];
  }
  double GetDinucNSyn31_CCCT(){
     return dinucNSyn_stat[2][5][7];
  }
  double GetDinucNSyn31_CCGA(){
     return dinucNSyn_stat[2][5][8];
  }
  double GetDinucNSyn31_CCGC(){
     return dinucNSyn_stat[2][5][9];
  }
  double GetDinucNSyn31_CCGG(){
     return dinucNSyn_stat[2][5][10];
  }
  double GetDinucNSyn31_CCGT(){
     return dinucNSyn_stat[2][5][11];
  }
  double GetDinucNSyn31_CCTA(){
     return dinucNSyn_stat[2][5][12];
  }
  double GetDinucNSyn31_CCTC(){
     return dinucNSyn_stat[2][5][13];
  }
  double GetDinucNSyn31_CCTG(){
     return dinucNSyn_stat[2][5][14];
  }
  double GetDinucNSyn31_CCTT(){
     return dinucNSyn_stat[2][5][15];
  }
  double GetDinucNSyn31_CGAA(){
     return dinucNSyn_stat[2][6][0];
  }
  double GetDinucNSyn31_CGAC(){
     return dinucNSyn_stat[2][6][1];
  }
  double GetDinucNSyn31_CGAG(){
     return dinucNSyn_stat[2][6][2];
  }
  double GetDinucNSyn31_CGAT(){
     return dinucNSyn_stat[2][6][3];
  }
  double GetDinucNSyn31_CGCA(){
     return dinucNSyn_stat[2][6][4];
  }
  double GetDinucNSyn31_CGCC(){
     return dinucNSyn_stat[2][6][5];
  }
  double GetDinucNSyn31_CGCG(){
     return dinucNSyn_stat[2][6][6];
  }
  double GetDinucNSyn31_CGCT(){
     return dinucNSyn_stat[2][6][7];
  }
  double GetDinucNSyn31_CGGA(){
     return dinucNSyn_stat[2][6][8];
  }
  double GetDinucNSyn31_CGGC(){
     return dinucNSyn_stat[2][6][9];
  }
  double GetDinucNSyn31_CGGG(){
     return dinucNSyn_stat[2][6][10];
  }
  double GetDinucNSyn31_CGGT(){
     return dinucNSyn_stat[2][6][11];
  }
  double GetDinucNSyn31_CGTA(){
     return dinucNSyn_stat[2][6][12];
  }
  double GetDinucNSyn31_CGTC(){
     return dinucNSyn_stat[2][6][13];
  }
  double GetDinucNSyn31_CGTG(){
     return dinucNSyn_stat[2][6][14];
  }
  double GetDinucNSyn31_CGTT(){
     return dinucNSyn_stat[2][6][15];
  }
  double GetDinucNSyn31_CTAA(){
     return dinucNSyn_stat[2][7][0];
  }
  double GetDinucNSyn31_CTAC(){
     return dinucNSyn_stat[2][7][1];
  }
  double GetDinucNSyn31_CTAG(){
     return dinucNSyn_stat[2][7][2];
  }
  double GetDinucNSyn31_CTAT(){
     return dinucNSyn_stat[2][7][3];
  }
  double GetDinucNSyn31_CTCA(){
     return dinucNSyn_stat[2][7][4];
  }
  double GetDinucNSyn31_CTCC(){
     return dinucNSyn_stat[2][7][5];
  }
  double GetDinucNSyn31_CTCG(){
     return dinucNSyn_stat[2][7][6];
  }
  double GetDinucNSyn31_CTCT(){
     return dinucNSyn_stat[2][7][7];
  }
  double GetDinucNSyn31_CTGA(){
     return dinucNSyn_stat[2][7][8];
  }
  double GetDinucNSyn31_CTGC(){
     return dinucNSyn_stat[2][7][9];
  }
  double GetDinucNSyn31_CTGG(){
     return dinucNSyn_stat[2][7][10];
  }
  double GetDinucNSyn31_CTGT(){
     return dinucNSyn_stat[2][7][11];
  }
  double GetDinucNSyn31_CTTA(){
     return dinucNSyn_stat[2][7][12];
  }
  double GetDinucNSyn31_CTTC(){
     return dinucNSyn_stat[2][7][13];
  }
  double GetDinucNSyn31_CTTG(){
     return dinucNSyn_stat[2][7][14];
  }
  double GetDinucNSyn31_CTTT(){
     return dinucNSyn_stat[2][7][15];
  }
  double GetDinucNSyn31_GAAA(){
     return dinucNSyn_stat[2][8][0];
  }
  double GetDinucNSyn31_GAAC(){
     return dinucNSyn_stat[2][8][1];
  }
  double GetDinucNSyn31_GAAG(){
     return dinucNSyn_stat[2][8][2];
  }
  double GetDinucNSyn31_GAAT(){
     return dinucNSyn_stat[2][8][3];
  }
  double GetDinucNSyn31_GACA(){
     return dinucNSyn_stat[2][8][4];
  }
  double GetDinucNSyn31_GACC(){
     return dinucNSyn_stat[2][8][5];
  }
  double GetDinucNSyn31_GACG(){
     return dinucNSyn_stat[2][8][6];
  }
  double GetDinucNSyn31_GACT(){
     return dinucNSyn_stat[2][8][7];
  }
  double GetDinucNSyn31_GAGA(){
     return dinucNSyn_stat[2][8][8];
  }
  double GetDinucNSyn31_GAGC(){
     return dinucNSyn_stat[2][8][9];
  }
  double GetDinucNSyn31_GAGG(){
     return dinucNSyn_stat[2][8][10];
  }
  double GetDinucNSyn31_GAGT(){
     return dinucNSyn_stat[2][8][11];
  }
  double GetDinucNSyn31_GATA(){
     return dinucNSyn_stat[2][8][12];
  }
  double GetDinucNSyn31_GATC(){
     return dinucNSyn_stat[2][8][13];
  }
  double GetDinucNSyn31_GATG(){
     return dinucNSyn_stat[2][8][14];
  }
  double GetDinucNSyn31_GATT(){
     return dinucNSyn_stat[2][8][15];
  }
  double GetDinucNSyn31_GCAA(){
     return dinucNSyn_stat[2][9][0];
  }
  double GetDinucNSyn31_GCAC(){
     return dinucNSyn_stat[2][9][1];
  }
  double GetDinucNSyn31_GCAG(){
     return dinucNSyn_stat[2][9][2];
  }
  double GetDinucNSyn31_GCAT(){
     return dinucNSyn_stat[2][9][3];
  }
  double GetDinucNSyn31_GCCA(){
     return dinucNSyn_stat[2][9][4];
  }
  double GetDinucNSyn31_GCCC(){
     return dinucNSyn_stat[2][9][5];
  }
  double GetDinucNSyn31_GCCG(){
     return dinucNSyn_stat[2][9][6];
  }
  double GetDinucNSyn31_GCCT(){
     return dinucNSyn_stat[2][9][7];
  }
  double GetDinucNSyn31_GCGA(){
     return dinucNSyn_stat[2][9][8];
  }
  double GetDinucNSyn31_GCGC(){
     return dinucNSyn_stat[2][9][9];
  }
  double GetDinucNSyn31_GCGG(){
     return dinucNSyn_stat[2][9][10];
  }
  double GetDinucNSyn31_GCGT(){
     return dinucNSyn_stat[2][9][11];
  }
  double GetDinucNSyn31_GCTA(){
     return dinucNSyn_stat[2][9][12];
  }
  double GetDinucNSyn31_GCTC(){
     return dinucNSyn_stat[2][9][13];
  }
  double GetDinucNSyn31_GCTG(){
     return dinucNSyn_stat[2][9][14];
  }
  double GetDinucNSyn31_GCTT(){
     return dinucNSyn_stat[2][9][15];
  }
  double GetDinucNSyn31_GGAA(){
     return dinucNSyn_stat[2][10][0];
  }
  double GetDinucNSyn31_GGAC(){
     return dinucNSyn_stat[2][10][1];
  }
  double GetDinucNSyn31_GGAG(){
     return dinucNSyn_stat[2][10][2];
  }
  double GetDinucNSyn31_GGAT(){
     return dinucNSyn_stat[2][10][3];
  }
  double GetDinucNSyn31_GGCA(){
     return dinucNSyn_stat[2][10][4];
  }
  double GetDinucNSyn31_GGCC(){
     return dinucNSyn_stat[2][10][5];
  }
  double GetDinucNSyn31_GGCG(){
     return dinucNSyn_stat[2][10][6];
  }
  double GetDinucNSyn31_GGCT(){
     return dinucNSyn_stat[2][10][7];
  }
  double GetDinucNSyn31_GGGA(){
     return dinucNSyn_stat[2][10][8];
  }
  double GetDinucNSyn31_GGGC(){
     return dinucNSyn_stat[2][10][9];
  }
  double GetDinucNSyn31_GGGG(){
     return dinucNSyn_stat[2][10][10];
  }
  double GetDinucNSyn31_GGGT(){
     return dinucNSyn_stat[2][10][11];
  }
  double GetDinucNSyn31_GGTA(){
     return dinucNSyn_stat[2][10][12];
  }
  double GetDinucNSyn31_GGTC(){
     return dinucNSyn_stat[2][10][13];
  }
  double GetDinucNSyn31_GGTG(){
     return dinucNSyn_stat[2][10][14];
  }
  double GetDinucNSyn31_GGTT(){
     return dinucNSyn_stat[2][10][15];
  }
  double GetDinucNSyn31_GTAA(){
     return dinucNSyn_stat[2][11][0];
  }
  double GetDinucNSyn31_GTAC(){
     return dinucNSyn_stat[2][11][1];
  }
  double GetDinucNSyn31_GTAG(){
     return dinucNSyn_stat[2][11][2];
  }
  double GetDinucNSyn31_GTAT(){
     return dinucNSyn_stat[2][11][3];
  }
  double GetDinucNSyn31_GTCA(){
     return dinucNSyn_stat[2][11][4];
  }
  double GetDinucNSyn31_GTCC(){
     return dinucNSyn_stat[2][11][5];
  }
  double GetDinucNSyn31_GTCG(){
     return dinucNSyn_stat[2][11][6];
  }
  double GetDinucNSyn31_GTCT(){
     return dinucNSyn_stat[2][11][7];
  }
  double GetDinucNSyn31_GTGA(){
     return dinucNSyn_stat[2][11][8];
  }
  double GetDinucNSyn31_GTGC(){
     return dinucNSyn_stat[2][11][9];
  }
  double GetDinucNSyn31_GTGG(){
     return dinucNSyn_stat[2][11][10];
  }
  double GetDinucNSyn31_GTGT(){
     return dinucNSyn_stat[2][11][11];
  }
  double GetDinucNSyn31_GTTA(){
     return dinucNSyn_stat[2][11][12];
  }
  double GetDinucNSyn31_GTTC(){
     return dinucNSyn_stat[2][11][13];
  }
  double GetDinucNSyn31_GTTG(){
     return dinucNSyn_stat[2][11][14];
  }
  double GetDinucNSyn31_GTTT(){
     return dinucNSyn_stat[2][11][15];
  }
  double GetDinucNSyn31_TAAA(){
     return dinucNSyn_stat[2][12][0];
  }
  double GetDinucNSyn31_TAAC(){
     return dinucNSyn_stat[2][12][1];
  }
  double GetDinucNSyn31_TAAG(){
     return dinucNSyn_stat[2][12][2];
  }
  double GetDinucNSyn31_TAAT(){
     return dinucNSyn_stat[2][12][3];
  }
  double GetDinucNSyn31_TACA(){
     return dinucNSyn_stat[2][12][4];
  }
  double GetDinucNSyn31_TACC(){
     return dinucNSyn_stat[2][12][5];
  }
  double GetDinucNSyn31_TACG(){
     return dinucNSyn_stat[2][12][6];
  }
  double GetDinucNSyn31_TACT(){
     return dinucNSyn_stat[2][12][7];
  }
  double GetDinucNSyn31_TAGA(){
     return dinucNSyn_stat[2][12][8];
  }
  double GetDinucNSyn31_TAGC(){
     return dinucNSyn_stat[2][12][9];
  }
  double GetDinucNSyn31_TAGG(){
     return dinucNSyn_stat[2][12][10];
  }
  double GetDinucNSyn31_TAGT(){
     return dinucNSyn_stat[2][12][11];
  }
  double GetDinucNSyn31_TATA(){
     return dinucNSyn_stat[2][12][12];
  }
  double GetDinucNSyn31_TATC(){
     return dinucNSyn_stat[2][12][13];
  }
  double GetDinucNSyn31_TATG(){
     return dinucNSyn_stat[2][12][14];
  }
  double GetDinucNSyn31_TATT(){
     return dinucNSyn_stat[2][12][15];
  }
  double GetDinucNSyn31_TCAA(){
     return dinucNSyn_stat[2][13][0];
  }
  double GetDinucNSyn31_TCAC(){
     return dinucNSyn_stat[2][13][1];
  }
  double GetDinucNSyn31_TCAG(){
     return dinucNSyn_stat[2][13][2];
  }
  double GetDinucNSyn31_TCAT(){
     return dinucNSyn_stat[2][13][3];
  }
  double GetDinucNSyn31_TCCA(){
     return dinucNSyn_stat[2][13][4];
  }
  double GetDinucNSyn31_TCCC(){
     return dinucNSyn_stat[2][13][5];
  }
  double GetDinucNSyn31_TCCG(){
     return dinucNSyn_stat[2][13][6];
  }
  double GetDinucNSyn31_TCCT(){
     return dinucNSyn_stat[2][13][7];
  }
  double GetDinucNSyn31_TCGA(){
     return dinucNSyn_stat[2][13][8];
  }
  double GetDinucNSyn31_TCGC(){
     return dinucNSyn_stat[2][13][9];
  }
  double GetDinucNSyn31_TCGG(){
     return dinucNSyn_stat[2][13][10];
  }
  double GetDinucNSyn31_TCGT(){
     return dinucNSyn_stat[2][13][11];
  }
  double GetDinucNSyn31_TCTA(){
     return dinucNSyn_stat[2][13][12];
  }
  double GetDinucNSyn31_TCTC(){
     return dinucNSyn_stat[2][13][13];
  }
  double GetDinucNSyn31_TCTG(){
     return dinucNSyn_stat[2][13][14];
  }
  double GetDinucNSyn31_TCTT(){
     return dinucNSyn_stat[2][13][15];
  }
  double GetDinucNSyn31_TGAA(){
     return dinucNSyn_stat[2][14][0];
  }
  double GetDinucNSyn31_TGAC(){
     return dinucNSyn_stat[2][14][1];
  }
  double GetDinucNSyn31_TGAG(){
     return dinucNSyn_stat[2][14][2];
  }
  double GetDinucNSyn31_TGAT(){
     return dinucNSyn_stat[2][14][3];
  }
  double GetDinucNSyn31_TGCA(){
     return dinucNSyn_stat[2][14][4];
  }
  double GetDinucNSyn31_TGCC(){
     return dinucNSyn_stat[2][14][5];
  }
  double GetDinucNSyn31_TGCG(){
     return dinucNSyn_stat[2][14][6];
  }
  double GetDinucNSyn31_TGCT(){
     return dinucNSyn_stat[2][14][7];
  }
  double GetDinucNSyn31_TGGA(){
     return dinucNSyn_stat[2][14][8];
  }
  double GetDinucNSyn31_TGGC(){
     return dinucNSyn_stat[2][14][9];
  }
  double GetDinucNSyn31_TGGG(){
     return dinucNSyn_stat[2][14][10];
  }
  double GetDinucNSyn31_TGGT(){
     return dinucNSyn_stat[2][14][11];
  }
  double GetDinucNSyn31_TGTA(){
     return dinucNSyn_stat[2][14][12];
  }
  double GetDinucNSyn31_TGTC(){
     return dinucNSyn_stat[2][14][13];
  }
  double GetDinucNSyn31_TGTG(){
     return dinucNSyn_stat[2][14][14];
  }
  double GetDinucNSyn31_TGTT(){
     return dinucNSyn_stat[2][14][15];
  }
  double GetDinucNSyn31_TTAA(){
     return dinucNSyn_stat[2][15][0];
  }
  double GetDinucNSyn31_TTAC(){
     return dinucNSyn_stat[2][15][1];
  }
  double GetDinucNSyn31_TTAG(){
     return dinucNSyn_stat[2][15][2];
  }
  double GetDinucNSyn31_TTAT(){
     return dinucNSyn_stat[2][15][3];
  }
  double GetDinucNSyn31_TTCA(){
     return dinucNSyn_stat[2][15][4];
  }
  double GetDinucNSyn31_TTCC(){
     return dinucNSyn_stat[2][15][5];
  }
  double GetDinucNSyn31_TTCG(){
     return dinucNSyn_stat[2][15][6];
  }
  double GetDinucNSyn31_TTCT(){
     return dinucNSyn_stat[2][15][7];
  }
  double GetDinucNSyn31_TTGA(){
     return dinucNSyn_stat[2][15][8];
  }
  double GetDinucNSyn31_TTGC(){
     return dinucNSyn_stat[2][15][9];
  }
  double GetDinucNSyn31_TTGG(){
     return dinucNSyn_stat[2][15][10];
  }
  double GetDinucNSyn31_TTGT(){
     return dinucNSyn_stat[2][15][11];
  }
  double GetDinucNSyn31_TTTA(){
     return dinucNSyn_stat[2][15][12];
  }
  double GetDinucNSyn31_TTTC(){
     return dinucNSyn_stat[2][15][13];
  }
  double GetDinucNSyn31_TTTG(){
     return dinucNSyn_stat[2][15][14];
  }
  double GetDinucNSyn31_TTTT(){
     return dinucNSyn_stat[2][15][15];
  }

  double GetNsub(){
       return (Nsub/lparam->Nsite_codon);
  }
  double GetNSynsub(){
       return (Nsynsub/lparam->Nsite_codon);
  }


    double GetssNnsynsub0(){
     return ssNnsynsub[0];
    }
    double GetssNnsynsub1(){
     return ssNnsynsub[1];
    }
    double GetssNnsynsub2(){
     return ssNnsynsub[2];
    }
    double GetssNnsynsub3(){
     return ssNnsynsub[3];
    }
    double GetssNnsynsub4(){
     return ssNnsynsub[4];
    }
    double GetssNnsynsub5(){
     return ssNnsynsub[5];
    }
    double GetssNnsynsub6(){
     return ssNnsynsub[6];
    }
    double GetssNnsynsub7(){
     return ssNnsynsub[7];
    }
    double GetssNnsynsub8(){
     return ssNnsynsub[8];
    }
    double GetssNnsynsub9(){
     return ssNnsynsub[9];
    }
    double GetssNnsynsub10(){
     return ssNnsynsub[10];
    }
    double GetssNnsynsub11(){
     return ssNnsynsub[11];
    }
    double GetssNnsynsub12(){
     return ssNnsynsub[12];
    }
    double GetssNnsynsub13(){
     return ssNnsynsub[13];
    }
    double GetssNnsynsub14(){
     return ssNnsynsub[14];
    }
    double GetssNnsynsub15(){
     return ssNnsynsub[15];
    }
    double GetssNnsynsub16(){
     return ssNnsynsub[16];
    }
    double GetssNnsynsub17(){
     return ssNnsynsub[17];
    }
    double GetssNnsynsub18(){
     return ssNnsynsub[18];
    }
    double GetssNnsynsub19(){
     return ssNnsynsub[19];
    }
    double GetssNnsynsub20(){
     return ssNnsynsub[20];
    }
    double GetssNnsynsub21(){
     return ssNnsynsub[21];
    }
    double GetssNnsynsub22(){
     return ssNnsynsub[22];
    }
    double GetssNnsynsub23(){
     return ssNnsynsub[23];
    }
    double GetssNnsynsub24(){
     return ssNnsynsub[24];
    }
    double GetssNnsynsub25(){
     return ssNnsynsub[25];
    }
    double GetssNnsynsub26(){
     return ssNnsynsub[26];
    }
    double GetssNnsynsub27(){
     return ssNnsynsub[27];
    }
    double GetssNnsynsub28(){
     return ssNnsynsub[28];
    }
    double GetssNnsynsub29(){
     return ssNnsynsub[29];
    }
    double GetssNnsynsub30(){
     return ssNnsynsub[30];
    }
    double GetssNnsynsub31(){
     return ssNnsynsub[31];
    }
    double GetssNnsynsub32(){
     return ssNnsynsub[32];
    }
    double GetssNnsynsub33(){
     return ssNnsynsub[33];
    }
    double GetssNnsynsub34(){
     return ssNnsynsub[34];
    }
    double GetssNnsynsub35(){
     return ssNnsynsub[35];
    }
    double GetssNnsynsub36(){
     return ssNnsynsub[36];
    }
    double GetssNnsynsub37(){
     return ssNnsynsub[37];
    }
    double GetssNnsynsub38(){
     return ssNnsynsub[38];
    }
    double GetssNnsynsub39(){
     return ssNnsynsub[39];
    }
    double GetssNnsynsub40(){
     return ssNnsynsub[40];
    }
    double GetssNnsynsub41(){
     return ssNnsynsub[41];
    }
    double GetssNnsynsub42(){
     return ssNnsynsub[42];
    }
    double GetssNnsynsub43(){
     return ssNnsynsub[43];
    }
    double GetssNnsynsub44(){
     return ssNnsynsub[44];
    }
    double GetssNnsynsub45(){
     return ssNnsynsub[45];
    }
    double GetssNnsynsub46(){
     return ssNnsynsub[46];
    }
    double GetssNnsynsub47(){
     return ssNnsynsub[47];
    }
    double GetssNnsynsub48(){
     return ssNnsynsub[48];
    }
    double GetssNnsynsub49(){
     return ssNnsynsub[49];
    }
    double GetssNnsynsub50(){
     return ssNnsynsub[50];
    }
    double GetssNnsynsub51(){
     return ssNnsynsub[51];
    }
    double GetssNnsynsub52(){
     return ssNnsynsub[52];
    }
    double GetssNnsynsub53(){
     return ssNnsynsub[53];
    }
    double GetssNnsynsub54(){
     return ssNnsynsub[54];
    }
    double GetssNnsynsub55(){
     return ssNnsynsub[55];
    }
    double GetssNnsynsub56(){
     return ssNnsynsub[56];
    }
    double GetssNnsynsub57(){
     return ssNnsynsub[57];
    }
    double GetssNnsynsub58(){
     return ssNnsynsub[58];
    }
    double GetssNnsynsub59(){
     return ssNnsynsub[59];
    }
    double GetssNnsynsub60(){
     return ssNnsynsub[60];
    }
    double GetssNnsynsub61(){
     return ssNnsynsub[61];
    }
    double GetssNnsynsub62(){
     return ssNnsynsub[62];
    }
    double GetssNnsynsub63(){
     return ssNnsynsub[63];
    }
    double GetssNnsynsub64(){
     return ssNnsynsub[64];
    }
    double GetssNnsynsub65(){
     return ssNnsynsub[65];
    }
    double GetssNnsynsub66(){
     return ssNnsynsub[66];
    }
    double GetssNnsynsub67(){
     return ssNnsynsub[67];
    }
    double GetssNnsynsub68(){
     return ssNnsynsub[68];
    }
    double GetssNnsynsub69(){
     return ssNnsynsub[69];
    }
    double GetssNnsynsub70(){
     return ssNnsynsub[70];
    }
    double GetssNnsynsub71(){
     return ssNnsynsub[71];
    }
    double GetssNnsynsub72(){
     return ssNnsynsub[72];
    }
    double GetssNnsynsub73(){
     return ssNnsynsub[73];
    }
    double GetssNnsynsub74(){
     return ssNnsynsub[74];
    }
    double GetssNnsynsub75(){
     return ssNnsynsub[75];
    }
    double GetssNnsynsub76(){
     return ssNnsynsub[76];
    }
    double GetssNnsynsub77(){
     return ssNnsynsub[77];
    }
    double GetssNnsynsub78(){
     return ssNnsynsub[78];
    }
    double GetssNnsynsub79(){
     return ssNnsynsub[79];
    }
    double GetssNnsynsub80(){
     return ssNnsynsub[80];
    }
    double GetssNnsynsub81(){
     return ssNnsynsub[81];
    }
    double GetssNnsynsub82(){
     return ssNnsynsub[82];
    }
    double GetssNnsynsub83(){
     return ssNnsynsub[83];
    }
    double GetssNnsynsub84(){
     return ssNnsynsub[84];
    }
    double GetssNnsynsub85(){
     return ssNnsynsub[85];
    }
    double GetssNnsynsub86(){
     return ssNnsynsub[86];
    }
    double GetssNnsynsub87(){
     return ssNnsynsub[87];
    }
    double GetssNnsynsub88(){
     return ssNnsynsub[88];
    }
    double GetssNnsynsub89(){
     return ssNnsynsub[89];
    }
    double GetssNnsynsub90(){
     return ssNnsynsub[90];
    }
    double GetssNnsynsub91(){
     return ssNnsynsub[91];
    }
    double GetssNnsynsub92(){
     return ssNnsynsub[92];
    }
    double GetssNnsynsub93(){
     return ssNnsynsub[93];
    }
    double GetssNnsynsub94(){
     return ssNnsynsub[94];
    }
    double GetssNnsynsub95(){
     return ssNnsynsub[95];
    }
    double GetssNnsynsub96(){
     return ssNnsynsub[96];
    }
    double GetssNnsynsub97(){
     return ssNnsynsub[97];
    }
    double GetssNnsynsub98(){
     return ssNnsynsub[98];
    }
    double GetssNnsynsub99(){
     return ssNnsynsub[99];
    }
    double GetssNnsynsub100(){
     return ssNnsynsub[100];
    }
    double GetssNnsynsub101(){
     return ssNnsynsub[101];
    }
    double GetssNnsynsub102(){
     return ssNnsynsub[102];
    }
    double GetssNnsynsub103(){
     return ssNnsynsub[103];
    }
    double GetssNnsynsub104(){
     return ssNnsynsub[104];
    }
    double GetssNnsynsub105(){
     return ssNnsynsub[105];
    }
    double GetssNnsynsub106(){
     return ssNnsynsub[106];
    }
    double GetssNnsynsub107(){
     return ssNnsynsub[107];
    }
    double GetssNnsynsub108(){
     return ssNnsynsub[108];
    }
    double GetssNnsynsub109(){
     return ssNnsynsub[109];
    }
    double GetssNnsynsub110(){
     return ssNnsynsub[110];
    }
    double GetssNnsynsub111(){
     return ssNnsynsub[111];
    }
    double GetssNnsynsub112(){
     return ssNnsynsub[112];
    }
    double GetssNnsynsub113(){
     return ssNnsynsub[113];
    }
    double GetssNnsynsub114(){
     return ssNnsynsub[114];
    }
    double GetssNnsynsub115(){
     return ssNnsynsub[115];
    }
    double GetssNnsynsub116(){
     return ssNnsynsub[116];
    }
    double GetssNnsynsub117(){
     return ssNnsynsub[117];
    }
    double GetssNnsynsub118(){
     return ssNnsynsub[118];
    }
    double GetssNnsynsub119(){
     return ssNnsynsub[119];
    }
    double GetssNnsynsub120(){
     return ssNnsynsub[120];
    }
    double GetssNnsynsub121(){
     return ssNnsynsub[121];
    }
    double GetssNnsynsub122(){
     return ssNnsynsub[122];
    }
    double GetssNnsynsub123(){
     return ssNnsynsub[123];
    }
    double GetssNnsynsub124(){
     return ssNnsynsub[124];
    }
    double GetssNnsynsub125(){
     return ssNnsynsub[125];
    }
    double GetssNnsynsub126(){
     return ssNnsynsub[126];
    }
    double GetssNnsynsub127(){
     return ssNnsynsub[127];
    }
    double GetssNnsynsub128(){
     return ssNnsynsub[128];
    }
    double GetssNnsynsub129(){
     return ssNnsynsub[129];
    }
    double GetssNnsynsub130(){
     return ssNnsynsub[130];
    }
    double GetssNnsynsub131(){
     return ssNnsynsub[131];
    }
    double GetssNnsynsub132(){
     return ssNnsynsub[132];
    }
    double GetssNnsynsub133(){
     return ssNnsynsub[133];
    }
    double GetssNnsynsub134(){
     return ssNnsynsub[134];
    }
    double GetssNnsynsub135(){
     return ssNnsynsub[135];
    }
    double GetssNnsynsub136(){
     return ssNnsynsub[136];
    }
    double GetssNnsynsub137(){
     return ssNnsynsub[137];
    }
    double GetssNnsynsub138(){
     return ssNnsynsub[138];
    }
    double GetssNnsynsub139(){
     return ssNnsynsub[139];
    }
    double GetssNnsynsub140(){
     return ssNnsynsub[140];
    }
    double GetssNnsynsub141(){
     return ssNnsynsub[141];
    }
    double GetssNnsynsub142(){
     return ssNnsynsub[142];
    }
    double GetssNnsynsub143(){
     return ssNnsynsub[143];
    }
    double GetssNnsynsub144(){
     return ssNnsynsub[144];
    }
    double GetssNnsynsub145(){
     return ssNnsynsub[145];
    }
    double GetssNnsynsub146(){
     return ssNnsynsub[146];
    }
    double GetssNnsynsub147(){
     return ssNnsynsub[147];
    }
    double GetssNnsynsub148(){
     return ssNnsynsub[148];
    }
    double GetssNnsynsub149(){
     return ssNnsynsub[149];
    }
    double GetssNnsynsub150(){
     return ssNnsynsub[150];
    }
    double GetssNnsynsub151(){
     return ssNnsynsub[151];
    }
    double GetssNnsynsub152(){
     return ssNnsynsub[152];
    }
    double GetssNnsynsub153(){
     return ssNnsynsub[153];
    }
    double GetssNnsynsub154(){
     return ssNnsynsub[154];
    }
    double GetssNnsynsub155(){
     return ssNnsynsub[155];
    }
    double GetssNnsynsub156(){
     return ssNnsynsub[156];
    }
    double GetssNnsynsub157(){
     return ssNnsynsub[157];
    }
    double GetssNnsynsub158(){
     return ssNnsynsub[158];
    }
    double GetssNnsynsub159(){
     return ssNnsynsub[159];
    }
    double GetssNnsynsub160(){
     return ssNnsynsub[160];
    }
    double GetssNnsynsub161(){
     return ssNnsynsub[161];
    }
    double GetssNnsynsub162(){
     return ssNnsynsub[162];
    }
    double GetssNnsynsub163(){
     return ssNnsynsub[163];
    }
    double GetssNnsynsub164(){
     return ssNnsynsub[164];
    }
    double GetssNnsynsub165(){
     return ssNnsynsub[165];
    }
    double GetssNnsynsub166(){
     return ssNnsynsub[166];
    }
    double GetssNnsynsub167(){
     return ssNnsynsub[167];
    }
    double GetssNnsynsub168(){
     return ssNnsynsub[168];
    }
    double GetssNnsynsub169(){
     return ssNnsynsub[169];
    }
    double GetssNnsynsub170(){
     return ssNnsynsub[170];
    }
    double GetssNnsynsub171(){
     return ssNnsynsub[171];
    }
    double GetssNnsynsub172(){
     return ssNnsynsub[172];
    }
    double GetssNnsynsub173(){
     return ssNnsynsub[173];
    }
    double GetssNnsynsub174(){
     return ssNnsynsub[174];
    }
    double GetssNnsynsub175(){
     return ssNnsynsub[175];
    }
    double GetssNnsynsub176(){
     return ssNnsynsub[176];
    }
    double GetssNnsynsub177(){
     return ssNnsynsub[177];
    }
    double GetssNnsynsub178(){
     return ssNnsynsub[178];
    }
    double GetssNnsynsub179(){
     return ssNnsynsub[179];
    }
    double GetssNnsynsub180(){
     return ssNnsynsub[180];
    }
    double GetssNnsynsub181(){
     return ssNnsynsub[181];
    }
    double GetssNnsynsub182(){
     return ssNnsynsub[182];
    }
    double GetssNnsynsub183(){
     return ssNnsynsub[183];
    }
    double GetssNnsynsub184(){
     return ssNnsynsub[184];
    }
    double GetssNnsynsub185(){
     return ssNnsynsub[185];
    }
    double GetssNnsynsub186(){
     return ssNnsynsub[186];
    }
    double GetssNnsynsub187(){
     return ssNnsynsub[187];
    }
    double GetssNnsynsub188(){
     return ssNnsynsub[188];
    }
    double GetssNnsynsub189(){
     return ssNnsynsub[189];
    }
    double GetssNnsynsub190(){
     return ssNnsynsub[190];
    }
    double GetssNnsynsub191(){
     return ssNnsynsub[191];
    }
    double GetssNnsynsub192(){
     return ssNnsynsub[192];
    }
    double GetssNnsynsub193(){
     return ssNnsynsub[193];
    }
    double GetssNnsynsub194(){
     return ssNnsynsub[194];
    }
    double GetssNnsynsub195(){
     return ssNnsynsub[195];
    }
    double GetssNnsynsub196(){
     return ssNnsynsub[196];
    }
    double GetssNnsynsub197(){
     return ssNnsynsub[197];
    }
    double GetssNnsynsub198(){
     return ssNnsynsub[198];
    }
    double GetssNnsynsub199(){
     return ssNnsynsub[199];
    }
    double GetssNnsynsub200(){
     return ssNnsynsub[200];
    }
    double GetssNnsynsub201(){
     return ssNnsynsub[201];
    }
    double GetssNnsynsub202(){
     return ssNnsynsub[202];
    }
    double GetssNnsynsub203(){
     return ssNnsynsub[203];
    }
    double GetssNnsynsub204(){
     return ssNnsynsub[204];
    }
    double GetssNnsynsub205(){
     return ssNnsynsub[205];
    }
    double GetssNnsynsub206(){
     return ssNnsynsub[206];
    }
    double GetssNnsynsub207(){
     return ssNnsynsub[207];
    }
    double GetssNnsynsub208(){
     return ssNnsynsub[208];
    }
    double GetssNnsynsub209(){
     return ssNnsynsub[209];
    }
    double GetssNnsynsub210(){
     return ssNnsynsub[210];
    }
    double GetssNnsynsub211(){
     return ssNnsynsub[211];
    }
    double GetssNnsynsub212(){
     return ssNnsynsub[212];
    }
    double GetssNnsynsub213(){
     return ssNnsynsub[213];
    }
    double GetssNnsynsub214(){
     return ssNnsynsub[214];
    }
    double GetssNnsynsub215(){
     return ssNnsynsub[215];
    }
    double GetssNnsynsub216(){
     return ssNnsynsub[216];
    }
    double GetssNnsynsub217(){
     return ssNnsynsub[217];
    }
    double GetssNnsynsub218(){
     return ssNnsynsub[218];
    }
    double GetssNnsynsub219(){
     return ssNnsynsub[219];
    }
    double GetssNnsynsub220(){
     return ssNnsynsub[220];
    }
    double GetssNnsynsub221(){
     return ssNnsynsub[221];
    }
    double GetssNnsynsub222(){
     return ssNnsynsub[222];
    }
    double GetssNnsynsub223(){
     return ssNnsynsub[223];
    }
    double GetssNnsynsub224(){
     return ssNnsynsub[224];
    }
    double GetssNnsynsub225(){
     return ssNnsynsub[225];
    }
    double GetssNnsynsub226(){
     return ssNnsynsub[226];
    }
    double GetssNnsynsub227(){
     return ssNnsynsub[227];
    }
    double GetssNnsynsub228(){
     return ssNnsynsub[228];
    }
    double GetssNnsynsub229(){
     return ssNnsynsub[229];
    }
    double GetssNnsynsub230(){
     return ssNnsynsub[230];
    }
    double GetssNnsynsub231(){
     return ssNnsynsub[231];
    }
    double GetssNnsynsub232(){
     return ssNnsynsub[232];
    }
    double GetssNnsynsub233(){
     return ssNnsynsub[233];
    }
    double GetssNnsynsub234(){
     return ssNnsynsub[234];
    }
    double GetssNnsynsub235(){
     return ssNnsynsub[235];
    }
    double GetssNnsynsub236(){
     return ssNnsynsub[236];
    }
    double GetssNnsynsub237(){
     return ssNnsynsub[237];
    }
    double GetssNnsynsub238(){
     return ssNnsynsub[238];
    }
    double GetssNnsynsub239(){
     return ssNnsynsub[239];
    }
    double GetssNnsynsub240(){
     return ssNnsynsub[240];
    }
    double GetssNnsynsub241(){
     return ssNnsynsub[241];
    }
    double GetssNnsynsub242(){
     return ssNnsynsub[242];
    }
    double GetssNnsynsub243(){
     return ssNnsynsub[243];
    }
    double GetssNnsynsub244(){
     return ssNnsynsub[244];
    }
    double GetssNnsynsub245(){
     return ssNnsynsub[245];
    }
    double GetssNnsynsub246(){
     return ssNnsynsub[246];
    }
    double GetssNnsynsub247(){
     return ssNnsynsub[247];
    }
    double GetssNnsynsub248(){
     return ssNnsynsub[248];
    }
    double GetssNnsynsub249(){
     return ssNnsynsub[249];
    }
    double GetssNnsynsub250(){
     return ssNnsynsub[250];
    }
    double GetssNnsynsub251(){
     return ssNnsynsub[251];
    }
    double GetssNnsynsub252(){
     return ssNnsynsub[252];
    }
    double GetssNnsynsub253(){
     return ssNnsynsub[253];
    }
    double GetssNnsynsub254(){
     return ssNnsynsub[254];
    }
    double GetssNnsynsub255(){
     return ssNnsynsub[255];
    }
    double GetssNnsynsub256(){
     return ssNnsynsub[256];
    }
    double GetssNnsynsub257(){
     return ssNnsynsub[257];
    }
    double GetssNnsynsub258(){
     return ssNnsynsub[258];
    }
    double GetssNnsynsub259(){
     return ssNnsynsub[259];
    }
    double GetssNnsynsub260(){
     return ssNnsynsub[260];
    }
    double GetssNnsynsub261(){
     return ssNnsynsub[261];
    }
    double GetssNnsynsub262(){
     return ssNnsynsub[262];
    }
    double GetssNnsynsub263(){
     return ssNnsynsub[263];
    }
    double GetssNnsynsub264(){
     return ssNnsynsub[264];
    }
    double GetssNnsynsub265(){
     return ssNnsynsub[265];
    }
    double GetssNnsynsub266(){
     return ssNnsynsub[266];
    }
    double GetssNnsynsub267(){
     return ssNnsynsub[267];
    }
    double GetssNnsynsub268(){
     return ssNnsynsub[268];
    }
    double GetssNnsynsub269(){
     return ssNnsynsub[269];
    }
    double GetssNnsynsub270(){
     return ssNnsynsub[270];
    }
    double GetssNnsynsub271(){
     return ssNnsynsub[271];
    }
    double GetssNnsynsub272(){
     return ssNnsynsub[272];
    }
    double GetssNnsynsub273(){
     return ssNnsynsub[273];
    }
    double GetssNnsynsub274(){
     return ssNnsynsub[274];
    }
    double GetssNnsynsub275(){
     return ssNnsynsub[275];
    }
    double GetssNnsynsub276(){
     return ssNnsynsub[276];
    }
    double GetssNnsynsub277(){
     return ssNnsynsub[277];
    }
    double GetssNnsynsub278(){
     return ssNnsynsub[278];
    }
    double GetssNnsynsub279(){
     return ssNnsynsub[279];
    }
    double GetssNnsynsub280(){
     return ssNnsynsub[280];
    }
    double GetssNnsynsub281(){
     return ssNnsynsub[281];
    }
    double GetssNnsynsub282(){
     return ssNnsynsub[282];
    }
    double GetssNnsynsub283(){
     return ssNnsynsub[283];
    }
    double GetssNnsynsub284(){
     return ssNnsynsub[284];
    }
    double GetssNnsynsub285(){
     return ssNnsynsub[285];
    }
    double GetssNnsynsub286(){
     return ssNnsynsub[286];
    }
    double GetssNnsynsub287(){
     return ssNnsynsub[287];
    }
    double GetssNnsynsub288(){
     return ssNnsynsub[288];
    }
    double GetssNnsynsub289(){
     return ssNnsynsub[289];
    }
    double GetssNnsynsub290(){
     return ssNnsynsub[290];
    }
    double GetssNnsynsub291(){
     return ssNnsynsub[291];
    }
    double GetssNnsynsub292(){
     return ssNnsynsub[292];
    }
    double GetssNnsynsub293(){
     return ssNnsynsub[293];
    }
    double GetssNnsynsub294(){
     return ssNnsynsub[294];
    }
    double GetssNnsynsub295(){
     return ssNnsynsub[295];
    }
    double GetssNnsynsub296(){
     return ssNnsynsub[296];
    }
    double GetssNnsynsub297(){
     return ssNnsynsub[297];
    }
    double GetssNnsynsub298(){
     return ssNnsynsub[298];
    }
    double GetssNnsynsub299(){
     return ssNnsynsub[299];
    }
    double GetssNnsynsub300(){
     return ssNnsynsub[300];
    }
    double GetssNnsynsub301(){
     return ssNnsynsub[301];
    }
    double GetssNnsynsub302(){
     return ssNnsynsub[302];
    }
    double GetssNnsynsub303(){
     return ssNnsynsub[303];
    }
    double GetssNnsynsub304(){
     return ssNnsynsub[304];
    }
    double GetssNnsynsub305(){
     return ssNnsynsub[305];
    }
    double GetssNnsynsub306(){
     return ssNnsynsub[306];
    }
    double GetssNnsynsub307(){
     return ssNnsynsub[307];
    }
    double GetssNnsynsub308(){
     return ssNnsynsub[308];
    }
    double GetssNnsynsub309(){
     return ssNnsynsub[309];
    }
    double GetssNnsynsub310(){
     return ssNnsynsub[310];
    }
    double GetssNnsynsub311(){
     return ssNnsynsub[311];
    }
    double GetssNnsynsub312(){
     return ssNnsynsub[312];
    }
    double GetssNnsynsub313(){
     return ssNnsynsub[313];
    }
    double GetssNnsynsub314(){
     return ssNnsynsub[314];
    }
    double GetssNnsynsub315(){
     return ssNnsynsub[315];
    }
    double GetssNnsynsub316(){
     return ssNnsynsub[316];
    }
    double GetssNnsynsub317(){
     return ssNnsynsub[317];
    }
    double GetssNnsynsub318(){
     return ssNnsynsub[318];
    }
    double GetssNnsynsub319(){
     return ssNnsynsub[319];
    }
    double GetssNnsynsub320(){
     return ssNnsynsub[320];
    }
    double GetssNnsynsub321(){
     return ssNnsynsub[321];
    }
    double GetssNnsynsub322(){
     return ssNnsynsub[322];
    }
    double GetssNnsynsub323(){
     return ssNnsynsub[323];
    }
    double GetssNnsynsub324(){
     return ssNnsynsub[324];
    }
    double GetssNnsynsub325(){
     return ssNnsynsub[325];
    }
    double GetssNnsynsub326(){
     return ssNnsynsub[326];
    }
    double GetssNnsynsub327(){
     return ssNnsynsub[327];
    }
    double GetssNnsynsub328(){
     return ssNnsynsub[328];
    }
    double GetssNnsynsub329(){
     return ssNnsynsub[329];
    }
    double GetssNnsynsub330(){
     return ssNnsynsub[330];
    }
    double GetssNnsynsub331(){
     return ssNnsynsub[331];
    }
    double GetssNnsynsub332(){
     return ssNnsynsub[332];
    }
    double GetssNnsynsub333(){
     return ssNnsynsub[333];
    }
    double GetssNnsynsub334(){
     return ssNnsynsub[334];
    }
    double GetssNnsynsub335(){
     return ssNnsynsub[335];
    }
    double GetssNnsynsub336(){
     return ssNnsynsub[336];
    }
    double GetssNnsynsub337(){
     return ssNnsynsub[337];
    }
    double GetssNnsynsub338(){
     return ssNnsynsub[338];
    }
    double GetssNnsynsub339(){
     return ssNnsynsub[339];
    }
    double GetssNnsynsub340(){
     return ssNnsynsub[340];
    }
    double GetssNnsynsub341(){
     return ssNnsynsub[341];
    }
    double GetssNnsynsub342(){
     return ssNnsynsub[342];
    }
    double GetssNnsynsub343(){
     return ssNnsynsub[343];
    }
    double GetssNnsynsub344(){
     return ssNnsynsub[344];
    }
    double GetssNnsynsub345(){
     return ssNnsynsub[345];
    }
    double GetssNnsynsub346(){
     return ssNnsynsub[346];
    }
    double GetssNnsynsub347(){
     return ssNnsynsub[347];
    }
    double GetssNnsynsub348(){
     return ssNnsynsub[348];
    }
    double GetssNnsynsub349(){
     return ssNnsynsub[349];
    }
    double GetssNnsynsub350(){
     return ssNnsynsub[350];
    }
    double GetssNnsynsub351(){
     return ssNnsynsub[351];
    }
    double GetssNnsynsub352(){
     return ssNnsynsub[352];
    }
    double GetssNnsynsub353(){
     return ssNnsynsub[353];
    }
    double GetssNnsynsub354(){
     return ssNnsynsub[354];
    }
    double GetssNnsynsub355(){
     return ssNnsynsub[355];
    }
    double GetssNnsynsub356(){
     return ssNnsynsub[356];
    }
    double GetssNnsynsub357(){
     return ssNnsynsub[357];
    }
    double GetssNnsynsub358(){
     return ssNnsynsub[358];
    }
    double GetssNnsynsub359(){
     return ssNnsynsub[359];
    }
    double GetssNnsynsub360(){
     return ssNnsynsub[360];
    }
    double GetssNnsynsub361(){
     return ssNnsynsub[361];
    }
    double GetssNnsynsub362(){
     return ssNnsynsub[362];
    }
    double GetssNnsynsub363(){
     return ssNnsynsub[363];
    }
    double GetssNnsynsub364(){
     return ssNnsynsub[364];
    }
    double GetssNnsynsub365(){
     return ssNnsynsub[365];
    }
    double GetssNnsynsub366(){
     return ssNnsynsub[366];
    }
    double GetssNnsynsub367(){
     return ssNnsynsub[367];
    }
    double GetssNnsynsub368(){
     return ssNnsynsub[368];
    }
    double GetssNnsynsub369(){
     return ssNnsynsub[369];
    }
    double GetssNnsynsub370(){
     return ssNnsynsub[370];
    }
    double GetssNnsynsub371(){
     return ssNnsynsub[371];
    }
    double GetssNnsynsub372(){
     return ssNnsynsub[372];
    }
    double GetssNnsynsub373(){
     return ssNnsynsub[373];
    }
    double GetssNnsynsub374(){
     return ssNnsynsub[374];
    }
    double GetssNnsynsub375(){
     return ssNnsynsub[375];
    }
    double GetssNnsynsub376(){
     return ssNnsynsub[376];
    }
    double GetssNnsynsub377(){
     return ssNnsynsub[377];
    }
    double GetssNnsynsub378(){
     return ssNnsynsub[378];
    }
    double GetssNnsynsub379(){
     return ssNnsynsub[379];
    }
    double GetssNnsynsub380(){
     return ssNnsynsub[380];
    }
    double GetssNnsynsub381(){
     return ssNnsynsub[381];
    }
    double GetssNnsynsub382(){
     return ssNnsynsub[382];
    }
    double GetssNnsynsub383(){
     return ssNnsynsub[383];
    }
    double GetssNnsynsub384(){
     return ssNnsynsub[384];
    }
    double GetssNnsynsub385(){
     return ssNnsynsub[385];
    }
    double GetssNnsynsub386(){
     return ssNnsynsub[386];
    }
    double GetssNnsynsub387(){
     return ssNnsynsub[387];
    }
    double GetssNnsynsub388(){
     return ssNnsynsub[388];
    }
    double GetssNnsynsub389(){
     return ssNnsynsub[389];
    }
    double GetssNnsynsub390(){
     return ssNnsynsub[390];
    }
    double GetssNnsynsub391(){
     return ssNnsynsub[391];
    }
    double GetssNnsynsub392(){
     return ssNnsynsub[392];
    }
    double GetssNnsynsub393(){
     return ssNnsynsub[393];
    }
    double GetssNnsynsub394(){
     return ssNnsynsub[394];
    }
    double GetssNnsynsub395(){
     return ssNnsynsub[395];
    }
    double GetssNnsynsub396(){
     return ssNnsynsub[396];
    }
    double GetssNnsynsub397(){
     return ssNnsynsub[397];
    }
    double GetssNnsynsub398(){
     return ssNnsynsub[398];
    }
    double GetssNnsynsub399(){
     return ssNnsynsub[399];
    }
    double GetssNnsynsub400(){
     return ssNnsynsub[400];
    }
    double GetssNnsynsub401(){
     return ssNnsynsub[401];
    }
    double GetssNnsynsub402(){
     return ssNnsynsub[402];
    }
    double GetssNnsynsub403(){
     return ssNnsynsub[403];
    }
    double GetssNnsynsub404(){
     return ssNnsynsub[404];
    }
    double GetssNnsynsub405(){
     return ssNnsynsub[405];
    }
    double GetssNnsynsub406(){
     return ssNnsynsub[406];
    }
    double GetssNnsynsub407(){
     return ssNnsynsub[407];
    }
    double GetssNnsynsub408(){
     return ssNnsynsub[408];
    }
    double GetssNnsynsub409(){
     return ssNnsynsub[409];
    }
    double GetssNnsynsub410(){
     return ssNnsynsub[410];
    }
    double GetssNnsynsub411(){
     return ssNnsynsub[411];
    }
    double GetssNnsynsub412(){
     return ssNnsynsub[412];
    }
    double GetssNnsynsub413(){
     return ssNnsynsub[413];
    }
    double GetssNnsynsub414(){
     return ssNnsynsub[414];
    }
    double GetssNnsynsub415(){
     return ssNnsynsub[415];
    }
    double GetssNnsynsub416(){
     return ssNnsynsub[416];
    }
    double GetssNnsynsub417(){
     return ssNnsynsub[417];
    }
    double GetssNnsynsub418(){
     return ssNnsynsub[418];
    }
    double GetssNnsynsub419(){
     return ssNnsynsub[419];
    }
    double GetssNnsynsub420(){
     return ssNnsynsub[420];
    }
    double GetssNnsynsub421(){
     return ssNnsynsub[421];
    }
    double GetssNnsynsub422(){
     return ssNnsynsub[422];
    }
    double GetssNnsynsub423(){
     return ssNnsynsub[423];
    }
    double GetssNnsynsub424(){
     return ssNnsynsub[424];
    }
    double GetssNnsynsub425(){
     return ssNnsynsub[425];
    }
    double GetssNnsynsub426(){
     return ssNnsynsub[426];
    }
    double GetssNnsynsub427(){
     return ssNnsynsub[427];
    }
    double GetssNnsynsub428(){
     return ssNnsynsub[428];
    }
    double GetssNnsynsub429(){
     return ssNnsynsub[429];
    }
    double GetssNnsynsub430(){
     return ssNnsynsub[430];
    }
    double GetssNnsynsub431(){
     return ssNnsynsub[431];
    }
    double GetssNnsynsub432(){
     return ssNnsynsub[432];
    }
    double GetssNnsynsub433(){
     return ssNnsynsub[433];
    }
    double GetssNnsynsub434(){
     return ssNnsynsub[434];
    }
    double GetssNnsynsub435(){
     return ssNnsynsub[435];
    }
    double GetssNnsynsub436(){
     return ssNnsynsub[436];
    }
    double GetssNnsynsub437(){
     return ssNnsynsub[437];
    }
    double GetssNnsynsub438(){
     return ssNnsynsub[438];
    }
    double GetssNnsynsub439(){
     return ssNnsynsub[439];
    }
    double GetssNnsynsub440(){
     return ssNnsynsub[440];
    }
    double GetssNnsynsub441(){
     return ssNnsynsub[441];
    }
    double GetssNnsynsub442(){
     return ssNnsynsub[442];
    }
    double GetssNnsynsub443(){
     return ssNnsynsub[443];
    }
    double GetssNnsynsub444(){
     return ssNnsynsub[444];
    }
    double GetssNnsynsub445(){
     return ssNnsynsub[445];
    }
    double GetssNnsynsub446(){
     return ssNnsynsub[446];
    }
    double GetssNnsynsub447(){
     return ssNnsynsub[447];
    }
    double GetssNnsynsub448(){
     return ssNnsynsub[448];
    }
    double GetssNnsynsub449(){
     return ssNnsynsub[449];
    }
    double GetssNnsynsub450(){
     return ssNnsynsub[450];
    }
    double GetssNnsynsub451(){
     return ssNnsynsub[451];
    }
    double GetssNnsynsub452(){
     return ssNnsynsub[452];
    }
    double GetssNnsynsub453(){
     return ssNnsynsub[453];
    }
    double GetssNnsynsub454(){
     return ssNnsynsub[454];
    }
    double GetssNnsynsub455(){
     return ssNnsynsub[455];
    }
    double GetssNnsynsub456(){
     return ssNnsynsub[456];
    }
    double GetssNnsynsub457(){
     return ssNnsynsub[457];
    }
    double GetssNnsynsub458(){
     return ssNnsynsub[458];
    }
    double GetssNnsynsub459(){
     return ssNnsynsub[459];
    }
    double GetssNnsynsub460(){
     return ssNnsynsub[460];
    }
    double GetssNnsynsub461(){
     return ssNnsynsub[461];
    }
    double GetssNnsynsub462(){
     return ssNnsynsub[462];
    }
    double GetssNnsynsub463(){
     return ssNnsynsub[463];
    }
    double GetssNnsynsub464(){
     return ssNnsynsub[464];
    }
    double GetssNnsynsub465(){
     return ssNnsynsub[465];
    }
    double GetssNnsynsub466(){
     return ssNnsynsub[466];
    }
    double GetssNnsynsub467(){
     return ssNnsynsub[467];
    }
    double GetssNnsynsub468(){
     return ssNnsynsub[468];
    }
    double GetssNnsynsub469(){
     return ssNnsynsub[469];
    }
    double GetssNnsynsub470(){
     return ssNnsynsub[470];
    }
    double GetssNnsynsub471(){
     return ssNnsynsub[471];
    }
    double GetssNnsynsub472(){
     return ssNnsynsub[472];
    }
    double GetssNnsynsub473(){
     return ssNnsynsub[473];
    }
    double GetssNnsynsub474(){
     return ssNnsynsub[474];
    }
    double GetssNnsynsub475(){
     return ssNnsynsub[475];
    }
    double GetssNnsynsub476(){
     return ssNnsynsub[476];
    }
    double GetssNnsynsub477(){
     return ssNnsynsub[477];
    }
    double GetssNnsynsub478(){
     return ssNnsynsub[478];
    }
    double GetssNnsynsub479(){
     return ssNnsynsub[479];
    }
    double GetssNnsynsub480(){
     return ssNnsynsub[480];
    }
    double GetssNnsynsub481(){
     return ssNnsynsub[481];
    }
    double GetssNnsynsub482(){
     return ssNnsynsub[482];
    }
    double GetssNnsynsub483(){
     return ssNnsynsub[483];
    }
    double GetssNnsynsub484(){
     return ssNnsynsub[484];
    }
    double GetssNnsynsub485(){
     return ssNnsynsub[485];
    }
    double GetssNnsynsub486(){
     return ssNnsynsub[486];
    }
    double GetssNnsynsub487(){
     return ssNnsynsub[487];
    }
    double GetssNnsynsub488(){
     return ssNnsynsub[488];
    }
    double GetssNnsynsub489(){
     return ssNnsynsub[489];
    }
    double GetssNnsynsub490(){
     return ssNnsynsub[490];
    }
    double GetssNnsynsub491(){
     return ssNnsynsub[491];
    }
    double GetssNnsynsub492(){
     return ssNnsynsub[492];
    }
    double GetssNnsynsub493(){
     return ssNnsynsub[493];
    }
    double GetssNnsynsub494(){
     return ssNnsynsub[494];
    }
    double GetssNnsynsub495(){
     return ssNnsynsub[495];
    }
    double GetssNnsynsub496(){
     return ssNnsynsub[496];
    }
    double GetssNnsynsub497(){
     return ssNnsynsub[497];
    }
    double GetssNnsynsub498(){
     return ssNnsynsub[498];
    }
    double GetssNnsynsub499(){
     return ssNnsynsub[499];
    }
    double GetssNnsynsub500(){
     return ssNnsynsub[500];
    }
    double GetssNnsynsub501(){
     return ssNnsynsub[501];
    }
    double GetssNnsynsub502(){
     return ssNnsynsub[502];
    }
    double GetssNnsynsub503(){
     return ssNnsynsub[503];
    }
    double GetssNnsynsub504(){
     return ssNnsynsub[504];
    }
    double GetssNnsynsub505(){
     return ssNnsynsub[505];
    }
    double GetssNnsynsub506(){
     return ssNnsynsub[506];
    }
    double GetssNnsynsub507(){
     return ssNnsynsub[507];
    }
    double GetssNnsynsub508(){
     return ssNnsynsub[508];
    }
    double GetssNnsynsub509(){
     return ssNnsynsub[509];
    }
    double GetssNnsynsub510(){
     return ssNnsynsub[510];
    }
    double GetssNnsynsub511(){
     return ssNnsynsub[511];
    }
    double GetssNnsynsub512(){
     return ssNnsynsub[512];
    }
    double GetssNnsynsub513(){
     return ssNnsynsub[513];
    }
    double GetssNnsynsub514(){
     return ssNnsynsub[514];
    }
    double GetssNnsynsub515(){
     return ssNnsynsub[515];
    }
    double GetssNnsynsub516(){
     return ssNnsynsub[516];
    }
    double GetssNnsynsub517(){
     return ssNnsynsub[517];
    }
    double GetssNnsynsub518(){
     return ssNnsynsub[518];
    }
    double GetssNnsynsub519(){
     return ssNnsynsub[519];
    }
    double GetssNnsynsub520(){
     return ssNnsynsub[520];
    }
    double GetssNnsynsub521(){
     return ssNnsynsub[521];
    }
    double GetssNnsynsub522(){
     return ssNnsynsub[522];
    }
    double GetssNnsynsub523(){
     return ssNnsynsub[523];
    }
    double GetssNnsynsub524(){
     return ssNnsynsub[524];
    }
    double GetssNnsynsub525(){
     return ssNnsynsub[525];
    }
    double GetssNnsynsub526(){
     return ssNnsynsub[526];
    }
    double GetssNnsynsub527(){
     return ssNnsynsub[527];
    }
    double GetssNnsynsub528(){
     return ssNnsynsub[528];
    }
    double GetssNnsynsub529(){
     return ssNnsynsub[529];
    }
    double GetssNnsynsub530(){
     return ssNnsynsub[530];
    }
    double GetssNnsynsub531(){
     return ssNnsynsub[531];
    }
    double GetssNnsynsub532(){
     return ssNnsynsub[532];
    }
    double GetssNnsynsub533(){
     return ssNnsynsub[533];
    }
    double GetssNnsynsub534(){
     return ssNnsynsub[534];
    }
    double GetssNnsynsub535(){
     return ssNnsynsub[535];
    }
    double GetssNnsynsub536(){
     return ssNnsynsub[536];
    }
    double GetssNnsynsub537(){
     return ssNnsynsub[537];
    }
    double GetssNnsynsub538(){
     return ssNnsynsub[538];
    }
    double GetssNnsynsub539(){
     return ssNnsynsub[539];
    }
    double GetssNnsynsub540(){
     return ssNnsynsub[540];
    }
    double GetssNnsynsub541(){
     return ssNnsynsub[541];
    }
    double GetssNnsynsub542(){
     return ssNnsynsub[542];
    }
    double GetssNnsynsub543(){
     return ssNnsynsub[543];
    }
    double GetssNnsynsub544(){
     return ssNnsynsub[544];
    }
    double GetssNnsynsub545(){
     return ssNnsynsub[545];
    }
    double GetssNnsynsub546(){
     return ssNnsynsub[546];
    }
    double GetssNnsynsub547(){
     return ssNnsynsub[547];
    }
    double GetssNnsynsub548(){
     return ssNnsynsub[548];
    }
    double GetssNnsynsub549(){
     return ssNnsynsub[549];
    }
    double GetssNnsynsub550(){
     return ssNnsynsub[550];
    }
    double GetssNnsynsub551(){
     return ssNnsynsub[551];
    }
    double GetssNnsynsub552(){
     return ssNnsynsub[552];
    }
    double GetssNnsynsub553(){
     return ssNnsynsub[553];
    }
    double GetssNnsynsub554(){
     return ssNnsynsub[554];
    }
    double GetssNnsynsub555(){
     return ssNnsynsub[555];
    }
    double GetssNnsynsub556(){
     return ssNnsynsub[556];
    }
    double GetssNnsynsub557(){
     return ssNnsynsub[557];
    }
    double GetssNnsynsub558(){
     return ssNnsynsub[558];
    }
    double GetssNnsynsub559(){
     return ssNnsynsub[559];
    }
    double GetssNnsynsub560(){
     return ssNnsynsub[560];
    }
    double GetssNnsynsub561(){
     return ssNnsynsub[561];
    }
    double GetssNnsynsub562(){
     return ssNnsynsub[562];
    }
    double GetssNnsynsub563(){
     return ssNnsynsub[563];
    }
    double GetssNnsynsub564(){
     return ssNnsynsub[564];
    }
    double GetssNnsynsub565(){
     return ssNnsynsub[565];
    }
    double GetssNnsynsub566(){
     return ssNnsynsub[566];
    }
    double GetssNnsynsub567(){
     return ssNnsynsub[567];
    }
    double GetssNnsynsub568(){
     return ssNnsynsub[568];
    }
    double GetssNnsynsub569(){
     return ssNnsynsub[569];
    }
    double GetssNnsynsub570(){
     return ssNnsynsub[570];
    }
    double GetssNnsynsub571(){
     return ssNnsynsub[571];
    }
    double GetssNnsynsub572(){
     return ssNnsynsub[572];
    }
    double GetssNnsynsub573(){
     return ssNnsynsub[573];
    }
    double GetssNnsynsub574(){
     return ssNnsynsub[574];
    }
    double GetssNnsynsub575(){
     return ssNnsynsub[575];
    }
    double GetssNnsynsub576(){
     return ssNnsynsub[576];
    }
    double GetssNnsynsub577(){
     return ssNnsynsub[577];
    }
    double GetssNnsynsub578(){
     return ssNnsynsub[578];
    }
    double GetssNnsynsub579(){
     return ssNnsynsub[579];
    }
    double GetssNnsynsub580(){
     return ssNnsynsub[580];
    }
    double GetssNnsynsub581(){
     return ssNnsynsub[581];
    }
    double GetssNnsynsub582(){
     return ssNnsynsub[582];
    }
    double GetssNnsynsub583(){
     return ssNnsynsub[583];
    }
    double GetssNnsynsub584(){
     return ssNnsynsub[584];
    }
    double GetssNnsynsub585(){
     return ssNnsynsub[585];
    }
    double GetssNnsynsub586(){
     return ssNnsynsub[586];
    }
    double GetssNnsynsub587(){
     return ssNnsynsub[587];
    }
    double GetssNnsynsub588(){
     return ssNnsynsub[588];
    }
    double GetssNnsynsub589(){
     return ssNnsynsub[589];
    }
    double GetssNnsynsub590(){
     return ssNnsynsub[590];
    }
    double GetssNnsynsub591(){
     return ssNnsynsub[591];
    }
    double GetssNnsynsub592(){
     return ssNnsynsub[592];
    }
    double GetssNnsynsub593(){
     return ssNnsynsub[593];
    }
    double GetssNnsynsub594(){
     return ssNnsynsub[594];
    }
    double GetssNnsynsub595(){
     return ssNnsynsub[595];
    }
    double GetssNnsynsub596(){
     return ssNnsynsub[596];
    }
    double GetssNnsynsub597(){
     return ssNnsynsub[597];
    }
    double GetssNnsynsub598(){
     return ssNnsynsub[598];
    }
    double GetssNnsynsub599(){
     return ssNnsynsub[599];
    }
    double GetssNnsynsub600(){
     return ssNnsynsub[600];
    }
    double GetssNnsynsub601(){
     return ssNnsynsub[601];
    }
    double GetssNnsynsub602(){
     return ssNnsynsub[602];
    }
    double GetssNnsynsub603(){
     return ssNnsynsub[603];
    }
    double GetssNnsynsub604(){
     return ssNnsynsub[604];
    }
    double GetssNnsynsub605(){
     return ssNnsynsub[605];
    }
    double GetssNnsynsub606(){
     return ssNnsynsub[606];
    }
    double GetssNnsynsub607(){
     return ssNnsynsub[607];
    }
    double GetssNnsynsub608(){
     return ssNnsynsub[608];
    }
    double GetssNnsynsub609(){
     return ssNnsynsub[609];
    }
    double GetssNnsynsub610(){
     return ssNnsynsub[610];
    }
    double GetssNnsynsub611(){
     return ssNnsynsub[611];
    }
    double GetssNnsynsub612(){
     return ssNnsynsub[612];
    }
    double GetssNnsynsub613(){
     return ssNnsynsub[613];
    }
    double GetssNnsynsub614(){
     return ssNnsynsub[614];
    }
    double GetssNnsynsub615(){
     return ssNnsynsub[615];
    }
    double GetssNnsynsub616(){
     return ssNnsynsub[616];
    }
    double GetssNnsynsub617(){
     return ssNnsynsub[617];
    }
    double GetssNnsynsub618(){
     return ssNnsynsub[618];
    }
    double GetssNnsynsub619(){
     return ssNnsynsub[619];
    }
    double GetssNnsynsub620(){
     return ssNnsynsub[620];
    }
    double GetssNnsynsub621(){
     return ssNnsynsub[621];
    }
    double GetssNnsynsub622(){
     return ssNnsynsub[622];
    }
    double GetssNnsynsub623(){
     return ssNnsynsub[623];
    }
    double GetssNnsynsub624(){
     return ssNnsynsub[624];
    }
    double GetssNnsynsub625(){
     return ssNnsynsub[625];
    }
    double GetssNnsynsub626(){
     return ssNnsynsub[626];
    }
    double GetssNnsynsub627(){
     return ssNnsynsub[627];
    }
    double GetssNnsynsub628(){
     return ssNnsynsub[628];
    }
    double GetssNnsynsub629(){
     return ssNnsynsub[629];
    }
    double GetssNnsynsub630(){
     return ssNnsynsub[630];
    }
    double GetssNnsynsub631(){
     return ssNnsynsub[631];
    }
    double GetssNnsynsub632(){
     return ssNnsynsub[632];
    }
    double GetssNnsynsub633(){
     return ssNnsynsub[633];
    }
    double GetssNnsynsub634(){
     return ssNnsynsub[634];
    }
    double GetssNnsynsub635(){
     return ssNnsynsub[635];
    }
    double GetssNnsynsub636(){
     return ssNnsynsub[636];
    }
    double GetssNnsynsub637(){
     return ssNnsynsub[637];
    }
    double GetssNnsynsub638(){
     return ssNnsynsub[638];
    }
    double GetssNnsynsub639(){
     return ssNnsynsub[639];
    }
    double GetssNnsynsub640(){
     return ssNnsynsub[640];
    }
    double GetssNnsynsub641(){
     return ssNnsynsub[641];
    }
    double GetssNnsynsub642(){
     return ssNnsynsub[642];
    }
    double GetssNnsynsub643(){
     return ssNnsynsub[643];
    }
    double GetssNnsynsub644(){
     return ssNnsynsub[644];
    }
    double GetssNnsynsub645(){
     return ssNnsynsub[645];
    }
    double GetssNnsynsub646(){
     return ssNnsynsub[646];
    }
    double GetssNnsynsub647(){
     return ssNnsynsub[647];
    }
    double GetssNnsynsub648(){
     return ssNnsynsub[648];
    }
    double GetssNnsynsub649(){
     return ssNnsynsub[649];
    }
    double GetssNnsynsub650(){
     return ssNnsynsub[650];
    }
    double GetssNnsynsub651(){
     return ssNnsynsub[651];
    }
    double GetssNnsynsub652(){
     return ssNnsynsub[652];
    }
    double GetssNnsynsub653(){
     return ssNnsynsub[653];
    }
    double GetssNnsynsub654(){
     return ssNnsynsub[654];
    }
    double GetssNnsynsub655(){
     return ssNnsynsub[655];
    }
    double GetssNnsynsub656(){
     return ssNnsynsub[656];
    }
    double GetssNnsynsub657(){
     return ssNnsynsub[657];
    }
    double GetssNnsynsub658(){
     return ssNnsynsub[658];
    }
    double GetssNnsynsub659(){
     return ssNnsynsub[659];
    }
    double GetssNnsynsub660(){
     return ssNnsynsub[660];
    }
    double GetssNnsynsub661(){
     return ssNnsynsub[661];
    }
    double GetssNnsynsub662(){
     return ssNnsynsub[662];
    }
    double GetssNnsynsub663(){
     return ssNnsynsub[663];
    }
    double GetssNnsynsub664(){
     return ssNnsynsub[664];
    }
    double GetssNnsynsub665(){
     return ssNnsynsub[665];
    }
    double GetssNnsynsub666(){
     return ssNnsynsub[666];
    }
    double GetssNnsynsub667(){
     return ssNnsynsub[667];
    }
    double GetssNnsynsub668(){
     return ssNnsynsub[668];
    }
    double GetssNnsynsub669(){
     return ssNnsynsub[669];
    }
    double GetssNnsynsub670(){
     return ssNnsynsub[670];
    }
    double GetssNnsynsub671(){
     return ssNnsynsub[671];
    }
    double GetssNnsynsub672(){
     return ssNnsynsub[672];
    }
    double GetssNnsynsub673(){
     return ssNnsynsub[673];
    }
    double GetssNnsynsub674(){
     return ssNnsynsub[674];
    }
    double GetssNnsynsub675(){
     return ssNnsynsub[675];
    }
    double GetssNnsynsub676(){
     return ssNnsynsub[676];
    }
    double GetssNnsynsub677(){
     return ssNnsynsub[677];
    }
    double GetssNnsynsub678(){
     return ssNnsynsub[678];
    }
    double GetssNnsynsub679(){
     return ssNnsynsub[679];
    }
    double GetssNnsynsub680(){
     return ssNnsynsub[680];
    }
    double GetssNnsynsub681(){
     return ssNnsynsub[681];
    }
    double GetssNnsynsub682(){
     return ssNnsynsub[682];
    }
    double GetssNnsynsub683(){
     return ssNnsynsub[683];
    }
    double GetssNnsynsub684(){
     return ssNnsynsub[684];
    }
    double GetssNnsynsub685(){
     return ssNnsynsub[685];
    }
    double GetssNnsynsub686(){
     return ssNnsynsub[686];
    }
    double GetssNnsynsub687(){
     return ssNnsynsub[687];
    }
    double GetssNnsynsub688(){
     return ssNnsynsub[688];
    }
    double GetssNnsynsub689(){
     return ssNnsynsub[689];
    }
    double GetssNnsynsub690(){
     return ssNnsynsub[690];
    }
    double GetssNnsynsub691(){
     return ssNnsynsub[691];
    }
    double GetssNnsynsub692(){
     return ssNnsynsub[692];
    }
    double GetssNnsynsub693(){
     return ssNnsynsub[693];
    }
    double GetssNnsynsub694(){
     return ssNnsynsub[694];
    }
    double GetssNnsynsub695(){
     return ssNnsynsub[695];
    }
    double GetssNnsynsub696(){
     return ssNnsynsub[696];
    }
    double GetssNnsynsub697(){
     return ssNnsynsub[697];
    }
    double GetssNnsynsub698(){
     return ssNnsynsub[698];
    }
    double GetssNnsynsub699(){
     return ssNnsynsub[699];
    }
    double GetssNnsynsub700(){
     return ssNnsynsub[700];
    }
    double GetssNnsynsub701(){
     return ssNnsynsub[701];
    }
    double GetssNnsynsub702(){
     return ssNnsynsub[702];
    }
    double GetssNnsynsub703(){
     return ssNnsynsub[703];
    }
    double GetssNnsynsub704(){
     return ssNnsynsub[704];
    }
    double GetssNnsynsub705(){
     return ssNnsynsub[705];
    }
    double GetssNnsynsub706(){
     return ssNnsynsub[706];
    }
    double GetssNnsynsub707(){
     return ssNnsynsub[707];
    }
    double GetssNnsynsub708(){
     return ssNnsynsub[708];
    }
    double GetssNnsynsub709(){
     return ssNnsynsub[709];
    }
    double GetssNnsynsub710(){
     return ssNnsynsub[710];
    }
    double GetssNnsynsub711(){
     return ssNnsynsub[711];
    }
    double GetssNnsynsub712(){
     return ssNnsynsub[712];
    }
    double GetssNnsynsub713(){
     return ssNnsynsub[713];
    }
    double GetssNnsynsub714(){
     return ssNnsynsub[714];
    }
    double GetssNnsynsub715(){
     return ssNnsynsub[715];
    }
    double GetssNnsynsub716(){
     return ssNnsynsub[716];
    }
    double GetssNnsynsub717(){
     return ssNnsynsub[717];
    }
    double GetssNnsynsub718(){
     return ssNnsynsub[718];
    }
    double GetssNnsynsub719(){
     return ssNnsynsub[719];
    }
    double GetssNnsynsub720(){
     return ssNnsynsub[720];
    }
    double GetssNnsynsub721(){
     return ssNnsynsub[721];
    }
    double GetssNnsynsub722(){
     return ssNnsynsub[722];
    }
    double GetssNnsynsub723(){
     return ssNnsynsub[723];
    }
    double GetssNnsynsub724(){
     return ssNnsynsub[724];
    }
    double GetssNnsynsub725(){
     return ssNnsynsub[725];
    }
    double GetssNnsynsub726(){
     return ssNnsynsub[726];
    }
    double GetssNnsynsub727(){
     return ssNnsynsub[727];
    }
    double GetssNnsynsub728(){
     return ssNnsynsub[728];
    }
    double GetssNnsynsub729(){
     return ssNnsynsub[729];
    }
    double GetssNnsynsub730(){
     return ssNnsynsub[730];
    }
    double GetssNnsynsub731(){
     return ssNnsynsub[731];
    }
    double GetssNnsynsub732(){
     return ssNnsynsub[732];
    }
    double GetssNnsynsub733(){
     return ssNnsynsub[733];
    }
    double GetssNnsynsub734(){
     return ssNnsynsub[734];
    }
    double GetssNnsynsub735(){
     return ssNnsynsub[735];
    }
    double GetssNnsynsub736(){
     return ssNnsynsub[736];
    }
    double GetssNnsynsub737(){
     return ssNnsynsub[737];
    }
    double GetssNnsynsub738(){
     return ssNnsynsub[738];
    }
    double GetssNnsynsub739(){
     return ssNnsynsub[739];
    }
    double GetssNnsynsub740(){
     return ssNnsynsub[740];
    }
    double GetssNnsynsub741(){
     return ssNnsynsub[741];
    }
    double GetssNnsynsub742(){
     return ssNnsynsub[742];
    }
    double GetssNnsynsub743(){
     return ssNnsynsub[743];
    }
    double GetssNnsynsub744(){
     return ssNnsynsub[744];
    }
    double GetssNnsynsub745(){
     return ssNnsynsub[745];
    }
    double GetssNnsynsub746(){
     return ssNnsynsub[746];
    }
    double GetssNnsynsub747(){
     return ssNnsynsub[747];
    }
    double GetssNnsynsub748(){
     return ssNnsynsub[748];
    }
    double GetssNnsynsub749(){
     return ssNnsynsub[749];
    }
    double GetssNnsynsub750(){
     return ssNnsynsub[750];
    }
    double GetssNnsynsub751(){
     return ssNnsynsub[751];
    }
    double GetssNnsynsub752(){
     return ssNnsynsub[752];
    }
    double GetssNnsynsub753(){
     return ssNnsynsub[753];
    }
    double GetssNnsynsub754(){
     return ssNnsynsub[754];
    }
    double GetssNnsynsub755(){
     return ssNnsynsub[755];
    }
    double GetssNnsynsub756(){
     return ssNnsynsub[756];
    }
    double GetssNnsynsub757(){
     return ssNnsynsub[757];
    }
    double GetssNnsynsub758(){
     return ssNnsynsub[758];
    }
    double GetssNnsynsub759(){
     return ssNnsynsub[759];
    }
    double GetssNnsynsub760(){
     return ssNnsynsub[760];
    }
    double GetssNnsynsub761(){
     return ssNnsynsub[761];
    }
    double GetssNnsynsub762(){
     return ssNnsynsub[762];
    }
    double GetssNnsynsub763(){
     return ssNnsynsub[763];
    }
    double GetssNnsynsub764(){
     return ssNnsynsub[764];
    }
    double GetssNnsynsub765(){
     return ssNnsynsub[765];
    }
    double GetssNnsynsub766(){
     return ssNnsynsub[766];
    }
    double GetssNnsynsub767(){
     return ssNnsynsub[767];
    }
    double GetssNnsynsub768(){
     return ssNnsynsub[768];
    }
    double GetssNnsynsub769(){
     return ssNnsynsub[769];
    }
    double GetssNnsynsub770(){
     return ssNnsynsub[770];
    }
    double GetssNnsynsub771(){
     return ssNnsynsub[771];
    }
    double GetssNnsynsub772(){
     return ssNnsynsub[772];
    }
    double GetssNnsynsub773(){
     return ssNnsynsub[773];
    }
    double GetssNnsynsub774(){
     return ssNnsynsub[774];
    }
    double GetssNnsynsub775(){
     return ssNnsynsub[775];
    }
    double GetssNnsynsub776(){
     return ssNnsynsub[776];
    }
    double GetssNnsynsub777(){
     return ssNnsynsub[777];
    }
    double GetssNnsynsub778(){
     return ssNnsynsub[778];
    }
    double GetssNnsynsub779(){
     return ssNnsynsub[779];
    }
    double GetssNnsynsub780(){
     return ssNnsynsub[780];
    }
    double GetssNnsynsub781(){
     return ssNnsynsub[781];
    }
    double GetssNnsynsub782(){
     return ssNnsynsub[782];
    }
    double GetssNnsynsub783(){
     return ssNnsynsub[783];
    }
    double GetssNnsynsub784(){
     return ssNnsynsub[784];
    }
    double GetssNnsynsub785(){
     return ssNnsynsub[785];
    }
    double GetssNnsynsub786(){
     return ssNnsynsub[786];
    }
    double GetssNnsynsub787(){
     return ssNnsynsub[787];
    }
    double GetssNnsynsub788(){
     return ssNnsynsub[788];
    }
    double GetssNnsynsub789(){
     return ssNnsynsub[789];
    }
    double GetssNnsynsub790(){
     return ssNnsynsub[790];
    }
    double GetssNnsynsub791(){
     return ssNnsynsub[791];
    }
    double GetssNnsynsub792(){
     return ssNnsynsub[792];
    }
    double GetssNnsynsub793(){
     return ssNnsynsub[793];
    }
    double GetssNnsynsub794(){
     return ssNnsynsub[794];
    }
    double GetssNnsynsub795(){
     return ssNnsynsub[795];
    }
    double GetssNnsynsub796(){
     return ssNnsynsub[796];
    }
    double GetssNnsynsub797(){
     return ssNnsynsub[797];
    }
    double GetssNnsynsub798(){
     return ssNnsynsub[798];
    }
    double GetssNnsynsub799(){
     return ssNnsynsub[799];
    }
    double GetssNnsynsub800(){
     return ssNnsynsub[800];
    }
    double GetssNnsynsub801(){
     return ssNnsynsub[801];
    }
    double GetssNnsynsub802(){
     return ssNnsynsub[802];
    }
    double GetssNnsynsub803(){
     return ssNnsynsub[803];
    }
    double GetssNnsynsub804(){
     return ssNnsynsub[804];
    }
    double GetssNnsynsub805(){
     return ssNnsynsub[805];
    }
    double GetssNnsynsub806(){
     return ssNnsynsub[806];
    }
    double GetssNnsynsub807(){
     return ssNnsynsub[807];
    }
    double GetssNnsynsub808(){
     return ssNnsynsub[808];
    }
    double GetssNnsynsub809(){
     return ssNnsynsub[809];
    }
    double GetssNnsynsub810(){
     return ssNnsynsub[810];
    }
    double GetssNnsynsub811(){
     return ssNnsynsub[811];
    }
    double GetssNnsynsub812(){
     return ssNnsynsub[812];
    }
    double GetssNnsynsub813(){
     return ssNnsynsub[813];
    }
    double GetssNnsynsub814(){
     return ssNnsynsub[814];
    }
    double GetssNnsynsub815(){
     return ssNnsynsub[815];
    }
    double GetssNnsynsub816(){
     return ssNnsynsub[816];
    }
    double GetssNnsynsub817(){
     return ssNnsynsub[817];
    }
    double GetssNnsynsub818(){
     return ssNnsynsub[818];
    }
    double GetssNnsynsub819(){
     return ssNnsynsub[819];
    }
    double GetssNnsynsub820(){
     return ssNnsynsub[820];
    }
    double GetssNnsynsub821(){
     return ssNnsynsub[821];
    }
    double GetssNnsynsub822(){
     return ssNnsynsub[822];
    }
    double GetssNnsynsub823(){
     return ssNnsynsub[823];
    }
    double GetssNnsynsub824(){
     return ssNnsynsub[824];
    }
    double GetssNnsynsub825(){
     return ssNnsynsub[825];
    }
    double GetssNnsynsub826(){
     return ssNnsynsub[826];
    }
    double GetssNnsynsub827(){
     return ssNnsynsub[827];
    }
    double GetssNnsynsub828(){
     return ssNnsynsub[828];
    }
    double GetssNnsynsub829(){
     return ssNnsynsub[829];
    }
    double GetssNnsynsub830(){
     return ssNnsynsub[830];
    }
    double GetssNnsynsub831(){
     return ssNnsynsub[831];
    }
    double GetssNnsynsub832(){
     return ssNnsynsub[832];
    }
    double GetssNnsynsub833(){
     return ssNnsynsub[833];
    }
    double GetssNnsynsub834(){
     return ssNnsynsub[834];
    }
    double GetssNnsynsub835(){
     return ssNnsynsub[835];
    }
    double GetssNnsynsub836(){
     return ssNnsynsub[836];
    }
    double GetssNnsynsub837(){
     return ssNnsynsub[837];
    }
    double GetssNnsynsub838(){
     return ssNnsynsub[838];
    }
    double GetssNnsynsub839(){
     return ssNnsynsub[839];
    }
    double GetssNnsynsub840(){
     return ssNnsynsub[840];
    }
    double GetssNnsynsub841(){
     return ssNnsynsub[841];
    }
    double GetssNnsynsub842(){
     return ssNnsynsub[842];
    }
    double GetssNnsynsub843(){
     return ssNnsynsub[843];
    }
    double GetssNnsynsub844(){
     return ssNnsynsub[844];
    }
    double GetssNnsynsub845(){
     return ssNnsynsub[845];
    }
    double GetssNnsynsub846(){
     return ssNnsynsub[846];
    }
    double GetssNnsynsub847(){
     return ssNnsynsub[847];
    }
    double GetssNnsynsub848(){
     return ssNnsynsub[848];
    }
    double GetssNnsynsub849(){
     return ssNnsynsub[849];
    }
    double GetssNnsynsub850(){
     return ssNnsynsub[850];
    }
    double GetssNnsynsub851(){
     return ssNnsynsub[851];
    }
    double GetssNnsynsub852(){
     return ssNnsynsub[852];
    }
    double GetssNnsynsub853(){
     return ssNnsynsub[853];
    }
    double GetssNnsynsub854(){
     return ssNnsynsub[854];
    }
    double GetssNnsynsub855(){
     return ssNnsynsub[855];
    }
    double GetssNnsynsub856(){
     return ssNnsynsub[856];
    }
    double GetssNnsynsub857(){
     return ssNnsynsub[857];
    }
    double GetssNnsynsub858(){
     return ssNnsynsub[858];
    }
    double GetssNnsynsub859(){
     return ssNnsynsub[859];
    }
    double GetssNnsynsub860(){
     return ssNnsynsub[860];
    }
    double GetssNnsynsub861(){
     return ssNnsynsub[861];
    }
    double GetssNnsynsub862(){
     return ssNnsynsub[862];
    }
    double GetssNnsynsub863(){
     return ssNnsynsub[863];
    }
    double GetssNnsynsub864(){
     return ssNnsynsub[864];
    }
    double GetssNnsynsub865(){
     return ssNnsynsub[865];
    }
    double GetssNnsynsub866(){
     return ssNnsynsub[866];
    }
    double GetssNnsynsub867(){
     return ssNnsynsub[867];
    }
    double GetssNnsynsub868(){
     return ssNnsynsub[868];
    }
    double GetssNnsynsub869(){
     return ssNnsynsub[869];
    }
    double GetssNnsynsub870(){
     return ssNnsynsub[870];
    }
    double GetssNnsynsub871(){
     return ssNnsynsub[871];
    }
    double GetssNnsynsub872(){
     return ssNnsynsub[872];
    }
    double GetssNnsynsub873(){
     return ssNnsynsub[873];
    }
    double GetssNnsynsub874(){
     return ssNnsynsub[874];
    }
    double GetssNnsynsub875(){
     return ssNnsynsub[875];
    }
    double GetssNnsynsub876(){
     return ssNnsynsub[876];
    }
    double GetssNnsynsub877(){
     return ssNnsynsub[877];
    }
    double GetssNnsynsub878(){
     return ssNnsynsub[878];
    }
    double GetssNnsynsub879(){
     return ssNnsynsub[879];
    }
    double GetssNnsynsub880(){
     return ssNnsynsub[880];
    }
    double GetssNnsynsub881(){
     return ssNnsynsub[881];
    }
    double GetssNnsynsub882(){
     return ssNnsynsub[882];
    }
    double GetssNnsynsub883(){
     return ssNnsynsub[883];
    }
    double GetssNnsynsub884(){
     return ssNnsynsub[884];
    }
    double GetssNnsynsub885(){
     return ssNnsynsub[885];
    }
    double GetssNnsynsub886(){
     return ssNnsynsub[886];
    }
    double GetssNnsynsub887(){
     return ssNnsynsub[887];
    }
    double GetssNnsynsub888(){
     return ssNnsynsub[888];
    }
    double GetssNnsynsub889(){
     return ssNnsynsub[889];
    }
    double GetssNnsynsub890(){
     return ssNnsynsub[890];
    }
    double GetssNnsynsub891(){
     return ssNnsynsub[891];
    }
    double GetssNnsynsub892(){
     return ssNnsynsub[892];
    }
    double GetssNnsynsub893(){
     return ssNnsynsub[893];
    }
    double GetssNnsynsub894(){
     return ssNnsynsub[894];
    }
    double GetssNnsynsub895(){
     return ssNnsynsub[895];
    }
    double GetssNnsynsub896(){
     return ssNnsynsub[896];
    }
    double GetssNnsynsub897(){
     return ssNnsynsub[897];
    }
    double GetssNnsynsub898(){
     return ssNnsynsub[898];
    }
    double GetssNnsynsub899(){
     return ssNnsynsub[899];
    }
    double GetssNnsynsub900(){
     return ssNnsynsub[900];
    }
    double GetssNnsynsub901(){
     return ssNnsynsub[901];
    }
    double GetssNnsynsub902(){
     return ssNnsynsub[902];
    }
    double GetssNnsynsub903(){
     return ssNnsynsub[903];
    }
    double GetssNnsynsub904(){
     return ssNnsynsub[904];
    }
    double GetssNnsynsub905(){
     return ssNnsynsub[905];
    }
    double GetssNnsynsub906(){
     return ssNnsynsub[906];
    }
    double GetssNnsynsub907(){
     return ssNnsynsub[907];
    }
    double GetssNnsynsub908(){
     return ssNnsynsub[908];
    }
    double GetssNnsynsub909(){
     return ssNnsynsub[909];
    }
    double GetssNnsynsub910(){
     return ssNnsynsub[910];
    }
    double GetssNnsynsub911(){
     return ssNnsynsub[911];
    }
    double GetssNnsynsub912(){
     return ssNnsynsub[912];
    }
    double GetssNnsynsub913(){
     return ssNnsynsub[913];
    }
    double GetssNnsynsub914(){
     return ssNnsynsub[914];
    }
    double GetssNnsynsub915(){
     return ssNnsynsub[915];
    }
    double GetssNnsynsub916(){
     return ssNnsynsub[916];
    }
    double GetssNnsynsub917(){
     return ssNnsynsub[917];
    }
    double GetssNnsynsub918(){
     return ssNnsynsub[918];
    }
    double GetssNnsynsub919(){
     return ssNnsynsub[919];
    }
    double GetssNnsynsub920(){
     return ssNnsynsub[920];
    }
    double GetssNnsynsub921(){
     return ssNnsynsub[921];
    }
    double GetssNnsynsub922(){
     return ssNnsynsub[922];
    }
    double GetssNnsynsub923(){
     return ssNnsynsub[923];
    }
    double GetssNnsynsub924(){
     return ssNnsynsub[924];
    }
    double GetssNnsynsub925(){
     return ssNnsynsub[925];
    }
    double GetssNnsynsub926(){
     return ssNnsynsub[926];
    }
    double GetssNnsynsub927(){
     return ssNnsynsub[927];
    }
    double GetssNnsynsub928(){
     return ssNnsynsub[928];
    }
    double GetssNnsynsub929(){
     return ssNnsynsub[929];
    }
    double GetssNnsynsub930(){
     return ssNnsynsub[930];
    }
    double GetssNnsynsub931(){
     return ssNnsynsub[931];
    }
    double GetssNnsynsub932(){
     return ssNnsynsub[932];
    }
    double GetssNnsynsub933(){
     return ssNnsynsub[933];
    }
    double GetssNnsynsub934(){
     return ssNnsynsub[934];
    }
    double GetssNnsynsub935(){
     return ssNnsynsub[935];
    }
    double GetssNnsynsub936(){
     return ssNnsynsub[936];
    }
    double GetssNnsynsub937(){
     return ssNnsynsub[937];
    }
    double GetssNnsynsub938(){
     return ssNnsynsub[938];
    }
    double GetssNnsynsub939(){
     return ssNnsynsub[939];
    }
    double GetssNnsynsub940(){
     return ssNnsynsub[940];
    }
    double GetssNnsynsub941(){
     return ssNnsynsub[941];
    }
    double GetssNnsynsub942(){
     return ssNnsynsub[942];
    }
    double GetssNnsynsub943(){
     return ssNnsynsub[943];
    }
    double GetssNnsynsub944(){
     return ssNnsynsub[944];
    }
    double GetssNnsynsub945(){
     return ssNnsynsub[945];
    }
    double GetssNnsynsub946(){
     return ssNnsynsub[946];
    }
    double GetssNnsynsub947(){
     return ssNnsynsub[947];
    }
    double GetssNnsynsub948(){
     return ssNnsynsub[948];
    }
    double GetssNnsynsub949(){
     return ssNnsynsub[949];
    }
    double GetssNnsynsub950(){
     return ssNnsynsub[950];
    }
    double GetssNnsynsub951(){
     return ssNnsynsub[951];
    }
    double GetssNnsynsub952(){
     return ssNnsynsub[952];
    }
    double GetssNnsynsub953(){
     return ssNnsynsub[953];
    }
    double GetssNnsynsub954(){
     return ssNnsynsub[954];
    }
    double GetssNnsynsub955(){
     return ssNnsynsub[955];
    }
    double GetssNnsynsub956(){
     return ssNnsynsub[956];
    }
    double GetssNnsynsub957(){
     return ssNnsynsub[957];
    }
    double GetssNnsynsub958(){
     return ssNnsynsub[958];
    }
    double GetssNnsynsub959(){
     return ssNnsynsub[959];
    }
    double GetssNnsynsub960(){
     return ssNnsynsub[960];
    }
    double GetssNnsynsub961(){
     return ssNnsynsub[961];
    }
    double GetssNnsynsub962(){
     return ssNnsynsub[962];
    }
    double GetssNnsynsub963(){
     return ssNnsynsub[963];
    }
    double GetssNnsynsub964(){
     return ssNnsynsub[964];
    }
    double GetssNnsynsub965(){
     return ssNnsynsub[965];
    }
    double GetssNnsynsub966(){
     return ssNnsynsub[966];
    }
    double GetssNnsynsub967(){
     return ssNnsynsub[967];
    }
    double GetssNnsynsub968(){
     return ssNnsynsub[968];
    }
    double GetssNnsynsub969(){
     return ssNnsynsub[969];
    }
    double GetssNnsynsub970(){
     return ssNnsynsub[970];
    }
    double GetssNnsynsub971(){
     return ssNnsynsub[971];
    }
    double GetssNnsynsub972(){
     return ssNnsynsub[972];
    }
    double GetssNnsynsub973(){
     return ssNnsynsub[973];
    }
    double GetssNnsynsub974(){
     return ssNnsynsub[974];
    }
    double GetssNnsynsub975(){
     return ssNnsynsub[975];
    }
    double GetssNnsynsub976(){
     return ssNnsynsub[976];
    }
    double GetssNnsynsub977(){
     return ssNnsynsub[977];
    }
    double GetssNnsynsub978(){
     return ssNnsynsub[978];
    }
    double GetssNnsynsub979(){
     return ssNnsynsub[979];
    }
    double GetssNnsynsub980(){
     return ssNnsynsub[980];
    }
    double GetssNnsynsub981(){
     return ssNnsynsub[981];
    }
    double GetssNnsynsub982(){
     return ssNnsynsub[982];
    }
    double GetssNnsynsub983(){
     return ssNnsynsub[983];
    }
    double GetssNnsynsub984(){
     return ssNnsynsub[984];
    }
    double GetssNnsynsub985(){
     return ssNnsynsub[985];
    }
    double GetssNnsynsub986(){
     return ssNnsynsub[986];
    }
    double GetssNnsynsub987(){
     return ssNnsynsub[987];
    }
    double GetssNnsynsub988(){
     return ssNnsynsub[988];
    }
    double GetssNnsynsub989(){
     return ssNnsynsub[989];
    }
    double GetssNnsynsub990(){
     return ssNnsynsub[990];
    }
    double GetssNnsynsub991(){
     return ssNnsynsub[991];
    }
    double GetssNnsynsub992(){
     return ssNnsynsub[992];
    }
    double GetssNnsynsub993(){
     return ssNnsynsub[993];
    }
    double GetssNnsynsub994(){
     return ssNnsynsub[994];
    }
    double GetssNnsynsub995(){
     return ssNnsynsub[995];
    }
    double GetssNnsynsub996(){
     return ssNnsynsub[996];
    }
    double GetssNnsynsub997(){
     return ssNnsynsub[997];
    }
    double GetssNnsynsub998(){
     return ssNnsynsub[998];
    }
    double GetssNnsynsub999(){
     return ssNnsynsub[999];
    }
    double GetssNnsynsub1000(){
     return ssNnsynsub[1000];
    }
    double GetssNnsynsub1001(){
     return ssNnsynsub[1001];
    }
    double GetssNnsynsub1002(){
     return ssNnsynsub[1002];
    }
    double GetssNnsynsub1003(){
     return ssNnsynsub[1003];
    }
    double GetssNnsynsub1004(){
     return ssNnsynsub[1004];
    }
    double GetssNnsynsub1005(){
     return ssNnsynsub[1005];
    }
    double GetssNnsynsub1006(){
     return ssNnsynsub[1006];
    }
    double GetssNnsynsub1007(){
     return ssNnsynsub[1007];
    }
    double GetssNnsynsub1008(){
     return ssNnsynsub[1008];
    }
    double GetssNnsynsub1009(){
     return ssNnsynsub[1009];
    }
    double GetssNnsynsub1010(){
     return ssNnsynsub[1010];
    }
    double GetssNnsynsub1011(){
     return ssNnsynsub[1011];
    }
    double GetssNnsynsub1012(){
     return ssNnsynsub[1012];
    }
    double GetssNnsynsub1013(){
     return ssNnsynsub[1013];
    }
    double GetssNnsynsub1014(){
     return ssNnsynsub[1014];
    }
    double GetssNnsynsub1015(){
     return ssNnsynsub[1015];
    }
    double GetssNnsynsub1016(){
     return ssNnsynsub[1016];
    }
    double GetssNnsynsub1017(){
     return ssNnsynsub[1017];
    }
    double GetssNnsynsub1018(){
     return ssNnsynsub[1018];
    }
    double GetssNnsynsub1019(){
     return ssNnsynsub[1019];
    }
    double GetssNnsynsub1020(){
     return ssNnsynsub[1020];
    }
    double GetssNnsynsub1021(){
     return ssNnsynsub[1021];
    }
    double GetssNnsynsub1022(){
     return ssNnsynsub[1022];
    }
    double GetssNnsynsub1023(){
     return ssNnsynsub[1023];
    }
    double GetssNnsynsub1024(){
     return ssNnsynsub[1024];
    }
    double GetssNnsynsub1025(){
     return ssNnsynsub[1025];
    }
    double GetssNnsynsub1026(){
     return ssNnsynsub[1026];
    }
    double GetssNnsynsub1027(){
     return ssNnsynsub[1027];
    }
    double GetssNnsynsub1028(){
     return ssNnsynsub[1028];
    }
    double GetssNnsynsub1029(){
     return ssNnsynsub[1029];
    }
    double GetssNnsynsub1030(){
     return ssNnsynsub[1030];
    }
    double GetssNnsynsub1031(){
     return ssNnsynsub[1031];
    }
    double GetssNnsynsub1032(){
     return ssNnsynsub[1032];
    }
    double GetssNnsynsub1033(){
     return ssNnsynsub[1033];
    }
    double GetssNnsynsub1034(){
     return ssNnsynsub[1034];
    }
    double GetssNnsynsub1035(){
     return ssNnsynsub[1035];
    }
    double GetssNnsynsub1036(){
     return ssNnsynsub[1036];
    }
    double GetssNnsynsub1037(){
     return ssNnsynsub[1037];
    }
    double GetssNnsynsub1038(){
     return ssNnsynsub[1038];
    }
    double GetssNnsynsub1039(){
     return ssNnsynsub[1039];
    }
    double GetssNnsynsub1040(){
     return ssNnsynsub[1040];
    }
    double GetssNnsynsub1041(){
     return ssNnsynsub[1041];
    }
    double GetssNnsynsub1042(){
     return ssNnsynsub[1042];
    }
    double GetssNnsynsub1043(){
     return ssNnsynsub[1043];
    }
    double GetssNnsynsub1044(){
     return ssNnsynsub[1044];
    }
    double GetssNnsynsub1045(){
     return ssNnsynsub[1045];
    }
    double GetssNnsynsub1046(){
     return ssNnsynsub[1046];
    }
    double GetssNnsynsub1047(){
     return ssNnsynsub[1047];
    }
    double GetssNnsynsub1048(){
     return ssNnsynsub[1048];
    }
    double GetssNnsynsub1049(){
     return ssNnsynsub[1049];
    }
    double GetssNnsynsub1050(){
     return ssNnsynsub[1050];
    }
    double GetssNnsynsub1051(){
     return ssNnsynsub[1051];
    }
    double GetssNnsynsub1052(){
     return ssNnsynsub[1052];
    }
    double GetssNnsynsub1053(){
     return ssNnsynsub[1053];
    }
    double GetssNnsynsub1054(){
     return ssNnsynsub[1054];
    }
    double GetssNnsynsub1055(){
     return ssNnsynsub[1055];
    }
    double GetssNnsynsub1056(){
     return ssNnsynsub[1056];
    }
    double GetssNnsynsub1057(){
     return ssNnsynsub[1057];
    }
    double GetssNnsynsub1058(){
     return ssNnsynsub[1058];
    }
    double GetssNnsynsub1059(){
     return ssNnsynsub[1059];
    }
    double GetssNnsynsub1060(){
     return ssNnsynsub[1060];
    }
    double GetssNnsynsub1061(){
     return ssNnsynsub[1061];
    }
    double GetssNnsynsub1062(){
     return ssNnsynsub[1062];
    }
    double GetssNnsynsub1063(){
     return ssNnsynsub[1063];
    }
    double GetssNnsynsub1064(){
     return ssNnsynsub[1064];
    }
    double GetssNnsynsub1065(){
     return ssNnsynsub[1065];
    }
    double GetssNnsynsub1066(){
     return ssNnsynsub[1066];
    }
    double GetssNnsynsub1067(){
     return ssNnsynsub[1067];
    }
    double GetssNnsynsub1068(){
     return ssNnsynsub[1068];
    }
    double GetssNnsynsub1069(){
     return ssNnsynsub[1069];
    }
    double GetssNnsynsub1070(){
     return ssNnsynsub[1070];
    }
    double GetssNnsynsub1071(){
     return ssNnsynsub[1071];
    }
    double GetssNnsynsub1072(){
     return ssNnsynsub[1072];
    }
    double GetssNnsynsub1073(){
     return ssNnsynsub[1073];
    }
    double GetssNnsynsub1074(){
     return ssNnsynsub[1074];
    }
    double GetssNnsynsub1075(){
     return ssNnsynsub[1075];
    }
    double GetssNnsynsub1076(){
     return ssNnsynsub[1076];
    }
    double GetssNnsynsub1077(){
     return ssNnsynsub[1077];
    }
    double GetssNnsynsub1078(){
     return ssNnsynsub[1078];
    }
    double GetssNnsynsub1079(){
     return ssNnsynsub[1079];
    }
    double GetssNnsynsub1080(){
     return ssNnsynsub[1080];
    }
    double GetssNnsynsub1081(){
     return ssNnsynsub[1081];
    }
    double GetssNnsynsub1082(){
     return ssNnsynsub[1082];
    }
    double GetssNnsynsub1083(){
     return ssNnsynsub[1083];
    }
    double GetssNnsynsub1084(){
     return ssNnsynsub[1084];
    }
    double GetssNnsynsub1085(){
     return ssNnsynsub[1085];
    }
    double GetssNnsynsub1086(){
     return ssNnsynsub[1086];
    }
    double GetssNnsynsub1087(){
     return ssNnsynsub[1087];
    }
    double GetssNnsynsub1088(){
     return ssNnsynsub[1088];
    }
    double GetssNnsynsub1089(){
     return ssNnsynsub[1089];
    }
    double GetssNnsynsub1090(){
     return ssNnsynsub[1090];
    }
    double GetssNnsynsub1091(){
     return ssNnsynsub[1091];
    }
    double GetssNnsynsub1092(){
     return ssNnsynsub[1092];
    }
    double GetssNnsynsub1093(){
     return ssNnsynsub[1093];
    }
    double GetssNnsynsub1094(){
     return ssNnsynsub[1094];
    }
    double GetssNnsynsub1095(){
     return ssNnsynsub[1095];
    }
    double GetssNnsynsub1096(){
     return ssNnsynsub[1096];
    }
    double GetssNnsynsub1097(){
     return ssNnsynsub[1097];
    }
    double GetssNnsynsub1098(){
     return ssNnsynsub[1098];
    }
    double GetssNnsynsub1099(){
     return ssNnsynsub[1099];
    }
    double GetssNnsynsub1100(){
     return ssNnsynsub[1100];
    }
    double GetssNnsynsub1101(){
     return ssNnsynsub[1101];
    }
    double GetssNnsynsub1102(){
     return ssNnsynsub[1102];
    }
    double GetssNnsynsub1103(){
     return ssNnsynsub[1103];
    }
    double GetssNnsynsub1104(){
     return ssNnsynsub[1104];
    }
    double GetssNnsynsub1105(){
     return ssNnsynsub[1105];
    }
    double GetssNnsynsub1106(){
     return ssNnsynsub[1106];
    }
    double GetssNnsynsub1107(){
     return ssNnsynsub[1107];
    }
    double GetssNnsynsub1108(){
     return ssNnsynsub[1108];
    }
    double GetssNnsynsub1109(){
     return ssNnsynsub[1109];
    }
    double GetssNnsynsub1110(){
     return ssNnsynsub[1110];
    }
    double GetssNnsynsub1111(){
     return ssNnsynsub[1111];
    }
    double GetssNnsynsub1112(){
     return ssNnsynsub[1112];
    }
    double GetssNnsynsub1113(){
     return ssNnsynsub[1113];
    }
    double GetssNnsynsub1114(){
     return ssNnsynsub[1114];
    }
    double GetssNnsynsub1115(){
     return ssNnsynsub[1115];
    }
    double GetssNnsynsub1116(){
     return ssNnsynsub[1116];
    }
    double GetssNnsynsub1117(){
     return ssNnsynsub[1117];
    }
    double GetssNnsynsub1118(){
     return ssNnsynsub[1118];
    }
    double GetssNnsynsub1119(){
     return ssNnsynsub[1119];
    }
    double GetssNnsynsub1120(){
     return ssNnsynsub[1120];
    }
    double GetssNnsynsub1121(){
     return ssNnsynsub[1121];
    }
    double GetssNnsynsub1122(){
     return ssNnsynsub[1122];
    }
    double GetssNnsynsub1123(){
     return ssNnsynsub[1123];
    }
    double GetssNnsynsub1124(){
     return ssNnsynsub[1124];
    }
    double GetssNnsynsub1125(){
     return ssNnsynsub[1125];
    }
    double GetssNnsynsub1126(){
     return ssNnsynsub[1126];
    }
    double GetssNnsynsub1127(){
     return ssNnsynsub[1127];
    }
    double GetssNnsynsub1128(){
     return ssNnsynsub[1128];
    }
    double GetssNnsynsub1129(){
     return ssNnsynsub[1129];
    }
    double GetssNnsynsub1130(){
     return ssNnsynsub[1130];
    }
    double GetssNnsynsub1131(){
     return ssNnsynsub[1131];
    }
    double GetssNnsynsub1132(){
     return ssNnsynsub[1132];
    }
    double GetssNnsynsub1133(){
     return ssNnsynsub[1133];
    }
    double GetssNnsynsub1134(){
     return ssNnsynsub[1134];
    }
    double GetssNnsynsub1135(){
     return ssNnsynsub[1135];
    }
    double GetssNnsynsub1136(){
     return ssNnsynsub[1136];
    }
    double GetssNnsynsub1137(){
     return ssNnsynsub[1137];
    }
    double GetssNnsynsub1138(){
     return ssNnsynsub[1138];
    }
    double GetssNnsynsub1139(){
     return ssNnsynsub[1139];
    }
    double GetssNnsynsub1140(){
     return ssNnsynsub[1140];
    }
    double GetssNnsynsub1141(){
     return ssNnsynsub[1141];
    }
    double GetssNnsynsub1142(){
     return ssNnsynsub[1142];
    }
    double GetssNnsynsub1143(){
     return ssNnsynsub[1143];
    }
    double GetssNnsynsub1144(){
     return ssNnsynsub[1144];
    }
    double GetssNnsynsub1145(){
     return ssNnsynsub[1145];
    }
    double GetssNnsynsub1146(){
     return ssNnsynsub[1146];
    }
    double GetssNnsynsub1147(){
     return ssNnsynsub[1147];
    }
    double GetssNnsynsub1148(){
     return ssNnsynsub[1148];
    }
    double GetssNnsynsub1149(){
     return ssNnsynsub[1149];
    }
    double GetssNnsynsub1150(){
     return ssNnsynsub[1150];
    }
    double GetssNnsynsub1151(){
     return ssNnsynsub[1151];
    }
    double GetssNnsynsub1152(){
     return ssNnsynsub[1152];
    }
    double GetssNnsynsub1153(){
     return ssNnsynsub[1153];
    }
    double GetssNnsynsub1154(){
     return ssNnsynsub[1154];
    }
    double GetssNnsynsub1155(){
     return ssNnsynsub[1155];
    }
    double GetssNnsynsub1156(){
     return ssNnsynsub[1156];
    }
    double GetssNnsynsub1157(){
     return ssNnsynsub[1157];
    }
    double GetssNnsynsub1158(){
     return ssNnsynsub[1158];
    }
    double GetssNnsynsub1159(){
     return ssNnsynsub[1159];
    }
    double GetssNnsynsub1160(){
     return ssNnsynsub[1160];
    }
    double GetssNnsynsub1161(){
     return ssNnsynsub[1161];
    }
    double GetssNnsynsub1162(){
     return ssNnsynsub[1162];
    }
    double GetssNnsynsub1163(){
     return ssNnsynsub[1163];
    }
    double GetssNnsynsub1164(){
     return ssNnsynsub[1164];
    }
    double GetssNnsynsub1165(){
     return ssNnsynsub[1165];
    }
    double GetssNnsynsub1166(){
     return ssNnsynsub[1166];
    }
    double GetssNnsynsub1167(){
     return ssNnsynsub[1167];
    }
    double GetssNnsynsub1168(){
     return ssNnsynsub[1168];
    }
    double GetssNnsynsub1169(){
     return ssNnsynsub[1169];
    }
    double GetssNnsynsub1170(){
     return ssNnsynsub[1170];
    }
    double GetssNnsynsub1171(){
     return ssNnsynsub[1171];
    }
    double GetssNnsynsub1172(){
     return ssNnsynsub[1172];
    }
    double GetssNnsynsub1173(){
     return ssNnsynsub[1173];
    }
    double GetssNnsynsub1174(){
     return ssNnsynsub[1174];
    }
    double GetssNnsynsub1175(){
     return ssNnsynsub[1175];
    }
    double GetssNnsynsub1176(){
     return ssNnsynsub[1176];
    }
    double GetssNnsynsub1177(){
     return ssNnsynsub[1177];
    }
    double GetssNnsynsub1178(){
     return ssNnsynsub[1178];
    }
    double GetssNnsynsub1179(){
     return ssNnsynsub[1179];
    }
    double GetssNnsynsub1180(){
     return ssNnsynsub[1180];
    }
    double GetssNnsynsub1181(){
     return ssNnsynsub[1181];
    }
    double GetssNnsynsub1182(){
     return ssNnsynsub[1182];
    }
    double GetssNnsynsub1183(){
     return ssNnsynsub[1183];
    }
    double GetssNnsynsub1184(){
     return ssNnsynsub[1184];
    }
    double GetssNnsynsub1185(){
     return ssNnsynsub[1185];
    }
    double GetssNnsynsub1186(){
     return ssNnsynsub[1186];
    }
    double GetssNnsynsub1187(){
     return ssNnsynsub[1187];
    }
    double GetssNnsynsub1188(){
     return ssNnsynsub[1188];
    }
    double GetssNnsynsub1189(){
     return ssNnsynsub[1189];
    }
    double GetssNnsynsub1190(){
     return ssNnsynsub[1190];
    }
    double GetssNnsynsub1191(){
     return ssNnsynsub[1191];
    }
    double GetssNnsynsub1192(){
     return ssNnsynsub[1192];
    }
    double GetssNnsynsub1193(){
     return ssNnsynsub[1193];
    }
    double GetssNnsynsub1194(){
     return ssNnsynsub[1194];
    }
    double GetssNnsynsub1195(){
     return ssNnsynsub[1195];
    }
    double GetssNnsynsub1196(){
     return ssNnsynsub[1196];
    }
    double GetssNnsynsub1197(){
     return ssNnsynsub[1197];
    }
    double GetssNnsynsub1198(){
     return ssNnsynsub[1198];
    }
    double GetssNnsynsub1199(){
     return ssNnsynsub[1199];
    }
    double GetssNnsynsub1200(){
     return ssNnsynsub[1200];
    }
    double GetssNnsynsub1201(){
     return ssNnsynsub[1201];
    }
    double GetssNnsynsub1202(){
     return ssNnsynsub[1202];
    }
    double GetssNnsynsub1203(){
     return ssNnsynsub[1203];
    }
    double GetssNnsynsub1204(){
     return ssNnsynsub[1204];
    }
    double GetssNnsynsub1205(){
     return ssNnsynsub[1205];
    }
    double GetssNnsynsub1206(){
     return ssNnsynsub[1206];
    }
    double GetssNnsynsub1207(){
     return ssNnsynsub[1207];
    }
    double GetssNnsynsub1208(){
     return ssNnsynsub[1208];
    }
    double GetssNnsynsub1209(){
     return ssNnsynsub[1209];
    }
    double GetssNnsynsub1210(){
     return ssNnsynsub[1210];
    }
    double GetssNnsynsub1211(){
     return ssNnsynsub[1211];
    }
    double GetssNnsynsub1212(){
     return ssNnsynsub[1212];
    }
    double GetssNnsynsub1213(){
     return ssNnsynsub[1213];
    }
    double GetssNnsynsub1214(){
     return ssNnsynsub[1214];
    }
    double GetssNnsynsub1215(){
     return ssNnsynsub[1215];
    }
    double GetssNnsynsub1216(){
     return ssNnsynsub[1216];
    }
    double GetssNnsynsub1217(){
     return ssNnsynsub[1217];
    }
    double GetssNnsynsub1218(){
     return ssNnsynsub[1218];
    }
    double GetssNnsynsub1219(){
     return ssNnsynsub[1219];
    }
    double GetssNnsynsub1220(){
     return ssNnsynsub[1220];
    }
    double GetssNnsynsub1221(){
     return ssNnsynsub[1221];
    }
    double GetssNnsynsub1222(){
     return ssNnsynsub[1222];
    }
    double GetssNnsynsub1223(){
     return ssNnsynsub[1223];
    }
    double GetssNnsynsub1224(){
     return ssNnsynsub[1224];
    }
    double GetssNnsynsub1225(){
     return ssNnsynsub[1225];
    }
    double GetssNnsynsub1226(){
     return ssNnsynsub[1226];
    }
    double GetssNnsynsub1227(){
     return ssNnsynsub[1227];
    }
    double GetssNnsynsub1228(){
     return ssNnsynsub[1228];
    }
    double GetssNnsynsub1229(){
     return ssNnsynsub[1229];
    }
    double GetssNnsynsub1230(){
     return ssNnsynsub[1230];
    }
    double GetssNnsynsub1231(){
     return ssNnsynsub[1231];
    }
    double GetssNnsynsub1232(){
     return ssNnsynsub[1232];
    }
    double GetssNnsynsub1233(){
     return ssNnsynsub[1233];
    }
    double GetssNnsynsub1234(){
     return ssNnsynsub[1234];
    }
    double GetssNnsynsub1235(){
     return ssNnsynsub[1235];
    }
    double GetssNnsynsub1236(){
     return ssNnsynsub[1236];
    }
    double GetssNnsynsub1237(){
     return ssNnsynsub[1237];
    }
    double GetssNnsynsub1238(){
     return ssNnsynsub[1238];
    }
    double GetssNnsynsub1239(){
     return ssNnsynsub[1239];
    }
    double GetssNnsynsub1240(){
     return ssNnsynsub[1240];
    }
    double GetssNnsynsub1241(){
     return ssNnsynsub[1241];
    }
    double GetssNnsynsub1242(){
     return ssNnsynsub[1242];
    }
    double GetssNnsynsub1243(){
     return ssNnsynsub[1243];
    }
    double GetssNnsynsub1244(){
     return ssNnsynsub[1244];
    }
    double GetssNnsynsub1245(){
     return ssNnsynsub[1245];
    }
    double GetssNnsynsub1246(){
     return ssNnsynsub[1246];
    }
    double GetssNnsynsub1247(){
     return ssNnsynsub[1247];
    }
    double GetssNnsynsub1248(){
     return ssNnsynsub[1248];
    }
    double GetssNnsynsub1249(){
     return ssNnsynsub[1249];
    }
    double GetssNnsynsub1250(){
     return ssNnsynsub[1250];
    }
    double GetssNnsynsub1251(){
     return ssNnsynsub[1251];
    }
    double GetssNnsynsub1252(){
     return ssNnsynsub[1252];
    }
    double GetssNnsynsub1253(){
     return ssNnsynsub[1253];
    }
    double GetssNnsynsub1254(){
     return ssNnsynsub[1254];
    }
    double GetssNnsynsub1255(){
     return ssNnsynsub[1255];
    }
    double GetssNnsynsub1256(){
     return ssNnsynsub[1256];
    }
    double GetssNnsynsub1257(){
     return ssNnsynsub[1257];
    }
    double GetssNnsynsub1258(){
     return ssNnsynsub[1258];
    }
    double GetssNnsynsub1259(){
     return ssNnsynsub[1259];
    }
    double GetssNnsynsub1260(){
     return ssNnsynsub[1260];
    }
    double GetssNnsynsub1261(){
     return ssNnsynsub[1261];
    }
    double GetssNnsynsub1262(){
     return ssNnsynsub[1262];
    }
    double GetssNnsynsub1263(){
     return ssNnsynsub[1263];
    }
    double GetssNnsynsub1264(){
     return ssNnsynsub[1264];
    }
    double GetssNnsynsub1265(){
     return ssNnsynsub[1265];
    }
    double GetssNnsynsub1266(){
     return ssNnsynsub[1266];
    }
    double GetssNnsynsub1267(){
     return ssNnsynsub[1267];
    }
    double GetssNnsynsub1268(){
     return ssNnsynsub[1268];
    }
    double GetssNnsynsub1269(){
     return ssNnsynsub[1269];
    }
    double GetssNnsynsub1270(){
     return ssNnsynsub[1270];
    }
    double GetssNnsynsub1271(){
     return ssNnsynsub[1271];
    }
    double GetssNnsynsub1272(){
     return ssNnsynsub[1272];
    }
    double GetssNnsynsub1273(){
     return ssNnsynsub[1273];
    }
    double GetssNnsynsub1274(){
     return ssNnsynsub[1274];
    }
    double GetssNnsynsub1275(){
     return ssNnsynsub[1275];
    }
    double GetssNnsynsub1276(){
     return ssNnsynsub[1276];
    }
    double GetssNnsynsub1277(){
     return ssNnsynsub[1277];
    }
    double GetssNnsynsub1278(){
     return ssNnsynsub[1278];
    }
    double GetssNnsynsub1279(){
     return ssNnsynsub[1279];
    }
    double GetssNnsynsub1280(){
     return ssNnsynsub[1280];
    }
    double GetssNnsynsub1281(){
     return ssNnsynsub[1281];
    }
    double GetssNnsynsub1282(){
     return ssNnsynsub[1282];
    }
    double GetssNnsynsub1283(){
     return ssNnsynsub[1283];
    }
    double GetssNnsynsub1284(){
     return ssNnsynsub[1284];
    }
    double GetssNnsynsub1285(){
     return ssNnsynsub[1285];
    }
    double GetssNnsynsub1286(){
     return ssNnsynsub[1286];
    }
    double GetssNnsynsub1287(){
     return ssNnsynsub[1287];
    }
    double GetssNnsynsub1288(){
     return ssNnsynsub[1288];
    }
    double GetssNnsynsub1289(){
     return ssNnsynsub[1289];
    }
    double GetssNnsynsub1290(){
     return ssNnsynsub[1290];
    }
    double GetssNnsynsub1291(){
     return ssNnsynsub[1291];
    }
    double GetssNnsynsub1292(){
     return ssNnsynsub[1292];
    }
    double GetssNnsynsub1293(){
     return ssNnsynsub[1293];
    }
    double GetssNnsynsub1294(){
     return ssNnsynsub[1294];
    }
    double GetssNnsynsub1295(){
     return ssNnsynsub[1295];
    }
    double GetssNnsynsub1296(){
     return ssNnsynsub[1296];
    }
    double GetssNnsynsub1297(){
     return ssNnsynsub[1297];
    }
    double GetssNnsynsub1298(){
     return ssNnsynsub[1298];
    }
    double GetssNnsynsub1299(){
     return ssNnsynsub[1299];
    }
    double GetssNnsynsub1300(){
     return ssNnsynsub[1300];
    }
    double GetssNnsynsub1301(){
     return ssNnsynsub[1301];
    }
    double GetssNnsynsub1302(){
     return ssNnsynsub[1302];
    }
    double GetssNnsynsub1303(){
     return ssNnsynsub[1303];
    }
    double GetssNnsynsub1304(){
     return ssNnsynsub[1304];
    }
    double GetssNnsynsub1305(){
     return ssNnsynsub[1305];
    }
    double GetssNnsynsub1306(){
     return ssNnsynsub[1306];
    }
    double GetssNnsynsub1307(){
     return ssNnsynsub[1307];
    }
    double GetssNnsynsub1308(){
     return ssNnsynsub[1308];
    }
    double GetssNnsynsub1309(){
     return ssNnsynsub[1309];
    }
    double GetssNnsynsub1310(){
     return ssNnsynsub[1310];
    }
    double GetssNnsynsub1311(){
     return ssNnsynsub[1311];
    }
    double GetssNnsynsub1312(){
     return ssNnsynsub[1312];
    }
    double GetssNnsynsub1313(){
     return ssNnsynsub[1313];
    }
    double GetssNnsynsub1314(){
     return ssNnsynsub[1314];
    }
    double GetssNnsynsub1315(){
     return ssNnsynsub[1315];
    }
    double GetssNnsynsub1316(){
     return ssNnsynsub[1316];
    }
    double GetssNnsynsub1317(){
     return ssNnsynsub[1317];
    }
    double GetssNnsynsub1318(){
     return ssNnsynsub[1318];
    }
    double GetssNnsynsub1319(){
     return ssNnsynsub[1319];
    }
    double GetssNnsynsub1320(){
     return ssNnsynsub[1320];
    }
    double GetssNnsynsub1321(){
     return ssNnsynsub[1321];
    }
    double GetssNnsynsub1322(){
     return ssNnsynsub[1322];
    }
    double GetssNnsynsub1323(){
     return ssNnsynsub[1323];
    }
    double GetssNnsynsub1324(){
     return ssNnsynsub[1324];
    }
    double GetssNnsynsub1325(){
     return ssNnsynsub[1325];
    }
    double GetssNnsynsub1326(){
     return ssNnsynsub[1326];
    }
    double GetssNnsynsub1327(){
     return ssNnsynsub[1327];
    }
    double GetssNnsynsub1328(){
     return ssNnsynsub[1328];
    }
    double GetssNnsynsub1329(){
     return ssNnsynsub[1329];
    }
    double GetssNnsynsub1330(){
     return ssNnsynsub[1330];
    }
    double GetssNnsynsub1331(){
     return ssNnsynsub[1331];
    }
    double GetssNnsynsub1332(){
     return ssNnsynsub[1332];
    }
    double GetssNnsynsub1333(){
     return ssNnsynsub[1333];
    }
    double GetssNnsynsub1334(){
     return ssNnsynsub[1334];
    }
    double GetssNnsynsub1335(){
     return ssNnsynsub[1335];
    }
    double GetssNnsynsub1336(){
     return ssNnsynsub[1336];
    }
    double GetssNnsynsub1337(){
     return ssNnsynsub[1337];
    }
    double GetssNnsynsub1338(){
     return ssNnsynsub[1338];
    }
    double GetssNnsynsub1339(){
     return ssNnsynsub[1339];
    }
    double GetssNnsynsub1340(){
     return ssNnsynsub[1340];
    }
    double GetssNnsynsub1341(){
     return ssNnsynsub[1341];
    }
    double GetssNnsynsub1342(){
     return ssNnsynsub[1342];
    }
    double GetssNnsynsub1343(){
     return ssNnsynsub[1343];
    }
    double GetssNnsynsub1344(){
     return ssNnsynsub[1344];
    }
    double GetssNnsynsub1345(){
     return ssNnsynsub[1345];
    }
    double GetssNnsynsub1346(){
     return ssNnsynsub[1346];
    }
    double GetssNnsynsub1347(){
     return ssNnsynsub[1347];
    }
    double GetssNnsynsub1348(){
     return ssNnsynsub[1348];
    }
    double GetssNnsynsub1349(){
     return ssNnsynsub[1349];
    }
    double GetssNnsynsub1350(){
     return ssNnsynsub[1350];
    }
    double GetssNnsynsub1351(){
     return ssNnsynsub[1351];
    }
    double GetssNnsynsub1352(){
     return ssNnsynsub[1352];
    }
    double GetssNnsynsub1353(){
     return ssNnsynsub[1353];
    }
    double GetssNnsynsub1354(){
     return ssNnsynsub[1354];
    }
    double GetssNnsynsub1355(){
     return ssNnsynsub[1355];
    }
    double GetssNnsynsub1356(){
     return ssNnsynsub[1356];
    }
    double GetssNnsynsub1357(){
     return ssNnsynsub[1357];
    }
    double GetssNnsynsub1358(){
     return ssNnsynsub[1358];
    }
    double GetssNnsynsub1359(){
     return ssNnsynsub[1359];
    }
    double GetssNnsynsub1360(){
     return ssNnsynsub[1360];
    }
    double GetssNnsynsub1361(){
     return ssNnsynsub[1361];
    }
    double GetssNnsynsub1362(){
     return ssNnsynsub[1362];
    }
    double GetssNnsynsub1363(){
     return ssNnsynsub[1363];
    }
    double GetssNnsynsub1364(){
     return ssNnsynsub[1364];
    }
    double GetssNnsynsub1365(){
     return ssNnsynsub[1365];
    }
    double GetssNnsynsub1366(){
     return ssNnsynsub[1366];
    }
    double GetssNnsynsub1367(){
     return ssNnsynsub[1367];
    }
    double GetssNnsynsub1368(){
     return ssNnsynsub[1368];
    }
    double GetssNnsynsub1369(){
     return ssNnsynsub[1369];
    }
    double GetssNnsynsub1370(){
     return ssNnsynsub[1370];
    }
    double GetssNnsynsub1371(){
     return ssNnsynsub[1371];
    }
    double GetssNnsynsub1372(){
     return ssNnsynsub[1372];
    }
    double GetssNnsynsub1373(){
     return ssNnsynsub[1373];
    }
    double GetssNnsynsub1374(){
     return ssNnsynsub[1374];
    }
    double GetssNnsynsub1375(){
     return ssNnsynsub[1375];
    }
    double GetssNnsynsub1376(){
     return ssNnsynsub[1376];
    }
    double GetssNnsynsub1377(){
     return ssNnsynsub[1377];
    }
    double GetssNnsynsub1378(){
     return ssNnsynsub[1378];
    }
    double GetssNnsynsub1379(){
     return ssNnsynsub[1379];
    }
    double GetssNnsynsub1380(){
     return ssNnsynsub[1380];
    }
    double GetssNnsynsub1381(){
     return ssNnsynsub[1381];
    }
    double GetssNnsynsub1382(){
     return ssNnsynsub[1382];
    }
    double GetssNnsynsub1383(){
     return ssNnsynsub[1383];
    }
    double GetssNnsynsub1384(){
     return ssNnsynsub[1384];
    }
    double GetssNnsynsub1385(){
     return ssNnsynsub[1385];
    }
    double GetssNnsynsub1386(){
     return ssNnsynsub[1386];
    }
    double GetssNnsynsub1387(){
     return ssNnsynsub[1387];
    }
    double GetssNnsynsub1388(){
     return ssNnsynsub[1388];
    }
    double GetssNnsynsub1389(){
     return ssNnsynsub[1389];
    }
    double GetssNnsynsub1390(){
     return ssNnsynsub[1390];
    }
    double GetssNnsynsub1391(){
     return ssNnsynsub[1391];
    }
    double GetssNnsynsub1392(){
     return ssNnsynsub[1392];
    }
    double GetssNnsynsub1393(){
     return ssNnsynsub[1393];
    }
    double GetssNnsynsub1394(){
     return ssNnsynsub[1394];
    }
    double GetssNnsynsub1395(){
     return ssNnsynsub[1395];
    }
    double GetssNnsynsub1396(){
     return ssNnsynsub[1396];
    }
    double GetssNnsynsub1397(){
     return ssNnsynsub[1397];
    }
    double GetssNnsynsub1398(){
     return ssNnsynsub[1398];
    }
    double GetssNnsynsub1399(){
     return ssNnsynsub[1399];
    }
    double GetssNnsynsub1400(){
     return ssNnsynsub[1400];
    }
    double GetssNnsynsub1401(){
     return ssNnsynsub[1401];
    }
    double GetssNnsynsub1402(){
     return ssNnsynsub[1402];
    }
    double GetssNnsynsub1403(){
     return ssNnsynsub[1403];
    }
    double GetssNnsynsub1404(){
     return ssNnsynsub[1404];
    }
    double GetssNnsynsub1405(){
     return ssNnsynsub[1405];
    }
    double GetssNnsynsub1406(){
     return ssNnsynsub[1406];
    }
    double GetssNnsynsub1407(){
     return ssNnsynsub[1407];
    }
    double GetssNnsynsub1408(){
     return ssNnsynsub[1408];
    }
    double GetssNnsynsub1409(){
     return ssNnsynsub[1409];
    }
    double GetssNnsynsub1410(){
     return ssNnsynsub[1410];
    }
    double GetssNnsynsub1411(){
     return ssNnsynsub[1411];
    }
    double GetssNnsynsub1412(){
     return ssNnsynsub[1412];
    }
    double GetssNnsynsub1413(){
     return ssNnsynsub[1413];
    }
    double GetssNnsynsub1414(){
     return ssNnsynsub[1414];
    }
    double GetssNnsynsub1415(){
     return ssNnsynsub[1415];
    }
    double GetssNnsynsub1416(){
     return ssNnsynsub[1416];
    }
    double GetssNnsynsub1417(){
     return ssNnsynsub[1417];
    }
    double GetssNnsynsub1418(){
     return ssNnsynsub[1418];
    }
    double GetssNnsynsub1419(){
     return ssNnsynsub[1419];
    }
    double GetssNnsynsub1420(){
     return ssNnsynsub[1420];
    }
    double GetssNnsynsub1421(){
     return ssNnsynsub[1421];
    }
    double GetssNnsynsub1422(){
     return ssNnsynsub[1422];
    }
    double GetssNnsynsub1423(){
     return ssNnsynsub[1423];
    }
    double GetssNnsynsub1424(){
     return ssNnsynsub[1424];
    }
    double GetssNnsynsub1425(){
     return ssNnsynsub[1425];
    }
    double GetssNnsynsub1426(){
     return ssNnsynsub[1426];
    }
    double GetssNnsynsub1427(){
     return ssNnsynsub[1427];
    }
    double GetssNnsynsub1428(){
     return ssNnsynsub[1428];
    }
    double GetssNnsynsub1429(){
     return ssNnsynsub[1429];
    }
    double GetssNnsynsub1430(){
     return ssNnsynsub[1430];
    }
    double GetssNnsynsub1431(){
     return ssNnsynsub[1431];
    }
    double GetssNnsynsub1432(){
     return ssNnsynsub[1432];
    }
    double GetssNnsynsub1433(){
     return ssNnsynsub[1433];
    }
    double GetssNnsynsub1434(){
     return ssNnsynsub[1434];
    }
    double GetssNnsynsub1435(){
     return ssNnsynsub[1435];
    }
    double GetssNnsynsub1436(){
     return ssNnsynsub[1436];
    }
    double GetssNnsynsub1437(){
     return ssNnsynsub[1437];
    }
    double GetssNnsynsub1438(){
     return ssNnsynsub[1438];
    }
    double GetssNnsynsub1439(){
     return ssNnsynsub[1439];
    }
    double GetssNnsynsub1440(){
     return ssNnsynsub[1440];
    }
    double GetssNnsynsub1441(){
     return ssNnsynsub[1441];
    }
    double GetssNnsynsub1442(){
     return ssNnsynsub[1442];
    }
    double GetssNnsynsub1443(){
     return ssNnsynsub[1443];
    }
    double GetssNnsynsub1444(){
     return ssNnsynsub[1444];
    }
    double GetssNnsynsub1445(){
     return ssNnsynsub[1445];
    }
    double GetssNnsynsub1446(){
     return ssNnsynsub[1446];
    }
    double GetssNnsynsub1447(){
     return ssNnsynsub[1447];
    }
    double GetssNnsynsub1448(){
     return ssNnsynsub[1448];
    }
    double GetssNnsynsub1449(){
     return ssNnsynsub[1449];
    }
    double GetssNnsynsub1450(){
     return ssNnsynsub[1450];
    }
    double GetssNnsynsub1451(){
     return ssNnsynsub[1451];
    }
    double GetssNnsynsub1452(){
     return ssNnsynsub[1452];
    }
    double GetssNnsynsub1453(){
     return ssNnsynsub[1453];
    }
    double GetssNnsynsub1454(){
     return ssNnsynsub[1454];
    }
    double GetssNnsynsub1455(){
     return ssNnsynsub[1455];
    }
    double GetssNnsynsub1456(){
     return ssNnsynsub[1456];
    }
    double GetssNnsynsub1457(){
     return ssNnsynsub[1457];
    }
    double GetssNnsynsub1458(){
     return ssNnsynsub[1458];
    }
    double GetssNnsynsub1459(){
     return ssNnsynsub[1459];
    }
    double GetssNnsynsub1460(){
     return ssNnsynsub[1460];
    }
    double GetssNnsynsub1461(){
     return ssNnsynsub[1461];
    }
    double GetssNnsynsub1462(){
     return ssNnsynsub[1462];
    }
    double GetssNnsynsub1463(){
     return ssNnsynsub[1463];
    }
    double GetssNnsynsub1464(){
     return ssNnsynsub[1464];
    }
    double GetssNnsynsub1465(){
     return ssNnsynsub[1465];
    }
    double GetssNnsynsub1466(){
     return ssNnsynsub[1466];
    }
    double GetssNnsynsub1467(){
     return ssNnsynsub[1467];
    }
    double GetssNnsynsub1468(){
     return ssNnsynsub[1468];
    }
    double GetssNnsynsub1469(){
     return ssNnsynsub[1469];
    }
    double GetssNnsynsub1470(){
     return ssNnsynsub[1470];
    }
    double GetssNnsynsub1471(){
     return ssNnsynsub[1471];
    }
    double GetssNnsynsub1472(){
     return ssNnsynsub[1472];
    }
    double GetssNnsynsub1473(){
     return ssNnsynsub[1473];
    }
    double GetssNnsynsub1474(){
     return ssNnsynsub[1474];
    }
    double GetssNnsynsub1475(){
     return ssNnsynsub[1475];
    }
    double GetssNnsynsub1476(){
     return ssNnsynsub[1476];
    }
    double GetssNnsynsub1477(){
     return ssNnsynsub[1477];
    }
    double GetssNnsynsub1478(){
     return ssNnsynsub[1478];
    }
    double GetssNnsynsub1479(){
     return ssNnsynsub[1479];
    }
    double GetssNnsynsub1480(){
     return ssNnsynsub[1480];
    }
    double GetssNnsynsub1481(){
     return ssNnsynsub[1481];
    }
    double GetssNnsynsub1482(){
     return ssNnsynsub[1482];
    }
    double GetssNnsynsub1483(){
     return ssNnsynsub[1483];
    }
    double GetssNnsynsub1484(){
     return ssNnsynsub[1484];
    }
    double GetssNnsynsub1485(){
     return ssNnsynsub[1485];
    }
    double GetssNnsynsub1486(){
     return ssNnsynsub[1486];
    }
    double GetssNnsynsub1487(){
     return ssNnsynsub[1487];
    }
    double GetssNnsynsub1488(){
     return ssNnsynsub[1488];
    }
    double GetssNnsynsub1489(){
     return ssNnsynsub[1489];
    }
    double GetssNnsynsub1490(){
     return ssNnsynsub[1490];
    }
    double GetssNnsynsub1491(){
     return ssNnsynsub[1491];
    }
    double GetssNnsynsub1492(){
     return ssNnsynsub[1492];
    }
    double GetssNnsynsub1493(){
     return ssNnsynsub[1493];
    }
    double GetssNnsynsub1494(){
     return ssNnsynsub[1494];
    }
    double GetssNnsynsub1495(){
     return ssNnsynsub[1495];
    }
    double GetssNnsynsub1496(){
     return ssNnsynsub[1496];
    }
    double GetssNnsynsub1497(){
     return ssNnsynsub[1497];
    }
    double GetssNnsynsub1498(){
     return ssNnsynsub[1498];
    }
    double GetssNnsynsub1499(){
     return ssNnsynsub[1499];
    }
    double GetssNnsynsub1500(){
     return ssNnsynsub[1500];
    }
    double GetssNnsynsub1501(){
     return ssNnsynsub[1501];
    }
    double GetssNnsynsub1502(){
     return ssNnsynsub[1502];
    }
    double GetssNnsynsub1503(){
     return ssNnsynsub[1503];
    }
    double GetssNnsynsub1504(){
     return ssNnsynsub[1504];
    }
    double GetssNnsynsub1505(){
     return ssNnsynsub[1505];
    }
    double GetssNnsynsub1506(){
     return ssNnsynsub[1506];
    }
    double GetssNnsynsub1507(){
     return ssNnsynsub[1507];
    }
    double GetssNnsynsub1508(){
     return ssNnsynsub[1508];
    }
    double GetssNnsynsub1509(){
     return ssNnsynsub[1509];
    }
    double GetssNnsynsub1510(){
     return ssNnsynsub[1510];
    }
    double GetssNnsynsub1511(){
     return ssNnsynsub[1511];
    }
    double GetssNnsynsub1512(){
     return ssNnsynsub[1512];
    }
    double GetssNnsynsub1513(){
     return ssNnsynsub[1513];
    }
    double GetssNnsynsub1514(){
     return ssNnsynsub[1514];
    }
    double GetssNnsynsub1515(){
     return ssNnsynsub[1515];
    }
    double GetssNnsynsub1516(){
     return ssNnsynsub[1516];
    }
    double GetssNnsynsub1517(){
     return ssNnsynsub[1517];
    }
    double GetssNnsynsub1518(){
     return ssNnsynsub[1518];
    }
    double GetssNnsynsub1519(){
     return ssNnsynsub[1519];
    }
    double GetssNnsynsub1520(){
     return ssNnsynsub[1520];
    }
    double GetssNnsynsub1521(){
     return ssNnsynsub[1521];
    }
    double GetssNnsynsub1522(){
     return ssNnsynsub[1522];
    }
    double GetssNnsynsub1523(){
     return ssNnsynsub[1523];
    }
    double GetssNnsynsub1524(){
     return ssNnsynsub[1524];
    }
    double GetssNnsynsub1525(){
     return ssNnsynsub[1525];
    }
    double GetssNnsynsub1526(){
     return ssNnsynsub[1526];
    }
    double GetssNnsynsub1527(){
     return ssNnsynsub[1527];
    }
    double GetssNnsynsub1528(){
     return ssNnsynsub[1528];
    }
    double GetssNnsynsub1529(){
     return ssNnsynsub[1529];
    }
    double GetssNnsynsub1530(){
     return ssNnsynsub[1530];
    }
    double GetssNnsynsub1531(){
     return ssNnsynsub[1531];
    }
    double GetssNnsynsub1532(){
     return ssNnsynsub[1532];
    }
    double GetssNnsynsub1533(){
     return ssNnsynsub[1533];
    }
    double GetssNnsynsub1534(){
     return ssNnsynsub[1534];
    }
    double GetssNnsynsub1535(){
     return ssNnsynsub[1535];
    }
    double GetssNnsynsub1536(){
     return ssNnsynsub[1536];
    }
    double GetssNnsynsub1537(){
     return ssNnsynsub[1537];
    }
    double GetssNnsynsub1538(){
     return ssNnsynsub[1538];
    }
    double GetssNnsynsub1539(){
     return ssNnsynsub[1539];
    }
    double GetssNnsynsub1540(){
     return ssNnsynsub[1540];
    }
    double GetssNnsynsub1541(){
     return ssNnsynsub[1541];
    }
    double GetssNnsynsub1542(){
     return ssNnsynsub[1542];
    }
    double GetssNnsynsub1543(){
     return ssNnsynsub[1543];
    }
    double GetssNnsynsub1544(){
     return ssNnsynsub[1544];
    }
    double GetssNnsynsub1545(){
     return ssNnsynsub[1545];
    }
    double GetssNnsynsub1546(){
     return ssNnsynsub[1546];
    }
    double GetssNnsynsub1547(){
     return ssNnsynsub[1547];
    }
    double GetssNnsynsub1548(){
     return ssNnsynsub[1548];
    }
    double GetssNnsynsub1549(){
     return ssNnsynsub[1549];
    }
    double GetssNnsynsub1550(){
     return ssNnsynsub[1550];
    }
    double GetssNnsynsub1551(){
     return ssNnsynsub[1551];
    }
    double GetssNnsynsub1552(){
     return ssNnsynsub[1552];
    }
    double GetssNnsynsub1553(){
     return ssNnsynsub[1553];
    }
    double GetssNnsynsub1554(){
     return ssNnsynsub[1554];
    }
    double GetssNnsynsub1555(){
     return ssNnsynsub[1555];
    }
    double GetssNnsynsub1556(){
     return ssNnsynsub[1556];
    }
    double GetssNnsynsub1557(){
     return ssNnsynsub[1557];
    }
    double GetssNnsynsub1558(){
     return ssNnsynsub[1558];
    }
    double GetssNnsynsub1559(){
     return ssNnsynsub[1559];
    }
    double GetssNnsynsub1560(){
     return ssNnsynsub[1560];
    }
    double GetssNnsynsub1561(){
     return ssNnsynsub[1561];
    }
    double GetssNnsynsub1562(){
     return ssNnsynsub[1562];
    }
    double GetssNnsynsub1563(){
     return ssNnsynsub[1563];
    }
    double GetssNnsynsub1564(){
     return ssNnsynsub[1564];
    }
    double GetssNnsynsub1565(){
     return ssNnsynsub[1565];
    }
    double GetssNnsynsub1566(){
     return ssNnsynsub[1566];
    }
    double GetssNnsynsub1567(){
     return ssNnsynsub[1567];
    }
    double GetssNnsynsub1568(){
     return ssNnsynsub[1568];
    }
    double GetssNnsynsub1569(){
     return ssNnsynsub[1569];
    }
    double GetssNnsynsub1570(){
     return ssNnsynsub[1570];
    }
    double GetssNnsynsub1571(){
     return ssNnsynsub[1571];
    }
    double GetssNnsynsub1572(){
     return ssNnsynsub[1572];
    }
    double GetssNnsynsub1573(){
     return ssNnsynsub[1573];
    }
    double GetssNnsynsub1574(){
     return ssNnsynsub[1574];
    }
    double GetssNnsynsub1575(){
     return ssNnsynsub[1575];
    }
    double GetssNnsynsub1576(){
     return ssNnsynsub[1576];
    }
    double GetssNnsynsub1577(){
     return ssNnsynsub[1577];
    }
    double GetssNnsynsub1578(){
     return ssNnsynsub[1578];
    }
    double GetssNnsynsub1579(){
     return ssNnsynsub[1579];
    }
    double GetssNnsynsub1580(){
     return ssNnsynsub[1580];
    }
    double GetssNnsynsub1581(){
     return ssNnsynsub[1581];
    }
    double GetssNnsynsub1582(){
     return ssNnsynsub[1582];
    }
    double GetssNnsynsub1583(){
     return ssNnsynsub[1583];
    }
    double GetssNnsynsub1584(){
     return ssNnsynsub[1584];
    }
    double GetssNnsynsub1585(){
     return ssNnsynsub[1585];
    }
    double GetssNnsynsub1586(){
     return ssNnsynsub[1586];
    }
    double GetssNnsynsub1587(){
     return ssNnsynsub[1587];
    }
    double GetssNnsynsub1588(){
     return ssNnsynsub[1588];
    }
    double GetssNnsynsub1589(){
     return ssNnsynsub[1589];
    }
    double GetssNnsynsub1590(){
     return ssNnsynsub[1590];
    }
    double GetssNnsynsub1591(){
     return ssNnsynsub[1591];
    }
    double GetssNnsynsub1592(){
     return ssNnsynsub[1592];
    }
    double GetssNnsynsub1593(){
     return ssNnsynsub[1593];
    }
    double GetssNnsynsub1594(){
     return ssNnsynsub[1594];
    }
    double GetssNnsynsub1595(){
     return ssNnsynsub[1595];
    }
    double GetssNnsynsub1596(){
     return ssNnsynsub[1596];
    }
    double GetssNnsynsub1597(){
     return ssNnsynsub[1597];
    }
    double GetssNnsynsub1598(){
     return ssNnsynsub[1598];
    }
    double GetssNnsynsub1599(){
     return ssNnsynsub[1599];
    }
    double GetssNnsynsub1600(){
     return ssNnsynsub[1600];
    }
    double GetssNnsynsub1601(){
     return ssNnsynsub[1601];
    }
    double GetssNnsynsub1602(){
     return ssNnsynsub[1602];
    }
    double GetssNnsynsub1603(){
     return ssNnsynsub[1603];
    }
    double GetssNnsynsub1604(){
     return ssNnsynsub[1604];
    }
    double GetssNnsynsub1605(){
     return ssNnsynsub[1605];
    }
    double GetssNnsynsub1606(){
     return ssNnsynsub[1606];
    }
    double GetssNnsynsub1607(){
     return ssNnsynsub[1607];
    }
    double GetssNnsynsub1608(){
     return ssNnsynsub[1608];
    }
    double GetssNnsynsub1609(){
     return ssNnsynsub[1609];
    }
    double GetssNnsynsub1610(){
     return ssNnsynsub[1610];
    }
    double GetssNnsynsub1611(){
     return ssNnsynsub[1611];
    }
    double GetssNnsynsub1612(){
     return ssNnsynsub[1612];
    }
    double GetssNnsynsub1613(){
     return ssNnsynsub[1613];
    }
    double GetssNnsynsub1614(){
     return ssNnsynsub[1614];
    }
    double GetssNnsynsub1615(){
     return ssNnsynsub[1615];
    }
    double GetssNnsynsub1616(){
     return ssNnsynsub[1616];
    }
    double GetssNnsynsub1617(){
     return ssNnsynsub[1617];
    }
    double GetssNnsynsub1618(){
     return ssNnsynsub[1618];
    }
    double GetssNnsynsub1619(){
     return ssNnsynsub[1619];
    }
    double GetssNnsynsub1620(){
     return ssNnsynsub[1620];
    }
    double GetssNnsynsub1621(){
     return ssNnsynsub[1621];
    }
    double GetssNnsynsub1622(){
     return ssNnsynsub[1622];
    }
    double GetssNnsynsub1623(){
     return ssNnsynsub[1623];
    }
    double GetssNnsynsub1624(){
     return ssNnsynsub[1624];
    }
    double GetssNnsynsub1625(){
     return ssNnsynsub[1625];
    }
    double GetssNnsynsub1626(){
     return ssNnsynsub[1626];
    }
    double GetssNnsynsub1627(){
     return ssNnsynsub[1627];
    }
    double GetssNnsynsub1628(){
     return ssNnsynsub[1628];
    }
    double GetssNnsynsub1629(){
     return ssNnsynsub[1629];
    }
    double GetssNnsynsub1630(){
     return ssNnsynsub[1630];
    }
    double GetssNnsynsub1631(){
     return ssNnsynsub[1631];
    }
    double GetssNnsynsub1632(){
     return ssNnsynsub[1632];
    }
    double GetssNnsynsub1633(){
     return ssNnsynsub[1633];
    }
    double GetssNnsynsub1634(){
     return ssNnsynsub[1634];
    }
    double GetssNnsynsub1635(){
     return ssNnsynsub[1635];
    }
    double GetssNnsynsub1636(){
     return ssNnsynsub[1636];
    }
    double GetssNnsynsub1637(){
     return ssNnsynsub[1637];
    }
    double GetssNnsynsub1638(){
     return ssNnsynsub[1638];
    }
    double GetssNnsynsub1639(){
     return ssNnsynsub[1639];
    }
    double GetssNnsynsub1640(){
     return ssNnsynsub[1640];
    }
    double GetssNnsynsub1641(){
     return ssNnsynsub[1641];
    }
    double GetssNnsynsub1642(){
     return ssNnsynsub[1642];
    }
    double GetssNnsynsub1643(){
     return ssNnsynsub[1643];
    }
    double GetssNnsynsub1644(){
     return ssNnsynsub[1644];
    }
    double GetssNnsynsub1645(){
     return ssNnsynsub[1645];
    }
    double GetssNnsynsub1646(){
     return ssNnsynsub[1646];
    }
    double GetssNnsynsub1647(){
     return ssNnsynsub[1647];
    }
    double GetssNnsynsub1648(){
     return ssNnsynsub[1648];
    }
    double GetssNnsynsub1649(){
     return ssNnsynsub[1649];
    }
    double GetssNnsynsub1650(){
     return ssNnsynsub[1650];
    }
    double GetssNnsynsub1651(){
     return ssNnsynsub[1651];
    }
    double GetssNnsynsub1652(){
     return ssNnsynsub[1652];
    }
    double GetssNnsynsub1653(){
     return ssNnsynsub[1653];
    }
    double GetssNnsynsub1654(){
     return ssNnsynsub[1654];
    }
    double GetssNnsynsub1655(){
     return ssNnsynsub[1655];
    }
    double GetssNnsynsub1656(){
     return ssNnsynsub[1656];
    }
    double GetssNnsynsub1657(){
     return ssNnsynsub[1657];
    }
    double GetssNnsynsub1658(){
     return ssNnsynsub[1658];
    }
    double GetssNnsynsub1659(){
     return ssNnsynsub[1659];
    }
    double GetssNnsynsub1660(){
     return ssNnsynsub[1660];
    }
    double GetssNnsynsub1661(){
     return ssNnsynsub[1661];
    }
    double GetssNnsynsub1662(){
     return ssNnsynsub[1662];
    }
    double GetssNnsynsub1663(){
     return ssNnsynsub[1663];
    }
    double GetssNnsynsub1664(){
     return ssNnsynsub[1664];
    }
    double GetssNnsynsub1665(){
     return ssNnsynsub[1665];
    }
    double GetssNnsynsub1666(){
     return ssNnsynsub[1666];
    }
    double GetssNnsynsub1667(){
     return ssNnsynsub[1667];
    }
    double GetssNnsynsub1668(){
     return ssNnsynsub[1668];
    }
    double GetssNnsynsub1669(){
     return ssNnsynsub[1669];
    }
    double GetssNnsynsub1670(){
     return ssNnsynsub[1670];
    }
    double GetssNnsynsub1671(){
     return ssNnsynsub[1671];
    }
    double GetssNnsynsub1672(){
     return ssNnsynsub[1672];
    }
    double GetssNnsynsub1673(){
     return ssNnsynsub[1673];
    }
    double GetssNnsynsub1674(){
     return ssNnsynsub[1674];
    }
    double GetssNnsynsub1675(){
     return ssNnsynsub[1675];
    }
    double GetssNnsynsub1676(){
     return ssNnsynsub[1676];
    }
    double GetssNnsynsub1677(){
     return ssNnsynsub[1677];
    }
    double GetssNnsynsub1678(){
     return ssNnsynsub[1678];
    }
    double GetssNnsynsub1679(){
     return ssNnsynsub[1679];
    }
    double GetssNnsynsub1680(){
     return ssNnsynsub[1680];
    }
    double GetssNnsynsub1681(){
     return ssNnsynsub[1681];
    }
    double GetssNnsynsub1682(){
     return ssNnsynsub[1682];
    }
    double GetssNnsynsub1683(){
     return ssNnsynsub[1683];
    }
    double GetssNnsynsub1684(){
     return ssNnsynsub[1684];
    }
    double GetssNnsynsub1685(){
     return ssNnsynsub[1685];
    }
    double GetssNnsynsub1686(){
     return ssNnsynsub[1686];
    }
    double GetssNnsynsub1687(){
     return ssNnsynsub[1687];
    }
    double GetssNnsynsub1688(){
     return ssNnsynsub[1688];
    }
    double GetssNnsynsub1689(){
     return ssNnsynsub[1689];
    }
    double GetssNnsynsub1690(){
     return ssNnsynsub[1690];
    }
    double GetssNnsynsub1691(){
     return ssNnsynsub[1691];
    }
    double GetssNnsynsub1692(){
     return ssNnsynsub[1692];
    }
    double GetssNnsynsub1693(){
     return ssNnsynsub[1693];
    }
    double GetssNnsynsub1694(){
     return ssNnsynsub[1694];
    }
    double GetssNnsynsub1695(){
     return ssNnsynsub[1695];
    }
    double GetssNnsynsub1696(){
     return ssNnsynsub[1696];
    }
    double GetssNnsynsub1697(){
     return ssNnsynsub[1697];
    }
    double GetssNnsynsub1698(){
     return ssNnsynsub[1698];
    }
    double GetssNnsynsub1699(){
     return ssNnsynsub[1699];
    }
    double GetssNnsynsub1700(){
     return ssNnsynsub[1700];
    }
    double GetssNnsynsub1701(){
     return ssNnsynsub[1701];
    }
    double GetssNnsynsub1702(){
     return ssNnsynsub[1702];
    }
    double GetssNnsynsub1703(){
     return ssNnsynsub[1703];
    }
    double GetssNnsynsub1704(){
     return ssNnsynsub[1704];
    }
    double GetssNnsynsub1705(){
     return ssNnsynsub[1705];
    }
    double GetssNnsynsub1706(){
     return ssNnsynsub[1706];
    }
    double GetssNnsynsub1707(){
     return ssNnsynsub[1707];
    }
    double GetssNnsynsub1708(){
     return ssNnsynsub[1708];
    }
    double GetssNnsynsub1709(){
     return ssNnsynsub[1709];
    }
    double GetssNnsynsub1710(){
     return ssNnsynsub[1710];
    }
    double GetssNnsynsub1711(){
     return ssNnsynsub[1711];
    }
    double GetssNnsynsub1712(){
     return ssNnsynsub[1712];
    }
    double GetssNnsynsub1713(){
     return ssNnsynsub[1713];
    }
    double GetssNnsynsub1714(){
     return ssNnsynsub[1714];
    }
    double GetssNnsynsub1715(){
     return ssNnsynsub[1715];
    }
    double GetssNnsynsub1716(){
     return ssNnsynsub[1716];
    }
    double GetssNnsynsub1717(){
     return ssNnsynsub[1717];
    }
    double GetssNnsynsub1718(){
     return ssNnsynsub[1718];
    }
    double GetssNnsynsub1719(){
     return ssNnsynsub[1719];
    }
    double GetssNnsynsub1720(){
     return ssNnsynsub[1720];
    }
    double GetssNnsynsub1721(){
     return ssNnsynsub[1721];
    }
    double GetssNnsynsub1722(){
     return ssNnsynsub[1722];
    }
    double GetssNnsynsub1723(){
     return ssNnsynsub[1723];
    }
    double GetssNnsynsub1724(){
     return ssNnsynsub[1724];
    }
    double GetssNnsynsub1725(){
     return ssNnsynsub[1725];
    }
    double GetssNnsynsub1726(){
     return ssNnsynsub[1726];
    }
    double GetssNnsynsub1727(){
     return ssNnsynsub[1727];
    }
    double GetssNnsynsub1728(){
     return ssNnsynsub[1728];
    }
    double GetssNnsynsub1729(){
     return ssNnsynsub[1729];
    }
    double GetssNnsynsub1730(){
     return ssNnsynsub[1730];
    }
    double GetssNnsynsub1731(){
     return ssNnsynsub[1731];
    }
    double GetssNnsynsub1732(){
     return ssNnsynsub[1732];
    }
    double GetssNnsynsub1733(){
     return ssNnsynsub[1733];
    }
    double GetssNnsynsub1734(){
     return ssNnsynsub[1734];
    }
    double GetssNnsynsub1735(){
     return ssNnsynsub[1735];
    }
    double GetssNnsynsub1736(){
     return ssNnsynsub[1736];
    }
    double GetssNnsynsub1737(){
     return ssNnsynsub[1737];
    }
    double GetssNnsynsub1738(){
     return ssNnsynsub[1738];
    }
    double GetssNnsynsub1739(){
     return ssNnsynsub[1739];
    }
    double GetssNnsynsub1740(){
     return ssNnsynsub[1740];
    }
    double GetssNnsynsub1741(){
     return ssNnsynsub[1741];
    }
    double GetssNnsynsub1742(){
     return ssNnsynsub[1742];
    }
    double GetssNnsynsub1743(){
     return ssNnsynsub[1743];
    }
    double GetssNnsynsub1744(){
     return ssNnsynsub[1744];
    }
    double GetssNnsynsub1745(){
     return ssNnsynsub[1745];
    }
    double GetssNnsynsub1746(){
     return ssNnsynsub[1746];
    }
    double GetssNnsynsub1747(){
     return ssNnsynsub[1747];
    }
    double GetssNnsynsub1748(){
     return ssNnsynsub[1748];
    }
    double GetssNnsynsub1749(){
     return ssNnsynsub[1749];
    }
    double GetssNnsynsub1750(){
     return ssNnsynsub[1750];
    }
    double GetssNnsynsub1751(){
     return ssNnsynsub[1751];
    }
    double GetssNnsynsub1752(){
     return ssNnsynsub[1752];
    }
    double GetssNnsynsub1753(){
     return ssNnsynsub[1753];
    }
    double GetssNnsynsub1754(){
     return ssNnsynsub[1754];
    }
    double GetssNnsynsub1755(){
     return ssNnsynsub[1755];
    }
    double GetssNnsynsub1756(){
     return ssNnsynsub[1756];
    }
    double GetssNnsynsub1757(){
     return ssNnsynsub[1757];
    }
    double GetssNnsynsub1758(){
     return ssNnsynsub[1758];
    }
    double GetssNnsynsub1759(){
     return ssNnsynsub[1759];
    }
    double GetssNnsynsub1760(){
     return ssNnsynsub[1760];
    }
    double GetssNnsynsub1761(){
     return ssNnsynsub[1761];
    }
    double GetssNnsynsub1762(){
     return ssNnsynsub[1762];
    }
    double GetssNnsynsub1763(){
     return ssNnsynsub[1763];
    }
    double GetssNnsynsub1764(){
     return ssNnsynsub[1764];
    }
    double GetssNnsynsub1765(){
     return ssNnsynsub[1765];
    }
    double GetssNnsynsub1766(){
     return ssNnsynsub[1766];
    }
    double GetssNnsynsub1767(){
     return ssNnsynsub[1767];
    }
    double GetssNnsynsub1768(){
     return ssNnsynsub[1768];
    }
    double GetssNnsynsub1769(){
     return ssNnsynsub[1769];
    }
    double GetssNnsynsub1770(){
     return ssNnsynsub[1770];
    }
    double GetssNnsynsub1771(){
     return ssNnsynsub[1771];
    }
    double GetssNnsynsub1772(){
     return ssNnsynsub[1772];
    }
    double GetssNnsynsub1773(){
     return ssNnsynsub[1773];
    }
    double GetssNnsynsub1774(){
     return ssNnsynsub[1774];
    }
    double GetssNnsynsub1775(){
     return ssNnsynsub[1775];
    }
    double GetssNnsynsub1776(){
     return ssNnsynsub[1776];
    }
    double GetssNnsynsub1777(){
     return ssNnsynsub[1777];
    }
    double GetssNnsynsub1778(){
     return ssNnsynsub[1778];
    }
    double GetssNnsynsub1779(){
     return ssNnsynsub[1779];
    }
    double GetssNnsynsub1780(){
     return ssNnsynsub[1780];
    }
    double GetssNnsynsub1781(){
     return ssNnsynsub[1781];
    }
    double GetssNnsynsub1782(){
     return ssNnsynsub[1782];
    }
    double GetssNnsynsub1783(){
     return ssNnsynsub[1783];
    }
    double GetssNnsynsub1784(){
     return ssNnsynsub[1784];
    }
    double GetssNnsynsub1785(){
     return ssNnsynsub[1785];
    }
    double GetssNnsynsub1786(){
     return ssNnsynsub[1786];
    }
    double GetssNnsynsub1787(){
     return ssNnsynsub[1787];
    }
    double GetssNnsynsub1788(){
     return ssNnsynsub[1788];
    }
    double GetssNnsynsub1789(){
     return ssNnsynsub[1789];
    }
    double GetssNnsynsub1790(){
     return ssNnsynsub[1790];
    }
    double GetssNnsynsub1791(){
     return ssNnsynsub[1791];
    }
    double GetssNnsynsub1792(){
     return ssNnsynsub[1792];
    }
    double GetssNnsynsub1793(){
     return ssNnsynsub[1793];
    }
    double GetssNnsynsub1794(){
     return ssNnsynsub[1794];
    }
    double GetssNnsynsub1795(){
     return ssNnsynsub[1795];
    }
    double GetssNnsynsub1796(){
     return ssNnsynsub[1796];
    }
    double GetssNnsynsub1797(){
     return ssNnsynsub[1797];
    }
    double GetssNnsynsub1798(){
     return ssNnsynsub[1798];
    }
    double GetssNnsynsub1799(){
     return ssNnsynsub[1799];
    }
    double GetssNnsynsub1800(){
     return ssNnsynsub[1800];
    }
    double GetssNnsynsub1801(){
     return ssNnsynsub[1801];
    }
    double GetssNnsynsub1802(){
     return ssNnsynsub[1802];
    }
    double GetssNnsynsub1803(){
     return ssNnsynsub[1803];
    }
    double GetssNnsynsub1804(){
     return ssNnsynsub[1804];
    }
    double GetssNnsynsub1805(){
     return ssNnsynsub[1805];
    }
    double GetssNnsynsub1806(){
     return ssNnsynsub[1806];
    }
    double GetssNnsynsub1807(){
     return ssNnsynsub[1807];
    }
    double GetssNnsynsub1808(){
     return ssNnsynsub[1808];
    }
    double GetssNnsynsub1809(){
     return ssNnsynsub[1809];
    }
    double GetssNnsynsub1810(){
     return ssNnsynsub[1810];
    }
    double GetssNnsynsub1811(){
     return ssNnsynsub[1811];
    }
    double GetssNnsynsub1812(){
     return ssNnsynsub[1812];
    }
    double GetssNnsynsub1813(){
     return ssNnsynsub[1813];
    }
    double GetssNnsynsub1814(){
     return ssNnsynsub[1814];
    }
    double GetssNnsynsub1815(){
     return ssNnsynsub[1815];
    }
    double GetssNnsynsub1816(){
     return ssNnsynsub[1816];
    }
    double GetssNnsynsub1817(){
     return ssNnsynsub[1817];
    }
    double GetssNnsynsub1818(){
     return ssNnsynsub[1818];
    }
    double GetssNnsynsub1819(){
     return ssNnsynsub[1819];
    }
    double GetssNnsynsub1820(){
     return ssNnsynsub[1820];
    }
    double GetssNnsynsub1821(){
     return ssNnsynsub[1821];
    }
    double GetssNnsynsub1822(){
     return ssNnsynsub[1822];
    }
    double GetssNnsynsub1823(){
     return ssNnsynsub[1823];
    }
    double GetssNnsynsub1824(){
     return ssNnsynsub[1824];
    }
    double GetssNnsynsub1825(){
     return ssNnsynsub[1825];
    }
    double GetssNnsynsub1826(){
     return ssNnsynsub[1826];
    }
    double GetssNnsynsub1827(){
     return ssNnsynsub[1827];
    }
    double GetssNnsynsub1828(){
     return ssNnsynsub[1828];
    }
    double GetssNnsynsub1829(){
     return ssNnsynsub[1829];
    }
    double GetssNnsynsub1830(){
     return ssNnsynsub[1830];
    }
    double GetssNnsynsub1831(){
     return ssNnsynsub[1831];
    }
    double GetssNnsynsub1832(){
     return ssNnsynsub[1832];
    }
    double GetssNnsynsub1833(){
     return ssNnsynsub[1833];
    }
    double GetssNnsynsub1834(){
     return ssNnsynsub[1834];
    }
    double GetssNnsynsub1835(){
     return ssNnsynsub[1835];
    }
    double GetssNnsynsub1836(){
     return ssNnsynsub[1836];
    }
    double GetssNnsynsub1837(){
     return ssNnsynsub[1837];
    }
    double GetssNnsynsub1838(){
     return ssNnsynsub[1838];
    }
    double GetssNnsynsub1839(){
     return ssNnsynsub[1839];
    }
    double GetssNnsynsub1840(){
     return ssNnsynsub[1840];
    }
    double GetssNnsynsub1841(){
     return ssNnsynsub[1841];
    }
    double GetssNnsynsub1842(){
     return ssNnsynsub[1842];
    }
    double GetssNnsynsub1843(){
     return ssNnsynsub[1843];
    }
    double GetssNnsynsub1844(){
     return ssNnsynsub[1844];
    }
    double GetssNnsynsub1845(){
     return ssNnsynsub[1845];
    }
    double GetssNnsynsub1846(){
     return ssNnsynsub[1846];
    }
    double GetssNnsynsub1847(){
     return ssNnsynsub[1847];
    }
    double GetssNnsynsub1848(){
     return ssNnsynsub[1848];
    }
    double GetssNnsynsub1849(){
     return ssNnsynsub[1849];
    }
    double GetssNnsynsub1850(){
     return ssNnsynsub[1850];
    }
    double GetssNnsynsub1851(){
     return ssNnsynsub[1851];
    }
    double GetssNnsynsub1852(){
     return ssNnsynsub[1852];
    }
    double GetssNnsynsub1853(){
     return ssNnsynsub[1853];
    }
    double GetssNnsynsub1854(){
     return ssNnsynsub[1854];
    }
    double GetssNnsynsub1855(){
     return ssNnsynsub[1855];
    }
    double GetssNnsynsub1856(){
     return ssNnsynsub[1856];
    }
    double GetssNnsynsub1857(){
     return ssNnsynsub[1857];
    }
    double GetssNnsynsub1858(){
     return ssNnsynsub[1858];
    }
    double GetssNnsynsub1859(){
     return ssNnsynsub[1859];
    }
    double GetssNnsynsub1860(){
     return ssNnsynsub[1860];
    }
    double GetssNnsynsub1861(){
     return ssNnsynsub[1861];
    }
    double GetssNnsynsub1862(){
     return ssNnsynsub[1862];
    }
    double GetssNnsynsub1863(){
     return ssNnsynsub[1863];
    }
    double GetssNnsynsub1864(){
     return ssNnsynsub[1864];
    }
    double GetssNnsynsub1865(){
     return ssNnsynsub[1865];
    }
    double GetssNnsynsub1866(){
     return ssNnsynsub[1866];
    }
    double GetssNnsynsub1867(){
     return ssNnsynsub[1867];
    }
    double GetssNnsynsub1868(){
     return ssNnsynsub[1868];
    }
    double GetssNnsynsub1869(){
     return ssNnsynsub[1869];
    }
    double GetssNnsynsub1870(){
     return ssNnsynsub[1870];
    }
    double GetssNnsynsub1871(){
     return ssNnsynsub[1871];
    }
    double GetssNnsynsub1872(){
     return ssNnsynsub[1872];
    }
    double GetssNnsynsub1873(){
     return ssNnsynsub[1873];
    }
    double GetssNnsynsub1874(){
     return ssNnsynsub[1874];
    }
    double GetssNnsynsub1875(){
     return ssNnsynsub[1875];
    }
    double GetssNnsynsub1876(){
     return ssNnsynsub[1876];
    }
    double GetssNnsynsub1877(){
     return ssNnsynsub[1877];
    }
    double GetssNnsynsub1878(){
     return ssNnsynsub[1878];
    }
    double GetssNnsynsub1879(){
     return ssNnsynsub[1879];
    }
    double GetssNnsynsub1880(){
     return ssNnsynsub[1880];
    }
    double GetssNnsynsub1881(){
     return ssNnsynsub[1881];
    }
    double GetssNnsynsub1882(){
     return ssNnsynsub[1882];
    }
    double GetssNnsynsub1883(){
     return ssNnsynsub[1883];
    }
    double GetssNnsynsub1884(){
     return ssNnsynsub[1884];
    }
    double GetssNnsynsub1885(){
     return ssNnsynsub[1885];
    }
    double GetssNnsynsub1886(){
     return ssNnsynsub[1886];
    }
    double GetssNnsynsub1887(){
     return ssNnsynsub[1887];
    }
    double GetssNnsynsub1888(){
     return ssNnsynsub[1888];
    }
    double GetssNnsynsub1889(){
     return ssNnsynsub[1889];
    }
    double GetssNnsynsub1890(){
     return ssNnsynsub[1890];
    }
    double GetssNnsynsub1891(){
     return ssNnsynsub[1891];
    }
    double GetssNnsynsub1892(){
     return ssNnsynsub[1892];
    }
    double GetssNnsynsub1893(){
     return ssNnsynsub[1893];
    }
    double GetssNnsynsub1894(){
     return ssNnsynsub[1894];
    }
    double GetssNnsynsub1895(){
     return ssNnsynsub[1895];
    }
    double GetssNnsynsub1896(){
     return ssNnsynsub[1896];
    }
    double GetssNnsynsub1897(){
     return ssNnsynsub[1897];
    }
    double GetssNnsynsub1898(){
     return ssNnsynsub[1898];
    }
    double GetssNnsynsub1899(){
     return ssNnsynsub[1899];
    }
    double GetssNnsynsub1900(){
     return ssNnsynsub[1900];
    }
    double GetssNnsynsub1901(){
     return ssNnsynsub[1901];
    }
    double GetssNnsynsub1902(){
     return ssNnsynsub[1902];
    }
    double GetssNnsynsub1903(){
     return ssNnsynsub[1903];
    }
    double GetssNnsynsub1904(){
     return ssNnsynsub[1904];
    }
    double GetssNnsynsub1905(){
     return ssNnsynsub[1905];
    }
    double GetssNnsynsub1906(){
     return ssNnsynsub[1906];
    }
    double GetssNnsynsub1907(){
     return ssNnsynsub[1907];
    }
    double GetssNnsynsub1908(){
     return ssNnsynsub[1908];
    }
    double GetssNnsynsub1909(){
     return ssNnsynsub[1909];
    }
    double GetssNnsynsub1910(){
     return ssNnsynsub[1910];
    }
    double GetssNnsynsub1911(){
     return ssNnsynsub[1911];
    }
    double GetssNnsynsub1912(){
     return ssNnsynsub[1912];
    }
    double GetssNnsynsub1913(){
     return ssNnsynsub[1913];
    }
    double GetssNnsynsub1914(){
     return ssNnsynsub[1914];
    }
    double GetssNnsynsub1915(){
     return ssNnsynsub[1915];
    }
    double GetssNnsynsub1916(){
     return ssNnsynsub[1916];
    }
    double GetssNnsynsub1917(){
     return ssNnsynsub[1917];
    }
    double GetssNnsynsub1918(){
     return ssNnsynsub[1918];
    }
    double GetssNnsynsub1919(){
     return ssNnsynsub[1919];
    }
    double GetssNnsynsub1920(){
     return ssNnsynsub[1920];
    }
    double GetssNnsynsub1921(){
     return ssNnsynsub[1921];
    }
    double GetssNnsynsub1922(){
     return ssNnsynsub[1922];
    }
    double GetssNnsynsub1923(){
     return ssNnsynsub[1923];
    }
    double GetssNnsynsub1924(){
     return ssNnsynsub[1924];
    }
    double GetssNnsynsub1925(){
     return ssNnsynsub[1925];
    }
    double GetssNnsynsub1926(){
     return ssNnsynsub[1926];
    }
    double GetssNnsynsub1927(){
     return ssNnsynsub[1927];
    }
    double GetssNnsynsub1928(){
     return ssNnsynsub[1928];
    }
    double GetssNnsynsub1929(){
     return ssNnsynsub[1929];
    }
    double GetssNnsynsub1930(){
     return ssNnsynsub[1930];
    }
    double GetssNnsynsub1931(){
     return ssNnsynsub[1931];
    }
    double GetssNnsynsub1932(){
     return ssNnsynsub[1932];
    }
    double GetssNnsynsub1933(){
     return ssNnsynsub[1933];
    }
    double GetssNnsynsub1934(){
     return ssNnsynsub[1934];
    }
    double GetssNnsynsub1935(){
     return ssNnsynsub[1935];
    }
    double GetssNnsynsub1936(){
     return ssNnsynsub[1936];
    }
    double GetssNnsynsub1937(){
     return ssNnsynsub[1937];
    }
    double GetssNnsynsub1938(){
     return ssNnsynsub[1938];
    }
    double GetssNnsynsub1939(){
     return ssNnsynsub[1939];
    }
    double GetssNnsynsub1940(){
     return ssNnsynsub[1940];
    }
    double GetssNnsynsub1941(){
     return ssNnsynsub[1941];
    }
    double GetssNnsynsub1942(){
     return ssNnsynsub[1942];
    }
    double GetssNnsynsub1943(){
     return ssNnsynsub[1943];
    }
    double GetssNnsynsub1944(){
     return ssNnsynsub[1944];
    }
    double GetssNnsynsub1945(){
     return ssNnsynsub[1945];
    }
    double GetssNnsynsub1946(){
     return ssNnsynsub[1946];
    }
    double GetssNnsynsub1947(){
     return ssNnsynsub[1947];
    }
    double GetssNnsynsub1948(){
     return ssNnsynsub[1948];
    }
    double GetssNnsynsub1949(){
     return ssNnsynsub[1949];
    }
    double GetssNnsynsub1950(){
     return ssNnsynsub[1950];
    }
    double GetssNnsynsub1951(){
     return ssNnsynsub[1951];
    }
    double GetssNnsynsub1952(){
     return ssNnsynsub[1952];
    }
    double GetssNnsynsub1953(){
     return ssNnsynsub[1953];
    }
    double GetssNnsynsub1954(){
     return ssNnsynsub[1954];
    }
    double GetssNnsynsub1955(){
     return ssNnsynsub[1955];
    }
    double GetssNnsynsub1956(){
     return ssNnsynsub[1956];
    }
    double GetssNnsynsub1957(){
     return ssNnsynsub[1957];
    }
    double GetssNnsynsub1958(){
     return ssNnsynsub[1958];
    }
    double GetssNnsynsub1959(){
     return ssNnsynsub[1959];
    }
    double GetssNnsynsub1960(){
     return ssNnsynsub[1960];
    }
    double GetssNnsynsub1961(){
     return ssNnsynsub[1961];
    }
    double GetssNnsynsub1962(){
     return ssNnsynsub[1962];
    }
    double GetssNnsynsub1963(){
     return ssNnsynsub[1963];
    }
    double GetssNnsynsub1964(){
     return ssNnsynsub[1964];
    }
    double GetssNnsynsub1965(){
     return ssNnsynsub[1965];
    }
    double GetssNnsynsub1966(){
     return ssNnsynsub[1966];
    }
    double GetssNnsynsub1967(){
     return ssNnsynsub[1967];
    }
    double GetssNnsynsub1968(){
     return ssNnsynsub[1968];
    }
    double GetssNnsynsub1969(){
     return ssNnsynsub[1969];
    }
    double GetssNnsynsub1970(){
     return ssNnsynsub[1970];
    }
    double GetssNnsynsub1971(){
     return ssNnsynsub[1971];
    }
    double GetssNnsynsub1972(){
     return ssNnsynsub[1972];
    }
    double GetssNnsynsub1973(){
     return ssNnsynsub[1973];
    }
    double GetssNnsynsub1974(){
     return ssNnsynsub[1974];
    }
    double GetssNnsynsub1975(){
     return ssNnsynsub[1975];
    }
    double GetssNnsynsub1976(){
     return ssNnsynsub[1976];
    }
    double GetssNnsynsub1977(){
     return ssNnsynsub[1977];
    }
    double GetssNnsynsub1978(){
     return ssNnsynsub[1978];
    }
    double GetssNnsynsub1979(){
     return ssNnsynsub[1979];
    }
    double GetssNnsynsub1980(){
     return ssNnsynsub[1980];
    }
    double GetssNnsynsub1981(){
     return ssNnsynsub[1981];
    }
    double GetssNnsynsub1982(){
     return ssNnsynsub[1982];
    }
    double GetssNnsynsub1983(){
     return ssNnsynsub[1983];
    }
    double GetssNnsynsub1984(){
     return ssNnsynsub[1984];
    }
    double GetssNnsynsub1985(){
     return ssNnsynsub[1985];
    }
    double GetssNnsynsub1986(){
     return ssNnsynsub[1986];
    }
    double GetssNnsynsub1987(){
     return ssNnsynsub[1987];
    }
    double GetssNnsynsub1988(){
     return ssNnsynsub[1988];
    }
    double GetssNnsynsub1989(){
     return ssNnsynsub[1989];
    }
    double GetssNnsynsub1990(){
     return ssNnsynsub[1990];
    }
    double GetssNnsynsub1991(){
     return ssNnsynsub[1991];
    }
    double GetssNnsynsub1992(){
     return ssNnsynsub[1992];
    }
    double GetssNnsynsub1993(){
     return ssNnsynsub[1993];
    }
    double GetssNnsynsub1994(){
     return ssNnsynsub[1994];
    }
    double GetssNnsynsub1995(){
     return ssNnsynsub[1995];
    }
    double GetssNnsynsub1996(){
     return ssNnsynsub[1996];
    }
    double GetssNnsynsub1997(){
     return ssNnsynsub[1997];
    }
    double GetssNnsynsub1998(){
     return ssNnsynsub[1998];
    }
    double GetssNnsynsub1999(){
     return ssNnsynsub[1999];
    }
    double GetssNnsynsub2000(){
     return ssNnsynsub[2000];
    }
    double GetssNnsynsub2001(){
     return ssNnsynsub[2001];
    }
    double GetssNnsynsub2002(){
     return ssNnsynsub[2002];
    }
    double GetssNnsynsub2003(){
     return ssNnsynsub[2003];
    }
    double GetssNnsynsub2004(){
     return ssNnsynsub[2004];
    }
    double GetssNnsynsub2005(){
     return ssNnsynsub[2005];
    }
    double GetssNnsynsub2006(){
     return ssNnsynsub[2006];
    }
    double GetssNnsynsub2007(){
     return ssNnsynsub[2007];
    }
    double GetssNnsynsub2008(){
     return ssNnsynsub[2008];
    }
    double GetssNnsynsub2009(){
     return ssNnsynsub[2009];
    }
    double GetssNnsynsub2010(){
     return ssNnsynsub[2010];
    }
    double GetssNnsynsub2011(){
     return ssNnsynsub[2011];
    }
    double GetssNnsynsub2012(){
     return ssNnsynsub[2012];
    }
    double GetssNnsynsub2013(){
     return ssNnsynsub[2013];
    }
    double GetssNnsynsub2014(){
     return ssNnsynsub[2014];
    }
    double GetssNnsynsub2015(){
     return ssNnsynsub[2015];
    }
    double GetssNnsynsub2016(){
     return ssNnsynsub[2016];
    }
    double GetssNnsynsub2017(){
     return ssNnsynsub[2017];
    }
    double GetssNnsynsub2018(){
     return ssNnsynsub[2018];
    }
    double GetssNnsynsub2019(){
     return ssNnsynsub[2019];
    }
    double GetssNnsynsub2020(){
     return ssNnsynsub[2020];
    }
    double GetssNnsynsub2021(){
     return ssNnsynsub[2021];
    }
    double GetssNnsynsub2022(){
     return ssNnsynsub[2022];
    }
    double GetssNnsynsub2023(){
     return ssNnsynsub[2023];
    }
    double GetssNnsynsub2024(){
     return ssNnsynsub[2024];
    }
    double GetssNnsynsub2025(){
     return ssNnsynsub[2025];
    }
    double GetssNnsynsub2026(){
     return ssNnsynsub[2026];
    }
    double GetssNnsynsub2027(){
     return ssNnsynsub[2027];
    }
    double GetssNnsynsub2028(){
     return ssNnsynsub[2028];
    }
    double GetssNnsynsub2029(){
     return ssNnsynsub[2029];
    }
    double GetssNnsynsub2030(){
     return ssNnsynsub[2030];
    }
    double GetssNnsynsub2031(){
     return ssNnsynsub[2031];
    }
    double GetssNnsynsub2032(){
     return ssNnsynsub[2032];
    }
    double GetssNnsynsub2033(){
     return ssNnsynsub[2033];
    }
    double GetssNnsynsub2034(){
     return ssNnsynsub[2034];
    }
    double GetssNnsynsub2035(){
     return ssNnsynsub[2035];
    }
    double GetssNnsynsub2036(){
     return ssNnsynsub[2036];
    }
    double GetssNnsynsub2037(){
     return ssNnsynsub[2037];
    }
    double GetssNnsynsub2038(){
     return ssNnsynsub[2038];
    }
    double GetssNnsynsub2039(){
     return ssNnsynsub[2039];
    }
    double GetssNnsynsub2040(){
     return ssNnsynsub[2040];
    }
    double GetssNnsynsub2041(){
     return ssNnsynsub[2041];
    }
    double GetssNnsynsub2042(){
     return ssNnsynsub[2042];
    }
    double GetssNnsynsub2043(){
     return ssNnsynsub[2043];
    }
    double GetssNnsynsub2044(){
     return ssNnsynsub[2044];
    }
    double GetssNnsynsub2045(){
     return ssNnsynsub[2045];
    }
    double GetssNnsynsub2046(){
     return ssNnsynsub[2046];
    }
    double GetssNnsynsub2047(){
     return ssNnsynsub[2047];
    }
    double GetssNnsynsub2048(){
     return ssNnsynsub[2048];
    }
    double GetssNnsynsub2049(){
     return ssNnsynsub[2049];
    }
    double GetssNnsynsub2050(){
     return ssNnsynsub[2050];
    }
    double GetssNnsynsub2051(){
     return ssNnsynsub[2051];
    }
    double GetssNnsynsub2052(){
     return ssNnsynsub[2052];
    }
    double GetssNnsynsub2053(){
     return ssNnsynsub[2053];
    }
    double GetssNnsynsub2054(){
     return ssNnsynsub[2054];
    }
    double GetssNnsynsub2055(){
     return ssNnsynsub[2055];
    }
    double GetssNnsynsub2056(){
     return ssNnsynsub[2056];
    }
    double GetssNnsynsub2057(){
     return ssNnsynsub[2057];
    }
    double GetssNnsynsub2058(){
     return ssNnsynsub[2058];
    }
    double GetssNnsynsub2059(){
     return ssNnsynsub[2059];
    }
    double GetssNnsynsub2060(){
     return ssNnsynsub[2060];
    }
    double GetssNnsynsub2061(){
     return ssNnsynsub[2061];
    }
    double GetssNnsynsub2062(){
     return ssNnsynsub[2062];
    }
    double GetssNnsynsub2063(){
     return ssNnsynsub[2063];
    }
    double GetssNnsynsub2064(){
     return ssNnsynsub[2064];
    }
    double GetssNnsynsub2065(){
     return ssNnsynsub[2065];
    }
    double GetssNnsynsub2066(){
     return ssNnsynsub[2066];
    }
    double GetssNnsynsub2067(){
     return ssNnsynsub[2067];
    }
    double GetssNnsynsub2068(){
     return ssNnsynsub[2068];
    }
    double GetssNnsynsub2069(){
     return ssNnsynsub[2069];
    }
    double GetssNnsynsub2070(){
     return ssNnsynsub[2070];
    }
    double GetssNnsynsub2071(){
     return ssNnsynsub[2071];
    }
    double GetssNnsynsub2072(){
     return ssNnsynsub[2072];
    }
    double GetssNnsynsub2073(){
     return ssNnsynsub[2073];
    }
    double GetssNnsynsub2074(){
     return ssNnsynsub[2074];
    }
    double GetssNnsynsub2075(){
     return ssNnsynsub[2075];
    }
    double GetssNnsynsub2076(){
     return ssNnsynsub[2076];
    }
    double GetssNnsynsub2077(){
     return ssNnsynsub[2077];
    }
    double GetssNnsynsub2078(){
     return ssNnsynsub[2078];
    }
    double GetssNnsynsub2079(){
     return ssNnsynsub[2079];
    }
    double GetssNnsynsub2080(){
     return ssNnsynsub[2080];
    }
    double GetssNnsynsub2081(){
     return ssNnsynsub[2081];
    }
    double GetssNnsynsub2082(){
     return ssNnsynsub[2082];
    }
    double GetssNnsynsub2083(){
     return ssNnsynsub[2083];
    }
    double GetssNnsynsub2084(){
     return ssNnsynsub[2084];
    }
    double GetssNnsynsub2085(){
     return ssNnsynsub[2085];
    }
    double GetssNnsynsub2086(){
     return ssNnsynsub[2086];
    }
    double GetssNnsynsub2087(){
     return ssNnsynsub[2087];
    }
    double GetssNnsynsub2088(){
     return ssNnsynsub[2088];
    }
    double GetssNnsynsub2089(){
     return ssNnsynsub[2089];
    }
    double GetssNnsynsub2090(){
     return ssNnsynsub[2090];
    }
    double GetssNnsynsub2091(){
     return ssNnsynsub[2091];
    }
    double GetssNnsynsub2092(){
     return ssNnsynsub[2092];
    }
    double GetssNnsynsub2093(){
     return ssNnsynsub[2093];
    }
    double GetssNnsynsub2094(){
     return ssNnsynsub[2094];
    }
    double GetssNnsynsub2095(){
     return ssNnsynsub[2095];
    }
    double GetssNnsynsub2096(){
     return ssNnsynsub[2096];
    }
    double GetssNnsynsub2097(){
     return ssNnsynsub[2097];
    }
    double GetssNnsynsub2098(){
     return ssNnsynsub[2098];
    }
    double GetssNnsynsub2099(){
     return ssNnsynsub[2099];
    }
    double GetssNnsynsub2100(){
     return ssNnsynsub[2100];
    }
    double GetssNnsynsub2101(){
     return ssNnsynsub[2101];
    }
    double GetssNnsynsub2102(){
     return ssNnsynsub[2102];
    }
    double GetssNnsynsub2103(){
     return ssNnsynsub[2103];
    }
    double GetssNnsynsub2104(){
     return ssNnsynsub[2104];
    }
    double GetssNnsynsub2105(){
     return ssNnsynsub[2105];
    }
    double GetssNnsynsub2106(){
     return ssNnsynsub[2106];
    }
    double GetssNnsynsub2107(){
     return ssNnsynsub[2107];
    }
    double GetssNnsynsub2108(){
     return ssNnsynsub[2108];
    }
    double GetssNnsynsub2109(){
     return ssNnsynsub[2109];
    }
    double GetssNnsynsub2110(){
     return ssNnsynsub[2110];
    }
    double GetssNnsynsub2111(){
     return ssNnsynsub[2111];
    }
    double GetssNnsynsub2112(){
     return ssNnsynsub[2112];
    }
    double GetssNnsynsub2113(){
     return ssNnsynsub[2113];
    }
    double GetssNnsynsub2114(){
     return ssNnsynsub[2114];
    }
    double GetssNnsynsub2115(){
     return ssNnsynsub[2115];
    }
    double GetssNnsynsub2116(){
     return ssNnsynsub[2116];
    }
    double GetssNnsynsub2117(){
     return ssNnsynsub[2117];
    }
    double GetssNnsynsub2118(){
     return ssNnsynsub[2118];
    }
    double GetssNnsynsub2119(){
     return ssNnsynsub[2119];
    }
    double GetssNnsynsub2120(){
     return ssNnsynsub[2120];
    }
    double GetssNnsynsub2121(){
     return ssNnsynsub[2121];
    }
    double GetssNnsynsub2122(){
     return ssNnsynsub[2122];
    }
    double GetssNnsynsub2123(){
     return ssNnsynsub[2123];
    }
    double GetssNnsynsub2124(){
     return ssNnsynsub[2124];
    }
    double GetssNnsynsub2125(){
     return ssNnsynsub[2125];
    }
    double GetssNnsynsub2126(){
     return ssNnsynsub[2126];
    }
    double GetssNnsynsub2127(){
     return ssNnsynsub[2127];
    }
    double GetssNnsynsub2128(){
     return ssNnsynsub[2128];
    }
    double GetssNnsynsub2129(){
     return ssNnsynsub[2129];
    }
    double GetssNnsynsub2130(){
     return ssNnsynsub[2130];
    }
    double GetssNnsynsub2131(){
     return ssNnsynsub[2131];
    }
    double GetssNnsynsub2132(){
     return ssNnsynsub[2132];
    }
    double GetssNnsynsub2133(){
     return ssNnsynsub[2133];
    }
    double GetssNnsynsub2134(){
     return ssNnsynsub[2134];
    }
    double GetssNnsynsub2135(){
     return ssNnsynsub[2135];
    }
    double GetssNnsynsub2136(){
     return ssNnsynsub[2136];
    }
    double GetssNnsynsub2137(){
     return ssNnsynsub[2137];
    }
    double GetssNnsynsub2138(){
     return ssNnsynsub[2138];
    }
    double GetssNnsynsub2139(){
     return ssNnsynsub[2139];
    }
    double GetssNnsynsub2140(){
     return ssNnsynsub[2140];
    }
    double GetssNnsynsub2141(){
     return ssNnsynsub[2141];
    }
    double GetssNnsynsub2142(){
     return ssNnsynsub[2142];
    }
    double GetssNnsynsub2143(){
     return ssNnsynsub[2143];
    }
    double GetssNnsynsub2144(){
     return ssNnsynsub[2144];
    }
    double GetssNnsynsub2145(){
     return ssNnsynsub[2145];
    }
    double GetssNnsynsub2146(){
     return ssNnsynsub[2146];
    }
    double GetssNnsynsub2147(){
     return ssNnsynsub[2147];
    }
    double GetssNnsynsub2148(){
     return ssNnsynsub[2148];
    }
    double GetssNnsynsub2149(){
     return ssNnsynsub[2149];
    }
    double GetssNnsynsub2150(){
     return ssNnsynsub[2150];
    }
    double GetssNnsynsub2151(){
     return ssNnsynsub[2151];
    }
    double GetssNnsynsub2152(){
     return ssNnsynsub[2152];
    }
    double GetssNnsynsub2153(){
     return ssNnsynsub[2153];
    }
    double GetssNnsynsub2154(){
     return ssNnsynsub[2154];
    }
    double GetssNnsynsub2155(){
     return ssNnsynsub[2155];
    }
    double GetssNnsynsub2156(){
     return ssNnsynsub[2156];
    }
    double GetssNnsynsub2157(){
     return ssNnsynsub[2157];
    }
    double GetssNnsynsub2158(){
     return ssNnsynsub[2158];
    }
    double GetssNnsynsub2159(){
     return ssNnsynsub[2159];
    }
    double GetssNnsynsub2160(){
     return ssNnsynsub[2160];
    }
    double GetssNnsynsub2161(){
     return ssNnsynsub[2161];
    }
    double GetssNnsynsub2162(){
     return ssNnsynsub[2162];
    }
    double GetssNnsynsub2163(){
     return ssNnsynsub[2163];
    }
    double GetssNnsynsub2164(){
     return ssNnsynsub[2164];
    }
    double GetssNnsynsub2165(){
     return ssNnsynsub[2165];
    }
    double GetssNnsynsub2166(){
     return ssNnsynsub[2166];
    }
    double GetssNnsynsub2167(){
     return ssNnsynsub[2167];
    }
    double GetssNnsynsub2168(){
     return ssNnsynsub[2168];
    }
    double GetssNnsynsub2169(){
     return ssNnsynsub[2169];
    }
    double GetssNnsynsub2170(){
     return ssNnsynsub[2170];
    }
    double GetssNnsynsub2171(){
     return ssNnsynsub[2171];
    }
    double GetssNnsynsub2172(){
     return ssNnsynsub[2172];
    }
    double GetssNnsynsub2173(){
     return ssNnsynsub[2173];
    }
    double GetssNnsynsub2174(){
     return ssNnsynsub[2174];
    }
    double GetssNnsynsub2175(){
     return ssNnsynsub[2175];
    }
    double GetssNnsynsub2176(){
     return ssNnsynsub[2176];
    }
    double GetssNnsynsub2177(){
     return ssNnsynsub[2177];
    }
    double GetssNnsynsub2178(){
     return ssNnsynsub[2178];
    }
    double GetssNnsynsub2179(){
     return ssNnsynsub[2179];
    }
    double GetssNnsynsub2180(){
     return ssNnsynsub[2180];
    }
    double GetssNnsynsub2181(){
     return ssNnsynsub[2181];
    }
    double GetssNnsynsub2182(){
     return ssNnsynsub[2182];
    }
    double GetssNnsynsub2183(){
     return ssNnsynsub[2183];
    }
    double GetssNnsynsub2184(){
     return ssNnsynsub[2184];
    }
    double GetssNnsynsub2185(){
     return ssNnsynsub[2185];
    }
    double GetssNnsynsub2186(){
     return ssNnsynsub[2186];
    }
    double GetssNnsynsub2187(){
     return ssNnsynsub[2187];
    }
    double GetssNnsynsub2188(){
     return ssNnsynsub[2188];
    }
    double GetssNnsynsub2189(){
     return ssNnsynsub[2189];
    }
    double GetssNnsynsub2190(){
     return ssNnsynsub[2190];
    }
    double GetssNnsynsub2191(){
     return ssNnsynsub[2191];
    }
    double GetssNnsynsub2192(){
     return ssNnsynsub[2192];
    }
    double GetssNnsynsub2193(){
     return ssNnsynsub[2193];
    }
    double GetssNnsynsub2194(){
     return ssNnsynsub[2194];
    }
    double GetssNnsynsub2195(){
     return ssNnsynsub[2195];
    }
    double GetssNnsynsub2196(){
     return ssNnsynsub[2196];
    }
    double GetssNnsynsub2197(){
     return ssNnsynsub[2197];
    }
    double GetssNnsynsub2198(){
     return ssNnsynsub[2198];
    }
    double GetssNnsynsub2199(){
     return ssNnsynsub[2199];
    }
    double GetssNnsynsub2200(){
     return ssNnsynsub[2200];
    }
    double GetssNnsynsub2201(){
     return ssNnsynsub[2201];
    }
    double GetssNnsynsub2202(){
     return ssNnsynsub[2202];
    }
    double GetssNnsynsub2203(){
     return ssNnsynsub[2203];
    }
    double GetssNnsynsub2204(){
     return ssNnsynsub[2204];
    }
    double GetssNnsynsub2205(){
     return ssNnsynsub[2205];
    }
    double GetssNnsynsub2206(){
     return ssNnsynsub[2206];
    }
    double GetssNnsynsub2207(){
     return ssNnsynsub[2207];
    }
    double GetssNnsynsub2208(){
     return ssNnsynsub[2208];
    }
    double GetssNnsynsub2209(){
     return ssNnsynsub[2209];
    }
    double GetssNnsynsub2210(){
     return ssNnsynsub[2210];
    }
    double GetssNnsynsub2211(){
     return ssNnsynsub[2211];
    }
    double GetssNnsynsub2212(){
     return ssNnsynsub[2212];
    }
    double GetssNnsynsub2213(){
     return ssNnsynsub[2213];
    }
    double GetssNnsynsub2214(){
     return ssNnsynsub[2214];
    }
    double GetssNnsynsub2215(){
     return ssNnsynsub[2215];
    }
    double GetssNnsynsub2216(){
     return ssNnsynsub[2216];
    }
    double GetssNnsynsub2217(){
     return ssNnsynsub[2217];
    }
    double GetssNnsynsub2218(){
     return ssNnsynsub[2218];
    }
    double GetssNnsynsub2219(){
     return ssNnsynsub[2219];
    }
    double GetssNnsynsub2220(){
     return ssNnsynsub[2220];
    }
    double GetssNnsynsub2221(){
     return ssNnsynsub[2221];
    }
    double GetssNnsynsub2222(){
     return ssNnsynsub[2222];
    }
    double GetssNnsynsub2223(){
     return ssNnsynsub[2223];
    }
    double GetssNnsynsub2224(){
     return ssNnsynsub[2224];
    }
    double GetssNnsynsub2225(){
     return ssNnsynsub[2225];
    }
    double GetssNnsynsub2226(){
     return ssNnsynsub[2226];
    }
    double GetssNnsynsub2227(){
     return ssNnsynsub[2227];
    }
    double GetssNnsynsub2228(){
     return ssNnsynsub[2228];
    }
    double GetssNnsynsub2229(){
     return ssNnsynsub[2229];
    }
    double GetssNnsynsub2230(){
     return ssNnsynsub[2230];
    }
    double GetssNnsynsub2231(){
     return ssNnsynsub[2231];
    }
    double GetssNnsynsub2232(){
     return ssNnsynsub[2232];
    }
    double GetssNnsynsub2233(){
     return ssNnsynsub[2233];
    }
    double GetssNnsynsub2234(){
     return ssNnsynsub[2234];
    }
    double GetssNnsynsub2235(){
     return ssNnsynsub[2235];
    }
    double GetssNnsynsub2236(){
     return ssNnsynsub[2236];
    }
    double GetssNnsynsub2237(){
     return ssNnsynsub[2237];
    }
    double GetssNnsynsub2238(){
     return ssNnsynsub[2238];
    }
    double GetssNnsynsub2239(){
     return ssNnsynsub[2239];
    }
    double GetssNnsynsub2240(){
     return ssNnsynsub[2240];
    }
    double GetssNnsynsub2241(){
     return ssNnsynsub[2241];
    }
    double GetssNnsynsub2242(){
     return ssNnsynsub[2242];
    }
    double GetssNnsynsub2243(){
     return ssNnsynsub[2243];
    }
    double GetssNnsynsub2244(){
     return ssNnsynsub[2244];
    }
    double GetssNnsynsub2245(){
     return ssNnsynsub[2245];
    }
    double GetssNnsynsub2246(){
     return ssNnsynsub[2246];
    }
    double GetssNnsynsub2247(){
     return ssNnsynsub[2247];
    }
    double GetssNnsynsub2248(){
     return ssNnsynsub[2248];
    }
    double GetssNnsynsub2249(){
     return ssNnsynsub[2249];
    }
    double GetssNnsynsub2250(){
     return ssNnsynsub[2250];
    }
    double GetssNnsynsub2251(){
     return ssNnsynsub[2251];
    }
    double GetssNnsynsub2252(){
     return ssNnsynsub[2252];
    }
    double GetssNnsynsub2253(){
     return ssNnsynsub[2253];
    }
    double GetssNnsynsub2254(){
     return ssNnsynsub[2254];
    }
    double GetssNnsynsub2255(){
     return ssNnsynsub[2255];
    }
    double GetssNnsynsub2256(){
     return ssNnsynsub[2256];
    }
    double GetssNnsynsub2257(){
     return ssNnsynsub[2257];
    }
    double GetssNnsynsub2258(){
     return ssNnsynsub[2258];
    }
    double GetssNnsynsub2259(){
     return ssNnsynsub[2259];
    }
    double GetssNnsynsub2260(){
     return ssNnsynsub[2260];
    }
    double GetssNnsynsub2261(){
     return ssNnsynsub[2261];
    }
    double GetssNnsynsub2262(){
     return ssNnsynsub[2262];
    }
    double GetssNnsynsub2263(){
     return ssNnsynsub[2263];
    }
    double GetssNnsynsub2264(){
     return ssNnsynsub[2264];
    }
    double GetssNnsynsub2265(){
     return ssNnsynsub[2265];
    }
    double GetssNnsynsub2266(){
     return ssNnsynsub[2266];
    }
    double GetssNnsynsub2267(){
     return ssNnsynsub[2267];
    }
    double GetssNnsynsub2268(){
     return ssNnsynsub[2268];
    }
    double GetssNnsynsub2269(){
     return ssNnsynsub[2269];
    }
    double GetssNnsynsub2270(){
     return ssNnsynsub[2270];
    }
    double GetssNnsynsub2271(){
     return ssNnsynsub[2271];
    }
    double GetssNnsynsub2272(){
     return ssNnsynsub[2272];
    }
    double GetssNnsynsub2273(){
     return ssNnsynsub[2273];
    }
    double GetssNnsynsub2274(){
     return ssNnsynsub[2274];
    }
    double GetssNnsynsub2275(){
     return ssNnsynsub[2275];
    }
    double GetssNnsynsub2276(){
     return ssNnsynsub[2276];
    }
    double GetssNnsynsub2277(){
     return ssNnsynsub[2277];
    }
    double GetssNnsynsub2278(){
     return ssNnsynsub[2278];
    }
    double GetssNnsynsub2279(){
     return ssNnsynsub[2279];
    }
    double GetssNnsynsub2280(){
     return ssNnsynsub[2280];
    }
    double GetssNnsynsub2281(){
     return ssNnsynsub[2281];
    }
    double GetssNnsynsub2282(){
     return ssNnsynsub[2282];
    }
    double GetssNnsynsub2283(){
     return ssNnsynsub[2283];
    }
    double GetssNnsynsub2284(){
     return ssNnsynsub[2284];
    }
    double GetssNnsynsub2285(){
     return ssNnsynsub[2285];
    }
    double GetssNnsynsub2286(){
     return ssNnsynsub[2286];
    }
    double GetssNnsynsub2287(){
     return ssNnsynsub[2287];
    }
    double GetssNnsynsub2288(){
     return ssNnsynsub[2288];
    }
    double GetssNnsynsub2289(){
     return ssNnsynsub[2289];
    }
    double GetssNnsynsub2290(){
     return ssNnsynsub[2290];
    }
    double GetssNnsynsub2291(){
     return ssNnsynsub[2291];
    }
    double GetssNnsynsub2292(){
     return ssNnsynsub[2292];
    }
    double GetssNnsynsub2293(){
     return ssNnsynsub[2293];
    }
    double GetssNnsynsub2294(){
     return ssNnsynsub[2294];
    }
    double GetssNnsynsub2295(){
     return ssNnsynsub[2295];
    }
    double GetssNnsynsub2296(){
     return ssNnsynsub[2296];
    }
    double GetssNnsynsub2297(){
     return ssNnsynsub[2297];
    }
    double GetssNnsynsub2298(){
     return ssNnsynsub[2298];
    }
    double GetssNnsynsub2299(){
     return ssNnsynsub[2299];
    }
    double GetssNnsynsub2300(){
     return ssNnsynsub[2300];
    }
    double GetssNnsynsub2301(){
     return ssNnsynsub[2301];
    }
    double GetssNnsynsub2302(){
     return ssNnsynsub[2302];
    }
    double GetssNnsynsub2303(){
     return ssNnsynsub[2303];
    }
    double GetssNnsynsub2304(){
     return ssNnsynsub[2304];
    }
    double GetssNnsynsub2305(){
     return ssNnsynsub[2305];
    }
    double GetssNnsynsub2306(){
     return ssNnsynsub[2306];
    }
    double GetssNnsynsub2307(){
     return ssNnsynsub[2307];
    }
    double GetssNnsynsub2308(){
     return ssNnsynsub[2308];
    }
    double GetssNnsynsub2309(){
     return ssNnsynsub[2309];
    }
    double GetssNnsynsub2310(){
     return ssNnsynsub[2310];
    }
    double GetssNnsynsub2311(){
     return ssNnsynsub[2311];
    }
    double GetssNnsynsub2312(){
     return ssNnsynsub[2312];
    }
    double GetssNnsynsub2313(){
     return ssNnsynsub[2313];
    }
    double GetssNnsynsub2314(){
     return ssNnsynsub[2314];
    }
    double GetssNnsynsub2315(){
     return ssNnsynsub[2315];
    }
    double GetssNnsynsub2316(){
     return ssNnsynsub[2316];
    }
    double GetssNnsynsub2317(){
     return ssNnsynsub[2317];
    }
    double GetssNnsynsub2318(){
     return ssNnsynsub[2318];
    }
    double GetssNnsynsub2319(){
     return ssNnsynsub[2319];
    }
    double GetssNnsynsub2320(){
     return ssNnsynsub[2320];
    }
    double GetssNnsynsub2321(){
     return ssNnsynsub[2321];
    }
    double GetssNnsynsub2322(){
     return ssNnsynsub[2322];
    }
    double GetssNnsynsub2323(){
     return ssNnsynsub[2323];
    }
    double GetssNnsynsub2324(){
     return ssNnsynsub[2324];
    }
    double GetssNnsynsub2325(){
     return ssNnsynsub[2325];
    }
    double GetssNnsynsub2326(){
     return ssNnsynsub[2326];
    }
    double GetssNnsynsub2327(){
     return ssNnsynsub[2327];
    }
    double GetssNnsynsub2328(){
     return ssNnsynsub[2328];
    }
    double GetssNnsynsub2329(){
     return ssNnsynsub[2329];
    }
    double GetssNnsynsub2330(){
     return ssNnsynsub[2330];
    }
    double GetssNnsynsub2331(){
     return ssNnsynsub[2331];
    }
    double GetssNnsynsub2332(){
     return ssNnsynsub[2332];
    }
    double GetssNnsynsub2333(){
     return ssNnsynsub[2333];
    }
    double GetssNnsynsub2334(){
     return ssNnsynsub[2334];
    }
    double GetssNnsynsub2335(){
     return ssNnsynsub[2335];
    }
    double GetssNnsynsub2336(){
     return ssNnsynsub[2336];
    }
    double GetssNnsynsub2337(){
     return ssNnsynsub[2337];
    }
    double GetssNnsynsub2338(){
     return ssNnsynsub[2338];
    }
    double GetssNnsynsub2339(){
     return ssNnsynsub[2339];
    }
    double GetssNnsynsub2340(){
     return ssNnsynsub[2340];
    }
    double GetssNnsynsub2341(){
     return ssNnsynsub[2341];
    }
    double GetssNnsynsub2342(){
     return ssNnsynsub[2342];
    }
    double GetssNnsynsub2343(){
     return ssNnsynsub[2343];
    }
    double GetssNnsynsub2344(){
     return ssNnsynsub[2344];
    }
    double GetssNnsynsub2345(){
     return ssNnsynsub[2345];
    }
    double GetssNnsynsub2346(){
     return ssNnsynsub[2346];
    }
    double GetssNnsynsub2347(){
     return ssNnsynsub[2347];
    }
    double GetssNnsynsub2348(){
     return ssNnsynsub[2348];
    }
    double GetssNnsynsub2349(){
     return ssNnsynsub[2349];
    }
    double GetssNnsynsub2350(){
     return ssNnsynsub[2350];
    }
    double GetssNnsynsub2351(){
     return ssNnsynsub[2351];
    }
    double GetssNnsynsub2352(){
     return ssNnsynsub[2352];
    }
    double GetssNnsynsub2353(){
     return ssNnsynsub[2353];
    }
    double GetssNnsynsub2354(){
     return ssNnsynsub[2354];
    }
    double GetssNnsynsub2355(){
     return ssNnsynsub[2355];
    }
    double GetssNnsynsub2356(){
     return ssNnsynsub[2356];
    }
    double GetssNnsynsub2357(){
     return ssNnsynsub[2357];
    }
    double GetssNnsynsub2358(){
     return ssNnsynsub[2358];
    }
    double GetssNnsynsub2359(){
     return ssNnsynsub[2359];
    }
    double GetssNnsynsub2360(){
     return ssNnsynsub[2360];
    }
    double GetssNnsynsub2361(){
     return ssNnsynsub[2361];
    }
    double GetssNnsynsub2362(){
     return ssNnsynsub[2362];
    }
    double GetssNnsynsub2363(){
     return ssNnsynsub[2363];
    }
    double GetssNnsynsub2364(){
     return ssNnsynsub[2364];
    }
    double GetssNnsynsub2365(){
     return ssNnsynsub[2365];
    }
    double GetssNnsynsub2366(){
     return ssNnsynsub[2366];
    }
    double GetssNnsynsub2367(){
     return ssNnsynsub[2367];
    }
    double GetssNnsynsub2368(){
     return ssNnsynsub[2368];
    }
    double GetssNnsynsub2369(){
     return ssNnsynsub[2369];
    }
    double GetssNnsynsub2370(){
     return ssNnsynsub[2370];
    }
    double GetssNnsynsub2371(){
     return ssNnsynsub[2371];
    }
    double GetssNnsynsub2372(){
     return ssNnsynsub[2372];
    }
    double GetssNnsynsub2373(){
     return ssNnsynsub[2373];
    }
    double GetssNnsynsub2374(){
     return ssNnsynsub[2374];
    }
    double GetssNnsynsub2375(){
     return ssNnsynsub[2375];
    }
    double GetssNnsynsub2376(){
     return ssNnsynsub[2376];
    }
    double GetssNnsynsub2377(){
     return ssNnsynsub[2377];
    }
    double GetssNnsynsub2378(){
     return ssNnsynsub[2378];
    }
    double GetssNnsynsub2379(){
     return ssNnsynsub[2379];
    }
    double GetssNnsynsub2380(){
     return ssNnsynsub[2380];
    }
    double GetssNnsynsub2381(){
     return ssNnsynsub[2381];
    }
    double GetssNnsynsub2382(){
     return ssNnsynsub[2382];
    }
    double GetssNnsynsub2383(){
     return ssNnsynsub[2383];
    }
    double GetssNnsynsub2384(){
     return ssNnsynsub[2384];
    }
    double GetssNnsynsub2385(){
     return ssNnsynsub[2385];
    }
    double GetssNnsynsub2386(){
     return ssNnsynsub[2386];
    }
    double GetssNnsynsub2387(){
     return ssNnsynsub[2387];
    }
    double GetssNnsynsub2388(){
     return ssNnsynsub[2388];
    }
    double GetssNnsynsub2389(){
     return ssNnsynsub[2389];
    }
    double GetssNnsynsub2390(){
     return ssNnsynsub[2390];
    }
    double GetssNnsynsub2391(){
     return ssNnsynsub[2391];
    }
    double GetssNnsynsub2392(){
     return ssNnsynsub[2392];
    }
    double GetssNnsynsub2393(){
     return ssNnsynsub[2393];
    }
    double GetssNnsynsub2394(){
     return ssNnsynsub[2394];
    }
    double GetssNnsynsub2395(){
     return ssNnsynsub[2395];
    }
    double GetssNnsynsub2396(){
     return ssNnsynsub[2396];
    }
    double GetssNnsynsub2397(){
     return ssNnsynsub[2397];
    }
    double GetssNnsynsub2398(){
     return ssNnsynsub[2398];
    }
    double GetssNnsynsub2399(){
     return ssNnsynsub[2399];
    }
    double GetssNnsynsub2400(){
     return ssNnsynsub[2400];
    }
    double GetssNnsynsub2401(){
     return ssNnsynsub[2401];
    }
    double GetssNnsynsub2402(){
     return ssNnsynsub[2402];
    }
    double GetssNnsynsub2403(){
     return ssNnsynsub[2403];
    }
    double GetssNnsynsub2404(){
     return ssNnsynsub[2404];
    }
    double GetssNnsynsub2405(){
     return ssNnsynsub[2405];
    }
    double GetssNnsynsub2406(){
     return ssNnsynsub[2406];
    }
    double GetssNnsynsub2407(){
     return ssNnsynsub[2407];
    }
    double GetssNnsynsub2408(){
     return ssNnsynsub[2408];
    }
    double GetssNnsynsub2409(){
     return ssNnsynsub[2409];
    }
    double GetssNnsynsub2410(){
     return ssNnsynsub[2410];
    }
    double GetssNnsynsub2411(){
     return ssNnsynsub[2411];
    }
    double GetssNnsynsub2412(){
     return ssNnsynsub[2412];
    }
    double GetssNnsynsub2413(){
     return ssNnsynsub[2413];
    }
    double GetssNnsynsub2414(){
     return ssNnsynsub[2414];
    }
    double GetssNnsynsub2415(){
     return ssNnsynsub[2415];
    }
    double GetssNnsynsub2416(){
     return ssNnsynsub[2416];
    }
    double GetssNnsynsub2417(){
     return ssNnsynsub[2417];
    }
    double GetssNnsynsub2418(){
     return ssNnsynsub[2418];
    }
    double GetssNnsynsub2419(){
     return ssNnsynsub[2419];
    }
    double GetssNnsynsub2420(){
     return ssNnsynsub[2420];
    }
    double GetssNnsynsub2421(){
     return ssNnsynsub[2421];
    }
    double GetssNnsynsub2422(){
     return ssNnsynsub[2422];
    }
    double GetssNnsynsub2423(){
     return ssNnsynsub[2423];
    }
    double GetssNnsynsub2424(){
     return ssNnsynsub[2424];
    }
    double GetssNnsynsub2425(){
     return ssNnsynsub[2425];
    }
    double GetssNnsynsub2426(){
     return ssNnsynsub[2426];
    }
    double GetssNnsynsub2427(){
     return ssNnsynsub[2427];
    }
    double GetssNnsynsub2428(){
     return ssNnsynsub[2428];
    }
    double GetssNnsynsub2429(){
     return ssNnsynsub[2429];
    }
    double GetssNnsynsub2430(){
     return ssNnsynsub[2430];
    }
    double GetssNnsynsub2431(){
     return ssNnsynsub[2431];
    }
    double GetssNnsynsub2432(){
     return ssNnsynsub[2432];
    }
    double GetssNnsynsub2433(){
     return ssNnsynsub[2433];
    }
    double GetssNnsynsub2434(){
     return ssNnsynsub[2434];
    }
    double GetssNnsynsub2435(){
     return ssNnsynsub[2435];
    }
    double GetssNnsynsub2436(){
     return ssNnsynsub[2436];
    }
    double GetssNnsynsub2437(){
     return ssNnsynsub[2437];
    }
    double GetssNnsynsub2438(){
     return ssNnsynsub[2438];
    }
    double GetssNnsynsub2439(){
     return ssNnsynsub[2439];
    }
    double GetssNnsynsub2440(){
     return ssNnsynsub[2440];
    }
    double GetssNnsynsub2441(){
     return ssNnsynsub[2441];
    }
    double GetssNnsynsub2442(){
     return ssNnsynsub[2442];
    }
    double GetssNnsynsub2443(){
     return ssNnsynsub[2443];
    }
    double GetssNnsynsub2444(){
     return ssNnsynsub[2444];
    }
    double GetssNnsynsub2445(){
     return ssNnsynsub[2445];
    }
    double GetssNnsynsub2446(){
     return ssNnsynsub[2446];
    }
    double GetssNnsynsub2447(){
     return ssNnsynsub[2447];
    }
    double GetssNnsynsub2448(){
     return ssNnsynsub[2448];
    }
    double GetssNnsynsub2449(){
     return ssNnsynsub[2449];
    }
    double GetssNnsynsub2450(){
     return ssNnsynsub[2450];
    }
    double GetssNnsynsub2451(){
     return ssNnsynsub[2451];
    }
    double GetssNnsynsub2452(){
     return ssNnsynsub[2452];
    }
    double GetssNnsynsub2453(){
     return ssNnsynsub[2453];
    }
    double GetssNnsynsub2454(){
     return ssNnsynsub[2454];
    }
    double GetssNnsynsub2455(){
     return ssNnsynsub[2455];
    }
    double GetssNnsynsub2456(){
     return ssNnsynsub[2456];
    }
    double GetssNnsynsub2457(){
     return ssNnsynsub[2457];
    }
    double GetssNnsynsub2458(){
     return ssNnsynsub[2458];
    }
    double GetssNnsynsub2459(){
     return ssNnsynsub[2459];
    }
    double GetssNnsynsub2460(){
     return ssNnsynsub[2460];
    }
    double GetssNnsynsub2461(){
     return ssNnsynsub[2461];
    }
    double GetssNnsynsub2462(){
     return ssNnsynsub[2462];
    }
    double GetssNnsynsub2463(){
     return ssNnsynsub[2463];
    }
    double GetssNnsynsub2464(){
     return ssNnsynsub[2464];
    }
    double GetssNnsynsub2465(){
     return ssNnsynsub[2465];
    }
    double GetssNnsynsub2466(){
     return ssNnsynsub[2466];
    }
    double GetssNnsynsub2467(){
     return ssNnsynsub[2467];
    }
    double GetssNnsynsub2468(){
     return ssNnsynsub[2468];
    }
    double GetssNnsynsub2469(){
     return ssNnsynsub[2469];
    }
    double GetssNnsynsub2470(){
     return ssNnsynsub[2470];
    }
    double GetssNnsynsub2471(){
     return ssNnsynsub[2471];
    }
    double GetssNnsynsub2472(){
     return ssNnsynsub[2472];
    }
    double GetssNnsynsub2473(){
     return ssNnsynsub[2473];
    }
    double GetssNnsynsub2474(){
     return ssNnsynsub[2474];
    }
    double GetssNnsynsub2475(){
     return ssNnsynsub[2475];
    }
    double GetssNnsynsub2476(){
     return ssNnsynsub[2476];
    }
    double GetssNnsynsub2477(){
     return ssNnsynsub[2477];
    }
    double GetssNnsynsub2478(){
     return ssNnsynsub[2478];
    }
    double GetssNnsynsub2479(){
     return ssNnsynsub[2479];
    }
    double GetssNnsynsub2480(){
     return ssNnsynsub[2480];
    }
    double GetssNnsynsub2481(){
     return ssNnsynsub[2481];
    }
    double GetssNnsynsub2482(){
     return ssNnsynsub[2482];
    }
    double GetssNnsynsub2483(){
     return ssNnsynsub[2483];
    }
    double GetssNnsynsub2484(){
     return ssNnsynsub[2484];
    }
    double GetssNnsynsub2485(){
     return ssNnsynsub[2485];
    }
    double GetssNnsynsub2486(){
     return ssNnsynsub[2486];
    }
    double GetssNnsynsub2487(){
     return ssNnsynsub[2487];
    }
    double GetssNnsynsub2488(){
     return ssNnsynsub[2488];
    }
    double GetssNnsynsub2489(){
     return ssNnsynsub[2489];
    }
    double GetssNnsynsub2490(){
     return ssNnsynsub[2490];
    }
    double GetssNnsynsub2491(){
     return ssNnsynsub[2491];
    }
    double GetssNnsynsub2492(){
     return ssNnsynsub[2492];
    }
    double GetssNnsynsub2493(){
     return ssNnsynsub[2493];
    }
    double GetssNnsynsub2494(){
     return ssNnsynsub[2494];
    }
    double GetssNnsynsub2495(){
     return ssNnsynsub[2495];
    }
    double GetssNnsynsub2496(){
     return ssNnsynsub[2496];
    }
    double GetssNnsynsub2497(){
     return ssNnsynsub[2497];
    }
    double GetssNnsynsub2498(){
     return ssNnsynsub[2498];
    }
    double GetssNnsynsub2499(){
     return ssNnsynsub[2499];
    }
    double GetssNnsynsub2500(){
     return ssNnsynsub[2500];
    }
    double GetssNnsynsub2501(){
     return ssNnsynsub[2501];
    }
    double GetssNnsynsub2502(){
     return ssNnsynsub[2502];
    }
    double GetssNnsynsub2503(){
     return ssNnsynsub[2503];
    }
    double GetssNnsynsub2504(){
     return ssNnsynsub[2504];
    }
    double GetssNnsynsub2505(){
     return ssNnsynsub[2505];
    }
    double GetssNnsynsub2506(){
     return ssNnsynsub[2506];
    }
    double GetssNnsynsub2507(){
     return ssNnsynsub[2507];
    }
    double GetssNnsynsub2508(){
     return ssNnsynsub[2508];
    }
    double GetssNnsynsub2509(){
     return ssNnsynsub[2509];
    }
    double GetssNnsynsub2510(){
     return ssNnsynsub[2510];
    }
    double GetssNnsynsub2511(){
     return ssNnsynsub[2511];
    }
    double GetssNnsynsub2512(){
     return ssNnsynsub[2512];
    }
    double GetssNnsynsub2513(){
     return ssNnsynsub[2513];
    }
    double GetssNnsynsub2514(){
     return ssNnsynsub[2514];
    }
    double GetssNnsynsub2515(){
     return ssNnsynsub[2515];
    }
    double GetssNnsynsub2516(){
     return ssNnsynsub[2516];
    }
    double GetssNnsynsub2517(){
     return ssNnsynsub[2517];
    }
    double GetssNnsynsub2518(){
     return ssNnsynsub[2518];
    }
    double GetssNnsynsub2519(){
     return ssNnsynsub[2519];
    }
    double GetssNnsynsub2520(){
     return ssNnsynsub[2520];
    }
    double GetssNnsynsub2521(){
     return ssNnsynsub[2521];
    }
    double GetssNnsynsub2522(){
     return ssNnsynsub[2522];
    }
    double GetssNnsynsub2523(){
     return ssNnsynsub[2523];
    }
    double GetssNnsynsub2524(){
     return ssNnsynsub[2524];
    }
    double GetssNnsynsub2525(){
     return ssNnsynsub[2525];
    }
    double GetssNnsynsub2526(){
     return ssNnsynsub[2526];
    }
    double GetssNnsynsub2527(){
     return ssNnsynsub[2527];
    }
    double GetssNnsynsub2528(){
     return ssNnsynsub[2528];
    }
    double GetssNnsynsub2529(){
     return ssNnsynsub[2529];
    }
    double GetssNnsynsub2530(){
     return ssNnsynsub[2530];
    }
    double GetssNnsynsub2531(){
     return ssNnsynsub[2531];
    }
    double GetssNnsynsub2532(){
     return ssNnsynsub[2532];
    }
    double GetssNnsynsub2533(){
     return ssNnsynsub[2533];
    }
    double GetssNnsynsub2534(){
     return ssNnsynsub[2534];
    }
    double GetssNnsynsub2535(){
     return ssNnsynsub[2535];
    }
    double GetssNnsynsub2536(){
     return ssNnsynsub[2536];
    }
    double GetssNnsynsub2537(){
     return ssNnsynsub[2537];
    }
    double GetssNnsynsub2538(){
     return ssNnsynsub[2538];
    }
    double GetssNnsynsub2539(){
     return ssNnsynsub[2539];
    }
    double GetssNnsynsub2540(){
     return ssNnsynsub[2540];
    }
    double GetssNnsynsub2541(){
     return ssNnsynsub[2541];
    }
    double GetssNnsynsub2542(){
     return ssNnsynsub[2542];
    }
    double GetssNnsynsub2543(){
     return ssNnsynsub[2543];
    }
    double GetssNnsynsub2544(){
     return ssNnsynsub[2544];
    }
    double GetssNnsynsub2545(){
     return ssNnsynsub[2545];
    }
    double GetssNnsynsub2546(){
     return ssNnsynsub[2546];
    }
    double GetssNnsynsub2547(){
     return ssNnsynsub[2547];
    }
    double GetssNnsynsub2548(){
     return ssNnsynsub[2548];
    }
    double GetssNnsynsub2549(){
     return ssNnsynsub[2549];
    }
    double GetssNnsynsub2550(){
     return ssNnsynsub[2550];
    }
    double GetssNnsynsub2551(){
     return ssNnsynsub[2551];
    }
    double GetssNnsynsub2552(){
     return ssNnsynsub[2552];
    }
    double GetssNnsynsub2553(){
     return ssNnsynsub[2553];
    }
    double GetssNnsynsub2554(){
     return ssNnsynsub[2554];
    }
    double GetssNnsynsub2555(){
     return ssNnsynsub[2555];
    }
    double GetssNnsynsub2556(){
     return ssNnsynsub[2556];
    }
    double GetssNnsynsub2557(){
     return ssNnsynsub[2557];
    }
    double GetssNnsynsub2558(){
     return ssNnsynsub[2558];
    }
    double GetssNnsynsub2559(){
     return ssNnsynsub[2559];
    }
    double GetssNnsynsub2560(){
     return ssNnsynsub[2560];
    }
    double GetssNnsynsub2561(){
     return ssNnsynsub[2561];
    }
    double GetssNnsynsub2562(){
     return ssNnsynsub[2562];
    }
    double GetssNnsynsub2563(){
     return ssNnsynsub[2563];
    }
    double GetssNnsynsub2564(){
     return ssNnsynsub[2564];
    }
    double GetssNnsynsub2565(){
     return ssNnsynsub[2565];
    }
    double GetssNnsynsub2566(){
     return ssNnsynsub[2566];
    }
    double GetssNnsynsub2567(){
     return ssNnsynsub[2567];
    }
    double GetssNnsynsub2568(){
     return ssNnsynsub[2568];
    }
    double GetssNnsynsub2569(){
     return ssNnsynsub[2569];
    }
    double GetssNnsynsub2570(){
     return ssNnsynsub[2570];
    }
    double GetssNnsynsub2571(){
     return ssNnsynsub[2571];
    }
    double GetssNnsynsub2572(){
     return ssNnsynsub[2572];
    }
    double GetssNnsynsub2573(){
     return ssNnsynsub[2573];
    }
    double GetssNnsynsub2574(){
     return ssNnsynsub[2574];
    }
    double GetssNnsynsub2575(){
     return ssNnsynsub[2575];
    }
    double GetssNnsynsub2576(){
     return ssNnsynsub[2576];
    }
    double GetssNnsynsub2577(){
     return ssNnsynsub[2577];
    }
    double GetssNnsynsub2578(){
     return ssNnsynsub[2578];
    }
    double GetssNnsynsub2579(){
     return ssNnsynsub[2579];
    }
    double GetssNnsynsub2580(){
     return ssNnsynsub[2580];
    }
    double GetssNnsynsub2581(){
     return ssNnsynsub[2581];
    }
    double GetssNnsynsub2582(){
     return ssNnsynsub[2582];
    }
    double GetssNnsynsub2583(){
     return ssNnsynsub[2583];
    }
    double GetssNnsynsub2584(){
     return ssNnsynsub[2584];
    }
    double GetssNnsynsub2585(){
     return ssNnsynsub[2585];
    }
    double GetssNnsynsub2586(){
     return ssNnsynsub[2586];
    }
    double GetssNnsynsub2587(){
     return ssNnsynsub[2587];
    }
    double GetssNnsynsub2588(){
     return ssNnsynsub[2588];
    }
    double GetssNnsynsub2589(){
     return ssNnsynsub[2589];
    }
    double GetssNnsynsub2590(){
     return ssNnsynsub[2590];
    }
    double GetssNnsynsub2591(){
     return ssNnsynsub[2591];
    }
    double GetssNnsynsub2592(){
     return ssNnsynsub[2592];
    }
    double GetssNnsynsub2593(){
     return ssNnsynsub[2593];
    }
    double GetssNnsynsub2594(){
     return ssNnsynsub[2594];
    }
    double GetssNnsynsub2595(){
     return ssNnsynsub[2595];
    }
    double GetssNnsynsub2596(){
     return ssNnsynsub[2596];
    }
    double GetssNnsynsub2597(){
     return ssNnsynsub[2597];
    }
    double GetssNnsynsub2598(){
     return ssNnsynsub[2598];
    }
    double GetssNnsynsub2599(){
     return ssNnsynsub[2599];
    }
    double GetssNnsynsub2600(){
     return ssNnsynsub[2600];
    }
    double GetssNnsynsub2601(){
     return ssNnsynsub[2601];
    }
    double GetssNnsynsub2602(){
     return ssNnsynsub[2602];
    }
    double GetssNnsynsub2603(){
     return ssNnsynsub[2603];
    }
    double GetssNnsynsub2604(){
     return ssNnsynsub[2604];
    }
    double GetssNnsynsub2605(){
     return ssNnsynsub[2605];
    }
    double GetssNnsynsub2606(){
     return ssNnsynsub[2606];
    }
    double GetssNnsynsub2607(){
     return ssNnsynsub[2607];
    }
    double GetssNnsynsub2608(){
     return ssNnsynsub[2608];
    }
    double GetssNnsynsub2609(){
     return ssNnsynsub[2609];
    }
    double GetssNnsynsub2610(){
     return ssNnsynsub[2610];
    }
    double GetssNnsynsub2611(){
     return ssNnsynsub[2611];
    }
    double GetssNnsynsub2612(){
     return ssNnsynsub[2612];
    }
    double GetssNnsynsub2613(){
     return ssNnsynsub[2613];
    }
    double GetssNnsynsub2614(){
     return ssNnsynsub[2614];
    }
    double GetssNnsynsub2615(){
     return ssNnsynsub[2615];
    }
    double GetssNnsynsub2616(){
     return ssNnsynsub[2616];
    }
    double GetssNnsynsub2617(){
     return ssNnsynsub[2617];
    }
    double GetssNnsynsub2618(){
     return ssNnsynsub[2618];
    }
    double GetssNnsynsub2619(){
     return ssNnsynsub[2619];
    }
    double GetssNnsynsub2620(){
     return ssNnsynsub[2620];
    }
    double GetssNnsynsub2621(){
     return ssNnsynsub[2621];
    }
    double GetssNnsynsub2622(){
     return ssNnsynsub[2622];
    }
    double GetssNnsynsub2623(){
     return ssNnsynsub[2623];
    }
    double GetssNnsynsub2624(){
     return ssNnsynsub[2624];
    }
    double GetssNnsynsub2625(){
     return ssNnsynsub[2625];
    }
    double GetssNnsynsub2626(){
     return ssNnsynsub[2626];
    }
    double GetssNnsynsub2627(){
     return ssNnsynsub[2627];
    }
    double GetssNnsynsub2628(){
     return ssNnsynsub[2628];
    }
    double GetssNnsynsub2629(){
     return ssNnsynsub[2629];
    }
    double GetssNnsynsub2630(){
     return ssNnsynsub[2630];
    }
    double GetssNnsynsub2631(){
     return ssNnsynsub[2631];
    }
    double GetssNnsynsub2632(){
     return ssNnsynsub[2632];
    }
    double GetssNnsynsub2633(){
     return ssNnsynsub[2633];
    }
    double GetssNnsynsub2634(){
     return ssNnsynsub[2634];
    }
    double GetssNnsynsub2635(){
     return ssNnsynsub[2635];
    }
    double GetssNnsynsub2636(){
     return ssNnsynsub[2636];
    }
    double GetssNnsynsub2637(){
     return ssNnsynsub[2637];
    }
    double GetssNnsynsub2638(){
     return ssNnsynsub[2638];
    }
    double GetssNnsynsub2639(){
     return ssNnsynsub[2639];
    }
    double GetssNnsynsub2640(){
     return ssNnsynsub[2640];
    }
    double GetssNnsynsub2641(){
     return ssNnsynsub[2641];
    }
    double GetssNnsynsub2642(){
     return ssNnsynsub[2642];
    }
    double GetssNnsynsub2643(){
     return ssNnsynsub[2643];
    }
    double GetssNnsynsub2644(){
     return ssNnsynsub[2644];
    }
    double GetssNnsynsub2645(){
     return ssNnsynsub[2645];
    }
    double GetssNnsynsub2646(){
     return ssNnsynsub[2646];
    }
    double GetssNnsynsub2647(){
     return ssNnsynsub[2647];
    }
    double GetssNnsynsub2648(){
     return ssNnsynsub[2648];
    }
    double GetssNnsynsub2649(){
     return ssNnsynsub[2649];
    }
    double GetssNnsynsub2650(){
     return ssNnsynsub[2650];
    }
    double GetssNnsynsub2651(){
     return ssNnsynsub[2651];
    }
    double GetssNnsynsub2652(){
     return ssNnsynsub[2652];
    }
    double GetssNnsynsub2653(){
     return ssNnsynsub[2653];
    }
    double GetssNnsynsub2654(){
     return ssNnsynsub[2654];
    }
    double GetssNnsynsub2655(){
     return ssNnsynsub[2655];
    }
    double GetssNnsynsub2656(){
     return ssNnsynsub[2656];
    }
    double GetssNnsynsub2657(){
     return ssNnsynsub[2657];
    }
    double GetssNnsynsub2658(){
     return ssNnsynsub[2658];
    }
    double GetssNnsynsub2659(){
     return ssNnsynsub[2659];
    }
    double GetssNnsynsub2660(){
     return ssNnsynsub[2660];
    }
    double GetssNnsynsub2661(){
     return ssNnsynsub[2661];
    }
    double GetssNnsynsub2662(){
     return ssNnsynsub[2662];
    }
    double GetssNnsynsub2663(){
     return ssNnsynsub[2663];
    }
    double GetssNnsynsub2664(){
     return ssNnsynsub[2664];
    }
    double GetssNnsynsub2665(){
     return ssNnsynsub[2665];
    }
    double GetssNnsynsub2666(){
     return ssNnsynsub[2666];
    }
    double GetssNnsynsub2667(){
     return ssNnsynsub[2667];
    }
    double GetssNnsynsub2668(){
     return ssNnsynsub[2668];
    }
    double GetssNnsynsub2669(){
     return ssNnsynsub[2669];
    }
    double GetssNnsynsub2670(){
     return ssNnsynsub[2670];
    }
    double GetssNnsynsub2671(){
     return ssNnsynsub[2671];
    }
    double GetssNnsynsub2672(){
     return ssNnsynsub[2672];
    }
    double GetssNnsynsub2673(){
     return ssNnsynsub[2673];
    }
    double GetssNnsynsub2674(){
     return ssNnsynsub[2674];
    }
    double GetssNnsynsub2675(){
     return ssNnsynsub[2675];
    }
    double GetssNnsynsub2676(){
     return ssNnsynsub[2676];
    }
    double GetssNnsynsub2677(){
     return ssNnsynsub[2677];
    }
    double GetssNnsynsub2678(){
     return ssNnsynsub[2678];
    }
    double GetssNnsynsub2679(){
     return ssNnsynsub[2679];
    }
    double GetssNnsynsub2680(){
     return ssNnsynsub[2680];
    }
    double GetssNnsynsub2681(){
     return ssNnsynsub[2681];
    }
    double GetssNnsynsub2682(){
     return ssNnsynsub[2682];
    }
    double GetssNnsynsub2683(){
     return ssNnsynsub[2683];
    }
    double GetssNnsynsub2684(){
     return ssNnsynsub[2684];
    }
    double GetssNnsynsub2685(){
     return ssNnsynsub[2685];
    }
    double GetssNnsynsub2686(){
     return ssNnsynsub[2686];
    }
    double GetssNnsynsub2687(){
     return ssNnsynsub[2687];
    }
    double GetssNnsynsub2688(){
     return ssNnsynsub[2688];
    }
    double GetssNnsynsub2689(){
     return ssNnsynsub[2689];
    }
    double GetssNnsynsub2690(){
     return ssNnsynsub[2690];
    }
    double GetssNnsynsub2691(){
     return ssNnsynsub[2691];
    }
    double GetssNnsynsub2692(){
     return ssNnsynsub[2692];
    }
    double GetssNnsynsub2693(){
     return ssNnsynsub[2693];
    }
    double GetssNnsynsub2694(){
     return ssNnsynsub[2694];
    }
    double GetssNnsynsub2695(){
     return ssNnsynsub[2695];
    }
    double GetssNnsynsub2696(){
     return ssNnsynsub[2696];
    }
    double GetssNnsynsub2697(){
     return ssNnsynsub[2697];
    }
    double GetssNnsynsub2698(){
     return ssNnsynsub[2698];
    }
    double GetssNnsynsub2699(){
     return ssNnsynsub[2699];
    }
    double GetssNnsynsub2700(){
     return ssNnsynsub[2700];
    }
    double GetssNnsynsub2701(){
     return ssNnsynsub[2701];
    }
    double GetssNnsynsub2702(){
     return ssNnsynsub[2702];
    }
    double GetssNnsynsub2703(){
     return ssNnsynsub[2703];
    }
    double GetssNnsynsub2704(){
     return ssNnsynsub[2704];
    }
    double GetssNnsynsub2705(){
     return ssNnsynsub[2705];
    }
    double GetssNnsynsub2706(){
     return ssNnsynsub[2706];
    }
    double GetssNnsynsub2707(){
     return ssNnsynsub[2707];
    }
    double GetssNnsynsub2708(){
     return ssNnsynsub[2708];
    }
    double GetssNnsynsub2709(){
     return ssNnsynsub[2709];
    }
    double GetssNnsynsub2710(){
     return ssNnsynsub[2710];
    }
    double GetssNnsynsub2711(){
     return ssNnsynsub[2711];
    }
    double GetssNnsynsub2712(){
     return ssNnsynsub[2712];
    }
    double GetssNnsynsub2713(){
     return ssNnsynsub[2713];
    }
    double GetssNnsynsub2714(){
     return ssNnsynsub[2714];
    }
    double GetssNnsynsub2715(){
     return ssNnsynsub[2715];
    }
    double GetssNnsynsub2716(){
     return ssNnsynsub[2716];
    }
    double GetssNnsynsub2717(){
     return ssNnsynsub[2717];
    }
    double GetssNnsynsub2718(){
     return ssNnsynsub[2718];
    }
    double GetssNnsynsub2719(){
     return ssNnsynsub[2719];
    }
    double GetssNnsynsub2720(){
     return ssNnsynsub[2720];
    }
    double GetssNnsynsub2721(){
     return ssNnsynsub[2721];
    }
    double GetssNnsynsub2722(){
     return ssNnsynsub[2722];
    }
    double GetssNnsynsub2723(){
     return ssNnsynsub[2723];
    }
    double GetssNnsynsub2724(){
     return ssNnsynsub[2724];
    }
    double GetssNnsynsub2725(){
     return ssNnsynsub[2725];
    }
    double GetssNnsynsub2726(){
     return ssNnsynsub[2726];
    }
    double GetssNnsynsub2727(){
     return ssNnsynsub[2727];
    }
    double GetssNnsynsub2728(){
     return ssNnsynsub[2728];
    }
    double GetssNnsynsub2729(){
     return ssNnsynsub[2729];
    }
    double GetssNnsynsub2730(){
     return ssNnsynsub[2730];
    }
    double GetssNnsynsub2731(){
     return ssNnsynsub[2731];
    }
    double GetssNnsynsub2732(){
     return ssNnsynsub[2732];
    }
    double GetssNnsynsub2733(){
     return ssNnsynsub[2733];
    }
    double GetssNnsynsub2734(){
     return ssNnsynsub[2734];
    }
    double GetssNnsynsub2735(){
     return ssNnsynsub[2735];
    }
    double GetssNnsynsub2736(){
     return ssNnsynsub[2736];
    }
    double GetssNnsynsub2737(){
     return ssNnsynsub[2737];
    }
    double GetssNnsynsub2738(){
     return ssNnsynsub[2738];
    }
    double GetssNnsynsub2739(){
     return ssNnsynsub[2739];
    }
    double GetssNnsynsub2740(){
     return ssNnsynsub[2740];
    }
    double GetssNnsynsub2741(){
     return ssNnsynsub[2741];
    }
    double GetssNnsynsub2742(){
     return ssNnsynsub[2742];
    }
    double GetssNnsynsub2743(){
     return ssNnsynsub[2743];
    }
    double GetssNnsynsub2744(){
     return ssNnsynsub[2744];
    }
    double GetssNnsynsub2745(){
     return ssNnsynsub[2745];
    }
    double GetssNnsynsub2746(){
     return ssNnsynsub[2746];
    }
    double GetssNnsynsub2747(){
     return ssNnsynsub[2747];
    }
    double GetssNnsynsub2748(){
     return ssNnsynsub[2748];
    }
    double GetssNnsynsub2749(){
     return ssNnsynsub[2749];
    }
    double GetssNnsynsub2750(){
     return ssNnsynsub[2750];
    }
    double GetssNnsynsub2751(){
     return ssNnsynsub[2751];
    }
    double GetssNnsynsub2752(){
     return ssNnsynsub[2752];
    }
    double GetssNnsynsub2753(){
     return ssNnsynsub[2753];
    }
    double GetssNnsynsub2754(){
     return ssNnsynsub[2754];
    }
    double GetssNnsynsub2755(){
     return ssNnsynsub[2755];
    }
    double GetssNnsynsub2756(){
     return ssNnsynsub[2756];
    }
    double GetssNnsynsub2757(){
     return ssNnsynsub[2757];
    }
    double GetssNnsynsub2758(){
     return ssNnsynsub[2758];
    }
    double GetssNnsynsub2759(){
     return ssNnsynsub[2759];
    }
    double GetssNnsynsub2760(){
     return ssNnsynsub[2760];
    }
    double GetssNnsynsub2761(){
     return ssNnsynsub[2761];
    }
    double GetssNnsynsub2762(){
     return ssNnsynsub[2762];
    }
    double GetssNnsynsub2763(){
     return ssNnsynsub[2763];
    }
    double GetssNnsynsub2764(){
     return ssNnsynsub[2764];
    }
    double GetssNnsynsub2765(){
     return ssNnsynsub[2765];
    }
    double GetssNnsynsub2766(){
     return ssNnsynsub[2766];
    }
    double GetssNnsynsub2767(){
     return ssNnsynsub[2767];
    }
    double GetssNnsynsub2768(){
     return ssNnsynsub[2768];
    }
    double GetssNnsynsub2769(){
     return ssNnsynsub[2769];
    }
    double GetssNnsynsub2770(){
     return ssNnsynsub[2770];
    }
    double GetssNnsynsub2771(){
     return ssNnsynsub[2771];
    }
    double GetssNnsynsub2772(){
     return ssNnsynsub[2772];
    }
    double GetssNnsynsub2773(){
     return ssNnsynsub[2773];
    }
    double GetssNnsynsub2774(){
     return ssNnsynsub[2774];
    }
    double GetssNnsynsub2775(){
     return ssNnsynsub[2775];
    }
    double GetssNnsynsub2776(){
     return ssNnsynsub[2776];
    }
    double GetssNnsynsub2777(){
     return ssNnsynsub[2777];
    }
    double GetssNnsynsub2778(){
     return ssNnsynsub[2778];
    }
    double GetssNnsynsub2779(){
     return ssNnsynsub[2779];
    }
    double GetssNnsynsub2780(){
     return ssNnsynsub[2780];
    }
    double GetssNnsynsub2781(){
     return ssNnsynsub[2781];
    }
    double GetssNnsynsub2782(){
     return ssNnsynsub[2782];
    }
    double GetssNnsynsub2783(){
     return ssNnsynsub[2783];
    }
    double GetssNnsynsub2784(){
     return ssNnsynsub[2784];
    }
    double GetssNnsynsub2785(){
     return ssNnsynsub[2785];
    }
    double GetssNnsynsub2786(){
     return ssNnsynsub[2786];
    }
    double GetssNnsynsub2787(){
     return ssNnsynsub[2787];
    }
    double GetssNnsynsub2788(){
     return ssNnsynsub[2788];
    }
    double GetssNnsynsub2789(){
     return ssNnsynsub[2789];
    }
    double GetssNnsynsub2790(){
     return ssNnsynsub[2790];
    }
    double GetssNnsynsub2791(){
     return ssNnsynsub[2791];
    }
    double GetssNnsynsub2792(){
     return ssNnsynsub[2792];
    }
    double GetssNnsynsub2793(){
     return ssNnsynsub[2793];
    }
    double GetssNnsynsub2794(){
     return ssNnsynsub[2794];
    }
    double GetssNnsynsub2795(){
     return ssNnsynsub[2795];
    }
    double GetssNnsynsub2796(){
     return ssNnsynsub[2796];
    }
    double GetssNnsynsub2797(){
     return ssNnsynsub[2797];
    }
    double GetssNnsynsub2798(){
     return ssNnsynsub[2798];
    }
    double GetssNnsynsub2799(){
     return ssNnsynsub[2799];
    }
    double GetssNnsynsub2800(){
     return ssNnsynsub[2800];
    }
    double GetssNnsynsub2801(){
     return ssNnsynsub[2801];
    }
    double GetssNnsynsub2802(){
     return ssNnsynsub[2802];
    }
    double GetssNnsynsub2803(){
     return ssNnsynsub[2803];
    }
    double GetssNnsynsub2804(){
     return ssNnsynsub[2804];
    }
    double GetssNnsynsub2805(){
     return ssNnsynsub[2805];
    }
    double GetssNnsynsub2806(){
     return ssNnsynsub[2806];
    }
    double GetssNnsynsub2807(){
     return ssNnsynsub[2807];
    }
    double GetssNnsynsub2808(){
     return ssNnsynsub[2808];
    }
    double GetssNnsynsub2809(){
     return ssNnsynsub[2809];
    }
    double GetssNnsynsub2810(){
     return ssNnsynsub[2810];
    }
    double GetssNnsynsub2811(){
     return ssNnsynsub[2811];
    }
    double GetssNnsynsub2812(){
     return ssNnsynsub[2812];
    }
    double GetssNnsynsub2813(){
     return ssNnsynsub[2813];
    }
    double GetssNnsynsub2814(){
     return ssNnsynsub[2814];
    }
    double GetssNnsynsub2815(){
     return ssNnsynsub[2815];
    }
    double GetssNnsynsub2816(){
     return ssNnsynsub[2816];
    }
    double GetssNnsynsub2817(){
     return ssNnsynsub[2817];
    }
    double GetssNnsynsub2818(){
     return ssNnsynsub[2818];
    }
    double GetssNnsynsub2819(){
     return ssNnsynsub[2819];
    }
    double GetssNnsynsub2820(){
     return ssNnsynsub[2820];
    }
    double GetssNnsynsub2821(){
     return ssNnsynsub[2821];
    }
    double GetssNnsynsub2822(){
     return ssNnsynsub[2822];
    }
    double GetssNnsynsub2823(){
     return ssNnsynsub[2823];
    }
    double GetssNnsynsub2824(){
     return ssNnsynsub[2824];
    }
    double GetssNnsynsub2825(){
     return ssNnsynsub[2825];
    }
    double GetssNnsynsub2826(){
     return ssNnsynsub[2826];
    }
    double GetssNnsynsub2827(){
     return ssNnsynsub[2827];
    }
    double GetssNnsynsub2828(){
     return ssNnsynsub[2828];
    }
    double GetssNnsynsub2829(){
     return ssNnsynsub[2829];
    }
    double GetssNnsynsub2830(){
     return ssNnsynsub[2830];
    }
    double GetssNnsynsub2831(){
     return ssNnsynsub[2831];
    }
    double GetssNnsynsub2832(){
     return ssNnsynsub[2832];
    }
    double GetssNnsynsub2833(){
     return ssNnsynsub[2833];
    }
    double GetssNnsynsub2834(){
     return ssNnsynsub[2834];
    }
    double GetssNnsynsub2835(){
     return ssNnsynsub[2835];
    }
    double GetssNnsynsub2836(){
     return ssNnsynsub[2836];
    }
    double GetssNnsynsub2837(){
     return ssNnsynsub[2837];
    }
    double GetssNnsynsub2838(){
     return ssNnsynsub[2838];
    }
    double GetssNnsynsub2839(){
     return ssNnsynsub[2839];
    }
    double GetssNnsynsub2840(){
     return ssNnsynsub[2840];
    }
    double GetssNnsynsub2841(){
     return ssNnsynsub[2841];
    }
    double GetssNnsynsub2842(){
     return ssNnsynsub[2842];
    }
    double GetssNnsynsub2843(){
     return ssNnsynsub[2843];
    }
    double GetssNnsynsub2844(){
     return ssNnsynsub[2844];
    }
    double GetssNnsynsub2845(){
     return ssNnsynsub[2845];
    }
    double GetssNnsynsub2846(){
     return ssNnsynsub[2846];
    }
    double GetssNnsynsub2847(){
     return ssNnsynsub[2847];
    }
    double GetssNnsynsub2848(){
     return ssNnsynsub[2848];
    }
    double GetssNnsynsub2849(){
     return ssNnsynsub[2849];
    }
    double GetssNnsynsub2850(){
     return ssNnsynsub[2850];
    }
    double GetssNnsynsub2851(){
     return ssNnsynsub[2851];
    }
    double GetssNnsynsub2852(){
     return ssNnsynsub[2852];
    }
    double GetssNnsynsub2853(){
     return ssNnsynsub[2853];
    }
    double GetssNnsynsub2854(){
     return ssNnsynsub[2854];
    }
    double GetssNnsynsub2855(){
     return ssNnsynsub[2855];
    }
    double GetssNnsynsub2856(){
     return ssNnsynsub[2856];
    }
    double GetssNnsynsub2857(){
     return ssNnsynsub[2857];
    }
    double GetssNnsynsub2858(){
     return ssNnsynsub[2858];
    }
    double GetssNnsynsub2859(){
     return ssNnsynsub[2859];
    }
    double GetssNnsynsub2860(){
     return ssNnsynsub[2860];
    }
    double GetssNnsynsub2861(){
     return ssNnsynsub[2861];
    }
    double GetssNnsynsub2862(){
     return ssNnsynsub[2862];
    }
    double GetssNnsynsub2863(){
     return ssNnsynsub[2863];
    }
    double GetssNnsynsub2864(){
     return ssNnsynsub[2864];
    }
    double GetssNnsynsub2865(){
     return ssNnsynsub[2865];
    }
    double GetssNnsynsub2866(){
     return ssNnsynsub[2866];
    }
    double GetssNnsynsub2867(){
     return ssNnsynsub[2867];
    }
    double GetssNnsynsub2868(){
     return ssNnsynsub[2868];
    }
    double GetssNnsynsub2869(){
     return ssNnsynsub[2869];
    }
    double GetssNnsynsub2870(){
     return ssNnsynsub[2870];
    }
    double GetssNnsynsub2871(){
     return ssNnsynsub[2871];
    }
    double GetssNnsynsub2872(){
     return ssNnsynsub[2872];
    }
    double GetssNnsynsub2873(){
     return ssNnsynsub[2873];
    }
    double GetssNnsynsub2874(){
     return ssNnsynsub[2874];
    }
    double GetssNnsynsub2875(){
     return ssNnsynsub[2875];
    }
    double GetssNnsynsub2876(){
     return ssNnsynsub[2876];
    }
    double GetssNnsynsub2877(){
     return ssNnsynsub[2877];
    }
    double GetssNnsynsub2878(){
     return ssNnsynsub[2878];
    }
    double GetssNnsynsub2879(){
     return ssNnsynsub[2879];
    }
    double GetssNnsynsub2880(){
     return ssNnsynsub[2880];
    }
    double GetssNnsynsub2881(){
     return ssNnsynsub[2881];
    }
    double GetssNnsynsub2882(){
     return ssNnsynsub[2882];
    }
    double GetssNnsynsub2883(){
     return ssNnsynsub[2883];
    }
    double GetssNnsynsub2884(){
     return ssNnsynsub[2884];
    }
    double GetssNnsynsub2885(){
     return ssNnsynsub[2885];
    }
    double GetssNnsynsub2886(){
     return ssNnsynsub[2886];
    }
    double GetssNnsynsub2887(){
     return ssNnsynsub[2887];
    }
    double GetssNnsynsub2888(){
     return ssNnsynsub[2888];
    }
    double GetssNnsynsub2889(){
     return ssNnsynsub[2889];
    }
    double GetssNnsynsub2890(){
     return ssNnsynsub[2890];
    }
    double GetssNnsynsub2891(){
     return ssNnsynsub[2891];
    }
    double GetssNnsynsub2892(){
     return ssNnsynsub[2892];
    }
    double GetssNnsynsub2893(){
     return ssNnsynsub[2893];
    }
    double GetssNnsynsub2894(){
     return ssNnsynsub[2894];
    }
    double GetssNnsynsub2895(){
     return ssNnsynsub[2895];
    }
    double GetssNnsynsub2896(){
     return ssNnsynsub[2896];
    }
    double GetssNnsynsub2897(){
     return ssNnsynsub[2897];
    }
    double GetssNnsynsub2898(){
     return ssNnsynsub[2898];
    }
    double GetssNnsynsub2899(){
     return ssNnsynsub[2899];
    }
    double GetssNnsynsub2900(){
     return ssNnsynsub[2900];
    }
    double GetssNnsynsub2901(){
     return ssNnsynsub[2901];
    }
    double GetssNnsynsub2902(){
     return ssNnsynsub[2902];
    }
    double GetssNnsynsub2903(){
     return ssNnsynsub[2903];
    }
    double GetssNnsynsub2904(){
     return ssNnsynsub[2904];
    }
    double GetssNnsynsub2905(){
     return ssNnsynsub[2905];
    }
    double GetssNnsynsub2906(){
     return ssNnsynsub[2906];
    }
    double GetssNnsynsub2907(){
     return ssNnsynsub[2907];
    }
    double GetssNnsynsub2908(){
     return ssNnsynsub[2908];
    }
    double GetssNnsynsub2909(){
     return ssNnsynsub[2909];
    }
    double GetssNnsynsub2910(){
     return ssNnsynsub[2910];
    }
    double GetssNnsynsub2911(){
     return ssNnsynsub[2911];
    }
    double GetssNnsynsub2912(){
     return ssNnsynsub[2912];
    }
    double GetssNnsynsub2913(){
     return ssNnsynsub[2913];
    }
    double GetssNnsynsub2914(){
     return ssNnsynsub[2914];
    }
    double GetssNnsynsub2915(){
     return ssNnsynsub[2915];
    }
    double GetssNnsynsub2916(){
     return ssNnsynsub[2916];
    }
    double GetssNnsynsub2917(){
     return ssNnsynsub[2917];
    }
    double GetssNnsynsub2918(){
     return ssNnsynsub[2918];
    }
    double GetssNnsynsub2919(){
     return ssNnsynsub[2919];
    }
    double GetssNnsynsub2920(){
     return ssNnsynsub[2920];
    }
    double GetssNnsynsub2921(){
     return ssNnsynsub[2921];
    }
    double GetssNnsynsub2922(){
     return ssNnsynsub[2922];
    }
    double GetssNnsynsub2923(){
     return ssNnsynsub[2923];
    }
    double GetssNnsynsub2924(){
     return ssNnsynsub[2924];
    }
    double GetssNnsynsub2925(){
     return ssNnsynsub[2925];
    }
    double GetssNnsynsub2926(){
     return ssNnsynsub[2926];
    }
    double GetssNnsynsub2927(){
     return ssNnsynsub[2927];
    }
    double GetssNnsynsub2928(){
     return ssNnsynsub[2928];
    }
    double GetssNnsynsub2929(){
     return ssNnsynsub[2929];
    }
    double GetssNnsynsub2930(){
     return ssNnsynsub[2930];
    }
    double GetssNnsynsub2931(){
     return ssNnsynsub[2931];
    }
    double GetssNnsynsub2932(){
     return ssNnsynsub[2932];
    }
    double GetssNnsynsub2933(){
     return ssNnsynsub[2933];
    }
    double GetssNnsynsub2934(){
     return ssNnsynsub[2934];
    }
    double GetssNnsynsub2935(){
     return ssNnsynsub[2935];
    }
    double GetssNnsynsub2936(){
     return ssNnsynsub[2936];
    }
    double GetssNnsynsub2937(){
     return ssNnsynsub[2937];
    }
    double GetssNnsynsub2938(){
     return ssNnsynsub[2938];
    }
    double GetssNnsynsub2939(){
     return ssNnsynsub[2939];
    }
    double GetssNnsynsub2940(){
     return ssNnsynsub[2940];
    }
    double GetssNnsynsub2941(){
     return ssNnsynsub[2941];
    }
    double GetssNnsynsub2942(){
     return ssNnsynsub[2942];
    }
    double GetssNnsynsub2943(){
     return ssNnsynsub[2943];
    }
    double GetssNnsynsub2944(){
     return ssNnsynsub[2944];
    }
    double GetssNnsynsub2945(){
     return ssNnsynsub[2945];
    }
    double GetssNnsynsub2946(){
     return ssNnsynsub[2946];
    }
    double GetssNnsynsub2947(){
     return ssNnsynsub[2947];
    }
    double GetssNnsynsub2948(){
     return ssNnsynsub[2948];
    }
    double GetssNnsynsub2949(){
     return ssNnsynsub[2949];
    }
    double GetssNnsynsub2950(){
     return ssNnsynsub[2950];
    }
    double GetssNnsynsub2951(){
     return ssNnsynsub[2951];
    }
    double GetssNnsynsub2952(){
     return ssNnsynsub[2952];
    }
    double GetssNnsynsub2953(){
     return ssNnsynsub[2953];
    }
    double GetssNnsynsub2954(){
     return ssNnsynsub[2954];
    }
    double GetssNnsynsub2955(){
     return ssNnsynsub[2955];
    }
    double GetssNnsynsub2956(){
     return ssNnsynsub[2956];
    }
    double GetssNnsynsub2957(){
     return ssNnsynsub[2957];
    }
    double GetssNnsynsub2958(){
     return ssNnsynsub[2958];
    }
    double GetssNnsynsub2959(){
     return ssNnsynsub[2959];
    }
    double GetssNnsynsub2960(){
     return ssNnsynsub[2960];
    }
    double GetssNnsynsub2961(){
     return ssNnsynsub[2961];
    }
    double GetssNnsynsub2962(){
     return ssNnsynsub[2962];
    }
    double GetssNnsynsub2963(){
     return ssNnsynsub[2963];
    }
    double GetssNnsynsub2964(){
     return ssNnsynsub[2964];
    }
    double GetssNnsynsub2965(){
     return ssNnsynsub[2965];
    }
    double GetssNnsynsub2966(){
     return ssNnsynsub[2966];
    }
    double GetssNnsynsub2967(){
     return ssNnsynsub[2967];
    }
    double GetssNnsynsub2968(){
     return ssNnsynsub[2968];
    }
    double GetssNnsynsub2969(){
     return ssNnsynsub[2969];
    }
    double GetssNnsynsub2970(){
     return ssNnsynsub[2970];
    }
    double GetssNnsynsub2971(){
     return ssNnsynsub[2971];
    }
    double GetssNnsynsub2972(){
     return ssNnsynsub[2972];
    }
    double GetssNnsynsub2973(){
     return ssNnsynsub[2973];
    }
    double GetssNnsynsub2974(){
     return ssNnsynsub[2974];
    }
    double GetssNnsynsub2975(){
     return ssNnsynsub[2975];
    }
    double GetssNnsynsub2976(){
     return ssNnsynsub[2976];
    }
    double GetssNnsynsub2977(){
     return ssNnsynsub[2977];
    }
    double GetssNnsynsub2978(){
     return ssNnsynsub[2978];
    }
    double GetssNnsynsub2979(){
     return ssNnsynsub[2979];
    }
    double GetssNnsynsub2980(){
     return ssNnsynsub[2980];
    }
    double GetssNnsynsub2981(){
     return ssNnsynsub[2981];
    }
    double GetssNnsynsub2982(){
     return ssNnsynsub[2982];
    }
    double GetssNnsynsub2983(){
     return ssNnsynsub[2983];
    }
    double GetssNnsynsub2984(){
     return ssNnsynsub[2984];
    }
    double GetssNnsynsub2985(){
     return ssNnsynsub[2985];
    }
    double GetssNnsynsub2986(){
     return ssNnsynsub[2986];
    }
    double GetssNnsynsub2987(){
     return ssNnsynsub[2987];
    }
    double GetssNnsynsub2988(){
     return ssNnsynsub[2988];
    }
    double GetssNnsynsub2989(){
     return ssNnsynsub[2989];
    }
    double GetssNnsynsub2990(){
     return ssNnsynsub[2990];
    }
    double GetssNnsynsub2991(){
     return ssNnsynsub[2991];
    }
    double GetssNnsynsub2992(){
     return ssNnsynsub[2992];
    }
    double GetssNnsynsub2993(){
     return ssNnsynsub[2993];
    }
    double GetssNnsynsub2994(){
     return ssNnsynsub[2994];
    }
    double GetssNnsynsub2995(){
     return ssNnsynsub[2995];
    }
    double GetssNnsynsub2996(){
     return ssNnsynsub[2996];
    }
    double GetssNnsynsub2997(){
     return ssNnsynsub[2997];
    }
    double GetssNnsynsub2998(){
     return ssNnsynsub[2998];
    }
    double GetssNnsynsub2999(){
     return ssNnsynsub[2999];
    }


};

#endif // EVOLHISTSTATISTICS_H
