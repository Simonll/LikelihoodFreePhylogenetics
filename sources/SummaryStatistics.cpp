#include "SummaryStatistics.h"


SummaryStatistics::SummaryStatistics(LocalParameters * lparam){
// speudo count are used to ensure the possibitlity of using log standardization

        this->lparam = lparam;

        //Statistics containers
        nuc_usage = new double [lparam->Nnucp];
        nuc1_usage = new double [lparam->Nnucp];
        nuc2_usage = new double [lparam->Nnucp];
        nuc3_usage = new double [lparam->Nnucp];
        aa_usage = new double [lparam->Nstate_aa];
        aa_usagewo_nr = new double [lparam->Nstate_aa];
        codon_usage = new double [lparam->Nstate_codon];

        for (int nuc = 0 ; nuc <lparam->Nnucp; nuc++) {
            nuc_usage[nuc] = 1.0;
            nuc1_usage[nuc] = 1.0;
            nuc2_usage[nuc] = 1.0;
            nuc3_usage[nuc] = 1.0;
        }

        for (int aa =0 ; aa < lparam->Nstate_aa ; aa++) {
            aa_usage[aa] = 1.0;
            aa_usagewo_nr[aa] = 1.0;
        }

        for (int codon =0 ; codon < lparam->Nstate_codon ; codon++) {
            codon_usage[codon] = 1.0;
        }

        dinuc_usage = new double* [lparam->Nnucp];
        dinuc12_usage = new double* [lparam->Nnucp];
        dinuc23_usage = new double* [lparam->Nnucp];
        dinuc31_usage = new double* [lparam->Nnucp];

        for(int nuc = 0 ; nuc < lparam->Nnucp ; nuc++) {
            dinuc_usage[nuc] = new double[lparam->Nnucp];
            dinuc12_usage[nuc] = new double[lparam->Nnucp];
            dinuc23_usage[nuc] = new double[lparam->Nnucp];
            dinuc31_usage[nuc] = new double[lparam->Nnucp];
        }

        for(int nuc1 = 0 ; nuc1 < lparam->Nnucp ; nuc1++) {
            for(int nuc2 = 0 ; nuc2 < lparam->Nnucp ; nuc2++) {
                dinuc_usage[nuc1][nuc2] = 1.0;
                dinuc12_usage[nuc1][nuc2] = 1.0;
                dinuc23_usage[nuc1][nuc2] = 1.0;
                dinuc31_usage[nuc1][nuc2] = 1.0;
            }
        }

        diaa_usage = new double*[lparam->Nstate_aa];

        for (int aa = 0 ; aa < lparam->Nstate_aa; aa++){
            diaa_usage[aa] = new double [lparam->Nstate_aa];

        }

        for (int aa1 = 0 ; aa1 < lparam->Nstate_aa; aa1++){
            for (int aa2 = 0 ; aa2 < lparam->Nstate_aa; aa2++){
                diaa_usage[aa1][aa2] = 1.0;

            }
        }

        dicodon_usage = new double* [lparam->Nstate_codon];

        for(int codon = 0 ; codon < lparam->Nstate_codon; codon++) {
            dicodon_usage[codon] = new double [lparam->Nstate_codon];

        }

        for(int codon1 = 0 ; codon1 < lparam->Nstate_codon; codon1++) {
            for(int codon2 = 0 ; codon2 < lparam->Nstate_codon; codon2++) {
                dicodon_usage[codon1][codon2] = 1.0;

            }
        }

        nuc_meandiff= new double [2];
        nuc1_meandiff= new double [2];
        nuc2_meandiff= new double [2];
        nuc3_meandiff= new double [2];
        codon_meandiff = new double [2];
        aa_meandiff = new double [2];
        aa_wonR_meandiff = new double [2];
        CGNAGR = new double[4];

        for (int i = 0 ; i < 2; i++)  {
            nuc_meandiff[i] = 0.0;
            nuc1_meandiff[i] = 0.0;
            nuc2_meandiff[i] = 0.0;
            nuc3_meandiff[i] = 0.0;
            codon_meandiff[i] = 0.0;
            aa_meandiff[i] = 0.0;
            aa_wonR_meandiff[i] = 0.0;

        }

        for (int i = 0 ; i < 4 ; i++ ) {
            CGNAGR[i] = 0.0;

        }

        nuc_pairwise = new int [7];
        nuc1_pairwise = new int [7];
        nuc2_pairwise = new int [7];
        nuc3_pairwise = new int [7];

        for (int i = 0 ; i < 7 ; i++) {
            nuc_pairwise[i] = 0;
            nuc1_pairwise[i] = 0;
            nuc2_pairwise[i] = 0;
            nuc3_pairwise[i] = 0;
        }

        aa_pairwise = new int [190];

        for (int i = 0 ; i < 190; i ++) {
            aa_pairwise[i] = 0 ;
        }



        dinucCpG_pairwise = new int [4];

        for (int i = 0 ; i <4; i++){
            dinucCpG_pairwise[i] = 0;
        }


        nuc_site_comphet = 0.0;
        nuc1_site_comphet = 0.0;
        nuc2_site_comphet = 0.0;
        nuc3_site_comphet = 0.0;
        nuc_taxa_comphet = 0.0;
        nuc1_taxa_comphet = 0.0;
        nuc2_taxa_comphet = 0.0;
        nuc3_taxa_comphet = 0.0;
        codon_site_comphet = 0.0;
        codon_taxa_comphet = 0.0;
        aa_site_comphet = 0.0;
        aa_taxa_comphet = 0.0;




        //ss.computeSummaries();

        codon_bool = false;
        dinuc_bool = false;
        dinuc12_bool = false;
        dinuc23_bool = false;
        dinuc31_bool = false;
        aa_bool = false;
        aa_wo_nr = false;
        dicodon_bool = false;
        diaa_bool = false;
        nuc_bool = false;
        nuc1_bool = false;
        nuc2_bool = false;
        nuc3_bool = false;
        nuc_meandiff_bool = false;
        nuc1_meandiff_bool = false;
        nuc2_meandiff_bool = false;
        nuc3_meandiff_bool = false;
        codon_meandiff_bool = false;
        aa_meandiff_bool = false;
        aa_wonR_meandiff_bool = false;
        CGNAGR_bool = false;
        nuc_pairwise_bool = false;
        nuc1_pairwise_bool = false;
        nuc2_pairwise_bool = false;
        nuc3_pairwise_bool = false;
        aa_pairwise_bool = false;
        dinucCpG_pairwise_bool = false;
        nuc_site_comphet_bool = false;
        nuc1_site_comphet_bool = false;
        nuc2_site_comphet_bool = false;
        nuc3_site_comphet_bool = false;
        nuc_taxa_comphet_bool = false;
        nuc1_taxa_comphet_bool = false;
        nuc2_taxa_comphet_bool = false;
        nuc3_taxa_comphet_bool = false;
        codon_site_comphet_bool = false;
        codon_taxa_comphet_bool = false;
        aa_site_comphet_bool = false;
        aa_taxa_comphet_bool = false;

        MapFunctions();


}

SummaryStatistics::SummaryStatistics(LocalData *ldata){
// speudo count are used to ensure the possibitlity of using log standardization

       this->ldata = ldata;

        //Statistics containers
       nuc_usage = new double [ldata->Nnucp];
       nuc1_usage = new double [ldata->Nnucp];
       nuc2_usage = new double [ldata->Nnucp];
       nuc3_usage = new double [ldata->Nnucp];
       aa_usage = new double [ldata->Nstate_aa];
       aa_usagewo_nr = new double [ldata->Nstate_aa];

       codon_usage = new double [ldata->Nstate_codon];

       for (int nuc = 0 ; nuc <ldata->Nnucp; nuc++) {
           nuc_usage[nuc] = 1.0;
           nuc1_usage[nuc] = 1.0;
           nuc2_usage[nuc] = 1.0;
           nuc3_usage[nuc] = 1.0;
       }

       for (int aa =0 ; aa < ldata->Nstate_aa ; aa++) {
           aa_usage[aa] = 1.0;
           aa_usagewo_nr[aa] = 1.0;
       }

       for (int codon =0 ; codon < ldata->Nstate_codon ; codon++) {
           codon_usage[codon] = 1.0;
       }

       dinuc_usage = new double* [ldata->Nnucp];
       dinuc12_usage = new double* [ldata->Nnucp];
       dinuc23_usage = new double* [ldata->Nnucp];
       dinuc31_usage = new double* [ldata->Nnucp];

       for(int nuc = 0 ; nuc < ldata->Nnucp ; nuc++) {
           dinuc_usage[nuc] = new double[ldata->Nnucp];
           dinuc12_usage[nuc] = new double[ldata->Nnucp];
           dinuc23_usage[nuc] = new double[ldata->Nnucp];
           dinuc31_usage[nuc] = new double[ldata->Nnucp];
       }

       for(int nuc1 = 0 ; nuc1 < ldata->Nnucp ; nuc1++) {
           for(int nuc2 = 0 ; nuc2 < ldata->Nnucp ; nuc2++) {
               dinuc_usage[nuc1][nuc2] = 1.0;
               dinuc12_usage[nuc1][nuc2] = 1.0;
               dinuc23_usage[nuc1][nuc2] = 1.0;
               dinuc31_usage[nuc1][nuc2] = 1.0;
           }
       }

       diaa_usage = new double*[ldata->Nstate_aa];

       for (int aa = 0 ; aa < ldata->Nstate_aa; aa++){
           diaa_usage[aa] = new double [ldata->Nstate_aa];

       }

       for (int aa1 = 0 ; aa1 < ldata->Nstate_aa; aa1++){
           for (int aa2 = 0 ; aa2 < ldata->Nstate_aa; aa2++){
               diaa_usage[aa1][aa2] = 1.0;

           }
       }

       dicodon_usage = new double* [ldata->Nstate_codon];

       for(int codon = 0 ; codon < ldata->Nstate_codon; codon++) {
           dicodon_usage[codon] = new double [ldata->Nstate_codon];

       }

       for(int codon1 = 0 ; codon1 < ldata->Nstate_codon; codon1++) {
           for(int codon2 = 0 ; codon2 < ldata->Nstate_codon; codon2++) {
               dicodon_usage[codon1][codon2] = 1.0;

           }
       }

       nuc_meandiff= new double [2];
       nuc1_meandiff= new double [2];
       nuc2_meandiff= new double [2];
       nuc3_meandiff= new double [2];
       codon_meandiff = new double [2];
       aa_meandiff = new double [2];
       aa_wonR_meandiff = new double [2];
       CGNAGR = new double [4];
       for (int i = 0 ; i < 2; i++)  {
           nuc_meandiff[i] = 0.0;
           nuc1_meandiff[i] = 0.0;
           nuc2_meandiff[i] = 0.0;
           nuc3_meandiff[i] = 0.0;
           codon_meandiff[i] = 0.0;
           aa_meandiff[i] = 0.0;
           aa_wonR_meandiff[i] = 0.0;

       }

       for (int i = 0 ; i < 4 ; i ++ ) {
            CGNAGR[i] = 0.0;

       }
       nuc_pairwise = new int [7];
       nuc1_pairwise = new int [7];
       nuc2_pairwise = new int [7];
       nuc3_pairwise = new int [7];

       for (int i = 0 ; i < 7 ; i++) {
           nuc_pairwise[i] = 0;
           nuc1_pairwise[i] = 0;
           nuc2_pairwise[i] = 0;
           nuc3_pairwise[i] = 0;
       }

       aa_pairwise = new int [190];

       for (int i = 0 ; i < 190; i ++) {
           aa_pairwise[i] = 0 ;
       }



       dinucCpG_pairwise = new int [4];

       for (int i = 0 ; i <4; i++){
           dinucCpG_pairwise[i] = 0;
       }


       nuc_site_comphet = 0.0;
       nuc1_site_comphet = 0.0;
       nuc2_site_comphet = 0.0;
       nuc3_site_comphet = 0.0;
       nuc_taxa_comphet = 0.0;
       nuc1_taxa_comphet = 0.0;
       nuc2_taxa_comphet = 0.0;
       nuc3_taxa_comphet = 0.0;
       codon_site_comphet = 0.0;
       codon_taxa_comphet = 0.0;
       aa_site_comphet = 0.0;
       aa_taxa_comphet = 0.0;

        cerr << "AAA22\n";


       //ss.computeSummaries();

       codon_bool = false;
       dinuc_bool = false;
       dinuc12_bool = false;
       dinuc23_bool = false;
       dinuc31_bool = false;
       aa_bool = false;
       aa_wo_nr = false;
       dicodon_bool = false;
       diaa_bool = false;
       nuc_bool = false;
       nuc1_bool = false;
       nuc2_bool = false;
       nuc3_bool = false;
       nuc_meandiff_bool = false;
       nuc1_meandiff_bool = false;
       nuc2_meandiff_bool = false;
       nuc3_meandiff_bool = false;
       codon_meandiff_bool = false;
       aa_meandiff_bool = false;
       aa_wonR_meandiff_bool = false;
       CGNAGR_bool = false;
       nuc_pairwise_bool = false;
       nuc1_pairwise_bool = false;
       nuc2_pairwise_bool = false;
       nuc3_pairwise_bool = false;
       aa_pairwise_bool = false;
       dinucCpG_pairwise_bool = false;
       nuc_site_comphet_bool = false;
       nuc1_site_comphet_bool = false;
       nuc2_site_comphet_bool = false;
       nuc3_site_comphet_bool = false;
       nuc_taxa_comphet_bool = false;
       nuc1_taxa_comphet_bool = false;
       nuc2_taxa_comphet_bool = false;
       nuc3_taxa_comphet_bool = false;
       codon_site_comphet_bool = false;
       codon_taxa_comphet_bool = false;
       aa_site_comphet_bool = false;
       aa_taxa_comphet_bool = false;
        cerr << "AAA23\n";
       MapFunctions();

        cerr << "AAA24\n";
}

void SummaryStatistics::MapFunctions(){

        GetSummariesMap["pwAC"] = &SummaryStatistics::GetpwAC;
        GetSummariesMap["pwAG"] = &SummaryStatistics::GetpwAG;
        GetSummariesMap["pwAT"] = &SummaryStatistics::GetpwAT;
        GetSummariesMap["pwCG"] = &SummaryStatistics::GetpwCG;
        GetSummariesMap["pwCT"] = &SummaryStatistics::GetpwCT;
        GetSummariesMap["pwGT"] = &SummaryStatistics::GetpwGT;

        GetSummariesMap["pw1AC"] = &SummaryStatistics::Getpw1AC;
        GetSummariesMap["pw1AG"] = &SummaryStatistics::Getpw1AG;
        GetSummariesMap["pw1AT"] = &SummaryStatistics::Getpw1AT;
        GetSummariesMap["pw1CG"] = &SummaryStatistics::Getpw1CG;
        GetSummariesMap["pw1CT"] = &SummaryStatistics::Getpw1CT;
        GetSummariesMap["pw1GT"] = &SummaryStatistics::Getpw1GT;

        GetSummariesMap["pw2AC"] = &SummaryStatistics::Getpw2AC;
        GetSummariesMap["pw2AG"] = &SummaryStatistics::Getpw2AG;
        GetSummariesMap["pw2AT"] = &SummaryStatistics::Getpw2AT;
        GetSummariesMap["pw2CG"] = &SummaryStatistics::Getpw2CG;
        GetSummariesMap["pw2CT"] = &SummaryStatistics::Getpw2CT;
        GetSummariesMap["pw2GT"] = &SummaryStatistics::Getpw2GT;

        GetSummariesMap["pw3AC"] = &SummaryStatistics::Getpw3AC;
        GetSummariesMap["pw3AG"] = &SummaryStatistics::Getpw3AG;
        GetSummariesMap["pw3AT"] = &SummaryStatistics::Getpw3AT;
        GetSummariesMap["pw3CG"] = &SummaryStatistics::Getpw3CG;
        GetSummariesMap["pw3CT"] = &SummaryStatistics::Getpw3CT;
        GetSummariesMap["pw3GT"] = &SummaryStatistics::Getpw3GT;

        GetSummariesMap["pwAA"]  = &SummaryStatistics::GetpwAA;

        GetSummariesMap["pwCpGTpG"]  = &SummaryStatistics::GetdinucCpG_TpG;
        GetSummariesMap["pwCpGCpA"]  = &SummaryStatistics::GetdinucCpG_CpA;
        GetSummariesMap["pwApGTpG"]  = &SummaryStatistics::GetdinucApG_TpG;


        GetSummariesMap["mAA"]  = &SummaryStatistics::GetAAmean;
        GetSummariesMap["vAA"]  = &SummaryStatistics::GetAAvar;

        GetSummariesMap["mAAwonR"]  = &SummaryStatistics::GetAAmean_wonR;
        GetSummariesMap["vAAwonR"]  = &SummaryStatistics::GetAAvar_wonR;

        GetSummariesMap["mCodon"]  = &SummaryStatistics::GetCodonmean;
        GetSummariesMap["vCodon"]  = &SummaryStatistics::GetCodonvar;

        GetSummariesMap["mNuc"]  = &SummaryStatistics::GetNucmean;
        GetSummariesMap["vNuc"]  = &SummaryStatistics::GetNucvar;

        GetSummariesMap["mNuc1"]  = &SummaryStatistics::GetNuc1mean;
        GetSummariesMap["vNuc1"]  = &SummaryStatistics::GetNuc1var;


        GetSummariesMap["mNuc2"]  = &SummaryStatistics::GetNuc2mean;
        GetSummariesMap["vNuc2"]  = &SummaryStatistics::GetNuc2var;


        GetSummariesMap["mNuc3"]  = &SummaryStatistics::GetNuc2mean;
        GetSummariesMap["vNuc3"]  = &SummaryStatistics::GetNuc2var;

        GetSummariesMap["CGNAGRcgnconst"]  = &SummaryStatistics::GetCGNAGRcgnconst;
        GetSummariesMap["CGNAGRagrconst"]  = &SummaryStatistics::GetCGNAGRagrconst;
        GetSummariesMap["CGNAGRcgnvar"]  = &SummaryStatistics::GetCGNAGRcgnvar;
        GetSummariesMap["CGNAGRagrvar"]  = &SummaryStatistics::GetCGNAGRagrvar;


        GetSummariesMap["nucA"] = &SummaryStatistics::GetNucA;
        GetSummariesMap["nucC"] = &SummaryStatistics::GetNucC;
        GetSummariesMap["nucG"] = &SummaryStatistics::GetNucG;
        GetSummariesMap["nucT"] = &SummaryStatistics::GetNucT;

        GetSummariesMap["nuc1A"] = &SummaryStatistics::GetNuc1A;
        GetSummariesMap["nuc1C"] = &SummaryStatistics::GetNuc1C;
        GetSummariesMap["nuc1G"] = &SummaryStatistics::GetNuc1G;
        GetSummariesMap["nuc1T"] = &SummaryStatistics::GetNuc1T;

        GetSummariesMap["nuc2A"] = &SummaryStatistics::GetNuc2A;
        GetSummariesMap["nuc2C"] = &SummaryStatistics::GetNuc2C;
        GetSummariesMap["nuc2G"] = &SummaryStatistics::GetNuc2G;
        GetSummariesMap["nuc2T"] = &SummaryStatistics::GetNuc2T;

        GetSummariesMap["nuc3A"] = &SummaryStatistics::GetNuc3A;
        GetSummariesMap["nuc3C"] = &SummaryStatistics::GetNuc3C;
        GetSummariesMap["nuc3G"] = &SummaryStatistics::GetNuc3G;
        GetSummariesMap["nuc3T"] = &SummaryStatistics::GetNuc3T;

        GetSummariesMap["dinucAA"] = &SummaryStatistics::GetDinucAA;
        GetSummariesMap["dinucAC"] = &SummaryStatistics::GetDinucAC;
        GetSummariesMap["dinucAG"] = &SummaryStatistics::GetDinucAG;
        GetSummariesMap["dinucAT"] = &SummaryStatistics::GetDinucAT;

        GetSummariesMap["dinucCA"] = &SummaryStatistics::GetDinucCA;
        GetSummariesMap["dinucCC"] = &SummaryStatistics::GetDinucCC;
        GetSummariesMap["dinucCG"] = &SummaryStatistics::GetDinucCG;
        GetSummariesMap["dinucCT"] = &SummaryStatistics::GetDinucCT;

        GetSummariesMap["dinucGA"] = &SummaryStatistics::GetDinucGA;
        GetSummariesMap["dinucGC"] = &SummaryStatistics::GetDinucGC;
        GetSummariesMap["dinucGG"] = &SummaryStatistics::GetDinucGG;
        GetSummariesMap["dinucGT"] = &SummaryStatistics::GetDinucGT;

        GetSummariesMap["dinucTA"] = &SummaryStatistics::GetDinucGA;
        GetSummariesMap["dinucTC"] = &SummaryStatistics::GetDinucGC;
        GetSummariesMap["dinucTG"] = &SummaryStatistics::GetDinucGG;
        GetSummariesMap["dinucTT"] = &SummaryStatistics::GetDinucGT;

        GetSummariesMap["dinuc12AA"] = &SummaryStatistics::GetDinuc12AA;
        GetSummariesMap["dinuc12AC"] = &SummaryStatistics::GetDinuc12AC;
        GetSummariesMap["dinuc12AG"] = &SummaryStatistics::GetDinuc12AG;
        GetSummariesMap["dinuc12AT"] = &SummaryStatistics::GetDinuc12AT;

        GetSummariesMap["dinuc12CA"] = &SummaryStatistics::GetDinuc12CA;
        GetSummariesMap["dinuc12CC"] = &SummaryStatistics::GetDinuc12CC;
        GetSummariesMap["dinuc12CG"] = &SummaryStatistics::GetDinuc12CG;
        GetSummariesMap["dinuc12CT"] = &SummaryStatistics::GetDinuc12CT;

        GetSummariesMap["dinuc12GA"] = &SummaryStatistics::GetDinuc12GA;
        GetSummariesMap["dinuc12GC"] = &SummaryStatistics::GetDinuc12GC;
        GetSummariesMap["dinuc12GG"] = &SummaryStatistics::GetDinuc12GG;
        GetSummariesMap["dinuc12GT"] = &SummaryStatistics::GetDinuc12GT;

        GetSummariesMap["dinuc12TA"] = &SummaryStatistics::GetDinuc12GA;
        GetSummariesMap["dinuc12TC"] = &SummaryStatistics::GetDinuc12GC;
        GetSummariesMap["dinuc12TG"] = &SummaryStatistics::GetDinuc12GG;
        GetSummariesMap["dinuc12TT"] = &SummaryStatistics::GetDinuc12GT;


        GetSummariesMap["dinuc23AA"] = &SummaryStatistics::GetDinuc23AA;
        GetSummariesMap["dinuc23AC"] = &SummaryStatistics::GetDinuc23AC;
        GetSummariesMap["dinuc23AG"] = &SummaryStatistics::GetDinuc23AG;
        GetSummariesMap["dinuc23AT"] = &SummaryStatistics::GetDinuc23AT;

        GetSummariesMap["dinuc23CA"] = &SummaryStatistics::GetDinuc23CA;
        GetSummariesMap["dinuc23CC"] = &SummaryStatistics::GetDinuc23CC;
        GetSummariesMap["dinuc23CG"] = &SummaryStatistics::GetDinuc23CG;
        GetSummariesMap["dinuc23CT"] = &SummaryStatistics::GetDinuc23CT;

        GetSummariesMap["dinuc23GA"] = &SummaryStatistics::GetDinuc23GA;
        GetSummariesMap["dinuc23GC"] = &SummaryStatistics::GetDinuc23GC;
        GetSummariesMap["dinuc23GG"] = &SummaryStatistics::GetDinuc23GG;
        GetSummariesMap["dinuc23GT"] = &SummaryStatistics::GetDinuc23GT;

        GetSummariesMap["dinuc23TA"] = &SummaryStatistics::GetDinuc23GA;
        GetSummariesMap["dinuc23TC"] = &SummaryStatistics::GetDinuc23GC;
        GetSummariesMap["dinuc23TG"] = &SummaryStatistics::GetDinuc23GG;
        GetSummariesMap["dinuc23TT"] = &SummaryStatistics::GetDinuc23GT;


        GetSummariesMap["dinuc31AA"] = &SummaryStatistics::GetDinuc31AA;
        GetSummariesMap["dinuc31AC"] = &SummaryStatistics::GetDinuc31AC;
        GetSummariesMap["dinuc31AG"] = &SummaryStatistics::GetDinuc31AG;
        GetSummariesMap["dinuc31AT"] = &SummaryStatistics::GetDinuc31AT;

        GetSummariesMap["dinuc31CA"] = &SummaryStatistics::GetDinuc31CA;
        GetSummariesMap["dinuc31CC"] = &SummaryStatistics::GetDinuc31CC;
        GetSummariesMap["dinuc31CG"] = &SummaryStatistics::GetDinuc31CG;
        GetSummariesMap["dinuc31CT"] = &SummaryStatistics::GetDinuc31CT;

        GetSummariesMap["dinuc31GA"] = &SummaryStatistics::GetDinuc31GA;
        GetSummariesMap["dinuc31GC"] = &SummaryStatistics::GetDinuc31GC;
        GetSummariesMap["dinuc31GG"] = &SummaryStatistics::GetDinuc31GG;
        GetSummariesMap["dinuc31GT"] = &SummaryStatistics::GetDinuc31GT;

        GetSummariesMap["dinuc31TA"] = &SummaryStatistics::GetDinuc31GA;
        GetSummariesMap["dinuc31TC"] = &SummaryStatistics::GetDinuc31GC;
        GetSummariesMap["dinuc31TG"] = &SummaryStatistics::GetDinuc31GG;
        GetSummariesMap["dinuc31TT"] = &SummaryStatistics::GetDinuc31GT;


        GetSummariesMap["A"] = &SummaryStatistics::GetA;
        GetSummariesMap["C"] = &SummaryStatistics::GetC;
        GetSummariesMap["D"] = &SummaryStatistics::GetD;
        GetSummariesMap["E"] = &SummaryStatistics::GetE;
        GetSummariesMap["F"] = &SummaryStatistics::GetF;
        GetSummariesMap["G"] = &SummaryStatistics::GetG;
        GetSummariesMap["H"] = &SummaryStatistics::GetH;
        GetSummariesMap["I"] = &SummaryStatistics::GetI;
        GetSummariesMap["K"] = &SummaryStatistics::GetK;
        GetSummariesMap["L"] = &SummaryStatistics::GetL;
        GetSummariesMap["M"] = &SummaryStatistics::GetM;
        GetSummariesMap["N"] = &SummaryStatistics::GetN;
        GetSummariesMap["P"] = &SummaryStatistics::GetP;
        GetSummariesMap["Q"] = &SummaryStatistics::GetQ;
        GetSummariesMap["R"] = &SummaryStatistics::GetR;
        GetSummariesMap["S"] = &SummaryStatistics::GetS;
        GetSummariesMap["T"] = &SummaryStatistics::GetT;
        GetSummariesMap["V"] = &SummaryStatistics::GetV;
        GetSummariesMap["W"] = &SummaryStatistics::GetW;
        GetSummariesMap["Y"] = &SummaryStatistics::GetY;

//        GetSummariesMap["AwonR"] = &SummaryStatistics::GetAwonR;
//        GetSummariesMap["CwonR"] = &SummaryStatistics::GetCwonR;
//        GetSummariesMap["DwonR"] = &SummaryStatistics::GetDwonR;
//        GetSummariesMap["EwonR"] = &SummaryStatistics::GetEwonR;
//        GetSummariesMap["FwonR"] = &SummaryStatistics::GetFwonR;
//        GetSummariesMap["GwonR"] = &SummaryStatistics::GetGwonR;
//        GetSummariesMap["HwonR"] = &SummaryStatistics::GetHwonR;
//        GetSummariesMap["IwonR"] = &SummaryStatistics::GetIwonR;
//        GetSummariesMap["KwonR"] = &SummaryStatistics::GetKwonR;
//        GetSummariesMap["LwonR"] = &SummaryStatistics::GetLwonR;
//        GetSummariesMap["MwonR"] = &SummaryStatistics::GetMwonR;
//        GetSummariesMap["NwonR"] = &SummaryStatistics::GetNwonR;
//        GetSummariesMap["PwonR"] = &SummaryStatistics::GetPwonR;
//        GetSummariesMap["QwonR"] = &SummaryStatistics::GetQwonR;
//        GetSummariesMap["RwonR"] = &SummaryStatistics::GetRwonR;
//        GetSummariesMap["SwonR"] = &SummaryStatistics::GetSwonR;
//        GetSummariesMap["TwonR"] = &SummaryStatistics::GetTwonR;
//        GetSummariesMap["VwonR"] = &SummaryStatistics::GetVwonR;
//        GetSummariesMap["WwonR"] = &SummaryStatistics::GetWwonR;
//        GetSummariesMap["YwonR"] = &SummaryStatistics::GetYwonR;


        GetSummariesMap["TTT"] = &SummaryStatistics::GetTTT;
        GetSummariesMap["TTC"] = &SummaryStatistics::GetTTC;
        GetSummariesMap["TTA"] = &SummaryStatistics::GetTTA;
        GetSummariesMap["TTG"] = &SummaryStatistics::GetTTG;
        GetSummariesMap["TCT"] = &SummaryStatistics::GetTCT;
        GetSummariesMap["TCC"] = &SummaryStatistics::GetTCC;
        GetSummariesMap["TCA"] = &SummaryStatistics::GetTCA;
        GetSummariesMap["TCG"] = &SummaryStatistics::GetTCG;
        GetSummariesMap["TAT"] = &SummaryStatistics::GetTAT;
        GetSummariesMap["TAC"] = &SummaryStatistics::GetTAC;
        GetSummariesMap["TAA"] = &SummaryStatistics::GetTAA;
        GetSummariesMap["TAG"] = &SummaryStatistics::GetTAG;
        GetSummariesMap["TGT"] = &SummaryStatistics::GetTGT;
        GetSummariesMap["TGC"] = &SummaryStatistics::GetTGC;
        GetSummariesMap["TGA"] = &SummaryStatistics::GetTGA;
        GetSummariesMap["TGG"] = &SummaryStatistics::GetTGG;
        GetSummariesMap["CTT"] = &SummaryStatistics::GetCTT;
        GetSummariesMap["CTC"] = &SummaryStatistics::GetCTC;
        GetSummariesMap["CTA"] = &SummaryStatistics::GetCTA;
        GetSummariesMap["CTG"] = &SummaryStatistics::GetCTG;
        GetSummariesMap["CCT"] = &SummaryStatistics::GetCCT;
        GetSummariesMap["CCC"] = &SummaryStatistics::GetCCC;
        GetSummariesMap["CCA"] = &SummaryStatistics::GetCCA;
        GetSummariesMap["CCG"] = &SummaryStatistics::GetCCG;
        GetSummariesMap["CAT"] = &SummaryStatistics::GetCAT;
        GetSummariesMap["CAC"] = &SummaryStatistics::GetCAC;
        GetSummariesMap["CAA"] = &SummaryStatistics::GetCAA;
        GetSummariesMap["CAG"] = &SummaryStatistics::GetCAG;
        GetSummariesMap["CGT"] = &SummaryStatistics::GetCGT;
        GetSummariesMap["CGC"] = &SummaryStatistics::GetCGC;
        GetSummariesMap["CGA"] = &SummaryStatistics::GetCGA;
        GetSummariesMap["CGG"] = &SummaryStatistics::GetCGG;
        GetSummariesMap["ATT"] = &SummaryStatistics::GetATT;
        GetSummariesMap["ATC"] = &SummaryStatistics::GetATC;
        GetSummariesMap["ATA"] = &SummaryStatistics::GetATA;
        GetSummariesMap["ATG"] = &SummaryStatistics::GetATG;
        GetSummariesMap["ACT"] = &SummaryStatistics::GetACT;
        GetSummariesMap["ACC"] = &SummaryStatistics::GetACC;
        GetSummariesMap["ACA"] = &SummaryStatistics::GetACA;
        GetSummariesMap["ACG"] = &SummaryStatistics::GetACG;
        GetSummariesMap["AAT"] = &SummaryStatistics::GetAAT;
        GetSummariesMap["AAC"] = &SummaryStatistics::GetAAC;
        GetSummariesMap["AAA"] = &SummaryStatistics::GetAAA;
        GetSummariesMap["AAG"] = &SummaryStatistics::GetAAG;
        GetSummariesMap["AGT"] = &SummaryStatistics::GetAGT;
        GetSummariesMap["AGC"] = &SummaryStatistics::GetAGC;
        GetSummariesMap["AGA"] = &SummaryStatistics::GetAGA;
        GetSummariesMap["AGG"] = &SummaryStatistics::GetAGG;
        GetSummariesMap["GTT"] = &SummaryStatistics::GetGTT;
        GetSummariesMap["GTC"] = &SummaryStatistics::GetGTC;
        GetSummariesMap["GTA"] = &SummaryStatistics::GetGTA;
        GetSummariesMap["GTG"] = &SummaryStatistics::GetGTG;
        GetSummariesMap["GCT"] = &SummaryStatistics::GetGCT;
        GetSummariesMap["GCC"] = &SummaryStatistics::GetGCC;
        GetSummariesMap["GCA"] = &SummaryStatistics::GetGCA;
        GetSummariesMap["GCG"] = &SummaryStatistics::GetGCG;
        GetSummariesMap["GAT"] = &SummaryStatistics::GetGAT;
        GetSummariesMap["GAC"] = &SummaryStatistics::GetGAC;
        GetSummariesMap["GAA"] = &SummaryStatistics::GetGAA;
        GetSummariesMap["GAG"] = &SummaryStatistics::GetGAG;
        GetSummariesMap["GGT"] = &SummaryStatistics::GetGGT;
        GetSummariesMap["GGC"] = &SummaryStatistics::GetGGC;
        GetSummariesMap["GGA"] = &SummaryStatistics::GetGGA;
        GetSummariesMap["GGG"] = &SummaryStatistics::GetGGG;

        GetSummariesMap["nucsitecomphet"] = &SummaryStatistics::Getnuc_site_comphet;
        GetSummariesMap["nuc1sitecomphet"] = &SummaryStatistics::Getnuc1_site_comphet;
        GetSummariesMap["nuc2sitecomphet"] = &SummaryStatistics::Getnuc2_site_comphet;
        GetSummariesMap["nuc3sitecomphet"] = &SummaryStatistics::Getnuc3_site_comphet;
        GetSummariesMap["nuc_taxacomphet"] = &SummaryStatistics::Getnuc_taxa_comphet;
        GetSummariesMap["nuc1taxacomphet"] = &SummaryStatistics::Getnuc1_taxa_comphet;
        GetSummariesMap["nuc2taxacomphet"] = &SummaryStatistics::Getnuc2_taxa_comphet;
        GetSummariesMap["nuc3taxacomphet"] = &SummaryStatistics::Getnuc3_taxa_comphet;
        GetSummariesMap["codonsitecomphet"] = &SummaryStatistics::Getcodon_site_comphet;
        GetSummariesMap["codontaxacomphet"] = &SummaryStatistics::Getcodon_taxa_comphet;
        GetSummariesMap["aasitecomphet"] = &SummaryStatistics::Getaa_site_comphet;
        GetSummariesMap["aataxacomphet"] = &SummaryStatistics::Getaa_taxa_comphet;

}

SummaryStatistics::~SummaryStatistics()
{

    for(int nuc = 0 ; nuc < lparam->Nnucp ; nuc++) {
        delete [] dinuc_usage[nuc];
        delete [] dinuc12_usage[nuc];
        delete [] dinuc23_usage[nuc];
        delete [] dinuc31_usage[nuc];
    }

        delete [] dinuc_usage;
        delete [] dinuc12_usage;
        delete [] dinuc23_usage;
        delete [] dinuc31_usage;
        delete [] dicodon_usage;
        delete [] diaa_usage;
        delete [] aa_usage;
        delete [] aa_usagewo_nr;
        delete [] CGNAGR;
        delete [] nuc_usage;
        delete [] nuc1_usage;
        delete [] nuc2_usage;
        delete [] nuc3_usage;
        delete [] nuc_meandiff;
        delete [] nuc1_meandiff;
        delete [] nuc2_meandiff;
        delete [] nuc3_meandiff;
        delete [] codon_meandiff;
        delete [] aa_meandiff;
        delete [] aa_wonR_meandiff;
        delete [] nuc_pairwise;
        delete [] nuc1_pairwise;
        delete [] nuc2_pairwise;
        delete [] nuc3_pairwise;
        delete [] aa_pairwise;


}

//void SummaryStatistics::computeSummaries(CodonSequenceAlignment* codondata) {

//}




void SummaryStatistics::computeSummaries(int** CurrentNodeLeafCodonSequence){
    lparam->summariesSimulatedData.clear();
    lparam->summariesSimulatedData.shrink_to_fit();

    codon_bool = false;
    dinuc_bool = false;
    dinuc12_bool = false;
    dinuc23_bool = false;
    dinuc31_bool = false;
    aa_bool = false;
    aa_wo_nr = false;
    dicodon_bool = false;
    diaa_bool = false;
    nuc_bool = false;
    nuc1_bool = false;
    nuc2_bool = false;
    nuc3_bool = false;
    nuc_meandiff_bool = false;
    nuc1_meandiff_bool = false;
    nuc2_meandiff_bool = false;
    nuc3_meandiff_bool = false;
    codon_meandiff_bool = false;
    aa_meandiff_bool = false;
    aa_wonR_meandiff_bool = false;
    CGNAGR_bool = false;
    nuc_pairwise_bool = false;
    nuc1_pairwise_bool = false;
    nuc2_pairwise_bool = false;
    nuc3_pairwise_bool = false;
    aa_pairwise_bool = false;
    dinucCpG_pairwise_bool = false;
    nuc_site_comphet_bool = false;
    nuc1_site_comphet_bool = false;
    nuc2_site_comphet_bool = false;
    nuc3_site_comphet_bool = false;
    nuc_taxa_comphet_bool = false;
    nuc1_taxa_comphet_bool = false;
    nuc2_taxa_comphet_bool = false;
    nuc3_taxa_comphet_bool = false;
    codon_site_comphet_bool = false;
    codon_taxa_comphet_bool = false;
    aa_site_comphet_bool = false;
    aa_taxa_comphet_bool = false;

    CodonSequenceAlignment* simulation = new CodonSequenceAlignment(lparam->codondata,CurrentNodeLeafCodonSequence);


    string* arrSummaries = new string[lparam->NusedSummaries];
    for (unsigned int i_summary = 0 ; i_summary < lparam->NSummaries ; i_summary++){
        auto it = lparam->mapUsedSummaries.find(lparam->listSummaries[i_summary]);
        if(it != lparam->mapUsedSummaries.end() && it->second != -1) {
            arrSummaries[it->second] = it->first;
        }
    }


    for(unsigned int i_summary = 0 ; i_summary < lparam->NusedSummaries;i_summary++){
        auto it = GetSummariesMap.find(arrSummaries[i_summary]);
        if (it != GetSummariesMap.end()){
            funcpt f = GetSummariesMap[arrSummaries[i_summary]];
            double s = (this->*f)(simulation);
            //double s = log2(it.*second(simulation));
            if (s < lparam->TOOSMALL || isinf(s)) {
                s = lparam->TOOSMALL;
            }
            s = log2(s);
            lparam->summariesSimulatedData.push_back(s);
        }
    }
    delete simulation;
    delete [] arrSummaries;
}

void SummaryStatistics::computeSummaries(){

    lparam->summariesRealData.clear();
    lparam->summariesRealData.shrink_to_fit();

    codon_bool = false;
    dinuc_bool = false;
    dinuc12_bool = false;
    dinuc23_bool = false;
    dinuc31_bool = false;
    aa_bool = false;
    aa_wo_nr = false;
    dicodon_bool = false;
    diaa_bool = false;
    nuc_bool = false;
    nuc1_bool = false;
    nuc2_bool = false;
    nuc3_bool = false;
    nuc_meandiff_bool = false;
    nuc1_meandiff_bool = false;
    nuc2_meandiff_bool = false;
    nuc3_meandiff_bool = false;
    codon_meandiff_bool = false;
    aa_meandiff_bool = false;
    aa_wonR_meandiff_bool = false;
    CGNAGR_bool =false ;
    nuc_pairwise_bool = false;
    nuc1_pairwise_bool = false;
    nuc2_pairwise_bool = false;
    nuc3_pairwise_bool = false;
    aa_pairwise_bool = false;
    dinucCpG_pairwise_bool = false;
    nuc_site_comphet_bool = false;
    nuc1_site_comphet_bool = false;
    nuc2_site_comphet_bool = false;
    nuc3_site_comphet_bool = false;
    nuc_taxa_comphet_bool = false;
    nuc1_taxa_comphet_bool = false;
    nuc2_taxa_comphet_bool = false;
    nuc3_taxa_comphet_bool = false;
    codon_site_comphet_bool = false;
    codon_taxa_comphet_bool = false;
    aa_site_comphet_bool = false;
    aa_taxa_comphet_bool = false;

    string* arrSummaries = new string[lparam->NusedSummaries];
    for (unsigned int i_summary = 0 ; i_summary < lparam->NSummaries ; i_summary++){
        auto it = lparam->mapUsedSummaries.find(lparam->listSummaries[i_summary]);
        if(it != lparam->mapUsedSummaries.end() && it->second != -1) {
            arrSummaries[it->second] = it->first;
        }
    }


    for(unsigned int i_summary = 0 ; i_summary < lparam->NusedSummaries;i_summary++){
        auto it = GetSummariesMap.find(arrSummaries[i_summary]);
        if (it != GetSummariesMap.end()){
            funcpt f = GetSummariesMap[arrSummaries[i_summary]];

            double s = (this->*f)(lparam->codondata);
            //double s = log2(it.*second(simulation));
            if (s < lparam->TOOSMALL || isinf(s)) {
                s = lparam->TOOSMALL;
            }
            s = log2(s);
            lparam->summariesRealData.push_back(s);
        }
    }

    delete [] arrSummaries;
}



void SummaryStatistics::computeSummariesFromData(){
    ldata->summariesRealData.clear();
    ldata->summariesRealData.shrink_to_fit();

    codon_bool = false;
    dinuc_bool = false;
    dinuc12_bool = false;
    dinuc23_bool = false;
    dinuc31_bool = false;
    aa_bool = false;
    aa_wo_nr = false;
    dicodon_bool = false;
    diaa_bool = false;
    nuc_bool = false;
    nuc1_bool = false;
    nuc2_bool = false;
    nuc3_bool = false;
    nuc_meandiff_bool = false;
    nuc1_meandiff_bool = false;
    nuc2_meandiff_bool = false;
    nuc3_meandiff_bool = false;
    codon_meandiff_bool = false;
    aa_meandiff_bool = false;
    aa_wonR_meandiff_bool = false;
    CGNAGR_bool = false ;
    nuc_pairwise_bool = false;
    nuc1_pairwise_bool = false;
    nuc2_pairwise_bool = false;
    nuc3_pairwise_bool = false;
    aa_pairwise_bool = false;
    dinucCpG_pairwise_bool = false;
    nuc_site_comphet_bool = false;
    nuc1_site_comphet_bool = false;
    nuc2_site_comphet_bool = false;
    nuc3_site_comphet_bool = false;
    nuc_taxa_comphet_bool = false;
    nuc1_taxa_comphet_bool = false;
    nuc2_taxa_comphet_bool = false;
    nuc3_taxa_comphet_bool = false;
    codon_site_comphet_bool = false;
    codon_taxa_comphet_bool = false;
    aa_site_comphet_bool = false;
    aa_taxa_comphet_bool = false;

    cerr << "DFDSF0:q\n";

    string* arrSummaries = new string[ldata->NusedSummaries];
    for (unsigned int i_summary = 0 ; i_summary < ldata->NSummaries ; i_summary++){
        auto it = ldata->mapUsedSummaries.find(ldata->listSummaries[i_summary]);
        if(it != ldata->mapUsedSummaries.end() && it->second != -1) {
            arrSummaries[it->second] = it->first;
        }
    }

    cerr << "DFDSF1\n";

    for(unsigned int i_summary = 0 ; i_summary < ldata->NusedSummaries;i_summary++){
        auto it = GetSummariesMap.find(arrSummaries[i_summary]);
        if (it != GetSummariesMap.end()){
            funcpt f = GetSummariesMap[arrSummaries[i_summary]];

            double s = (this->*f)(ldata->codondata);
            //double s = log2(it.*second(simulation));
            if (s < ldata->TOOSMALL || isinf(s)) {
                s = ldata->TOOSMALL;
            }
            s = log2(s);
            ldata->summariesRealData.push_back(s);
        }
    }

    delete [] arrSummaries;
}







