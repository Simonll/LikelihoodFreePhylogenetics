#ifndef GLOBALPARAMETERS_H
#define GLOBALPARAMETERS_H

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
#include "BiologicalSequences.h"

class GlobalParameters
{
public:


    //const
    static constexpr double TOOSMALL = 1e-30;
    static constexpr double TOOLARGE = 500;
    static constexpr double TOOLARGENEGATIVE = -500;


    static const int NSummaries = 805;
    static const int NParam = 29;
    static const int NEvoStats = 1358;
    static const int NSiteSpecificEvoStats = 2;
    static const int NDistances = 3; 
    static const int NTransformations = 3; 
    static const int NCodes = 3; 

    const string listCodes[NCodes] = {
        "Universal","MtMam","MtInv"

    };  

    const string listTransformtations[NTransformations] = {
        "none","log2","log10"
    };

    const string listDistances[NDistances] = {
        "Euclidian","dist1","normalized"

    };

    const string listParam[NParam] = 
    {"chainID","root","lambda","lambda_CpG","lambda_TpA","lambdaTG","lambdaCA","lambda_TBL", "lambda_omega",
                                      "nucsA", "nucsC", "nucsG","nucsT",
                                      "nucrrAC","nucrrAG","nucrrAT","nucrrCA","nucrrCG","nucrrCT","nucrrGA","nucrrGC","nucrrGT","nucrrTA","nucrrTC","nucrrTG",
                                      "wR_CHQW","lambda_CpG_GpG","lambda_GpT","fitCpG"
    };



    const string listSummaries[NSummaries] =
    {
        "pwAC","pwAG","pwAT","pwCG","pwCT","pwGT",
        "pw1AC","pw1AG","pw1AT","pw1CG","pw1CT","pw1GT",
        "pw2AC","pw2AG","pw2AT","pw2CG","pw2CT","pw2GT",
        "pw3AC","pw3AG","pw3AT","pw3CG","pw3CT","pw3GT",
        "pwAC10","pwAG10","pwAT10","pwCG10","pwCT10","pwGT10","K80nuc10",
        "pwAC30","pwAG30","pwAT30","pwCG30","pwCT30","pwGT30","K80nuc30",
        "pwAC50","pwAG50","pwAT50","pwCG50","pwCT50","pwGT50","K80nuc50",
        "pwAA10","pwAA30","pwAA50","K80aa10",
        "pwAA","pwCpGTpG","pwCpGCpA","pwApGTpG",
        "nucA","nucC","nucG","nucT",
        "nuc1A","nuc1C","nuc1G","nuc1T",
        "nuc2A","nuc2C","nuc2G","nuc2T",
        "nuc3A","nuc3C","nuc3G","nuc3T",
        "dinucAA","dinucAC","dinucAG","dinucAT","dinucCA","dinucCC","dinucCG","dinucCT","dinucGA","dinucGC","dinucGG","dinucGT","dinucTA","dinucTC","dinucTG","dinucTT",
        "dinuc12AA","dinuc12AC","dinuc12AG","dinuc12AT","dinuc12CA","dinuc12CC","dinuc12CG","dinuc12CT","dinuc12GA","dinuc12GC","dinuc12GG","dinuc12GT","dinuc12TA","dinuc12TC","dinuc12TG","dinuc12TT",
        "dinuc23AA","dinuc23AC","dinuc23AG","dinuc23AT","dinuc23CA","dinuc23CC","dinuc23CG","dinuc23CT","dinuc23GA","dinuc23GC","dinuc23GG","dinuc23GT","dinuc23TA","dinuc23TC","dinuc23TG","dinuc23TT",
        "dinuc31AA","dinuc31AC","dinuc31AG","dinuc31AT","dinuc31CA","dinuc31CC","dinuc31CG","dinuc31CT","dinuc31GA","dinuc31GC","dinuc31GG","dinuc31GT","dinuc31TA","dinuc31TC","dinuc31TG","dinuc31TT",
        "TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT",
        "CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC",
        "GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG",
        "fTTT","fTTC","fTTA","fTTG","fTCT","fTCC","fTCA","fTCG","fTAT","fTAC","fTAA","fTAG","fTGT","fTGC","fTGA","fTGG","fCTT","fCTC","fCTA","fCTG","fCCT","fCCC","fCCA","fCCG","fCAT",
        "fCAC","fCAA","fCAG","fCGT","fCGC","fCGA","fCGG","fATT","fATC","fATA","fATG","fACT","fACC","fACA","fACG","fAAT","fAAC","fAAA","fAAG","fAGT","fAGC","fAGA","fAGG","fGTT","fGTC",
        "fGTA","fGTG","fGCT","fGCC","fGCA","fGCG","fGAT","fGAC","fGAA","fGAG","fGGT","fGGC","fGGA","fGGG",
        "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y",
        "nucsitecomphet","nuc1sitecomphet","nuc2sitecomphet","nuc3sitecomphet","nuctaxacomphet","nuc1taxacomphet","nuc2taxacomphet","nuc3taxacomphet","codonsitecomphet","codontaxacomphet","aasitecomphet","aataxacomphet",
        "mNuc","vNuc","mNuc1","vNuc1","mNuc2","vNuc2","mNuc3","vNuc3","mCodon","vCodon","mAA","vAA","mAAwonR","vAAwonR",
        "CGNAGRcgnconst","CGNAGRagrconst","CGNAGRcgnvar","CGNAGRagrvar",
        "AwonR","CwonR", "DwonR", "EwonR", "FwonR", "GwonR", "HwonR", "IwonR", "KwonR","LwonR", "MwonR", "NwonR", "PwonR", "QwonR", "RwonR", "SwonR", "TwonR", "VwonR", "WwonR", "YwonR",
        "TTTwonR","TTCwonR","TTAwonR","TTGwonR","TCTwonR","TCCwonR","TCAwonR","TCGwonR","TATwonR","TACwonR","TAAwonR","TAGwonR","TGTwonR","TGCwonR","TGAwonR","TGGwonR","CTTwonR","CTCwonR","CTAwonR","CTGwonR","CCTwonR","CCCwonR","CCAwonR","CCGwonR","CATwonR",
        "CACwonR","CAAwonR","CAGwonR","CGTwonR","CGCwonR","CGAwonR","CGGwonR","ATTwonR","ATCwonR","ATAwonR","ATGwonR","ACTwonR","ACCwonR","ACAwonR","ACGwonR","AATwonR","AACwonR","AAAwonR","AAGwonR","AGTwonR","AGCwonR","AGAwonR","AGGwonR","GTTwonR","GTCwonR",
        "GTAwonR","GTGwonR","GCTwonR","GCCwonR","GCAwonR","GCGwonR","GATwonR","GACwonR","GAAwonR","GAGwonR","GGTwonR","GGCwonR","GGAwonR","GGGwonR",
        "diaaAA","diaaAC","diaaAD","diaaAE","diaaAF","diaaAG","diaaAH","diaaAI","diaaAK","diaaAL","diaaAM","diaaAN","diaaAP","diaaAQ","diaaAR","diaaAS","diaaAT","diaaAV","diaaAW","diaaAY","diaaCA","diaaCC","diaaCD","diaaCE","diaaCF","diaaCG","diaaCH","diaaCI","diaaCK","diaaCL","diaaCM","diaaCN","diaaCP","diaaCQ","diaaCR","diaaCS","diaaCT","diaaCV","diaaCW","diaaCY","diaaDA","diaaDC","diaaDD","diaaDE","diaaDF","diaaDG","diaaDH","diaaDI","diaaDK","diaaDL","diaaDM","diaaDN","diaaDP","diaaDQ","diaaDR","diaaDS","diaaDT","diaaDV","diaaDW","diaaDY","diaaEA","diaaEC","diaaED","diaaEE","diaaEF","diaaEG","diaaEH","diaaEI","diaaEK","diaaEL","diaaEM","diaaEN","diaaEP","diaaEQ","diaaER","diaaES","diaaET","diaaEV","diaaEW","diaaEY","diaaFA","diaaFC","diaaFD","diaaFE","diaaFF","diaaFG","diaaFH","diaaFI","diaaFK","diaaFL","diaaFM","diaaFN","diaaFP","diaaFQ","diaaFR","diaaFS","diaaFT","diaaFV","diaaFW","diaaFY","diaaGA","diaaGC","diaaGD","diaaGE","diaaGF","diaaGG","diaaGH","diaaGI","diaaGK","diaaGL","diaaGM","diaaGN","diaaGP","diaaGQ","diaaGR","diaaGS","diaaGT","diaaGV","diaaGW","diaaGY","diaaHA","diaaHC","diaaHD","diaaHE","diaaHF","diaaHG","diaaHH","diaaHI","diaaHK","diaaHL","diaaHM","diaaHN","diaaHP","diaaHQ","diaaHR","diaaHS","diaaHT","diaaHV","diaaHW","diaaHY","diaaIA","diaaIC","diaaID","diaaIE","diaaIF","diaaIG","diaaIH","diaaII","diaaIK","diaaIL","diaaIM","diaaIN","diaaIP","diaaIQ","diaaIR","diaaIS","diaaIT","diaaIV","diaaIW","diaaIY","diaaKA","diaaKC","diaaKD","diaaKE","diaaKF","diaaKG","diaaKH","diaaKI","diaaKK","diaaKL","diaaKM",
        "diaaKN","diaaKP","diaaKQ","diaaKR","diaaKS","diaaKT","diaaKV","diaaKW","diaaKY","diaaLA","diaaLC","diaaLD","diaaLE","diaaLF","diaaLG","diaaLH","diaaLI","diaaLK","diaaLL","diaaLM","diaaLN","diaaLP","diaaLQ","diaaLR","diaaLS","diaaLT","diaaLV","diaaLW","diaaLY","diaaMA","diaaMC","diaaMD","diaaME","diaaMF","diaaMG","diaaMH","diaaMI","diaaMK","diaaML","diaaMM","diaaMN","diaaMP","diaaMQ","diaaMR","diaaMS","diaaMT","diaaMV","diaaMW","diaaMY","diaaNA","diaaNC","diaaND","diaaNE","diaaNF","diaaNG","diaaNH","diaaNI","diaaNK","diaaNL","diaaNM","diaaNN","diaaNP","diaaNQ","diaaNR","diaaNS","diaaNT","diaaNV","diaaNW","diaaNY","diaaPA","diaaPC","diaaPD","diaaPE","diaaPF","diaaPG","diaaPH","diaaPI","diaaPK","diaaPL","diaaPM","diaaPN","diaaPP","diaaPQ","diaaPR","diaaPS","diaaPT","diaaPV","diaaPW","diaaPY","diaaQA","diaaQC","diaaQD","diaaQE","diaaQF","diaaQG","diaaQH","diaaQI","diaaQK","diaaQL","diaaQM","diaaQN","diaaQP","diaaQQ","diaaQR","diaaQS","diaaQT","diaaQV","diaaQW","diaaQY","diaaRA","diaaRC","diaaRD","diaaRE","diaaRF","diaaRG","diaaRH","diaaRI","diaaRK","diaaRL","diaaRM","diaaRN","diaaRP","diaaRQ","diaaRR","diaaRS","diaaRT","diaaRV","diaaRW","diaaRY","diaaSA","diaaSC","diaaSD","diaaSE","diaaSF","diaaSG","diaaSH","diaaSI","diaaSK","diaaSL","diaaSM","diaaSN","diaaSP","diaaSQ","diaaSR","diaaSS","diaaST","diaaSV","diaaSW","diaaSY","diaaTA","diaaTC","diaaTD","diaaTE","diaaTF","diaaTG","diaaTH","diaaTI","diaaTK","diaaTL","diaaTM","diaaTN","diaaTP","diaaTQ","diaaTR","diaaTS","diaaTT","diaaTV","diaaTW","diaaTY","diaaVA","diaaVC",
        "diaaVD","diaaVE","diaaVF","diaaVG","diaaVH","diaaVI","diaaVK","diaaVL","diaaVM","diaaVN","diaaVP","diaaVQ","diaaVR","diaaVS","diaaVT","diaaVV","diaaVW","diaaVY","diaaWA","diaaWC","diaaWD","diaaWE","diaaWF","diaaWG","diaaWH","diaaWI","diaaWK","diaaWL","diaaWM","diaaWN","diaaWP","diaaWQ","diaaWR","diaaWS","diaaWT","diaaWV","diaaWW","diaaWY","diaaYA","diaaYC","diaaYD","diaaYE","diaaYF","diaaYG","diaaYH","diaaYI","diaaYK","diaaYL","diaaYM","diaaYN","diaaYP","diaaYQ","diaaYR","diaaYS","diaaYT","diaaYV","diaaYW","diaaYY",
        "RSCUentropy","GC","GC1","GC2","GC3","Codonfentropy","AAentropy","CGNonAGR","CHQWonR","LMVonAPST"

    };

    const string listEvoStats[NEvoStats] =
    {
        "gtnrAA","gtnrAC","gtnrAG","gtnrAT","gtnrCA","gtnrCC","gtnrCG","gtnrCT","gtnrGA","gtnrGC","gtnrGG","gtnrGT","gtnrTA","gtnrTC","gtnrTG","gtnrTT",
        "gtnrSAA","gtnrSAC","gtnrSAG","gtnrSAT","gtnrSCA","gtnrSCC","gtnrSCG","gtnrSCT","gtnrSGA","gtnrSGC","gtnrSGG","gtnrSGT","gtnrSTA","gtnrSTC","gtnrSTG","gtnrSTT",
        "gtnrSNAA","gtnrSNAC","gtnrSNAG","gtnrSNAT","gtnrSNCA","gtnrSNCC","gtnrSNCG","gtnrSNCT","gtnrSNGA","gtnrSNGC","gtnrSNGG","gtnrSNGT","gtnrSNTA","gtnrSNTC","gtnrSNTG","gtnrSNTT",

        "gtnr1AA","gtnr1AC","gtnr1AG","gtnr1AT","gtnr1CA","gtnr1CC","gtnr1CG","gtnr1CT","gtnr1GA","gtnr1GC","gtnr1GG","gtnr1GT","gtnr1TA","gtnr1TC","gtnr1TG","gtnr1TT",
        "gtnr2AA","gtnr2AC","gtnr2AG","gtnr2AT","gtnr2CA","gtnr2CC","gtnr2CG","gtnr2CT","gtnr2GA","gtnr2GC","gtnr2GG","gtnr2GT","gtnr2TA","gtnr2TC","gtnr2TG","gtnr2TT",
        "gtnr3AA","gtnr3AC","gtnr3AG","gtnr3AT","gtnr3CA","gtnr3CC","gtnr3CG","gtnr3CT","gtnr3GA","gtnr3GC","gtnr3GG","gtnr3GT","gtnr3TA","gtnr3TC","gtnr3TG","gtnr3TT",

        "gtnr1SAA","gtnr1SAC","gtnr1SAG","gtnr1SAT","gtnr1SCA","gtnr1SCC","gtnr1SCG","gtnr1SCT","gtnr1SGA","gtnr1SGC","gtnr1SGG","gtnr1SGT","gtnr1STA","gtnr1STC","gtnr1STG","gtnr1STT",
        "gtnr2SAA","gtnr2SAC","gtnr2SAG","gtnr2SAT","gtnr2SCA","gtnr2SCC","gtnr2SCG","gtnr2SCT","gtnr2SGA","gtnr2SGC","gtnr2SGG","gtnr2SGT","gtnr2STA","gtnr2STC","gtnr2STG","gtnr2STT",
        "gtnr3SAA","gtnr3SAC","gtnr3SAG","gtnr3SAT","gtnr3SCA","gtnr3SCC","gtnr3SCG","gtnr3SCT","gtnr3SGA","gtnr3SGC","gtnr3SGG","gtnr3SGT","gtnr3STA","gtnr3STC","gtnr3STG","gtnr3STT",

        "gtnr1NSAA","gtnr1NSAC","gtnr1NSAG","gtnr1NSAT","gtnr1NSCA","gtnr1NSCC","gtnr1NSCG","gtnr1NSCT","gtnr1NSGA","gtnr1NSGC","gtnr1NSGG","gtnr1NSGT","gtnr1NSTA","gtnr1NSTC","gtnr1NSTG","gtnr1NSTT",
        "gtnr2NSAA","gtnr2NSAC","gtnr2NSAG","gtnr2NSAT","gtnr2NSCA","gtnr2NSCC","gtnr2NSCG","gtnr2NSCT","gtnr2NSGA","gtnr2NSGC","gtnr2NSGG","gtnr2NSGT","gtnr2NSTA","gtnr2NSTC","gtnr2NSTG","gtnr2NSTT",
        "gtnr3NSAA","gtnr3NSAC","gtnr3NSAG","gtnr3NSAT","gtnr3NSCA","gtnr3NSCC","gtnr3NSCG","gtnr3NSCT","gtnr3NSGA","gtnr3NSGC","gtnr3NSGG","gtnr3NSGT","gtnr3NSTA","gtnr3NSTC","gtnr3NSTG","gtnr3NSTT",

        "dinucAAAC","dinucAAAG","dinucAAAT","dinucAACA","dinucAAGA","dinucAATA","dinucACAA","dinucACAG","dinucACAT","dinucACCC","dinucACGC","dinucACTC","dinucAGAA","dinucAGAC","dinucAGAT","dinucAGCG","dinucAGGG","dinucAGTG","dinucATAA","dinucATAC","dinucATAG","dinucATCT","dinucATGT","dinucATTT","dinucCAAA","dinucCACC","dinucCACG","dinucCACT","dinucCAGA","dinucCATA","dinucCCAC","dinucCCCA","dinucCCCG","dinucCCCT","dinucCCGC","dinucCCTC","dinucCGAG","dinucCGCA","dinucCGCC","dinucCGCT","dinucCGGG","dinucCGTG","dinucCTAT","dinucCTCA","dinucCTCC","dinucCTCG","dinucCTGT","dinucCTTT","dinucGAAA","dinucGACA","dinucGAGC","dinucGAGG","dinucGAGT","dinucGATA","dinucGCAC","dinucGCCC","dinucGCGA","dinucGCGG","dinucGCGT","dinucGCTC","dinucGGAG","dinucGGCG","dinucGGGA","dinucGGGC","dinucGGGT","dinucGGTG","dinucGTAT","dinucGTCT","dinucGTGA","dinucGTGC","dinucGTGG","dinucGTTT","dinucTAAA","dinucTACA","dinucTAGA","dinucTATC","dinucTATG","dinucTATT","dinucTCAC","dinucTCCC","dinucTCGC","dinucTCTA","dinucTCTG","dinucTCTT","dinucTGAG","dinucTGCG","dinucTGGG","dinucTGTA","dinucTGTC","dinucTGTT","dinucTTAT","dinucTTCT","dinucTTGT","dinucTTTA","dinucTTTC","dinucTTTG",
        "dinuc12AAAC","dinuc12AAAG","dinuc12AAAT","dinuc12AACA","dinuc12AAGA","dinuc12AATA","dinuc12ACAA","dinuc12ACAG","dinuc12ACAT","dinuc12ACCC","dinuc12ACGC","dinuc12ACTC","dinuc12AGAA","dinuc12AGAC","dinuc12AGAT","dinuc12AGCG","dinuc12AGGG","dinuc12AGTG","dinuc12ATAA","dinuc12ATAC","dinuc12ATAG","dinuc12ATCT","dinuc12ATGT","dinuc12ATTT","dinuc12CAAA","dinuc12CACC","dinuc12CACG","dinuc12CACT","dinuc12CAGA","dinuc12CATA","dinuc12CCAC","dinuc12CCCA","dinuc12CCCG","dinuc12CCCT","dinuc12CCGC","dinuc12CCTC","dinuc12CGAG","dinuc12CGCA","dinuc12CGCC","dinuc12CGCT","dinuc12CGGG","dinuc12CGTG","dinuc12CTAT","dinuc12CTCA","dinuc12CTCC","dinuc12CTCG","dinuc12CTGT","dinuc12CTTT","dinuc12GAAA","dinuc12GACA","dinuc12GAGC","dinuc12GAGG","dinuc12GAGT","dinuc12GATA","dinuc12GCAC","dinuc12GCCC","dinuc12GCGA","dinuc12GCGG","dinuc12GCGT","dinuc12GCTC","dinuc12GGAG","dinuc12GGCG","dinuc12GGGA","dinuc12GGGC","dinuc12GGGT","dinuc12GGTG","dinuc12GTAT","dinuc12GTCT","dinuc12GTGA","dinuc12GTGC","dinuc12GTGG","dinuc12GTTT","dinuc12TAAA","dinuc12TACA","dinuc12TAGA","dinuc12TATC","dinuc12TATG","dinuc12TATT","dinuc12TCAC","dinuc12TCCC","dinuc12TCGC","dinuc12TCTA","dinuc12TCTG","dinuc12TCTT","dinuc12TGAG","dinuc12TGCG","dinuc12TGGG","dinuc12TGTA","dinuc12TGTC","dinuc12TGTT","dinuc12TTAT","dinuc12TTCT","dinuc12TTGT","dinuc12TTTA","dinuc12TTTC","dinuc12TTTG",
        "dinuc23AAAC","dinuc23AAAG","dinuc23AAAT","dinuc23AACA","dinuc23AAGA","dinuc23AATA","dinuc23ACAA","dinuc23ACAG","dinuc23ACAT","dinuc23ACCC","dinuc23ACGC","dinuc23ACTC","dinuc23AGAA","dinuc23AGAC","dinuc23AGAT","dinuc23AGCG","dinuc23AGGG","dinuc23AGTG","dinuc23ATAA","dinuc23ATAC","dinuc23ATAG","dinuc23ATCT","dinuc23ATGT","dinuc23ATTT","dinuc23CAAA","dinuc23CACC","dinuc23CACG","dinuc23CACT","dinuc23CAGA","dinuc23CATA","dinuc23CCAC","dinuc23CCCA","dinuc23CCCG","dinuc23CCCT","dinuc23CCGC","dinuc23CCTC","dinuc23CGAG","dinuc23CGCA","dinuc23CGCC","dinuc23CGCT","dinuc23CGGG","dinuc23CGTG","dinuc23CTAT","dinuc23CTCA","dinuc23CTCC","dinuc23CTCG","dinuc23CTGT","dinuc23CTTT","dinuc23GAAA","dinuc23GACA","dinuc23GAGC","dinuc23GAGG","dinuc23GAGT","dinuc23GATA","dinuc23GCAC","dinuc23GCCC","dinuc23GCGA","dinuc23GCGG","dinuc23GCGT","dinuc23GCTC","dinuc23GGAG","dinuc23GGCG","dinuc23GGGA","dinuc23GGGC","dinuc23GGGT","dinuc23GGTG","dinuc23GTAT","dinuc23GTCT","dinuc23GTGA","dinuc23GTGC","dinuc23GTGG","dinuc23GTTT","dinuc23TAAA","dinuc23TACA","dinuc23TAGA","dinuc23TATC","dinuc23TATG","dinuc23TATT","dinuc23TCAC","dinuc23TCCC","dinuc23TCGC","dinuc23TCTA","dinuc23TCTG","dinuc23TCTT","dinuc23TGAG","dinuc23TGCG","dinuc23TGGG","dinuc23TGTA","dinuc23TGTC","dinuc23TGTT","dinuc23TTAT","dinuc23TTCT","dinuc23TTGT","dinuc23TTTA","dinuc23TTTC","dinuc23TTTG",
        "dinuc31AAAC","dinuc31AAAG","dinuc31AAAT","dinuc31AACA","dinuc31AAGA","dinuc31AATA","dinuc31ACAA","dinuc31ACAG","dinuc31ACAT","dinuc31ACCC","dinuc31ACGC","dinuc31ACTC","dinuc31AGAA","dinuc31AGAC","dinuc31AGAT","dinuc31AGCG","dinuc31AGGG","dinuc31AGTG","dinuc31ATAA","dinuc31ATAC","dinuc31ATAG","dinuc31ATCT","dinuc31ATGT","dinuc31ATTT","dinuc31CAAA","dinuc31CACC","dinuc31CACG","dinuc31CACT","dinuc31CAGA","dinuc31CATA","dinuc31CCAC","dinuc31CCCA","dinuc31CCCG","dinuc31CCCT","dinuc31CCGC","dinuc31CCTC","dinuc31CGAG","dinuc31CGCA","dinuc31CGCC","dinuc31CGCT","dinuc31CGGG","dinuc31CGTG","dinuc31CTAT","dinuc31CTCA","dinuc31CTCC","dinuc31CTCG","dinuc31CTGT","dinuc31CTTT","dinuc31GAAA","dinuc31GACA","dinuc31GAGC","dinuc31GAGG","dinuc31GAGT","dinuc31GATA","dinuc31GCAC","dinuc31GCCC","dinuc31GCGA","dinuc31GCGG","dinuc31GCGT","dinuc31GCTC","dinuc31GGAG","dinuc31GGCG","dinuc31GGGA","dinuc31GGGC","dinuc31GGGT","dinuc31GGTG","dinuc31GTAT","dinuc31GTCT","dinuc31GTGA","dinuc31GTGC","dinuc31GTGG","dinuc31GTTT","dinuc31TAAA","dinuc31TACA","dinuc31TAGA","dinuc31TATC","dinuc31TATG","dinuc31TATT","dinuc31TCAC","dinuc31TCCC","dinuc31TCGC","dinuc31TCTA","dinuc31TCTG","dinuc31TCTT","dinuc31TGAG","dinuc31TGCG","dinuc31TGGG","dinuc31TGTA","dinuc31TGTC","dinuc31TGTT","dinuc31TTAT","dinuc31TTCT","dinuc31TTGT","dinuc31TTTA","dinuc31TTTC","dinuc31TTTG",

        "dinucSAAAC","dinucSAAAG","dinucSAAAT","dinucSAACA","dinucSAAGA","dinucSAATA","dinucSACAA","dinucSACAG","dinucSACAT","dinucSACCC","dinucSACGC","dinucSACTC","dinucSAGAA","dinucSAGAC","dinucSAGAT","dinucSAGCG","dinucSAGGG","dinucSAGTG","dinucSATAA","dinucSATAC","dinucSATAG","dinucSATCT","dinucSATGT","dinucSATTT","dinucSCAAA","dinucSCACC","dinucSCACG","dinucSCACT","dinucSCAGA","dinucSCATA","dinucSCCAC","dinucSCCCA","dinucSCCCG","dinucSCCCT","dinucSCCGC","dinucSCCTC","dinucSCGAG","dinucSCGCA","dinucSCGCC","dinucSCGCT","dinucSCGGG","dinucSCGTG","dinucSCTAT","dinucSCTCA","dinucSCTCC","dinucSCTCG","dinucSCTGT","dinucSCTTT","dinucSGAAA","dinucSGACA","dinucSGAGC","dinucSGAGG","dinucSGAGT","dinucSGATA","dinucSGCAC","dinucSGCCC","dinucSGCGA","dinucSGCGG","dinucSGCGT","dinucSGCTC","dinucSGGAG","dinucSGGCG","dinucSGGGA","dinucSGGGC","dinucSGGGT","dinucSGGTG","dinucSGTAT","dinucSGTCT","dinucSGTGA","dinucSGTGC","dinucSGTGG","dinucSGTTT","dinucSTAAA","dinucSTACA","dinucSTAGA","dinucSTATC","dinucSTATG","dinucSTATT","dinucSTCAC","dinucSTCCC","dinucSTCGC","dinucSTCTA","dinucSTCTG","dinucSTCTT","dinucSTGAG","dinucSTGCG","dinucSTGGG","dinucSTGTA","dinucSTGTC","dinucSTGTT","dinucSTTAT","dinucSTTCT","dinucSTTGT","dinucSTTTA","dinucSTTTC","dinucSTTTG",
        "dinuc12SAAAC","dinuc12SAAAG","dinuc12SAAAT","dinuc12SAACA","dinuc12SAAGA","dinuc12SAATA","dinuc12SACAA","dinuc12SACAG","dinuc12SACAT","dinuc12SACCC","dinuc12SACGC","dinuc12SACTC","dinuc12SAGAA","dinuc12SAGAC","dinuc12SAGAT","dinuc12SAGCG","dinuc12SAGGG","dinuc12SAGTG","dinuc12SATAA","dinuc12SATAC","dinuc12SATAG","dinuc12SATCT","dinuc12SATGT","dinuc12SATTT","dinuc12SCAAA","dinuc12SCACC","dinuc12SCACG","dinuc12SCACT","dinuc12SCAGA","dinuc12SCATA","dinuc12SCCAC","dinuc12SCCCA","dinuc12SCCCG","dinuc12SCCCT","dinuc12SCCGC","dinuc12SCCTC","dinuc12SCGAG","dinuc12SCGCA","dinuc12SCGCC","dinuc12SCGCT","dinuc12SCGGG","dinuc12SCGTG","dinuc12SCTAT","dinuc12SCTCA","dinuc12SCTCC","dinuc12SCTCG","dinuc12SCTGT","dinuc12SCTTT","dinuc12SGAAA","dinuc12SGACA","dinuc12SGAGC","dinuc12SGAGG","dinuc12SGAGT","dinuc12SGATA","dinuc12SGCAC","dinuc12SGCCC","dinuc12SGCGA","dinuc12SGCGG","dinuc12SGCGT","dinuc12SGCTC","dinuc12SGGAG","dinuc12SGGCG","dinuc12SGGGA","dinuc12SGGGC","dinuc12SGGGT","dinuc12SGGTG","dinuc12SGTAT","dinuc12SGTCT","dinuc12SGTGA","dinuc12SGTGC","dinuc12SGTGG","dinuc12SGTTT","dinuc12STAAA","dinuc12STACA","dinuc12STAGA","dinuc12STATC","dinuc12STATG","dinuc12STATT","dinuc12STCAC","dinuc12STCCC","dinuc12STCGC","dinuc12STCTA","dinuc12STCTG","dinuc12STCTT","dinuc12STGAG","dinuc12STGCG","dinuc12STGGG","dinuc12STGTA","dinuc12STGTC","dinuc12STGTT","dinuc12STTAT","dinuc12STTCT","dinuc12STTGT","dinuc12STTTA","dinuc12STTTC","dinuc12STTTG",
        "dinuc23SAAAC","dinuc23SAAAG","dinuc23SAAAT","dinuc23SAACA","dinuc23SAAGA","dinuc23SAATA","dinuc23SACAA","dinuc23SACAG","dinuc23SACAT","dinuc23SACCC","dinuc23SACGC","dinuc23SACTC","dinuc23SAGAA","dinuc23SAGAC","dinuc23SAGAT","dinuc23SAGCG","dinuc23SAGGG","dinuc23SAGTG","dinuc23SATAA","dinuc23SATAC","dinuc23SATAG","dinuc23SATCT","dinuc23SATGT","dinuc23SATTT","dinuc23SCAAA","dinuc23SCACC","dinuc23SCACG","dinuc23SCACT","dinuc23SCAGA","dinuc23SCATA","dinuc23SCCAC","dinuc23SCCCA","dinuc23SCCCG","dinuc23SCCCT","dinuc23SCCGC","dinuc23SCCTC","dinuc23SCGAG","dinuc23SCGCA","dinuc23SCGCC","dinuc23SCGCT","dinuc23SCGGG","dinuc23SCGTG","dinuc23SCTAT","dinuc23SCTCA","dinuc23SCTCC","dinuc23SCTCG","dinuc23SCTGT","dinuc23SCTTT","dinuc23SGAAA","dinuc23SGACA","dinuc23SGAGC","dinuc23SGAGG","dinuc23SGAGT","dinuc23SGATA","dinuc23SGCAC","dinuc23SGCCC","dinuc23SGCGA","dinuc23SGCGG","dinuc23SGCGT","dinuc23SGCTC","dinuc23SGGAG","dinuc23SGGCG","dinuc23SGGGA","dinuc23SGGGC","dinuc23SGGGT","dinuc23SGGTG","dinuc23SGTAT","dinuc23SGTCT","dinuc23SGTGA","dinuc23SGTGC","dinuc23SGTGG","dinuc23SGTTT","dinuc23STAAA","dinuc23STACA","dinuc23STAGA","dinuc23STATC","dinuc23STATG","dinuc23STATT","dinuc23STCAC","dinuc23STCCC","dinuc23STCGC","dinuc23STCTA","dinuc23STCTG","dinuc23STCTT","dinuc23STGAG","dinuc23STGCG","dinuc23STGGG","dinuc23STGTA","dinuc23STGTC","dinuc23STGTT","dinuc23STTAT","dinuc23STTCT","dinuc23STTGT","dinuc23STTTA","dinuc23STTTC","dinuc23STTTG",
        "dinuc31SAAAC","dinuc31SAAAG","dinuc31SAAAT","dinuc31SAACA","dinuc31SAAGA","dinuc31SAATA","dinuc31SACAA","dinuc31SACAG","dinuc31SACAT","dinuc31SACCC","dinuc31SACGC","dinuc31SACTC","dinuc31SAGAA","dinuc31SAGAC","dinuc31SAGAT","dinuc31SAGCG","dinuc31SAGGG","dinuc31SAGTG","dinuc31SATAA","dinuc31SATAC","dinuc31SATAG","dinuc31SATCT","dinuc31SATGT","dinuc31SATTT","dinuc31SCAAA","dinuc31SCACC","dinuc31SCACG","dinuc31SCACT","dinuc31SCAGA","dinuc31SCATA","dinuc31SCCAC","dinuc31SCCCA","dinuc31SCCCG","dinuc31SCCCT","dinuc31SCCGC","dinuc31SCCTC","dinuc31SCGAG","dinuc31SCGCA","dinuc31SCGCC","dinuc31SCGCT","dinuc31SCGGG","dinuc31SCGTG","dinuc31SCTAT","dinuc31SCTCA","dinuc31SCTCC","dinuc31SCTCG","dinuc31SCTGT","dinuc31SCTTT","dinuc31SGAAA","dinuc31SGACA","dinuc31SGAGC","dinuc31SGAGG","dinuc31SGAGT","dinuc31SGATA","dinuc31SGCAC","dinuc31SGCCC","dinuc31SGCGA","dinuc31SGCGG","dinuc31SGCGT","dinuc31SGCTC","dinuc31SGGAG","dinuc31SGGCG","dinuc31SGGGA","dinuc31SGGGC","dinuc31SGGGT","dinuc31SGGTG","dinuc31SGTAT","dinuc31SGTCT","dinuc31SGTGA","dinuc31SGTGC","dinuc31SGTGG","dinuc31SGTTT","dinuc31STAAA","dinuc31STACA","dinuc31STAGA","dinuc31STATC","dinuc31STATG","dinuc31STATT","dinuc31STCAC","dinuc31STCCC","dinuc31STCGC","dinuc31STCTA","dinuc31STCTG","dinuc31STCTT","dinuc31STGAG","dinuc31STGCG","dinuc31STGGG","dinuc31STGTA","dinuc31STGTC","dinuc31STGTT","dinuc31STTAT","dinuc31STTCT","dinuc31STTGT","dinuc31STTTA","dinuc31STTTC","dinuc31STTTG",

        "dinucNSAAAC","dinucNSAAAG","dinucNSAAAT","dinucNSAACA","dinucNSAAGA","dinucNSAATA","dinucNSACAA","dinucNSACAG","dinucNSACAT","dinucNSACCC","dinucNSACGC","dinucNSACTC","dinucNSAGAA","dinucNSAGAC","dinucNSAGAT","dinucNSAGCG","dinucNSAGGG","dinucNSAGTG","dinucNSATAA","dinucNSATAC","dinucNSATAG","dinucNSATCT","dinucNSATGT","dinucNSATTT","dinucNSCAAA","dinucNSCACC","dinucNSCACG","dinucNSCACT","dinucNSCAGA","dinucNSCATA","dinucNSCCAC","dinucNSCCCA","dinucNSCCCG","dinucNSCCCT","dinucNSCCGC","dinucNSCCTC","dinucNSCGAG","dinucNSCGCA","dinucNSCGCC","dinucNSCGCT","dinucNSCGGG","dinucNSCGTG","dinucNSCTAT","dinucNSCTCA","dinucNSCTCC","dinucNSCTCG","dinucNSCTGT","dinucNSCTTT","dinucNSGAAA","dinucNSGACA","dinucNSGAGC","dinucNSGAGG","dinucNSGAGT","dinucNSGATA","dinucNSGCAC","dinucNSGCCC","dinucNSGCGA","dinucNSGCGG","dinucNSGCGT","dinucNSGCTC","dinucNSGGAG","dinucNSGGCG","dinucNSGGGA","dinucNSGGGC","dinucNSGGGT","dinucNSGGTG","dinucNSGTAT","dinucNSGTCT","dinucNSGTGA","dinucNSGTGC","dinucNSGTGG","dinucNSGTTT","dinucNSTAAA","dinucNSTACA","dinucNSTAGA","dinucNSTATC","dinucNSTATG","dinucNSTATT","dinucNSTCAC","dinucNSTCCC","dinucNSTCGC","dinucNSTCTA","dinucNSTCTG","dinucNSTCTT","dinucNSTGAG","dinucNSTGCG","dinucNSTGGG","dinucNSTGTA","dinucNSTGTC","dinucNSTGTT","dinucNSTTAT","dinucNSTTCT","dinucNSTTGT","dinucNSTTTA","dinucNSTTTC","dinucNSTTTG",
        "dinuc12NSAAAC","dinuc12NSAAAG","dinuc12NSAAAT","dinuc12NSAACA","dinuc12NSAAGA","dinuc12NSAATA","dinuc12NSACAA","dinuc12NSACAG","dinuc12NSACAT","dinuc12NSACCC","dinuc12NSACGC","dinuc12NSACTC","dinuc12NSAGAA","dinuc12NSAGAC","dinuc12NSAGAT","dinuc12NSAGCG","dinuc12NSAGGG","dinuc12NSAGTG","dinuc12NSATAA","dinuc12NSATAC","dinuc12NSATAG","dinuc12NSATCT","dinuc12NSATGT","dinuc12NSATTT","dinuc12NSCAAA","dinuc12NSCACC","dinuc12NSCACG","dinuc12NSCACT","dinuc12NSCAGA","dinuc12NSCATA","dinuc12NSCCAC","dinuc12NSCCCA","dinuc12NSCCCG","dinuc12NSCCCT","dinuc12NSCCGC","dinuc12NSCCTC","dinuc12NSCGAG","dinuc12NSCGCA","dinuc12NSCGCC","dinuc12NSCGCT","dinuc12NSCGGG","dinuc12NSCGTG","dinuc12NSCTAT","dinuc12NSCTCA","dinuc12NSCTCC","dinuc12NSCTCG","dinuc12NSCTGT","dinuc12NSCTTT","dinuc12NSGAAA","dinuc12NSGACA","dinuc12NSGAGC","dinuc12NSGAGG","dinuc12NSGAGT","dinuc12NSGATA","dinuc12NSGCAC","dinuc12NSGCCC","dinuc12NSGCGA","dinuc12NSGCGG","dinuc12NSGCGT","dinuc12NSGCTC","dinuc12NSGGAG","dinuc12NSGGCG","dinuc12NSGGGA","dinuc12NSGGGC","dinuc12NSGGGT","dinuc12NSGGTG","dinuc12NSGTAT","dinuc12NSGTCT","dinuc12NSGTGA","dinuc12NSGTGC","dinuc12NSGTGG","dinuc12NSGTTT","dinuc12NSTAAA","dinuc12NSTACA","dinuc12NSTAGA","dinuc12NSTATC","dinuc12NSTATG","dinuc12NSTATT","dinuc12NSTCAC","dinuc12NSTCCC","dinuc12NSTCGC","dinuc12NSTCTA","dinuc12NSTCTG","dinuc12NSTCTT","dinuc12NSTGAG","dinuc12NSTGCG","dinuc12NSTGGG","dinuc12NSTGTA","dinuc12NSTGTC","dinuc12NSTGTT","dinuc12NSTTAT","dinuc12NSTTCT","dinuc12NSTTGT","dinuc12NSTTTA","dinuc12NSTTTC","dinuc12NSTTTG",
        "dinuc23NSAAAC","dinuc23NSAAAG","dinuc23NSAAAT","dinuc23NSAACA","dinuc23NSAAGA","dinuc23NSAATA","dinuc23NSACAA","dinuc23NSACAG","dinuc23NSACAT","dinuc23NSACCC","dinuc23NSACGC","dinuc23NSACTC","dinuc23NSAGAA","dinuc23NSAGAC","dinuc23NSAGAT","dinuc23NSAGCG","dinuc23NSAGGG","dinuc23NSAGTG","dinuc23NSATAA","dinuc23NSATAC","dinuc23NSATAG","dinuc23NSATCT","dinuc23NSATGT","dinuc23NSATTT","dinuc23NSCAAA","dinuc23NSCACC","dinuc23NSCACG","dinuc23NSCACT","dinuc23NSCAGA","dinuc23NSCATA","dinuc23NSCCAC","dinuc23NSCCCA","dinuc23NSCCCG","dinuc23NSCCCT","dinuc23NSCCGC","dinuc23NSCCTC","dinuc23NSCGAG","dinuc23NSCGCA","dinuc23NSCGCC","dinuc23NSCGCT","dinuc23NSCGGG","dinuc23NSCGTG","dinuc23NSCTAT","dinuc23NSCTCA","dinuc23NSCTCC","dinuc23NSCTCG","dinuc23NSCTGT","dinuc23NSCTTT","dinuc23NSGAAA","dinuc23NSGACA","dinuc23NSGAGC","dinuc23NSGAGG","dinuc23NSGAGT","dinuc23NSGATA","dinuc23NSGCAC","dinuc23NSGCCC","dinuc23NSGCGA","dinuc23NSGCGG","dinuc23NSGCGT","dinuc23NSGCTC","dinuc23NSGGAG","dinuc23NSGGCG","dinuc23NSGGGA","dinuc23NSGGGC","dinuc23NSGGGT","dinuc23NSGGTG","dinuc23NSGTAT","dinuc23NSGTCT","dinuc23NSGTGA","dinuc23NSGTGC","dinuc23NSGTGG","dinuc23NSGTTT","dinuc23NSTAAA","dinuc23NSTACA","dinuc23NSTAGA","dinuc23NSTATC","dinuc23NSTATG","dinuc23NSTATT","dinuc23NSTCAC","dinuc23NSTCCC","dinuc23NSTCGC","dinuc23NSTCTA","dinuc23NSTCTG","dinuc23NSTCTT","dinuc23NSTGAG","dinuc23NSTGCG","dinuc23NSTGGG","dinuc23NSTGTA","dinuc23NSTGTC","dinuc23NSTGTT","dinuc23NSTTAT","dinuc23NSTTCT","dinuc23NSTTGT","dinuc23NSTTTA","dinuc23NSTTTC","dinuc23NSTTTG",
        "dinuc31NSAAAC","dinuc31NSAAAG","dinuc31NSAAAT","dinuc31NSAACA","dinuc31NSAAGA","dinuc31NSAATA","dinuc31NSACAA","dinuc31NSACAG","dinuc31NSACAT","dinuc31NSACCC","dinuc31NSACGC","dinuc31NSACTC","dinuc31NSAGAA","dinuc31NSAGAC","dinuc31NSAGAT","dinuc31NSAGCG","dinuc31NSAGGG","dinuc31NSAGTG","dinuc31NSATAA","dinuc31NSATAC","dinuc31NSATAG","dinuc31NSATCT","dinuc31NSATGT","dinuc31NSATTT","dinuc31NSCAAA","dinuc31NSCACC","dinuc31NSCACG","dinuc31NSCACT","dinuc31NSCAGA","dinuc31NSCATA","dinuc31NSCCAC","dinuc31NSCCCA","dinuc31NSCCCG","dinuc31NSCCCT","dinuc31NSCCGC","dinuc31NSCCTC","dinuc31NSCGAG","dinuc31NSCGCA","dinuc31NSCGCC","dinuc31NSCGCT","dinuc31NSCGGG","dinuc31NSCGTG","dinuc31NSCTAT","dinuc31NSCTCA","dinuc31NSCTCC","dinuc31NSCTCG","dinuc31NSCTGT","dinuc31NSCTTT","dinuc31NSGAAA","dinuc31NSGACA","dinuc31NSGAGC","dinuc31NSGAGG","dinuc31NSGAGT","dinuc31NSGATA","dinuc31NSGCAC","dinuc31NSGCCC","dinuc31NSGCGA","dinuc31NSGCGG","dinuc31NSGCGT","dinuc31NSGCTC","dinuc31NSGGAG","dinuc31NSGGCG","dinuc31NSGGGA","dinuc31NSGGGC","dinuc31NSGGGT","dinuc31NSGGTG","dinuc31NSGTAT","dinuc31NSGTCT","dinuc31NSGTGA","dinuc31NSGTGC","dinuc31NSGTGG","dinuc31NSGTTT","dinuc31NSTAAA","dinuc31NSTACA","dinuc31NSTAGA","dinuc31NSTATC","dinuc31NSTATG","dinuc31NSTATT","dinuc31NSTCAC","dinuc31NSTCCC","dinuc31NSTCGC","dinuc31NSTCTA","dinuc31NSTCTG","dinuc31NSTCTT","dinuc31NSTGAG","dinuc31NSTGCG","dinuc31NSTGGG","dinuc31NSTGTA","dinuc31NSTGTC","dinuc31NSTGTT","dinuc31NSTTAT","dinuc31NSTTCT","dinuc31NSTTGT","dinuc31NSTTTA","dinuc31NSTTTC","dinuc31NSTTTG",

        "Nsub","Nsynsub","MutRateStart","SubRateStart","MutRateNonSynStart","SubRateNonSynStart","MutRateSynStart","SubRateSynStart","MutRateEnd","SubRateEnd","MutRateNonSynEnd","SubRateNonSynEnd","MutRateSynEnd","SubRateSynEnd"

    };

    const string listSiteSpecificEvoStats[NSiteSpecificEvoStats] = {"ssNsub","ssNsynsub"};


    std::map<string,int> mapUsedParam;
    std::map<string,int> mapUsedSummaries;
    std::map<string,int> mapUsedAncSummaries;
    std::map<string,int> mapUsedAccessorySummaries;
    std::map<string,int> mapUsedEvoStats;
    std::map<string,int> mapUsedEvoAncStats;
    std::map<string,int> mapUsedSiteSpecificEvoStats;

    int OutPartialDistance,verbose, NusedAccessorySummaries,NusedSiteSpecificEvoStats, NusedEvoStats, NusedEvoAncStats, NusedParam, NusedSummaries, NusedAncSummaries, Ngenes, chainPointStart, chainPointEnd, chainPointEvery,Nthread, Niter, Nrun, threshold, Nrep, Nsite_codon, Ntaxa;

    std::vector<string> listGenes;
    std::vector<string> listChains;
    std::vector<int>    listPoints;
    std::vector<string> listOutput;
    std::vector<string> listSpecies;

    string model, controlfile, localcontrolfile, output, distance, transformation;

    void readInstructions();

    GlobalParameters(string model, string controlfile);
    GlobalParameters(); 
    virtual ~GlobalParameters();

protected:

private:
};

#endif // GLOBALPARAMETERS_H
