#include "SiteInterSubMatrix.h"

SiteInterSubMatrix::SiteInterSubMatrix(LocalParameters* lparam)
{

    this->lparam = lparam;

    //init SiteInterSubMatrix containers

    submatrixTreeSim = new double**[lparam->refTree->GetNnode()];
    mutmatrixTreeSim = new double**[lparam->refTree->GetNnode()];
    selmatrixTreeSim = new double**[lparam->refTree->GetNnode()];
    for (int node = 0 ; node < lparam->refTree->GetNnode() ; node ++ )
    {
        submatrixTreeSim[node] =  new double*[lparam->Nsite_nuc];
        mutmatrixTreeSim[node] = new  double*[lparam->Nsite_nuc];
        selmatrixTreeSim[node] = new double*[lparam->Nsite_nuc];
        for (int site_nuc = 0 ; site_nuc < lparam->Nsite_nuc ; site_nuc ++ )
        {
            submatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
            mutmatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
            selmatrixTreeSim[node][site_nuc] = new double[lparam->Nnucp];
        }
    }

    TotalSubRate = new double[lparam->refTree->GetNnode()];
    TotalMutRate = new double[lparam->refTree->GetNnode()];
    TotalSubRateNonSyn = new double[lparam->refTree->GetNnode()];
    TotalMutRateNonSyn = new double[lparam->refTree->GetNnode()];
    TotalSubRateSyn = new double[lparam->refTree->GetNnode()];
    TotalMutRateSyn = new double[lparam->refTree->GetNnode()];


}

SiteInterSubMatrix::~SiteInterSubMatrix()
{
    //dtor
}


void SiteInterSubMatrix::resetSubMatrix()
{



    for (int node = 0 ; node < lparam->refTree->GetNnode() ; node ++ )
    {

        TotalMutRate[node] = 0.0;
        TotalSubRate[node] = 0.0;
        TotalMutRateNonSyn[node] = 0.0;
        TotalSubRateNonSyn[node] = 0.0;
        TotalMutRateSyn[node] = 0.0;
        TotalSubRateSyn[node] = 0.0;

        for (int site_nuc = 0 ; site_nuc < lparam->Nsite_nuc ; site_nuc ++ )
        {
            for(int nuc = 0 ; nuc < Nnuc; nuc++)
            {
                submatrixTreeSim[node][site_nuc][nuc] = 0.0 ;
                mutmatrixTreeSim[node][site_nuc][nuc] = 0.0;
                selmatrixTreeSim[node][site_nuc][nuc] = 0.0;
            }
        }
    }

}




int  SiteInterSubMatrix::testCpGcontext(int inNnodeIndex, int insite, int innucFrom, int innucTo,int**CurrentNodeNucSequence)
{
// if insite is one before end seq and  cur_insite = C and cur_insite+1 = G and C goes to T then CpGHyp
// if insite is one after start seq and  cur_insite = G and cur_insite-1 = G and G goes to A then CpGHyp

    if (insite < lparam->Nsite_codon*3-1 && (CurrentNodeNucSequence[inNnodeIndex][insite] == 1 && CurrentNodeNucSequence[inNnodeIndex][insite+1] == 2) && innucFrom == 1 && innucTo == 3)
    {
        return 1;
    }
    else if (insite > 0 && (CurrentNodeNucSequence[inNnodeIndex][insite] == 2 && CurrentNodeNucSequence[inNnodeIndex][insite-1] == 1) && innucFrom == 2 && innucTo == 0)
    {
        return 2;
    }
    else
    {
        return -1;
    }
}

int  SiteInterSubMatrix::testTpAcontext(int inNnodeIndex, int insite, int innucFrom, int innucTo,int**CurrentNodeNucSequence)
{
// if insite is one before end seq and  cur_insite = T and cur_insite+1 = A and T goes to C then TpAHyp
// if insite is one after start seq and  cur_insite = A and cur_insite-1 = T and A goes to G then TpAHyp

    if (insite < lparam->Nsite_codon*3-1 && (CurrentNodeNucSequence[inNnodeIndex][insite] == 3 && CurrentNodeNucSequence[inNnodeIndex][insite+1] == 0) && innucFrom == 3 && innucTo == 1)
    {
        return 1; //TA->CA
    }
    else if (insite > 0 && (CurrentNodeNucSequence[inNnodeIndex][insite] == 0 && CurrentNodeNucSequence[inNnodeIndex][insite-1] == 3) && innucFrom == 0 && innucTo == 2)
    {
        return 2; //TA->TG
    }
    else
    {
        return -1;
    }

}

void SiteInterSubMatrix::UpdateSubMatrixTreeSim(int NodeIndex, int site_codon,int**CurrentNodeNucSequence)
{
    bool whole = false;
    if (site_codon == -1)
    {
        whole = true;

    }
    double deltaTotalSub = 0;
    double deltaTotalMut = 0;
    double deltaTotalSubNonSyn = 0;
    double deltaTotalMutNonSyn = 0;
    double deltaTotalSubSyn = 0;
    double deltaTotalMutSyn = 0;

// if -1, loop over all sites
    int site_codon_start = 0 ;
    int site_codon_end = lparam->Nsite_codon;
// else loop over site_codon-1 to site_codon+1: 3 codons to take CpG into account, which is in fact the worst case.
    if (!whole)
    {
        if (site_codon < lparam->Nsite_codon-2)
        {
            site_codon_end = site_codon + 2;
        }
        if (site_codon > 0)
        {
            site_codon_start =  site_codon - 1 ;
        }

        deltaTotalSub -= GetSubRate(NodeIndex,site_codon);
        deltaTotalMut -= GetMutRate(NodeIndex,site_codon);
        deltaTotalSubNonSyn -= GetSubRateNonSyn(NodeIndex,site_codon,CurrentNodeNucSequence);
        deltaTotalMutNonSyn -= GetMutRateNonSyn(NodeIndex,site_codon,CurrentNodeNucSequence);
        deltaTotalSubSyn -= GetSubRateSyn(NodeIndex,site_codon,CurrentNodeNucSequence);
        deltaTotalMutSyn -= GetMutRateSyn(NodeIndex,site_codon,CurrentNodeNucSequence);

    }

    int* nucposFrom = new int[3];
    int* nucposTo = new int[3];
    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ; site_codon_i++)
    {

        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
        for (int codonPos = 0; codonPos < 3; codonPos++)
        {
            nucposFrom[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
            nucposTo[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
        }

        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each nucleotide codon postions [0,2] we will be computing adjacent nucleotide.
        {
            int site_nuc =  site_nuc_start+codonPos;
            for (int nucTo = 0; nucTo < 4; nucTo++)
            {
                double  S = 0.0;
                double  MutRate = 0.0;
                double  SubRate = 0.0;
                double  MutRateNonSyn = 0.0;
                double  SubRateNonSyn = 0.0;
                double  MutRateSyn = 0.0;
                double  SubRateSyn = 0.0;

                if (nucposFrom[codonPos] != nucTo)
                {
                    nucposTo[codonPos] = nucTo;
                    int codonFrom = lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0], nucposFrom[1], nucposFrom[2]);
                    int codonTo = lparam->codonstatespace->GetCodonFromDNA(nucposTo[0], nucposTo[1], nucposTo[2]);
                    if(!lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1], nucposTo[2]))
                    {
                        double m;

                        auto it = lparam->gtrMap.find(NodeIndex);
                        if (it == lparam->gtrMap.end())
                        {
                            cerr << "Error when looking to the gtr along the tree\n";
                            exit(0);

                        }
                        if(it->second == 1)
                        {
                            //              cerr << "GTR " << 1 << "\n";
                            m = lparam->gtnr1[nucposFrom[codonPos]][nucposTo[codonPos]];
                        }
                        else if(it->second == 2)
                        {
                            //               cerr << "GTR " << 2 << "\n";
                            m = lparam->gtnr2[nucposFrom[codonPos]][nucposTo[codonPos]];
                        }
                        else
                        {
                            m = lparam->gtnr[nucposFrom[codonPos]][nucposTo[codonPos]];
                        }

                        int CpGcont = testCpGcontext(NodeIndex,site_nuc, nucposFrom[codonPos], nucposTo[codonPos],CurrentNodeNucSequence);
                        int TpAcont = testTpAcontext(NodeIndex,site_nuc, nucposFrom[codonPos], nucposTo[codonPos],CurrentNodeNucSequence);

                        if (CpGcont == 1 || CpGcont == 2)
                        {
                            m *=  lparam->lambda_CpG;
                        }

                        if (TpAcont == 1 || TpAcont == 2)
                        {
                            m *=  lparam->lambda_TpA;
                        }


                        if (m < lparam->TOOSMALL)
                        {
                            MutRate = lparam->TOOSMALL;
                        }
                        else
                        {
                            MutRate =  m;
                        }


                        MutRate *= lparam->lambda_TBL;
                        if (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
                        {
                            S = log(lparam->ssaaprofiles[lparam->alloc[site_codon_i]][lparam->codonstatespace->Translation(codonTo)] /
                                    lparam->ssaaprofiles[lparam->alloc[site_codon_i]][lparam->codonstatespace->Translation(codonFrom)] *
                                    lparam->codonprofile[codonTo] /
                                    lparam->codonprofile[codonFrom]);


                            SubRate = MutRate * lparam->lambda_omega *  lparam->omega;
                            MutRateNonSyn = MutRate;
                            SubRateNonSyn = SubRate;

                        }
                        else
                        {
                            S = log(lparam->codonprofile[codonTo] / lparam->codonprofile[codonFrom]);
                            SubRate = MutRate;
                            MutRateSyn = MutRate;
                            SubRateSyn = SubRate;

                        }

                        //NUMERICAL SECURITY
                        if (fabs(S) < lparam->TOOSMALL)
                        {
                            SubRate /=  (1.0 - (S / 2));
                        }
                        else if (S > lparam->TOOLARGE)
                        {
                            SubRate *= S;

                        }
                        else if (S < lparam->TOOLARGENEGATIVE)
                        {
                            SubRate = 0.0;
                        }
                        else
                        {
                            SubRate *=  (S / (1.0 -exp(-S)));
                        }

                        if (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
                        {

                            SubRateNonSyn = SubRate;
                        }
                        else
                        {

                            SubRateSyn = SubRate;
                        }

                        //NUMERICAL SECURITY linked to substitution process
                        if (SubRate<0)
                        {
                            cerr << "negative entry in matrix\n";
                            cerr << "S: " << S << "\n";
                            exit(1);
                        }
                        if (isinf(SubRate))
                        {
                            cerr << "isinf\n";
                            cerr << site_nuc << "\t" << nucTo << "\t" << MutRate << "\t" << S << "\n";
                            exit(1);
                        }
                        if (isnan(SubRate))
                        {
                            cerr << "isnan\n";
                            cerr << site_nuc << "\t" << nucTo << "\t" << MutRate << "\t" << S << "\n";
                            exit(1);
                        }
                        if(lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1], nucposTo[2]))
                        {
                            cerr << nucposTo[0] << " " <<  nucposTo[1] << " " << nucposTo[2] << " " <<  MutRate << " " <<  SubRate <<  " " << submatrixTreeSim[0][site_nuc][nucTo] << "\n";
                            exit(0);
                        }
                    }
                    nucposTo[codonPos]  = nucposFrom[codonPos];
                }

                deltaTotalSub += SubRate;
                deltaTotalMut += MutRate;
                deltaTotalSubNonSyn += SubRateNonSyn;
                deltaTotalMutNonSyn += MutRateNonSyn;
                deltaTotalSubSyn += SubRateSyn;
                deltaTotalMutSyn += MutRateSyn;

                mutmatrixTreeSim[NodeIndex][site_nuc][nucTo] = MutRate;
                selmatrixTreeSim[NodeIndex][site_nuc][nucTo] = S;
                submatrixTreeSim[NodeIndex][site_nuc][nucTo] = SubRate;

//              cerr << "MutRate" << submatrix_mut[inNodeIndex][site_nuc][nucTo] << "\n";
                if (lparam->opt == 3)
                {
                    submatrixTreeSim[NodeIndex][site_nuc][nucTo] = MutRate;
                }

            }
        }
    }

    if (whole)
    {
        TotalSubRate[NodeIndex] = deltaTotalSub;
        TotalMutRate[NodeIndex] = deltaTotalMut;
        TotalSubRateNonSyn[NodeIndex] = deltaTotalSubNonSyn;
        TotalMutRateNonSyn[NodeIndex] = deltaTotalMutNonSyn;
        TotalSubRateSyn[NodeIndex] = deltaTotalSubSyn;
        TotalMutRateSyn[NodeIndex] = deltaTotalMutSyn;
    }
    else
    {
        TotalSubRate[NodeIndex] += deltaTotalSub;
        TotalMutRate[NodeIndex] += deltaTotalMut;
        TotalSubRateNonSyn[NodeIndex] += deltaTotalSubNonSyn;
        TotalMutRateNonSyn[NodeIndex] += deltaTotalMutNonSyn;
        TotalSubRateSyn[NodeIndex] += deltaTotalSubSyn;
        TotalMutRateSyn[NodeIndex] += deltaTotalMutSyn;
    }

    if (lparam->opt == 3)
    {
        TotalSubRate[NodeIndex] = TotalMutRate[NodeIndex];
    }

    delete [] nucposFrom;
    delete [] nucposTo;

}

void SiteInterSubMatrix::transfertTotalRate(int sourceNodeIndex,int sinkNodeIndex)
{
    TotalSubRate[sinkNodeIndex] = TotalSubRate[sourceNodeIndex] ;
    TotalMutRate[sinkNodeIndex] = TotalMutRate[sourceNodeIndex] ;

}

void SiteInterSubMatrix::transfertNodeMatrix(int sourceNodeIndex, int sinkNodeIndex,int site_nuc)
{
    for (int i = 0 ; i < 4; i++)
    {
        submatrixTreeSim[sinkNodeIndex][site_nuc][i] = submatrixTreeSim[sourceNodeIndex][site_nuc][i];
        mutmatrixTreeSim[sinkNodeIndex][site_nuc][i] = mutmatrixTreeSim[sourceNodeIndex][site_nuc][i];
        selmatrixTreeSim[sinkNodeIndex][site_nuc][i] = selmatrixTreeSim[sourceNodeIndex][site_nuc][i];
    }
}


int  SiteInterSubMatrix::testContextDinuc(int NodeIndex, int site_nuc, int* context, int nucTo, int** CurrentNodeNucSequence)
{

    int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
    int MutPos = -1;

    //Find dinucleotide mutating position for the actual partern to be tested
    if (context[0] == context[2] && context[1] == nucFrom && context[3] == nucTo)
    {
        MutPos = 1;

    }
    else if (context[0] == nucFrom && context[1] == context[3] && context[2] == nucTo)
    {
        MutPos = 0;

    }

    //Test if the mutating patern is recovered from the actual sequence
    if (MutPos == 0 && site_nuc < lparam->Nsite_nuc-1)
    {
        if (CurrentNodeNucSequence[NodeIndex][site_nuc] == context[MutPos] && CurrentNodeNucSequence[NodeIndex][site_nuc+1] == context[1])
        {
            return 1;

        }

    }
    else if (MutPos == 1 && site_nuc > 0)
    {
        if (CurrentNodeNucSequence[NodeIndex][site_nuc] == context[MutPos] && CurrentNodeNucSequence[NodeIndex][site_nuc-1] == context[0])
        {
            return 1;

        }

    }
    else
    {

        return -1;

    }
}


void SiteInterSubMatrix::findCodonContext(int NodeIndex,  int site_nuc,int nucFrom, int nucTo, int &pos1From, int &pos2From, int &pos3From, int &pos1To, int &pos2To, int &pos3To,int**CurrentNodeNucSequence )
{
    // this function aims at determining adjacent codon from one mutation

    int codonPos = (site_nuc % 3);
    if (codonPos == 0)
    {
        pos1From = nucFrom;
        pos1To = nucTo;

        pos2From = CurrentNodeNucSequence[NodeIndex][site_nuc + 1];
        pos2To = CurrentNodeNucSequence[NodeIndex][site_nuc + 1] ;

        pos3From = CurrentNodeNucSequence[NodeIndex][NodeIndex + 2];
        pos3To = CurrentNodeNucSequence[NodeIndex][site_nuc + 2];

    }
    else if (codonPos == 1)
    {
        pos1From = CurrentNodeNucSequence[NodeIndex][site_nuc - 1];
        pos1To = CurrentNodeNucSequence[NodeIndex][site_nuc - 1];

        pos2From = nucFrom;
        pos2To = nucTo;

        pos3From = CurrentNodeNucSequence[NodeIndex][site_nuc + 1];
        pos3To = CurrentNodeNucSequence[NodeIndex][site_nuc + 1];

    }
    else if (codonPos == 2)
    {
        pos1From = CurrentNodeNucSequence[NodeIndex][site_nuc - 2];
        pos1To = CurrentNodeNucSequence[NodeIndex][site_nuc - 2];

        pos2From = CurrentNodeNucSequence[NodeIndex][site_nuc - 1];
        pos2To = CurrentNodeNucSequence[NodeIndex][site_nuc - 1];

        pos3From = nucFrom;
        pos3To = nucTo;
    }
}

//
//double SiteInterSubMatrix::GetMutRateSyn(int NodeIndex,int** CurrentNodeNucSequence)
//{
//    int nucFrom, pos1From, pos2From, pos3From, pos1To, pos2To, pos3To, codonFrom, codonTo;
//    double sum  = 0;
//    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++)
//    {
//        nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
//        for (int nucTo = 0; nucTo < 4; nucTo++)
//        {
//            findCodonContext(NodeIndex, site_nuc, nucFrom, nucTo, pos1From, pos2From, pos3From, pos1To, pos2To, pos3To,CurrentNodeNucSequence);
//            if(!lparam->codonstatespace->CheckStop(pos1To, pos2To, pos3To))
//            {
//                codonFrom = lparam->codonstatespace->GetCodonFromDNA(pos1From, pos2From, pos3From);
//                codonTo = lparam->codonstatespace->GetCodonFromDNA(pos1To, pos2To, pos3To) ;
//                if (lparam->codonstatespace->Synonymous(codonFrom,codonTo))
//                {
//                    sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
//                }
//            }
//        }
//    }
//    return sum;
//}
//
//double SiteInterSubMatrix::GetSubRateCpG(int NodeIndex,int** CurrentNodeNucSequence)
//{
//    double sum  = 0.0;
//    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++)
//    {
//        int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
//        for (int nucTo = 0; nucTo < 4; nucTo++)
//        {
//            if(testCpGcontext(NodeIndex, site_nuc, nucFrom, nucTo,CurrentNodeNucSequence) > 0)
//            {
//                sum += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
//            }
//        }
//    }
//    return sum;
//}
//
//double SiteInterSubMatrix::GetMutRateCpG(int NodeIndex,int** CurrentNodeNucSequence)
//{
//    double sum  = 0.0;
//    for (int site_nuc = 0; site_nuc < lparam->Nsite_nuc; site_nuc++)
//    {
//        int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
//        for (int nucTo = 0; nucTo < 4; nucTo++)
//        {
//            if(testCpGcontext(NodeIndex, site_nuc, nucFrom, nucTo,CurrentNodeNucSequence)> 0)
//            {
//                sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
//            }
//        }
//    }
//    return sum;
//}


double SiteInterSubMatrix::GetMutRate(int NodeIndex, int site_codon)
{
    double sum  = 0.0;
    int site_codon_start = 0 ;
    int site_codon_end = lparam->Nsite_codon;
    if (site_codon > -1)
    {
        if (site_codon < lparam->Nsite_codon-1)
        {
            site_codon_end = site_codon + 2;
        }
        if (site_codon > 0)
        {
            site_codon_start =  site_codon - 1 ;
        }
    }
    int site_nuc_start = site_codon_start * 3;
    int site_nuc_end =  site_codon_end * 3;
    for (int site_nuc = site_nuc_start ; site_nuc < site_nuc_end; site_nuc++)
    {
        for (int nucTo = 0; nucTo < 4; nucTo++)
        {
            sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];
        }
    }
    //TotalMut[NodeIndex] = sum;
    return sum;
}


double SiteInterSubMatrix::GetSubRate(int NodeIndex, int site_codon)
{
    double sum  = 0;
    int site_codon_start = 0 ;
    int site_codon_end =  lparam->Nsite_codon;
    if (site_codon > -1)
    {
        if (site_codon < lparam->Nsite_codon-1)
        {
            site_codon_end = site_codon + 2;
        }
        if (site_codon > 0)
        {
            site_codon_start =  site_codon - 1 ;
        }
    }
    int site_nuc_start = site_codon_start * 3;
    int site_nuc_end =  site_codon_end * 3;
    for (int site_nuc = site_nuc_start ; site_nuc < site_nuc_end; site_nuc++)
    {
        for (int nucTo = 0; nucTo < 4; nucTo++)
        {
            sum += submatrixTreeSim[NodeIndex][site_nuc][nucTo];
        }
    }
    //TotalSub[NodeIndex] = sum;
    return sum;
}

double SiteInterSubMatrix::GetSubRateNonSyn(int NodeIndex, int site_codon,int** CurrentNodeNucSequence)
{

    double sum  = 0.0;
    int site_codon_start = 0 ;
    int site_codon_end = lparam->Nsite_codon;


    if (site_codon > -1)
    {
        if (site_codon < lparam->Nsite_codon-1)
        {
            site_codon_end = site_codon + 2;
        }
        if (site_codon > 0)
        {
            site_codon_start =  site_codon - 1;
        }
    }


    int* nucposFrom = new int[3];
    int* nucposTo = new int[3];

    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ; site_codon_i++)
    {

        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
        for (int codonPos = 0; codonPos < 3; codonPos++)
        {
            nucposFrom[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
            nucposTo[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
        }

        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each nucleotide codon postions [0,2] we will be computing adjacent nucleotide.
        {
            int site_nuc =  site_nuc_start+codonPos;
            for (int nucTo = 0; nucTo < 4; nucTo++)
            {
                if (nucposFrom[codonPos] != nucTo)
                {
                    nucposTo[codonPos] = nucTo;

                    if(lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1], nucposTo[2]))
                    {

                        int codonFrom = lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0], nucposFrom[1], nucposFrom[2]);
                        int codonTo = lparam->codonstatespace->GetCodonFromDNA(nucposTo[0], nucposTo[1], nucposTo[2]);

                        if (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
                        {

                            sum += submatrixTreeSim[NodeIndex][site_nuc][nucTo];

                        }

                    }

                }

            }
        }

    }

    delete [] nucposFrom;
    delete [] nucposTo;


    return sum;
}

double SiteInterSubMatrix::GetSubRateSyn(int NodeIndex, int site_codon,int** CurrentNodeNucSequence)
{
    double sum  = 0.0;
    int site_codon_start = 0 ;
    int site_codon_end = lparam->Nsite_codon;


    if (site_codon > -1)
    {
        if (site_codon < lparam->Nsite_codon-1)
        {
            site_codon_end = site_codon + 2;
        }
        if (site_codon > 0)
        {
            site_codon_start =  site_codon - 1;
        }
    }


    int* nucposFrom = new int[3];
    int* nucposTo = new int[3];

    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ; site_codon_i++)
    {

        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
        for (int codonPos = 0; codonPos < 3; codonPos++)
        {
            nucposFrom[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
            nucposTo[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
        }

        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each nucleotide codon postions [0,2] we will be computing adjacent nucleotide.
        {
            int site_nuc =  site_nuc_start+codonPos;
            for (int nucTo = 0; nucTo < 4; nucTo++)
            {
                if (nucposFrom[codonPos] != nucTo)
                {
                    nucposTo[codonPos] = nucTo;

                    if(lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1], nucposTo[2]))
                    {

                        int codonFrom = lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0], nucposFrom[1], nucposFrom[2]);
                        int codonTo = lparam->codonstatespace->GetCodonFromDNA(nucposTo[0], nucposTo[1], nucposTo[2]);

                        if (lparam->codonstatespace->Synonymous(codonFrom,codonTo))
                        {

                            sum += submatrixTreeSim[NodeIndex][site_nuc][nucTo];

                        }



                    }

                }

            }
        }

    }

    delete [] nucposFrom;
    delete [] nucposTo;


    return sum;
}

double SiteInterSubMatrix::GetMutRateNonSyn(int NodeIndex, int site_codon,int** CurrentNodeNucSequence)
{
    double sum  = 0.0;
    int site_codon_start = 0 ;
    int site_codon_end = lparam->Nsite_codon;


    if (site_codon > -1)
    {
        if (site_codon < lparam->Nsite_codon-1)
        {
            site_codon_end = site_codon + 2;
        }
        if (site_codon > 0)
        {
            site_codon_start =  site_codon - 1;
        }
    }


    int* nucposFrom = new int[3];
    int* nucposTo = new int[3];

    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ; site_codon_i++)
    {

        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
        for (int codonPos = 0; codonPos < 3; codonPos++)
        {
            nucposFrom[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
            nucposTo[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
        }

        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each nucleotide codon postions [0,2] we will be computing adjacent nucleotide.
        {
            int site_nuc =  site_nuc_start+codonPos;
            for (int nucTo = 0; nucTo < 4; nucTo++)
            {
                if (nucposFrom[codonPos] != nucTo)
                {
                    nucposTo[codonPos] = nucTo;

                    if(lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1], nucposTo[2]))
                    {

                        int codonFrom = lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0], nucposFrom[1], nucposFrom[2]);
                        int codonTo = lparam->codonstatespace->GetCodonFromDNA(nucposTo[0], nucposTo[1], nucposTo[2]);

                        if (!lparam->codonstatespace->Synonymous(codonFrom,codonTo))
                        {

                            sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];

                        }

                    }

                }

            }
        }

    }

    delete [] nucposFrom;
    delete [] nucposTo;


    return sum;
}

double SiteInterSubMatrix::GetMutRateSyn(int NodeIndex, int site_codon,int** CurrentNodeNucSequence)
{
    double sum  = 0.0;
    int site_codon_start = 0 ;
    int site_codon_end = lparam->Nsite_codon;


    if (site_codon > -1)
    {
        if (site_codon < lparam->Nsite_codon-1)
        {
            site_codon_end = site_codon + 2;
        }
        if (site_codon > 0)
        {
            site_codon_start =  site_codon - 1;
        }
    }



    int* nucposFrom = new int[3];
    int* nucposTo = new int[3];

    for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ; site_codon_i++)
    {

        int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
        for (int codonPos = 0; codonPos < 3; codonPos++)
        {
            nucposFrom[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
            nucposTo[codonPos] = CurrentNodeNucSequence[NodeIndex][site_nuc_start+codonPos];
        }

        for ( int codonPos = 0; codonPos < 3; codonPos++ )   // for each nucleotide codon postions [0,2] we will be computing adjacent nucleotide.
        {
            int site_nuc =  site_nuc_start+codonPos;
            for (int nucTo = 0; nucTo < 4; nucTo++)
            {
                if (nucposFrom[codonPos] != nucTo)
                {
                    nucposTo[codonPos] = nucTo;

                    if(lparam->codonstatespace->CheckStop(nucposTo[0], nucposTo[1], nucposTo[2]))
                    {

                        int codonFrom = lparam->codonstatespace->GetCodonFromDNA(nucposFrom[0], nucposFrom[1], nucposFrom[2]);
                        int codonTo = lparam->codonstatespace->GetCodonFromDNA(nucposTo[0], nucposTo[1], nucposTo[2]);

                        if (lparam->codonstatespace->Synonymous(codonFrom,codonTo))
                        {

                            sum += mutmatrixTreeSim[NodeIndex][site_nuc][nucTo];

                        }

                    }

                }

            }
        }

    }

    delete [] nucposFrom;
    delete [] nucposTo;


    return sum;
}



