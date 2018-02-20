#include "TreeSimulator.h"

TreeSimulator::TreeSimulator(LocalParameters* lparam, SiteInterSubMatrix* submatrix, AncestralSequence* ancestralseq)
{
    this->lparam = lparam;
    this->submatrix = submatrix;
    this->ancestralseq = ancestralseq;
    this->treeEvoStats = new EvolHistStatistics(this->lparam);
    this->rootBranchEvoStats = new EvolHistStatistics(this->lparam);


    CurrentNodeCodonSequence = new int*[lparam->refTree->GetNnode()];
    CurrentNodeNucSequence = new int*[lparam->refTree->GetNnode()];

    for (int node = 0 ; node < lparam->refTree->GetNnode(); node++)
    {
        CurrentNodeCodonSequence[node]= new int[lparam->Nsite_codon];
        CurrentNodeNucSequence[node]= new int[lparam->Nsite_codon*3];
    }

    CurrentLeafNodeCodonSequences= new int * [lparam->Ntaxa];
    CurrentLeafNodeNucSequence= new int * [lparam->Ntaxa];

    for (int taxa = 0 ; taxa < lparam->Ntaxa ; taxa ++ )
    {
        CurrentLeafNodeCodonSequences[taxa] = new int [lparam->Nsite_codon];
        CurrentLeafNodeNucSequence[taxa] = new int [lparam->Nsite_nuc];
    }


    CurrentAncestralCodonSequence = new int**[11];
    for (int point_i = 0 ; point_i < 11 ; point_i ++)
    {
        CurrentAncestralCodonSequence[point_i] = new int *[1];
        CurrentAncestralCodonSequence[point_i][0] = new int [lparam->Nsite_codon];
    }


}

TreeSimulator::~TreeSimulator()
{
    //dtor
}

void TreeSimulator::resetSimulator()
{

    int verbose = lparam->verbose;

    if(verbose)
    {
        cerr << "resetSimulator1\n";
    }
    if(verbose)
    {
        cerr << lparam->refTree->GetNnode() << "resetSimulator1.1\n";
    }
    for (int node = 0 ; node < lparam->refTree->GetNnode(); node++)
    {
        for(int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon++)
        {

            CurrentNodeCodonSequence[node][site_codon] = -2 ;

            for(int j = 0 ; j < 3; j ++)
            {
                CurrentNodeNucSequence[node][site_codon*3+j] = -2;

            }
        }
    }

    if(verbose)
    {
        cerr << "resetSimulator2\n";
    }
    // reset nodeleaf sequences
    for (int taxa = 0 ; taxa < lparam->Ntaxa ; taxa ++ )
    {
        for(int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon++)
        {
            CurrentLeafNodeCodonSequences[taxa][site_codon] = -2 ;
            for(int j = 0 ; j < 3; j ++)
            {
                CurrentLeafNodeNucSequence[taxa][site_codon*3+j] = -2;
            }
        }
    }

    if(verbose)
    {
        cerr << "resetSimulator3\n";
    }

    for (int point_i = 0 ; point_i < 11 ; point_i ++)
    {

        for(int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon++)
        {
            CurrentAncestralCodonSequence[point_i][0][site_codon] = -2;
        }
    }


}

void TreeSimulator::GetNewSimulatedCodonAlignment()
{

    int verbose = lparam->verbose;

    submatrix->resetSubMatrix();
    if(verbose)
    {
        cerr << "submatrix->resetSubMatrix()\n";
    }

    rootBranchEvoStats->resetEvoStats();
    if(verbose)
    {
        cerr << "rootBranchEvoStats->resetEvoStats()\n";
    }

    treeEvoStats->resetEvoStats();
    if(verbose)
    {
        cerr << "treeEvoStats->resetEvoStats()\n";
    }

    resetSimulator();
    if(verbose)
    {
        cerr << "resetSimulator()\n";
    }

    ancestralseq->GetNewStationaryCodonSequence();
    if(verbose)
    {
        cerr << "ancestralseq->GetNewStationaryCodonSequence()\n";
    }

    SetAncestralSequence();
    if(verbose)
    {
        cerr << "SetAncestralSequence()\n";
    }

    //launch recursive simulation on a phylogenetic tree
    ComputeRecursiveSimulation(lparam->refTree->GetRoot());
    if(verbose)
    {
        cerr << "ComputeRecursiveSimulation\n";
    }
    //register mappingstats

    resetEvoStatVectors();
    if(verbose)
    {
        cerr << "resetEvoStatVectors()\n";
    }

    rootBranchEvoStats->GetEvoAncStats();
    if(verbose)
    {
        cerr << "ancestralSequenceEvohist->GetEvoAncStats()\n";
    }

    treeEvoStats->GetEvoStats();
    if(verbose)
    {
        cerr << "evohist->GetEvoStats()\n";
    }

    treeEvoStats->GetSiteSpecificEvoStats();
    if(verbose)
    {
        cerr << "evohist->GetSiteSpecificEvoStats()\n";
    }
}


void TreeSimulator::resetEvoStatVectors()
{
    lparam->ancevostats.clear();
    lparam->ancevostats.shrink_to_fit();

    lparam->evostats.clear();
    lparam->evostats.shrink_to_fit();

    lparam->sitespecificevostats.clear();
    lparam->sitespecificevostats.shrink_to_fit();

}

void TreeSimulator::SetAncestralSequence()
{
    int NodeIndex = lparam->refTree->GetRoot()->GetNode()->GetIndex();
    for (int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon++)
    {
        this->CurrentNodeCodonSequence[NodeIndex][site_codon] = this->ancestralseq->GetCurrentAncestralCodonSequence(site_codon);
        for(int j = 0 ; j < 3; j++)
        {
            this->CurrentNodeNucSequence[NodeIndex][site_codon*3+j] = this->ancestralseq->GetCurrentAncestralNucSequence(site_codon*3+j);
        }
    }


}

void TreeSimulator::RegisterSubTreeSim(int NodeIndex, int site_nuc, int nucTo)
{

    int verbose = lparam->verbose;

    int nucFrom = CurrentNodeNucSequence[NodeIndex][site_nuc];
    int CodonPos = site_nuc % 3;
    //3 codons length context

    int* nucposFrom = new int[9];
    int* nucposTo = new int[9];


    int* codonFrom = new int[3];
    int* codonTo = new int[3];

    int site_codon = (int) (site_nuc / 3);
    int site_codon_start = site_codon-1 ;
    int site_codon_end = site_codon+2;


    ///
    /// GET 3 ADJACENT CODONS FOR CODON INTERFACE STATS
    {
        int codon_count = 0;
        for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ; site_codon_i++)
        {

            int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
            for (int codonPos = 0; codonPos < 3; codonPos++)
            {

                if (site_codon_i >= 0 && site_codon_i < lparam->Nsite_codon)
                {


                    nucposFrom[codon_count*3+codonPos] = CurrentNodeNucSequence[NodeIndex][site_codon_i*3+codonPos];
                    //nucposTo[codon_count*3+codonPos] = CurrentNodeNucSequence[NodeIndex][site_codon_i*3+codonPos];
                    nucposTo[codon_count*3+codonPos] = nucposFrom[codon_count*3+codonPos];

                }
                else
                {

                    nucposFrom[codon_count*3+codonPos] = -1;
                    nucposTo[codon_count*3+codonPos] = -1;

                }



            }
            if (site_codon_i >= 0 && site_codon_i < lparam->Nsite_codon)
            {
                codonTo[codon_count] = lparam->codonstatespace->GetCodonFromDNA(nucposTo[codon_count*3], nucposTo[codon_count*3+1], nucposTo[codon_count*3+2]);
                codonFrom[codon_count] = codonTo[codon_count];
            }
            else
            {
                codonTo[codon_count] = -1;
                codonFrom[codon_count] = -1;
            }
            codon_count++;
        }
    }

    int site_nuc_To = site_nuc%3+3;
    nucposTo[site_nuc_To] = nucTo;
    codonTo[1] = lparam->codonstatespace->GetCodonFromDNA(nucposTo[3], nucposTo[4], nucposTo[5]);




    if(lparam->codonstatespace->CheckStop(nucposTo[3], nucposTo[4], nucposTo[5]))
    {
        cerr << "error while registring a substitution\n";
        cerr << nucTo << "\n";
        cerr << nucposTo[3] << " " <<  nucposTo[4] << " " <<   nucposTo[5] << " " << codonTo[1]<<"\n";
        cerr << nucposFrom[3] << " " <<  nucposFrom[4] << " " <<   nucposFrom[5] << " "<< codonFrom[1] <<"\n";
        exit(0);
    }


    // Mapping statitics
    //Statistics on the ancestral sequence
    if(NodeIndex == lparam->refTree->GetRoot()->GetNode()->GetIndex())
    {
        rootBranchEvoStats->Nsub++;
        //rootBranchEvoStats->ssNsub[site_codon]++;
        if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
        {
            rootBranchEvoStats->Nsynsub++;
        }


        if (CodonPos == 0)
        {
            rootBranchEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
            rootBranchEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;
            rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
            rootBranchEvoStats->dinuc_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
            if (site_nuc>0)
            {
                rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                rootBranchEvoStats->dinuc_stat[2][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
            }

//            if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
//            {
//
//                //ancestralSequenceEvohist->ssNsynsub[site_codon]++;
//
////                rootBranchEvoStats->gtnrSyn_stat[3][nucFrom][nucTo]++;
////                rootBranchEvoStats->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;
////
////                rootBranchEvoStats->dinucSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                rootBranchEvoStats->dinucSyn_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                if (site_nuc>0)
////                {
////                    rootBranchEvoStats->dinucSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                    rootBranchEvoStats->dinucSyn_stat[2][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                }
//
//            }
//            else
//            {
//                rootBranchEvoStats->gtnrNSyn_stat[3][nucFrom][nucTo]++;
//                rootBranchEvoStats->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;
//
//                rootBranchEvoStats->dinucNSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                rootBranchEvoStats->dinucNSyn_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                if (site_nuc>0)
//                {
//                    rootBranchEvoStats->dinucNSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//                    rootBranchEvoStats->dinucNSyn_stat[2][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//                }
//
//            }

        }
        else if(CodonPos == 1)
        {

            rootBranchEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
            rootBranchEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;

            rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
            rootBranchEvoStats->dinuc_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

            rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
            rootBranchEvoStats->dinuc_stat[0][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //



//            if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
//            {
//
//                //ancestralSequenceEvohist->ssNsynsub[site_codon]++;
//
//                rootBranchEvoStats->gtnrSyn_stat[3][nucFrom][nucTo]++;
//                rootBranchEvoStats->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;
//
//                rootBranchEvoStats->dinucSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                rootBranchEvoStats->dinucSyn_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
//                rootBranchEvoStats->dinucSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//                rootBranchEvoStats->dinucSyn_stat[0][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
//            }
//            else
//            {
//                rootBranchEvoStats->gtnrNSyn_stat[3][nucFrom][nucTo]++;
//                rootBranchEvoStats->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;
//
//                rootBranchEvoStats->dinucNSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                rootBranchEvoStats->dinucNSyn_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
//                rootBranchEvoStats->dinucNSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//                rootBranchEvoStats->dinucNSyn_stat[0][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
//
//            }

        }
        else if (CodonPos == 2)
        {

            rootBranchEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
            rootBranchEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;

            rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
            rootBranchEvoStats->dinuc_stat[1][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

            if (site_codon < lparam->Nsite_codon-2)
            {
                rootBranchEvoStats->dinuc_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                rootBranchEvoStats->dinuc_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
            }


//            if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
//            {
//                //ancestralSequenceEvohist->ssNsynsub[site_codon]++;
//
//                rootBranchEvoStats->gtnrSyn_stat[3][nucFrom][nucTo]++;
//                rootBranchEvoStats->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;
//
//                rootBranchEvoStats->dinucSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//                rootBranchEvoStats->dinucSyn_stat[1][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
//                if (site_codon < lparam->Nsite_codon-2)
//                {
//                    rootBranchEvoStats->dinucSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                    rootBranchEvoStats->dinucSyn_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                }
//
//            }
//            else
//            {
//                rootBranchEvoStats->gtnrNSyn_stat[3][nucFrom][nucTo]++;
//                rootBranchEvoStats->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;
//
//                rootBranchEvoStats->dinucNSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//                rootBranchEvoStats->dinucNSyn_stat[1][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
//                if (site_codon < lparam->Nsite_codon-2)
//                {
//                    rootBranchEvoStats->dinucNSyn_stat[3][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                    rootBranchEvoStats->dinucNSyn_stat[CodonPos][rootBranchEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][rootBranchEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                }
//
//            }


        }



    }
    else
    {

        // Non-root events
        treeEvoStats->Nsub++;
        treeEvoStats->ssNsub[site_codon]++;

        treeEvoStats->gtnr_stat[3][nucFrom][nucTo]++;
        treeEvoStats->gtnr_stat[CodonPos][nucFrom][nucTo]++;

//            treeEvoStats->codon_stat[3][codonFrom[1]][codonTo[1]]++;
//            treeEvoStats->codon_stat[CodonPos][codonFrom[1]][codonTo[1]]++;
//            treeEvoStats->ssgtnr_stat[site_codon][3][nucFrom][nucTo]++;
//            treeEvoStats->ssgtnr_stat[site_codon][CodonPos][nucFrom][nucTo]++;
//            treeEvoStats->sscodon_stat[site_codon][3][codonFrom[1]][codonTo[1]]++;
//            treeEvoStats->sscodon_stat[site_codon][CodonPos][codonFrom[1]][codonTo[1]]++;

        if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
        {
            treeEvoStats->Nsynsub++;
            treeEvoStats->ssNsynsub[site_codon]++;
            treeEvoStats->gtnrSyn_stat[3][nucFrom][nucTo]++;
            treeEvoStats->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;
//          treeEvoStats->ssgtnrSyn_stat[site_codon][3][nucFrom][nucTo]++;
//          treeEvoStats->ssgtnrSyn_stat[site_codon][CodonPos][nucFrom][nucTo]++;
        }
        else
        {
//
            treeEvoStats->gtnrNSyn_stat[3][nucFrom][nucTo]++;
            treeEvoStats->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;
//          treeEvoStats->ssgtnrNSyn_stat[site_codon][3][nucFrom][nucTo]++;
//          treeEvoStats->ssgtnrNSyn_stat[site_codon][CodonPos][nucFrom][nucTo]++;

        }
        if(verbose)
        {
            cerr << "10\n";
        }
//
        if (CodonPos == 0)
        {

            treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
            treeEvoStats->dinuc_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

////                    evohist->ssdinuc_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                    evohist->ssdinuc_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
//            // SECOND FRAM FOR DINUC
            if (site_nuc>0)
            {
                treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                treeEvoStats->dinuc_stat[2][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

////                        evohist->ssdinuc_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                        evohist->ssdinuc_stat[site_codon][2][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
            }

//
            if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
            {
//
                treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                treeEvoStats->dinucSyn_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
////                        evohist->ssdinucSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                        evohist->ssdinucSyn_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
                if (site_nuc>0)
                {
                    treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    treeEvoStats->dinucSyn_stat[2][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
////                            evohist->ssdinucSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                            evohist->ssdinucSyn_stat[site_codon][2][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                }
//
            }
            else
            {
//
//
                treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                treeEvoStats->dinucNSyn_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
////                        evohist->ssdinucNSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                        evohist->ssdinucNSyn_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

                if (site_nuc>0)
                {
                    treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    treeEvoStats->dinucNSyn_stat[2][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
////                            evohist->ssdinucNSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                            evohist->ssdinucNSyn_stat[site_codon][2][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                }

            }
//
        }
        else if(CodonPos == 1)
        {

//
            treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
            treeEvoStats->dinuc_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
////                    evohist->ssdinuc_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                    evohist->ssdinuc_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
//
            treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
            treeEvoStats->dinuc_stat[0][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
//
////                    evohist->ssdinuc_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                    evohist->ssdinuc_stat[site_codon][0][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
            if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
            {
//
//
                treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                treeEvoStats->dinucSyn_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
////                        evohist->ssdinucSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                        evohist->ssdinucSyn_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
//
                treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                treeEvoStats->dinucSyn_stat[0][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
////                        evohist->ssdinucSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                        evohist->ssdinucSyn_stat[site_codon][0][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
            }
            else
            {
                treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                treeEvoStats->dinucNSyn_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

////                        evohist->ssdinucNSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                        evohist->ssdinucNSyn_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//
//
                treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                treeEvoStats->dinucNSyn_stat[0][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

//
////                       evohist->ssdinucNSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                       evohist->ssdinucNSyn_stat[site_codon][0][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
            }
//
        }
        else if (CodonPos == 2)
        {

            treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
            treeEvoStats->dinuc_stat[1][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

////                    evohist->ssdinuc_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                    evohist->ssdinuc_stat[site_codon][1][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//
            if (site_codon < lparam->Nsite_codon-2)
            {
                treeEvoStats->dinuc_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                treeEvoStats->dinuc_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

////                        evohist->ssdinuc_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                        evohist->ssdinuc_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
            }
//
//
            if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1]))
            {
                treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                treeEvoStats->dinucSyn_stat[1][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

////                        evohist->ssdinucSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
////                        evohist->ssdinucSyn_stat[site_codon][1][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                if (site_codon < lparam->Nsite_codon-2)
                {
                    treeEvoStats->dinucSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    treeEvoStats->dinucSyn_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

//                            evohist->ssdinucSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
//                            evohist->ssdinucSyn_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                }

            }
            else
            {

                treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                treeEvoStats->dinucNSyn_stat[1][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

//                //evohist->ssdinucNSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
//                //evohist->ssdinucNSyn_stat[site_codon][1][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                if (site_codon < lparam->Nsite_codon-2)
                {
                    treeEvoStats->dinucNSyn_stat[3][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    treeEvoStats->dinucNSyn_stat[CodonPos][treeEvoStats->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][treeEvoStats->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

////                            evohist->ssdinucNSyn_stat[site_codon][3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
////                            evohist->ssdinucNSyn_stat[site_codon][CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++;
                }

            }

        }

    }




    //The evolving sequence is updated here



    CurrentNodeNucSequence[NodeIndex][site_nuc] =  nucTo;
    CurrentNodeCodonSequence[NodeIndex][site_codon] = codonTo[1];

    //lparam->codonstatespace->GetCodonFromDNA(nucposTo[3], nucposTo[4], nucposTo[5]);

    delete [] nucposFrom;
    delete [] nucposTo;
    delete [] codonFrom;
    delete [] codonTo;

}






void TreeSimulator::ComputeRecursiveSimulation(Link* from)
{

    int verbose = lparam->verbose;

    int FromNodeIndex = from->GetNode()->GetIndex();

    ////
    // IF is ROOT
    ////
    if(from->isRoot())
    {

        if (verbose)
        {
            cerr << "CRS1\n";
        }
        submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);


        double rate = 0.0;
        double time = 0.0;
        double blength = 0.0;



        rate = submatrix->GetTotalSubRate(FromNodeIndex);
        //cerr << "Sub " << rate << "\n";
        if (lparam->model == "FMutSelSimu")
        {

            rate = submatrix->GetTotalMutRate(FromNodeIndex);
            //cerr << "Mut " << rate << "\n";
        }

        if (verbose)
        {
            cerr << "TreeSimulator::ComputeRecursiveSimulation1.1\n";
        }
        ////
        // IF is ROOT
        ////

//        rootBranchEvoStats->MutRate[0][0] = submatrix->GetTotalMutRate(FromNodeIndex);
//        rootBranchEvoStats->MutRate[0][1] = submatrix->GetTotalMutRateNonSyn(FromNodeIndex);
//        rootBranchEvoStats->MutRate[0][2] = submatrix->GetTotalMutRateSyn(FromNodeIndex);
//        rootBranchEvoStats->SubRate[0][0] = submatrix->GetTotalSubRate(FromNodeIndex);
//        rootBranchEvoStats->SubRate[0][1] = submatrix->GetTotalSubRateNonSyn(FromNodeIndex);
//        rootBranchEvoStats->SubRate[0][2] = submatrix->GetTotalSubRateSyn(FromNodeIndex);


        time = (lparam->rnd->sExpo()) /rate;

        blength = lparam->rootlength;

        if (verbose)
        {
            cerr << "TreeSimulator::ComputeRecursiveSimulation2\n";
            cerr << time << "\n" << blength << "\n";
            cerr << submatrix->GetTotalMutRate(FromNodeIndex) << "\n" << submatrix->GetTotalMutRateNonSyn(FromNodeIndex) << "\n" << submatrix->GetTotalMutRateSyn(FromNodeIndex) << "\n";
            cerr << submatrix->GetTotalSubRate(FromNodeIndex) << "\n" << submatrix->GetTotalSubRateNonSyn(FromNodeIndex) << "\n" << submatrix->GetTotalSubRateSyn(FromNodeIndex) << "\n";
        }

        double IntervalLength = blength/10;
        int interval = 0;
        GetSampleAncestralCodonSequence(FromNodeIndex, interval);
        interval ++;
        while (time < blength)
        {



            double u = lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);
            int site_nuc = 0;
            int nucTo = 0;
            double testcummul = submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);
//            if (lparam->model == "FMutSelSimu") {
//
//                testcummul = submatrix->GetMutRate(FromNodeIndex,site_nuc,nucTo);
//
//            }
            int k = 0 ;
            while (testcummul < u )
            {

                k++;
                site_nuc =(int) k/4 ;
                nucTo = k%4;
                //cerr << FromNodeIndex << " " << site_nuc << " " << nucTo << " " <<testcummul << " " << u << "\n";
                if(site_nuc >= lparam->Nsite_nuc)
                {
                    cerr << "error " << site_nuc << " > " << lparam->Nsite_nuc << "testcummul " << testcummul <<"\n";
                    exit(0);
                }

                testcummul +=  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);

//                if (lparam->model == "FMutSelSimu") {
//
//                    testcummul = submatrix->GetMutRate(FromNodeIndex,site_nuc,nucTo);
//
//                }



            }



            //cerr << submatrix[FromNodeIndex][site_nuc][nucTo] << "\n";
            int site_codon = int(site_nuc/3);
            if (verbose)
            {
                cerr << "TreeSimulator::ComputeRecursiveSimulation3\n";
            }
            submatrix->ComputePartialRates(FromNodeIndex, site_codon,CurrentNodeNucSequence);
            RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo) ;
            if (verbose)
            {
                cerr << "TreeSimulator::ComputeRecursiveSimulation4\n";
            }
            submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, site_codon,CurrentNodeNucSequence);
            if (verbose)
            {
                cerr << "TreeSimulator::ComputeRecursiveSimulation5\n";
            }
            //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, -1,CurrentNodeNucSequence);
            rate = submatrix->GetTotalSubRate(FromNodeIndex);
            //cerr << "Sub " << rate << "\n";
            if (lparam->model == "FMutSelSimu")
            {

                rate = submatrix->GetTotalMutRate(FromNodeIndex);
                //cerr << "Mut " << rate << "\n";
            }

            time += (lparam->rnd->sExpo()) /rate;

            if (verbose)
            {
                cerr << "TreeSimulator::ComputeRecursiveSimulation5.0 " << interval << "\n" ;
            }
            if (time>IntervalLength*interval)
            {
                GetSampleAncestralCodonSequence(FromNodeIndex,interval);
                interval ++;
            }

        }

        if (verbose)
        {
            cerr << "TreeSimulator::ComputeRecursiveSimulation5.1\n";
        }

//        rootBranchEvoStats->MutRate[1][0] = submatrix->GetTotalMutRate(FromNodeIndex);
//        rootBranchEvoStats->MutRate[1][1] = submatrix->GetTotalMutRateNonSyn(FromNodeIndex);
//        rootBranchEvoStats->MutRate[1][2] = submatrix->GetTotalMutRateSyn(FromNodeIndex);
//        rootBranchEvoStats->SubRate[1][0] = submatrix->GetTotalSubRate(FromNodeIndex);
//        rootBranchEvoStats->SubRate[1][1] = submatrix->GetTotalSubRateNonSyn(FromNodeIndex);
//        rootBranchEvoStats->SubRate[1][2] = submatrix->GetTotalSubRateSyn(FromNodeIndex);

        for (Link* link = from->Next(); link != from; link = link->Next())
        {
            int OutNodeIndex = link->Out()->GetNode()->GetIndex();
            submatrix->transfertTotalRate(FromNodeIndex,OutNodeIndex);


            for (int site_codon = 0 ; site_codon < lparam->Nsite_codon ; site_codon++)
            {
                CurrentNodeCodonSequence[OutNodeIndex][site_codon] = CurrentNodeCodonSequence[FromNodeIndex][site_codon];
                for (int j = 0; j < 3; j ++)
                {
                    CurrentNodeNucSequence[OutNodeIndex][site_codon*3+j] = CurrentNodeNucSequence[FromNodeIndex][site_codon*3+j];
                    submatrix->transfertNodeMatrix(FromNodeIndex,OutNodeIndex,site_codon*3+j);


                }
            }

            ComputeRecursiveSimulation(link->Out());

        }
    }
    else if(!from->isLeaf())
    {

        double rate = 0.0;
        double time = 0.0;
        double blength = 0.0;

        //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
        rate = submatrix->GetTotalSubRate(FromNodeIndex);

        //cerr << "Sub " << rate << "\n";
        if (lparam->model == "FMutSelSimu")
        {
            rate = submatrix->GetTotalMutRate(FromNodeIndex);
            //cerr << "Mut " << rate << "\n";

        }

        time = (lparam->rnd->sExpo()) /rate;

        blength = atof(from->GetBranch()->GetName().c_str());
        if (verbose)
        {
            cerr << "CRS5\n";
        }

        while (time < blength)
        {
            double u = lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);
            int site_nuc = 0;
            int nucTo = 0;
            double testcummul = submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);
            int k = 0 ;
            while (testcummul < u )
            {
                k++;
                site_nuc =(int) k/4 ;
                nucTo = k%4;
                if(site_nuc >= lparam->Nsite_nuc)
                {
                    cerr << "error " << site_nuc << " >= " << lparam->Nsite_nuc << "\n";
                    exit(0);
                }
                testcummul +=  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);

            }
            if (verbose)
            {
                cerr << "CRS6\n";
            }

            int site_codon = int(site_nuc/3);

            if (verbose)
            {
                cerr << "CRS7\n";
            }
            submatrix->ComputePartialRates(FromNodeIndex, site_codon,CurrentNodeNucSequence);
            RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo) ;

            if (verbose)
            {
                cerr << "CRS8\n";
            }

            submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, site_codon,CurrentNodeNucSequence);
            //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
            rate = submatrix->GetTotalSubRate(FromNodeIndex);

            if (lparam->model == "FMutSelSimu")
            {
                rate = submatrix->GetTotalMutRate(FromNodeIndex);

            }

            time += (lparam->rnd->sExpo()) /rate;
            //cerr << rate << "\n";

        }


        for (Link* link = from->Next(); link != from; link = link->Next())
        {
            int OutNodeIndex = link->Out()->GetNode()->GetIndex();
            submatrix->transfertTotalRate(FromNodeIndex,OutNodeIndex);
            for (int site_codon = 0 ; site_codon < lparam->Nsite_codon ; site_codon++)
            {
                CurrentNodeCodonSequence[OutNodeIndex][site_codon] = CurrentNodeCodonSequence[FromNodeIndex][site_codon];
                for (int j = 0; j < 3; j ++)
                {
                    CurrentNodeNucSequence[OutNodeIndex][site_codon*3+j] =  CurrentNodeNucSequence[FromNodeIndex][site_codon*3+j];
                    submatrix->transfertNodeMatrix(FromNodeIndex,OutNodeIndex,site_codon*3+j);
                }
            }

            ComputeRecursiveSimulation(link->Out());

        }
    }
    else if (from->isLeaf())
    {


        double rate = 0.0;
        double time = 0.0;
        double blength = 0.0;
        submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
        rate = submatrix->GetTotalSubRate(FromNodeIndex);
        //cerr << rate << "\n";
        if (lparam->model == "FMutSelSimu")
        {
            rate = submatrix->GetTotalMutRate(FromNodeIndex);
        }

        time = (lparam->rnd->sExpo()) /rate;

        blength = atof(from->GetBranch()->GetName().c_str());


        while (time < blength)
        {

            double u = lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);

            int site_nuc = 0;

            int nucTo = 0;
            double testcummul =  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);
            int k = 0 ;
            while (testcummul < u )
            {
                k++;
                site_nuc =(int) k/4 ;
                nucTo = k%4;
                if(site_nuc >= lparam->Nsite_nuc)
                {
                    cerr << "error " << site_nuc << " >= " << lparam->Nsite_nuc << "\n";
                    exit(0);
                }
                testcummul +=  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);

            }
            //cerr << submatrix[FromNodeIndex][site_nuc][nucTo] << "\n";
            int site_codon = int(site_nuc/3);

            submatrix->ComputePartialRates(FromNodeIndex, site_codon,CurrentNodeNucSequence);
            RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo) ;
            submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,site_codon,CurrentNodeNucSequence);
            //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
            rate = submatrix->GetTotalSubRate(FromNodeIndex);

            if (lparam->model == "FMutSelSimu")
            {
                rate = submatrix->GetTotalMutRate(FromNodeIndex);


            }

            time += (lparam->rnd->sExpo()) /rate;
            //cerr << rate << "\n";
        }


        for (int site_codon = 0 ; site_codon < lparam->Nsite_codon ; site_codon++)
        {
            int state = lparam->codondata->GetState(FromNodeIndex,site_codon);
            if (state != unknown)
            {
                CurrentLeafNodeCodonSequences[FromNodeIndex][site_codon] = CurrentNodeCodonSequence[FromNodeIndex][site_codon];
                for (int j = 0; j < 3; j ++)
                {
                    CurrentLeafNodeNucSequence[FromNodeIndex][site_codon*3+j] = CurrentNodeNucSequence[FromNodeIndex][site_codon*3+j];
                }
            }
            else
            {
                CurrentLeafNodeCodonSequences[FromNodeIndex][site_codon] = -1;
                for (int j = 0; j < 3; j ++)
                {
                    CurrentLeafNodeNucSequence[FromNodeIndex][site_codon*3+j] = -1;
                }
            }
        }

    }




}

void TreeSimulator::GetSampleAncestralCodonSequence(int FromNodeIndex,int interval)
{

    for(int site_codon = 0 ; site_codon < lparam->Nsite_codon ; site_codon++)
    {
        CurrentAncestralCodonSequence[interval][0][site_codon] = CurrentNodeCodonSequence[FromNodeIndex][site_codon];

    }

}


void TreeSimulator::WriteAncestral()
{
    // Write stationary
    ofstream ancestral_os ((lparam->output + ".ancestral").c_str(),APPEND);
    for (int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon ++)
    {
        ancestral_os << lparam->codonstatespace->GetState(CurrentNodeCodonSequence[0][site_codon]);
    }
    ancestral_os << "\n";
    ancestral_os.close();
}





