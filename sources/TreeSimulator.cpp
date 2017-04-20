#include "TreeSimulator.h"

TreeSimulator::TreeSimulator(LocalParameters* lparam, SiteInterSubMatrix* submatrix, AncestralSequence* ancestralseq)
{
    this->lparam = lparam;
    this->submatrix = submatrix;
    this->ancestralseq = ancestralseq;
    this->evohist = new EvolHistStatistics(this->lparam);
    this->ancestralSequenceEvohist = new EvolHistStatistics(this->lparam);

        CurrentNodeCodonSequence = new int*[lparam->refTree->GetNnode()];
        CurrentNodeNucSequence = new int*[lparam->refTree->GetNnode()];

        for (int node = 0 ; node < lparam->refTree->GetNnode(); node++) {
            CurrentNodeCodonSequence[node]= new int[lparam->Nsite_codon];
            CurrentNodeNucSequence[node]= new int[lparam->Nsite_codon*3];
        }
        CurrentLeafNodeCodonSequences= new int * [lparam->Ntaxa];
        CurrentLeafNodeNucSequence= new int * [lparam->Ntaxa];
        for (int taxa = 0 ; taxa < lparam->Ntaxa ; taxa ++ ) {
            CurrentLeafNodeCodonSequences[taxa] = new int [lparam->Nsite_codon];
            CurrentLeafNodeNucSequence[taxa] = new int [lparam->Nsite_nuc];
        }


}

TreeSimulator::~TreeSimulator()
{
    //dtor
}

void TreeSimulator::resetSimulator(){
// reset node sequeneces
    for (int node = 0 ; node < lparam->refTree->GetNnode(); node++) {
            for(int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon++){
                CurrentNodeCodonSequence[node][site_codon] = -2 ;
                for(int j = 0 ; j < 3; j ++){
                    CurrentNodeNucSequence[node][site_codon*3+j] = -2;
                }
            }
    }
    // reset nodeleaf sequences
    for (int taxa = 0 ; taxa < lparam->Ntaxa ; taxa ++ ) {
        for(int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon++){
            CurrentLeafNodeCodonSequences[taxa][site_codon] = -2 ;
            for(int j = 0 ; j < 3; j ++){
                CurrentLeafNodeNucSequence[taxa][site_codon*3+j] = -2;
            }
        }
    }
}

void TreeSimulator::GetNewSimulatedCodonAlignment(){

    // reset  submatrix
    submatrix->resetSubMatrix();

    evohist->resetMappingStat();
    ancestralSequenceEvohist->resetMappingStat();

    resetSimulator();

    ancestralseq->GetNewStationaryCodonSequence();


    // copy stationary sequence to node sequence
    int NodeIndex = lparam->refTree->GetRoot()->GetNode()->GetIndex();

    //exit(0);




    for (int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon++){
        CurrentNodeCodonSequence[NodeIndex][site_codon] = ancestralseq->GetCurrentAncestralCodonSequence(site_codon);
        for(int j = 0 ; j < 3; j++) {
            CurrentNodeNucSequence[NodeIndex][site_codon*3+j] = ancestralseq->GetCurrentAncestralNucSequence(site_codon*3+j);
        }
    }

    //lauchn recursive simulation on a phylogenetic tree
    ComputeRecursiveSimulation(lparam->refTree->GetRoot());

    //register mappingstats
    lparam->mappingstats.clear();
    lparam->mappingstats.shrink_to_fit();
    ancestralSequenceEvohist->GetMapAncStats();

    evohist->GetMapStats();


}


void TreeSimulator::RegisterSubTreeSim(int NodeIndex, int site_nuc, int nucTo) {

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

        {
            int codon_count = 0;
            for (int site_codon_i = site_codon_start; site_codon_i < site_codon_end ; site_codon_i++) {

                int site_nuc_start = (site_codon_i * 3);  // site_codon to site_nuc
                for (int codonPos = 0; codonPos < 3; codonPos++) {

                    if (site_codon_i >= 0 && site_codon_i < lparam->Nsite_codon) {


                        nucposFrom[codon_count*3+codonPos] = CurrentNodeNucSequence[NodeIndex][site_codon_i*3+codonPos];
                        //nucposTo[codon_count*3+codonPos] = CurrentNodeNucSequence[NodeIndex][site_codon_i*3+codonPos];
                        nucposTo[codon_count*3+codonPos] = nucposFrom[codon_count*3+codonPos];

                    } else {

                        nucposFrom[codon_count*3+codonPos] = -1;
                        nucposTo[codon_count*3+codonPos] = -1;

                    }



                }
                if (site_codon_i >= 0 && site_codon_i < lparam->Nsite_codon) {
                    codonTo[codon_count] = lparam->codonstatespace->GetCodonFromDNA(nucposTo[codon_count*3], nucposTo[codon_count*3+1], nucposTo[codon_count*3+2]);
                    codonFrom[codon_count] = codonTo[codon_count];
                } else {
                    codonTo[codon_count] = -1;
                    codonFrom[codon_count] = -1;
                }
                codon_count++;
            }
        }

        int site_nuc_To = site_nuc%3+3;
        nucposTo[site_nuc_To] = nucTo;
        codonTo[1] = lparam->codonstatespace->GetCodonFromDNA(nucposTo[3], nucposTo[4], nucposTo[5]);





        if(lparam->codonstatespace->CheckStop(nucposTo[3], nucposTo[4], nucposTo[5])){
                cerr << "error while registring a substitution\n";
                cerr << nucTo << "\n";
                cerr << nucposTo[3] << " " <<  nucposTo[4] << " " <<   nucposTo[5] << " " << codonTo[1]<<"\n";
                cerr << nucposFrom[3] << " " <<  nucposFrom[4] << " " <<   nucposFrom[5] << " "<< codonFrom[1] <<"\n";
                exit(0);
        }


        // Mapping statitics
         //Statistics on the ancestral sequence
         if(NodeIndex == lparam->refTree->GetRoot()->GetNode()->GetIndex()) {
                ancestralSequenceEvohist->Nsub++;

                if (CodonPos == 0){
                    ancestralSequenceEvohist->gtnr_stat[3][nucFrom][nucTo]++;
                    ancestralSequenceEvohist->gtnr_stat[CodonPos][nucFrom][nucTo]++;
                    ancestralSequenceEvohist->dinuc_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    ancestralSequenceEvohist->dinuc_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    if (site_nuc>0){
                        ancestralSequenceEvohist->dinuc_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        ancestralSequenceEvohist->dinuc_stat[2][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    }

                    if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1])){
                        ancestralSequenceEvohist->Nsynsub++;
                        ancestralSequenceEvohist->ssNsynsub[site_codon]++;

                        ancestralSequenceEvohist->gtnrSyn_stat[3][nucFrom][nucTo]++;
                        ancestralSequenceEvohist->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;

                        ancestralSequenceEvohist->dinucSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        ancestralSequenceEvohist->dinucSyn_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        if (site_nuc>0){
                            ancestralSequenceEvohist->dinucSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                            ancestralSequenceEvohist->dinucSyn_stat[2][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        }

                    } else {
                        ancestralSequenceEvohist->ssNnsynsub[site_codon]++;
                        ancestralSequenceEvohist->gtnrNSyn_stat[3][nucFrom][nucTo]++;
                        ancestralSequenceEvohist->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;

                        ancestralSequenceEvohist->dinucNSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        ancestralSequenceEvohist->dinucNSyn_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        if (site_nuc>0){
                            ancestralSequenceEvohist->dinucNSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                            ancestralSequenceEvohist->dinucNSyn_stat[2][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        }

                    }

                } else if(CodonPos == 1){

                    ancestralSequenceEvohist->gtnr_stat[3][nucFrom][nucTo]++;
                    ancestralSequenceEvohist->gtnr_stat[CodonPos][nucFrom][nucTo]++;

                    ancestralSequenceEvohist->dinuc_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    ancestralSequenceEvohist->dinuc_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

                    ancestralSequenceEvohist->dinuc_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    ancestralSequenceEvohist->dinuc_stat[0][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //



                    if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1])){
                        ancestralSequenceEvohist->Nsynsub++;
                        ancestralSequenceEvohist->ssNsynsub[site_codon]++;

                        ancestralSequenceEvohist->gtnrSyn_stat[3][nucFrom][nucTo]++;
                        ancestralSequenceEvohist->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;

                        ancestralSequenceEvohist->dinucSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        ancestralSequenceEvohist->dinucSyn_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

                        ancestralSequenceEvohist->dinucSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        ancestralSequenceEvohist->dinucSyn_stat[0][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                    } else {
                        ancestralSequenceEvohist->ssNnsynsub[site_codon]++;
                        ancestralSequenceEvohist->gtnrNSyn_stat[3][nucFrom][nucTo]++;
                        ancestralSequenceEvohist->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;

                        ancestralSequenceEvohist->dinucNSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        ancestralSequenceEvohist->dinucNSyn_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

                        ancestralSequenceEvohist->dinucNSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        ancestralSequenceEvohist->dinucNSyn_stat[0][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //


                    }

                } else if (CodonPos == 2){

                    ancestralSequenceEvohist->gtnr_stat[3][nucFrom][nucTo]++;
                    ancestralSequenceEvohist->gtnr_stat[CodonPos][nucFrom][nucTo]++;

                    ancestralSequenceEvohist->dinuc_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    ancestralSequenceEvohist->dinuc_stat[1][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                    if (site_codon < lparam->Nsite_codon-2){
                        ancestralSequenceEvohist->dinuc_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        ancestralSequenceEvohist->dinuc_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    }


                    if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1])){
                        ancestralSequenceEvohist->Nsynsub++;
                        ancestralSequenceEvohist->ssNsynsub[site_codon]++;

                        ancestralSequenceEvohist->gtnrSyn_stat[3][nucFrom][nucTo]++;
                        ancestralSequenceEvohist->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;

                        ancestralSequenceEvohist->dinucSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        ancestralSequenceEvohist->dinucSyn_stat[1][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                        if (site_codon < lparam->Nsite_codon-2){
                            ancestralSequenceEvohist->dinucSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                            ancestralSequenceEvohist->dinucSyn_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        }

                    } else {
                        ancestralSequenceEvohist->ssNnsynsub[site_codon]++;
                        ancestralSequenceEvohist->gtnrNSyn_stat[3][nucFrom][nucTo]++;
                        ancestralSequenceEvohist->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;

                        ancestralSequenceEvohist->dinucNSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        ancestralSequenceEvohist->dinucNSyn_stat[1][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                        if (site_codon < lparam->Nsite_codon-2){
                            ancestralSequenceEvohist->dinucNSyn_stat[3][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                            ancestralSequenceEvohist->dinucNSyn_stat[CodonPos][ancestralSequenceEvohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][ancestralSequenceEvohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        }

                    }


                }



        } else {

            // Non-root event
            evohist->Nsub++;

            if (CodonPos == 0){

                    evohist->gtnr_stat[3][nucFrom][nucTo]++;
                    evohist->gtnr_stat[CodonPos][nucFrom][nucTo]++;

                    evohist->dinuc_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    evohist->dinuc_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    if (site_nuc>0){
                        evohist->dinuc_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        evohist->dinuc_stat[2][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    }


                    if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1])){
                        evohist->Nsynsub++;
                        evohist->ssNsynsub[site_codon]++;

                        evohist->gtnrSyn_stat[3][nucFrom][nucTo]++;
                        evohist->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;

                        evohist->dinucSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        evohist->dinucSyn_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        if (site_nuc>0){
                            evohist->dinucSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                            evohist->dinucSyn_stat[2][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        }

                    } else {
                        evohist->ssNnsynsub[site_codon]++;
                        evohist->gtnrNSyn_stat[3][nucFrom][nucTo]++;
                        evohist->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;

                        evohist->dinucNSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        evohist->dinucNSyn_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        if (site_nuc>0){
                            evohist->dinucNSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                            evohist->dinucNSyn_stat[2][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        }

                    }

                } else if(CodonPos == 1){

                    evohist->gtnr_stat[3][nucFrom][nucTo]++;
                    evohist->gtnr_stat[CodonPos][nucFrom][nucTo]++;

                    evohist->dinuc_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    evohist->dinuc_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

                    evohist->dinuc_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    evohist->dinuc_stat[0][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //



                    if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1])){
                        evohist->Nsynsub++;
                        evohist->ssNsynsub[site_codon]++;

                        evohist->gtnrSyn_stat[3][nucFrom][nucTo]++;
                        evohist->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;

                        evohist->dinucSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        evohist->dinucSyn_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

                        evohist->dinucSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        evohist->dinucSyn_stat[0][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                    } else {
                        evohist->ssNnsynsub[site_codon]++;
                        evohist->gtnrNSyn_stat[3][nucFrom][nucTo]++;
                        evohist->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;

                        evohist->dinucNSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        evohist->dinucNSyn_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //

                        evohist->dinucNSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        evohist->dinucNSyn_stat[0][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //


                    }

                } else if (CodonPos == 2){

                    evohist->gtnr_stat[3][nucFrom][nucTo]++;
                    evohist->gtnr_stat[CodonPos][nucFrom][nucTo]++;

                    evohist->dinuc_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                    evohist->dinuc_stat[1][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                    if (site_codon < lparam->Nsite_codon-2){
                        evohist->dinuc_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        evohist->dinuc_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                    }


                    if(lparam->codonstatespace->Synonymous(codonFrom[1],codonTo[1])){
                        evohist->Nsynsub++;
                        evohist->ssNsynsub[site_codon]++;

                        evohist->gtnrSyn_stat[3][nucFrom][nucTo]++;
                        evohist->gtnrSyn_stat[CodonPos][nucFrom][nucTo]++;

                        evohist->dinucSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        evohist->dinucSyn_stat[1][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                        if (site_codon < lparam->Nsite_codon-2){
                            evohist->dinucSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                            evohist->dinucSyn_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                        }

                    } else {
                        evohist->ssNnsynsub[site_codon]++;
                        evohist->gtnrNSyn_stat[3][nucFrom][nucTo]++;
                        evohist->gtnrNSyn_stat[CodonPos][nucFrom][nucTo]++;

                        evohist->dinucNSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //
                        evohist->dinucNSyn_stat[1][evohist->GetDinucContext(nucposFrom[site_nuc_To-1],nucposFrom[site_nuc_To])][evohist->GetDinucContext(nucposTo[site_nuc_To-1],nucposTo[site_nuc_To])]++; //

                        if (site_codon < lparam->Nsite_codon-2){
                            evohist->dinucNSyn_stat[3][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
                            evohist->dinucNSyn_stat[CodonPos][evohist->GetDinucContext(nucposFrom[site_nuc_To],nucposFrom[site_nuc_To+1])][evohist->GetDinucContext(nucposTo[site_nuc_To],nucposTo[site_nuc_To+1])]++; //
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






void TreeSimulator::ComputeRecursiveSimulation(Link* from){


    int FromNodeIndex = from->GetNode()->GetIndex();

    if(from->isRoot()){


        submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);


        double rate = 0.0;
        double time = 0.0;
        double blength = 0.0;


        rate = submatrix->GetTotalSubRate(FromNodeIndex);
        //cerr << "Sub " << rate << "\n";
        if (lparam->model == "FMutSelSimu") {

            rate = submatrix->GetTotalMutRate(FromNodeIndex);
            //cerr << "Mut " << rate << "\n";
        }

        time = (lparam->rnd->sExpo()) /rate;

        blength = lparam->rootlength;

        while (time < blength) {
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
            while (testcummul < u ) {

                k++;
                site_nuc =(int) k/4 ;
                nucTo = k%4;
                //cerr << FromNodeIndex << " " << site_nuc << " " << nucTo << " " <<testcummul << " " << u << "\n";
                if(site_nuc >= lparam->Nsite_nuc) {cerr << "error " << site_nuc << " > " << lparam->Nsite_nuc << "testcummul " << testcummul <<"\n"; exit(0);}

                testcummul +=  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);

//                if (lparam->model == "FMutSelSimu") {
//
//                    testcummul = submatrix->GetMutRate(FromNodeIndex,site_nuc,nucTo);
//
//                }

            }


            //cerr << submatrix[FromNodeIndex][site_nuc][nucTo] << "\n";
            int site_codon = int(site_nuc/3);
            RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo) ;
            submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, site_codon,CurrentNodeNucSequence);
            //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, -1,CurrentNodeNucSequence);
            rate = submatrix->GetTotalSubRate(FromNodeIndex);
            //cerr << "Sub " << rate << "\n";
            if (lparam->model == "FMutSelSimu") {

                rate = submatrix->GetTotalMutRate(FromNodeIndex);
                //cerr << "Mut " << rate << "\n";
            }

            time += (lparam->rnd->sExpo()) /rate;


        }

	    for (Link* link = from->Next(); link != from; link = link->Next()) {
            int OutNodeIndex = link->Out()->GetNode()->GetIndex();
            submatrix->transfertTotalRate(FromNodeIndex,OutNodeIndex);


            for (int site_codon = 0 ; site_codon < lparam->Nsite_codon ; site_codon++){
                CurrentNodeCodonSequence[OutNodeIndex][site_codon] = CurrentNodeCodonSequence[FromNodeIndex][site_codon];
                for (int j = 0; j < 3; j ++){
                    CurrentNodeNucSequence[OutNodeIndex][site_codon*3+j] = CurrentNodeNucSequence[FromNodeIndex][site_codon*3+j];
                    submatrix->transfertNodeMatrix(FromNodeIndex,OutNodeIndex,site_codon*3+j);


                }
            }

            ComputeRecursiveSimulation(link->Out());

        }
   } else if(!from->isLeaf()) {

        double rate = 0.0;
        double time = 0.0;
        double blength = 0.0;

        //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
        rate = submatrix->GetTotalSubRate(FromNodeIndex);

        //cerr << "Sub " << rate << "\n";
        if (lparam->model == "FMutSelSimu") {
            rate = submatrix->GetTotalMutRate(FromNodeIndex);
            //cerr << "Mut " << rate << "\n";

        }

        time = (lparam->rnd->sExpo()) /rate;

        blength = atof(from->GetBranch()->GetName().c_str());

        while (time < blength) {
            double u = lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);
            int site_nuc = 0;
            int nucTo = 0;
            double testcummul = submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);
            int k = 0 ;
            while (testcummul < u ) {
                k++;
                site_nuc =(int) k/4 ;
                nucTo = k%4;
                if(site_nuc >= lparam->Nsite_nuc) {cerr << "error " << site_nuc << " >= " << lparam->Nsite_nuc << "\n"; exit(0);}
                testcummul +=  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);

            }
            int site_codon = int(site_nuc/3);

            RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo) ;

            submatrix->UpdateSubMatrixTreeSim(FromNodeIndex, site_codon,CurrentNodeNucSequence);
            //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
            rate = submatrix->GetTotalSubRate(FromNodeIndex);

            if (lparam->model == "FMutSelSimu") {
                rate = submatrix->GetTotalMutRate(FromNodeIndex);


            }

            time += (lparam->rnd->sExpo()) /rate;
            //cerr << rate << "\n";

        }


        for (Link* link = from->Next(); link != from; link = link->Next()) {
            int OutNodeIndex = link->Out()->GetNode()->GetIndex();
            submatrix->transfertTotalRate(FromNodeIndex,OutNodeIndex);
            for (int site_codon = 0 ; site_codon < lparam->Nsite_codon ; site_codon++){
                CurrentNodeCodonSequence[OutNodeIndex][site_codon] = CurrentNodeCodonSequence[FromNodeIndex][site_codon];
                for (int j = 0; j < 3; j ++){
                    CurrentNodeNucSequence[OutNodeIndex][site_codon*3+j] =  CurrentNodeNucSequence[FromNodeIndex][site_codon*3+j];
                    submatrix->transfertNodeMatrix(FromNodeIndex,OutNodeIndex,site_codon*3+j);
                }
            }

            ComputeRecursiveSimulation(link->Out());

        }
     } else if (from->isLeaf()){


        double rate = 0.0;
        double time = 0.0;
        double blength = 0.0;
        submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
        rate = submatrix->GetTotalSubRate(FromNodeIndex);
        //cerr << rate << "\n";
        if (lparam->model == "FMutSelSimu") {
            rate = submatrix->GetTotalMutRate(FromNodeIndex);
        }

        time = (lparam->rnd->sExpo()) /rate;

        blength = atof(from->GetBranch()->GetName().c_str());


        while (time < blength) {

			double u = lparam->rnd->Uniform() * submatrix->GetTotalSubRate(FromNodeIndex);

           	int site_nuc = 0;

            int nucTo = 0;
			double testcummul =  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);
			int k = 0 ;
			while (testcummul < u ) {
                k++;
                site_nuc =(int) k/4 ;
                nucTo = k%4;
                if(site_nuc >= lparam->Nsite_nuc) {cerr << "error " << site_nuc << " >= " << lparam->Nsite_nuc << "\n"; exit(0);}
                testcummul +=  submatrix->GetSubRate(FromNodeIndex,site_nuc,nucTo);

			}
			//cerr << submatrix[FromNodeIndex][site_nuc][nucTo] << "\n";
			int site_codon = int(site_nuc/3);


            RegisterSubTreeSim(FromNodeIndex, site_nuc, nucTo) ;
	        submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,site_codon,CurrentNodeNucSequence);
            //submatrix->UpdateSubMatrixTreeSim(FromNodeIndex,-1,CurrentNodeNucSequence);
            rate = submatrix->GetTotalSubRate(FromNodeIndex);

            if (lparam->model == "FMutSelSimu") {
                rate = submatrix->GetTotalMutRate(FromNodeIndex);


            }

            time += (lparam->rnd->sExpo()) /rate;
            //cerr << rate << "\n";
        }


        for (int site_codon = 0 ; site_codon < lparam->Nsite_codon ; site_codon++){
            int state = lparam->codondata->GetState(FromNodeIndex,site_codon);
            if (state != unknown) {
                CurrentLeafNodeCodonSequences[FromNodeIndex][site_codon] = CurrentNodeCodonSequence[FromNodeIndex][site_codon];
                for (int j = 0; j < 3; j ++){
                    CurrentLeafNodeNucSequence[FromNodeIndex][site_codon*3+j] = CurrentNodeNucSequence[FromNodeIndex][site_codon*3+j];
                }
            } else {
                CurrentLeafNodeCodonSequences[FromNodeIndex][site_codon] = -1;
                for (int j = 0; j < 3; j ++){
                    CurrentLeafNodeNucSequence[FromNodeIndex][site_codon*3+j] = -1;
                }
            }
        }

    }




}




void TreeSimulator::WriteAncestral(){
     // Write stationary
     ofstream ancestral_os ((lparam->output + ".ancestral").c_str(),APPEND);
     for (int site_codon = 0 ; site_codon < lparam->Nsite_codon; site_codon ++) {
        ancestral_os << lparam->codonstatespace->GetState(CurrentNodeCodonSequence[0][site_codon]);
     }
     ancestral_os << "\n";
     ancestral_os.close();
}





