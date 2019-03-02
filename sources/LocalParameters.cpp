/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/
#include "LocalParameters.h"


LocalParameters::LocalParameters(GlobalParameters* gparam)
{

    // input info
    //this->gparam = gparam;
    this->localcontrolfile = gparam->localcontrolfile;
    this->verbose = gparam->verbose;
    this->output = gparam->output;
    this->model = gparam->model;
    this->distance = gparam->distance;
    this->transformation = gparam->transformation;
    this->TOOSMALL = gparam->TOOSMALL;
    this->TOOLARGE = gparam->TOOLARGE;
    this->TOOLARGENEGATIVE = gparam->TOOLARGENEGATIVE;
    this->Ntaxa = gparam->Ntaxa;
    this->Nsite_codon = gparam->Nsite_codon;
    this->Nsite_nuc = 3 * this->Nsite_codon;


    if ((int) gparam->listSpecies.size() > 0 && this->Ntaxa >0 && (int) gparam->listSpecies.size() == this->Ntaxa)
    {
        this->listSpecies  = new string[this->Ntaxa];
        for (int i = 0 ; i < this->Ntaxa; i ++)
        {
            this->listSpecies[i] = gparam->listSpecies[i];
        }
    }


    this->NSummaries = gparam->NSummaries;
    this->NParam = gparam->NParam;
    this->NEvoStats = gparam->NEvoStats;
    this->NSiteSpecificEvoStats = gparam->NSiteSpecificEvoStats;

    if(gparam->verbose)
    {
        cerr << "LocalParam1\n";
    }

    this->listParam = new string[this->NParam];
    for (int param_i = 0; param_i < this->NParam; param_i++)
    {
        this->listParam[param_i] = gparam->listParam[param_i];
    }

    if(gparam->verbose)
    {
        cerr << "LocalParam2\n";
    }


    this->listSummaries = new string[this->NSummaries];
    for (int summary_i = 0; summary_i < this->NSummaries; summary_i++)
    {
        this->listSummaries[summary_i] = gparam->listSummaries[summary_i];
    }

    if(gparam->verbose)
    {
        cerr << "LocalParam3\n";
    }


    this->listEvoStats = new string[this->NEvoStats];
    for (int EvoStats_i = 0; EvoStats_i < this->NEvoStats; EvoStats_i++)
    {
        this->listEvoStats[EvoStats_i] = gparam->listEvoStats[EvoStats_i];
    }

    if(gparam->verbose)
    {
        cerr << "LocalParam4\n";
    }


    this->listSiteSpecificEvoStats = new string[this->NSiteSpecificEvoStats];
    for (int EvoStats_i = 0; EvoStats_i < this->NSiteSpecificEvoStats; EvoStats_i++)
    {
        this->listSiteSpecificEvoStats[EvoStats_i] = gparam->listSiteSpecificEvoStats[EvoStats_i];
    }

    if(gparam->verbose)
    {
        cerr << "LocalParam5\n";
    }

    this->NusedEvoStats = gparam->NusedEvoStats;
    this->NusedSiteSpecificEvoStats = gparam->NusedSiteSpecificEvoStats;
    this->NusedEvoAncStats = gparam->NusedEvoAncStats;
    this->NusedParam = gparam->NusedParam;
    this->NusedSummaries = gparam->NusedSummaries;
    this->NusedAncSummaries = gparam->NusedAncSummaries;
    this->NusedAccessorySummaries = gparam->NusedAccessorySummaries;
    this->Ngenes = gparam->Ngenes;

    this->mapUsedParam.insert(gparam->mapUsedParam.begin(),gparam->mapUsedParam.end());
    this->mapUsedSummaries.insert(gparam->mapUsedSummaries.begin(),gparam->mapUsedSummaries.end());
    this->mapUsedAncSummaries.insert(gparam->mapUsedAncSummaries.begin(),gparam->mapUsedAncSummaries.end());
    this->mapUsedAccessorySummaries.insert(gparam->mapUsedAccessorySummaries.begin(),gparam->mapUsedAccessorySummaries.end());
    this->mapUsedEvoStats.insert(gparam->mapUsedEvoStats.begin(),gparam->mapUsedEvoStats.end());
    this->mapUsedSiteSpecificEvoStats.insert(gparam->mapUsedSiteSpecificEvoStats.begin(),gparam->mapUsedSiteSpecificEvoStats.end());
    this->mapUsedEvoAncStats.insert(gparam->mapUsedEvoAncStats.begin(),gparam->mapUsedEvoAncStats.end());

    this->randomseed = -1;
    this->rnd = new Random(randomseed);
    this->MutationNormFactor = this->MutationNormFactor1 = this->MutationNormFactor2 =  1;
    this->getrate = false;
    this->getrate1 = false; 
    this->getrate2 = false; 

    //default parameters values
    this->percentFromOutGroup = 0.5;
    this->omega = 1.0;
    this->lambda_TBL = 1.0;
    this->lambda_omega = 1.0;
    this->lambda_CpG = 1.0;
    this->lambda_tvCpG = 1.0;
    this->lambda_tvTpA = 1.0;
    this->lambda_tstvCpG = 1.0; 
    this->lambda_tstvTpA = 1.0;
    this->lambda_TpA = 1.0;
    this->lambda_R = 0.5;
    this->fitCpG = 0.5;
    this->fitTpA = 0.5;
    this->fitGC = 0.5;
    this->AAadj = new double[20];
    for (int i = 0 ; i < 20 ; i ++)
    {
        AAadj[i] = 1; 
    }
    this->Aadj = 1.0;
    this->Cadj = 1.0;
    this->Dadj = 1.0;
    this->Eadj = 1.0;
    this->Fadj = 1.0;
    this->Gadj = 1.0;
    this->Hadj = 1.0;
    this->Iadj = 1.0;
    this->Kadj = 1.0;
    this->Ladj = 1.0;
    this->Madj = 1.0;
    this->Nadj = 1.0;
    this->Padj = 1.0;
    this->Qadj = 1.0;
    this->Radj = 1.0;
    this->Sadj = 1.0;
    this->Tadj = 1.0;
    this->Vadj = 1.0;
    this->Wadj = 1.0;
    this->Yadj = 1.0;
    
    
    //parameter switch
    this->fixNsite = 0; 
    this->fixfitGC = 1;
    this->fixfitCpG = 1;
    this->fixfitTpA = 1;
    this->fixlambda_TBL = 1;
    this->fixlambda_omega = 1;
    this->fixlambda_CpG = 1;
    this->fixlambda_tvCpG = 1;
    this->fixlambda_tvTpA = 1;
    this->fixlambda_tstvCpG = 1;
    this->fixlambda_tstvTpA = 1;
    this->fixlambda_TpA = 1;
    this->fixlambda_R = 1;
    this->fixAAadj = 1; 
    this->fixCODONadj = 1;
    this->fixgtr = 1;
    this->fixgtr2 = -1;
    this->fixgtr1 = -1;
    this->fixgtnr = 1;
    this->fixhky = 1;
    this->fixkappa = 1;
    this->fixstat = 1;
    this->fixtr = 1;
    this->fixts = 1;
    this->fixrr = 1;
    this->fixss = 1;
    this->rooted = 0;
    this->fixroot = -1;
    this->randomseed = -1;
    this->rootlength = 10.0;
    this->Ninterval = (int) this->rootlength;
    this->sampleAncSeq = 0;
    this->Nsite_codon = 1000;
    this->Nsite_nuc = 3*this->Nsite_codon;
    this->startPoint = 1;
    this->endPoint = 1;
    this->everyPoint = 1;
    this->Nrep = 1;
    this->taxa_gtr1_a = "";
    this->taxa_gtr1_b = "";
    this->taxa_gtr2_a = "";
    this->taxa_gtr2_b = "";
    this->gtr1NodeIndex = -1;
    this->gtr2NodeIndex = -1;
    this->taxa_a = "";
    this->taxa_b = "";
    this->tofasta = 0;
    this->isdata = false;
    this->iscodon = false;
    this->lambda_TBL_prior = "log2Unif";
    this->lambda_CpG_prior = "log10Unif";
    this->lambda_TpA_prior = "log10Unif";
    this->lambda_omega_prior = "log2Unif";

    this->nucrrnr = new double*[this->Nnucp];
    this->nucrrnr1 = new double*[this->Nnucp];
    this->nucrrnr2 = new double*[this->Nnucp];
    this->gtnr = new double*[this->Nnucp];
    this->gtnr1 = new double*[this->Nnucp];
    this->gtnr2 = new double*[this->Nnucp];
    for (int i = 0; i < this->Nnucp; i++)
    {
        this->gtnr[i] = new double[this->Nnucp];
        this->nucrrnr[i] = new double[this->Nnucp];
        this->gtnr1[i] = new double[this->Nnucp];
        this->nucrrnr1[i] = new double[this->Nnucp];
        this->gtnr2[i] = new double[this->Nnucp];
        this->nucrrnr2[i] = new double[this->Nnucp];
    }
    for (int i = 0; i < this->Nnucp; i++)
    {
        for (int j = 0; j < this->Nnucp; j++)
        {
            this->gtnr[i][j] = 0.0;
            this->nucrrnr[i][j] = 0.0;
            this->gtnr1[i][j] = 0.0;
            this->nucrrnr1[i][j] = 0.0;
            this->gtnr2[i][j] = 0.0;
            this->nucrrnr2[i][j] = 0.0;
        }
    }
    this->nucp = new double [this->Nnucp];
    this->nucp1 = new double [this->Nnucp];
    this->nucp2 = new double [this->Nnucp];
    for (int i = 0; i < this->Nnucp; i++)
    {
        this->nucp[i] = 0.0;
        this->nucp1[i] = 0.0;
        this->nucp2[i] = 0.0;
    }
    this->nucrr = new double [this->Nnucrr];
    this->nucrr1 = new double [this->Nnucrr];
    this->nucrr2 = new double [this->Nnucrr];
    for (int i= 0; i < this->Nnucrr; i++)
    {
        this->nucrr[i] = 0.0;
        this->nucrr1[i] = 0.0;
        this->nucrr2[i] = 0.0;
    }

    readLocalInstructions();

    if (this->isdata)
    {
        cerr << "alignment found\n";
        this->dnadata = new FileSequenceAlignment((data).c_str(), 0, 0);
        if (this->iscodon)
        {
            if(code == "Universal" )
            {
                cerr << "Universal\n";
                this->codonstatespace =  new CodonStateSpace(Universal);
                this->codondata  = new CodonSequenceAlignment(dnadata, true, Universal);
                this->taxonset = this->codondata->GetTaxonSet();
            }
            else if (code == "MtMam")
            {
                cerr << "MtMam\n";
                this->codonstatespace =  new CodonStateSpace(MtMam);
                this->codondata  = new CodonSequenceAlignment(dnadata, true, MtMam);
                this->taxonset = this->codondata->GetTaxonSet();
            }
            else if (code == "MtInv")
            {
                cerr << "MtInv\n";
                this->codonstatespace =  new CodonStateSpace(MtInv);
                this->codondata  = new CodonSequenceAlignment(dnadata, true, MtInv);
                this->taxonset = this->codondata->GetTaxonSet();
            }
            else
            {
                cerr << "wrong genetic code\n";
            }
        }
    }
    else
    {
        cerr << "no alignment found\n";
        if (this->Ntaxa > 0 && this->Nsite_nuc > 0 && this->listSpecies >0)
        {
            this->dnadata = new FileSequenceAlignment(this->Ntaxa, this->Nsite_nuc,this->listSpecies);

            if (this->iscodon)
            {
                if(code == "Universal" )
                {
                    cerr << "Universal\n";
                    this->codonstatespace =  new CodonStateSpace(Universal);
                    this->codondata  = new CodonSequenceAlignment(dnadata, true, Universal);
                    this->taxonset = this->codondata->GetTaxonSet();
                }
                else if (code == "MtMam")
                {
                    cerr << "MtMam\n";
                    this->codonstatespace =  new CodonStateSpace(MtMam);
                    this->codondata  = new CodonSequenceAlignment(dnadata, true, MtMam);
                    this->taxonset = this->codondata->GetTaxonSet();
                }
                else if (code == "MtInv")
                {
                    cerr << "MtInv\n";
                    this->codonstatespace =  new CodonStateSpace(MtInv);
                    this->codondata  = new CodonSequenceAlignment(dnadata, true, MtInv);
                    this->taxonset = this->codondata->GetTaxonSet();
                }
                else
                {
                    cerr << "wrong genetic code\n";
                }
            }
        }


    }
    this->Nsite_codon  = this->codondata->GetNsite();
    this->Nsite_nuc    = this->Nsite_codon * 3;
    this->Ntaxa        = this->codondata->GetNtaxa();
    this->Nstate_codon = this->codondata->GetNstate();
    if (this->Nsite_codon >4999)
    {
        cerr << "number of sites too large >=5000: " << this->Nsite_codon  <<"\n";
        exit(0);
    }

    this->codonprofile = new double [this->Nstate_codon];
    for (int i = 0 ; i < this->Nstate_codon; i++)
    {
        this->codonprofile[i] = 0.0;
    }

    this->ssaaprofiles = new double*[Nsite_codon];
    for(int i =0 ; i < this->Nsite_codon; i++)
    {
        this->ssaaprofiles[i] = new double [Nstate_aa];
    }
    for(int i = 0; i < this->Nsite_codon; i++)
    {
        for(int j =0 ; j < this->Nstate_aa; j++)
        {
            this->ssaaprofiles[i][j] = 0.0;
        }
    }

    this->alloc = new int [Nsite_codon];
    for(int i =0 ; i < this->Nsite_codon; i++)
    {
        this->alloc[i] = -1;
    }

    this->newlink = new Link();
    this->newnext = new Link();
    this->branchToInGroup = new Branch();
    this->branchToOutGroup = new Branch();
    this->newnode = new Node();

    cerr << "Nsite_codon : " << this->Nsite_codon << "\n";
    cerr << "Nsite_nuc : " << this->Nsite_nuc << "\n";
    cerr << "Ntaxa : " <<  this->Ntaxa << "\n";
    cerr << "Nstate_codon : " << this->Nstate_codon << "\n";

    int taxa_a_val = this->taxonset->GetTaxonIndexWithIncompleteName(this->taxa_a);
    int taxa_b_val = this->taxonset->GetTaxonIndexWithIncompleteName(this->taxa_b);
    
    if (taxa_a_val == -1 && taxa_b_val == -1)
    {
        cerr << "both taxon used for determining root position are absent from the sequence aligment\n"; 
        cerr << taxa_a << " " << taxa_b << "\n"; 
        exit(0);
    } 
    else if (taxa_a_val == -1)
    {
        cerr << "the taxon used for determining root position is absent from the sequence aligment\n";
        cerr << taxa_a << "\n"; 
        exit(0);
    }
    else if (taxa_b_val == -1)
    {
        cerr << "the taxon used for determining root position is absent from the sequence aligment\n";
        cerr << taxa_b << "\n"; 
        exit(0);
    }

    cerr << "Rooted at " << this->taxa_a << " " << this->taxa_b  << "\n";
}





LocalParameters::~LocalParameters()
{
    //dtor
}

void LocalParameters::GetGTR1()
{
    double sum =  0.0 ;
    for (int nuc1 = 0 ; nuc1 < this->Nnucp; nuc1++)
    {
        sum += this->nucp1[nuc1];
    }

    for (int nuc1 = 0 ; nuc1 < this->Nnucp; nuc1++)
    {
        this->nucp1[nuc1]/=sum;
    }

    //nucrrnr1[0][0]; //AA
    nucrrnr1[0][1] =  nucrr1[0]; //AC
    nucrrnr1[0][2] =  nucrr1[1]; //AG
    nucrrnr1[0][3] =  nucrr1[2]; //AT
    nucrrnr1[1][0] =  nucrr1[0]; //CA
    //nucrrnr1[1][1]; //CG
    nucrrnr1[1][2] =  nucrr1[3]; //CG
    nucrrnr1[1][3] =  nucrr1[4]; //CT
    nucrrnr1[2][0] =  nucrr1[1]; //GA
    nucrrnr1[2][1] =  nucrr1[3]; //GC
    //nucrrnr1[2][2]; //GG
    nucrrnr1[2][3] =  nucrr1[5]; //GT
    nucrrnr1[3][0] =  nucrr1[2]; //TA
    nucrrnr1[3][1] =  nucrr1[4]; //TC
    nucrrnr1[3][2] =  nucrr1[5]; //TG
    //nucrrnr1[3][3]; //TT


    sum = 0.0 ;
    for (int nuc1=0; nuc1<4-1; nuc1++)
    {
        for (int nuc2=nuc1+1; nuc2<4; nuc2++)
        {
            sum += this->nucrrnr1[nuc1][nuc2];
        }
    }
    for (int nuc1=0; nuc1<4; nuc1++)
    {
        for (int nuc2=0; nuc2<4; nuc2++)
        {
            this->nucrrnr1[nuc1][nuc2] /= sum;
        }
    }

    this->getrate1 = false;
    double norm = 0.0;
    for (int i=0; i<4-1; i++)
    {
        for (int j=i+1; j<4; j++)
        {
            //norm += nucp[i] * nucp[j] * nucrr[GetNucRRIndex(i,j)];
            norm += this->nucp1[i] * this->nucp1[j] * this->nucrrnr1[i][j];
            //cerr << this->nucp1[i] << " " << this->nucp1[j] << " " << this->nucrrnr1[i][j] << " " << norm << "\n";
        }
    }
    // 2 for the symetry of the matrix??, and 3 for the number of codon positons
    getrate1 = true;
    MutationNormFactor1 = 2 * (norm * 3);

    for (int i=0; i<4-1; i++)
    {
        for (int j=i+1; j<4; j++)
        {
            this->gtnr1[i][j] = (this->nucp1[j] * this->nucrrnr1[i][j]) / MutationNormFactor1 ;
        }
    }
}


void LocalParameters::GetGTR2()
{
    double sum =  0.0 ;
    for (int nuc1 = 0 ; nuc1 < this->Nnucp; nuc1++)
    {
        sum += this->nucp2[nuc1];
    }

    for (int nuc1 = 0 ; nuc1 < this->Nnucp; nuc1++)
    {
        this->nucp2[nuc1]/=sum;
    }
    //nucrrnr2[0][0]; //AA
    nucrrnr2[0][1] =  nucrr2[0]; //AC
    nucrrnr2[0][2] =  nucrr2[1]; //AG
    nucrrnr2[0][3] =  nucrr2[2]; //AT
    nucrrnr2[1][0] =  nucrr2[0]; //CA
    //nucrrnr2[1][1]; //CG
    nucrrnr2[1][2] =  nucrr2[3]; //CG
    nucrrnr2[1][3] =  nucrr2[4]; //CT
    nucrrnr2[2][0] =  nucrr2[1]; //GA
    nucrrnr2[2][1] =  nucrr2[3]; //GC
    //nucrrnr2[2][2]; //GG
    nucrrnr2[2][3] =  nucrr2[5]; //GT
    nucrrnr2[3][0] =  nucrr2[2]; //TA
    nucrrnr2[3][1] =  nucrr2[4]; //TC
    nucrrnr2[3][2] =  nucrr2[5]; //TG
    //nucrrnr2[3][3]; //TT

    sum = 0.0 ;
    for (int nuc1=0; nuc1<4-1; nuc1++)
    {
        for (int nuc2=nuc1+1; nuc2<4; nuc2++)
        {
            sum += this->nucrrnr2[nuc1][nuc2];
        }
    }

    for (int nuc1=0; nuc1<4; nuc1++)
    {
        for (int nuc2=0; nuc2<4; nuc2++)
        {
            this->nucrrnr2[nuc1][nuc2] /= sum;
        }
    }

    this->getrate1 = false;

    double norm = 0.0;
    for (int i=0; i<4-1; i++)
    {
        for (int j=i+1; j<4; j++)
        {
            //norm += nucp[i] * nucp[j] * nucrr[GetNucRRIndex(i,j)];
            norm += this->nucp2[i] * this->nucp2[j] * this->nucrrnr2[i][j];
            //cerr << this->nucp1[i] << " " << this->nucp1[j] << " " << this->nucrrnr1[i][j] << " " << norm << "\n";
        }
    }
    // 2 for the symetry of the matrix??, and 3 for the number of codon positons
    this->getrate2 = true;
    this->MutationNormFactor2 = 2 * (norm * 3);

    for (int i=0; i<4-1; i++)
    {
        for (int j=i+1; j<4; j++)
        {
            this->gtnr2[i][j] = (this->nucp2[j] * this->nucrrnr2[i][j]) / MutationNormFactor1 ;
        }
    }
}

void LocalParameters::readLocalInstructions()
{
    string line = localcontrolfile ;
    istringstream iss(line);
    string s;

    while(iss >> s)
    {
        if ( s =="-d")
        {
            iss >> s;
            this->data = s;
            cerr << "data " << data<< "\n";
            this->isdata = true;
        }
        else if (s =="-chain")
        {
            iss >> s;
            this->chain = s;
            cerr << "chain " << this->chain<< "\n";
        }
        else if (s == "-posterior")
        {
            iss >> s;
            this->posteriorfile = s;
            cerr << "posterior " << this->posteriorfile << "\n";
        }
        else if (s == "-start")
        {
            iss >> s;
            this->startPoint = atoi(s.c_str());
            cerr << "start " << this->startPoint<< "\n";
        }
        else if (s == "-until")
        {
            iss >> s;
            this->endPoint = atoi(s.c_str());
            cerr << "until " << this->endPoint<< "\n";
        }
        else if (s == "-every")
        {
            iss >> s ;
            this->everyPoint = atoi(s.c_str());
            cerr << "every " << this->everyPoint << "\n";
        }
        else if (s == "-rep")
        {
            iss >> s;
            this->Nrep = atoi(s.c_str());
            cerr << "Nreplicate " << this->Nrep << "\n";
        }
        else if (s =="-iscodon")
        {
            this->iscodon = true;
        }
        else if (s =="-code")
        {
            iss >> s;
            this->code = s;
            cerr << "code " << this->code << "\n";
        }
        else if (s =="-root")
        {
            rooted = 1;
            iss >> s;
            this->taxa_a = s;
            iss >> s;
            this->taxa_b = s;
        }
        else if (s =="-fixroot")
        {
            fixroot = 1;
            iss >> s;
            this->taxa_a = s;
            iss >> s;
            this->taxa_b = s;
            iss >> s;
            this->percentFromOutGroup = atof(s.c_str());
//                    for (unsigned int param_i = 0 ; param_i < gparam->NParam; param_i++){
//                        if("root" == gparam->listParam[param_i]){
//                                gparam->listUsedParam.push_back(param_i);
//                                break;
//                        }
//                    }



            cerr << ((this->fixroot ==0)? "randomroot " : "fixroot ") << this->taxa_a   << " " << this->taxa_b  << " " << this->percentFromOutGroup << "\n";
        }
        else if (s =="-freeroot")
        {
            this->fixroot = 0;
            iss >> s;
            this->taxa_a = s;
            iss >> s;
            this->taxa_b = s;
            cerr << ((this->fixroot ==0)? "freeroot " : "fixroot ") << this->taxa_a  << " " << this->taxa_b << "\n";
        }
        else if (s =="-fitCpG")
        {
            iss >> s;
            this->fitCpG = atof(s.c_str());
            this->fixfitCpG = 1;
            cerr << "fixfitCpG " << this->fitCpG << "\n";
        }
        else if (s =="-fitGC")
        {
            iss >> s;
            this->fitGC = atof(s.c_str());
            this->fixfitGC = 1;
            cerr << "fixfitGC " << this->fitGC << "\n";
        }
        else if (s == "-freefitGC")
        {
            this->fixfitGC = 0;
            cerr << "freefitGC\n";
        }
        else if (s == "-freefitCpG")
        {
            this->fixfitCpG = 0;
            cerr << "free fitCpG\n";
        }
        else if (s == "-freefitTpA")
        {
            this->fixfitTpA = 0;
            cerr << "free fitTpA\n";
        }
        else if (s =="-lambdaTBL")
        {
            iss >> s;
            this->lambda_TBL = atof(s.c_str());
            this->fixlambda_TBL = 1;
            cerr << "fix TBL " << this->lambda_TBL << "\n";
        }
        else if (s == "-freelambdaTBL")
        {
            this->fixlambda_TBL = 0;
            cerr << "free TBL\n";
        }
        else if (s == "-priorlambdaTBL")
        {
            iss >> s;
            this->lambda_TBL_prior = s;
            cerr << "prior TBL " << this->lambda_TBL_prior << "\n";
        }
        else if(s == "-freelambdaR")
        {
            this->fixlambda_R = 0;
            cerr << "freelambda_R\n";
        }
        else if (s =="-lambdaR")
        {
            iss >> s;
            this->lambda_R = atof(s.c_str());
            this->fixlambda_R = 1;
            cerr << "fix lambda_R " << this->lambda_R << "\n";
        }
        else if (s =="-lambdaOmega")
        {
            iss >> s;
            this->lambda_omega = atof(s.c_str());
            this->fixlambda_omega = 1;
            cerr << "fix lambdaomega " << this->lambda_omega << "\n";
        }
        else if (s =="-omega")
        {
            iss >> s;
            this->omega = atof(s.c_str());
            this->fixomega = 1;
            cerr << "fix omega " << this->omega << "\n";
        }
        else if (s == "-freelambdaomega")
        {
            this->fixlambda_omega = 0;
            cerr << "freelambdaomega\n";
        }
        else if (s == "-priorlambdaomega")
        {
            iss >> s;
            this->lambda_omega_prior = s;
            cerr << "prior omega " << this->lambda_omega_prior << "\n";
        }
        else if (s =="-lambdaCpG" || s == "-fixlambdaCpG")
        {
            iss >> s;
            this->lambda_CpG = atof(s.c_str());
            this->fixlambda_CpG = 1;
            cerr << "fix lambda CpG " << this->lambda_CpG << "\n";
        }
        else if (s =="-fixgtr1")
        {
            this->fixgtr1 = 1;
            iss >> s ;
            this->taxa_gtr1_a = s;
            iss >> s ;
            this->taxa_gtr1_b = s;
            iss >> s ;
            this->taxa_gtr1_c = s;
            iss >> s;
            this->nucp1[0] = atof(s.c_str());
            iss >> s;
            this->nucp1[1] = atof(s.c_str());
            iss >> s;
            this->nucp1[2] = atof(s.c_str());
            iss >> s;
            this->nucp1[3] = atof(s.c_str());
            iss >> s;
            this->nucrr1[0] = atof(s.c_str());
            iss >> s;
            this->nucrr1[1] = atof(s.c_str());
            iss >> s;
            this->nucrr1[2] = atof(s.c_str());
            iss >> s;
            this->nucrr1[3] = atof(s.c_str());
            iss >> s;
            this->nucrr1[4] = atof(s.c_str());
            iss >> s;
            this->nucrr1[5] = atof(s.c_str());
            cerr << " fixgtr1 " << this->fixgtr1 << "\n";
            cerr << this->taxa_gtr1_a << "\t" << this->taxa_gtr1_b <<"\t" << this->taxa_gtr1_c << "\n";
            for (int v = 0 ; v < Nnuc ; v ++ )
            {
                cerr << this->nucp1[v] << "\t";
            }
            cerr << "\n";
            for (int v = 0 ; v < Nnucrr ; v ++ )
            {
                cerr << this->nucrr1[v] << "\t";
            }
            cerr << "\n";
        }
        else if (s =="-fixgtr2")
        {
            this->fixgtr2 = 1;
            iss >> s ;
            this->taxa_gtr2_a = s;
            iss >> s ;
            this->taxa_gtr2_b = s;
            iss >> s ;
            this->taxa_gtr2_c = s;
            iss >> s;
            this->nucp2[0] = atof(s.c_str());
            iss >> s;
            this->nucp2[1] = atof(s.c_str());
            iss >> s;
            this->nucp2[2] = atof(s.c_str());
            iss >> s;
            this->nucp2[3] = atof(s.c_str());
            iss >> s;
            this->nucrr2[0] = atof(s.c_str());
            iss >> s;
            this->nucrr2[1] = atof(s.c_str());
            iss >> s;
            this->nucrr2[2] = atof(s.c_str());
            iss >> s;
            this->nucrr2[3] = atof(s.c_str());
            iss >> s;
            this->nucrr2[4] = atof(s.c_str());
            iss >> s;
            this->nucrr2[5] = atof(s.c_str());


            cerr << " fixgtr2 " << this->fixgtr2 << "\n";
            cerr << this->taxa_gtr2_a << "\t" << this->taxa_gtr2_b << "\t" << this->taxa_gtr2_c << "\n";
            for (int v = 0 ; v < Nnuc ; v ++ )
            {
                cerr << this->nucp2[v] << "\t";
            }
            cerr << "\n";
            for (int v = 0 ; v < Nnucrr ; v ++ )
            {
                cerr << this->nucrr2[v] << "\t";
            }
            cerr << "\n";
        }
        else if (s =="-lambdaTpA" || s == "-fixlambdaTpA")
        {
            iss >> s;
            this->lambda_TpA = atof(s.c_str());
            this->fixlambda_TpA = 1;
            cerr << "fix lambda TpA " << this->fixlambda_TpA << "\n";
        }
        else if (s == "-lambdatstvCpG")
        {
            iss >> s;
            this->lambda_tstvCpG = atof(s.c_str());
            this->fixlambda_tstvCpG = 1;
            cerr << "fix lambda tstvCpG " << this->lambda_tstvCpG << "\n";
        }
        else if (s == "-lambdatstvTpA")
        {
            iss >> s;
            this->lambda_tstvTpA = atof(s.c_str());
            this->fixlambda_tstvTpA = 1;
            cerr << "fix lambda tstvTpA " << this->lambda_tstvTpA << "\n";
        }
        else if (s == "-freelambdatvCpG")
        {
            this->fixlambda_tvCpG = 0; 
            cerr << "freelambdatvCpG\n";
        }
        else if (s == "-freelambdatvTpA")
        {
            this->fixlambda_tvTpA = 0; 
            cerr << "freelambdatvTpA\n";
        }
        else if (s == "-freelambdatstvCpG")
        {
            this->fixlambda_tstvCpG = 0; 
            cerr << "freelambdatstvCpG\n";
        }
        else if (s == "-freelambdatstvTpA")
        {
            this->fixlambda_tstvTpA = 0; 
            cerr << "freelambdatstvCpG\n";
        }
        else if (s == "-freelambdaCpG" || s == "-freelambdatsCpG")
        {
            this->fixlambda_CpG = 0;
            cerr << "freelambdatsCpG\n";

        }
        else if (s == "-priorlambdaCpG")
        {
            iss >> s;
            this->lambda_CpG_prior = s;
            cerr << "prior CpG " << this->lambda_CpG_prior << "\n";

        }
        else if (s == "-freelambdaTpA" || s == "-freelambdatsTpA")
        {
            this->fixlambda_TpA = 0;
            cerr << "freelambdatsTpA\n";

        }
        else if (s == "-priorlambdaTpA")
        {
            iss >> s;
            this->lambda_TpA_prior = s;
            cerr << "prior TpA " << this->lambda_TpA_prior << "\n";

        }
        else if (s =="-freegtnr")
        {
            this->fixgtnr = 0;
            cerr << "fixgtnr\n";

        }
        else if (s =="-freehky")
        {
            this->fixhky = 0;
            cerr << "freehky\n";

        }
        else if (s =="-freeAAadj")
        {
            this->fixAAadj = 0;
            cerr << "freeAAadj\n";
        }
        else if (s =="-freegtr")
        {
            this->fixgtr = 0;
            cerr << "freegtr\n";

        }
        else if (s =="-freestat")
        {
            this->fixstat = 0;
            cerr << "freestat\n";

        }
        else if (s =="-freets")
        {
            this->fixts = 0;
            cerr << "freets\n";

        }
        else if (s =="-freetr")
        {
            this->fixtr = 0;
            cerr << "freetr\n";

        }
        else if (s =="-freekappa")
        {
            this->fixkappa = 0;
            cerr << "freekappa\n";

        }
        else if (s =="-freerr")
        {
            this->fixrr = 0;
            cerr << "freerr\n";

        }
        else if (s =="-rootlength")
        {
            iss >> s;
            this->rootlength = atof(s.c_str());
            cerr << "rootlength " << this->rootlength <<"\n";

        }      
        else if (s =="-tofasta" || s =="-tophylip")
        {
            this->tofasta = 1;

        }
        else if (s == "-randomfix")
        {
            this->randomseed = 12345;
            cerr << "randomseed " << this->randomseed << "\n";

        }
        else if (s == "-verbose")
        {
            iss >> s;
            this->verbose = 1;
            cerr << "verbose " << this->verbose << "\n";

        }
        else if (s == "-fixNsite")
        {
            iss >> s;
            this->fixNsite = 1;
            iss >> s;
            this->Nsite_codon = atoi(s.c_str());
            this->Nsite_nuc = this->Nsite_codon * 3;
            cerr << "fixNsite " << this->Nsite_codon << " " << this->Nsite_nuc << "\n";
        }
    }
}

void LocalParameters::SetRootLCA()
{
    outgroupLink = refTree->GetLCA(taxa_a,taxa_b);
    refTree->RootAt(outgroupLink);
    refTree->SetIndices();
}

void LocalParameters::SetTree()
{
    if (this->fixroot == 1 )
    {
        SetRootBetweenInAndOutGroup();
    }
    else if (this->fixroot == 0)
    {
        percentFromOutGroup = rnd->Uniform();
        SetRootBetweenInAndOutGroup();
    }
    else if (this->rooted == 1)
    {
        SetRootLCA();
    }
}

void LocalParameters::SetRootBetweenInAndOutGroup()
{

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup1\n";

    
        
    if (refTree->GetLCA(taxa_a,taxa_b)->isRoot())
    {    
        if (!refTree->CheckRootDegree())
        {
            cerr << "error: root should be of degree tree\n";
            cerr << "i.e. tree should have following format: (A,B,C) and not (A,(B,C));\n";
            exit(1);
        }
        
        if (verbose)
        {
            cerr << "LocalParameters::SetRootBetweenInAndOutGroup1 isRoot\n";
        }
        
        if (refTree->GetLCA(taxa_a,taxa_b)->Next()->Out()->GetNode()->GetName() != taxa_a && refTree->GetLCA(taxa_a,taxa_b)->Next()->Out()->GetNode()->GetName() != taxa_b)
        {
            outgroupLink = refTree->GetLCA(taxa_a,taxa_b)->Next()->Out();
        }
        else if (refTree->GetLCA(taxa_a,taxa_b)->Next()->Next()->Out()->GetNode()->GetName() != taxa_a && refTree->GetLCA(taxa_a,taxa_b)->Next()->Next()->Out()->GetNode()->GetName() != taxa_b)
        {
            outgroupLink = refTree->GetLCA(taxa_a,taxa_b)->Next()->Next()->Out();
        }
        else if (refTree->GetLCA(taxa_a,taxa_b)->Next()->Next()->Next()->Out()->GetNode()->GetName() != taxa_a && refTree->GetLCA(taxa_a,taxa_b)->Next()->Next()->Next()->Out()->GetNode()->GetName() != taxa_b)
        {
            outgroupLink = refTree->GetLCA(taxa_a,taxa_b)->Next()->Next()->Next()->Out();
        }

    }
    else
    {
        outgroupLink = refTree->GetLCA(taxa_a,taxa_b);
    }

        if (verbose)
            cerr << "LocalParameters::SetRootBetweenInAndOutGroup1 isRoot\n";
        //cerr << "You have to redefine the root position\n"; 
        //cerr << "The system will exit by now\n";
        //exit(0);

    branchLengthBetweenInAndOutGroup = atof(outgroupLink->GetBranch()->GetName().c_str());
    
    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup2\n";
    branchToOutGroup = outgroupLink->GetBranch();

    newnode->SetName("1");

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup3\n";
    newnode->SetIndex(refTree->GetNnode());

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup4\n";

    // set linkToIngroup

    outgroupLink->Out()->SetBranch(branchToInGroup);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup5\n";
    // set newlink

    newlink->SetBranch(branchToInGroup);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup6\n";
    newlink->SetNode(newnode);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup7\n";
    newlink->SetOut(outgroupLink->Out());

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup8\n";
    newlink->SetNext(newnext);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup9\n";
    // set newnext

    newnext->SetBranch(branchToOutGroup);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup10\n";
    newnext->SetNode(newnode);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup11\n";
    newnext->SetOut(outgroupLink);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup12\n";
    newnext->SetNext(newlink);

    // set outgroupLink

    outgroupLink->SetOut(newlink);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup13\n";
    refTree->RootAt(newnext);

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup14\n";
    refTree->SetIndices();

    double branchLengthToOutGroup = branchLengthBetweenInAndOutGroup * percentFromOutGroup;

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup15\n";
    branchToInGroup->SetName(std::to_string(branchLengthBetweenInAndOutGroup-branchLengthToOutGroup));

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup16\n";
    branchToOutGroup->SetName(std::to_string(branchLengthToOutGroup));

    if (verbose)
        cerr << "LocalParameters::SetRootBetweenInAndOutGroup17\n";

    ofstream os ((this->output+".inputparam").c_str(),APPEND);
    refTree->Print(os);
    os.close();
}

std::vector<double> LocalParameters::GetCurrentDistances()
{
    std::vector<double> cur_dist;
    double dist = 0.0 ;
    //squared Euclidian
    if (this->distance == "Euclidian")
    {
        for (int i_summary = 0 ; i_summary < this->NusedSummaries; i_summary++)
        {
            double sqdisc  = this->summariesRealData[i_summary] - summariesSimulatedData[i_summary];
            sqdisc *= sqdisc ;
            dist += sqdisc;
            cur_dist.push_back(sqdisc);
        }
        cur_dist.push_back(dist);
    }
    else if (this->distance == "normalized")
    {
        for (int i_summary = 0 ; i_summary < this->NusedSummaries; i_summary++)
        {
            double sqdisc  = (this->summariesRealData[i_summary] - summariesSimulatedData[i_summary])/this->summariesRealData[i_summary];
            sqdisc *= sqdisc ;
            dist += sqdisc;
            cur_dist.push_back(sqdisc);
        }
        cur_dist.push_back(dist);
    } else if


    (this->distance == "dist1")
    {
        for (int i_summary = 0 ; i_summary < this->NusedSummaries; i_summary++)
        {
            double sqdisc  = this->summariesRealData[i_summary] - summariesSimulatedData[i_summary];
            dist += sqdisc;
            cur_dist.push_back(sqdisc);
        }
        cur_dist.push_back(dist);
    }
    return cur_dist;
}

void LocalParameters::SetBranchesLengthsBetweenInAndOutGroup()
{
    double branchLengthToOutGroup = branchLengthBetweenInAndOutGroup * percentFromOutGroup;
    branchToInGroup->SetName(std::to_string(branchLengthBetweenInAndOutGroup-branchLengthToOutGroup));
    branchToOutGroup->SetName(std::to_string(branchLengthToOutGroup));
}


void LocalParameters::SetCurrentParametersFromPosterior(std::vector<std::vector<double>>posterior, int it)
{
    int k = 0; 
    string* arrParam = new string[this->NusedParam];
    for (int param_i = 0 ; param_i < this->NParam ; param_i++)
    {
        auto it_ = this->mapUsedParam.find(this->listParam[param_i]);
        if(it_ != this->mapUsedParam.end())
        {
            if(it_->second != -1)
            {
                arrParam[k] = it_->first;
                k++;
            }
        }
    }
    for (int param_i = 0 ; param_i < this->NusedParam ; param_i++)
    {
        if (arrParam[param_i] == "chainID")
        {
            this->MCMCpointID = static_cast<int> (posterior[it][param_i]);
            cerr << "chainID " << this->MCMCpointID << "\n";
        }
        else if (arrParam[param_i] == "root")
        {
            this->percentFromOutGroup = posterior[it][param_i];
            if (percentFromOutGroup >= 1.0 )
            {
                cerr << "percentFromOutGroup >= 1\n";
                this->percentFromOutGroup = 0.9999999999;
            }
            else if (percentFromOutGroup <= 0.0)
            {
                cerr << "percentFromOutGroup <= 0\n";
                this->percentFromOutGroup = 0.000000001;
            }
            this->SetBranchesLengthsBetweenInAndOutGroup();
        }
        else if(arrParam[param_i]  == "fitGC")
        {
            this->fitGC = posterior[it][param_i];
            if (this->fitGC >= 1.0)
            {
                this->fitGC = 0.9999999999;
            } 
            else if (this->fitGC <= 0.0)
            {
                this->fitGC = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "fitCpG")
        {
            this->fitCpG = posterior[it][param_i];
            if (fitCpG >= 1.0)
            {
                this->fitCpG = 0.9999999999;
            } 
            else if (fitCpG <= 0.0)
            {
                this->fitCpG = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "fitTpA")
        {
            this->fitTpA = posterior[it][param_i];
            if (fitTpA >= 1.0)
            {
                this->fitTpA = 0.9999999999;
            } 
            else if (fitTpA <= 0.0)
            {
                this->fitTpA = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "lambda_CpG")
        {
            this->lambda_CpG = posterior[it][param_i];
            if (this->lambda_CpG <= 0)
            {
                this->lambda_CpG = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "lambda_TpA")
        {
            this->lambda_TpA = posterior[it][param_i];
            if (this->lambda_TpA <= 0)
            {
                this->lambda_TpA = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "lambda_tvTpA")
        {
            this->lambda_tvTpA = posterior[it][param_i];
            if (this->lambda_tvTpA <= 0)
            {
                this->lambda_tvTpA = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "lambda_tvCpG")
        {
            this->lambda_tvCpG = posterior[it][param_i];
            if (this->lambda_tvCpG <= 0)
            {
                this->lambda_tvCpG = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "lambda_tstvTpA")
        {
            this->lambda_tstvTpA = posterior[it][param_i];
            if (this->lambda_tstvTpA <= 0)
            {
                this->lambda_tstvTpA = 0.0000000001;
            }
        }
        else if(arrParam[param_i]  == "lambda_tstvCpG")
        {
            this->lambda_tstvCpG = posterior[it][param_i];
            if (this->lambda_tstvCpG <= 0)
            {
                this->lambda_tstvCpG = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "lambda_TBL")
        {
            this->lambda_TBL = posterior[it][param_i];
            if (this->lambda_TBL <= 0)
            {
                this->lambda_TBL = 0.0000000001;
            }
//            for(int node = 0; node < refTree->GetNnode(); node++){
//                this->muBranch[node] = this->lambda_TBL;
//            }
        }
        else if (arrParam[param_i]  == "lambda_omega")
        {
            this->lambda_omega = posterior[it][param_i];
            if (this->lambda_omega <= 0)
            {
                this->lambda_omega = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "lambda_R")
        {
            this->lambda_R = posterior[it][param_i];
            if (this->lambda_R <= 0)
            {
                this->lambda_R = 0.0000000001;
            }
            if (this->lambda_R >= 1)
            {
                this->lambda_R = 0.9999999999;
            }


        }
        else if (arrParam[param_i]  == "nucsA")
        {
            this->nucp[0] = posterior[it][param_i];
            if (this->nucp[0] <= 0)
            {
                this->nucp[0] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucsC")
        {
            this->nucp[1] = posterior[it][param_i];
            if (this->nucp[1] <= 0)
            {
                this->nucp[1] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucsG")
        {
            this->nucp[2] = posterior[it][param_i];
            if (this->nucp[2] <= 0)
            {
                this->nucp[2] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucsT")
        {
            this->nucp[3] = posterior[it][param_i];
            if (this->nucp[3] <= 0)
            {
                this->nucp[3] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrAC")
        {
            this->nucrrnr[0][1] = posterior[it][param_i];
            if (this->nucrrnr[0][1] <= 0)
            {
                this->nucrrnr[0][1] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrAG")
        {
            this->nucrrnr[0][2] = posterior[it][param_i];
            if (this->nucrrnr[0][2] <= 0)
            {
                this->nucrrnr[0][2] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrAT")
        {
            this->nucrrnr[0][3] = posterior[it][param_i];
            if (this->nucrrnr[0][3] <= 0)
            {
                this->nucrrnr[0][3] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrCA")
        {
            this->nucrrnr[1][0] = posterior[it][param_i];
            if (this->nucrrnr[1][0] <= 0)
            {
                this->nucrrnr[1][0] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrCG")
        {
            this->nucrrnr[1][2] = posterior[it][param_i];
            if (this->nucrrnr[1][2] <= 0)
            {
                this->nucrrnr[1][2] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrCT")
        {
            this->nucrrnr[1][3] = posterior[it][param_i];
            if (this->nucrrnr[1][3] <= 0)
            {
                this->nucrrnr[1][3] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrGA")
        {
            this->nucrrnr[2][0] = posterior[it][param_i];
            if (this->nucrrnr[2][0] <= 0)
            {
                this->nucrrnr[2][0] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrGC")
        {
            this->nucrrnr[2][1] = posterior[it][param_i];
            if (this->nucrrnr[2][1] <= 0)
            {
                this->nucrrnr[2][1] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrGT")
        {
            this->nucrrnr[2][3] = posterior[it][param_i];
            if (this->nucrrnr[2][3] <= 0)
            {
                this->nucrrnr[2][3] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrTA")
        {
            this->nucrrnr[3][0] = posterior[it][param_i];
            if (this->nucrrnr[3][0] <= 0)
            {
                this->nucrrnr[3][0] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrTC")
        {
            this->nucrrnr[3][1] = posterior[it][param_i];
            if (this->nucrrnr[3][1] <= 0)
            {
                this->nucrrnr[3][1] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "nucrrTG")
        {
            this->nucrrnr[3][2] = posterior[it][param_i];
            if (this->nucrrnr[3][2] <= 0)
            {
                this->nucrrnr[3][2] = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Aadj")
        {
            this->Aadj = posterior[it][param_i];
            if (this->Aadj <= 0)
            {
                this->Aadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Cadj")
        {
            this->Cadj = posterior[it][param_i];
            if (this->Cadj <= 0)
            {
                this->Cadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Dadj")
        {
            this->Dadj = posterior[it][param_i];
            if (this->Dadj <= 0)
            {
                this->Dadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Eadj")
        {
            this->Eadj = posterior[it][param_i];
            if (this->Eadj <= 0)
            {
                this->Eadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Fadj")
        {
            this->Fadj = posterior[it][param_i];
            if (this->Fadj <= 0)
            {
                this->Fadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Gadj")
        {
            this->Gadj = posterior[it][param_i];
            if (this->Gadj <= 0)
            {
                this->Gadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Hadj")
        {
            this->Hadj = posterior[it][param_i];
            if (this->Hadj <= 0)
            {
                this->Hadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Iadj")
        {
            this->Iadj = posterior[it][param_i];
            if (this->Iadj <= 0)
            {
                this->Iadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Kadj")
        {
            this->Kadj = posterior[it][param_i];
            if (this->Kadj <= 0)
            {
                this->Kadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Ladj")
        {
            this->Ladj = posterior[it][param_i];
            if (this->Ladj <= 0)
            {
                this->Ladj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Madj")
        {
            this->Madj = posterior[it][param_i];
            if (this->Madj <= 0)
            {
                this->Madj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Nadj")
        {
            this->Nadj = posterior[it][param_i];
            if (this->Nadj <= 0)
            {
                this->Nadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Padj")
        {
            this->Padj = posterior[it][param_i];
            if (this->Padj <= 0)
            {
                this->Padj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Qadj")
        {
            this->Qadj = posterior[it][param_i];
            if (this->Qadj <= 0)
            {
                this->Qadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Radj")
        {
            this->Radj = posterior[it][param_i];
            if (this->Radj <= 0)
            {
                this->Radj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Sadj")
        {
            this->Sadj = posterior[it][param_i];
            if (this->Sadj <= 0)
            {
                this->Sadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Tadj")
        {
            this->Tadj = posterior[it][param_i];
            if (this->Tadj <= 0)
            {
                this->Tadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Vadj")
        {
            this->Vadj = posterior[it][param_i];
            if (this->Vadj <= 0)
            {
                this->Vadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Wadj")
        {
            this->Wadj = posterior[it][param_i];
            if (this->Wadj <= 0)
            {
                this->Wadj = 0.0000000001;
            }
        }
        else if (arrParam[param_i]  == "Yadj")
        {
            this->Yadj = posterior[it][param_i];
            if (this->Yadj <= 0)
            {
                this->Yadj = 0.0000000001;
            }
        }
    }
    delete [] arrParam;


    double sum =  0.0 ;
    for (int nuc1 = 0 ; nuc1 < this->Nnucp; nuc1++)
    {
        sum += this->nucp[nuc1];
    }
    for (int nuc1 = 0 ; nuc1 < this->Nnucp; nuc1++)
    {
        this->nucp[nuc1]/=sum;
    }

    sum = 0.0 ;
    for (int nuc1=0; nuc1<4-1; nuc1++)
    {
        for (int nuc2=nuc1+1; nuc2<4; nuc2++)
        {
            sum += this->nucrrnr[nuc1][nuc2];
        }
    }

    for (int nuc1=0; nuc1<4; nuc1++)
    {
        for (int nuc2=0; nuc2<4; nuc2++)
        {
            this->nucrrnr[nuc1][nuc2] /= sum;
        }
    }

    this->getrate = false;

    this->gtnr[0][1] = GetGTR(0,1); //ac
    this->gtnr[0][2] = GetGTR(0,2); //ag
    this->gtnr[0][3] = GetGTR(0,3); //at
    this->gtnr[1][0] = GetGTR(1,0); //ca
    //gtnr[1][1] = 0.0; //cc
    this->gtnr[1][2] = GetGTR(1,2); //cg
    this->gtnr[1][3] = GetGTR(1,3); //ct
    this->gtnr[2][0] = GetGTR(2,0); //ga
    this->gtnr[2][1] = GetGTR(2,1); //gc
    //gtnr[2][2] = 0.0; //gg
    this->gtnr[2][3] = GetGTR(2,3); //gt
    this->gtnr[3][0] = GetGTR(3,0); //ta
    this->gtnr[3][1] = GetGTR(3,1); //tc
    this->gtnr[3][2] = GetGTR(3,2); //tg
    //gtnr[3][3] = 0.0; //tt
}


std::vector<double> LocalParameters::GetCurrentWeights()
{
    return weights;
}

std::vector<double> LocalParameters::GetCurrentEvoStats()
{
    return evostats;
}

std::vector<double> LocalParameters::GetCurrentAncEvoStats()
{
    return ancevostats;
}

std::vector<double> LocalParameters::GetCurrentSiteSpecificEvoStats()
{
    return sitespecificevostats;
}

std::vector<double> LocalParameters::GetCurrentSummaries()
{
    return summariesSimulatedData;
}

std::vector<double> LocalParameters::GetCurrentAccessorySummaries()
{
    return accessorysummariesSimulatedData;
}


std::vector<double> LocalParameters::GetCurrentParameters()
{
    string* arrParam = new string[this->NusedParam];
    int k = 0 ; 
    for (int param_i = 0 ; param_i < this->NParam ; param_i++)
    {
        auto it = this->mapUsedParam.find(this->listParam[param_i]);
        if(it != this->mapUsedParam.end() )
        {
            if(it->second != -1)
            {
                arrParam[k] = it->first;
                k++;
            }
        }
    }

    std::vector<double> cur_param;
    for (int param_i = 0 ; param_i < this->NusedParam; param_i++)
    {
        if (arrParam[param_i] == "root")
        {
            cur_param.push_back(this->percentFromOutGroup);
        }
        else if(arrParam[param_i]  == "fitGC")
        {
            cur_param.push_back(this->fitGC);
        }
        else if(arrParam[param_i]  == "fitCpG")
        {
            cur_param.push_back(this->fitCpG);
        }
        else if(arrParam[param_i]  == "fitTpA")
        {
            cur_param.push_back(this->fitTpA);
        }
        else if(arrParam[param_i]  == "lambda_CpG")
        {
            cur_param.push_back(this->lambda_CpG);
        }
        else if(arrParam[param_i]  == "lambda_TpA")
        {
            cur_param.push_back(this->lambda_TpA);
        }
        else if(arrParam[param_i]  == "lambda_tvTpA")
        {
            cur_param.push_back(this->lambda_tvTpA);
        }
        else if(arrParam[param_i]  == "lambda_tvCpG")
        {
            cur_param.push_back(this->lambda_tvCpG);
        }
        else if(arrParam[param_i]  == "lambda_tstvTpA")
        {
            cur_param.push_back(this->lambda_tstvTpA);
        }
        else if(arrParam[param_i]  == "lambda_tstvCpG")
        {
            cur_param.push_back(this->lambda_tstvCpG);
        }
        else if (arrParam[param_i]  == "lambda_TBL")
        {
            cur_param.push_back(this->lambda_TBL);
        }
        else if (arrParam[param_i]  == "lambda_omega")
        {
            cur_param.push_back(this->lambda_omega);
        }
        else if (arrParam[param_i]  == "lambda_R")
        {
            cur_param.push_back(this->lambda_R);
        }
        else if (arrParam[param_i]  == "nucsA")
        {
            cur_param.push_back(this->nucp[0]);
        }
        else if (arrParam[param_i]  == "nucsC")
        {
            cur_param.push_back(this->nucp[1]);
        }
        else if (arrParam[param_i]  == "nucsG")
        {
            cur_param.push_back(this->nucp[2]);
        }
        else if (arrParam[param_i]  == "nucsT")
        {
            cur_param.push_back(this->nucp[3]);
        }
        else if (arrParam[param_i]  == "nucrrAC")
        {
            cur_param.push_back(this->nucrrnr[0][1]);
        }
        else if (arrParam[param_i]  == "nucrrAG")
        {
            cur_param.push_back(this->nucrrnr[0][2]);
        }
        else if (arrParam[param_i]  == "nucrrAT")
        {
            cur_param.push_back(this->nucrrnr[0][3]);
        }
        else if (arrParam[param_i]  == "nucrrCA")
        {
            cur_param.push_back(this->nucrrnr[1][0]);
        }
        else if (arrParam[param_i]  == "nucrrCG")
        {
            cur_param.push_back(this->nucrrnr[1][2]);
        }
        else if (arrParam[param_i]  == "nucrrCT")
        {
            cur_param.push_back(this->nucrrnr[1][3]);
        }
        else if (arrParam[param_i]  == "nucrrGA")
        {
            cur_param.push_back(this->nucrrnr[2][0]);
        }
        else if (arrParam[param_i]  == "nucrrGC")
        {
            cur_param.push_back(this->nucrrnr[2][1]);
        }
        else if (arrParam[param_i]  == "nucrrGT")
        {
            cur_param.push_back(this->nucrrnr[2][3]);
        }
        else if (arrParam[param_i]  == "nucrrTA")
        {
            cur_param.push_back(this->nucrrnr[3][0]);
        }
        else if (arrParam[param_i]  == "nucrrTC")
        {
            cur_param.push_back(this->nucrrnr[3][1]);
        }
        else if (arrParam[param_i]  == "nucrrTG")
        {
            cur_param.push_back(this->nucrrnr[3][2]);
        }
        else if (arrParam[param_i]  == "Aadj")
        {
            cur_param.push_back(this->Aadj);
        }
        else if (arrParam[param_i]  == "Cadj")
        {
            cur_param.push_back(this->Cadj);
        }
        else if (arrParam[param_i]  == "Dadj")
        {
            cur_param.push_back(this->Dadj);
        }
        else if (arrParam[param_i]  == "Eadj")
        {
            cur_param.push_back(this->Eadj);
        }
        else if (arrParam[param_i]  == "Fadj")
        {
            cur_param.push_back(this->Fadj);
        }
        else if (arrParam[param_i]  == "Gadj")
        {
            cur_param.push_back(this->Gadj);
        }
        else if (arrParam[param_i]  == "Hadj")
        {
            cur_param.push_back(this->Hadj);
        }
        else if (arrParam[param_i]  == "Iadj")
        {
            cur_param.push_back(this->Iadj);
        }
        else if (arrParam[param_i]  == "Kadj")
        {
            cur_param.push_back(this->Kadj);
        }
        else if (arrParam[param_i]  == "Ladj")
        {
            cur_param.push_back(this->Ladj);
        }
        else if (arrParam[param_i]  == "Madj")
        {
            cur_param.push_back(this->Madj);
        }
        else if (arrParam[param_i]  == "Nadj")
        {
            cur_param.push_back(this->Nadj);
        }
        else if (arrParam[param_i]  == "Padj")
        {
            cur_param.push_back(this->Padj);
        }
        else if (arrParam[param_i]  == "Qadj")
        {
            cur_param.push_back(this->Qadj);
        }
        else if (arrParam[param_i]  == "Radj")
        {
            cur_param.push_back(this->Radj);
        }
        else if (arrParam[param_i]  == "Sadj")
        {
            cur_param.push_back(this->Sadj);
        }
        else if (arrParam[param_i]  == "Tadj")
        {
            cur_param.push_back(this->Tadj);
        }
        else if (arrParam[param_i]  == "Vadj")
        {
            cur_param.push_back(this->Vadj);
        }
        else if (arrParam[param_i]  == "Wadj")
        {
            cur_param.push_back(this->Wadj);
        }
        else if (arrParam[param_i]  == "Yadj")
        {
            cur_param.push_back(this->Yadj);
        }
    }
    delete [] arrParam;
    return cur_param;
}



void LocalParameters::writeParam(ofstream& os)
{

    os << "\n";
    os << fixlambda_TBL << "\t" << lambda_TBL << "\n";
    os << fixlambda_omega << "\t" << lambda_omega << "\n";
    os << fixlambda_CpG << "\t" <<lambda_CpG << "\n";
    os << fixlambda_TpA << "\t" <<lambda_TpA << "\n";
    os << fixgtr << "\t" << fixgtr1 << "\t" << fixgtr2 << "\n";
    os << fixroot << "\n";
    os << rootlength << "\n";
    os << "#####\n";

    refTree->Print(os);
    os << "\n";

    for(int i = 0 ; i <Nnucp; i++)
    {

        os << nucp[i] << "\t";
    }
    os << "\n";

    for(int i = 0; i <Nnucrr; i++)
    {

        os << nucrr[i] << "\t";
    }
    os << "\n";
    os << omega << "\n";
    for(int i = 0; i <Nstate_codon; i++)
    {
        if (i < Nstate_codon-1)
        {
            os << codonprofile[i] << "\t";
        }
        else
        {
            os << codonprofile[i] << "\n";
        }
    }
    os << "\n";
    for(int j = 0; j < Nsite_codon; j++)
    {
        for(int i = 0; i <Nstate_aa; i++ )
        {
            if (i < Nstate_aa -1)
            {
                os << ssaaprofiles[j][i] << "\t";
            }
            else
            {
                os << ssaaprofiles[j][i] << "\n";
            }
        }
    }
    for(int i = 0; i < Nsite_codon; i++)
    {
        if (i < Nsite_codon -1)
        {
            os << alloc[i] << "\t";
        }
        else
        {
            os << alloc[i] << "\n";
        }
    }

    if (this->model == "FMutSelSimu" || this->model == "FMutSel0Simu")
    {

        os << "stationary distribution\n";
        double* stat = new double[Nstate_codon];
        double Z = 0.0 ;

        for (int state = 0; state < Nstate_codon; state++)
        {

            stat[state] =
                nucp[codonstatespace->GetCodonPosition(0, state)] *
                nucp[codonstatespace->GetCodonPosition(1, state)] *
                nucp[codonstatespace->GetCodonPosition(2, state)] *
                codonprofile[state]*
                ssaaprofiles[alloc[0]][codonstatespace->Translation(state)];

            Z += stat[state];
        }

        for (int state = 0; state < Nstate_codon; state++)
        {

            stat[state] /= Z;
        }

        for(int i = 0; i <Nstate_codon; i++)
        {

            if (i < Nstate_codon -1)
            {
                os << stat[i] << "\t";
            }
            else
            {
                os << stat[i] << "\n";
            }
        }
    }

}



void LocalParameters::readFMutSelCodeML()
{

    fstream is((this->chain + ".chain").c_str());
    if (!is)
    {
        cerr << "error: did not find " << this->chain << ".chain\n";
        exit(1);
    }


    this->refTree = new Tree(is);
    this->refTree->RegisterWith(taxonset, 0);

    for (int k=0; k<this->Nnucp; k++)
    {
        is >> this->nucp[k];
        //cerr << nucp[k] << "\t";
    }

    for (int k=0; k<this->Nnucrr; k++)
    {
        is >> this->nucrr[k];
        //cerr << nucrr[k] << "\t";
    }

    //cerr << "\n\n\n";


    //nucrrnr[0][0]; //AA
    this->nucrrnr[0][1] =  nucrr[0]; //AC
    this->nucrrnr[0][2] =  nucrr[1]; //AG
    this->nucrrnr[0][3] =  nucrr[2]; //AT
    this->nucrrnr[1][0] =  nucrr[0]; //CA
    //nucrrnr[1][1]; //CG
    this->nucrrnr[1][2] =  nucrr[3]; //CG
    this->nucrrnr[1][3] =  nucrr[4]; //CT
    this->nucrrnr[2][0] =  nucrr[1]; //GA
    this->nucrrnr[2][1] =  nucrr[3]; //GC
    //nucrrnr[2][2]; //GG
    this->nucrrnr[2][3] =  nucrr[5]; //GT
    this->nucrrnr[3][0] =  nucrr[2]; //TA
    this->nucrrnr[3][1] =  nucrr[4]; //TC
    this->nucrrnr[3][2] =  nucrr[5]; //TG
    //nucrrnr[3][3]; //TT

    double sum =  0.0 ;
    for (int nuc1 = 0 ; nuc1 < Nnucp; nuc1++)
    {
        sum += this->nucp[nuc1];
    }
    for (int nuc1 = 0 ; nuc1 < Nnucp; nuc1++)
    {
        this->nucp[nuc1]/=sum;
    }

    sum = 0.0 ;
    for (int nuc1=0; nuc1<4-1; nuc1++)
    {
        for (int nuc2=nuc1+1; nuc2<4; nuc2++)
        {
            sum += this->nucrrnr[nuc1][nuc2];
        }
    }

    for (int nuc1=0; nuc1<4; nuc1++)
    {
        for (int nuc2=0; nuc2<4; nuc2++)
        {
            this->nucrrnr[nuc1][nuc2] /= sum;
        }
    }



    is >> omega;


    for (int k=0; k<this->Nstate_codon; k++)
    {
        is >>  codonprofile[k];
    }

    for (int l=0; l<Nstate_aa; l++)
    {
        is >> ssaaprofiles[0][l];
    }

    for (int k=0; k<this->Nsite_codon; k++)
    {
        alloc[k] = 0;
    }


//    if (this->model == "FMutSel0Simu") {
//
//        double *aa_vec = new double[20];
//
//        for (int l = 0; l< this->Nstate_aa; l++){
//            aa_vec[l] = 0.0;
//        }
//
//
//        double Z = 0.0;
//        for (int k=0; k<this->Nstate_codon; k++){
//            double tmp;
//            is >>  tmp;
//            aa_vec[codonstatespace->Translation(k)] = tmp;
//            this->codonprofile[k] = (double) 1/this->Nstate_codon;
//            Z += tmp;
//            //cerr << this->codonprofile[k] << "\t";
//
//        }
//        for (int l=0; l<this->Nstate_aa; l++){
//            aa_vec[l] /= Z;
//        }
//
//        for (int k=0; k<this->Nsite_codon; k++){
//            for (int l=0; l<this->Nstate_aa; l++){
//                this->ssaaprofiles[k][l] = aa_vec[l];
//                //cerr << this->ssaaprofiles[k][l] << "\t";
//            }
//        }
//        delete [] aa_vec;
//
//    } else if (this->model == "FMutSelSimu") {
//
//        double Z = 0.0;
//        for (int k=0; k<this->Nstate_codon; k++){
//            is >>  this->codonprofile[k];
//            Z += this->codonprofile[k];
//
//        }
//        for (int l=0; l<this->Nstate_codon; l++){
//            this->codonprofile[l] /= Z;
//        }
//
//        for (int k=0; k<this->Nsite_codon; k++){
//            for (int l=0; l<this->Nstate_aa; l++){
//                this->ssaaprofiles[k][l] = (double) 1/this->Nstate_aa;
//            }
//        }
//    }


    for (int k=0; k<this->Nsite_codon; k++)
    {
        this->alloc[k] = 0;
    }

    this->getrate = false;

    this->gtnr[0][0] = 0.0; //AA
    this->gtnr[0][1] = GetGTR(0,1); //ac
    this->gtnr[0][2] = GetGTR(0,2); //ag
    this->gtnr[0][3] = GetGTR(0,3); //at
    this->gtnr[1][0] = GetGTR(1,0); //ca
    this->gtnr[1][1] = 0.0; //cc
    this->gtnr[1][2] = GetGTR(1,2); //cg
    this->gtnr[1][3] = GetGTR(1,3); //ct
    this->gtnr[2][0] = GetGTR(2,0); //ga
    this->gtnr[2][1] = GetGTR(2,1); //gc
    this->gtnr[2][2] = 0.0; //gg
    this->gtnr[2][3] = GetGTR(2,3); //gt
    this->gtnr[3][0] = GetGTR(3,0); //ta
    this->gtnr[3][1] = GetGTR(3,1); //tc
    this->gtnr[3][2] = GetGTR(3,2); //tg
    this->gtnr[3][3] = 0.0; //tt

    SetTreeStuff();

}

void LocalParameters::readChainCodonMutSelSBDP(int pt_i)
{
    this->MCMCpointID = pt_i;
//set parameters : posterior specific
    ifstream is((this->chain + ".chain").c_str());
    if (!is)
    {
        cerr << "error: did not find " << this->chain << ".chain\n";
        exit(1);
    }

    //cerr << this->chain;

    int j = 0;
    string tmp = "" ;
    while (j < pt_i)
    {

        is >> tmp; // tree
        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> tmp; //nucp
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> tmp; //nucrr
        }

        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >> tmp; //codon
        }
        is >> tmp; // omega
        is >> tmp; // kappa
        is >> tmp; // Ncomponents
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }
        for (int k=0; k<this->Nsite_codon; k++)
        {
            for (int l=0; l<this->Nstate_aa; l++)
            {
                is >> tmp; // ssprofiles
            }
        }
        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> tmp; // alloc
        }
        j++;

    }

    if(j == pt_i)
    {
        refTree = new Tree(is);
        if (verbose)
            cerr << "READCHAIN1\n"; 
        refTree->RegisterWith(taxonset, 0);
        if (verbose)
            cerr << "READCHAIN2\n"; 

        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> nucp[k];
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> nucrr[k];
        }

        //nucrrnr[0][0]; //AA
        nucrrnr[0][1] =  nucrr[0]; //AC
        nucrrnr[0][2] =  nucrr[1]; //AG
        nucrrnr[0][3] =  nucrr[2]; //AT
        nucrrnr[1][0] =  nucrr[0]; //CA
        //nucrrnr[1][1]; //CG
        nucrrnr[1][2] =  nucrr[3]; //CG
        nucrrnr[1][3] =  nucrr[4]; //CT
        nucrrnr[2][0] =  nucrr[1]; //GA
        nucrrnr[2][1] =  nucrr[3]; //GC
        //nucrrnr[2][2]; //GG
        nucrrnr[2][3] =  nucrr[5]; //GT
        nucrrnr[3][0] =  nucrr[2]; //TA
        nucrrnr[3][1] =  nucrr[4]; //TC
        nucrrnr[3][2] =  nucrr[5]; //TG
        //nucrrnr[3][3]; //TT

        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >>  codonprofile[k];
        }
        is >> omega;
        is >> tmp;
        is >> tmp;
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }

        for (int k=0; k<this->Nsite_codon; k++)
        {
            for (int l=0; l<Nstate_aa; l++)
            {
                is >> ssaaprofiles[k][l];
            }
        }

        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> alloc[k];
        }
    }
    is.close();

    this->getrate = false;

    gtnr[0][0] = 0.0; //AA
    gtnr[0][1] = GetGTR(0,1); //ac
    gtnr[0][2] = GetGTR(0,2); //ag
    gtnr[0][3] = GetGTR(0,3); //at
    gtnr[1][0] = GetGTR(1,0); //ca
    gtnr[1][1] = 0.0; //cc
    gtnr[1][2] = GetGTR(1,2); //cg
    gtnr[1][3] = GetGTR(1,3); //ct
    gtnr[2][0] = GetGTR(2,0); //ga
    gtnr[2][1] = GetGTR(2,1); //gc
    gtnr[2][2] = 0.0; //gg
    gtnr[2][3] = GetGTR(2,3); //gt
    gtnr[3][0] = GetGTR(3,0); //ta
    gtnr[3][1] = GetGTR(3,1); //tc
    gtnr[3][2] = GetGTR(3,2); //tg
    gtnr[3][3] = 0.0; //tt

    SetTreeStuff();

    //cerr << "Nnode : " << refTree->GetNnode() << "\n";
}



void LocalParameters::readChainCodonMutSelSBDP()
{

//set parameters : posterior specific
    ifstream is((this->chain + ".chain").c_str());
    if (!is)
    {
        cerr << "error: did not find " << this->chain << ".chain\n";
        exit(1);
    }

    //cerr << this->chain;

    int j = 0;
    string tmp = "" ;
    while (j < this->startPoint)
    {

        is >> tmp; // tree
        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> tmp; //nucp
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> tmp; //nucrr
        }

        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >> tmp; //codon
        }
        is >> tmp; // omega
        is >> tmp; // kappa
        is >> tmp; // Ncomponents
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }
        for (int k=0; k<this->Nsite_codon; k++)
        {
            for (int l=0; l<this->Nstate_aa; l++)
            {
                is >> tmp; // ssprofiles
            }
        }
        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> tmp; // alloc
        }
        j++;

    }

    if(j == this->startPoint)
    {
        refTree = new Tree(is);

        refTree->RegisterWith(taxonset, 0);


        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> nucp[k];
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> nucrr[k];
        }

        //nucrrnr[0][0]; //AA
        nucrrnr[0][1] =  nucrr[0]; //AC
        nucrrnr[0][2] =  nucrr[1]; //AG
        nucrrnr[0][3] =  nucrr[2]; //AT
        nucrrnr[1][0] =  nucrr[0]; //CA
        //nucrrnr[1][1]; //CG
        nucrrnr[1][2] =  nucrr[3]; //CG
        nucrrnr[1][3] =  nucrr[4]; //CT
        nucrrnr[2][0] =  nucrr[1]; //GA
        nucrrnr[2][1] =  nucrr[3]; //GC
        //nucrrnr[2][2]; //GG
        nucrrnr[2][3] =  nucrr[5]; //GT
        nucrrnr[3][0] =  nucrr[2]; //TA
        nucrrnr[3][1] =  nucrr[4]; //TC
        nucrrnr[3][2] =  nucrr[5]; //TG
        //nucrrnr[3][3]; //TT

        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >>  codonprofile[k];
        }
        is >> omega;
        is >> tmp;
        is >> tmp;
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }

        for (int k=0; k<this->Nsite_codon; k++)
        {
            for (int l=0; l<Nstate_aa; l++)
            {
                is >> ssaaprofiles[k][l];
            }
        }

        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> alloc[k];
        }
    }
    is.close();

    this->getrate = false;

    gtnr[0][0] = 0.0; //AA
    gtnr[0][1] = GetGTR(0,1); //ac
    gtnr[0][2] = GetGTR(0,2); //ag
    gtnr[0][3] = GetGTR(0,3); //at
    gtnr[1][0] = GetGTR(1,0); //ca
    gtnr[1][1] = 0.0; //cc
    gtnr[1][2] = GetGTR(1,2); //cg
    gtnr[1][3] = GetGTR(1,3); //ct
    gtnr[2][0] = GetGTR(2,0); //ga
    gtnr[2][1] = GetGTR(2,1); //gc
    gtnr[2][2] = 0.0; //gg
    gtnr[2][3] = GetGTR(2,3); //gt
    gtnr[3][0] = GetGTR(3,0); //ta
    gtnr[3][1] = GetGTR(3,1); //tc
    gtnr[3][2] = GetGTR(3,2); //tg
    gtnr[3][3] = 0.0; //tt

    SetTreeStuff();

    //cerr << "Nnode : " << refTree->GetNnode() << "\n";
}

void LocalParameters::readChainCodonMutSelFinite(int it)
{
    this->MCMCpointID = it;
//set parameters : posterior specific
    ifstream is((this->chain + ".chain").c_str());
    if (!is)
    {
        cerr << "error: did not find " << this->chain << ".chain\n";
        exit(1);
    }

    //cerr << this->chain;

    int j = 0;
    string tmp = "" ;
    while (j < it)
    {

        is >> tmp; // tree
        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> tmp; //nucp
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> tmp; //nucrr
        }

        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >> tmp; //codon
        }
        is >> tmp; // omega
        is >> tmp; // kappa
        //is >> tmp; // Ncomponents
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }
        //for (int k=0; k<this->Nsite_codon; k++){
        for (int l=0; l<this->Nstate_aa; l++)
        {
            is >> tmp; // ssprofiles
        }
        //}
        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> tmp; // alloc
        }
        j++;

    }

    if(j == it)
    {
        refTree = new Tree(is);

        refTree->RegisterWith(taxonset, 0);


        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> nucp[k];
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> nucrr[k];
        }

        //nucrrnr[0][0]; //AA
        nucrrnr[0][1] =  nucrr[0]; //AC
        nucrrnr[0][2] =  nucrr[1]; //AG
        nucrrnr[0][3] =  nucrr[2]; //AT
        nucrrnr[1][0] =  nucrr[0]; //CA
        //nucrrnr[1][1]; //CG
        nucrrnr[1][2] =  nucrr[3]; //CG
        nucrrnr[1][3] =  nucrr[4]; //CT
        nucrrnr[2][0] =  nucrr[1]; //GA
        nucrrnr[2][1] =  nucrr[3]; //GC
        //nucrrnr[2][2]; //GG
        nucrrnr[2][3] =  nucrr[5]; //GT
        nucrrnr[3][0] =  nucrr[2]; //TA
        nucrrnr[3][1] =  nucrr[4]; //TC
        nucrrnr[3][2] =  nucrr[5]; //TG
        //nucrrnr[3][3]; //TT


        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < Nnucp; nuc1++)
        {
            sum += this->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < Nnucp; nuc1++)
        {
            this->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += this->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                this->nucrrnr[nuc1][nuc2] /= sum;
            }
        }



        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >>  codonprofile[k];
        }
        is >> omega;
        is >> tmp;
        //is >> tmp;
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }

        //for (int k=0; k<this->Nsite_codon; k++){
        for (int l=0; l<Nstate_aa; l++)
        {
            is >> ssaaprofiles[0][l];
        }
        //}

        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> alloc[k];
        }
    }
    is.close();

    this->getrate = false;

    gtnr[0][0] = 0.0; //AA
    gtnr[0][1] = GetGTR(0,1); //ac
    gtnr[0][2] = GetGTR(0,2); //ag
    gtnr[0][3] = GetGTR(0,3); //at
    gtnr[1][0] = GetGTR(1,0); //ca
    gtnr[1][1] = 0.0; //cc
    gtnr[1][2] = GetGTR(1,2); //cg
    gtnr[1][3] = GetGTR(1,3); //ct
    gtnr[2][0] = GetGTR(2,0); //ga
    gtnr[2][1] = GetGTR(2,1); //gc
    gtnr[2][2] = 0.0; //gg
    gtnr[2][3] = GetGTR(2,3); //gt
    gtnr[3][0] = GetGTR(3,0); //ta
    gtnr[3][1] = GetGTR(3,1); //tc
    gtnr[3][2] = GetGTR(3,2); //tg
    gtnr[3][3] = 0.0; //tt

    SetTreeStuff();

    //cerr << "Nnode : " << refTree->GetNnode() << "\n";
}


void LocalParameters::readChainCodonMutSelFinite()
{

//set parameters : posterior specific
    ifstream is((this->chain + ".chain").c_str());
    if (!is)
    {
        cerr << "error: did not find " << this->chain << ".chain\n";
        exit(1);
    }

    //cerr << this->chain;

    int j = 0;
    string tmp = "" ;
    while (j < this->startPoint)
    {

        is >> tmp; // tree
        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> tmp; //nucp
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> tmp; //nucrr
        }

        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >> tmp; //codon
        }
        is >> tmp; // omega
        is >> tmp; // kappa
        //is >> tmp; // Ncomponents
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }
        //for (int k=0; k<this->Nsite_codon; k++){
        for (int l=0; l<this->Nstate_aa; l++)
        {
            is >> tmp; // ssprofiles
        }
        //}
        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> tmp; // alloc
        }
        j++;

    }

    if(j == this->startPoint)
    {
        refTree = new Tree(is);

        refTree->RegisterWith(taxonset, 0);


        is >> tmp; // branchalpha
        is >> tmp; // branchbeta
        for (int k=0; k<this->Nnucp; k++)
        {
            is >> nucp[k];
        }
        for (int k=0; k<this->Nnucrr; k++)
        {
            is >> nucrr[k];
        }

        //nucrrnr[0][0]; //AA
        nucrrnr[0][1] =  nucrr[0]; //AC
        nucrrnr[0][2] =  nucrr[1]; //AG
        nucrrnr[0][3] =  nucrr[2]; //AT
        nucrrnr[1][0] =  nucrr[0]; //CA
        //nucrrnr[1][1]; //CG
        nucrrnr[1][2] =  nucrr[3]; //CG
        nucrrnr[1][3] =  nucrr[4]; //CT
        nucrrnr[2][0] =  nucrr[1]; //GA
        nucrrnr[2][1] =  nucrr[3]; //GC
        //nucrrnr[2][2]; //GG
        nucrrnr[2][3] =  nucrr[5]; //GT
        nucrrnr[3][0] =  nucrr[2]; //TA
        nucrrnr[3][1] =  nucrr[4]; //TC
        nucrrnr[3][2] =  nucrr[5]; //TG
        //nucrrnr[3][3]; //TT

        double sum =  0.0 ;
        for (int nuc1 = 0 ; nuc1 < Nnucp; nuc1++)
        {
            sum += this->nucp[nuc1];
        }
        for (int nuc1 = 0 ; nuc1 < Nnucp; nuc1++)
        {
            this->nucp[nuc1]/=sum;
        }

        sum = 0.0 ;
        for (int nuc1=0; nuc1<4-1; nuc1++)
        {
            for (int nuc2=nuc1+1; nuc2<4; nuc2++)
            {
                sum += this->nucrrnr[nuc1][nuc2];
            }
        }

        for (int nuc1=0; nuc1<4; nuc1++)
        {
            for (int nuc2=0; nuc2<4; nuc2++)
            {
                this->nucrrnr[nuc1][nuc2] /= sum;
            }
        }


        for (int k=0; k<this->Nstate_codon; k++)
        {
            is >>  codonprofile[k];
        }
        is >> omega;
        is >> tmp;
        //is >> tmp;
        for (int k=0; k<this->Nstate_aa; k++)
        {
            is >> tmp;
        }

        //for (int k=0; k<this->Nsite_codon; k++){
        for (int l=0; l<Nstate_aa; l++)
        {
            is >> ssaaprofiles[0][l];
        }
        //}

        for (int k=0; k<this->Nsite_codon; k++)
        {
            is >> alloc[k];
        }
    }
    is.close();

    this->getrate = false;

    gtnr[0][0] = 0.0; //AA
    gtnr[0][1] = GetGTR(0,1); //ac
    gtnr[0][2] = GetGTR(0,2); //ag
    gtnr[0][3] = GetGTR(0,3); //at
    gtnr[1][0] = GetGTR(1,0); //ca
    gtnr[1][1] = 0.0; //cc
    gtnr[1][2] = GetGTR(1,2); //cg
    gtnr[1][3] = GetGTR(1,3); //ct
    gtnr[2][0] = GetGTR(2,0); //ga
    gtnr[2][1] = GetGTR(2,1); //gc
    gtnr[2][2] = 0.0; //gg
    gtnr[2][3] = GetGTR(2,3); //gt
    gtnr[3][0] = GetGTR(3,0); //ta
    gtnr[3][1] = GetGTR(3,1); //tc
    gtnr[3][2] = GetGTR(3,2); //tg
    gtnr[3][3] = 0.0; //tt

    SetTreeStuff();

    //cerr << "Nnode : " << refTree->GetNnode() << "\n";
}

void LocalParameters::SetTreeStuff()
{
    SetTree();
    for (int node_i = 0 ; node_i < this->refTree->GetNnode(); node_i++)
    {
        this->gtrMap[node_i] = 0;

    }

    if (fixgtr1==1)
    {
        GetGTR1();

        Link* a = refTree->GetLCA(this->taxa_gtr1_a,taxa_gtr1_b);
        int notNodeIndex = refTree->GetLCA(this->taxa_gtr1_a,taxa_gtr1_c)->GetNode()->GetIndex();
        cerr << "notNodeIndex\t" << notNodeIndex << " nodeIndex" << a->GetNode()->GetIndex() <<  "\n";
        SetTreeStuffRecursively(a,notNodeIndex,1);

    }

    if (fixgtr2==1)
    {
        GetGTR2();

        Link* a = refTree->GetLCA(this->taxa_gtr2_a,taxa_gtr2_b);
        int notNodeIndex = refTree->GetLCA(this->taxa_gtr2_a,taxa_gtr2_c)->GetNode()->GetIndex();
        cerr << "notNodeIndex\t" << notNodeIndex << " nodeIndex" << a->GetNode()->GetIndex() <<  "\n";
        SetTreeStuffRecursively(a,notNodeIndex,2);

    }

}


void LocalParameters::SetTreeStuffRecursively(Link *from,int notNodeIndex,int gtrIndex)
{

    gtrMap[from->GetNode()->GetIndex()] = gtrIndex;

    auto it = gtrMap.find(from->GetNode()->GetIndex());
    if(it != gtrMap.end())
    {
        cerr << gtrIndex << "\t" <<from->GetNode()->GetIndex() << "\t" << it->second << "\n";
    }
    else
    {
        cerr << "Error while registring gtr along the tree\n";
        exit(0);
    }
    if(from->isRoot())
    {
        for (Link* link = from->Next(); link != from; link = link->Next())
        {
            if (link->Out()->GetNode()->GetIndex() != notNodeIndex)
            {
                SetTreeStuffRecursively(link->Out(),notNodeIndex,gtrIndex);
            }
        }
    }
    else if(!from->isLeaf())
    {
        for (Link* link = from->Next(); link != from; link = link->Next())
        {
            if (link->Out()->GetNode()->GetIndex() != notNodeIndex)
            {
                SetTreeStuffRecursively(link->Out(),notNodeIndex,gtrIndex);
            }
        }
    }
    else if (from->isLeaf())
    {


    }
}

void LocalParameters::toFasta(ofstream &os, int** curent_nodeleaf_sequence_codon)
{
    for(int taxa = 0 ; taxa < Ntaxa ; taxa++)
    {
        os << ">" <<codondata->taxset->GetTaxon(taxa) << "\n";
        for (int site_codon = 0 ; site_codon < Nsite_codon ; site_codon++)
        {
            os << codonstatespace->GetState(curent_nodeleaf_sequence_codon[taxa][site_codon]);
        }
        os << "\n";
    }
}

void LocalParameters::toAli(ofstream &os, int** curent_nodeleaf_sequence_codon)
{
    os << Ntaxa << "\t" << (3* Nsite_codon) << "\n";
    for(int taxa = 0 ; taxa < Ntaxa ; taxa++)
    {
        os <<codondata->taxset->GetTaxon(taxa) << ' ';
        for (int site_codon = 0 ; site_codon < Nsite_codon ; site_codon++)
        {
            os << codonstatespace->GetState(curent_nodeleaf_sequence_codon[taxa][site_codon]);
        }
        os << "\n";
    }
    os << "\n";
}


int LocalParameters::GetPointID()
{
    return MCMCpointID;
}

void LocalParameters::writeRealDataSummaries(ofstream&os, bool headers)
{
    string* arrSummaries = new string[NusedSummaries];
    for(int summary_i = 0 ; summary_i < NSummaries; summary_i++)
    {
        auto it = mapUsedSummaries.find(listSummaries[summary_i]);
        if(it != mapUsedSummaries.end() )
        {
            if(it->second != -1)
            {
                arrSummaries[it->second] = it->first;
            }

        }
    }

    if (headers)
    {
        for(int summary_i = 0 ; summary_i < NusedSummaries; summary_i++)
        {
            if (summary_i < NusedSummaries-1)
            {
                os << arrSummaries[summary_i]  << "\t";
            }
            else
            {
                os << arrSummaries[summary_i] << "\n";
            }
        }
    }

    for(int summary_i = 0 ; summary_i < NusedSummaries; summary_i++)
    {
        if (summary_i < NusedSummaries-1)
        {
            os << summariesRealData[summary_i]  << "\t";
        }
        else
        {
            os << summariesRealData[summary_i] << "\n";
        }
    }
    delete [] arrSummaries;

}

void LocalParameters::writeAncestralDataSummaries(ofstream&os, bool headers)
{
    // should be incorporated to populatio_t
    string* arrSummaries = new string[NusedAncSummaries];
    for(int summary_i = 0 ; summary_i < NSummaries; summary_i++)
    {
        auto it = mapUsedAncSummaries.find(listSummaries[summary_i]);
        if(it != mapUsedAncSummaries.end() )
        {
            if(it->second != -1)
            {
                arrSummaries[it->second] = it->first;
            }

        }
    }

    if (headers)
    {
        for(int summary_i = 0 ; summary_i < NusedAncSummaries; summary_i++)
        {
            if (summary_i < NusedAncSummaries-1)
            {
                os << arrSummaries[summary_i]  << "\t";
            }
            else
            {
                os << arrSummaries[summary_i] << "\n";
            }
        }
    }

    for(int summary_i = 0 ; summary_i < NusedAncSummaries; summary_i++)
    {
        if (summary_i < NusedAncSummaries-1)
        {
            os << summariesAncestralData[summary_i]  << "\t";
        }
        else
        {
            os << summariesAncestralData[summary_i] << "\n";
        }
    }
    delete [] arrSummaries;
}


