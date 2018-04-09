/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LikelihoodFreePhylogenetics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with LikelihoodFreePhylogenetics. If not, see <http://www.gnu.org/licenses/>.
*/

#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "GlobalParameters.h"
#include "LocalParameters.h"
#include "SummaryStatistics.h"
#include "SiteInterSubMatrix.h"
#include "PriorSampler.h"
#include "TreeSimulator.h"
#include "AncestralSequence.h"
#include "Posterior.h"

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



int main(int argc, char* argv[])
{
// Comments

    // program options
    string model = "";
    string controlfile = "";

    try
    {
        if(argc < 2)
        {
            throw(0);
        }
        int i = 1;
        while (i < argc)
        {
            string s = argv[i];
            if (s == "-v" || s == "--version" )
            {
                throw(0);
            }
            else if ( s =="-m")
            {
                i++;
                model = argv[i];
                i++;
                controlfile = argv[i];
            }
            i++;
        }// end while

    }// end try
    catch (...)
    {
        cerr << "\n";
        cerr << "version 1.0\n";
        cerr << "###########################\n";
        cerr << "-m < stats | show | CodonMutSelFiniteABC | CodonMutSelSBDPABC | CodonDegMutSelFiniteABC | CodonDegMutSelSBDPABC | CodonMutSelFinite | CodonMutSelSBDP | CodonMutSelFinitePPred | CodonMutSelSBDPPPred > <controlfile>\n";
        cerr << "###########################\n";
        cerr << "#SUMMARIES\n";
        cerr << "#ANCSUMMARIES\n";
        cerr << "#ACCSUMMARIES\n";
        cerr << "#PARAM\n";
        cerr << "#SSMAP\n";
        cerr << "#MAP\n";
        cerr << "#ANCESTRALMAP\n";
        cerr << "#DIST\n";
        cerr << "#TRANS\n";
        cerr << "#SPEUDODATA\n";
        cerr << "#NRUN\n";
        cerr << "#SAMPLING\n";
        cerr << "#LOCALPARAM\n";
        cerr << "###########################\n";
        exit(1);
    }

    cerr << "models\n";
    if (model == "stats")
    {

        GlobalParameters* gparam = new GlobalParameters(model, controlfile);
        //Posterior* post = new Posterior(gparam);

        LocalData* ldata = new  LocalData(gparam);
        ldata->readLocalData(1);
        SummaryStatistics* ss = new SummaryStatistics(ldata);
        ofstream realDataSummaries_os ((ldata->output+".stats").c_str());
        int k = 1 ;

        while (k < (int) gparam->listGenes.size())
        {
            cerr << "#########\n";
            cerr << k << "\n";
            cerr << "#########\n";

            ldata->readLocalData(k);

            ss->computeSummariesFromData();

            ldata->writeRealDataSummaries(realDataSummaries_os,k==1);
            k++;

        }
        realDataSummaries_os.close();
    }

    else if (model == "show")
    {
        GlobalParameters* gparam = new GlobalParameters(); 
        cerr << "#PARAMETERS\n";
        for (auto i : gparam->listParam)
            cerr << i << "\t"; 
        cerr << "\n";
        cerr << "#SUMMARY\n";
        for (auto i : gparam->listSummaries)
            cerr << i << "\t"; 
        cerr << "\n";
        cerr << "#EVOLSTATS\n";
        for (auto i : gparam->listEvoStats)
            cerr << i << "\t"; 
        cerr << "\n";
        cerr << "#DIST\n";
        for (auto i : gparam->listDistances)
            cerr << i << "\t"; 
        cerr << "\n";
        cerr << "#TRANS\n";
        for (auto i : gparam->listTransformtations)
            cerr << i << "\t"; 
        cerr << "\n";

    }

    else if (model == "CodonMutSelFiniteABC")
    {
        cerr << "CodonMutSelFiniteCpG\n";

        GlobalParameters* gparam = new GlobalParameters(model, controlfile);


        int Npoint = (int) (gparam->chainPointEnd-gparam->chainPointStart)/gparam->chainPointEvery;


        Posterior* post1 = new Posterior(gparam);
        Posterior* post2 = new Posterior(gparam);
        Posterior* post3 = new Posterior(gparam);

        cerr << gparam->chainPointStart << " " << gparam->chainPointEnd << " " << gparam->chainPointEvery << " " << Npoint << "\n";


        LocalParameters** lparam = new LocalParameters*[Npoint];
        SummaryStatistics** ss = new SummaryStatistics*[Npoint];
        PriorSampler** sampler = new PriorSampler*[Npoint];
        SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[Npoint];
        AncestralSequence** ancestraseq = new AncestralSequence*[Npoint];
        TreeSimulator** simulator = new TreeSimulator*[Npoint];


        omp_set_dynamic(0);
        omp_set_num_threads(gparam->Nthread);
        int pt_i;
        cerr << "Npoint\t" << Npoint << "\n";
        #pragma omp parallel for
        for (pt_i = gparam->chainPointStart ; pt_i < gparam->chainPointEnd; pt_i+=gparam->chainPointEvery)
        {
            //int thread_id = omp_get_thread_num();
            //for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
            int l = (int) (pt_i-gparam->chainPointStart)/gparam->chainPointEvery;
            lparam[l] = new LocalParameters(gparam);
            lparam[l]->readChainCodonMutSelFinite(pt_i);
            ss[l] = new SummaryStatistics(lparam[l]);
            ss[l]->computeSummaries();
            sampler[l] = new PriorSampler(lparam[l]);
            submatrix[l] = new SiteInterSubMatrix(lparam[l]);
            ancestraseq[l] = new AncestralSequence(lparam[l]);
            simulator[l] = new TreeSimulator(lparam[l],submatrix[l],ancestraseq[l]);


            if (l == 0)
            {

                post1->SetNsite(lparam[l]->Nsite_codon);
                post2->SetNsite(lparam[l]->Nsite_codon);
                post3->SetNsite(lparam[l]->Nsite_codon);


                ostringstream ost2;
                ost2 <<  gparam->output << ".realdata";
                ofstream realDataSummaries_os (ost2.str());
                lparam[l]->writeRealDataSummaries(realDataSummaries_os);
                realDataSummaries_os.close();


            }

        }

        ifstream monitor_is((gparam->output+".monitor").c_str());
        if(monitor_is)
        {
            monitor_is.close();
            post1->readMonitor(monitor_is);
            monitor_is.close();

            ifstream posterior_is((gparam->output+".post").c_str());
            post1->readPosterior(posterior_is);
            post2->readPosterior(posterior_is);
            post3->readPosterior(posterior_is);
        }

        cerr << "The simulation process started\n";

        while(post3->Niter < post3->Nrun)
        {
            cerr << ".";
            omp_set_dynamic(0);
            omp_set_num_threads(gparam->Nthread);
            #pragma omp parallel
            {

                #pragma omp for
                for (int l = 0 ; l < Npoint; l++)
                {

                    sampler[l]->sample();

                    simulator[l]->GetNewSimulatedCodonAlignment();

                    ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);

                    #pragma omp critical
                    {


                        if(post3->Niter < 100000 && post3->Nrun == 100000)
                        {

                            post1->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );

                            post2->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );

                            post3->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );


                            if (post3->Niter % post3->threshold == 0)
                            {


                                ofstream dist_os1((gparam->output+"-100K.post").c_str(),OUT);
                                post1->writeHeader(dist_os1);
                                post1->writePosterior(dist_os1);
                                dist_os1.close();

                                ofstream monitor_os1((gparam->output+"-100K.monitor").c_str(),OUT);
                                post1->writeMonitorPosterior(monitor_os1);
                                monitor_os1.close();


                            }



                        }

                        if(post3->Niter < 1000000 && post3->Nrun == 1000000)
                        {

                            post2->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );

                            post3->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );


                            if (post3->Niter % post3->threshold == 0)
                            {


                                ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
                                post2->writeHeader(dist_os2);
                                post2->writePosterior(dist_os2);
                                dist_os2.close();

                                ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
                                post2->writeMonitorPosterior(monitor_os2);
                                monitor_os2.close();

                            }

                        }


                        if(post3->Niter < 10000000  && post3->Nrun == 10000000 )
                        {

                            post3->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );


                            if (post3->Niter % post3->threshold == 0)
                            {

                                ofstream dist_os3((gparam->output+"-10M.post").c_str(),OUT);
                                post3->writeHeader(dist_os3);
                                post3->writePosterior(dist_os3);
                                dist_os3.close();

                                ofstream monitor_os3((gparam->output+"-10M.monitor").c_str(),OUT);
                                post3->writeMonitorPosterior(monitor_os3);
                                monitor_os3.close();



                            }

                        }



                    }
                }

            }

        }

        cerr << "End of the simulation process\n";
        exit(0);

    }
    else if (model == "CodonMutSelSBDPABC") 
    {
        cerr << "CodonMutSelSBDPCpG\n";

        GlobalParameters* gparam = new GlobalParameters(model, controlfile);
        cerr << gparam->chainPointStart << " " << gparam->chainPointEnd << " " << gparam->chainPointEvery << "\n";
        int Npoint = (int) (gparam->chainPointEnd-gparam->chainPointStart)/gparam->chainPointEvery;
//        Posterior* post1 = new Posterior(gparam);
//        Posterior* post2 = new Posterior(gparam);
        Posterior* post3 = new Posterior(gparam);

        LocalParameters** lparam = new LocalParameters*[Npoint];

        if(gparam->verbose)
        {
            cerr << "debug1\n";
        }

        if(gparam->verbose)
        {
            cerr << "debug2\n";
        }
        cerr << "Npoint\t" << Npoint << "\n";
        if(gparam->verbose)
        {
            cerr << "debug3\n";
        }
        SummaryStatistics** ss = new SummaryStatistics*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug4\n";
        }
        PriorSampler** sampler = new PriorSampler*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug5\n";
        }
        SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug6\n";
        }
        AncestralSequence** ancestraseq = new AncestralSequence*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug7\n";
        }
        TreeSimulator** simulator = new TreeSimulator*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug8\n";
        }




        omp_set_dynamic(0);
        omp_set_num_threads(gparam->Nthread);
        int pt_i;
        #pragma omp parallel for
        for (pt_i = gparam->chainPointStart ; pt_i < gparam->chainPointEnd; pt_i+=gparam->chainPointEvery)
        {
            //int thread_id = omp_get_thread_num();
            //for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
            int l = (int) (pt_i-gparam->chainPointStart)/gparam->chainPointEvery;
            if(gparam->verbose)
            {
                cerr << "debug9\n";
                cerr << (l) << "\n";
            }
            lparam[l] = new LocalParameters(gparam);
            if(gparam->verbose)
            {
                cerr << "debug10\n";
            }
            lparam[l]->readChainCodonMutSelSBDP(pt_i);
            if(gparam->verbose)
            {
                cerr << "debug11\n";
            }
            ss[l] = new SummaryStatistics(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug12\n";
            }
            ss[l]->computeSummaries();
            if(gparam->verbose)
            {
                cerr << "debug13\n";
            }
            sampler[l] = new PriorSampler(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug14\n";
            }
            submatrix[l] = new SiteInterSubMatrix(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug15\n";
            }
            ancestraseq[l] = new AncestralSequence(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug16\n";
            }
            simulator[l] = new TreeSimulator(lparam[l],submatrix[l],ancestraseq[l]);
            if(gparam->verbose)
            {
                cerr << "debug17\n";
            }


            if (l == 0)
            {

//                post1->SetNsite(lparam[l]->Nsite_codon);
//                post2->SetNsite(lparam[l]->Nsite_codon);
                post3->SetNsite(lparam[l]->Nsite_codon);


                ostringstream ost1;
                ost1 <<  gparam->output << ".inputparam";
                ofstream lparam_os (ost1.str());
                lparam[l]->writeParam(lparam_os);
                lparam_os.close();


                ostringstream ost2;
                ost2 <<  gparam->output << ".realdata";
                ofstream realDataSummaries_os (ost2.str());
                lparam[l]->writeRealDataSummaries(realDataSummaries_os);
                realDataSummaries_os.close();


            }

        }



        ifstream monitor_is((gparam->output+"-1M.monitor").c_str());
        ifstream monitor_is_100K((gparam->output+"-100K.monitor").c_str());
        if(monitor_is)
        {
            cerr << "monitor 1M\n";
            post3->readMonitor(monitor_is);
            monitor_is.close();

            ifstream posterior_is((gparam->output+"-1M.post").c_str());
            if (!posterior_is)
            {
                cerr << "error: did not find posteriorfile"<< "\n";
                exit(1);
            }

            post3->readPosterior(posterior_is);
            posterior_is.close();

            cerr << post3->Niter << " on " << post3->Nrun << "\n";

        }
        else if (monitor_is_100K)
        {
            cerr << "monitor 100K\n";
            post3->readMonitor(monitor_is_100K);
            monitor_is_100K.close();

            ifstream posterior_is((gparam->output+"-100K.post").c_str());
            if (!posterior_is)
            {
                cerr << "error: did not find posteriorfile"<< "\n";
                exit(1);
            }

            post3->readPosterior(posterior_is);
            posterior_is.close();

            cerr << post3->Niter << " on " << post3->Nrun << "\n";

        }

        cerr << "The simulation process started\n";


        while(post3->Niter < post3->Nrun)
        {
            cerr << ".";
            omp_set_dynamic(0);
            omp_set_num_threads(gparam->Nthread);
            #pragma omp parallel
            {


                #pragma omp for
                for (int l = 0 ; l < Npoint; l++)
                {


                    if(gparam->verbose)
                    {
                        cerr << "debug18\n";
                    }
                    sampler[l]->sample();

                    if(gparam->verbose)
                    {
                        cerr << "debug19\n";
                    }
                    simulator[l]->GetNewSimulatedCodonAlignment();
                    if(gparam->verbose)
                    {
                        cerr << "debug20\n";
                    }

                    /* for (int interval_i = 0 ; interval_i < lparam[l]->Ninterval; interval_i++)
                    {
                        ss[l]->computeSummariesAncestralSequence(simulator[l]->CurrentAncestralCodonSequence[interval_i]);
                    } */

                    if(gparam->verbose)
                    {
                        cerr << "debug20.1\n";
                    }

                    ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);
                    if(gparam->verbose)
                    {
                        cerr << "debug21\n";
                    }

                    #pragma omp critical
                    {
                        
                        if(post3->Niter < 100000 && post3->Nrun < 100000)
                        {

                            post3->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );

                            if (post3->Niter % post3->threshold == 0)
                            {


                                ofstream dist_os1((gparam->output+".post").c_str(),OUT);
                                post3->writeHeader(dist_os1);
                                post3->writePosterior(dist_os1);
                                dist_os1.close();

                                ofstream monitor_os1((gparam->output+".monitor").c_str(),OUT);
                                post3->writeMonitorPosterior(monitor_os1);
                                monitor_os1.close();


                            }



                        }


                        if(post3->Niter < 100000 && post3->Nrun == 100000)
                        {
                            post3->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );

                            if (post3->Niter % post3->threshold == 0)
                            {

                                ofstream dist_os1((gparam->output+"-100K.post").c_str(),OUT);
                                post3->writeHeader(dist_os1);
                                post3->writePosterior(dist_os1);
                                dist_os1.close();

                                ofstream monitor_os1((gparam->output+"-100K.monitor").c_str(),OUT);
                                post3->writeMonitorPosterior(monitor_os1);
                                monitor_os1.close();

                            }
                        }

                        if(post3->Niter < 1000000 && post3->Nrun == 1000000)
                        {
                            post3->registerNewSimulation(
                                lparam[l]->MCMCpointID,
                                lparam[l]->GetCurrentParameters(),
                                lparam[l]->GetCurrentSummaries(),
                                lparam[l]->GetCurrentAccessorySummaries(),
                                lparam[l]->GetCurrentAncEvoStats(),
                                lparam[l]->GetCurrentEvoStats(),
                                lparam[l]->GetCurrentSiteSpecificEvoStats(),
                                lparam[l]->GetCurrentDistances(),
                                lparam[l]->GetCurrentWeights()
                            );

                            if (post3->Niter % post3->threshold == 0)
                            {
                                ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
                                post3->writeHeader(dist_os2);
                                post3->writePosterior(dist_os2);
                                dist_os2.close();

                                ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
                                post3->writeMonitorPosterior(monitor_os2);
                                monitor_os2.close();

                            }

                            if (post3->Niter == 100000)
                            {
                                ofstream dist_os1((gparam->output+"-100K.post").c_str(),OUT);
                                post3->writeHeader(dist_os1);
                                post3->writePosterior(dist_os1);
                                dist_os1.close();

                                ofstream monitor_os1((gparam->output+"-100K.monitor").c_str(),OUT);
                                post3->writeMonitorPosterior(monitor_os1);
                                monitor_os1.close();
                            }
                        }
                    }
                }
            }
        }
        cerr << "End of the simulation process\n";
        exit(0);
    }
    else if (model == "CodonMutSelSBDPABC-v2")
    {
        cerr << "CodonMutSelSBDPCpG-v2\n";

        GlobalParameters* gparam = new GlobalParameters(model, controlfile);
        cerr << gparam->chainPointStart << " " << gparam->chainPointEnd << " " << gparam->chainPointEvery << "\n";
        int Npoint = (int) (gparam->chainPointEnd-gparam->chainPointStart)/gparam->chainPointEvery;
//        Posterior* post1 = new Posterior(gparam);
//        Posterior* post2 = new Posterior(gparam);
        Posterior* postMaster = new Posterior(gparam);

        Posterior** postSlave = new Posterior*[Npoint];

        LocalParameters** lparam = new LocalParameters*[Npoint];

        if(gparam->verbose)
        {
            cerr << "debug1\n";
        }

        if(gparam->verbose)
        {
            cerr << "debug2\n";
        }
        cerr << "Npoint\t" << Npoint << "\n";
        if(gparam->verbose)
        {
            cerr << "debug3\n";
        }
        SummaryStatistics** ss = new SummaryStatistics*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug4\n";
        }
        PriorSampler** sampler = new PriorSampler*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug5\n";
        }
        SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug6\n";
        }
        AncestralSequence** ancestraseq = new AncestralSequence*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug7\n";
        }
        TreeSimulator** simulator = new TreeSimulator*[Npoint];
        if(gparam->verbose)
        {
            cerr << "debug8\n";
        }




        omp_set_dynamic(0);
        omp_set_num_threads(gparam->Nthread);
        int pt_i;
        #pragma omp parallel for
        for (pt_i = gparam->chainPointStart ; pt_i < gparam->chainPointEnd; pt_i+=gparam->chainPointEvery)
        {


            //int thread_id = omp_get_thread_num();
            //for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
            int l = (int) (pt_i-gparam->chainPointStart)/gparam->chainPointEvery;
            if(gparam->verbose)
            {
                cerr << "debug9\n";
                cerr << (l) << "\n";
            }


            postSlave[l] = new Posterior(gparam);


            lparam[l] = new LocalParameters(gparam);
            if(gparam->verbose)
            {
                cerr << "debug10\n";
            }
            lparam[l]->readChainCodonMutSelSBDP(pt_i);
            if(gparam->verbose)
            {
                cerr << "debug11\n";
            }
            ss[l] = new SummaryStatistics(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug12\n";
            }
            ss[l]->computeSummaries();
            if(gparam->verbose)
            {
                cerr << "debug13\n";
            }
            sampler[l] = new PriorSampler(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug14\n";
            }
            submatrix[l] = new SiteInterSubMatrix(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug15\n";
            }
            ancestraseq[l] = new AncestralSequence(lparam[l]);
            if(gparam->verbose)
            {
                cerr << "debug16\n";
            }
            simulator[l] = new TreeSimulator(lparam[l],submatrix[l],ancestraseq[l]);
            if(gparam->verbose)
            {
                cerr << "debug17\n";
            }

            postSlave[l]->SetNsite(lparam[l]->Nsite_codon);

            if (l == 0)
            {


                ostringstream ost1;
                ost1 <<  gparam->output << ".inputparam";
                ofstream lparam_os (ost1.str());
                lparam[l]->writeParam(lparam_os);
                lparam_os.close();


                ostringstream ost2;
                ost2 <<  gparam->output << ".realdata";
                ofstream realDataSummaries_os (ost2.str());
                lparam[l]->writeRealDataSummaries(realDataSummaries_os);
                realDataSummaries_os.close();


            }

        }

        cerr << "The simulation process started\n";
        cerr << postMaster->Niter << " on " << postMaster->Nrun << "\n";

        ifstream monitor_is((gparam->output+"-1M.monitor").c_str());
        ifstream monitor_is_100K((gparam->output+"-100K.monitor").c_str());
        if(monitor_is)
        {

            postMaster->readMonitor(monitor_is);
            monitor_is.close();

            ifstream posterior_is((gparam->output+"-1M.post").c_str());
            if (!posterior_is)
            {
                cerr << "error: did not find posteriorfile"<< "\n";
                exit(1);
            }

            postMaster->readPosterior(posterior_is);
            posterior_is.close();

        }
        else if (monitor_is_100K)
        {

            postMaster->readMonitor(monitor_is_100K);
            monitor_is_100K.close();

            ifstream posterior_is((gparam->output+"-100K.post").c_str());
            if (!posterior_is)
            {
                cerr << "error: did not find posteriorfile"<< "\n";
                exit(1);
            }

            postMaster->readPosterior(posterior_is);
            posterior_is.close();

        }

        int runTodo = int((postMaster->Nrun-postMaster->Niter)/Npoint);
        int run = 0;
//        while(run < runTodo)
//        {

        omp_set_dynamic(0);
        omp_set_num_threads(gparam->Nthread);

        #pragma omp parallel for private(run) collapse(2)
        for (int l = 0 ; l < Npoint; l++)
        {
            for (run  = 0; run < runTodo; run++)

            {
                if(gparam->verbose)
                {
                    cerr << "debug18\n";
                }
                sampler[l]->sample();

                if(gparam->verbose)
                {
                    cerr << "debug19\n";
                }
                simulator[l]->GetNewSimulatedCodonAlignment();
                if(gparam->verbose)
                {
                    cerr << "debug20\n";
                }

                for (int interval_i = 0 ; interval_i < 11; interval_i++)
                {
                    ss[l]->computeSummariesAncestralSequence(simulator[l]->CurrentAncestralCodonSequence[interval_i]);
                }

                if(gparam->verbose)
                {
                    cerr << "debug20.1\n";
                }

                ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);
                if(gparam->verbose)
                {
                    cerr << "debug21\n";
                }

                postSlave[l]->slaveRegisterNewSimulation(
                    lparam[l]->MCMCpointID,
                    lparam[l]->GetCurrentParameters(),
                    lparam[l]->GetCurrentSummaries(),
                    lparam[l]->GetCurrentAccessorySummaries(),
                    lparam[l]->GetCurrentAncEvoStats(),
                    lparam[l]->GetCurrentEvoStats(),
                    lparam[l]->GetCurrentSiteSpecificEvoStats(),
                    lparam[l]->GetCurrentDistances(),
                    lparam[l]->GetCurrentWeights()
                );

                cerr << ".";
            }
//          }

//            run++;

        }

        for (int l = 0 ; l < Npoint; l++)
        {
            postMaster->slaveToMaster(postSlave[l]->population_t);

        }

        postMaster->sortPopulation();

        ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
        postMaster->writeHeader(dist_os2);
        postMaster->writePosterior(dist_os2);
        dist_os2.close();

//        ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
//        postMaster->writeMonitorPosterior(monitor_os2);
//        monitor_os2.close();

        cerr << "End of the simulation process\n";
        exit(0);

    }
    else if (model == "FMutSelSimu")
    {

        GlobalParameters* gparam = new GlobalParameters(model, controlfile);

        if(gparam->verbose)
        {
            cerr << "debug2\n";
        }
        LocalParameters* lparam =  new LocalParameters(gparam);
        lparam->readFMutSelCodeML();

        if(gparam->verbose)
        {
            cerr << "debug3\n";
        }
        Posterior* post = new Posterior(gparam);
        post->SetNsite(lparam->Nsite_codon);


        if(gparam->verbose)
        {
            cerr << "debug4\n";
        }
        SummaryStatistics* ss = new SummaryStatistics(lparam);

        if(gparam->verbose)
        {
            cerr << "debug5\n";
        }
        ss->computeSummaries();

        if(gparam->verbose)
        {
            cerr << "debug6\n";
        }
        PriorSampler* sampler = new PriorSampler(lparam);

        if(gparam->verbose)
        {
            cerr << "debug7\n";
        }
        SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);

        if(gparam->verbose)
        {
            cerr << "debug8\n";
        }
        AncestralSequence* ancestraseq = new AncestralSequence(lparam);


        TreeSimulator* simulator = new TreeSimulator(lparam,submatrix,ancestraseq);


        ofstream lparam_os ((gparam->output+".inputparam").c_str());
        lparam->writeParam(lparam_os);
        lparam_os.close();

        ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
        lparam->writeRealDataSummaries(realDataSummaries_os);
        realDataSummaries_os.close();

        while(post->Niter < post->Nrun)
        {
            if(gparam->verbose)
            {
                cerr << "debug11\n";
            }
            sampler->sample();

            if(gparam->verbose)
            {
                cerr << "debug12\n";
            }
            simulator->GetNewSimulatedCodonAlignment();

            if(gparam->verbose)
            {
                cerr << "debug13\n";
            }
            ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);


            ss->computeSummariesAncestralSequence(simulator->CurrentAncestralCodonSequence[10]);

            if (post->Niter == 0)
            {
                ofstream AncestralDataSummaries_os ((gparam->output+".ancestral").c_str(), OUT);
                bool headers = true;
                lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,headers);
                AncestralDataSummaries_os.close();
            }
            else

            {
                ofstream AncestralDataSummaries_os ((gparam->output+".ancestral").c_str(), APPEND);
                bool headers = false;
                lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,headers);
                AncestralDataSummaries_os.close();
            }

            if(gparam->verbose)
            {
                cerr << "debug14\n";
            }
            post->registerNewSimulation(
                1,
                lparam->GetCurrentParameters(),
                lparam->GetCurrentSummaries(),
                lparam->GetCurrentAccessorySummaries(),
                lparam->GetCurrentAncEvoStats(),
                lparam->GetCurrentEvoStats(),
                lparam->GetCurrentSiteSpecificEvoStats(),
                lparam->GetCurrentDistances(),
                lparam->GetCurrentWeights()
            );

            if (lparam->tofasta)
            {
                ostringstream oss;
                oss << gparam->output << "-" << post->Niter << ".phylip";
                string output = oss.str();
                ofstream ali_os((output).c_str(),OUT);
                lparam->toAli(ali_os,simulator->CurrentLeafNodeCodonSequences);
                ali_os.close();
            }
            cerr << ".";

        }
        cerr << "End of the simulation process\n";
        if(gparam->verbose)
        {
            cerr << "debug15\n";
        }
        ofstream dist_os((gparam->output+".post").c_str(),OUT);
        post->writeHeader(dist_os);
        post->writePosterior(dist_os);
        dist_os.close();


        ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
        post->writePosteriorPredictiveStatistics(ppp_os,lparam->summariesRealData);
        ppp_os.close();

        exit(0);

    }
    else if (model == "CodonMutSelSBDP")
    {
        cerr << "CodonMutSelSBDP" << "\n";
        GlobalParameters* gparam = new GlobalParameters(model, controlfile);
        Posterior* post = new Posterior(gparam);

        LocalParameters* lparam =  new LocalParameters(gparam);
        lparam->readChainCodonMutSelSBDP();


        ofstream lparam_os ((gparam->output+".inputparam").c_str());
        lparam->writeParam(lparam_os);
        lparam_os.close();

        SummaryStatistics* ss =new SummaryStatistics(lparam);
        ss->computeSummaries();

        ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
        lparam->writeRealDataSummaries(realDataSummaries_os);
        realDataSummaries_os.close();

        PriorSampler* sampler = new PriorSampler(lparam);
        SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
        AncestralSequence* ancestraseq = new AncestralSequence(lparam);
        TreeSimulator* simulator = new TreeSimulator(lparam,submatrix,ancestraseq);

        cerr << "Start simulating \n" ;
        while(lparam->startPoint < lparam->endPoint)
        {
            lparam->readChainCodonMutSelSBDP();

            int iter = 0 ;
            while(iter < post->Nrun)
            {



                sampler->sample();

                simulator->GetNewSimulatedCodonAlignment();

                ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

                if (lparam->sampleAncSeq)
                {
                    for (int sample_i = 0 ; sample_i < lparam->Ninterval; sample_i++){
                        ss->computeSummariesAncestralSequence(simulator->CurrentAncestralCodonSequence[sample_i]);

                        if (post->Niter == 0)
                        {
                            ofstream AncestralDataSummaries_os ((gparam->output+".ancestral").c_str(), OUT);
                            bool headers = true;
                            lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,headers);
                            AncestralDataSummaries_os.close();
                        }
                        else
                        {
                            ofstream AncestralDataSummaries_os ((gparam->output+".ancestral").c_str(), APPEND);
                            bool headers = false;
                            lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,headers);
                            AncestralDataSummaries_os.close();
                        }
                    }
                }

                post->registerNewSimulation(
                    lparam->startPoint,
                    lparam->GetCurrentParameters(),
                    lparam->GetCurrentSummaries(),
                    lparam->GetCurrentAccessorySummaries(),
                    lparam->GetCurrentAncEvoStats(),
                    lparam->GetCurrentEvoStats(),
                    lparam->GetCurrentSiteSpecificEvoStats(),
                    lparam->GetCurrentDistances(),
                    lparam->GetCurrentWeights()
                );

                if (lparam->tofasta)
                {
                    ostringstream oss;
                    oss << gparam->output << "-pt" <<lparam->startPoint <<"-" << iter << ".phylip";
                    string output = oss.str();
                    ofstream ali_os((output).c_str(),OUT);
                    lparam->toAli(ali_os,simulator->CurrentLeafNodeCodonSequences);
                    ali_os.close();
                }
                iter++;

            }

            int rep = 1;
            while (rep < lparam->everyPoint)
            {
                rep++;
                lparam->incrementStartPoint();

            }
            lparam->incrementStartPoint();
            cerr << ".";

        }
        cerr << "End of the simulation process\n";

        ofstream dist_os((gparam->output+".post").c_str(),OUT);
        post->writeHeader(dist_os);
        post->writePosterior(dist_os);
        dist_os.close();

        ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
        post->writePosteriorPredictiveStatistics(ppp_os,lparam->summariesRealData);
        ppp_os.close();

        exit(0);


    }
    else if (model == "CodonMutSelFinite")
    {
        cerr << "CodonMutSelFinite" << "\n";
        GlobalParameters* gparam = new GlobalParameters(model, controlfile);
        Posterior* post = new Posterior(gparam);
        if(gparam->verbose)
        {
            cerr << "debug1\n";
        }
        LocalParameters* lparam =  new LocalParameters(gparam);
        lparam->readChainCodonMutSelFinite();
        if(gparam->verbose)
        {
            cerr << "debug2\n";
        }
        ofstream lparam_os ((gparam->output+".inputparam").c_str());
        lparam->writeParam(lparam_os);
        lparam_os.close();
        if(gparam->verbose)
        {
            cerr << "debug3\n";
        }
        SummaryStatistics* ss =new SummaryStatistics(lparam);
        ss->computeSummaries();
        if(gparam->verbose)
        {
            cerr << "debug4\n";
        }
        ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
        lparam->writeRealDataSummaries(realDataSummaries_os);
        realDataSummaries_os.close();
        if(gparam->verbose)
        {
            cerr << "debug5\n";
        }
        PriorSampler* sampler = new PriorSampler(lparam);
        SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
        AncestralSequence* ancestraseq = new AncestralSequence(lparam);
        TreeSimulator* simulator = new TreeSimulator(lparam,submatrix,ancestraseq);
        if(gparam->verbose)
        {
            cerr << "debug6\n";
        }
        cerr << "Start simulating \n" ;
        while(lparam->startPoint < lparam->endPoint)
        {
            lparam->readChainCodonMutSelFinite();
            int iter = 0 ;
            while(iter < post->Nrun)
            {



                sampler->sample();
                if(gparam->verbose)
                {
                    cerr << "debug7\n";
                }
                simulator->GetNewSimulatedCodonAlignment();
                if(gparam->verbose)
                {
                    cerr << "debug8\n";
                }
                ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);
                if(gparam->verbose)
                {
                    cerr << "debug9\n";
                }
                post->registerNewSimulation(
                    lparam->startPoint,
                    lparam->GetCurrentParameters(),
                    lparam->GetCurrentSummaries(),
                    lparam->GetCurrentAccessorySummaries(),
                    lparam->GetCurrentAncEvoStats(),
                    lparam->GetCurrentEvoStats(),
                    lparam->GetCurrentSiteSpecificEvoStats(),
                    lparam->GetCurrentDistances(),
                    lparam->GetCurrentWeights()
                );
                if(gparam->verbose)
                {
                    cerr << "debug10\n";
                }
                if (lparam->tofasta)
                {
                    ostringstream oss;
                    oss << gparam->output << "-pt" <<lparam->startPoint <<"-" << iter << ".phylip";
                    string output = oss.str();
                    ofstream ali_os((output).c_str(),OUT);
                    lparam->toAli(ali_os,simulator->CurrentLeafNodeCodonSequences);
                    ali_os.close();
                }
                iter++;

            }

            int rep = 1;
            while (rep < lparam->everyPoint)
            {
                rep++;
                lparam->incrementStartPoint();

            }
            lparam->incrementStartPoint();
            cerr << ".";

        }
        cerr << "End of the simulation process\n";
        if(gparam->verbose)
        {
            cerr << "debug11\n";
        }
        ofstream dist_os((gparam->output+".post").c_str(),OUT);
        post->writeHeader(dist_os);
        post->writePosterior(dist_os);
        dist_os.close();
        if(gparam->verbose)
        {
            cerr << "debug12\n";
        }

        ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
        post->writePosteriorPredictiveStatistics(ppp_os,lparam->summariesRealData);
        ppp_os.close();

        exit(0);



    }
    else if (model == "CodonMutSelFinitePPPred" || model == "CodonMutSelSBDPPPred")
    {
        // the chain pointS are extract from the posterior file according to chainID
        cerr << model  << "\n";

        GlobalParameters* gparam = new GlobalParameters(model, controlfile);
        Posterior* post = new Posterior(gparam);
        LocalParameters* lparam =  new LocalParameters(gparam);

        if (model == "CodonMutSelSBDPPPred")
        {

            lparam->readChainCodonMutSelSBDP();

        }
        else if (model == "CodonMutSelFinitePPred")
        {

            lparam->readChainCodonMutSelFinite();

        }
        /* 
        ofstream lparam_os ((gparam->output+".inputparam").c_str());
        lparam->writeParam(lparam_os);
        lparam_os.close();
        */

        SummaryStatistics* ss = new SummaryStatistics(lparam);
        ss->computeSummaries();

        ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
        lparam->writeRealDataSummaries(realDataSummaries_os);
        realDataSummaries_os.close();

        if(gparam->verbose)
        {
            cerr << "debug1\n";
        }

        //PriorSampler* sampler = new PriorSampler(lparam);

        if(gparam->verbose)
        {
            cerr << "debug2\n";
        }

        SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);

        if(gparam->verbose)
        {
            cerr << "debug3\n";
        }

        AncestralSequence* ancestraseq = new AncestralSequence(lparam);

        if(gparam->verbose)
        {
            cerr << "debug4\n";
        }

        cerr << lparam->Nsite_codon << "\n";
        TreeSimulator* simulator = new TreeSimulator(lparam,submatrix,ancestraseq);

        if(gparam->verbose)
        {
            cerr << "debug5\n";
        }
        post->readPosterior(lparam->posteriorfile);

        cerr << "The simulation process started\n";
        if (!post->posterior.empty())
        {

            int it = 0 ;
            while(it < post->threshold)
            {
                int point = static_cast<int> (lparam->rnd->Uniform() * post->posterior.size());
                lparam->SetCurrentParametersFromPosterior(post->posterior,point);
                
                if (model == "CodonMutSelSBDPPPred")
                {
                    if(gparam->verbose)
                    {
                        cerr << "CodonMutSelSBDPPPred " <<"debug6\n";
                    }
                    lparam->readChainCodonMutSelSBDP(lparam->GetPointID());

                }
                else if (model == "CodonMutSelFinitePPred")
                {
                    if(gparam->verbose)
                    {
                        cerr << "CodonMutSelFinitePPred " <<"debug6\n";
                    }
                    lparam->readChainCodonMutSelFinite(lparam->GetPointID());

                }
                int rep = 0;
                while (rep < post->Nrun)
                {
                    rep++;
                    if(gparam->verbose)
                    {
                        cerr << rep << " " << point <<" debug7\n";
                    }
                    simulator->GetNewSimulatedCodonAlignment();

                    if(gparam->verbose)
                    {
                        cerr <<" debug7.1\n";
                    }
                    ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);
                    if(gparam->verbose)
                    {
                        cerr << " debug7.2\n";
                    }
/* 
                    for (int interval_i = 0 ; interval_i < lparam[l]->Ninterval; interval_i++)
                    {
                        ss[l]->computeSummariesAncestralSequence(simulator[l]->CurrentAncestralCodonSequence[interval_i]);
                        if (post->Niter == 0)
                        {
                            ofstream AncestralDataSummaries_os ((gparam->output+".ancestral").c_str(), OUT);
                            bool headers = true;
                            lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,headers);
                            AncestralDataSummaries_os.close();
                        }
                        else
                        {
                            ofstream AncestralDataSummaries_os ((gparam->output+".ancestral").c_str(), APPEND);
                            bool headers = false;
                            lparam->writeAncestralDataSummaries(AncestralDataSummaries_os,headers);
                            AncestralDataSummaries_os.close();
                        }                     
                    }                        
 */
                    post->registerNewSimulation(
                        lparam->GetPointID(),
                        lparam->GetCurrentParameters(),
                        lparam->GetCurrentSummaries(),
                        lparam->GetCurrentAccessorySummaries(),
                        lparam->GetCurrentAncEvoStats(),
                        lparam->GetCurrentEvoStats(),
                        lparam->GetCurrentSiteSpecificEvoStats(),
                        lparam->GetCurrentDistances(),
                        lparam->GetCurrentWeights()
                    );
                    it++;
                }

                cerr << ".";
            }

            cerr << "End of the simulation process\n";

            ofstream dist_os((gparam->output+".post").c_str(),OUT);
            post->writeHeader(dist_os);
            post->writePosterior(dist_os);
            dist_os.close();

            ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
            post->writePosteriorPredictiveStatistics(ppp_os,lparam->summariesRealData);
            ppp_os.close();
        }

    }
    // end main
}








