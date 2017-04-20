
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



int main(int argc, char* argv[]){
// Comments

    // program options
    string model = "";
    string controlfile = "";




        try {
            if(argc < 2) {
                throw(0);
            }
            int i = 1;
            while (i < argc) {
                string s = argv[i];
                if (s == "-v" || s == "--version" ) {
                   throw(0);
                } else if ( s =="-m") {
                    i++;
                    model = argv[i];
                    i++;
                    controlfile = argv[i];
                }
                i++;
            }// end while

        }// end try
        catch (...) {
                    cerr << "\n";
                    cerr << "version 1.0\n";
                    cerr << "###########################\n";
                    cerr << "-m < stats | CodonMutSelFiniteABC | CodonMutSelSBDPABC | CodonMutSelFinite | CodonMutSelSBDP | MAP100CodonMutSelFinitePPChecks | MAP100CodonMutSelSBDPPPChecks > <controlfile>\n";
                    cerr << "###########################\n";
                    exit(1);
        }

        cerr << "models\n";
        if (model == "stats") {

            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
            Posterior* post = new Posterior(gparam);

            LocalData* ldata = new  LocalData(gparam);
            ldata->readLocalData(1);
            SummaryStatistics* ss = new SummaryStatistics(ldata);
            ofstream realDataSummaries_os ((ldata->output+".stats").c_str());
            int k = 1 ;

            while (k < gparam->listGenes.size()){
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
//        else if (model == "FMSCpG" ) {
//            cerr << model <<"\n";
//
//            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
//            Posterior* post = new Posterior(gparam);
//
//            LocalParameters** lparam = new LocalParameters*[gparam->Nthread];
//            SummaryStatistics** ss = new SummaryStatistics*[gparam->Nthread];
//            PriorSampler** sampler = new PriorSampler*[gparam->Nthread];
//            SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[gparam->Nthread];
//            AncestralSequence** ancestraseq = new AncestralSequence*[gparam->Nthread];
//            TreeSimulator** simulator = new TreeSimulator*[gparam->Nthread];
//
//            omp_set_dynamic(0);
//            omp_set_num_threads(gparam->Nthread);
//            #pragma omp parallel for
//            for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++) {
//                int thread_id = omp_get_thread_num();
//            //for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
//
//                lparam[thread_id] = new LocalParameters(gparam);
//                lparam[thread_id]->readChainCodonMutSelFinite();
//                ss[thread_id] = new SummaryStatistics(lparam[thread_id]);
//                ss[thread_id]->computeSummaries();
//                sampler[thread_id] = new PriorSampler(lparam[thread_id]);
//                submatrix[thread_id] = new SiteInterSubMatrix(lparam[thread_id]);
//                ancestraseq[thread_id] = new AncestralSequence(lparam[thread_id]);
//                simulator[thread_id] = new TreeSimulator(lparam[thread_id],submatrix[thread_id],ancestraseq[thread_id]);
//
//
//                if (thread_id == 0) {
//
//                    ofstream lparam_os ((gparam->output+".inputparam").c_str());
//                    lparam[0]->writeParam(lparam_os);
//                    lparam_os.close();
//
//
//                    ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
//                    lparam[0]->writeRealDataSummaries(realDataSummaries_os);
//                    realDataSummaries_os.close();
//
//                    ifstream monitor_is((gparam->output+".monitor").c_str());
//                    if (!monitor_is){
//                        monitor_is.close();
//                        ofstream monitor_os((gparam->output+".monitor").c_str(),OUT);
//                        monitor_os.close();
//                    } else {
//                        post->readMonitorPosterior(monitor_is);
//                        monitor_is.close();
//                    }
//
//                }
//
//
//            }
//
//
//
//
//
//            cerr << "The simulation process started\n";
//            cerr << post->Niter << " " << post->Nrun << "\n";
//
//
//
//            while(post->Niter < post->Nrun) {
//                     cerr << ".";
//
//                    #pragma omp parallel
//                    {
//
////                        std::vector<double> param;
////                        std::vector<double> summaries;
////                        std::vector<double> mappingstats;
////                        std::vector<double> distances;
////                        std::vector<double> weights;
//
//                        #pragma omp for
//                        for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
//
//                            int thread_id = omp_get_thread_num();
//
//                            sampler[thread_id]->sample();
//
//                            simulator[thread_id]->GetNewSimulatedCodonAlignment();
//
//                            ss[thread_id]->computeSummaries(simulator[thread_id]->CurrentLeafNodeCodonSequences);
//
////                            std::vector<double> param = lparam[thread_id]->GetCurrentParameters();
////                            std::vector<double> summaries = lparam[thread_id]->GetCurrentSummaries();
////                            std::vector<double> mappingstats = lparam[thread_id]->GetCurrentMappingStats();
////                            std::vector<double> distances  = lparam[thread_id]->GetCurrentDistances();
////                            std::vector<double> weights = lparam[thread_id]->GetCurrentWeights();
//
//                            #pragma omp critical
//                            {
//                                post->registerNewSimulation(
//                                    1,
//                                    lparam[thread_id]->GetCurrentParameters(),
//                                    lparam[thread_id]->GetCurrentSummaries(),
//                                    lparam[thread_id]->GetCurrentMappingStats(),
//                                    lparam[thread_id]->GetCurrentDistances(),
//                                    lparam[thread_id]->GetCurrentWeights()
//                                    );
//
//                                if (post->thresholdAchieved()) {
//                                    ofstream emp_os((gparam->output+".emp").c_str(),OUT);
//                                    post->writeHeader(emp_os);
//                                    post->writePosterior(emp_os);
//                                    emp_os.close();
//                                }
//
//                                if (post->Niter == post->Nrun) {
//                                    cerr << "End of the simulation process\n";
//
//                                    ofstream dist_os((gparam->output+".post").c_str(),OUT);
//                                    post->writeHeader(dist_os);
//                                    post->writePosterior(dist_os);
//                                    dist_os.close();
//
//                                    ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
//                                    post->writePosteriorPredictivePvalues(ppp_os,lparam[0]->summariesRealData);
//                                    ppp_os.close();
//
//                                    exit(0);
//
//                                } else if (post->Niter % post->threshold == 0) {
//
//                                    ofstream dist_os((gparam->output+".post").c_str(),OUT);
//                                    post->writeHeader(dist_os);
//                                    post->writePosterior(dist_os);
//                                    dist_os.close();
//
//                                    ofstream monitor_os((gparam->output+".monitor").c_str(),OUT);
//                                    post->writeMonitorPosterior(monitor_os);
//                                    monitor_os.close();
//
//                                }
//
//                            }
//                        }
//
//                    }
//
//            }
//            exit(0);
//
//
//        }
        else if (model == "FMSCpGV2" || model == "CodonMutSelFiniteABC") {
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
            for (pt_i = gparam->chainPointStart ; pt_i < gparam->chainPointEnd; pt_i+=gparam->chainPointEvery) {
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


                if (l == 0) {
//                    ostringstream ost1;
//                    ost1 <<  gparam->output << ".inputparam";
//                    ofstream lparam_os (ost1.str());
//                    lparam[l]->writeParam(lparam_os);
//                    lparam_os.close();
//

                    ostringstream ost2;
                    ost2 <<  gparam->output << ".realdata";
                    ofstream realDataSummaries_os (ost2.str());
                    lparam[l]->writeRealDataSummaries(realDataSummaries_os);
                    realDataSummaries_os.close();


                }
//

            }

            cerr << "The simulation process started\n";




            while(post3->Niter < 1000000) {
                    cerr << ".";
                    omp_set_dynamic(0);
                    omp_set_num_threads(gparam->Nthread);
                    #pragma omp parallel
                    {


                        #pragma omp for
                        for (int l = 0 ; l < Npoint; l++){



                            sampler[l]->sample();

                            simulator[l]->GetNewSimulatedCodonAlignment();

                            ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);


                            #pragma omp critical
                            {



                                if(post3->Niter < 100000) {

                                    post1->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );

                                    post2->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );

                                    post3->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );


                                    if (post3->Niter % 9999 == 0) {

//                                        ofstream dist_os3((gparam->output+"-10M.post").c_str(),OUT);
//                                        post3->writeHeader(dist_os3);
//                                        post3->writePosterior(dist_os3);
//                                        dist_os3.close();
//
//                                        ofstream monitor_os3((gparam->output+"-10M.monitor").c_str(),OUT);
//                                        post3->writeMonitorPosterior(monitor_os3);
//                                        monitor_os3.close();

                                        ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
                                        post2->writeHeader(dist_os2);
                                        post2->writePosterior(dist_os2);
                                        dist_os2.close();

                                        ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
                                        post2->writeMonitorPosterior(monitor_os2);
                                        monitor_os2.close();

                                        ofstream dist_os1((gparam->output+"-100K.post").c_str(),OUT);
                                        post1->writeHeader(dist_os1);
                                        post1->writePosterior(dist_os1);
                                        dist_os1.close();

                                        ofstream monitor_os1((gparam->output+"-100K.monitor").c_str(),OUT);
                                        post1->writeMonitorPosterior(monitor_os1);
                                        monitor_os1.close();


                                    }



                                }

                                if(post3->Niter >= 100000 && post3->Niter < 1000000 ) {

                                    post2->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );

                                    post3->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );


                                    if (post3->Niter % 9999 == 0) {

//                                        ofstream dist_os3((gparam->output+"-10M.post").c_str(),OUT);
//                                        post3->writeHeader(dist_os3);
//                                        post3->writePosterior(dist_os3);
//                                        dist_os3.close();
//
//                                        ofstream monitor_os3((gparam->output+"-10M.monitor").c_str(),OUT);
//                                        post3->writeMonitorPosterior(monitor_os3);
//                                        monitor_os3.close();

                                        ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
                                        post2->writeHeader(dist_os2);
                                        post2->writePosterior(dist_os2);
                                        dist_os2.close();

                                        ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
                                        post2->writeMonitorPosterior(monitor_os2);
                                        monitor_os2.close();

                                    }

                                }


//                                if(post3->Niter >= 1000000 && post3->Niter < 10000000 ) {
//
//                                    post3->registerNewSimulation(
//                                        lparam[l]->chainID,
//                                        lparam[l]->GetCurrentParameters(),
//                                        lparam[l]->GetCurrentSummaries(),
//                                        lparam[l]->GetCurrentMappingStats(),
//                                        lparam[l]->GetCurrentDistances(),
//                                        lparam[l]->GetCurrentWeights()
//                                        );
//
//
//                                    if (post3->Niter % 9999 == 0) {
//
//                                        ofstream dist_os3((gparam->output+"-10M.post").c_str(),OUT);
//                                        post3->writeHeader(dist_os3);
//                                        post3->writePosterior(dist_os3);
//                                        dist_os3.close();
//
//                                        ofstream monitor_os3((gparam->output+"-10M.monitor").c_str(),OUT);
//                                        post3->writeMonitorPosterior(monitor_os3);
//                                        monitor_os3.close();
//
//                                        ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
//                                        post2->writeHeader(dist_os2);
//                                        post2->writePosterior(dist_os2);
//                                        dist_os2.close();
//
//                                        ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
//                                        post2->writeMonitorPosterior(monitor_os2);
//                                        monitor_os2.close();
//
//                                    }
//
//                                }





                            }
                        }

                    }

            }

            cerr << "End of the simulation process\n";
            exit(0);

        }


//        else if (model == "MutSelAACpG") {
//            cerr << "MutSelAACpG\n";
//
//            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
//
//
//
//            Posterior* post = new Posterior(gparam);
//
//            LocalParameters** lparam = new LocalParameters*[gparam->Nthread];
//            SummaryStatistics** ss = new SummaryStatistics*[gparam->Nthread];
//            PriorSampler** sampler = new PriorSampler*[gparam->Nthread];
//            SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[gparam->Nthread];
//            AncestralSequence** ancestraseq = new AncestralSequence*[gparam->Nthread];
//            TreeSimulator** simulator = new TreeSimulator*[gparam->Nthread];
//
//            omp_set_dynamic(0);
//            omp_set_num_threads(gparam->Nthread);
//            #pragma omp parallel for
//            for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++) {
//                int thread_id = omp_get_thread_num();
//            //for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
//
//                lparam[thread_id] = new LocalParameters(gparam);
//                lparam[thread_id]->readChainCodonMutSelSBDP();
//                ss[thread_id] = new SummaryStatistics(lparam[thread_id]);
//                ss[thread_id]->computeSummaries();
//                sampler[thread_id] = new PriorSampler(lparam[thread_id]);
//                submatrix[thread_id] = new SiteInterSubMatrix(lparam[thread_id]);
//                ancestraseq[thread_id] = new AncestralSequence(lparam[thread_id]);
//                simulator[thread_id] = new TreeSimulator(lparam[thread_id],submatrix[thread_id],ancestraseq[thread_id]);
//
//
//                if (thread_id == 0) {
//
//                    ofstream lparam_os ((gparam->output+".inputparam").c_str());
//                    lparam[0]->writeParam(lparam_os);
//                    lparam_os.close();
//
//
//                    ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
//                    lparam[0]->writeRealDataSummaries(realDataSummaries_os);
//                    realDataSummaries_os.close();
//
//                    ifstream monitor_is((gparam->output+".monitor").c_str());
//                    if (!monitor_is){
//                        monitor_is.close();
//                        ofstream monitor_os((gparam->output+".monitor").c_str(),OUT);
//                        monitor_os.close();
//                    } else {
//                        post->readMonitorPosterior(monitor_is);
//                        monitor_is.close();
//                    }
//
//                }
//
//
//            }
//
//
//
//
//
//            cerr << "The simulation process started\n";
//            cerr << post->Niter << " " << post->Nrun << "\n";
//
//
//
//            while(post->Niter < post->Nrun) {
//                     cerr << ".";
//
//                    #pragma omp parallel
//                    {
//
////                        std::vector<double> param;
////                        std::vector<double> summaries;
////                        std::vector<double> mappingstats;
////                        std::vector<double> distances;
////                        std::vector<double> weights;
//
//                        #pragma omp for
//                        for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
//
//                            int thread_id = omp_get_thread_num();
//
//                            sampler[thread_id]->sample();
//
//                            simulator[thread_id]->GetNewSimulatedCodonAlignment();
//
//                            ss[thread_id]->computeSummaries(simulator[thread_id]->CurrentLeafNodeCodonSequences);
//
////                            std::vector<double> param = lparam[thread_id]->GetCurrentParameters();
////                            std::vector<double> summaries = lparam[thread_id]->GetCurrentSummaries();
////                            std::vector<double> mappingstats = lparam[thread_id]->GetCurrentMappingStats();
////                            std::vector<double> distances  = lparam[thread_id]->GetCurrentDistances();
////                            std::vector<double> weights = lparam[thread_id]->GetCurrentWeights();
//
//                            #pragma omp critical
//                            {
//                                post->registerNewSimulation(
//                                    1,
//                                    lparam[thread_id]->GetCurrentParameters(),
//                                    lparam[thread_id]->GetCurrentSummaries(),
//                                    lparam[thread_id]->GetCurrentMappingStats(),
//                                    lparam[thread_id]->GetCurrentDistances(),
//                                    lparam[thread_id]->GetCurrentWeights()
//                                    );
//
//                                if (post->thresholdAchieved()) {
//                                    ofstream emp_os((gparam->output+".emp").c_str(),OUT);
//                                    post->writeHeader(emp_os);
//                                    post->writePosterior(emp_os);
//                                    emp_os.close();
//                                }
//
//                                if (post->Niter == post->Nrun) {
//                                    cerr << "End of the simulation process\n";
//
//                                    ofstream dist_os((gparam->output+".post").c_str(),OUT);
//                                    post->writeHeader(dist_os);
//                                    post->writePosterior(dist_os);
//                                    dist_os.close();
//
//                                    ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
//                                    post->writePosteriorPredictivePvalues(ppp_os,lparam[0]->summariesRealData);
//                                    ppp_os.close();
//
//                                    exit(0);
//
//                                } else if (post->Niter % post->threshold == 0) {
//
//                                    ofstream dist_os((gparam->output+".post").c_str(),OUT);
//                                    post->writeHeader(dist_os);
//                                    post->writePosterior(dist_os);
//                                    dist_os.close();
//
//                                    ofstream monitor_os((gparam->output+".monitor").c_str(),OUT);
//                                    post->writeMonitorPosterior(monitor_os);
//                                    monitor_os.close();
//
//                                }
//
//                            }
//                        }
//
//                    }
//
//            }
//            exit(0);


//        }

        else if (model == "MutSelAACpGV2" || model == "CodonMutSelSBDPABC") {
            cerr << "CodonMutSelSBDPCpG\n";

            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
            cerr << gparam->chainPointStart << " " << gparam->chainPointEnd << " " << gparam->chainPointEvery << "\n";

            Posterior* post1 = new Posterior(gparam);
            Posterior* post2 = new Posterior(gparam);
            Posterior* post3 = new Posterior(gparam);


            int Npoint = (int) (gparam->chainPointEnd-gparam->chainPointStart)/gparam->chainPointEvery;
            cerr << "Npoint" << Npoint << "\n";
            LocalParameters** lparam = new LocalParameters*[Npoint];
            SummaryStatistics** ss = new SummaryStatistics*[Npoint];
            PriorSampler** sampler = new PriorSampler*[Npoint];
            SiteInterSubMatrix** submatrix = new SiteInterSubMatrix*[Npoint];
            AncestralSequence** ancestraseq = new AncestralSequence*[Npoint];
            TreeSimulator** simulator = new TreeSimulator*[Npoint];




            omp_set_dynamic(0);
            omp_set_num_threads(gparam->Nthread);
            int pt_i;
            #pragma omp parallel for
            for (pt_i = gparam->chainPointStart ; pt_i < gparam->chainPointEnd; pt_i+=gparam->chainPointEvery) {
                //int thread_id = omp_get_thread_num();
                //for (int thread_i = 0 ; thread_i < gparam->Nthread; thread_i++){
                int l = (int) (pt_i-gparam->chainPointStart)/gparam->chainPointEvery;
                lparam[l] = new LocalParameters(gparam);
                lparam[l]->readChainCodonMutSelSBDP(pt_i);
                ss[l] = new SummaryStatistics(lparam[l]);
                ss[l]->computeSummaries();
                sampler[l] = new PriorSampler(lparam[l]);
                submatrix[l] = new SiteInterSubMatrix(lparam[l]);
                ancestraseq[l] = new AncestralSequence(lparam[l]);
                simulator[l] = new TreeSimulator(lparam[l],submatrix[l],ancestraseq[l]);




                if (l == 0) {
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




            while(post3->Niter < 1000000) {
                    cerr << ".";
                    omp_set_dynamic(0);
                    omp_set_num_threads(gparam->Nthread);
                    #pragma omp parallel
                    {


                        #pragma omp for
                        for (int l = 0 ; l < Npoint; l++){



                            sampler[l]->sample();

                            simulator[l]->GetNewSimulatedCodonAlignment();

                            ss[l]->computeSummaries(simulator[l]->CurrentLeafNodeCodonSequences);


                            #pragma omp critical
                            {



                                if(post3->Niter < 100000) {

                                    post1->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );

                                    post2->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );

                                    post3->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );


                                    if (post3->Niter % 9999 == 0) {

//                                        ofstream dist_os3((gparam->output+"-10M.post").c_str(),OUT);
//                                        post3->writeHeader(dist_os3);
//                                        post3->writePosterior(dist_os3);
//                                        dist_os3.close();
//
//                                        ofstream monitor_os3((gparam->output+"-10M.monitor").c_str(),OUT);
//                                        post3->writeMonitorPosterior(monitor_os3);
//                                        monitor_os3.close();

                                        ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
                                        post2->writeHeader(dist_os2);
                                        post2->writePosterior(dist_os2);
                                        dist_os2.close();

                                        ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
                                        post2->writeMonitorPosterior(monitor_os2);
                                        monitor_os2.close();

                                        ofstream dist_os1((gparam->output+"-100K.post").c_str(),OUT);
                                        post1->writeHeader(dist_os1);
                                        post1->writePosterior(dist_os1);
                                        dist_os1.close();

                                        ofstream monitor_os1((gparam->output+"-100K.monitor").c_str(),OUT);
                                        post1->writeMonitorPosterior(monitor_os1);
                                        monitor_os1.close();


                                    }



                                }

                                if(post3->Niter >= 100000 && post3->Niter < 1000000 ) {

                                    post2->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );

                                    post3->registerNewSimulation(
                                        lparam[l]->MCMCpointID,
                                        lparam[l]->GetCurrentParameters(),
                                        lparam[l]->GetCurrentSummaries(),
                                        lparam[l]->GetCurrentMappingStats(),
                                        lparam[l]->GetCurrentDistances(),
                                        lparam[l]->GetCurrentWeights()
                                        );


                                    if (post3->Niter % 9999 == 0) {

//                                        ofstream dist_os3((gparam->output+"-10M.post").c_str(),OUT);
//                                        post3->writeHeader(dist_os3);
//                                        post3->writePosterior(dist_os3);
//                                        dist_os3.close();
//
//                                        ofstream monitor_os3((gparam->output+"-10M.monitor").c_str(),OUT);
//                                        post3->writeMonitorPosterior(monitor_os3);
//                                        monitor_os3.close();

                                        ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
                                        post2->writeHeader(dist_os2);
                                        post2->writePosterior(dist_os2);
                                        dist_os2.close();

                                        ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
                                        post2->writeMonitorPosterior(monitor_os2);
                                        monitor_os2.close();

                                    }

                                }


//                                if(post3->Niter >= 1000000 && post3->Niter < 10000000 ) {
//
//                                    post3->registerNewSimulation(
//                                        lparam[l]->chainID,
//                                        lparam[l]->GetCurrentParameters(),
//                                        lparam[l]->GetCurrentSummaries(),
//                                        lparam[l]->GetCurrentMappingStats(),
//                                        lparam[l]->GetCurrentDistances(),
//                                        lparam[l]->GetCurrentWeights()
//                                        );
//
//
//                                    if (post3->Niter % 9999 == 0) {
//
//                                        ofstream dist_os3((gparam->output+"-10M.post").c_str(),OUT);
//                                        post3->writeHeader(dist_os3);
//                                        post3->writePosterior(dist_os3);
//                                        dist_os3.close();
//
//                                        ofstream monitor_os3((gparam->output+"-10M.monitor").c_str(),OUT);
//                                        post3->writeMonitorPosterior(monitor_os3);
//                                        monitor_os3.close();
//
//                                        ofstream dist_os2((gparam->output+"-1M.post").c_str(),OUT);
//                                        post2->writeHeader(dist_os2);
//                                        post2->writePosterior(dist_os2);
//                                        dist_os2.close();
//
//                                        ofstream monitor_os2((gparam->output+"-1M.monitor").c_str(),OUT);
//                                        post2->writeMonitorPosterior(monitor_os2);
//                                        monitor_os2.close();
//
//                                    }
//
//                                }





                            }
                        }

                    }

            }

            cerr << "End of the simulation process\n";
            exit(0);

        }
        else if (model == "FMutSelSimu" || model == "FMutSel0Simu" || model == "MGMutSelSimu"){
            cerr << "CodeML Simu" << "\n";

            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
            Posterior* post = new Posterior(gparam);

            LocalParameters* lparam =  new LocalParameters(gparam);
            lparam->readFMutSelCodeML();

            SummaryStatistics* ss = new SummaryStatistics(lparam);
            ss->computeSummaries();
            PriorSampler* sampler = new PriorSampler(lparam);
            SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
            AncestralSequence* ancestraseq = new AncestralSequence(lparam);
            TreeSimulator* simulator = new TreeSimulator(lparam,submatrix,ancestraseq);

            ofstream lparam_os ((gparam->output+".inputparam").c_str());
            lparam->writeParam(lparam_os);
            lparam_os.close();

            ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
            lparam->writeRealDataSummaries(realDataSummaries_os);
            realDataSummaries_os.close();

            while(post->Niter < post->Nrun) {

                    sampler->sample();

                    simulator->GetNewSimulatedCodonAlignment();

                    ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

                    post->registerNewSimulation(
                                    1,
                                    lparam->GetCurrentParameters(),
                                    lparam->GetCurrentSummaries(),
                                    lparam->GetCurrentMappingStats(),
                                    lparam->GetCurrentDistances(),
                                    lparam->GetCurrentWeights()
                                    );

                    if (lparam->tofasta){
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

            ofstream dist_os((gparam->output+".post").c_str(),OUT);
            post->writeHeader(dist_os);
            post->writePosterior(dist_os);
            dist_os.close();

            ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
            post->writePosteriorPredictivePvalues(ppp_os,lparam->summariesRealData);
            ppp_os.close();
            exit(0);

        }

        else if (model == "CodonMutSelSBDP"){
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
            while(lparam->startPoint < lparam->endPoint){
                    lparam->readChainCodonMutSelSBDP();

                    int iter = 0 ;
                    while(iter < post->Nrun) {



                        sampler->sample();

                        simulator->GetNewSimulatedCodonAlignment();

                        ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

                        post->registerNewSimulation(
                                    lparam->startPoint,
                                    lparam->GetCurrentParameters(),
                                    lparam->GetCurrentSummaries(),
                                    lparam->GetCurrentMappingStats(),
                                    lparam->GetCurrentDistances(),
                                    lparam->GetCurrentWeights()
                                    );

                        if (lparam->tofasta){
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
                    while (rep < lparam->everyPoint){
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

            //ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
            //post->writePosteriorPredictivePvalues(ppp_os,lparam->summariesRealData);
            //ppp_os.close();

            exit(0);


        }
        else if (model == "CodonMutSelFinite"){
            cerr << "CodonMutSelFinite" << "\n";
            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
            Posterior* post = new Posterior(gparam);

            LocalParameters* lparam =  new LocalParameters(gparam);
            lparam->readChainCodonMutSelFinite();

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
            while(lparam->startPoint < lparam->endPoint){
                    lparam->readChainCodonMutSelFinite();
                    int iter = 0 ;
                    while(iter < post->Nrun) {



                        sampler->sample();

                        simulator->GetNewSimulatedCodonAlignment();

                        ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

                        post->registerNewSimulation(
                                    lparam->startPoint,
                                    lparam->GetCurrentParameters(),
                                    lparam->GetCurrentSummaries(),
                                    lparam->GetCurrentMappingStats(),
                                    lparam->GetCurrentDistances(),
                                    lparam->GetCurrentWeights()
                                    );

                        if (lparam->tofasta){
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
                    while (rep < lparam->everyPoint){
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

            //ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
            //post->writePosteriorPredictivePvalues(ppp_os,lparam->summariesRealData);
            //ppp_os.close();

            exit(0);



        }

//        else if (model ==  "MutSelAACpGEmpirical") {
//
//            cerr << "MutSelAACpGEmpirical" << "\n";
//
//            gparam = new GlobalParameters(controlfile);
//
//            lparam = new LocalParameters*[1];
//            ss = new SummaryStatistics*[1];
//            sampler = new PriorSampler*[1];
//            simulator = new SiteInterSubMatrix*[1];
//
//            lparam[0] = new LocalParameters(gparam);
//            ss[0] = new SummaryStatistics(lparam[0]);
//            ss[0]->computeSummaries();
//            sampler[0] = new PriorSampler(lparam[0]);
//            simulator[0] = new SiteInterSubMatrix(lparam[0]);
//
//
//            while(lparam[0]->startPoint < lparam[0]->endPoint) {
//                int rep = 0;
//                while (rep < lparam[0]->Nrep) {
//                    lparam[0]->readChainMutSelAAC();
//                    lparam[0]->SetTree();
//                    sampler[0]->sample();
//                    simulator[0]->GetNewStationaryCodonSequence();
//                    simulator[0]->GetNewAncestralCodonSequence();
//                    simulator[0]->GetNewSimulatedCodonAlignment();
//                    ss[0]->computeSummaries(simulator[0]->CurrentLeafNodeCodonSequences);
//                    lparam[0]->registerNewSimulation();
//
//
//                    if (lparam[0]->tofasta) {
//                        ostringstream oss;
//                        oss << lparam[0]->output << "_" << lparam[0]->startPoint << "." << rep << ".fasta";
//                        ofstream fasta_os(oss.str(),OUT);
//                        lparam[0]->toAli(fasta_os,simulator[0]->CurrentLeafNodeCodonSequences);
//                        fasta_os.close();
//                        oss.str("");
//                        oss << lparam[0]->output << "_" << lparam[0]->startPoint << "." << rep << ".tre";
//                        ofstream tree_os(oss.str(),OUT);
//                        lparam[0]->refTree->Print(tree_os);
//                        tree_os.close();
//
//                    }
//                    rep++;
//                    cerr << ".";
//                }
//
//                lparam[0]->startPoint +=lparam[0]->everyPoint;
//
//            }
//
//            ofstream dist_os((gparam->output+".post").c_str(),OUT);
//            gparam->writeHeader(dist_os);
//            gparam->writePosterior(dist_os);
//            dist_os.close();
//
//            ofstream emp_os((gparam->output+".emp").c_str(),OUT);
//            gparam->writeHeader(emp_os);
//            gparam->writePosterior(emp_os);
//            emp_os.close();
//
//            ofstream ppp_os((gparam->output+".empppp").c_str(),OUT);
//            gparam->writePosteriorPredictivePvalues(ppp_os);
//            ppp_os.close();
//
//
//         else if (model == "MutSelAACpGppchecks" || model == "FMutSelppchecks") {
//            // the chain point is defined within the conf file
//
//            if (model == "MutSelAACpGppchecks") {
//
//                cerr << model  << "\n";
//
//            } else if (model == "FMutSelppchecks") {
//
//                cerr << model  << "\n";
//
//            }
//
//
//            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
//            Posterior* post = new Posterior(gparam);
//
//            LocalParameters* lparam =  new LocalParameters(gparam);
//            if (model == "MutSelAACpGppchecks") {
//
//                lparam->readChainCodonMutSelSBDP();
//
//            } else if (model == "FMutSelppchecks") {
//
//                lparam->readChainCodonMutSelFinite();
//
//            }
//
//
//            ofstream lparam_os ((gparam->output+".inputparam").c_str());
//            lparam->writeParam(lparam_os);
//            lparam_os.close();
//
//            SummaryStatistics* ss = new SummaryStatistics(lparam);
//            ss->computeSummaries();
//
//            ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
//            lparam->writeRealDataSummaries(realDataSummaries_os);
//            realDataSummaries_os.close();
//
//            PriorSampler* sampler = new PriorSampler(lparam);
//            SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
//            AncestralSequence* ancestraseq = new AncestralSequence(lparam);
//            TreeSimulator* simulator = new TreeSimulator(lparam,submatrix,ancestraseq);
//
//            post->readPosterior(lparam->posteriorfile);
//
//            cerr << "The simulation process started\n";
//
//            if (!post->posterior.empty()) {
//                cerr << post->posterior.size() << "\n";
//                unsigned int it = 0 ;
//                while(it < 1000) {
//                    lparam->SetCurrentParametersFromPosterior(post->posterior,0);
//                    //cerr << "AAA1\n";
//                    simulator->GetNewSimulatedCodonAlignment();
//
//                    ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);
//
//                    post->registerNewSimulation(
//                                it,
//                                lparam->GetCurrentParameters(),
//                                lparam->GetCurrentSummaries(),
//                                lparam->GetCurrentMappingStats(),
//                                lparam->GetCurrentDistances(),
//                                lparam->GetCurrentWeights()
//                                );
//                    it++;
//                    cerr << ".";
//                }
//
//            cerr << "End of the simulation process\n";
//
//            ofstream dist_os((gparam->output+".post").c_str(),OUT);
//            post->writeHeader(dist_os);
//            post->writePosterior(dist_os);
//            dist_os.close();
//
////            ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
////            post->writePosteriorPredictivePvalues(ppp_os);
////            ppp_os.close();
//
//
//            }
//
//        }

        else if (model == "MutSelAACpGppchecksV2" || model == "FMutSelppchecksV2" || model == "MAP100CodonMutSelFinitePPChecks" || model == "MAP100CodonMutSelSBDPPPChecks") {

            // the chain pointS are extract from the posterior file according to chainID
            cerr << model  << "\n";


            GlobalParameters* gparam = new GlobalParameters(model, controlfile);
            Posterior* post = new Posterior(gparam);

            LocalParameters* lparam =  new LocalParameters(gparam);

            if (model == "MutSelAACpGppchecksV2" || model == "MAP100CodonMutSelSBDPPPChecks") {

                lparam->readChainCodonMutSelSBDP();

            } else if (model == "FMutSelppchecksV2" || model == "MAP100CodonMutSelFinitePPChecks") {

                lparam->readChainCodonMutSelFinite();

            }


            ofstream lparam_os ((gparam->output+".inputparam").c_str());
            lparam->writeParam(lparam_os);
            lparam_os.close();

            SummaryStatistics* ss = new SummaryStatistics(lparam);
            ss->computeSummaries();

            ofstream realDataSummaries_os ((gparam->output+".realdata").c_str());
            lparam->writeRealDataSummaries(realDataSummaries_os);
            realDataSummaries_os.close();

            PriorSampler* sampler = new PriorSampler(lparam);
            SiteInterSubMatrix* submatrix = new SiteInterSubMatrix(lparam);
            AncestralSequence* ancestraseq = new AncestralSequence(lparam);
            TreeSimulator* simulator = new TreeSimulator(lparam,submatrix,ancestraseq);

            post->readPosterior(lparam->posteriorfile);

            cerr << "The simulation process started\n";

            if (!post->posterior.empty()) {

                unsigned int it = 0 ;
                while(it < 1000) {

                    int point = static_cast<int> (lparam->rnd->Uniform() * 99);

                    lparam->SetCurrentParametersFromPosterior(post->posterior,point);

                    if (model == "MutSelAACpGppchecksV2" || model == "MAP100CodonMutSelSBDPPPChecks") {

                        lparam->readChainCodonMutSelSBDP(lparam->GetPointID());

                    } else if (model == "FMutSelppchecksV2" || model == "MAP100CodonMutSelFinitePPChecks") {

                        lparam->readChainCodonMutSelFinite(lparam->GetPointID());

                    }



                    simulator->GetNewSimulatedCodonAlignment();

                    ss->computeSummaries(simulator->CurrentLeafNodeCodonSequences);

                    post->registerNewSimulation(
                                point,
                                lparam->GetCurrentParameters(),
                                lparam->GetCurrentSummaries(),
                                lparam->GetCurrentMappingStats(),
                                lparam->GetCurrentDistances(),
                                lparam->GetCurrentWeights()
                                );
                    it++;
                    cerr << ".";
                }

            cerr << "End of the simulation process\n";

            ofstream dist_os((gparam->output+".post").c_str(),OUT);
            post->writeHeader(dist_os);
            post->writePosterior(dist_os);
            dist_os.close();

//            ofstream ppp_os((gparam->output+".ppp").c_str(),OUT);
//            post->writePosteriorPredictivePvalues(ppp_os);
//            ppp_os.close();


            }

        }
        // end main

}








