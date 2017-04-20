#include "Posterior.h"

Posterior::Posterior(GlobalParameters* gparam)
{


    this->Niter = gparam->Niter;
    this->Nrun = gparam->Nrun;
    this->threshold = gparam->threshold;

    this->localcontrolfile = gparam->localcontrolfile;
    this->output = gparam->output;
    this->model = gparam->model;

    this->TOOSMALL = gparam->TOOSMALL;
    this->TOOLARGE = gparam->TOOLARGE;
    this->TOOLARGENEGATIVE = gparam->TOOLARGENEGATIVE;

    this->NSummaries = gparam->NSummaries;
    this->NParam = gparam->NParam;
    this->NMapStats = gparam->NMapStats;

    this->listParam = new string[this->NParam];
    for (int param_i = 0; param_i < this->NParam; param_i++){
        this->listParam[param_i] = gparam->listParam[param_i];

    }



    this->listSummaries = new string[this->NSummaries];
    for (int summary_i = 0; summary_i < this->NSummaries; summary_i++){
        this->listSummaries[summary_i] = gparam->listSummaries[summary_i];

    }


    this->listMapStats = new string[this->NMapStats];
    for (int MapStats_i = 0; MapStats_i < this->NMapStats; MapStats_i++){
        this->listMapStats[MapStats_i] = gparam->listMapStats[MapStats_i];

    }


    this->NusedMapStats = gparam->NusedMapStats;
    this->NusedMapAncStats = gparam->NusedMapAncStats;
    this->NusedParam = gparam->NusedParam;
    this->NusedSummaries = gparam->NusedSummaries;
    this->Ngenes = gparam->Ngenes;


    this->mapUsedParam.insert(gparam->mapUsedParam.begin(),gparam->mapUsedParam.end());
    this->mapUsedSummaries.insert(gparam->mapUsedSummaries.begin(),gparam->mapUsedSummaries.end());
    this->mapUsedMapStats.insert(gparam->mapUsedMapStats.begin(),gparam->mapUsedMapStats.end());
    this->mapUsedMapAncStats.insert(gparam->mapUsedMapAncStats.begin(),gparam->mapUsedMapAncStats.end());


    this->sorted = false;
    this->Naccepted = 0;

    this->randomseed = -1;
    this->rnd = new Random(randomseed);



    population_t.reserve(this->threshold);

    empVar = new double [this->NusedParam];
    empMean = new double [this->NusedParam];
    for (unsigned int i = 0 ; i < this->NusedParam; i++) {
        empVar[i] = 1.0;
        empMean[i] = 1.0;
    }

}

Posterior::~Posterior()
{
    //dtor
}



std::vector<std::vector<double>> Posterior::GetPartialDistances(){

    std::vector<std::vector<double>> X;

    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ ){
        std::vector<double> X_i;
        for (unsigned int distance_i = 0 ; distance_i < NusedSummaries+1; distance_i++){
            if (distance_i == 0) {
                X_i.push_back(1.0);
            } else {
                X_i.push_back(std::get<4>(population_t[simu_i])[distance_i]);
            }

        }
        X.push_back(X_i);
    }

    return X;

}

std::vector<std::vector<double>> Posterior::GetPartialDistancesT(){

    std::vector<std::vector<double>> X_T;


        for (int distance_i = 0 ; distance_i < NusedSummaries+1; distance_i++){
            std::vector<double> X_Ti;
            for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ ){
                if (distance_i == 0) {
                    X_Ti[simu_i] = 1.0;
                } else {
                    X_Ti.push_back(std::get<4>(population_t[simu_i])[distance_i]);
                }


            }
            X_T.push_back(X_Ti);
        }

    return X_T;

}



std::vector<std::vector<double>> Posterior::GetTheta() {

    std::vector<std::vector<double> > theta;

    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ ){

            theta.push_back(std::get<1>(population_t[simu_i]));



    }

    return theta;

}




double Posterior::GetEpanechnikov(double x, double y){

    return 1-(x*x)/(y*y);

}

std::vector<double> Posterior::GetWeights(){
    std::vector<double> W;

    double max_ = std::get<4>(population_t[population_t.size()-1])[NusedSummaries];


    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ ){
        W.push_back(GetEpanechnikov(std::get<4>(population_t[simu_i])[NusedSummaries],max_));

    }

    return W;
}


std::vector<std::vector<double>> Posterior::GetLocalWeights(){
    std::vector<std::vector<double>> W;

    double* max_ = new double[NusedSummaries];
    for (unsigned int distance_i = 0 ; distance_i < NusedSummaries; distance_i++){
        max_[distance_i] = std::get<4>(population_t[population_t.size()-1])[distance_i];
    }

    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ ){
        std::vector<double> W_i;
        for (unsigned int distance_i = 0 ; distance_i < NusedSummaries; distance_i++){
            W_i.push_back(GetEpanechnikov(std::get<4>(population_t[simu_i])[distance_i],max_[distance_i]));
        }

        W.push_back(W_i);
    }

    return W;
}


void Posterior::writePosterior(ofstream&os) {

    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ ){

        //write chainID
        os << std::get<0>(population_t[simu_i]) << "\t";

        //write parameters
        if (this->NusedParam > 0) {
            for (unsigned int param_i = 0 ; param_i < std::get<1>(population_t[simu_i]).size(); param_i++){
                os << std::get<1>(population_t[simu_i])[param_i] << "\t";
            }
        }
        //write summaries
        if (this->NusedSummaries > 0) {
            for (unsigned int summary_i = 0 ; summary_i < std::get<2>(population_t[simu_i]).size(); summary_i++){
                os << std::get<2>(population_t[simu_i])[summary_i] << "\t";
            }
        }

        //write distances (sum of square discrepancies)
        if (this->NusedSummaries > 0) {
            for (unsigned int distance_i = 0 ; distance_i < std::get<4>(population_t[simu_i]).size(); distance_i++){
                os << std::get<4>(population_t[simu_i])[distance_i] << "\t";
            }
        }

        //write mappingstats
        if (this->NusedMapStats > 0) {
            for (unsigned int mapping_i = 0 ; mapping_i < std::get<3>(population_t[simu_i]).size(); mapping_i++){
                if(mapping_i < std::get<3>(population_t[simu_i]).size()-1){
                    os << std::get<3>(population_t[simu_i])[mapping_i] << "\t";
                } else {
                    os << std::get<3>(population_t[simu_i])[mapping_i] << "\t";
                }
            }
        }
        os <<  "\n";
    }
}



void Posterior::readPosterior(string posteriorfile){
    ifstream is(posteriorfile.c_str());
    if (!is)       {
        cerr << "error: did not find " << posteriorfile << "\n";
        exit(1);
    }
    string line;
    std::getline(is, line);
    istringstream iss(line);
    string w;
    int k = 0;
    while(iss >> w){

        auto it = mapUsedParam.find(w);
        if (it != mapUsedParam.end()) {
            cerr << it->second << " " << k << "\n";
            it->second = k;
            k++;

        } else {

            cerr << "Undefined parameter " << w << "\n";
            exit(0);

        }

    }

    NusedParam = k;
    while(std::getline(is, line)) {
        if(!line.empty()) {
            //cerr << line << "\n";
            istringstream iss(line);
            double w;
            std::vector <double> cur_param;

            for (unsigned int param_i = 0; param_i < NusedParam; param_i++) {
                    iss >> w;
                    cur_param.push_back(w);
            }

            posterior.push_back(cur_param);
        }
    }
    is.close();
}

void Posterior::writePosteriorPredictivePvalues(ofstream& os, std::vector<double>summariesRealData){


    unsigned int pop_size = population_t.size();
    int k =0 ;
    if (NusedSummaries > 0) {
        string* arrSummaries = new string[this->NusedSummaries];
        for(unsigned int summary_i = 0 ; summary_i < this->NSummaries; summary_i++){
            auto it = mapUsedSummaries.find(this->listSummaries[summary_i]);
            if(it != mapUsedSummaries.end() ) {
                if(it->second != -1) {
                    arrSummaries[it->second] = it->first;

                }

            }

        }


        for(unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
            if(k == 0 ){
                os << arrSummaries[summary_i];
                k = 1;
            } else {
                os << "\t" << arrSummaries[summary_i];
            }

        }
        os << "\n";
    }





    double* ppp  = new double[this->NusedSummaries];
    double* mean = new double[this->NusedSummaries];
    double* var  = new double[this->NusedSummaries];
    double* sum1  = new double[this->NusedSummaries];
    double* sum2  = new double[this->NusedSummaries];
    double* sum3  = new double[this->NusedSummaries];

    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
        ppp[summary_i] = 0.0;
        mean[summary_i] = 0.0;
        var[summary_i] = 0.0;
        sum1[summary_i] = 0.0;
        sum2[summary_i] = 0.0;
        sum3[summary_i] = 0.0;
    }

    for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++ ){
        for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
            if(std::get<2>(population_t[simu_i])[summary_i] < summariesRealData[summary_i]){
                ppp[summary_i]+=1;
            }
            sum1[summary_i]+=exp2(std::get<2>(population_t[simu_i])[summary_i]);
        }
     }


    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
        mean[summary_i] = sum1[summary_i]/pop_size;

    }

    for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++ ){
        for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
            sum2[summary_i] += (exp2(std::get<2>(population_t[simu_i])[summary_i]-mean[summary_i])*exp2(std::get<2>(population_t[simu_i])[summary_i]-mean[summary_i]));
            sum3[summary_i] += exp2(std::get<2>(population_t[simu_i])[summary_i]-mean[summary_i]);
        }
    }

    //posterior predictive p-values (ppp)
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
        var[summary_i] = (sum2[summary_i]-(sum3[summary_i]*sum3[summary_i])/pop_size)/(pop_size-1);
        ppp[summary_i] /= pop_size;
        if(summary_i < this->NusedSummaries-1) {
            os << ppp[summary_i] << "\t";
        } else {
            os << ppp[summary_i] << "\n";
        }

    }

    //means
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
        if(summary_i < this->NusedSummaries-1) {
            os << mean[summary_i] << "\t";
        } else {
            os << mean[summary_i] << "\n";
        }
    }

    //variances
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
        if(summary_i < this->NusedSummaries-1) {
            os << var[summary_i] << "\t";
        } else {
            os << var[summary_i] << "\n";
        }
    }

    //z-scores
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++){
        if(summary_i < this->NusedSummaries-1) {
            os <<  (exp2(summariesRealData[summary_i])-mean[summary_i] / sqrt(var[summary_i])) << "\t";
        } else {
            os <<  (exp2(summariesRealData[summary_i])-mean[summary_i] / sqrt(var[summary_i])) << "\n";
        }
    }



    delete [] ppp;
    delete [] mean;
    delete [] var;
    delete [] sum1;
    delete [] sum2;
    delete [] sum3;


}


void Posterior::readMonitorPosterior(ifstream & is) {
    string line;
    std::getline(is, line);
    if(!line.empty()) {
        istringstream iss(line);
        iss >> this->Niter;

    } else {
        this->Niter = 0;
        cerr << "Monitor file is empty " << "\n";
    }

    std::getline(is, line);
    if(!line.empty()) {
        istringstream iss(line);
        iss >> this->Naccepted;

    } else {
        this->Naccepted = 0;
        cerr << "Monitor file is empty " << "\n";
    }
}

void Posterior::writeMonitorPosterior(ofstream& os){
    os << this->Niter << "\n";
    os << this->Naccepted << "\n";
    os << GetAcceptanceRate() << "\n";

}

double Posterior::GetAcceptanceRate(){
    return (double) Naccepted/Niter;
}




void Posterior::writeHeader(ofstream&os){

    // write parameters' header
    //for (unsigned int param_i = 0 ; param_i < listUsedParam.size() ; param_i++){
    int k = 0;




    if (NusedParam > 0) {
        string* arrParam = new string[NusedParam];
        for (unsigned int param_i = 0 ; param_i < NParam ; param_i++){
            auto it = mapUsedParam.find(listParam[param_i]);
            if(it != mapUsedParam.end() ) {
                if(it->second != -1) {
                    arrParam[it->second] = it->first;

                }

            }
            //os << listParam[listUsedParam[param_i]] << "\t";

        }

        for (unsigned int param_i = 0 ; param_i < NusedParam ; param_i++){
            if(k == 0) {
                os << arrParam[param_i];
                k = 1;
            } else {
                os << "\t" << arrParam[param_i];
            }
        }
        delete [] arrParam;
    }

    //write summaires' header
    //for(unsigned int summary_i = 0 ; summary_i < listUsedSummaries.size(); summary_i++){
    if (NusedSummaries > 0) {
        string* arrSummaries = new string[NusedSummaries];
        for(unsigned int summary_i = 0 ; summary_i < NSummaries; summary_i++){
            auto it = mapUsedSummaries.find(listSummaries[summary_i]);
            if(it != mapUsedSummaries.end() ) {
                if(it->second != -1) {
                    arrSummaries[it->second] = it->first;

                }

            }


            //os << listSummaries[listUsedSummaries[summary_i]] << "\t";

        }


        for(unsigned int summary_i = 0 ; summary_i < NusedSummaries; summary_i++){
            if(k == 0 ){
                os << arrSummaries[summary_i];
                k = 1;
            } else {
                os << "\t" << arrSummaries[summary_i];
            }

        }

        //write distance header
//    for(unsigned int summary_i = 0 ; summary_i < listUsedSummaries.size(); summary_i++){
//        os << "D_" << listSummaries[listUsedSummaries[summary_i]]  <<"\t";
//    }

        for(unsigned int summary_i = 0 ; summary_i < NusedSummaries; summary_i++){
            if (k == 0) {
                os << "D_" << arrSummaries[summary_i];
                k = 1;
            } else {
                os << "\t" << "D_" << arrSummaries[summary_i];
            }
        }

        os << "\tD_sum";
        delete [] arrSummaries;


    }
    //write mappingstats

    if (NusedMapAncStats > 0) {
        string* arrMapAnc = new string[NusedMapAncStats];
        for(unsigned int map_i = 0 ; map_i < NMapStats; map_i++){
            auto it = mapUsedMapAncStats.find(listMapStats[map_i]);
            if(it != mapUsedMapAncStats.end()) {
                if(it->second != -1) {
                    arrMapAnc[it->second] = it->first;

                }

            }
        }


        for(unsigned int map_i = 0; map_i < NusedMapAncStats; map_i++){
            if (k == 0) {
                os << "Manc_" << arrMapAnc[map_i];
                k = 1;
            } else {
                os << "\t" << "Manc_" << arrMapAnc[map_i];
            }


        }
        delete [] arrMapAnc;
    }

    if (NusedMapStats > 0) {
        string* arrMap = new string[NusedMapStats];
        for(unsigned int map_i = 0 ; map_i < NMapStats; map_i++){
            auto it = mapUsedMapStats.find(listMapStats[map_i]);
            if(it != mapUsedMapStats.end()) {
                if(it->second != -1) {
                    arrMap[it->second] = it->first;

                }
            }
        }


        for(unsigned int map_i = 0; map_i < NusedMapStats; map_i++){
            if (k == 0) {
                os << "M_" << arrMap[map_i];
                k = 1;
            } else {
                os << "\t" << "M_" << arrMap[map_i];
            }

        }
        delete [] arrMap;


    }

    os << "\n";

}

void Posterior::GetWeights(string kernel) {
        unsigned int pop_size = population_t.size();    if (kernel == "sNormal") {

        for (unsigned int param_i = 0 ; param_i < NusedParam ; param_i++){
            double new_weight = 0.0;
            for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++) {
                new_weight += std::get<5>(population_t[simu_i])[param_i]* (std::get<1>(population_t[simu_i])[param_i]  + 2.0 * (empVar[param_i]) * rnd->sNormal());

            }
            new_weight /= pop_size;
//            std::get<4>(population_t[simu_i])[param_i] = new_weight;
        }



    }


}

void Posterior::GetEmpVar() {

    unsigned int pop_size = population_t.size();


    double* mean = new double[NusedParam];
    double* sum1  = new double[NusedParam];
    double* sum2  = new double[NusedParam];
    double* sum3  = new double[NusedParam];

    for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++){
        mean[param_i] = 0.0;
        sum1[param_i] = 0.0;
        sum2[param_i] = 0.0;
        sum3[param_i] = 0.0;
    }

    for (unsigned int simu_i = 0 ; simu_i < NusedParam; simu_i++ ){
        for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++){
            sum1[param_i] += std::get<1>(population_t[simu_i])[param_i];
        }
     }


    for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++){
        mean[param_i] = sum1[param_i]/pop_size;
        empMean[param_i] = mean[param_i];
    }

    for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++ ){
        for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++){
            //(x - K) * (x - K)
            sum2[param_i] += (std::get<1>(population_t[simu_i])[param_i]-mean[param_i])*(std::get<1>(population_t[simu_i])[param_i]-mean[param_i]);
            //x - K
            sum3[param_i] += std::get<1>(population_t[simu_i])[param_i]-mean[param_i];
        }
    }


    for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++){
        // sum_sqr - (sum_ * sum_)/n)/(n - 1)
        empVar[param_i] = (sum2[param_i]-(sum3[param_i]*sum3[param_i])/pop_size)/(pop_size-1);
    }

    delete [] mean;
    delete [] sum1;
    delete [] sum2;
    delete [] sum3;

}



void Posterior::registerNewSimulation(int chainID, std::vector<double> param ,std::vector<double> summaries,std::vector<double> mappingstats ,std::vector<double> distances, std::vector<double> weights){

    if(population_t.empty()) {

        //cerr << "POPULATION IS EMPTY" << Naccepted <<  " "<< threshold <<  " \n";

        population_t.push_back(make_tuple(chainID,param,summaries,mappingstats,distances,weights));
        Naccepted++;

    } else if(population_t.size() < (unsigned) threshold){

        //cerr << "POPULATION IS LESS THAN THRESHOLD" << " " << population_t.size() << " "<< Naccepted  << " "<< threshold << " \n";

        population_t.push_back(make_tuple(chainID,param,summaries,mappingstats,distances,weights));
        Naccepted++;

    } else {

        if (!sorted) {

            std::sort(population_t.begin(), population_t.end(),
                [](const std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &left,
                   const std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &right)
                    {
                        return std::get<4>(left).back() < std::get<4>(right).back();
                    }
                );
            sorted = true;
        }

        std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> cur_tuple
        = make_tuple(chainID,param,summaries,mappingstats,distances,weights);

        auto it = std::lower_bound(population_t.begin(), population_t.end(), cur_tuple,
        [](const std::tuple<int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> &left,
           const std::tuple<int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> &right)
            {
                return std::get<4>(left).back() < std::get<4>(right).back();
            }
        );

        // compute acceptance rate
        if (it != population_t.end()) {
            //cerr << "INSERTED" << " " << population_t.size()  << " "<<  Naccepted << " "<< threshold << " \n";
            population_t.insert(it,cur_tuple);
            population_t.pop_back();
            population_t.shrink_to_fit();
            Naccepted++;
        }

    }
    Niter++;
}
