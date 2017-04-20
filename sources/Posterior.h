#ifndef POSTERIOR_H
#define POSTERIOR_H

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
#include "GlobalParameters.h"


class Posterior
{
    public:

        Posterior(GlobalParameters* gparam);
        virtual ~Posterior();


        double TOOSMALL;
        double TOOLARGE;
        double TOOLARGENEGATIVE;

        int NSummaries;
        int NParam;
        int NMapStats;

        string* listParam;
        string* listSummaries;
        string* listMapStats;

        std::map<string,int> mapUsedParam;
        std::map<string,int> mapUsedSummaries;
        std::map<string,int> mapUsedMapStats;
        std::map<string,int> mapUsedMapAncStats;

        int NusedMapStats;
        int NusedMapAncStats;
        int NusedParam;
        int NusedSummaries;
        int Ngenes;

        string localcontrolfile, output, model;

        std::vector<std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>> population_t;
        std::vector<std::vector<double>> posterior;
        double*empVar;
        double*empMean;

        Random* rnd;
        int randomseed;

        int Niter, Nrun, Naccepted, threshold, Nthread;
        bool sorted;
        // writter
        void writePosterior(ofstream&os);
        void writeMonitorPosterior(ofstream& os);
        void writePosteriorPredictivePvalues(ofstream& os, std::vector<double> realDataSummaries);
        void writeHeader(ofstream&os);


        //readers
        void readPosterior(string posteriorfile);
        void readMonitorPosterior();
        void readMonitorPosterior(ifstream & is);

        //Getters
        double GetAcceptanceRate();
        std::vector<std::vector<double>> GetPartialDistances();
        std::vector<std::vector<double>> GetPartialDistancesT();
        std::vector<double> GetTheta_i(int theta_i);
        std::vector<std::vector<double>> GetTheta();
        void GetEmpVar();
        void GetWeights(string kernel);
        std::vector<std::vector<double>> GetLocalWeights();
        std::vector<double> GetWeights();
        double GetEpanechnikov(double x, double y);


        //Setters
        void registerNewSimulation(int chainID, std::vector<double> param ,std::vector<double> summaries,std::vector<double> mappingstats, std::vector<double> distances,std::vector<double> weights);


        void incNiter(){
            this->Niter++;
        }


        int GetSize(){return population_t.size();}
        bool thresholdAchieved(){
            bool test = false;
            if (this->Niter == threshold){
                test = true;
            }
            return test;
        }


    protected:

    private:
};

#endif // POSTERIOR_H
