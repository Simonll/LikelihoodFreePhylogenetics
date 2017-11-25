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

    static const int chainIDGetter =  0;
    static const int paramGetter = 1;
    static const int summariesGetter = 2;
    static const int accsummariesGetter = 3;
    static const int evoancstatsGetter = 4;
    static const int evostatsGetter = 5;
    static const int ssevostatsGetter = 6;
    static const int distancesGetter = 7;
    static const int weightsGetter = 8;


    double TOOSMALL;
    double TOOLARGE;
    double TOOLARGENEGATIVE;

    int verbose;
    int NSummaries;
    int NParam;
    int NEvoStats;
    int NSiteSpecificEvoStats;
    int NAccessorySummaries;
    int NEvoAncStats;

    string* listParam;
    string* listSummaries;
    string* listEvoStats;
    string* listSiteSpecificEvoStats;

    std::map<string,int> mapUsedParam;
    std::map<string,int> mapUsedSummaries;
    std::map<string,int> mapUsedAccessorySummaries;
    std::map<string,int> mapUsedEvoStats;
    std::map<string,int> mapUsedSiteSpecificEvoStats;
    std::map<string,int> mapUsedEvoAncStats;

    int NusedEvoStats;
    int NusedSiteSpecificEvoStats;
    int NusedEvoAncStats;
    int NusedParam;
    int NusedSummaries;
    int NusedAccessorySummaries;
    int Ngenes;

    string localcontrolfile, output, model;

    std::vector<std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>> population_t;
    std::vector<std::vector<double>> posterior;
    double*empVar;
    double*empMean;

    Random* rnd;
    int randomseed;

    int OutPartialDistance, Niter, Nrun, Naccepted, threshold, Nthread, Nsite_codon;
    bool sorted;
    // writter
    void writePosterior(ofstream& os);
    void writeMonitorPosterior(ofstream& os);
    void writePosteriorPredictivePvalues(ofstream& os, std::vector<double> realDataSummaries);
    void writePosteriorPredictiveStatistics(ofstream& os, std::vector<double> realDataSummaries);
    void writeHeader(ofstream& os);


    //readers
    void readPosterior(string posteriorfile);
    void readPosterior(ifstream& is);
    void readMonitor();
    void readMonitor(ifstream& is);

    //Getters
    int PosteriorGetSize();
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
    void sortPopulation();
    void slaveToMaster(std::vector<std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>> population_i);
    //Setters

    void slaveRegisterNewSimulation(
        int chainID,
        std::vector<double> param,
        std::vector<double> summaries,
        std::vector<double> accsummaries,
        std::vector<double> ancevostat,
        std::vector<double> evostat,
        std::vector<double> ssevostat,
        std::vector<double> distances,
        std::vector<double> weights
    );

    void registerNewSimulation(
        int chainID,
        std::vector<double> param,
        std::vector<double> summaries,
        std::vector<double> accsummaries,
        std::vector<double> ancevostat,
        std::vector<double> evostat,
        std::vector<double> ssevostat,
        std::vector<double> distances,
        std::vector<double> weights
    );

    void registerOldSimulation(
        int chainID,
        std::vector<double> param,
        std::vector<double> summaries,
        std::vector<double> accsummaries,
        std::vector<double> ancevostat,
        std::vector<double> evostat,
        std::vector<double> ssevostat,
        std::vector<double> distances,
        std::vector<double> weights
    );

    void SetNsite(int i);




    int GetSize()
    {
        return population_t.size();
    }
    bool thresholdAchieved()
    {
        bool test = false;
        if (this->Niter == threshold)
        {
            test = true;
        }
        return test;
    }

protected:
private:
};

#endif // POSTERIOR_H
