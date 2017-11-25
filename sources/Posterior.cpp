#include "Posterior.h"

Posterior::Posterior(GlobalParameters* gparam)
{

    this->verbose = gparam->verbose;
    this->Niter = gparam->Niter;
    this->Nrun = gparam->Nrun;
    this->threshold = gparam->threshold;

    this->localcontrolfile = gparam->localcontrolfile;
    this->output = gparam->output;
    this->model = gparam->model;
    this->OutPartialDistance = gparam->OutPartialDistance;

    this->TOOSMALL = gparam->TOOSMALL;
    this->TOOLARGE = gparam->TOOLARGE;
    this->TOOLARGENEGATIVE = gparam->TOOLARGENEGATIVE;

    this->NSummaries = gparam->NSummaries;
    this->NParam = gparam->NParam;
    this->NEvoStats = gparam->NEvoStats;
    this->NSiteSpecificEvoStats = gparam->NSiteSpecificEvoStats;

    if(verbose)
    {
        cerr << "Posterior1\n";
    }

    this->listParam = new string[this->NParam];
    for (int param_i = 0; param_i < this->NParam; param_i++)
    {
        this->listParam[param_i] = gparam->listParam[param_i];

    }

    if(verbose)
    {
        cerr << "Posterior2\n";
    }

    this->listSummaries = new string[this->NSummaries];
    for (int summary_i = 0; summary_i < this->NSummaries; summary_i++)
    {
        this->listSummaries[summary_i] = gparam->listSummaries[summary_i];

    }

    if(verbose)
    {
        cerr << "Posterior3\n";
    }


    this->listEvoStats = new string[this->NEvoStats];
    for (int EvoStats_i = 0; EvoStats_i < this->NEvoStats; EvoStats_i++)
    {
        this->listEvoStats[EvoStats_i] = gparam->listEvoStats[EvoStats_i];

    }

    if(verbose)
    {
        cerr << "Posterior4\n";
    }

    this->listSiteSpecificEvoStats = new string[this->NSiteSpecificEvoStats];
    for (int EvoStats_i = 0; EvoStats_i < this->NSiteSpecificEvoStats; EvoStats_i++)
    {
        this->listSiteSpecificEvoStats[EvoStats_i] = gparam->listSiteSpecificEvoStats[EvoStats_i];

    }

    if(verbose)
    {
        cerr << "Posterior5\n";
    }

    this->NusedEvoStats = gparam->NusedEvoStats;
    this->NusedSiteSpecificEvoStats = gparam->NusedSiteSpecificEvoStats;
    this->NusedEvoAncStats = gparam->NusedEvoAncStats;
    this->NusedParam = gparam->NusedParam;
    this->NusedSummaries = gparam->NusedSummaries;
    this->NusedAccessorySummaries = gparam->NusedAccessorySummaries;
    this->Ngenes = gparam->Ngenes;

    if(verbose)
    {
        cerr << "Posterior6\n";
    }

    this->mapUsedParam.insert(gparam->mapUsedParam.begin(),gparam->mapUsedParam.end());
    this->mapUsedSummaries.insert(gparam->mapUsedSummaries.begin(),gparam->mapUsedSummaries.end());
    this->mapUsedAccessorySummaries.insert(gparam->mapUsedAccessorySummaries.begin(),gparam->mapUsedAccessorySummaries.end());
    this->mapUsedEvoStats.insert(gparam->mapUsedEvoStats.begin(),gparam->mapUsedEvoStats.end());
    this->mapUsedSiteSpecificEvoStats.insert(gparam->mapUsedSiteSpecificEvoStats.begin(),gparam->mapUsedSiteSpecificEvoStats.end());
    this->mapUsedEvoAncStats.insert(gparam->mapUsedEvoAncStats.begin(),gparam->mapUsedEvoAncStats.end());


    if(verbose)
    {
        cerr << "Posterior7\n";
    }

    this->sorted = false;
    this->Naccepted = 0;
    this->randomseed = -1;
    this->rnd = new Random(randomseed);



    population_t.reserve(this->Nrun);

    empVar = new double [this->NusedParam];
    empMean = new double [this->NusedParam];
    for (unsigned int i = 0 ; i < this->NusedParam; i++)
    {
        empVar[i] = 1.0;
        empMean[i] = 1.0;
    }

}

Posterior::~Posterior()
{
    //dtor
}

int Posterior::PosteriorGetSize()
{

    return posterior.size();
}


std::vector<std::vector<double>> Posterior::GetPartialDistances()
{

    std::vector<std::vector<double>> X;

    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ )
    {
        std::vector<double> X_i;
        for (unsigned int distance_i = 0 ; distance_i < NusedSummaries+1; distance_i++)
        {
            if (distance_i == 0)
            {
                X_i.push_back(1.0);
            }
            else
            {
                X_i.push_back(std::get<distancesGetter>(population_t[simu_i])[distance_i]);
            }

        }
        X.push_back(X_i);
    }

    return X;

}

std::vector<std::vector<double>> Posterior::GetPartialDistancesT()
{

    std::vector<std::vector<double>> X_T;


    for (int distance_i = 0 ; distance_i < NusedSummaries+1; distance_i++)
    {
        std::vector<double> X_Ti;
        for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ )
        {
            if (distance_i == 0)
            {
                X_Ti[simu_i] = 1.0;
            }
            else
            {
                X_Ti.push_back(std::get<distancesGetter>(population_t[simu_i])[distance_i]);
            }


        }
        X_T.push_back(X_Ti);
    }

    return X_T;

}



std::vector<std::vector<double>> Posterior::GetTheta()
{

    std::vector<std::vector<double> > theta;

    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ )
    {

        theta.push_back(std::get<paramGetter>(population_t[simu_i]));



    }

    return theta;

}




double Posterior::GetEpanechnikov(double x, double y)
{

    return 1-(x*x)/(y*y);

}

std::vector<double> Posterior::GetWeights()
{
    std::vector<double> W;

    double max_ = std::get<distancesGetter>(population_t[population_t.size()-1])[NusedSummaries];


    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ )
    {
        W.push_back(GetEpanechnikov(std::get<distancesGetter>(population_t[simu_i])[NusedSummaries],max_));

    }

    return W;
}


std::vector<std::vector<double>> Posterior::GetLocalWeights()
{
    std::vector<std::vector<double>> W;

    double* max_ = new double[NusedSummaries];
    for (unsigned int distance_i = 0 ; distance_i < NusedSummaries; distance_i++)
    {
        max_[distance_i] = std::get<distancesGetter>(population_t[population_t.size()-1])[distance_i];
    }

    for (unsigned int simu_i = 0 ; simu_i < population_t.size(); simu_i++ )
    {
        std::vector<double> W_i;
        for (unsigned int distance_i = 0 ; distance_i < NusedSummaries; distance_i++)
        {
            W_i.push_back(GetEpanechnikov(std::get<distancesGetter>(population_t[simu_i])[distance_i],max_[distance_i]));
        }

        W.push_back(W_i);
    }

    return W;
}


void Posterior::writePosterior(ofstream&os)
{

    for (unsigned int simu_i = 0 ; simu_i < this->threshold; simu_i++ )
    {

        //write chainID
        os << std::get<chainIDGetter>(population_t[simu_i]) << "\t";

        //write parameters
        if (this->NusedParam > 0)
        {
            for (unsigned int param_i = 0 ; param_i < std::get<paramGetter>(population_t[simu_i]).size(); param_i++)
            {
                os << std::get<paramGetter>(population_t[simu_i])[param_i] << "\t";
            }
        }
        //write summaries
        if (this->NusedSummaries > 0)
        {
            for (unsigned int summary_i = 0 ; summary_i < std::get<summariesGetter>(population_t[simu_i]).size(); summary_i++)
            {
                os << std::get<summariesGetter>(population_t[simu_i])[summary_i] << "\t";
            }
        }

        //write distances (sum of square discrepancies)
        if (this->NusedSummaries > 0)
        {
            if (this->OutPartialDistance)
            {
                for (unsigned int distance_i = 0 ; distance_i < std::get<distancesGetter>(population_t[simu_i]).size(); distance_i++)
                {
                    os << std::get<distancesGetter>(population_t[simu_i])[distance_i] << "\t";
                }
            }
            else
            {

                os << std::get<distancesGetter>(population_t[simu_i])[std::get<distancesGetter>(population_t[simu_i]).size()-1] << "\t";

            }
        }

        //write accessory summaries
        if (this->NusedAccessorySummaries > 0)
        {
            for (unsigned int summary_i = 0 ; summary_i < std::get<accsummariesGetter>(population_t[simu_i]).size(); summary_i++)
            {
                os << std::get<accsummariesGetter>(population_t[simu_i])[summary_i] << "\t";
            }
        }

        //write mappingstats
        if (this->NusedEvoAncStats > 0)
        {
            for (unsigned int mapping_i = 0 ; mapping_i < std::get<evoancstatsGetter>(population_t[simu_i]).size(); mapping_i++)
            {
                if(mapping_i < std::get<evoancstatsGetter>(population_t[simu_i]).size()-1)
                {
                    os << std::get<evoancstatsGetter>(population_t[simu_i])[mapping_i] << "\t";
                }
                else
                {
                    os << std::get<evoancstatsGetter>(population_t[simu_i])[mapping_i] << "\t";
                }
            }
        }

        //write mappingstats
        if (this->NusedEvoStats > 0)
        {
            for (unsigned int mapping_i = 0 ; mapping_i < std::get<evostatsGetter>(population_t[simu_i]).size(); mapping_i++)
            {
                if(mapping_i < std::get<evostatsGetter>(population_t[simu_i]).size()-1)
                {
                    os << std::get<evostatsGetter>(population_t[simu_i])[mapping_i] << "\t";
                }
                else
                {
                    os << std::get<evostatsGetter>(population_t[simu_i])[mapping_i] << "\t";
                }
            }
        }

        //write mappingstats
        if (this->NusedSiteSpecificEvoStats > 0)
        {
            for (unsigned int mapping_i = 0 ; mapping_i < std::get<ssevostatsGetter>(population_t[simu_i]).size(); mapping_i++)
            {
                if(mapping_i < std::get<ssevostatsGetter>(population_t[simu_i]).size()-1)
                {
                    os << std::get<ssevostatsGetter>(population_t[simu_i])[mapping_i] << "\t";
                }
                else
                {
                    os << std::get<ssevostatsGetter>(population_t[simu_i])[mapping_i] << "\t";
                }
            }
        }


        os <<  "\n";
    }
}



void Posterior::readPosterior(string posteriorfile)
{
    ifstream is(posteriorfile.c_str());
    if (!is)
    {
        cerr << "error: did not find " << posteriorfile << "\n";
        exit(1);
    }
    string line;
    std::getline(is, line);
    istringstream iss(line);
    string w;
    int k = 0;
    while(iss >> w)
    {

        auto it = mapUsedParam.find(w);
        if (it != mapUsedParam.end())
        {
            cerr << it->second << " " << k << "\n";
            it->second = k;
            k++;

        }
        else
        {

            cerr << "Undefined parameter " << w << "\n";
            exit(0);

        }

    }

    NusedParam = k;
    while(std::getline(is, line))
    {
        if(!line.empty())
        {
            //cerr << line << "\n";
            istringstream iss(line);
            double w;
            std::vector <double> cur_param;

            for (unsigned int param_i = 0; param_i < NusedParam; param_i++)
            {
                iss >> w;
                cur_param.push_back(w);
            }

            posterior.push_back(cur_param);
        }
    }
    is.close();
}

void Posterior::readPosterior(ifstream& is)
{
    is.clear();                 // clear fail and eof bits
    is.seekg(0, std::ios::beg);
    int verbose = 0;



    //check correspondence between control file and posterior file



    if(verbose)
    {
        cerr << "readPosterior NusedParam"<< this->NusedParam << "\n";
    }

    std::map<int,string> mapHeader;
    int mapHeaderIndex = 0;
    if (this->NusedParam > 0)
    {
        is.clear();                 // clear fail and eof bits
        is.seekg(0, std::ios::beg); // back to the start!
        string line;
        std::getline(is, line);
        istringstream iss(line);
        string w;
        int k = 0;
        string* arr = new string[this->NusedParam];
        while(iss >> w || k < this->NusedParam)
        {

            auto it = mapUsedParam.find(w);
            if (it != mapUsedParam.end())
            {
                if(it->second != -1)
                {
                    if(k!=it->second)
                    {
                        cerr << k << " "<< it->second << " " << w << " " << it->first << "\n";
                        exit(0);
                    }

                    if (w == "chainID")
                    {
                        arr[it->second] = it->first;
                        mapHeader[mapHeaderIndex] = "chainID";
                        k++;
                        mapHeaderIndex++;

                    }
                    else
                    {
                        arr[it->second] = it->first;
                        mapHeader[mapHeaderIndex] = "P";
                        k++;
                        mapHeaderIndex++;
                    }

                }
            }
        }
        for(int v = 0 ; v < this->NusedParam; v++)
        {
            cerr << arr[v] << "\n";
        }
        delete[] arr;
    }

    if(verbose)
    {
        cerr << "readPosterior NusedSummaries"<< this->NusedSummaries << "\n";
    }

    if (this->NusedSummaries > 0)
    {
        is.clear();                 // clear fail and eof bits
        is.seekg(0, std::ios::beg); // back to the start!
        string line;
        std::getline(is, line);
        istringstream iss(line);
        string w;
        int k = 0;
        string* arr = new string[this->NusedSummaries+1];
        while(iss >> w || k < this->NusedSummaries+1)
        {

            auto it = mapUsedSummaries.find(w);

            if (it != mapUsedSummaries.end())
            {
                if(it->second != -1)
                {
                    if(k!=it->second)
                    {
                        cerr << k << " "<< it->second << " " << w << " " << it->first << "\n";
                        exit(0);
                    }
                    arr[it->second] = it->first;
                    mapHeader[mapHeaderIndex] = "S";
                    k++;
                    mapHeaderIndex++;
                }
            }
            else if (w == "D_sum")
            {
                if (k != this->NusedSummaries)
                {
                    cerr << k << " "<< this->NusedSummaries << " " << w << " " << "\n";
                    exit(0);
                }
                arr[this->NusedSummaries] = "D_sum";
                mapHeader[mapHeaderIndex] = "D_sum";
                k++;
                mapHeaderIndex++;
            }


        }

        for(int v = 0 ; v < this->NusedSummaries+1; v++)
        {
            cerr << arr[v] << "\n";
        }
        delete[] arr;
    }


//    if (this->OutPartialDistance)
//    {
//        is.clear();                 // clear fail and eof bits
//        is.seekg(0, std::ios::beg); // back to the start!
//        string line;
//        std::getline(is, line);
//        istringstream iss(line);
//        string w;
//        int k = 0;
//        string* arr = new string[this->NusedSummaries];
//        while(iss >> w || k < this->NusedSummaries)
//        {
//
//
//            w = w.erase(0,2);
//            auto it = mapUsedSummaries.find(w);
//
//            if (it != mapUsedSummaries.end())
//            {
//                if(it->second != -1)
//                {
//                    if(k!=it->second)
//                    {
//                        cerr << k << " "<< it->second << " " << w << " " << it->first << "\n";
//                        exit(0);
//                    }
//                    arr[it->second] = it->first;
//                    mapHeader[mapHeaderIndex] = "D";
//                    k++;
//                    mapHeaderIndex++;
//                }
//            }
//
//        }
//        for(int v = 0 ; v < this->NusedSummaries; v++)
//        {
//            cerr << "D_"<< arr[v] << "\n";
//        }
//        delete[] arr;
//    }


    if(verbose)
    {
        cerr << "readPosterior NusedEvoStats"<< this->NusedEvoStats << "\n";
    }

    if (this->NusedEvoStats > 0)
    {
        is.clear();                 // clear fail and eof bits
        is.seekg(0, std::ios::beg); // back to the start!
        string line;
        std::getline(is, line);
        istringstream iss(line);
        string w;
        int k = 0;
        string* arr = new string[this->NusedEvoStats];
        while(iss >> w || k < this->NusedEvoStats)
        {
            w = w.erase(0,2);
            auto it = mapUsedEvoStats.find(w);
            if (it != mapUsedEvoStats.end())
            {
                if(it->second != -1)
                {
                    if(k!=it->second)
                    {
                        cerr << k << " "<< it->second << " " << w << " " << it->first << "\n";
                        exit(0);
                    }
                    arr[it->second] = it->first;
                    mapHeader[mapHeaderIndex] = "ES";
                    k++;
                    mapHeaderIndex++;
                }
            }
        }
        for(int v = 0 ; v < this->NusedEvoStats; v++)
        {
            cerr << arr[v] << "\n";
        }
        delete[] arr;
    }


    // get all the information



    is.clear();
    is.seekg(0, std::ios::beg);
    string line;
    std::getline(is, line);// skip header
    while(std::getline(is, line))
    {
        cerr << ".";
        if(!line.empty())
        {
            //cerr << line << "\n";


            istringstream iss_tmp(line);
            int chainID = 1;
            std::vector<double> cur_param;
            std::vector<double> cur_summaries;
            std::vector<double> cur_evostats;
            std::vector<double> cur_distances;
            mapHeaderIndex = 0;
            string  w;


            while(iss_tmp >> w)
            {
                auto it = mapHeader.find(mapHeaderIndex);
                if (it != mapHeader.end())
                {

                    if (it->second == "P")
                    {
                        cur_param.push_back(std::stof(w));

                    }
                    else if (it->second == "chainID")

                    {
                        chainID = std::stoi(w);
                    }
                    else if (it->second == "S")
                    {
                        cur_summaries.push_back(std::stof(w));

                    }
                    else if (it->second == "ES")
                    {
                        cur_evostats.push_back(std::stof(w));

                    }
                    else if (it->second == "D_sum")
                    {
                        cur_distances.push_back(std::stof(w));

                    }
                    mapHeaderIndex++;

                }
            }

            std::vector<double> tmp;
            registerNewSimulation(
                chainID,
                cur_param,
                cur_summaries,
                tmp,
                tmp,
                cur_evostats,
                tmp,
                cur_distances,
                tmp
            );

        }


    }


}

void Posterior::writePosteriorPredictivePvalues(ofstream& os, std::vector<double>summariesRealData)
{


    unsigned int pop_size = population_t.size();
    int k =0 ;
    if (NusedSummaries > 0)
    {
        string* arrSummaries = new string[this->NusedSummaries];
        for(unsigned int summary_i = 0 ; summary_i < this->NSummaries; summary_i++)
        {
            auto it = mapUsedSummaries.find(this->listSummaries[summary_i]);
            if(it != mapUsedSummaries.end() )
            {
                if(it->second != -1)
                {
                    arrSummaries[it->second] = it->first;

                }

            }

        }


        for(unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
        {
            if(k == 0 )
            {
                os << arrSummaries[summary_i];
                k = 1;
            }
            else
            {
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

    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        ppp[summary_i] = 0.0;
        mean[summary_i] = 0.0;
        var[summary_i] = 0.0;
        sum1[summary_i] = 0.0;
        sum2[summary_i] = 0.0;
        sum3[summary_i] = 0.0;
    }

    for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++ )
    {
        for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
        {
            if(std::get<summariesGetter>(population_t[simu_i])[summary_i] < summariesRealData[summary_i])
            {
                ppp[summary_i]+=1;
            }
            sum1[summary_i]+=exp2(std::get<summariesGetter>(population_t[simu_i])[summary_i]);
        }
    }


    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        mean[summary_i] = sum1[summary_i]/pop_size;

    }

    for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++ )
    {
        for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
        {
            sum2[summary_i] += (std::get<summariesGetter>(population_t[simu_i])[summary_i]-mean[summary_i])*(std::get<summariesGetter>(population_t[simu_i])[summary_i]-mean[summary_i]);
            sum3[summary_i] += (std::get<summariesGetter>(population_t[simu_i])[summary_i]-mean[summary_i]);
        }
    }

    //posterior predictive p-values (ppp)
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {

        var[summary_i] = (sum2[summary_i]-(sum3[summary_i]*sum3[summary_i])/pop_size)/(pop_size-1);
        ppp[summary_i] /= pop_size;

        if(summary_i < this->NusedSummaries-1)
        {
            os << ppp[summary_i] << "\t";
        }
        else
        {
            os << ppp[summary_i] << "\n";
        }

    }

    //means
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        if(summary_i < this->NusedSummaries-1)
        {
            os << mean[summary_i] << "\t";
        }
        else
        {
            os << mean[summary_i] << "\n";
        }
    }

    //variances
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        if(summary_i < this->NusedSummaries-1)
        {
            os << var[summary_i] << "\t";
        }
        else
        {
            os << var[summary_i] << "\n";
        }
    }

    //z-scores
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        if(summary_i < this->NusedSummaries-1)
        {
            os <<  (summariesRealData[summary_i])-mean[summary_i] / sqrt(var[summary_i]) << "\t";
        }
        else
        {
            os <<  (summariesRealData[summary_i])-mean[summary_i] / sqrt(var[summary_i]) << "\n";
        }
    }

    delete [] ppp;
    delete [] mean;
    delete [] var;
    delete [] sum1;
    delete [] sum2;
    delete [] sum3;


}

void Posterior::writePosteriorPredictiveStatistics(ofstream& os, std::vector<double>summariesRealData)
{


    unsigned int pop_size = population_t.size();
    int k =0 ;
    if (NusedSummaries > 0)
    {
        string* arrSummaries = new string[this->NusedSummaries];
        for(unsigned int summary_i = 0 ; summary_i < this->NSummaries; summary_i++)
        {
            auto it = mapUsedSummaries.find(this->listSummaries[summary_i]);
            if(it != mapUsedSummaries.end() )
            {
                if(it->second != -1)
                {
                    arrSummaries[it->second] = it->first;

                }

            }

        }


        for(unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
        {
            if(k == 0 )
            {
                os << arrSummaries[summary_i];
                k = 1;
            }
            else
            {
                os << "\t" << arrSummaries[summary_i];
            }

        }
        os << "\n";
    }



    int n = 0;



    double* ppp  = new double[this->NusedSummaries];
    double* mean = new double[this->NusedSummaries];
    double* var  = new double[this->NusedSummaries];
    double* a  = new double[this->NusedSummaries];
    double* b  = new double[this->NusedSummaries];
    double* K  = new double[this->NusedSummaries];

    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        ppp[summary_i] = 0.0;
        mean[summary_i] = 0.0;
        var[summary_i] = 0.0;
        a[summary_i] = 0.0;
        b[summary_i] = 0.0;
        K[summary_i] = std::get<summariesGetter>(population_t[0])[summary_i]; //const for shifted data
    }

    for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++ )
    {
        for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
        {


            if(std::get<summariesGetter>(population_t[simu_i])[summary_i] < summariesRealData[summary_i])
            {
                ppp[summary_i]+=1;
            }

            mean[summary_i]+=std::get<summariesGetter>(population_t[simu_i])[summary_i];
            a[summary_i] +=  std::get<summariesGetter>(population_t[simu_i])[summary_i] - K[summary_i];
            b[summary_i] +=  (std::get<summariesGetter>(population_t[simu_i])[summary_i] - K[summary_i]) * (std::get<summariesGetter>(population_t[simu_i])[summary_i] - K[summary_i]);

        }
    }


    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        mean[summary_i] /= pop_size;
        var[summary_i]  =  (b[summary_i]  - (a[summary_i]*a[summary_i])/pop_size)/(pop_size);

    }



    //posterior predictive p-values (ppp)
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {

        ppp[summary_i] /= pop_size;

        if(summary_i < this->NusedSummaries-1)
        {
            os << ppp[summary_i] << "\t";
        }
        else
        {
            os << ppp[summary_i] << "\n";
        }

    }

    //means
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        if(summary_i < this->NusedSummaries-1)
        {
            os << mean[summary_i] << "\t";
        }
        else
        {
            os << mean[summary_i] << "\n";
        }
    }

    //variances
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        if(summary_i < this->NusedSummaries-1)
        {
            os << var[summary_i] << "\t";
        }
        else
        {
            os << var[summary_i] << "\n";
        }
    }

    //z-scores
    for (unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
    {
        if(summary_i < this->NusedSummaries-1)
        {
            os <<  (summariesRealData[summary_i]-mean[summary_i]) / sqrt(var[summary_i]) << "\t";
        }
        else
        {
            os <<  (summariesRealData[summary_i]-mean[summary_i]) / sqrt(var[summary_i]) << "\n";
        }
    }

    delete [] ppp;
    delete [] mean;
    delete [] var;
    delete [] a;
    delete [] b;
    delete [] K;


}
void Posterior::readMonitor(ifstream & is)
{



    is.clear();                 // clear fail and eof bits
    is.seekg(0, std::ios::beg);



    string line;
    std::getline(is, line);
    if(!line.empty())
    {
        istringstream iss1(line);
        iss1 >> this->Niter;
        std::getline(is, line);
        istringstream iss2(line);
        iss2 >> this->Naccepted;
        cerr << "Niter "  << this->Niter << " " << "Naccepted " << this->Naccepted << "\n";

    }
    else
    {
        this->Niter = 0;
        cerr << "Monitor file is empty " << "\n";
    }


}

void Posterior::writeMonitorPosterior(ofstream& os)
{
    os << this->Niter << "\n";
    os << this->Naccepted << "\n";
    os << GetAcceptanceRate() << "\n";

}

double Posterior::GetAcceptanceRate()
{
    return (double) Naccepted/Niter;
}




void Posterior::writeHeader(ofstream&os)
{

    // write parameters' header
    //for (unsigned int param_i = 0 ; param_i < listUsedParam.size() ; param_i++){
    int k = 0;

    //MCMCpt should go there




    if(verbose)
    {
        cerr << "writeHeader1 "<< this->NusedParam << "\n";
    }
    if (this->NusedParam > 0)
    {
        string* arrParam = new string[this->NusedParam];
        for (unsigned int param_i = 0 ; param_i < this->NParam ; param_i++)
        {
            auto it = this->mapUsedParam.find(this->listParam[param_i]);
            if(it != this->mapUsedParam.end() )
            {
                if(it->second != -1)
                {
                    arrParam[it->second] = it->first;

                }

            }

        }

        for (unsigned int param_i = 0 ; param_i < this->NusedParam ; param_i++)
        {
            if(k == 0)
            {
                os << arrParam[param_i];
                k = 1;
            }
            else
            {
                os << "\t" << arrParam[param_i];
            }
        }
        delete [] arrParam;
    }

    if(verbose)
    {
        cerr << "writeHeader2 "<< this->NusedSummaries << "\n";
    }
    //write summaires' header
    if (this->NusedSummaries > 0)
    {
        string* arrSummaries = new string[this->NusedSummaries];
        for(unsigned int summary_i = 0 ; summary_i < this->NSummaries; summary_i++)
        {
            auto it = this->mapUsedSummaries.find(this->listSummaries[summary_i]);
            if(it != this->mapUsedSummaries.end() )
            {
                if(it->second != -1)
                {
                    arrSummaries[it->second] = it->first;

                }

            }

        }


        for(unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
        {
            if(k == 0 )
            {
                os << arrSummaries[summary_i];
                k = 1;
            }
            else
            {
                os << "\t" << arrSummaries[summary_i];
            }

        }

        //write distance header
        if(this->OutPartialDistance)
        {
            for(unsigned int summary_i = 0 ; summary_i < this->NusedSummaries; summary_i++)
            {
                if (k == 0)
                {
                    os << "D_" << arrSummaries[summary_i];
                    k = 1;
                }
                else
                {
                    os << "\t" << "D_" << arrSummaries[summary_i];
                }
            }
        }
        os << "\tD_sum";
        delete [] arrSummaries;


    }


    if (this->NusedAccessorySummaries > 0)
    {
        string* arrSummaries = new string[this->NusedAccessorySummaries];
        for(unsigned int summary_i = 0 ; summary_i < this->NSummaries; summary_i++)
        {
            auto it = this->mapUsedAccessorySummaries.find(this->listSummaries[summary_i]);
            if(it != this->mapUsedAccessorySummaries.end() )
            {
                if(it->second != -1)
                {
                    arrSummaries[it->second] = it->first;

                }

            }

        }


        for(unsigned int summary_i = 0 ; summary_i < this->NusedAccessorySummaries; summary_i++)
        {
            if(k == 0 )
            {
                os << "Acc_" << arrSummaries[summary_i];
                k = 1;
            }
            else
            {
                os << "\t" << "Acc_" <<arrSummaries[summary_i];
            }

        }

        delete [] arrSummaries;


    }


    //write mappingstats
    if(verbose)
    {
        cerr << "writeHeader3 "<< this->NusedEvoAncStats << "\n";
    }
    if (this->NusedEvoAncStats > 0)
    {
        string* arrMapAnc = new string[this->NusedEvoAncStats];
        for(unsigned int map_i = 0 ; map_i < this->NEvoStats; map_i++)
        {
            auto it = this->mapUsedEvoAncStats.find(this->listEvoStats[map_i]);
            if(it != this->mapUsedEvoAncStats.end())
            {
                if(it->second != -1)
                {
                    arrMapAnc[it->second] = it->first;

                }

            }
        }


        for(unsigned int map_i = 0; map_i < NusedEvoAncStats; map_i++)
        {
            if (k == 0)
            {
                os << "A_" << arrMapAnc[map_i];
                k = 1;
            }
            else
            {
                os << "\t" << "A_" << arrMapAnc[map_i];
            }


        }
        delete [] arrMapAnc;
    }


    if(verbose)
    {
        cerr << "writeHeader4 "<< this->NusedEvoStats << "\n";
    }
    if (this->NusedEvoStats > 0)
    {
        string* arrMap = new string[this->NusedEvoStats];
        for(unsigned int map_i = 0 ; map_i < this->NEvoStats; map_i++)
        {
            auto it = this->mapUsedEvoStats.find(this->listEvoStats[map_i]);
            if(it != this->mapUsedEvoStats.end())
            {
                if(it->second != -1)
                {
                    arrMap[it->second] = it->first;

                }
            }
        }


        for(unsigned int map_i = 0; map_i < this->NusedEvoStats; map_i++)
        {
            if (k == 0)
            {
                os << "T_" << arrMap[map_i];
                k = 1;
            }
            else
            {
                os << "\t" << "T_" << arrMap[map_i];
            }

        }
        delete [] arrMap;


    }

    if(verbose)
    {
        cerr << "writeHeader5 "<< this->NusedSiteSpecificEvoStats << "\n";
    }
    if (this->NusedSiteSpecificEvoStats > 0)
    {
        string* arrMap = new string [this->NusedSiteSpecificEvoStats];
        for(unsigned int map_i = 0 ; map_i < this->NSiteSpecificEvoStats; map_i++)
        {
            auto it = this->mapUsedSiteSpecificEvoStats.find(this->listSiteSpecificEvoStats[map_i]);
            if(it != this->mapUsedSiteSpecificEvoStats.end())
            {
                if(it->second != -1)
                {
                    arrMap[it->second] = it->first;

                }
            }
        }


        for(unsigned int map_i = 0; map_i < this->NusedSiteSpecificEvoStats; map_i++)
        {
            for (int bin_i = 0; bin_i < 100; bin_i++)
            {

                if (k == 0)
                {
                    os << "SS_" << bin_i << "_" << arrMap[map_i];
                    k = 1;
                }
                else
                {
                    os << "\t" << "SS_" << bin_i << "_" <<  arrMap[map_i];
                }
            }
        }
        delete [] arrMap;


    }

    os << "\n";

}

void Posterior::SetNsite(int i)
{
    this->Nsite_codon = i;


}

void Posterior::GetWeights(string kernel)
{
    unsigned int pop_size = population_t.size();
    if (kernel == "sNormal")
    {

        for (unsigned int param_i = 0 ; param_i < NusedParam ; param_i++)
        {
            double new_weight = 0.0;
            for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++)
            {
                new_weight += std::get<weightsGetter>(population_t[simu_i])[param_i] * (std::get<paramGetter>(population_t[simu_i])[param_i]  + 2.0 * (empVar[param_i]) * rnd->sNormal());

            }
            new_weight /= pop_size;
//            std::get<4>(population_t[simu_i])[param_i] = new_weight;
        }
    }
}

void Posterior::GetEmpVar()
{

    unsigned int pop_size = population_t.size();


    double* mean = new double[NusedParam];
    double* sum1  = new double[NusedParam];
    double* sum2  = new double[NusedParam];
    double* sum3  = new double[NusedParam];

    for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++)
    {
        mean[param_i] = 0.0;
        sum1[param_i] = 0.0;
        sum2[param_i] = 0.0;
        sum3[param_i] = 0.0;
    }

    for (unsigned int simu_i = 0 ; simu_i < NusedParam; simu_i++ )
    {
        for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++)
        {
            sum1[param_i] += std::get<paramGetter>(population_t[simu_i])[param_i];
        }
    }


    for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++)
    {
        mean[param_i] = sum1[param_i]/pop_size;
        empMean[param_i] = mean[param_i];
    }

    for (unsigned int simu_i = 0 ; simu_i < pop_size; simu_i++ )
    {
        for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++)
        {
            //(x - K) * (x - K)
            sum2[param_i] += (std::get<paramGetter>(population_t[simu_i])[param_i]-mean[param_i])*(std::get<paramGetter>(population_t[simu_i])[param_i]-mean[param_i]);
            //x - K
            sum3[param_i] += std::get<paramGetter>(population_t[simu_i])[param_i]-mean[param_i];
        }
    }


    for (unsigned int param_i = 0 ; param_i < NusedParam; param_i++)
    {
        // sum_sqr - (sum_ * sum_)/n)/(n - 1)
        empVar[param_i] = (sum2[param_i]-(sum3[param_i]*sum3[param_i])/pop_size)/(pop_size-1);
    }

    delete [] mean;
    delete [] sum1;
    delete [] sum2;
    delete [] sum3;

}


void Posterior::slaveToMaster(std::vector<std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>> population_i){

    population_t.insert(population_t.end(),population_i.begin(),population_i.end());

}



void Posterior::sortPopulation(){
    if (!sorted)
    {

        std::sort(population_t.begin(), population_t.end(),
                  [](const std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &left,
                     const std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &right)
        {
            return std::get<distancesGetter>(left).back() < std::get<distancesGetter>(right).back();
        }
                 );
        sorted = true;
    }


}

void Posterior::slaveRegisterNewSimulation(int chainID, std::vector<double> param,std::vector<double> summaries,std::vector<double> accsummaries,std::vector<double> evoancstat,std::vector<double> evostat,std::vector<double> ssevostat,std::vector<double> distances, std::vector<double> weights)
{

    population_t.push_back(make_tuple(chainID,param,summaries,accsummaries,evoancstat, evostat,ssevostat,distances,weights));
    Naccepted++;
    this->Niter++;

}



void Posterior::registerNewSimulation(int chainID, std::vector<double> param,std::vector<double> summaries,std::vector<double> accsummaries,std::vector<double> evoancstat,std::vector<double> evostat,std::vector<double> ssevostat,std::vector<double> distances, std::vector<double> weights)
{

    if(population_t.empty())
    {

        //cerr << "POPULATION IS EMPTY" << Naccepted <<  " "<< threshold <<  " \n";

        population_t.push_back(make_tuple(chainID,param,summaries,accsummaries,evoancstat, evostat,ssevostat,distances,weights));
        Naccepted++;

    }
    else if(population_t.size() < (unsigned) threshold)
    {

        //cerr << "POPULATION IS LESS THAN THRESHOLD" << " " << population_t.size() << " "<< Naccepted  << " "<< threshold << " \n";

        population_t.push_back(make_tuple(chainID,param,summaries,accsummaries,evoancstat, evostat,ssevostat,distances,weights));
        Naccepted++;

    }
    else
    {

        if (!sorted)
        {

            std::sort(population_t.begin(), population_t.end(),
                      [](const std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &left,
                         const std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &right)
            {
                return std::get<distancesGetter>(left).back() < std::get<distancesGetter>(right).back();
            }
                     );
            sorted = true;
        }

        std::tuple<int, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> cur_tuple
                = make_tuple(chainID,param,summaries,accsummaries,evoancstat, evostat,ssevostat,distances,weights);

        auto it = std::lower_bound(population_t.begin(), population_t.end(), cur_tuple,
                                   [](const std::tuple<int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &left,
                                      const std::tuple<int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> &right)
        {
            return std::get<distancesGetter>(left).back() < std::get<distancesGetter>(right).back();
        }
                                  );

        // compute acceptance rate
        if (it != population_t.end())
        {
            //cerr << "INSERTED" << " " << population_t.size()  << " "<<  Naccepted << " "<< threshold << " \n";
            population_t.insert(it,cur_tuple);
            population_t.pop_back();
            population_t.shrink_to_fit();
            Naccepted++;
        }

    }
    this->Niter++;
    //cerr << this->Niter << " ";
}
