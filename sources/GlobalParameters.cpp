#include "GlobalParameters.h"

GlobalParameters::GlobalParameters(string model, string controlfile)
{
    this->model = model;
    this->controlfile = controlfile;
    this->Nrun = 1;
    this->Nthread = 1;
    this->Niter = 0;
    this->threshold = 10000;
    this->verbose = 0;

    this->chainPointStart = 1;
    this->chainPointEnd = 101;
    this->chainPointEvery = 1;

    this->NusedEvoStats = 0;
    this->NusedSiteSpecificEvoStats = 0;
    this->NusedEvoAncStats = 0;
    this->NusedParam = 0;
    this->NusedSummaries = 0;
    this->NusedAccessorySummaries = 0;
    this->Ngenes = 0;

    this->distance = "Euclidian";
    this->transformation = "log2";
    this->OutPartialDistance = 0;

    cerr << "Constructing mapUsedParam\n";
    for (unsigned int i = 0 ; i < this->NParam; i++ )
    {
        mapUsedParam[listParam[i]] = -1;
    }

    cerr << "Constructing mapUsedSummaries\n";
    for (unsigned int i = 0 ; i <  this->NSummaries; i++ )
    {
        mapUsedSummaries[listSummaries[i]] = -1;
    }

    cerr << "Constructing mapUsedAccessorySummaries\n";
    for (unsigned int i = 0 ; i <  this->NSummaries; i++ )
    {
        mapUsedAccessorySummaries[listSummaries[i]] = -1;
    }

    cerr << "Constructing mapUsedEvoStats\n";
    for (unsigned int i = 0 ; i <  this->NEvoStats; i++ )
    {
        mapUsedEvoStats[listEvoStats[i]] = -1;
    }

    cerr << "Constructing mapUsedSiteSpecificEvoStats\n";
    for (unsigned int i = 0 ; i <  this->NSiteSpecificEvoStats; i++ )
    {
        mapUsedSiteSpecificEvoStats[listSiteSpecificEvoStats[i]] = -1;
    }

    cerr << "Constructing mapUsedEvoAncStats\n";
    for (unsigned int i = 0 ; i <  this->NEvoStats; i++ )
    {
        mapUsedEvoAncStats[listEvoStats[i]] = -1;
    }

    cerr << "Reading instructions\n";
    readInstructions();


}

GlobalParameters::~GlobalParameters()
{
    //dtor
}





void GlobalParameters::readInstructions()
{

    ifstream is(this->controlfile.c_str());
    if (!is)
    {
        cerr << "error: did not find " << this->controlfile << "\n";
        exit(1);
    }

    string line = "";
    while(std::getline(is, line))
    {
        //cerr << line << "\n";
        if(!line.empty() && line.substr(0,10) == "#SUMMARIES")
        {
            istringstream iss(line);
            string w;
            int k = 0 ;
            while(iss >> w)

            {

                auto it = mapUsedSummaries.find(w);
                if (it != mapUsedSummaries.end())
                {

                    it->second = k;
                    k++;

                }
                else if (w != "#SUMMARIES")
                {

                    cerr << "Undefined summary " << w << "\n";
                    exit(0);

                }

            }
            NusedSummaries = k;
            cerr << "#SUMMARIES " << NusedSummaries << "\n";

        }
        else if(!line.empty() && line.substr(0,13) == "#ACCSUMMARIES")
        {
            istringstream iss(line);
            string w;
            int k = 0 ;
            while(iss >> w)

            {

                auto it = mapUsedAccessorySummaries.find(w);
                if (it != mapUsedAccessorySummaries.end())
                {

                    it->second = k;
                    k++;

                }
                else if (w != "#ACCSUMMARIES")
                {

                    cerr << "Undefined summary " << w << "\n";
                    exit(0);

                }

            }
            NusedAccessorySummaries = k;
            cerr << "#ACCSUMMARIES " << NusedAccessorySummaries << "\n";

        }
        else if(!line.empty() && line.substr(0,6) == "#PARAM")
        {
            istringstream iss(line);
            string w;
            int k = 0;
            while(iss >> w)

            {

                auto it = mapUsedParam.find(w);
                if (it != mapUsedParam.end())
                {

                    it->second = k;
                    cerr << "UsedParam " << k << " " <<  w << "\n";
                    k++;

                }
                else if (w !="#PARAM")
                {


                    cerr << "Undefined parameter " << w << "\n";
                    exit(0);

                }

            }
            NusedParam = k;
            cerr << "#PARAM " << NusedParam << "\n";

        }
        else if(!line.empty() && line.substr(0,6) == "#SSMAP")
        {
            istringstream iss(line);
            string w;
            int k = 0;
            while(iss >> w)

            {

                auto it = mapUsedSiteSpecificEvoStats.find(w);
                if (it != mapUsedSiteSpecificEvoStats.end())
                {

                    it->second = k;
                    k++;

                }
                else if (w !="#SSMAP")
                {

                    cerr << "Undefined SiteSpecificEvoStats " << w << "\n";
                    exit(0);

                }

            }
            NusedSiteSpecificEvoStats = k;
            cerr << "#SSMAP " << NusedSiteSpecificEvoStats << "\n";

        }
        else if(!line.empty() && line.substr(0,4) == "#MAP")
        {
            istringstream iss(line);
            string w;
            int k = 0;
            while(iss >> w)

            {

                auto it = mapUsedEvoStats.find(w);
                if (it != mapUsedEvoStats.end())
                {

                    it->second = k;
                    k++;

                }
                else if (w !="#MAP")
                {

                    cerr << "Undefined EvoStats parameter " << w << "\n";
                    exit(0);

                }

            }
            NusedEvoStats = k;
            cerr << "#MAP " << NusedEvoStats << "\n";

        }
        else if(!line.empty() && line.substr(0,13) == "#ANCESTRALMAP")
        {
            istringstream iss(line);
            string w;
            int k = 0;
            while(iss >> w)

            {

                auto it = mapUsedEvoAncStats.find(w);
                if (it != mapUsedEvoAncStats.end())
                {

                    it->second = k;
                    k++;

                }
                else if (w !="#ANCESTRALMAP")
                {

                    cerr << "Undefined EvoAncStats parameter" << w << "\n";
                    exit(0);

                }

            }
            this->NusedEvoAncStats = k;
            cerr << "#ANCESTRALMAP " << NusedEvoAncStats << "\n";

        }
        else if (!line.empty() && line.substr(0,6) == "#GENES")
        {
            istringstream iss(line);
            string w;
            int k = 0 ;
            while(iss >> w)
            {
                if (w != "#GENES")
                {
                    k++;
                    this->listGenes.push_back(w);
                }

            }
            this->Ngenes = k;
            cerr << "#GENES\t" << this->Ngenes << "\n";
            for (int i = 0 ; i < this->listGenes.size(); i++)
            {
                cerr << this->listGenes[i] << "\t";
            }
            cerr << "\n";

        }
        else if (!line.empty() && line.substr(0,5) == "#TRANS")
        {
            istringstream iss(line);
            string w;
            int k = 0 ;
            while(iss >> w)
            {
                if (w != "#TRANS")
                {
                    k++;
                    this->transformation = w;
                }

            }

            cerr << "#TRANS\t" << this->transformation << "\n";


        }
        else if (!line.empty() && line.substr(0,8) == "#OUTDIST")
        {

            this->OutPartialDistance = 1;

            cerr << "#OUTDIST\t" << this->distance << "\n";


        }

        else if (!line.empty() && line.substr(0,11) == "#SPEUDODATA")
        {
            istringstream iss(line);
            string w;
            int k = 0 ;
            iss >> w;
            iss >> w;
            this->Ntaxa = atoi(w.c_str());
            iss >> w;
            this->Nsite_codon = atoi(w.c_str());
            while(iss >> w)
            {
                if (w != "#SPEUDODATA")
                {
                    k++;
                    this->listSpecies.push_back(w);
                }

            }

            if (this->Ntaxa != this->listSpecies.size())
            {
                cerr << "Error the number of taxa does not match the list species size\n";
                exit(0);
            }

            cerr << "#SPEUDODATA " << this->Ntaxa << " " << this->Nsite_codon << "\n";
            for (auto i : this->listSpecies)
            {
                cerr << i << " ";
            }
            cerr << "\n";



        }
        else if (!line.empty() && line.substr(0,5) == "#DIST")
        {
            istringstream iss(line);
            string w;
            int k = 0 ;
            while(iss >> w)
            {
                if (w != "#DIST")
                {
                    k++;
                    this->distance = w;
                }

            }

            cerr << "#DIST\t" << this->distance << "\n";


        }
        else if (!line.empty() && line.substr(0,7) == "#CHAINS")
        {
            istringstream iss(line);
            string w;
            while(iss >> w)
            {
                listChains.push_back(w);
                iss >> w;
                listPoints.push_back(atoi(w.c_str()));

            }

        }
        else if (!line.empty() && line.substr(0,15) == "#HYPERMUT-DINUC")
        {
            cerr << "### HYPERMUT-DINUC ###\n";
            istringstream iss(line);
            string w;
            while(iss >> w)
            {
                cerr << w << "\n";
            }

        }
        else if (!line.empty() && line.substr(0,5) == "#NRUN")
        {
            istringstream iss(line);
            string w;
            iss >> w;
            iss >> w;
            this->Nrun = atoi(w.c_str());
            iss >> w;
            this->threshold = atoi(w.c_str());
            cerr << "#Nrun " << this->Nrun << " threshold " << this->threshold << "\n";

        }
        else if (!line.empty() && line.substr(0,9) == "#NTHREADS")
        {
            istringstream iss(line);
            string w;
            iss >> w;
            iss >> w;
            this->Nthread = atoi(w.c_str());
            cerr << "#Nthread " << this->Nthread << "\n";

        }
        else if (!line.empty() && line.substr(0,7) == "#OUTPUT")
        {
            istringstream iss(line);
            string w;
            iss >> w;
            iss >> w;
            this->output = w;
            cerr << "#OUTPUT " << this->output << "\n";

        }
        else if (!line.empty() && line.substr(0,8) == "#VERBOSE")
        {
            this->verbose = 1;
            cerr << "#VERBOSE " << this->verbose << "\n";

        }
        else if (!line.empty() && line.substr(0,9) == "#SAMPLING")
        {
            istringstream iss(line);
            string w;
            iss >> w;
            iss >> w;
            this->chainPointStart = atoi(w.c_str());
            iss >> w;
            this->chainPointEvery = atoi(w.c_str());
            iss >> w;
            this->chainPointEnd = atoi(w.c_str());
            cerr << "#SAMPLING " <<  this->chainPointStart << " " << this->chainPointEvery << " " << this->chainPointEnd << "\n";

        }
        else if (!line.empty() && line.substr(0,10) == "#OLDPARAMS")
        {
            cerr << "### OLDPARAMS ### \n";
            localcontrolfile = line;

        }
        else if (!line.empty() && line.substr(0,11) == "#LOCALPARAM")
        {
            cerr << "### LOCALPARAM ### \n";
            localcontrolfile = line;

        }
        /////////////////////////
    }
    is.close();

}


