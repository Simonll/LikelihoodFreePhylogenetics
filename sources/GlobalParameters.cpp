#include "GlobalParameters.h"

GlobalParameters::GlobalParameters(string model, string controlfile)
{
    this->model = model;
    this->controlfile = controlfile;
    this->Nrun = 0;
    this->Nthread = 1;
    this->Niter = 0;
    this->threshold = 10000;

    this->chainPointStart = 1;
    this->chainPointEnd = 101;
    this->chainPointEvery = 1;

    this->NusedMapStats = 0;
    this->NusedMapAncStats = 0;
    this->NusedParam = 0;
    this->NusedSummaries = 0;
    this->Ngenes = 0;

    this->distance = "log2";

    cerr << "Constructing mapUsedParam\n";
    for (unsigned int i = 0 ; i < NParam; i++ ){
       mapUsedParam[listParam[i]] = -1;
    }

    cerr << "Constructing mapUsedSummaries\n";
    for (unsigned int i = 0 ; i < NSummaries; i++ ){
       mapUsedSummaries[listSummaries[i]] = -1;
    }

    cerr << "Constructing mapUsedMapStats\n";
    for (unsigned int i = 0 ; i < NMapStats; i++ ){
       mapUsedMapStats[listMapStats[i]] = -1;
    }

    cerr << "Constructing mapUsedMapAncStats\n";
    for (unsigned int i = 0 ; i < NMapStats; i++ ){
       mapUsedMapAncStats[listMapStats[i]] = -1;
    }

    cerr << "Reading instructions\n";
    readInstructions();


}

GlobalParameters::~GlobalParameters()
{
    //dtor
}





void GlobalParameters::readInstructions() {

    ifstream is(this->controlfile.c_str());
    if (!is)       {
        cerr << "error: did not find " << this->controlfile << "\n";
        exit(1);
    }

    string line = "";
    while(std::getline(is, line)) {
        //cerr << line << "\n";
       if(!line.empty() && line.substr(0,10) == "#SUMMARIES") {

           istringstream iss(line);
           string w;
           int k = 0 ;
           while(iss >> w)

           {
//               for (int i = 0; i < NSummaries; i++) {
//                   if (w == listSummaries[i]) {
//                       cerr << w << " " << i << "\n";
//                       listUsedSummaries.push_back(i);
//                   }
//               }

                auto it = mapUsedSummaries.find(w);
                if (it != mapUsedSummaries.end()) {

                    it->second = k;
                    //cerr << "UsedSummaries " << k << "\n";
                    k++;

                } else if (w != "#SUMMARIES"){


                    cerr << "Undefined summary " << w << "\n";
                    exit(0);

                }

        }
        NusedSummaries = k;
        cerr << "#SUMMARIES " << NusedSummaries << "\n";
       } else if(!line.empty() && line.substr(0,6) == "#PARAM") {

           istringstream iss(line);
           string w;
           int k = 0;
           while(iss >> w)

           {

                auto it = mapUsedParam.find(w);
                if (it != mapUsedParam.end()) {

                    it->second = k;
                    //cerr << "UsedParam " << k << "\n";
                    k++;

                } else if (w !="#PARAM"){


                    cerr << "Undefined parameter " << w << "\n";
                    exit(0);

                }

           }
           NusedParam = k;
           cerr << "#PARAM " << NusedParam << "\n";
       } else if(!line.empty() && line.substr(0,4) == "#MAP") {

           istringstream iss(line);
           string w;
           int k = 0;
           while(iss >> w)

           {

                auto it = mapUsedMapStats.find(w);
                if (it != mapUsedMapStats.end()) {

                    it->second = k;
                    //cerr << "UsedParam " << k << "\n";
                    k++;

                } else if (w !="#MAP"){

                    cerr << "Undefined mapping parameter " << w << "\n";
                    exit(0);

                }

           }
           NusedMapStats = k;
           cerr << "#MAP " << NusedMapStats << "\n";
       }  else if(!line.empty() && line.substr(0,13) == "#ANCESTRALMAP") {

           istringstream iss(line);
           string w;
           int k = 0;
           while(iss >> w)

           {

                auto it = mapUsedMapAncStats.find(w);
                if (it != mapUsedMapAncStats.end()) {

                    it->second = k;
                    //cerr << "UsedParam " << k << "\n";
                    k++;

                } else if (w !="#ANCESTRALMAP"){

                    cerr << "Undefined parameter" << w << "\n";
                    exit(0);

                }

           }
           this->NusedMapAncStats = k;
           cerr << "#ANCESTRALMAP " << NusedMapAncStats << "\n";
       } else if (!line.empty() && line.substr(0,6) == "#GENES") {

           istringstream iss(line);
           string w;
           int k = 0 ;
           while(iss >> w)
           {
                if (w != "#GENES") {
                    k++;
                    this->listGenes.push_back(w);
                }

           }
           this->Ngenes = k;
           cerr << "#GENES\t" << this->Ngenes << "\n";
           for (int i = 0 ; i < this->listGenes.size(); i++){
                cerr << this->listGenes[i] << "\t";
           }
           cerr << "\n";

       } else if (!line.empty() && line.substr(0,5) == "#DIST") {

           istringstream iss(line);
           string w;
           int k = 0 ;
           while(iss >> w)
           {
                if (w != "#DIST") {
                    k++;
                    this->distance = w;
                }

           }

           cerr << "#DIST\t" << this->distance << "\n";


       } else if (!line.empty() && line.substr(0,9) == "#SAMPLING") {

           istringstream iss(line);
           string w;
           int k = 0 ;
           while(iss >> w)
           {
                if (w != "#SAMPLING") {
                    k++;
                    if (k == 1) {
                        this->chainPointStart = atoi(w.c_str());
                    } else if (k == 2) {
                        this->chainPointEvery = atoi(w.c_str());
                    } else if (k == 3) {
                        this->chainPointEnd = atoi(w.c_str());
                    }
                }

           }

           cerr << "#SAMPLING\t" << this->chainPointStart << "\t" << this->chainPointEvery << "\t" << this->chainPointEnd << "\n";


       } else if (!line.empty() && line.substr(0,7) == "#CHAINS") {
           cerr << "### CHAINS ###\n";
           istringstream iss(line);
           string w;
           while(iss >> w)
           {
               listChains.push_back(w);
               iss >> w;
               listPoints.push_back(atoi(w.c_str()));

           }

       } else if (!line.empty() && line.substr(0,15) == "#HYPERMUT-DINUC") {
           cerr << "### HYPERMUT-DINUC ###\n";
           istringstream iss(line);
           string w;
           while(iss >> w)
           {
                cerr << w << "\n";
           }

       } else if (!line.empty() && line.substr(0,5) == "#NRUN") {
           cerr << "### NRUN ###\n";
           istringstream iss(line);
           string w;
           iss >> w;
           iss >> w;
           this->Nrun = atoi(w.c_str());
           iss >> w;
           this->threshold = atoi(w.c_str());
           cerr << "Nrun " << this->Nrun << " threshold " << this->threshold << "\n";

       } else if (!line.empty() && line.substr(0,9) == "#NTHREADS") {
           cerr << "### NTHREADS ###\n";
           istringstream iss(line);
           string w;
           iss >> w;
           iss >> w;
           this->Nthread = atoi(w.c_str());
           cerr << "Nthread " << this->Nthread << "\n";

       } else if (!line.empty() && line.substr(0,7) == "#OUTPUT") {
           istringstream iss(line);
           string w;
           iss >> w;
           iss >> w;
           this->output = w;
           cerr << "#OUTPUT " << this->output << "\n";

       } else if (!line.empty() && line.substr(0,10) == "#OLDPARAMS") {
           cerr << "### OLDPARAMS ### \n";
           localcontrolfile = line;

       }
   }
   is.close();

}


