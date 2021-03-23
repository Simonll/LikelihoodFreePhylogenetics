/*
LikelihoodFreePhylogenetics, Copyright (C) 2017, Simon Laurin-Lemay

LikelihoodFreePhylogenetics is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version. LikelihoodFreePhylogenetics is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details. You should have received a copy of the GNU
General Public License along with LikelihoodFreePhylogenetics. If not, see
<http://www.gnu.org/licenses/>.
*/
#include "Posterior.h"

#include <tuple>

Posterior::Posterior(GlobalParameters* gparam) {
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

  if (verbose) {
    std::cerr << "Posterior1\n";
  }

  this->listParam = new string[this->NParam];
  for (int param_i = 0; param_i < this->NParam; param_i++) {
    this->listParam[param_i] = gparam->listParam[param_i];
  }

  if (verbose) {
    std::cerr << "Posterior2\n";
  }

  this->listSummaries = new string[this->NSummaries];
  for (int summary_i = 0; summary_i < this->NSummaries; summary_i++) {
    this->listSummaries[summary_i] = gparam->listSummaries[summary_i];
  }

  if (verbose) {
    std::cerr << "Posterior3\n";
  }

  this->listEvoStats = new string[this->NEvoStats];
  for (int EvoStats_i = 0; EvoStats_i < this->NEvoStats; EvoStats_i++) {
    this->listEvoStats[EvoStats_i] = gparam->listEvoStats[EvoStats_i];
  }

  if (verbose) {
    std::cerr << "Posterior4\n";
  }

  this->listSiteSpecificEvoStats = new string[this->NSiteSpecificEvoStats];
  for (int EvoStats_i = 0; EvoStats_i < this->NSiteSpecificEvoStats;
       EvoStats_i++) {
    this->listSiteSpecificEvoStats[EvoStats_i] =
        gparam->listSiteSpecificEvoStats[EvoStats_i];
  }

  if (verbose) {
    std::cerr << "Posterior5\n";
  }

  this->NusedEvoStats = gparam->NusedEvoStats;
  this->NusedSiteSpecificEvoStats = gparam->NusedSiteSpecificEvoStats;
  this->NusedEvoAncStats = gparam->NusedEvoAncStats;
  this->NusedParam = gparam->NusedParam;
  this->NusedSummaries = gparam->NusedSummaries;
  this->NusedAccessorySummaries = gparam->NusedAccessorySummaries;
  this->Ngenes = gparam->Ngenes;

  if (verbose) {
    std::cerr << "Posterior6\n";
  }

  this->mapUsedParam.insert(gparam->mapUsedParam.begin(),
                            gparam->mapUsedParam.end());
  this->mapUsedSummaries.insert(gparam->mapUsedSummaries.begin(),
                                gparam->mapUsedSummaries.end());
  this->mapUsedAccessorySummaries.insert(
      gparam->mapUsedAccessorySummaries.begin(),
      gparam->mapUsedAccessorySummaries.end());
  this->mapUsedEvoStats.insert(gparam->mapUsedEvoStats.begin(),
                               gparam->mapUsedEvoStats.end());
  this->mapUsedSiteSpecificEvoStats.insert(
      gparam->mapUsedSiteSpecificEvoStats.begin(),
      gparam->mapUsedSiteSpecificEvoStats.end());
  this->mapUsedEvoAncStats.insert(gparam->mapUsedEvoAncStats.begin(),
                                  gparam->mapUsedEvoAncStats.end());

  if (verbose) {
    std::cerr << "Posterior7\n";
  }

  this->sorted = false;
  this->Naccepted = 0;
  this->randomseed = -1;
  this->rnd = new Random(randomseed);

  population_t.reserve(this->threshold);

  empVar = new double[this->NusedParam];
  empMean = new double[this->NusedParam];
  for (int i = 0; i < this->NusedParam; i++) {
    empVar[i] = 1.0;
    empMean[i] = 1.0;
  }
}

Posterior::~Posterior() {
  // dtor
}

int Posterior::PosteriorGetSize() { return posterior.size(); }

std::vector<std::vector<double>> Posterior::GetPartialDistances() {
  std::vector<std::vector<double>> X;
  for (int simu_i = 0; simu_i < static_cast<int>(population_t.size());
       simu_i++) {
    std::vector<double> X_i;
    for (int distance_i = 0; distance_i < NusedSummaries + 1; distance_i++) {
      if (distance_i == 0) {
        X_i.push_back(1.0);
      } else {
        X_i.push_back(
            std::get<distancesGetter>(population_t[simu_i])[distance_i]);
      }
    }
    X.push_back(X_i);
  }
  return X;
}

std::vector<std::vector<double>> Posterior::GetPartialDistancesT() {
  std::vector<std::vector<double>> X_T;
  for (int distance_i = 0; distance_i < NusedSummaries + 1; distance_i++) {
    std::vector<double> X_Ti;
    for (int simu_i = 0; simu_i < static_cast<int>(population_t.size());
         simu_i++) {
      if (distance_i == 0) {
        X_Ti[simu_i] = 1.0;
      } else {
        X_Ti.push_back(
            std::get<distancesGetter>(population_t[simu_i])[distance_i]);
      }
    }
    X_T.push_back(X_Ti);
  }
  return X_T;
}

std::vector<std::vector<double>> Posterior::GetTheta() {
  std::vector<std::vector<double>> theta;
  for (int simu_i = 0; simu_i < static_cast<int>(population_t.size());
       simu_i++) {
    theta.push_back(std::get<paramGetter>(population_t[simu_i]));
  }
  return theta;
}

double Posterior::GetEpanechnikov(double x, double y) {
  return 1 - (x * x) / (y * y);
}

std::vector<double> Posterior::GetWeights() {
  std::vector<double> W;
  double max_ = std::get<distancesGetter>(
      population_t[population_t.size() - 1])[NusedSummaries];
  for (int simu_i = 0; simu_i < static_cast<int>(population_t.size());
       simu_i++) {
    W.push_back(GetEpanechnikov(
        std::get<distancesGetter>(population_t[simu_i])[NusedSummaries], max_));
  }
  return W;
}

std::vector<std::vector<double>> Posterior::GetLocalWeights() {
  std::vector<std::vector<double>> W;

  double* max_ = new double[NusedSummaries];
  for (int distance_i = 0; distance_i < NusedSummaries; distance_i++) {
    max_[distance_i] = std::get<distancesGetter>(
        population_t[population_t.size() - 1])[distance_i];
  }
  for (int simu_i = 0; simu_i < static_cast<int>(population_t.size());
       simu_i++) {
    std::vector<double> W_i;
    for (int distance_i = 0; distance_i < NusedSummaries; distance_i++) {
      W_i.push_back(GetEpanechnikov(
          std::get<distancesGetter>(population_t[simu_i])[distance_i],
          max_[distance_i]));
    }
    W.push_back(W_i);
  }
  return W;
}

void Posterior::writePosterior(ofstream& os, int Nsimu) {
  for (int simu_i = 0; simu_i < Nsimu; simu_i++) {
    // write chainID
    os << std::get<chainIDGetter>(population_t[simu_i]) << "\t";

    // write parameters
    if (this->NusedParam > 0) {
      for (int param_i = 0;
           param_i <
           static_cast<int>(std::get<paramGetter>(population_t[simu_i]).size());
           param_i++) {
        os << std::get<paramGetter>(population_t[simu_i])[param_i] << "\t";
      }
    }
    // write summaries
    if (this->NusedSummaries > 0) {
      for (int summary_i = 0;
           summary_i <
           static_cast<int>(
               std::get<summariesGetter>(population_t[simu_i]).size());
           summary_i++) {
        os << std::get<summariesGetter>(population_t[simu_i])[summary_i]
           << "\t";
      }
    }
    // write distances (sum of square discrepancies)
    if (this->NusedSummaries > 0) {
      if (this->OutPartialDistance) {
        for (int distance_i = 0;
             distance_i <
             static_cast<int>(
                 std::get<distancesGetter>(population_t[simu_i]).size());
             distance_i++) {
          os << std::get<distancesGetter>(population_t[simu_i])[distance_i]
             << "\t";
        }
      } else {
        os << std::get<distancesGetter>(population_t[simu_i])
                  [std::get<distancesGetter>(population_t[simu_i]).size() - 1]
           << "\t";
      }
    }
    // write accessory summaries
    if (this->NusedAccessorySummaries > 0) {
      for (int summary_i = 0;
           summary_i <
           static_cast<int>(
               std::get<accsummariesGetter>(population_t[simu_i]).size());
           summary_i++) {
        os << std::get<accsummariesGetter>(population_t[simu_i])[summary_i]
           << "\t";
      }
    }
    // write mappingstats
    if (this->NusedEvoAncStats > 0) {
      for (int mapping_i = 0;
           mapping_i <
           static_cast<int>(
               std::get<evoancstatsGetter>(population_t[simu_i]).size());
           mapping_i++) {
        if (mapping_i <
            static_cast<int>(
                std::get<evoancstatsGetter>(population_t[simu_i]).size()) -
                1) {
          os << std::get<evoancstatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        } else {
          os << std::get<evoancstatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        }
      }
    }
    // write mappingstats
    if (this->NusedEvoStats > 0) {
      for (int mapping_i = 0;
           mapping_i <
           static_cast<int>(
               std::get<evostatsGetter>(population_t[simu_i]).size());
           mapping_i++) {
        if (mapping_i <
            static_cast<int>(
                std::get<evostatsGetter>(population_t[simu_i]).size()) -
                1) {
          os << std::get<evostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        } else {
          os << std::get<evostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        }
      }
    }
    // write mappingstats
    if (this->NusedSiteSpecificEvoStats > 0) {
      for (int mapping_i = 0;
           mapping_i <
           static_cast<int>(
               std::get<ssevostatsGetter>(population_t[simu_i]).size());
           mapping_i++) {
        if (mapping_i <
            static_cast<int>(
                std::get<ssevostatsGetter>(population_t[simu_i]).size()) -
                1) {
          os << std::get<ssevostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        } else {
          os << std::get<ssevostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        }
      }
    }
    os << "\n";
  }
}

void Posterior::writePosterior(ofstream& os) {
  for (int simu_i = 0; simu_i < this->threshold; simu_i++) {
    // write chainID
    os << std::get<chainIDGetter>(population_t[simu_i]) << "\t";

    // write parameters
    if (this->NusedParam > 0) {
      for (int param_i = 0;
           param_i <
           static_cast<int>(std::get<paramGetter>(population_t[simu_i]).size());
           param_i++) {
        os << std::get<paramGetter>(population_t[simu_i])[param_i] << "\t";
      }
    }
    // write summaries
    if (this->NusedSummaries > 0) {
      for (int summary_i = 0;
           summary_i <
           static_cast<int>(
               std::get<summariesGetter>(population_t[simu_i]).size());
           summary_i++) {
        os << std::get<summariesGetter>(population_t[simu_i])[summary_i]
           << "\t";
      }
    }
    // write distances (sum of square discrepancies)
    if (this->NusedSummaries > 0) {
      if (this->OutPartialDistance) {
        for (int distance_i = 0;
             distance_i <
             static_cast<int>(
                 std::get<distancesGetter>(population_t[simu_i]).size());
             distance_i++) {
          os << std::get<distancesGetter>(population_t[simu_i])[distance_i]
             << "\t";
        }
      } else {
        os << std::get<distancesGetter>(population_t[simu_i])
                  [std::get<distancesGetter>(population_t[simu_i]).size() - 1]
           << "\t";
      }
    }
    // write accessory summaries
    if (this->NusedAccessorySummaries > 0) {
      for (int summary_i = 0;
           summary_i <
           static_cast<int>(
               std::get<accsummariesGetter>(population_t[simu_i]).size());
           summary_i++) {
        os << std::get<accsummariesGetter>(population_t[simu_i])[summary_i]
           << "\t";
      }
    }
    // write mappingstats
    if (this->NusedEvoAncStats > 0) {
      for (int mapping_i = 0;
           mapping_i <
           static_cast<int>(
               std::get<evoancstatsGetter>(population_t[simu_i]).size());
           mapping_i++) {
        if (mapping_i <
            static_cast<int>(
                std::get<evoancstatsGetter>(population_t[simu_i]).size()) -
                1) {
          os << std::get<evoancstatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        } else {
          os << std::get<evoancstatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        }
      }
    }
    // write mappingstats
    if (this->NusedEvoStats > 0) {
      for (int mapping_i = 0;
           mapping_i <
           static_cast<int>(
               std::get<evostatsGetter>(population_t[simu_i]).size());
           mapping_i++) {
        if (mapping_i <
            static_cast<int>(
                std::get<evostatsGetter>(population_t[simu_i]).size()) -
                1) {
          os << std::get<evostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        } else {
          os << std::get<evostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        }
      }
    }
    // write mappingstats
    if (this->NusedSiteSpecificEvoStats > 0) {
      for (int mapping_i = 0;
           mapping_i <
           static_cast<int>(
               std::get<ssevostatsGetter>(population_t[simu_i]).size());
           mapping_i++) {
        if (mapping_i <
            static_cast<int>(
                std::get<ssevostatsGetter>(population_t[simu_i]).size()) -
                1) {
          os << std::get<ssevostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        } else {
          os << std::get<ssevostatsGetter>(population_t[simu_i])[mapping_i]
             << "\t";
        }
      }
    }
    os << "\n";
  }
}

void Posterior::readPosterior(string posteriorfile) {
  ifstream is(posteriorfile.c_str());
  if (!is) {
    std::cerr << "error: did not find " << posteriorfile << "\n";
    exit(1);
  }
  string line;
  std::getline(is, line);
  istringstream iss(line);
  string w;
  int k = 0;
  string* arrParam = new string[this->NusedParam];
  while (iss >> w) {
    auto it = mapUsedParam.find(w);
    if (it != mapUsedParam.end()) {
      if (it->second == -1) {
        std::cerr << "Undefined parameter " << w << "\n";
        std::cerr << "Not present in the posterior\n";
        exit(0);
      } else {
        std::cerr << it->first << " " << it->second << " " << w << "\n";
        if (k == this->NusedParam) {
          std::cerr << "Wrong number of parameters " << k << " "
                    << this->NusedParam << "\n";
          exit(0);
        }
        arrParam[k] = w;
        k++;
      }
    } else {
      std::cerr << "Undefined parameter " << w << "\n";
      exit(0);
    }
  }

  if (k != NusedParam) {
    std::cerr << "Different number of parameters between posterior and "
                 "configuration file\n";
    exit(0);
  }
  std::cerr << "NusedParam " << this->NusedParam << " " << k << "\n";

  // transform for global order
  std::map<string, int> map_neworder;
  k = 0;
  for (int param_i = 0; param_i < this->NParam; param_i++) {
    auto it_ = this->mapUsedParam.find(this->listParam[param_i]);
    if (it_ != this->mapUsedParam.end()) {
      if (it_->second != -1) {
        for (int param_j = 0; param_j < this->NusedParam; param_j++) {
          if (it_->first == arrParam[param_j]) {
            std::cerr << param_j << " " << arrParam[param_j] << " " << k
                      << "\n";
            map_neworder.insert({arrParam[param_j], k});
            k++;
          }
        }
      }
    }
  }
  while (std::getline(is, line)) {
    if (!line.empty()) {
      // std::cerr << line << "\n";
      istringstream iss(line);
      std::vector<double> cur_param;
      for (int param_i = 0; param_i < NusedParam; param_i++) {
        cur_param.push_back(-1);
      }
      double d = 0.0;
      for (int param_i = 0; param_i < NusedParam; param_i++) {
        auto it_ = map_neworder.find(arrParam[param_i]);
        if (it_ != map_neworder.end()) {
          iss >> d;
          cur_param[it_->second] = d;
        }
      }
      posterior.push_back(cur_param);
    }
  }
  delete[] arrParam;
  is.close();
}

void Posterior::readPosterior(ifstream& is) {
  is.clear();  // clear fail and eof bits
  is.seekg(0, std::ios::beg);
  int verbose = 0;
  // check correspondence between control file and posterior file
  if (verbose) {
    std::cerr << "readPosterior NusedParam" << this->NusedParam << "\n";
  }

  std::map<int, string> mapHeader;
  int mapHeaderIndex = 0;
  if (this->NusedParam > 0) {
    is.clear();                  // clear fail and eof bits
    is.seekg(0, std::ios::beg);  // back to the start!
    string line;
    std::getline(is, line);
    istringstream iss(line);
    string w;
    int k = 0;
    string* arr = new string[this->NusedParam];
    while (iss >> w || k < this->NusedParam) {
      auto it = mapUsedParam.find(w);
      if (it != mapUsedParam.end()) {
        if (it->second != -1) {
          if (k != it->second) {
            std::cerr << k << " " << it->second << " " << w << " " << it->first
                      << "\n";
            exit(0);
          }
          if (w == "chainID") {
            arr[it->second] = it->first;
            mapHeader[mapHeaderIndex] = "chainID";
            k++;
            mapHeaderIndex++;
          } else {
            arr[it->second] = it->first;
            mapHeader[mapHeaderIndex] = "P";
            k++;
            mapHeaderIndex++;
          }
        }
      }
    }
    for (int v = 0; v < this->NusedParam; v++) {
      std::cerr << "P " << arr[v] << "\n";
    }
    delete[] arr;
  }

  if (verbose) {
    std::cerr << "readPosterior NusedSummaries" << this->NusedSummaries << "\n";
  }

  if (this->NusedSummaries > 0) {
    is.clear();                  // clear fail and eof bits
    is.seekg(0, std::ios::beg);  // back to the start!
    string line;
    std::getline(is, line);
    istringstream iss(line);
    string w;
    int k = 0;
    string* arr = new string[this->NusedSummaries + 1];
    while (iss >> w || k < this->NusedSummaries + 1) {
      auto it = mapUsedSummaries.find(w);
      if (it != mapUsedSummaries.end()) {
        if (it->second != -1) {
          if (k != it->second) {
            std::cerr << k << " " << it->second << " " << w << " " << it->first
                      << "\n";
            exit(0);
          }
          arr[it->second] = it->first;
          mapHeader[mapHeaderIndex] = "S";
          k++;
          mapHeaderIndex++;
        }
      } else if (w == "D_sum") {
        if (k != this->NusedSummaries) {
          std::cerr << k << " " << this->NusedSummaries << " " << w << " "
                    << "\n";
          exit(0);
        }
        arr[this->NusedSummaries] = "D_sum";
        mapHeader[mapHeaderIndex] = "D_sum";
        k++;
        mapHeaderIndex++;
      }
    }
    for (int v = 0; v < this->NusedSummaries; v++) {
      std::cerr << "S " << arr[v] << "\n";
    }
    std::cerr << "D_sum " << arr[this->NusedSummaries] << "\n";
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
  //                        std::cerr << k << " "<< it->second << " " << w << "
  //                        " << it->first << "\n"; exit(0);
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
  //            std::cerr << "D_"<< arr[v] << "\n";
  //        }
  //        delete[] arr;
  //    }

  if (verbose) {
    std::cerr << "readPosterior NusedEvoStats" << this->NusedEvoStats << "\n";
  }

  if (this->NusedEvoStats > 0) {
    is.clear();                  // clear fail and eof bits
    is.seekg(0, std::ios::beg);  // back to the start!
    string line;
    std::getline(is, line);
    istringstream iss(line);
    string w;
    int k = 0;
    string* arr = new string[this->NusedEvoStats];
    while (iss >> w || k < this->NusedEvoStats) {
      w = w.erase(0, 2);
      auto it = mapUsedEvoStats.find(w);
      if (it != mapUsedEvoStats.end()) {
        if (it->second != -1) {
          if (k != it->second) {
            std::cerr << k << " " << it->second << " " << w << " " << it->first
                      << "\n";
            exit(0);
          }
          arr[it->second] = it->first;
          mapHeader[mapHeaderIndex] = "ES";
          k++;
          mapHeaderIndex++;
        }
      }
    }
    for (int v = 0; v < this->NusedEvoStats; v++) {
      std::cerr << "ES " << arr[v] << "\n";
    }
    delete[] arr;
  }

  // get all the information

  is.clear();
  is.seekg(0, std::ios::beg);
  string line;
  std::getline(is, line);  // skip header
  while (std::getline(is, line)) {
    std::cerr << ".";
    if (!line.empty()) {
      // std::cerr << line << "\n";

      istringstream iss_tmp(line);
      int chainID = 1;
      std::vector<double> cur_param;
      std::vector<double> cur_summaries;
      std::vector<double> cur_evostats;
      std::vector<double> cur_distances;
      mapHeaderIndex = 0;
      string w;

      while (iss_tmp >> w) {
        auto it = mapHeader.find(mapHeaderIndex);
        if (it != mapHeader.end()) {
          if (it->second == "P") {
            cur_param.push_back(std::stof(w));

          } else if (it->second == "chainID") {
            chainID = std::stoi(w);
          } else if (it->second == "S") {
            cur_summaries.push_back(std::stof(w));

          } else if (it->second == "ES") {
            cur_evostats.push_back(std::stof(w));

          } else if (it->second == "D_sum") {
            cur_distances.push_back(std::stof(w));
          }
          mapHeaderIndex++;
        }
      }

      std::vector<double> tmp;
      registerOldSimulation(chainID, cur_param, cur_summaries, tmp, tmp,
                            cur_evostats, tmp, cur_distances, tmp);
    }
  }
}

void Posterior::writePosteriorPredictiveStatistics(
    ofstream& os, std::vector<double> summariesRealData) {
  int pop_size = static_cast<int>(population_t.size());
  int k = 0;
  if (NusedSummaries > 0) {
    string* arrSummaries = new string[this->NusedSummaries];
    for (int summary_i = 0; summary_i < this->NSummaries; summary_i++) {
      auto it = mapUsedSummaries.find(this->listSummaries[summary_i]);
      if (it != mapUsedSummaries.end()) {
        if (it->second != -1) {
          arrSummaries[it->second] = it->first;
        }
      }
    }
    for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
      if (k == 0) {
        os << arrSummaries[summary_i];
        k = 1;
      } else {
        os << "\t" << arrSummaries[summary_i];
      }
    }
    os << "\n";
  }

  double* ppp = new double[this->NusedSummaries];
  double* mean = new double[this->NusedSummaries];
  double* var = new double[this->NusedSummaries];
  double* a = new double[this->NusedSummaries];
  double* b = new double[this->NusedSummaries];
  double* K = new double[this->NusedSummaries];

  for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
    ppp[summary_i] = 0.0;
    mean[summary_i] = 0.0;
    var[summary_i] = 0.0;
    a[summary_i] = 0.0;
    b[summary_i] = 0.0;
    K[summary_i] = std::get<summariesGetter>(
        population_t[0])[summary_i];  // const for shifted data
  }

  for (int simu_i = 0; simu_i < pop_size; simu_i++) {
    for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
      if (std::get<summariesGetter>(population_t[simu_i])[summary_i] <
          summariesRealData[summary_i]) {
        ppp[summary_i] += 1;
      }
      mean[summary_i] +=
          std::get<summariesGetter>(population_t[simu_i])[summary_i];
      a[summary_i] +=
          std::get<summariesGetter>(population_t[simu_i])[summary_i] -
          K[summary_i];
      b[summary_i] +=
          (std::get<summariesGetter>(population_t[simu_i])[summary_i] -
           K[summary_i]) *
          (std::get<summariesGetter>(population_t[simu_i])[summary_i] -
           K[summary_i]);
    }
  }
  for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
    mean[summary_i] /= pop_size;
    var[summary_i] =
        (b[summary_i] - (a[summary_i] * a[summary_i]) / pop_size) / (pop_size);
  }
  // posterior predictive p-values (ppp)
  for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
    ppp[summary_i] /= pop_size;
    if (summary_i < this->NusedSummaries - 1) {
      os << ppp[summary_i] << "\t";
    } else {
      os << ppp[summary_i] << "\n";
    }
  }

  // means
  for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
    if (summary_i < this->NusedSummaries - 1) {
      os << mean[summary_i] << "\t";
    } else {
      os << mean[summary_i] << "\n";
    }
  }

  // variances
  for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
    if (summary_i < this->NusedSummaries - 1) {
      os << var[summary_i] << "\t";
    } else {
      os << var[summary_i] << "\n";
    }
  }

  // z-scores
  for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
    if (summary_i < this->NusedSummaries - 1) {
      os << (summariesRealData[summary_i] - mean[summary_i]) /
                sqrt(var[summary_i])
         << "\t";
    } else {
      os << (summariesRealData[summary_i] - mean[summary_i]) /
                sqrt(var[summary_i])
         << "\n";
    }
  }
  delete[] ppp;
  delete[] mean;
  delete[] var;
  delete[] a;
  delete[] b;
  delete[] K;
}
void Posterior::readMonitor(ifstream& is) {
  is.clear();  // clear fail and eof bits
  is.seekg(0, std::ios::beg);

  string line;
  std::getline(is, line);
  if (!line.empty()) {
    istringstream iss1(line);
    iss1 >> this->Niter;
    std::getline(is, line);
    istringstream iss2(line);
    iss2 >> this->Naccepted;
    std::cerr << "Niter " << this->Niter << " "
              << "Naccepted " << this->Naccepted << "\n";

  } else {
    this->Niter = 0;
    std::cerr << "Monitor file is empty "
              << "\n";
  }
}

void Posterior::writeMonitorPosterior(ofstream& os) {
  os << this->Niter << "\n";
  os << this->Naccepted << "\n";
  os << GetAcceptanceRate() << "\n";
}

double Posterior::GetAcceptanceRate() {
  return static_cast<double>(Naccepted / Niter);
}

void Posterior::writeHeader(ofstream& os) {
  // write parameters' header
  int k = 0;
  int v = 0;
  if (verbose) {
    std::cerr << "writeHeader1 " << this->NusedParam << "\n";
  }
  if (this->NusedParam > 0) {
    string* arrParam = new string[this->NusedParam];
    for (int param_i = 0; param_i < this->NParam; param_i++) {
      auto it = this->mapUsedParam.find(this->listParam[param_i]);
      if (it != this->mapUsedParam.end()) {
        if (it->second != -1) {
          arrParam[v] = it->first;
          v++;
        }
      }
    }

    for (int param_i = 0; param_i < this->NusedParam; param_i++) {
      if (k == 0) {
        os << arrParam[param_i];
        k = 1;
      } else {
        os << "\t" << arrParam[param_i];
      }
    }
    delete[] arrParam;
  }

  if (verbose) {
    std::cerr << "writeHeader2 " << this->NusedSummaries << "\n";
  }
  // write summaires' header
  if (this->NusedSummaries > 0) {
    string* arrSummaries = new string[this->NusedSummaries];
    for (int summary_i = 0; summary_i < this->NSummaries; summary_i++) {
      auto it = this->mapUsedSummaries.find(this->listSummaries[summary_i]);
      if (it != this->mapUsedSummaries.end()) {
        if (it->second != -1) {
          arrSummaries[it->second] = it->first;
        }
      }
    }

    for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
      if (k == 0) {
        os << arrSummaries[summary_i];
        k = 1;
      } else {
        os << "\t" << arrSummaries[summary_i];
      }
    }

    // write distance header
    if (this->OutPartialDistance) {
      for (int summary_i = 0; summary_i < this->NusedSummaries; summary_i++) {
        if (k == 0) {
          os << "D_" << arrSummaries[summary_i];
          k = 1;
        } else {
          os << "\t"
             << "D_" << arrSummaries[summary_i];
        }
      }
    }
    os << "\tD_sum";
    delete[] arrSummaries;
  }

  if (this->NusedAccessorySummaries > 0) {
    string* arrSummaries = new string[this->NusedAccessorySummaries];
    for (int summary_i = 0; summary_i < this->NSummaries; summary_i++) {
      auto it =
          this->mapUsedAccessorySummaries.find(this->listSummaries[summary_i]);
      if (it != this->mapUsedAccessorySummaries.end()) {
        if (it->second != -1) {
          arrSummaries[it->second] = it->first;
        }
      }
    }

    for (int summary_i = 0; summary_i < this->NusedAccessorySummaries;
         summary_i++) {
      if (k == 0) {
        os << "Acc_" << arrSummaries[summary_i];
        k = 1;
      } else {
        os << "\t"
           << "Acc_" << arrSummaries[summary_i];
      }
    }
    delete[] arrSummaries;
  }

  // write mappingstats
  if (verbose) {
    std::cerr << "writeHeader3 " << this->NusedEvoAncStats << "\n";
  }
  if (this->NusedEvoAncStats > 0) {
    string* arrMapAnc = new string[this->NusedEvoAncStats];
    for (int map_i = 0; map_i < this->NEvoStats; map_i++) {
      auto it = this->mapUsedEvoAncStats.find(this->listEvoStats[map_i]);
      if (it != this->mapUsedEvoAncStats.end()) {
        if (it->second != -1) {
          arrMapAnc[it->second] = it->first;
        }
      }
    }

    for (int map_i = 0; map_i < NusedEvoAncStats; map_i++) {
      if (k == 0) {
        os << "A_" << arrMapAnc[map_i];
        k = 1;
      } else {
        os << "\t"
           << "A_" << arrMapAnc[map_i];
      }
    }
    delete[] arrMapAnc;
  }

  if (verbose) {
    std::cerr << "writeHeader4 " << this->NusedEvoStats << "\n";
  }
  if (this->NusedEvoStats > 0) {
    string* arrMap = new string[this->NusedEvoStats];
    for (int map_i = 0; map_i < this->NEvoStats; map_i++) {
      auto it = this->mapUsedEvoStats.find(this->listEvoStats[map_i]);
      if (it != this->mapUsedEvoStats.end()) {
        if (it->second != -1) {
          arrMap[it->second] = it->first;
        }
      }
    }

    for (int map_i = 0; map_i < this->NusedEvoStats; map_i++) {
      if (k == 0) {
        os << "T_" << arrMap[map_i];
        k = 1;
      } else {
        os << "\t"
           << "T_" << arrMap[map_i];
      }
    }
    delete[] arrMap;
  }

  if (verbose) {
    std::cerr << "writeHeader5 " << this->NusedSiteSpecificEvoStats << "\n";
  }
  if (this->NusedSiteSpecificEvoStats > 0) {
    string* arrMap = new string[this->NusedSiteSpecificEvoStats];
    for (int map_i = 0; map_i < this->NSiteSpecificEvoStats; map_i++) {
      auto it = this->mapUsedSiteSpecificEvoStats.find(
          this->listSiteSpecificEvoStats[map_i]);
      if (it != this->mapUsedSiteSpecificEvoStats.end()) {
        if (it->second != -1) {
          arrMap[it->second] = it->first;
        }
      }
    }

    for (int map_i = 0; map_i < this->NusedSiteSpecificEvoStats; map_i++) {
      for (int bin_i = 0; bin_i < 100; bin_i++) {
        if (k == 0) {
          os << "SS_" << bin_i << "_" << arrMap[map_i];
          k = 1;
        } else {
          os << "\t"
             << "SS_" << bin_i << "_" << arrMap[map_i];
        }
      }
    }
    delete[] arrMap;
  }
  os << "\n";
}

void Posterior::SetNsite(int i) { this->Nsite_codon = i; }

void Posterior::GetWeights(string kernel) {
  int pop_size = population_t.size();
  if (kernel == "sNormal") {
    for (int param_i = 0; param_i < NusedParam; param_i++) {
      double new_weight = 0.0;
      for (int simu_i = 0; simu_i < pop_size; simu_i++) {
        new_weight += std::get<weightsGetter>(population_t[simu_i])[param_i] *
                      (std::get<paramGetter>(population_t[simu_i])[param_i] +
                       2.0 * (empVar[param_i]) * rnd->sNormal());
      }
      new_weight /= pop_size;
      //            std::get<4>(population_t[simu_i])[param_i] = new_weight;
    }
  }
}

void Posterior::GetEmpVar() {
  int pop_size = population_t.size();

  double* mean = new double[NusedParam];
  double* sum1 = new double[NusedParam];
  double* sum2 = new double[NusedParam];
  double* sum3 = new double[NusedParam];

  for (int param_i = 0; param_i < NusedParam; param_i++) {
    mean[param_i] = 0.0;
    sum1[param_i] = 0.0;
    sum2[param_i] = 0.0;
    sum3[param_i] = 0.0;
  }

  for (int simu_i = 0; simu_i < NusedParam; simu_i++) {
    for (int param_i = 0; param_i < NusedParam; param_i++) {
      sum1[param_i] += std::get<paramGetter>(population_t[simu_i])[param_i];
    }
  }

  for (int param_i = 0; param_i < NusedParam; param_i++) {
    mean[param_i] = sum1[param_i] / pop_size;
    empMean[param_i] = mean[param_i];
  }
  for (int simu_i = 0; simu_i < pop_size; simu_i++) {
    for (int param_i = 0; param_i < NusedParam;
         param_i++) {  // (x - K) * (x - K)
      sum2[param_i] += (std::get<paramGetter>(population_t[simu_i])[param_i] -
                        mean[param_i]) *
                       (std::get<paramGetter>(population_t[simu_i])[param_i] -
                        mean[param_i]);
      // x - K
      sum3[param_i] +=
          std::get<paramGetter>(population_t[simu_i])[param_i] - mean[param_i];
    }
  }

  for (int param_i = 0; param_i < NusedParam; param_i++) {
    // sum_sqr - (sum_ * sum_)/n)/(n - 1)
    empVar[param_i] =
        (sum2[param_i] - (sum3[param_i] * sum3[param_i]) / pop_size) /
        (pop_size - 1);
  }

  delete[] mean;
  delete[] sum1;
  delete[] sum2;
  delete[] sum3;
}

void Posterior::slaveToMaster(
    std::vector<std::tuple<int, std::vector<double>, std::vector<double>,
                           std::vector<double>, std::vector<double>,
                           std::vector<double>, std::vector<double>,
                           std::vector<double>, std::vector<double>>>
        population_i) {
  population_t.insert(population_t.end(), population_i.begin(),
                      population_i.end());
}

void Posterior::sortPopulation() {
  if (!sorted) {
    std::sort(
        population_t.begin(), population_t.end(),
        [](const std::tuple<int, std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>>& left,
           const std::tuple<int, std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>>& right) {
          return std::get<distancesGetter>(left).back() <
                 std::get<distancesGetter>(right).back();
        });
    sorted = true;
  }
}

void Posterior::slaveRegisterNewSimulation(
    int chainID, std::vector<double> param, std::vector<double> summaries,
    std::vector<double> accsummaries, std::vector<double> evoancstat,
    std::vector<double> evostat, std::vector<double> ssevostat,
    std::vector<double> distances, std::vector<double> weights) {
  population_t.push_back(make_tuple(chainID, param, summaries, accsummaries,
                                    evoancstat, evostat, ssevostat, distances,
                                    weights));
  Naccepted++;
  this->Niter++;
}

void Posterior::registerNewSimulation(
    int chainID, std::vector<double> param, std::vector<double> summaries,
    std::vector<double> accsummaries, std::vector<double> evoancstat,
    std::vector<double> evostat, std::vector<double> ssevostat,
    std::vector<double> distances, std::vector<double> weights) {
  if (population_t.empty()) {
    // std::cerr << "POPULATION IS EMPTY" << Naccepted <<  " "<< threshold <<  "
    // \n";
    population_t.push_back(make_tuple(chainID, param, summaries, accsummaries,
                                      evoancstat, evostat, ssevostat, distances,
                                      weights));
    Naccepted++;
  } else if (static_cast<int>(population_t.size()) < threshold) {
    // std::cerr << "POPULATION IS LESS THAN THRESHOLD" << " " <<
    // population_t.size()
    // << " "<< Naccepted  << " "<< threshold << " \n";
    population_t.push_back(make_tuple(chainID, param, summaries, accsummaries,
                                      evoancstat, evostat, ssevostat, distances,
                                      weights));
    Naccepted++;
  } else {
    if (!sorted) {
      std::sort(
          population_t.begin(), population_t.end(),
          [](const std::tuple<int, std::vector<double>, std::vector<double>,
                              std::vector<double>, std::vector<double>,
                              std::vector<double>, std::vector<double>,
                              std::vector<double>, std::vector<double>>& left,
             const std::tuple<int, std::vector<double>, std::vector<double>,
                              std::vector<double>, std::vector<double>,
                              std::vector<double>, std::vector<double>,
                              std::vector<double>, std::vector<double>>&
                 right) {
            return std::get<distancesGetter>(left).back() <
                   std::get<distancesGetter>(right).back();
          });
      sorted = true;
    }
    std::tuple<int, std::vector<double>, std::vector<double>,
               std::vector<double>, std::vector<double>, std::vector<double>,
               std::vector<double>, std::vector<double>, std::vector<double>>
        cur_tuple =
            make_tuple(chainID, param, summaries, accsummaries, evoancstat,
                       evostat, ssevostat, distances, weights);

    auto it = std::lower_bound(
        population_t.begin(), population_t.end(), cur_tuple,
        [](const std::tuple<int, std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>>& left,
           const std::tuple<int, std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>,
                            std::vector<double>, std::vector<double>>& right) {
          return std::get<distancesGetter>(left).back() <
                 std::get<distancesGetter>(right).back();
        });
    // compute acceptance rate
    if (it != population_t.end()) {
      // std::cerr << "INSERTED" << " " << population_t.size()  << " "<<
      // Naccepted
      // << " "<< threshold << " \n";
      population_t.insert(it, cur_tuple);
      population_t.pop_back();
      population_t.shrink_to_fit();
      Naccepted++;
    }
  }
  this->Niter++;
  // std::cerr << this->Niter << " ";
}

void Posterior::registerOldSimulation(
    int chainID, std::vector<double> param, std::vector<double> summaries,
    std::vector<double> accsummaries, std::vector<double> evoancstat,
    std::vector<double> evostat, std::vector<double> ssevostat,
    std::vector<double> distances, std::vector<double> weights) {
  population_t.push_back(make_tuple(chainID, param, summaries, accsummaries,
                                    evoancstat, evostat, ssevostat, distances,
                                    weights));
}
