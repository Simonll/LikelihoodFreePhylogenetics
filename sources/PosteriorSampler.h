#ifndef POSTERIORSAMPLER_H
#define POSTERIORSAMPLER_H

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
#include "LocalParameters.h"

class PosteriorSampler
{
public:
    LocalParameters* lparam;
    PosteriorSampler();
    virtual ~PosteriorSampler();

protected:

    void sample();
    std::tuple<double,double> sampleFromAjustedDensity(int param,double inf, double sup);

private:
};

#endif // POSTERIORSAMPLER_H
