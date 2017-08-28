#ifndef SAMPLER_H
#define SAMPLER_H

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


class PriorSampler
{
public:
    LocalParameters* lparam;

    PriorSampler(LocalParameters* lparam);
    virtual ~PriorSampler();

    void sample();

protected:

    double Unif(double inf, double sup)
    {
        return lparam->rnd->Uniform()* (sup-inf) + inf;
    }

    double log2Unif()
    {
        return pow(2,lparam->rnd->Uniform()*2-1);
    }


    double log10Unif()
    {
        return pow(10,lparam->rnd->Uniform()*2-1);
    }
    
    double log50Unif()
    {
        return pow(50,lparam->rnd->Uniform()*2-1);
    }


    double logUnifTruncated()
    {
        double a = pow (100,lparam->rnd->Uniform()*2-1);
        // 1 < x <= 4
        if ( a <= 20.0 && a > 5.0 )
        {

            a /= 5;

            // 4 < x <= 20
        }
        else if ( a > 20 )
        {

            a /= 5;

        }
        return a;
    }


private:
};

#endif // SAMPLER_H
