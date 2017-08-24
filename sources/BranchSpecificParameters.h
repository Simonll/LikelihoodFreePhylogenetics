#ifndef BRANCHSPECIFICPARAMETERS_H
#define BRANCHSPECIFICPARAMETERS_H

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



class BranchSpecificParameters
{
public:
    static const int Nnucp = 4;
    static const int Nnucrr = 6;
    static const int  Ndinuc = 16;
    double* nucp, *nucrr; // 6 param
    double** nucrrnr; //12 param
    double** gtnr;

    BranchSpecificParameters();
    virtual ~BranchSpecificParameters();

    void SetLocalParaemters(double* nucp, double* nucrr);


protected:

private:
};

#endif // BRANCHSPECIFICPARAMETERS_H
