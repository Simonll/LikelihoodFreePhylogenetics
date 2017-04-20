#ifndef REGRESSION_H
#define REGRESSION_H

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

#include "GlobalParameters.h"
#include "linalg.h"
#include "Posterior.h"

class Regression
{
    public:


        GlobalParameters* gparam;
        Posterior* post;

        //local regression adjustment
        //(X_T dot W dot X)^{-1} dot X_T dot W dot Y
        std::vector<std::vector<double>> Y_HAT;
        std::vector<std::vector<double>> Yres; //Matrix of residual = Y - Y_HAT
        std::vector<std::vector<double>> B_HAT; //Matrix of regression coefficient




        Regression(GlobalParameters* gparam, Posterior* post);
        virtual ~Regression();

        //Getters
        double GetRSquare(int i);
        double GetAjustedRSquare(int i);
        double GetRSquarePred(int i);
        double GetS(int i);
        double** GetXWX_t();
        double** GetX_tWY();

        void ComputeMultipleRegression();
        void ComputeX_TWY();


    protected:

    private:
};

#endif // REGRESSION_H
