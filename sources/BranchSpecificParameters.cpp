#include "BranchSpecificParameters.h"

BranchSpecificParameters::BranchSpecificParameters()
{
        nucrrnr = new double*[Nnucp];
        gtnr = new double*[Nnucp];
        for (int i = 0; i < Nnucp; i++){
            gtnr[i] = new double[Nnucp];
            nucrrnr[i] = new double[Nnucp];
        }
        for (int i = 0; i < Nnucp; i++){
            for (int j = 0; j < Nnucp; j++){
                gtnr[i][j] = 0.0;
                nucrrnr[i][j] = 0.0;
            }
        }

        nucp = new double [Nnucp];
        for (int i = 0; i < Nnucp; i++){
            nucp[i] = 0.0;
        }

        nucrr = new double [Nnucrr];
        for (int i= 0; i < Nnucrr; i++){
            nucrr[i] = 0.0;
        }


}

BranchSpecificParameters::~BranchSpecificParameters()
{
    //dtor
}

void BranchSpecificParameters::SetLocalParaemters(double* nucp, double* nucrr){

        for (int i = 0; i < Nnucp; i++){
            this->nucp[i] = nucp[i];
        }


        for (int i= 0; i < Nnucrr; i++){
            this->nucrr[i] = nucrr[i];
        }

}
