#ifndef MODEL_NODSCOV2_FUN__H
#define MODEL_NODSCOV2_FUN__H

#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

extern double beta;
extern double epsilon;
extern double nu;
extern double mu;
extern double tau;


// Update ENV
RcppExport DataFrame Update_environment(
    DataFrame environment_tim1,
    DataFrame interaction_tim1,
    DataFrame status_tim1,
    double mu,
    double nu,
    double dt
    //double t
);



// List encountered

// FOI (lambda_c / lambda_e)

// Update STATUS

// UPDATE t --> GIVEN GLOBAL LISTS -> GIVE DATAFRAME FOR TIME t
#endif