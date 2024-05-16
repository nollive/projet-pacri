#ifndef MODEL_NODSCOV2_FUN__H
#define MODEL_NODSCOV2_FUN__H

#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

extern double beta;
extern double epsilon;
// extern double nu;
// extern double mu;
// extern double tau;
extern double deltat;



RcppExport Rcpp::DataFrame Get_t(
    Rcpp::List Global_list,
     int t
);


RcppExport Rcpp::NumericVector Update_status(
    Rcpp::DataFrame status_tim1,
    Rcpp::DataFrame lambda_ti
);


RcppExport Rcpp::NumericVector Update_environment(
    Rcpp::DataFrame environment_tim1,
    Rcpp::DataFrame localization_tim1,
    Rcpp::DataFrame status_tim1,
    Rcpp::DataFrame info_patient_HCW, //(id: id of the individual, info: "0" IF PATIENT, "1" IF HCW, room: room assigned to the individual, easier for patients...) 
    const double mu,
    const double nu,
    const int deltat
);


RcppExport Rcpp::List List_encountered(
    Rcpp::String id,
    Rcpp::DataFrame interaction_ti
);

RcppExport Rcpp::NumericVector Lambda_c (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame interaction_ti,
    Rcpp::DataFrame status_ti,
    const double beta,
    const double deltat
);

RcppExport Rcpp::NumericVector Lambda_e (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame localization_ti,
    Rcpp::DataFrame environment_ti,
    Rcpp::DataFrame info_patient_HCW, 
    const double epsilon,
    const double deltat
);

// Update ENV
// RcppExport DataFrame Update_environment(
//     DataFrame environment_tim1,
//     DataFrame interaction_tim1,
//     DataFrame status_tim1,
//     //DataFrame rooms,
//     double mu,
//     double nu,
//     double deltat
//     //double t
// );



// List encountered

// FOI (lambda_c / lambda_e)

// Update STATUS

// UPDATE t --> GIVEN GLOBAL LISTS -> GIVE DATAFRAME FOR TIME t
#endif