#ifndef MODEL_NODSCOV2_FUN__H
#define MODEL_NODSCOV2_FUN__H

#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// extern double beta;
// extern double epsilon;
// extern double nu;
// extern double mu;
// extern double tau;
extern double deltat;
extern double mIncub;
extern double sdIncub;
extern double maxPCRDetectability;
extern double m_incub_g;
extern double sd_incub_g;
extern double shape_incub_g;
extern double scale_incub_g;

extern double mInf;
extern double sdInf;
extern double m_inf_g;
extern double sd_inf_g;
extern double shape_inf_g;
extern double scale_inf_g;


RcppExport int Incub_period_gamma();

RcppExport int Incub_period_lognormal();

RcppExport int Inf_period_gamma();

RcppExport int Inf_period_lognormal();

RcppExport int Incub_period_uniform();

RcppExport int Inf_period_uniform();

    

RcppExport Rcpp::NumericVector Update_environment(
    Rcpp::DataFrame environment_tim1,
    Rcpp::DataFrame localization_ti,
    Rcpp::IntegerVector status_tim1,
    Rcpp::DataFrame admission, //(id: id of the individual, info: "0" IF PATIENT, "1" IF HCW, room: room assigned to the individual, easier for patients...) 
    const double mu,
    const double nu,
    const double deltat
);


RcppExport Rcpp::List List_encountered(
    Rcpp::String id,
    Rcpp::DataFrame interaction_ti
);

RcppExport Rcpp::List List_inf_encountered(
    Rcpp::String id,
    Rcpp::DataFrame interaction_ti,
    Rcpp::DataFrame global_status,
    Rcpp::DataFrame admission,
    int t
);

RcppExport Rcpp::NumericVector Lambda_c (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame interaction_ti,
    Rcpp::IntegerVector status_ti,
    const double beta,
    const double deltat
);

RcppExport Rcpp::NumericVector Lambda_e (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame localization_ti,
    Rcpp::DataFrame environment_ti,
    Rcpp::DataFrame admission, 
    const double epsilon,
    const double env_thresold,
    const double deltat
);

RcppExport Rcpp::DataFrame Get_t(
    Rcpp::List Global_list,
     int t
);

RcppExport Rcpp::IntegerVector Get_status_t(
    Rcpp::DataFrame global_status,
    int t
);

RcppExport int Get_status_j(
    Rcpp::String id,
    Rcpp::DataFrame global_status,
    Rcpp::DataFrame admission,
    int t
);

RcppExport Rcpp::String Sample_inf(
    Rcpp::String id,
    Rcpp::List list_inf_encountered,
    Rcpp::DataFrame admission,
    Rcpp::DataFrame localization_ti,
    double lambda_e_j,
    double lambda_c_j
);


RcppExport Rcpp::DataFrame Update_status_bis(
    Rcpp::DataFrame global_status,
    Rcpp::DataFrame lambda_ti,
    Rcpp::DataFrame admission,
    Rcpp::DataFrame interactions_ti,
    Rcpp::DataFrame localization_ti,
    int t
);

RcppExport int Get_loc_HCW(
    Rcpp::String id_HCW,
    Rcpp::DataFrame admission,
    Rcpp::DataFrame localization_ti
);

RcppExport int Get_loc_j(
    Rcpp::String id,
    Rcpp::DataFrame admission,
    Rcpp::DataFrame localization_ti
); 
#endif