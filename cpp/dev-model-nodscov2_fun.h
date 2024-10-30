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
extern const int z;
extern const double deltat;
extern const double mIncub;
extern const double sdIncub;
extern const double maxPCRDetectability;
extern const double m_incub_g;
extern const double sd_incub_g;
extern const double shape_incub_g;
extern const double scale_incub_g;

extern const double mInf;
extern const double sdInf;
extern const double m_inf_g;
extern const double sd_inf_g;
extern const double shape_inf_g;
extern const double scale_inf_g;


RcppExport int Incub_period_gamma();

RcppExport int Incub_period_lognormal();

RcppExport int Inf_period_gamma();

RcppExport int Inf_period_lognormal();

RcppExport int Incub_period_uniform();

RcppExport int Inf_period_uniform();

    

RcppExport Rcpp::NumericVector Update_environment(
    const Rcpp::CharacterVector& ids_ti,
    const Rcpp::IntegerVector& info_ti,
    const Rcpp::DataFrame& environment_tim1,
    const Rcpp::IntegerVector& location_ti,
    const Rcpp::IntegerVector& status_tim1,
    const double& mu,
    const double& nu,
    const double& deltat,
    const int& t
);


RcppExport Rcpp::List List_encountered(
    const Rcpp::String& id,
    const Rcpp::DataFrame& interaction_ti
);

RcppExport Rcpp::List List_inf_encountered(
    const Rcpp::String& id,
    const Rcpp::DataFrame& interaction_ti,
    const Rcpp::DataFrame& global_status,
    const int& t
);

RcppExport Rcpp::NumericVector Lambda_c (
    const Rcpp::CharacterVector& ids_ti, 
    const Rcpp::DataFrame& interaction_ti,
    const Rcpp::IntegerVector& status_ti,
    const double& beta_c,
    const double& deltat
);

RcppExport Rcpp::NumericVector Lambda_e (
    const Rcpp::IntegerVector& info_ti,
    const Rcpp::IntegerVector& location_ti,
    const Rcpp::DataFrame& environment_ti,
    const double& beta_e,
    const double& B,
    const String& env_model,
    const double& deltat
);

RcppExport Rcpp::DataFrame Get_t(
    const Rcpp::List& Global_list,
    const int& t
);

RcppExport Rcpp::IntegerVector Get_status_t(
    const Rcpp::DataFrame& global_status,
    const Rcpp::StringVector& ids,
    const int& t
);

RcppExport int Get_status_j(
    const Rcpp::String& id,
    const Rcpp::DataFrame& global_status,
    const int& t
);

RcppExport Rcpp::String Sample_inf(
    const Rcpp::String& id,
    const Rcpp::List& list_inf_encountered,
    const Rcpp::CharacterVector& ids_ti,
    const Rcpp::IntegerVector& location_ti,
    const double& lambda_e_j,
    const double& lambda_c_j
);

RcppExport Rcpp::DataFrame Update_status_bis(
    const Rcpp::DataFrame& global_status,
    const Rcpp::IntegerVector& status_tim1,
    const Rcpp::NumericVector& lambda_c_ti,
    const Rcpp::NumericVector& lambda_e_ti,
    const Rcpp::IntegerVector& info_ti,
    const Rcpp::CharacterVector& ids_ti,
    const Rcpp::DataFrame& interactions_ti,
    const Rcpp::IntegerVector& location_ti,
    const int& t
);


RcppExport int Get_loc_j(
    const Rcpp::String& id,
    const Rcpp::IntegerVector& location_ti,
    const Rcpp::CharacterVector& ids_ti
); 
#endif