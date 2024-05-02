#ifndef MODEL__H
#define MODEL__H

#include <Rcpp.h>
#include <iostream>


/////////////
// STRUCTS //  
/////////////

typedef struct {
    Rcpp::IntegerVector individual;
    Rcpp::NumericVector lambda;
} lambda_c;

typedef struct{
    Rcpp::IntegerVector individual;
    Rcpp::NumericVector lambda;
} lambda_e ;

typedef struct {
    Rcpp::IntegerVector from;
    Rcpp::IntegerVector to;
}  interaction ;

typedef struct {
    Rcpp::IntegerVector individual;
    Rcpp::IntegerVector position;
} localisation ;

typedef struct {
    Rcpp::NumericVector env;
    Rcpp::IntegerVector room;
} environment ;

typedef struct {
    Rcpp::IntegerVector individual;
    Rcpp::IntegerVector status;
} status ;


///////////
// Utils // 
///////////

Rcpp::List List_encountered(
    const int id,
    const interaction interactions
);



////////////
// UPDATE // 
////////////

interaction Update_interaction(
    Rcpp::DataFrame data,
    int ti,
    Rcpp::String time,
    const Rcpp::Datetime begin_date, // can change (apply offset of (begin_date) in data --> origin is 00-00-00 00:00:00)
    const int dt
);

environment Update_environment(
    environment environment_tim1,
    const localisation& localisation_tim1,
    const status& status_tim1,
    const double mu,
    const double nu,
    const int dt
);

status Update_status(
    const localisation localisation_ti,
    const environment environment_ti,
    const status status_tim1,
    const double beta,
    const double epsilon,
    const double nu,
    const double mu,
    const int dt
);



/////////
// FOI //
/////////

lambda_c Lambda_c (
    lambda_c lambda_c_tim1,
    const interaction interaction_ti,
    const status status_ti,
    const double beta,
    const double dt,
    const int ti
);

lambda_e Lambda_e (
    lambda_e lambda_e_tim1,
    const localisation localisation_ti,
    const environment environment_ti,
    const double epsilon,
    const double dt
);


#endif