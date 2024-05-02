#ifndef MODEL__H
#define MODEL__H

#include <Rcpp.h>
#include <iostream>


using namespace Rcpp;

// STRUCTURES 

struct lambda_c {
    Rcpp::IntegerVector individual;
    Rcpp::NumericVector lambda;
};


struct lambda_e {
    Rcpp::IntegerVector individual;
    Rcpp::NumericVector lambda;
};



struct interaction {
    Rcpp::IntegerVector from;
    Rcpp::IntegerVector to;
};



struct localisation {
    Rcpp::IntegerVector individual;
    Rcpp::IntegerVector position;
};

struct environment {
    Rcpp::NumericVector env;
    Rcpp::IntegerVector room;
};

struct status {
    Rcpp::IntegerVector individual;
    Rcpp::IntegerVector status;
    //double length_infection; total db for infection then splice?
};
// delete? 
struct id {
    std::string id;
};



// Fonctions

interaction Update_interaction(
    Rcpp::DataFrame data,
    int ti
);

environment Update_environment(
    environment env_prev,
    const localisation loc_current,
    const status status_prev,
    const double mu,
    const double nu,
    const int dt
);

status Update_status(
    const localisation loc_current,
    const environment env_current,
    const status infected_prev,
    const double alpha,
    const double beta,
    const double epsilon,
    const int dt
);

Rcpp::List List_encountered(
    const int id,
    const interaction interactions
);




#endif