#ifndef MODEL__H
#define MODEL__H

#include <Rcpp.h>
#include <iostream>


using namespace Rcpp;

// STRUCTURES 

struct lambda_c {
    id individual;
    double lambda;
}

struct graph {
    graph graph;
};

struct id {
    std::string id;
};

struct interaction {
    id from;
    id to;
};

struct interaction_loc {
    id individual;
    int position;
};

struct localisation {
    id individual;
    int position;
};

struct environment {
    double env;
    int position;
};

struct status {
    id individual;
    int status;
    //double length_infection; total db for infection then splice?
};

struct visited

// Fonctions

interaction Update_interaction(
    graph graph,
    int subdivision
);

interaction_loc Associate_interaction(
    interaction int_current
);

localisation Update_localisation(
    localisation loc_prev,
    interaction_loc int_loc_current
);

environment Update_environment(
    environment env_prev,
    const localisation loc_current,
    const status status_prev,
    const double mu,
    const double nu,
    const integer dt
);

status Update_status(
    const localisation loc_current,
    const environment env_current,
    const infected infected_prev,
    const double alpha,
    const double beta,
    const double epsilon,
    const integer dt
);




#endif