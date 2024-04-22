#ifndef MODEL__H
#define MODEL__H

#include <Rcpp.h>
#include <iostream>


using namespace Rcpp;

// STRUCTURES 

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
};

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
    localisation loc_current,
    status status_prev,
    double mu,
    double nu
);

status Update_status(
    localisation loc_current,
    environment env_current,
    infected infected_prev,
    double alpha,
    double beta,
    double epsilon,
    integer dt
);




#endif