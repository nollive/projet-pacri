#include <Rcpp.h>
#include "model.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]



//////////////////////////////////////////////
// [[Rcpp::export]]

interaction Update_interaction(graph graph, int subdivision){
    // Filtration du graphe
    return int_current;
};


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
){
    for (int r = 1; t < sizeof(env_prev); ++t) {
        if (r ==)

    }


}

status Update_status(
    localisation loc_current,
    environment env_current,
    infected infected_prev,
    double alpha,
    double beta,
    double epsilon)
    {



    }
