#include <Rcpp.h>
#include "model.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]



//////////////////////////////////////////////
// [[Rcpp::export]]

// OK
environment Update_environment_bis(
    environment env_prev,
    const localisation loc_prev,
    const status status_prev,
    const double mu,
    const double nu,
    const integer dt
){

        for (const auto& room : env_prev) {
            int n_inf = 0;
            for (const auto& ind : loc_prev){
                if (ind.location == room.position){
                    int identifiant = ind.id;
                    if (status_prev[identifiant].status == 1){
                        n_inf +=1;
                    }
                }
            }
            room.env = room.env * Rcpp:exp(- mu * dt) + n_inf * nu * dt;
        }

    return env_prev; 
}

lambda_c Lambda_c (
    lambda_c lambda_c_prev,
    const double beta,
    const double deltat,
    const interaction interaction_ti,
    const status status_ti,
    const int ti,
){
    for (const auto& ind : lambda_c_prev.individual){
        int nb_inf_r = 0;
        list liste_ind_r = Rencontre(ind, ti);
        for (const auto& r : liste_ind_r){
            if (status_ti[r].status == 1){
                nb_inf_r += 1;
            }
        }
        ind.lambda = beta * dt * nb_inf_r;
    }
    return lambda_c_prev;
};




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








status Update_status(
    localisation loc_current,
    environment env_current,
    infected infected_prev,
    double alpha,
    double beta,
    double epsilon)
    {



    }


// environment Update_environment(
//     environment env_prev,
//     localisation loc_current,
//     status status_prev,
//     double mu,
//     double nu
// ){
//     for (int i = 1; t < sizeof(status_prev); ++t) {
//         if (status_prev[i].status == 1) {


//         }

//     }


// }