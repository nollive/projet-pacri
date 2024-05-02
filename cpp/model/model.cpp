#include <Rcpp.h>
#include "model.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]



//////////////////////////////////////////////
// [[Rcpp::export]]

interaction Update_interaction(
    Rcpp::DataFrame data,
    int ti
) {
    interaction int_tim1;
    interaction int_ti;


    return int_ti;
};


// à voir sià changer (localisation avec prev et future)
environment Update_environment_ter(
    environment environment_tim1,
    const localisation& localisation_tim1,
    const status& status_tim1,
    const double mu,
    const double nu,
    const int dt
) {
    for (int r = 0; r < environment_tim1.room.size(); ++r) {
        int n_inf = 0;
        for (int j = 0; j < localisation_tim1.individual.size(); ++j) {
            if (localisation_tim1.position[j] == environment_tim1.room[r]) { // if ind j is in room r (at time ti-1)
                int identifiant = localisation_tim1.individual[j];
                if (status_tim1.status[identifiant] == 1) { // if ind j is infected (at time ti-1)
                    n_inf += 1;
                }
            }
        }
        environment_tim1.env[r] = environment_tim1.env[r] * exp(-mu * dt) + n_inf * nu * dt;
    }

    return environment_tim1;
};


Rcpp::List List_encountered(
    const int id,
    const interaction interactions
) {
    Rcpp::List list_id;
    for (int j = 0; j < interactions.from.size(); ++j) {
        if (interactions.from[j] == id) {
            list_id.push_back(interactions.to[j]);
        }
        if (interactions.to[j] == id) {
            list_id.push_back(interactions.from[j]);
        }
    }

    return list_id;
};


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
) {
    for (int j = 0; j < lambda_c_tim1.individual.size(); ++j){
        int identifiant = lambda_c_tim1.individual[j];
        int nb_inf_r = 0;
        Rcpp::List liste_ind_r = List_encountered(identifiant, interaction_ti); // 
        for (int i = 0; i < liste_ind_r.size(); ++i){
            int identifiant_r = liste_ind_r[i];
            if (status_ti.status[identifiant_r] == 1){
                nb_inf_r += 1;
            }
        }
        lambda_c_tim1.lambda[j] = beta * dt * nb_inf_r;
    }
    return lambda_c_tim1;
};

lambda_e Lambda_e (
    lambda_e lambda_e_prev,
    const localisation localisation_ti,
    const environment environment_ti,
    const double epsilon,
    const double dt
) {
    for (int j = 0; j < lambda_e_prev.individual.size(); ++j){
        int identifiant = lambda_e_prev.individual[j];
        int room  = localisation_ti.position[identifiant];
        int env = environment_ti.env[room];

        lambda_e_prev.lambda[j] = epsilon * dt * env;
    }
    return lambda_e_prev;
};

status Update_status(
    const localisation localisation_ti,
    const environment environment_ti,
    const status infected_tim1,
    const double alpha,
    const double beta,
    const double epsilon,
    const int dt
);


//////////////////// GARBAGE //////////////////////////






// interaction Update_interaction(graph graph, int subdivision){
//     // Filtration du graphe
//     return int_current;
// };


// interaction_loc Associate_interaction(
//     interaction int_current
// );

// localisation Update_localisation(
//     localisation loc_prev,
//     interaction_loc int_loc_current
// );



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

// environment Update_environment_bis(
//     environment env_prev,
//     const localisation loc_prev,
//     const status status_prev,
//     const double mu,
//     const double nu,
//     const int dt
// ){

//         for (const auto& room : env_prev.room) {
//             int n_inf = 0;
//             for (const auto& ind : loc_prev.individual){
//                 if (ind.position == room.room){
//                     int identifiant = ind.id;
//                     if (status_prev[identifiant].status == 1){
//                         n_inf +=1;
//                     }
//                 }
//             }
//             room.env = room.env * exp(- mu * dt) + n_inf * nu * dt;
//         }

//     return env_prev; 
// }