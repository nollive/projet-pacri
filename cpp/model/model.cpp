#include "model.h"
#include <Rcpp.h>

using namespace std;

// [[Rcpp::plugins(cpp11)]]

//////////////////////////////////////////////
// [[Rcpp::export]]

environment Update_environment(
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
}


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
}


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
}

lambda_e Lambda_e (
    lambda_e lambda_e_tim1,
    const localisation localisation_ti,
    const environment environment_ti,
    const double epsilon,
    const double dt
) {
    for (int j = 0; j < lambda_e_tim1.individual.size(); ++j){
        int identifiant = lambda_e_tim1.individual[j];
        int room  = localisation_ti.position[identifiant];
        int env = environment_ti.env[room];

        lambda_e_tim1.lambda[j] = epsilon * dt * env;
    }
    return lambda_e_tim1;
}

status Update_status(
    const localisation localisation_ti,
    const environment environment_ti,
    const status status_tim1,
    const double beta,
    const double epsilon,
    const double nu,
    const double mu,
    const int dt
){
    return status_tim1;
}


interaction Update_interaction(
    Rcpp::DataFrame data,
    int ti,
    Rcpp::String time,
    const Rcpp::Datetime begin_date, // can change (apply offset of (begin_date) in data --> origin is 00-00-00 00:00:00)
    const int dt
) {
    
    if (time == "tim1"){
        interaction int_tim1;
        int_tim1.from = 1;
        int_tim1.to = 1;

        return int_tim1;
        }

    if (time == "ti"){
        interaction int_ti;


        return int_ti;
    }
    interaction int_ti;
    return int_ti; // if time is not ti or tim1
}