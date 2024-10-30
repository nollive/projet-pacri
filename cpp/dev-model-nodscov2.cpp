#include "dev-model-nodscov2_fun.h"
#include <Rcpp.h>


#include <iostream>

// [[Rcpp::plugins(cpp11)]]

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List simulation(
    Rcpp::List global_interaction,
    Rcpp::List global_environment,
    Rcpp::List global_data,
    Rcpp::DataFrame global_status,
    double beta_c,
    double beta_e,
    double B,
    double nu,
    double mu,
    String env_model,
    double dt
    
) {
    // INITIALIZATION
    double tau = 60 * 60 * 24; //seconds in a day
    double deltat = dt/(tau);

    Rcpp::DataFrame environment_ti = Get_t(global_environment, 0);
    Rcpp::DataFrame environment_tim1 = Get_t(global_environment, 0);
    
    Rcpp::IntegerVector status_ti;
    Rcpp::IntegerVector status_tim1;

    Rcpp::DataFrame interaction_ti;
    Rcpp::DataFrame interaction_tim1;

    Rcpp::IntegerVector location_ti;
    Rcpp::IntegerVector location_tmi;

    Rcpp::DataFrame global_data_t;
    Rcpp::CharacterVector ids_ti;
    Rcpp::IntegerVector info_ti;
    

    ///////////////
    // R's t = 1 //
    ///////////////
    ids_ti = Get_t(global_data, 0)["id"];
    info_ti = Get_t(global_data, 0)["info"];
    location_ti = Get_t(global_data, 0)["location_ti"];
    status_tim1 = Get_status_t(global_status, ids_ti, 0); 
    // Shedding of the index patient
    environment_ti["env"] = Update_environment(ids_ti, info_ti, environment_tim1, location_ti, status_tim1, mu, nu, deltat, 0);
    global_environment[0] = environment_ti;
    // update status for time t = 1?
    
    //////////////////////////
    // SIMULATION (R's t>1) //
    //////////////////////////
    int sim_size = global_interaction.size();
    for (int t = 1; t < sim_size; t++){
        Rcpp::DataFrame global_data_t = Get_t(global_data, t);
        Rcpp::CharacterVector ids_ti = global_data_t["id"];
        Rcpp::IntegerVector info_ti = global_data_t["info"];
        Rcpp::IntegerVector location_ti = global_data_t["location_ti"];

        interaction_ti = Get_t(global_interaction, t);
        status_tim1 = Get_status_t(global_status, ids_ti, t-1);
        status_ti = Get_status_t(global_status, ids_ti, t);
        ////////////////////////////
        // Update the environment //
        ////////////////////////////
        environment_tim1 = Get_t(global_environment, t-1);
        environment_ti = clone(environment_tim1);
        environment_ti["env"] = Update_environment(ids_ti, info_ti, environment_tim1, location_ti, status_tim1, mu, nu, deltat, 0);
        global_environment[t] = environment_ti;

        ///////////////////
        // Update Lambda //
        ///////////////////
        Rcpp::DataFrame temp_global_data = clone(global_data_t);
        temp_global_data["lambda_e"] = Lambda_e(info_ti, location_ti, environment_ti, beta_e, B, env_model, deltat);
        temp_global_data["lambda_c"] = Lambda_c(ids_ti, interaction_ti, status_ti, beta_c, deltat);
        global_data[t] = temp_global_data;
        ///////////////////////
        // Update the status //
        ///////////////////////
        Rcpp::DataFrame temp = Update_status_bis(global_status, status_tim1, temp_global_data["lambda_c"], temp_global_data["lambda_e"], info_ti, ids_ti, interaction_ti, location_ti, t);
        global_status = clone(temp);
    }

    Rcpp::List res = Rcpp::List::create(_["global_status"] = global_status, _["global_data"] = global_data, _["global_environment"] = global_environment);
    return res;    
}