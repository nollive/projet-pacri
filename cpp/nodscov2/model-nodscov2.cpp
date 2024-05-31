#include "model-nodscov2_fun.h"
#include <Rcpp.h>


#include <iostream>

// [[Rcpp::plugins(cpp11)]]

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List simulation(
    Rcpp::List global_interaction,
    Rcpp::List global_localization,
    Rcpp::List global_environment,
    Rcpp::List global_lambda,
    Rcpp::DataFrame global_status,
    Rcpp::DataFrame admission,
    double beta,
    double epsilon,
    double nu,
    double mu,
    double env_thresold,
    double dt
    
) {
    // INITIALIZATION
    double tau = 60 * 60 * 24; //seconds in a day
    double deltat = dt/(tau);

    Rcpp::DataFrame environment_ti = Get_t(global_environment, 0);
    Rcpp::DataFrame environment_tim1 = Get_t(global_environment, 0);
    
    Rcpp::IntegerVector status_ti =  Get_status_t(global_status, 0);
    Rcpp::IntegerVector status_tim1;

    Rcpp::DataFrame interaction_ti;
    Rcpp::DataFrame interaction_tim1;
    
    Rcpp::DataFrame lambda_ti;
    Rcpp::DataFrame lambda_template;
    lambda_template = Get_t(global_lambda, 0);
    
    Rcpp::DataFrame localization_ti = Get_t(global_localization, 0);
    Rcpp::DataFrame localization_tim1;
    
    ///////////////
    // R's t = 1 //
    ///////////////
    // Shedding of the index patient
    environment_ti["env"] = Update_environment(environment_tim1, localization_ti , status_ti, admission, mu, nu, deltat);
    global_environment[0] = environment_ti;
    // update status for time t = 1?
    
    //////////////////////////
    // SIMULATION (R's t>1) //
    //////////////////////////
    int sim_size = global_interaction.size();
    for (int t = 1; t < sim_size; t++){
        interaction_ti = Get_t(global_interaction, t);
        localization_tim1 = Get_t(global_localization, t-1);
        status_tim1 = Get_status_t(global_status, t-1);
        ////////////////////////////
        // Update the environment //
        ////////////////////////////
        environment_tim1 = Get_t(global_environment, t-1);
        environment_ti = clone(environment_tim1);
        localization_ti = Get_t(global_localization, t);
        environment_ti["env"] = Update_environment(environment_ti, localization_ti, status_tim1, admission, mu, nu, deltat);
        global_environment[t] = environment_ti;

        ///////////////////
        // Update Lambda //
        ///////////////////
        lambda_ti = clone(lambda_template);
        lambda_ti["lambda_e"] = Lambda_e(lambda_template, localization_ti, environment_ti, admission, epsilon, env_thresold, deltat);
        lambda_ti["lambda_c"] = Lambda_c(lambda_template, interaction_ti, status_tim1, beta, deltat);
        global_lambda[t] = lambda_ti;
        
        ///////////////////////
        // Update the status //
        ///////////////////////
        Rcpp::DataFrame temp = Update_status_bis(global_status, lambda_ti,admission, interaction_ti, localization_ti, t);
        global_status = clone(temp);
    }

    Rcpp::List res = Rcpp::List::create(_["global_status"] = global_status, _["global_lambda"] = global_lambda, _["global_environment"] = global_environment);
    return res;    
}