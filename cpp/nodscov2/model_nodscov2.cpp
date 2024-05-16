#include "model-nodscov2_fun.h"
#include <Rcpp.h>


// [[Rcpp::plugins(cpp11)]]

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame simulation(
    Rcpp::List global_interaction,
    Rcpp::List global_localization,
    Rcpp::List global_environment,
    Rcpp::List global_status,
    Rcpp::List global_lambda,
    Rcpp::DataFrame info_patient_HCW,
    double beta,
    double epsilon,
    double nu,
    double mu,
    double dt
    
) {
    // INITIALIZATION
    double tau = 60 * 60 * 24; //seconds in a day
    double deltat = dt/(tau);

    Rcpp::DataFrame environment_ti;
    Rcpp::DataFrame environment_tim1;
    
    Rcpp::DataFrame status_ti;
    Rcpp::DataFrame status_tim1;

    Rcpp::DataFrame interaction_ti;
    Rcpp::DataFrame interaction_tim1;
    
    Rcpp::DataFrame lambda_ti;
    Rcpp::DataFrame lambda_template = Get_t(global_lambda, 0);

    Rcpp::DataFrame localization_ti;
    Rcpp::DataFrame localization_tim1;
    

    ////////////////
    // SIMULATION //
    ////////////////
    int sim_size = global_interaction.size();
    for (int t = 1; t < sim_size; t++){
        interaction_ti = Get_t(global_interaction, t);
        localization_tim1 = Get_t(global_localization, t-1);
        status_tim1 = Get_t(global_status, t-1);
        ////////////////////////////
        // Update the environment //
        ////////////////////////////
        environment_tim1 = Get_t(global_environment, t-1);
        environment_ti = clone(environment_tim1);
        environment_ti["env"] = Update_environment(environment_ti, localization_tim1, status_tim1, info_patient_HCW, mu, nu, deltat);
        global_environment[t] = environment_ti;

        ///////////////////
        // Update Lambda //
        ///////////////////
        lambda_ti = clone(lambda_template);
        lambda_ti["lambda_e"] = Lambda_e(lambda_template, localization_ti, environment_ti, info_patient_HCW, epsilon, deltat);
        lambda_ti["lambda_c"] = Lambda_c(lambda_template, interaction_ti, status_tim1, beta, deltat);
        
        
        ///////////////////////
        // Update the status //
        ///////////////////////
        status_tim1 = Get_t(global_status, t-1);
        status_ti = clone(status_tim1);
        status_ti["status"] = Update_status(lambda_ti, status_tim1);
        global_status[t] = status_ti;




    }


    return(info_patient_HCW);
}

    // Rcpp::List clusters,
    // Rcpp::DataFrame inferred_admission,
    // Rcpp::DataFrame HCW_interacting_id,
    // Rcpp::List double_rooms,
    // Rcpp::Datetime begin_date,
    // Rcpp::Datetime end_date,
    // int n_subdivisions,