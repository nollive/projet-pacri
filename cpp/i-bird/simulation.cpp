#include <Rcpp.h>
#include "model.h"
#include "simulation.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]



//////////////////////////////////////////////
// [[Rcpp::export]]

int main()
{
    Rcpp::DataFrame data;

    ///////////
    // CONST //
    ///////////
    const int dt = 30;
    const Rcpp::Datetime begin_date("2009-07-06 00:00:00");
    const Rcpp::Datetime end_date("2009-09-28 00:00:00");
    const int nb_ti = (end_date - begin_date) / 30;

    double mu;
    double nu;
    double epsilon;
    double beta;

    ////////////////////
    // INITIALIZATION //
    ////////////////////
    environment environment_ti;
    localisation localisation_ti;
    localisation localisation_tim1;
    status status_ti;
    status status_tim1;
    // lambda_c lambda_c;
    // lambda_e lambda_e;

    ////////////////
    // SIMULATION //
    ////////////////
    for (int ti = 0; ti < nb_ti; ti++){
        // GET INT //
        interaction interaction_ti = Update_interaction(data, ti, "ti", begin_date, dt);
        interaction interaction_tim1 = Update_interaction(data, ti, "tim1", begin_date, dt);
        // GET LOC //
        localisation localisation_ti ;
        localisation localisation_tim1 ;
        // UPDATE ENV //
        environment_ti = Update_environment(environment_ti, localisation_tim1, status_tim1, mu, nu, dt );
        // UPDATE STATUS //
        // I.E. UPDATE STATUS & RANDOM DRAW OF INFECTION/INFECTION TIME //
        status_ti = Update_status(localisation_ti, environment_ti, status_tim1, beta, epsilon, nu, mu, dt);
        // OPTIONNAL: STORE ti DATA IN GLOBAL DATA //
    };

}