#include "model-nodscov2_fun.h"
#include <Rcpp.h>
#include <iostream>

using namespace std;

// [[Rcpp::plugins(cpp11)]]
double mIncub = 1.63 * 24*60*2;
double sdIncub = 0.5 * 24*60*2;

double m_incub_g = 4.07 * 24*60*2;
double sd_incub_g = 2.12 * 24*60*2;
double shape_incub_g = pow(m_incub_g,2) / pow(sd_incub_g, 2);
double scale_incub_g = pow(sd_incub_g,2) / m_incub_g;




// R UNIQUE(X) FUNCTION
Rcpp::Environment base("package:base");
Function do_unique = base["unique"];
Function do_sample = base["sample"];

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Update_status(
    Rcpp::DataFrame status_tim1,
    Rcpp::DataFrame lambda_ti
) {
    Rcpp::NumericVector status_prev = status_tim1["status"];
    Rcpp::NumericVector lambda_c = lambda_ti["lambda_c"];
    Rcpp::NumericVector lambda_e = lambda_ti["lambda_c"];
    Rcpp::NumericVector status_ti = clone(status_prev);


    Rcpp::NumericVector FOI = (lambda_ti.nrows(), 1) - exp(- (lambda_c  + lambda_e));
    for (int j=0; j < lambda_ti.nrows(); j++){
        if (status_prev[j] == 0 && R::runif(0, 1) <= FOI[j]){
            status_ti[j] = 1;
        }
    };

    return status_ti;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::IntegerVector Get_status_t(
    Rcpp::DataFrame global_status,
    int t
) {
    Rcpp::IntegerVector t_inf = global_status["t_inf"];
    Rcpp::IntegerVector t_recover = global_status["t_recover"];
    Rcpp::IntegerVector status_t (t_inf.size()) ;
    for(int j = 0; j < t_inf.size(); j++){
        if (t_inf[j] != -1 && t >= t_inf[j] && t <= t_recover[j]){
            status_t[j] = 1;
        } else if (t_inf[j] != -1 && t > t_recover[j]){
            status_t[j] = 2;
        } else{
            status_t[j] = 0;
        }
    }

    return status_t;
};





//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame Update_status_bis(
    Rcpp::DataFrame global_status,
    Rcpp::DataFrame lambda_ti,
    Rcpp::DataFrame info_patient_HCW,
    Rcpp::DataFrame interactions_ti,
    Rcpp::DataFrame localization_ti,
    int t
) {
    Rcpp::NumericVector lambda_c = lambda_ti["lambda_c"];
    Rcpp::NumericVector lambda_e = lambda_ti["lambda_c"];
    
    Rcpp::IntegerVector status_tim1 = Get_status_t(global_status, t);

    Rcpp::CharacterVector ids = info_patient_HCW["id"];
    Rcpp::IntegerVector info_int = info_patient_HCW["info"];
    Rcpp::IntegerVector info_room = info_patient_HCW["room"];
    
    Rcpp::IntegerVector t_inf_tim1 = global_status["t_inf"];
    Rcpp::IntegerVector t_recover_tim1 = global_status["t_recover"];
    Rcpp::CharacterVector inf_by = global_status["inf_by"];
    Rcpp::IntegerVector inf_room_tim1 = global_status["inf_room"];
    
    Rcpp::DataFrame global_status_updated = clone(global_status);
    Rcpp::IntegerVector t_inf_ti = clone(t_inf_tim1);
    Rcpp::IntegerVector t_recover_ti = clone(t_recover_tim1);
    Rcpp::IntegerVector inf_room_ti = clone(inf_room_tim1);

    Rcpp::NumericVector FOI = (lambda_ti.nrows(), 1) - exp(- (lambda_c  + lambda_e));

    for (int j=0; j < lambda_ti.nrows(); j++){
        if (status_tim1[j] == 0 && R::runif(0, 1) <= FOI[j]){
            int room_j = -1;
            t_inf_ti[j] = t;
            t_recover_ti[j] = t + Incub_period_gamma();

            // // ROOM ?
            if (info_int[j] == 0){
                room_j= info_room[j];
            } else if (info_int[j] == 1){
                room_j = Get_loc_HCW(ids[j], info_patient_HCW, localization_ti);
            }
            inf_room_ti[j] = room_j;
            
            // Rcpp::CharacterVector ids = info_patient_HCW["id"];
            // Rcpp::List ind_inf = List_encountered(ids[j], interactions_ti);
            // Rcpp::CharacterVector sample();
            // Rcpp::String inf_by_j;

        }
    };
    global_status_updated["t_inf"] = t_inf_ti;
    global_status_updated["t_recover"] = t_recover_ti;
    global_status_updated["inf_room"] = inf_room_ti;

    return global_status_updated;
};
            

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame Get_t(
    Rcpp::List Global_list,
     int t
){
    // WARNING: c++ COUNTS STARTS AT 0
    Rcpp::DataFrame df(Global_list[t]);
    return df;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Update_environment(
    Rcpp::DataFrame environment_tim1,
    Rcpp::DataFrame localization_ti,
    Rcpp::IntegerVector status_tim1,
    Rcpp::DataFrame info_patient_HCW, //(id: id of the individual, info: "0" IF PATIENT, "1" IF HCW, room: room assigned to the individual, easier for patients...) 
    const double mu,
    const double nu,
    const double deltat
) {
    // THERE IS MULTIPLE WAYS TO ACHIEVE THIS
    // A. (naive) for loop on room ( for loop on individuals (check if patient/HCW and if the localization == room r then check if infected etc))
    // B. for loop on room (Patient assigned to this room is infected?  THEN for loop on localization(is the localization == room r ? and is the individual infected?))
    // C. FIRST loop on rooms/patients (is the patient infected? if so, update the env) THEN another for loop on localization (update the env for each row if the individual is infected)
    // C. FIRST loop on rooms/patients (is the patient infected? if so, update the env) THEN another loop on infected patients etc.
    // THIS FUNCTION USES THE C. METHOD
    // WARNING: IF THERE IS DOUBLE ROOMS, DONT DO THE DOUBLE EXPONENTIAL INACTIVATION OF E(t-1)
    
    Rcpp::IntegerVector rooms_environment = environment_tim1["id_room"];
    Rcpp::NumericVector temp_env = environment_tim1["env"];
    Rcpp::NumericVector env_ti = clone(temp_env);
    Rcpp::CharacterVector ids = info_patient_HCW["id"];
    Rcpp::IntegerVector info_patient_HCW_int = info_patient_HCW["info"];
    Rcpp::IntegerVector info_patient_HCW_room = info_patient_HCW["room"];
    Rcpp::IntegerVector localizations = localization_ti["localization"];
    Rcpp::CharacterVector id_HCW = localization_ti["id"];

    // EXPONENTIAL INACTIVATION 
    env_ti = env_ti * exp(-mu * deltat);
    
    // INFECTING INDIVIDUALS SHEDDING
    for (int j = 0; j < info_patient_HCW.nrows(); ++j) {
        // INFECTED PATIENTS SHEDDING IN THEIR ROOMS
        if (info_patient_HCW_int[j] == 0 && status_tim1[j] == 1){
            int room_j = info_patient_HCW_room[j];
            // Search for the index of the room in environment dataframe
            int index_room_j = -1;
            for (int k = 0; k < rooms_environment.size(); ++k) {
                if (rooms_environment[k] == room_j) {
                    index_room_j = k;
                    break;
                }
            }
            // UPDATE THE ENVIRONMENT FOR THE INFECTED PATIENT IN THAT ROOM
            if(index_room_j != -1){
                env_ti[index_room_j] += 1 * nu * deltat;}
        } 

        // INFECTED HCW SHEDDING IN THE ROOM HE WAS AT t-1
        if (info_patient_HCW_int[j] == 1 && status_tim1[j] == 1){
            // WARNING: LOCALIZATION INDEX != INDIVIDUAL INDEX ETC
            Rcpp::String id_j = ids[j];
            // Search for the room where was the HCW j at t-1
            int index_localization_j = -1;
            for (int k = 0; k < id_HCW.size(); k++){
                if (ids[j] == id_HCW[k]){
                   index_localization_j = k;
                   break; 
                }
            }
            int room_j = localizations[index_localization_j];
            // Search for the index associated to this room
            int index_room_j = -1;
            for (int k = 0; k < rooms_environment.size(); k++){
                if (rooms_environment[k] == room_j){
                   index_room_j = k;
                   break; 
                }
            }
            // UPDATE THE ENVIRONMENT FOR THE INFECTED HCW IN THAT ROOM
            if(index_room_j != -1){
                env_ti[index_room_j] += 1 * nu * deltat;
            }
        }
    }

  return env_ti;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List List_encountered(
    Rcpp::String id,
    Rcpp::DataFrame interaction_ti
) {
  Rcpp::List list_id;
  Rcpp::CharacterVector from = interaction_ti["from"];
  Rcpp::CharacterVector to = interaction_ti["to"];

  for (int j = 0; j < interaction_ti.nrows(); ++j) {
    if (from[j] == id) {
        Rcpp::String push = to[j];
        list_id.push_back(push);
    }
    if (to[j] == id) {
        Rcpp::String push = from[j];
        list_id.push_back(push);
    }
  }
  // Need to have unique list of ids
  list_id = do_unique(list_id);
  return list_id;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Lambda_c (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame interaction_ti,
    Rcpp::IntegerVector status_ti,
    const double beta,
    const double deltat
) {
    Rcpp::CharacterVector ids = lambda_tim1["id"];
    Rcpp::NumericVector temp_lambda_c = lambda_tim1["lambda_c"];
    Rcpp::NumericVector lambda_c_ti = clone(temp_lambda_c);
    Rcpp::String id;
    Rcpp::List list_ind_r;
    int nb_inf_r; 

    for (int j = 0; j < lambda_tim1.nrows(); ++j){
        id = ids[j];
        nb_inf_r = 0;
        list_ind_r = List_encountered(id, interaction_ti); // 
        for (int i = 0; i < list_ind_r.size(); ++i){
            Rcpp::String id_r = list_ind_r[i];
            // Search for the index of individual encountered in ids vector 
            int index_r = -1;
            for (int k = 0; k < ids.size(); ++k) {
                if (ids[k] == id_r) {
                    index_r = k;
                    break;
                }
            }
            // if individual r is infected & we found its index (for safety)
            if (index_r != -1 && status_ti[index_r] == 1) {
                nb_inf_r += 1;
            }
        }
        
        lambda_c_ti[j] = beta * deltat * nb_inf_r;
    }

    return lambda_c_ti;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Lambda_e (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame localization_ti,
    Rcpp::DataFrame environment_ti,
    Rcpp::DataFrame info_patient_HCW, 
    const double epsilon,
    const double deltat
) {
    Rcpp::CharacterVector ids_lambda = lambda_tim1["id"];
    Rcpp::NumericVector temp_lambda_e = lambda_tim1["lambda_e"];
    Rcpp::NumericVector lambda_e_ti = clone(temp_lambda_e);
    Rcpp::CharacterVector ids_localization = localization_ti["id"];
    Rcpp::IntegerVector info_patient_HCW_int = info_patient_HCW["info"]; // "0" IF PATIENT, "1" IF HCW
    Rcpp::IntegerVector info_patient_HCW_room = info_patient_HCW["room"];
    Rcpp::NumericVector environment = environment_ti["env"];
    Rcpp::IntegerVector rooms_environment = environment_ti["id_room"];
    Rcpp::IntegerVector localizations = localization_ti["localization"];

  for (int j = 0; j < lambda_tim1.nrows(); ++j){
    
    // TWO CASES (PATIENTS AND HCWS)
    // CASE 1. IF INDIVIDUAL j IS A PATIENT --> ENVIRONMENT ACCORDING TO ITS ROOM
    if (info_patient_HCW_int[j] == 0){
        int room_j = info_patient_HCW_room[j];
        // Search for the index of patient's room
        int index_room = -1;
            for (int k = 0; k < rooms_environment.size(); ++k) {
                if (rooms_environment[k] == room_j) {
                    index_room = k;
                    break;
                }
            }
        lambda_e_ti[j] = epsilon * deltat * environment[index_room];
    }
    
    // CASE 2. IF INDIVIDUAL j IS A HCW --> ENVIRONMENT ACCORDING TO ITS LOCALIZATION
    // WARNING, LOCALIZATION DF INDEX != LAMBDA DF INDEX ETC
    if (info_patient_HCW_int[j] == 1){
        // Search for the room where the HCW is located
        int index_localization = -1;
            for (int k = 0; k < ids_localization.size(); ++k) {
                if (ids_localization[k] == ids_lambda[j]) {
                    index_localization = k;
                    break;
                }
            }
        int room_j = localizations[index_localization];
        // Search for the index of this room
        int index_room = -1;
            for (int k = 0; k < rooms_environment.size(); ++k) {
                if (rooms_environment[k] == room_j) {
                    index_room = k;
                    break;
                }
            }
        lambda_e_ti[j] = epsilon * deltat * environment[index_room];
    }
  }

  return lambda_e_ti;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Incub_period_gamma() {
    // Incubation period -> Gamma distribution (shape,scale)
    double incubation_period_seconds = R::rgamma(shape_incub_g, scale_incub_g);
    int incubation_period_subdivisions = static_cast<int>(incubation_period_seconds / 30.0);
    
    return incubation_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Incub_period_lognormal() {
    // Incubation period -> Log-normal distribution (meanlog, sdlog)
    double incubation_period_seconds = R::rlnorm(mIncub, sd_incub_g);
    int incubation_period_subdivisions = static_cast<int>(incubation_period_seconds / 30.0);
    
    return incubation_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Get_loc_HCW(
    Rcpp::String id_HCW,
    Rcpp::DataFrame info_patient_HCW,
    Rcpp::DataFrame localization_ti
) {
    int localization;
    Rcpp::CharacterVector ids = info_patient_HCW["id"];
    Rcpp::CharacterVector ids_HCW = localization_ti["id"];
    Rcpp::IntegerVector localizations_HCW = localization_ti["localization"];
    Rcpp::IntegerVector status = info_patient_HCW["info"];
    
    int index_ind = -1;
    for (int k = 0; k < ids.size(); k++){
        if (id_HCW == ids[k]){
                index_ind = k;
                break;
                }
        }
    if (status[index_ind] == 1){
        int index_localization_j = -1;
        for (int k = 0; k < ids_HCW.size(); k++){
            if (ids[index_ind] == ids_HCW[k]){
            index_localization_j = k;
            break;
            }
        }
        int room_j = localizations_HCW[index_localization_j];
        if (room_j == -1){
            localization = -1;
        } else {
            localization = room_j;
            }
        }
    return localization;
};