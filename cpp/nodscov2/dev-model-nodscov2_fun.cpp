#include "model-nodscov2_fun.h"
#include <Rcpp.h>
#include <iostream>

using namespace std;

// [[Rcpp::plugins(cpp11)]]
const int z = 10;
const double mIncub = (1.63) * (24 * 60 * 2) * z;
const double sdIncub = (0.5) * (24 * 60 * 2) * z;

const double m_incub_g = (4.07) * (24 * 60 * 2) * z;
const double sd_incub_g = (2.12) * (24 * 60 * 2) * z;
const double shape_incub_g = pow(m_incub_g, 2) / pow(sd_incub_g, 2);
const double scale_incub_g = pow(sd_incub_g, 2) / m_incub_g;

const double mInf = (1.63) * (24 * 60 * 2) * z;
const double sdInf = (0.5) * (24 * 60 * 2) * z;

const double m_inf_g = (4.07) * (24 * 60 * 2) * z;
const double sd_inf_g = (2.12) * (24 * 60 * 2) * z;
const double shape_inf_g = pow(m_inf_g, 2) / pow(sd_inf_g, 2);
const double scale_inf_g = pow(sd_inf_g, 2) / m_inf_g;


// R UNIQUE(X) FUNCTION
Rcpp::Environment base("package:base");
Function do_unique = base["unique"];
Function do_sample = base["sample"];


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::IntegerVector Get_status_t(
    const Rcpp::DataFrame& global_status,
    const int& t
) {
    Rcpp::IntegerVector t_inf = global_status["t_inf"];
    Rcpp::IntegerVector t_incub = global_status["t_incub"];
    Rcpp::IntegerVector t_recover = global_status["t_recover"];
    Rcpp::IntegerVector status_t (t_inf.size()) ;
    
    for(int j = 0; j < t_inf.size(); j++){
        if (t_inf[j] != -1 && (t+1) >= t_inf[j] &&(t+1) <= t_incub[j]){
            status_t[j] = 1; // individual j is EXPOSED
        } else if (t_inf[j] != -1 && (t+1) > t_incub[j] && (t+1) <= t_recover[j]){ //cpp index begins at 0 & R's at 1, we chose to use R's index for time
            status_t[j] = 2; // individual j is INFECTIOUS
        } else if (t_inf[j] != -1 && (t+1) > t_recover[j]){
            status_t[j] = 3; // individual j is RECOVERED
        } else{
            status_t[j] = 0; // individual j is SUSCEPTIBLE
        }
    }

    return status_t;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Get_status_j(
    const Rcpp::String& id,
    const Rcpp::DataFrame& global_status,
    const Rcpp::DataFrame& admission,
    const int& t
) {
    int status_j = -1; // if returns -1 --> error (id not in admission)
    Rcpp::CharacterVector ids = admission["id"];
    Rcpp::IntegerVector t_inf = global_status["t_inf"];
    Rcpp::IntegerVector t_incub = global_status["t_incub"];
    Rcpp::IntegerVector t_recover = global_status["t_recover"];
    int index_j = -1;
    for (int k = 0; k < admission.nrows(); k++){
        if (id == ids[k]){
            index_j = k;
        }
    }
    // WE CHECK THE STATUS FOR ONLY INDIVIDUAL J (with index_j)
    if (index_j != -1){
        if (t_inf[index_j] != -1 && (t+1) >= t_inf[index_j] && (t+1) <= t_incub[index_j]){
            status_j = 1; // individual j is EXPOSED
        } else if (t_inf[index_j] != -1 && (t+1) > t_incub[index_j] && (t+1) <= t_recover[index_j]){ //cpp index begins at 0 & R's at 1, we chose to use R's index for time
            status_j = 2; // individual j is INFECTIOUS
        } else if (t_inf[index_j] != -1 && (t+1) > t_recover[index_j]){
            status_j = 3; // individual j is RECOVERED
        } else{
            status_j = 0; // individual j is SUSCEPTIBLE
        }
    }
    
    return status_j;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::String Sample_inf(
    const Rcpp::String& id,
    const Rcpp::List& list_inf_encountered,
    const Rcpp::DataFrame& admission,
    const Rcpp::DataFrame& localization_ti,
    const double& lambda_e_j,
    const double& lambda_c_j
) {
    if (lambda_e_j == 0){ // NO ENV --> CONTACT
        int n_ind = list_inf_encountered.size();
        double FOI_j = 1 - exp(-lambda_c_j);
        double ind_weight = (1 - exp(-lambda_c_j)) / (n_ind * FOI_j);
        Rcpp::NumericVector weights(n_ind, ind_weight);  // Initialize with ind_weight, CAUTION, R its rep(1,3) but cpp its v(3,1)

        Rcpp::CharacterVector elements(n_ind);
        for (int i = 0; i < elements.size(); i++) { //
            Rcpp::String id_patient_i = list_inf_encountered[i];
            elements[i] = id_patient_i;
        }
        
        Rcpp::String sampled = Rcpp::sample(elements, 1, true, weights)[0];
        Rcpp::String res("CONTACT-");
        res += sampled; // "sampled" is the id of the individual we've sampled
        return res;

    } else if (lambda_c_j == 0){ // NO CONTACTS -> ENV
        Rcpp::String env_string("ENVIRONMENT-");
        int room_j = Get_loc_j(id,localization_ti);
        Rcpp::String room_string = std::to_string(room_j);
        Rcpp::String res(env_string);
        res += room_string; 
        return res;

    } else { // CONTACT AND ENV
        int n_ind = list_inf_encountered.size(); 
        double FOI_j = 1 - exp(-lambda_c_j -lambda_e_j);
        double ind_weight = (1 - exp(-lambda_c_j)) / (n_ind * FOI_j); // lambda_c !=0 so n_ind!=0
        double env_weight = (1 - exp(-lambda_e_j))/ FOI_j;
        // WEIGHTS VECTOR
        Rcpp::NumericVector weights(n_ind, ind_weight);  // Initialize with ind_weight, CAUTION, R its rep(1,3) but cpp its v(3,1)
        weights.push_front(env_weight); //WEIGHT FOR ENVIRONMENT
        
        // ELEMENTS VECTOR
        Rcpp::CharacterVector elements(n_ind);
        for (int i = 0; i < elements.size(); i++) { //
            Rcpp::String id_patient_i = list_inf_encountered[i];
            elements[i] = id_patient_i;
        }
        Rcpp::String env_string("ENVIRONMENT-");
        elements.push_front(env_string);

        // SAMPLE
        Rcpp::String sampled = Rcpp::sample(elements, 1, true, weights)[0];
        
        // RETURN THE CAUSE OF INFECTION
        if (sampled == env_string){ // IF ENVIRONMENT
            int room_j = Get_loc_j(id, localization_ti);
            Rcpp::String room_string = std::to_string(room_j);
            Rcpp::String res(env_string);
            res += room_string; 
            return res;

        } else { // IF CLOSE CONTACT
            Rcpp::String res("CONTACT-");
            res += sampled; // "sampled" is the id of the individual we've sampled
            return res;
        };
    }
};

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame Update_status_bis(
    const Rcpp::DataFrame& global_status,
    const Rcpp::DataFrame& lambda_ti,
    const Rcpp::DataFrame& admission,
    const Rcpp::DataFrame& interactions_ti,
    const Rcpp::DataFrame& localization_ti,
    const int& t
) {
    Rcpp::NumericVector lambda_c = lambda_ti["lambda_c"];
    Rcpp::NumericVector lambda_e = lambda_ti["lambda_e"];
    
    Rcpp::IntegerVector status_tim1 = Get_status_t(global_status, t);

    Rcpp::CharacterVector ids = admission["id"];
    Rcpp::IntegerVector admission_int = admission["info"];
    
    Rcpp::IntegerVector t_inf_tim1 = global_status["t_inf"];
    Rcpp::IntegerVector t_incub_tim1 = global_status["t_incub"];
    Rcpp::IntegerVector t_recover_tim1 = global_status["t_recover"];
    Rcpp::CharacterVector inf_by_tim1 = global_status["inf_by"];
    Rcpp::IntegerVector inf_room_tim1 = global_status["inf_room"];
    
    Rcpp::DataFrame global_status_updated = clone(global_status);
    Rcpp::IntegerVector t_inf_ti = clone(t_inf_tim1);
    Rcpp::IntegerVector t_incub_ti = clone(t_incub_tim1);
    Rcpp::IntegerVector t_recover_ti = clone(t_recover_tim1);
    Rcpp::IntegerVector inf_room_ti = clone(inf_room_tim1);
    Rcpp::CharacterVector inf_by_ti = clone(inf_by_tim1);


    Rcpp::NumericVector FOI = (lambda_ti.nrows(), 1) - exp(- (lambda_c  + lambda_e));

    for (int j=0; j < lambda_ti.nrows(); j++){
        if (status_tim1[j] == 0 && R::runif(0, 1) <= FOI[j]){
            t_inf_ti[j] = (t+1); // C++ INDEX BEGINS AT 0 / R BEGINS AT 1
            t_incub_ti[j] = t_inf_ti[j] + Incub_period_uniform();
            t_recover_ti[j] = t_incub_ti[j] + Inf_period_uniform(); // C++ INDEX BEGINS AT 0 / R BEGINS AT 1

            // ROOM WHERE j IS INFECTED //
            inf_room_ti[j] = Get_loc_j(ids[j], localization_ti);
            
            // CAUSE OF INFECTION
            Rcpp::String id_j = ids[j];
            Rcpp::List list_inf_encountered = List_inf_encountered(id_j, interactions_ti, global_status, admission, t);
            double lambda_e_j = lambda_e[j];
            double lambda_c_j = lambda_c[j];
            Rcpp::String sampled_j = Sample_inf(ids[j], list_inf_encountered, admission, localization_ti, lambda_e_j, lambda_c_j);
            inf_by_ti[j] = sampled_j;
            }
    };
    global_status_updated["t_inf"] = t_inf_ti;
    global_status_updated["t_incub"] = t_incub_ti;
    global_status_updated["t_recover"] = t_recover_ti;
    global_status_updated["inf_room"] = inf_room_ti;
    global_status_updated["inf_by"] = inf_by_ti;

    return global_status_updated;
};
            

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame Get_t(
    const Rcpp::List& Global_list,
    const int& t
){
    // WARNING: c++ COUNTS STARTS AT 0
    Rcpp::DataFrame df(Global_list[t]);
    return df;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Update_environment(
    const Rcpp::DataFrame& environment_tim1,
    const Rcpp::DataFrame& localization_ti,
    const Rcpp::IntegerVector& status_tim1,
    const Rcpp::DataFrame& admission,
    const double& mu,
    const double& nu,
    const double& deltat
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
    Rcpp::CharacterVector ids = admission["id"];
    Rcpp::IntegerVector admission_int = admission["info"];
    
    double individual_weigth;

    // EXPONENTIAL INACTIVATION 
    env_ti = env_ti * exp(-mu * deltat);
    
    // INFECTING INDIVIDUALS SHEDDING
    for (int j = 0; j < admission.nrows(); ++j) {
        if (admission_int[j] == 0){
            individual_weigth = 1;
        } else if (admission_int[j] == 1){
            individual_weigth = 1; // mask?
        } else {
            individual_weigth = 1;
        }

        if (status_tim1[j] == 2){ // 2 = INFECTIOUS -> SHEDDING
        int room_j = Get_loc_j(ids[j], localization_ti);
            // Search for the index of the room in environment dataframe
            int index_room_j = -1;
            for (int k = 0; k < rooms_environment.size(); ++k) {
                if (rooms_environment[k] == room_j) {
                    index_room_j = k;
                    break;
                }
            }
            if(index_room_j != -1){
                env_ti[index_room_j] += individual_weigth * nu * deltat;}
        }

  return env_ti;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List List_encountered(
    const Rcpp::String& id,
    const Rcpp::DataFrame& interaction_ti
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
Rcpp::List List_inf_encountered(
    const Rcpp::String& id,
    const Rcpp::DataFrame& interaction_ti,
    const Rcpp::DataFrame& global_status,
    const Rcpp::DataFrame& admission,
    const int& t
) {
    Rcpp::List list_encountered = List_encountered(id, interaction_ti);
    Rcpp::List list_inf_encountered;
    for (int j = 0; j < list_encountered.size(); j++){
        Rcpp::String id_encountered = list_encountered[j];
        // 0 = Susceptible, 1 = Exposed, 2 = Infected, 3 = Recovered
        if (Get_status_j(id_encountered, global_status, admission, t) == 2){ //INFECTIOUS
            list_inf_encountered.push_back(id_encountered);
        }
    }

    return list_inf_encountered;
};




//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Lambda_c (
    const Rcpp::DataFrame& lambda_tim1,
    const Rcpp::DataFrame& interaction_ti,
    const Rcpp::IntegerVector& status_ti,
    const double& beta,
    const double& deltat
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
        list_ind_r = List_encountered(id, interaction_ti);
        // Dont use List_inf_encountered because we dont need admission for Lambda_c
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
            // if individual r is INFECTIOUS & we found its index (for safety)
            if (index_r != -1 && status_ti[index_r] == 2) { // 2 = INFECTIOUS
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
    const Rcpp::DataFrame& lambda_tim1,
    const Rcpp::DataFrame& localization_ti,
    const Rcpp::DataFrame& environment_ti,
    const Rcpp::DataFrame& admission, 
    const double& epsilon,
    const double& env_thresold,
    const double& deltat
) {
    Rcpp::CharacterVector ids_lambda = lambda_tim1["id"];
    Rcpp::NumericVector temp_lambda_e = lambda_tim1["lambda_e"];
    Rcpp::NumericVector lambda_e_ti = clone(temp_lambda_e);
    Rcpp::CharacterVector ids_localization = localization_ti["id"];
    Rcpp::IntegerVector localizations = localization_ti["localization"];
    
    Rcpp::IntegerVector admission_int = admission["info"]; // "0" IF PATIENT, "1" IF HCW
    Rcpp::IntegerVector admission_room = admission["room"];
    Rcpp::NumericVector environment = environment_ti["env"];
    Rcpp::IntegerVector rooms_environment = environment_ti["id_room"];
    

  for (int j = 0; j < lambda_tim1.nrows(); ++j){
    Rcpp::String id_j = ids_lambda[j];
    // TWO CASES (PATIENTS AND HCWS)
    // CASE 1. IF INDIVIDUAL j IS A PATIENT --> 
    if (admission_int[j] == 0){
        int room_j = Get_loc_j(id_j, localization_ti)
        // Search for the index of patient's room
        int index_room = -1;
            for (int k = 0; k < rooms_environment.size(); ++k) {
                if (rooms_environment[k] == room_j) {
                    index_room = k;
                    break;
                }
            }
        // VIRAL LOAD THRESOLD
        if (environment[index_room] > env_thresold){
            lambda_e_ti[j] = epsilon * deltat * environment[index_room];
        } else{
            lambda_e_ti[j] = 0;
        }
    }
    
    // CASE 2. IF INDIVIDUAL j IS A HCW --> ENVIRONMENT ACCORDING TO ITS LOCALIZATION
    // WARNING, LOCALIZATION DF INDEX != LAMBDA DF INDEX ETC
    if (admission_int[j] == 1){
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
         
        // VIRAL LOAD THRESOLD
        if (environment[index_room] > env_thresold){
            lambda_e_ti[j] = epsilon * deltat * environment[index_room];
        } else{
            lambda_e_ti[j] = 0;
        }
        
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
int Inf_period_gamma() {
    // Infection period -> Gamma distribution (shape,scale)
    double infection_period_seconds = R::rgamma(shape_inf_g, scale_inf_g);
    int infection_period_subdivisions = static_cast<int>(infection_period_seconds / 30.0);
    
    return infection_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Inf_period_lognormal() {
    // Infection period -> Log-normal distribution (meanlog, sdlog)
    double infection_period_seconds = R::rlnorm(mInf, sd_inf_g);
    int infection_period_subdivisions = static_cast<int>(infection_period_seconds / 30.0);
    
    return infection_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Inf_period_uniform() {
    // Infection period -> Uniform distribution
    double infection_period_days = R::runif(3, 7);
    int infection_period_subdivisions = static_cast<int>( (infection_period_days * 3600*24) / 30.0);
    
    return infection_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Incub_period_uniform() {
    // Incub period -> Uniform distribution
    double incubation_period_days = R::runif(1, 3);
    int incubation_period_subdivisions = static_cast<int>( (incubation_period_days * 3600*24) / 30.0);
    
    return incubation_period_subdivisions;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
int Get_loc_j(
    const Rcpp::String& id,
    const Rcpp::DataFrame& localization_ti
) {
    Rcpp::CharacterVector ids = localization_ti["id"];
    Rcpp::IntegerVector localizations = localization_ti["localization"];
    // Search for individual's index
    int index_j = -1;
    for (int k = 0; k < ids.size(); k++){
        if (id == ids[k]){
            index_j = k;
            break;
        }
    }
    return localizations[index_j];
}
