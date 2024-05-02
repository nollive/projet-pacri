#ifndef TEST_STRUCT__H
#define TEST_STRUCT__H

#include <Rcpp.h>

/////////////
// STRUCTS //  
/////////////

typedef struct {
  Rcpp::IntegerVector individual;
  Rcpp::NumericVector lambda;
} lambda_c;

typedef struct{
  Rcpp::IntegerVector individual;
  Rcpp::NumericVector lambda;
} lambda_e ;

typedef struct {
  Rcpp::IntegerVector from;
  Rcpp::IntegerVector to;
}  interaction ;

typedef struct {
  Rcpp::IntegerVector individual;
  Rcpp::IntegerVector position;
} localisation ;

typedef struct {
  Rcpp::NumericVector env;
  Rcpp::IntegerVector room;
} environment ;

typedef struct {
  Rcpp::IntegerVector individual;
  Rcpp::IntegerVector status;
} status ;


// Define custom converters for your structs
namespace Rcpp {

// Convert R list to lambda_c struct
template <>
lambda_c as(SEXP x) {
  Rcpp::List list(x);
  lambda_c result;
  result.individual = list["individual"];
  result.lambda = list["lambda"];
  return result;
}

// Convert lambda_c struct to R list
template <>
SEXP wrap(const lambda_c& lc) {
  Rcpp::List list;
  list["individual"] = lc.individual;
  list["lambda"] = lc.lambda;
  return list;
}

// Convert R list to lambda_e struct
template <>
lambda_e as(SEXP x) {
  Rcpp::List list(x);
  lambda_e result;
  result.individual = list["individual"];
  result.lambda = list["lambda"];
  return result;
}

// Convert lambda_e struct to R list
template <>
SEXP wrap(const lambda_e& le) {
  Rcpp::List list;
  list["individual"] = le.individual;
  list["lambda"] = le.lambda;
  return list;
}

// Convert R list to interaction struct
template <>
interaction as(SEXP x) {
  Rcpp::List list(x);
  interaction result;
  result.from = list["from"];
  result.to = list["to"];
  return result;
}

// Convert interaction struct to R list
template <>
SEXP wrap(const interaction& inter) {
  Rcpp::List list;
  list["from"] = inter.from;
  list["to"] = inter.to;
  return list;
}

// Convert R list to localisation struct
template <>
localisation as(SEXP x) {
  Rcpp::List list(x);
  localisation result;
  result.individual = list["individual"];
  result.position = list["position"];
  return result;
}

// Convert localisation struct to R list
template <>
SEXP wrap(const localisation& loc) {
  Rcpp::List list;
  list["individual"] = loc.individual;
  list["position"] = loc.position;
  return list;
}

// Convert R list to environment struct
template <>
environment as(SEXP x) {
  Rcpp::List list(x);
  environment result;
  result.env = list["env"];
  result.room = list["room"];
  return result;
}

// Convert environment struct to R list
template <>
SEXP wrap(const environment& env) {
  Rcpp::List list;
  list["env"] = env.env;
  list["room"] = env.room;
  return list;
}

// Convert R list to status struct
template <>
status as(SEXP x) {
  Rcpp::List list(x);
  status result;
  result.individual = list["individual"];
  result.status = list["status"];
  return result;
}

// Convert status struct to R list
template <>
SEXP wrap(const status& stat) {
  Rcpp::List list;
  list["individual"] = stat.individual;
  list["status"] = stat.status;
  return list;
}

} // namespace Rcpp


///////////
// Utils // 
///////////

Rcpp::List List_encountered(
    const int id,
    const interaction interactions
);

////////////
// UPDATE // 
////////////

interaction Update_interaction(
    Rcpp::DataFrame data,
    int ti,
    Rcpp::String time,
    const Rcpp::Datetime begin_date, // can change (apply offset of (begin_date) in data --> origin is 00-00-00 00:00:00)
    const int dt
);

environment Update_environment(
    environment environment_tim1,
    const localisation& localisation_tim1,
    const status& status_tim1,
    const double mu,
    const double nu,
    const int dt
);

status Update_status(
    const localisation localisation_ti,
    const environment environment_ti,
    const interaction interaction_ti,
    const status status_ti,
    const status status_tim1,
    lambda_c lambda_c_tim1,
    lambda_e lambda_e_tim1,
    const double alpha,
    const double beta,
    const double epsilon,
    const int dt,
    const int ti
);

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
);

lambda_e Lambda_e (
    lambda_e lambda_e_tim1,
    const localisation localisation_ti,
    const environment environment_ti,
    const double epsilon,
    const double dt
);

#endif