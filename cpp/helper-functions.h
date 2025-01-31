#ifndef HELPER_FUN__H
#define HELPER_FUN__H

#include <Rcpp.h>

using namespace Rcpp;

RcppExport DataFrame get_global_interaction(int t, DataFrame interactions);

RcppExport std::list<std::vector<std::string>> get_clusters(DataFrame df);

RcppExport CharacterMatrix get_global_location(
    std::vector<std::vector<std::vector<std::string>>> clusters, 
    DataFrame admission,
    int n_subdivisions,
    DataFrame rooms
);

RcppExport Datetime floor_date_hour(Datetime date_posix);

RcppExport DataFrame trim_interactions_agenda(
    DataFrame interactions, 
    DataFrame admission, 
    DataFrame agenda
);

#endif