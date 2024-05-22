#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct InteractionsWorker : public RcppParallel::Worker {
  // Input data
  const Rcpp::NumericVector date_posix;
  const Rcpp::NumericVector length;
  const Rcpp::CharacterVector from;
  const Rcpp::CharacterVector to;
  const double begin_date;
  const int n_subdivisions;
  
  // Output data
  Rcpp::CharacterVector output_from;
  Rcpp::CharacterVector output_to;
  Rcpp::NumericVector time;
  
  // Constructor
  InteractionsWorker(const Rcpp::NumericVector date_posix, const Rcpp::NumericVector length, const Rcpp::CharacterVector from, const Rcpp::CharacterVector to, const double begin_date, const int n_subdivisions)
    : date_posix(date_posix), length(length), from(from), to(to), begin_date(begin_date), n_subdivisions(n_subdivisions) {}
  
  // Parallel function
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      double t_begin = begin_date + (i) * 30;
      double t_end = begin_date + (i+1) * 30;
      for (int j = 0; j < date_posix.size(); ++j) {
        double interaction_start = date_posix[j];
        double interaction_end = date_posix[j] + length[j];
        if (interaction_start <= t_begin && interaction_end >= t_end) {
            // output_from.push_back(Rcpp::as<std::string>(from[j]));
            // output_to.push_back(Rcpp::as<std::string>(to[j])); 
            output_from.push_back(from[j]);
            output_to.push_back(to[j]); 
            time.push_back(i + 1);
        }
      }
    }
  }
};

// [[Rcpp::export]]
DataFrame parallelInteractions(Rcpp::NumericVector date_posix, Rcpp::NumericVector length, Rcpp::CharacterVector from, Rcpp::CharacterVector to, double begin_date, int n_subdivisions) {
  // Create worker
  // int num_workers = 2; // temporary
  // int chunk_size = (n_subdivisions + num_workers - 1) / num_workers; // distribute the work load
  InteractionsWorker worker(date_posix, length, from, to, begin_date, n_subdivisions);
  
  parallelFor(0, n_subdivisions, worker);
  
  // Combine results
  return DataFrame::create(_["from"] = wrap(worker.output_from),
                            _["to"] = wrap(worker.output_to),
                            _["time"] = wrap(worker.time));
}
