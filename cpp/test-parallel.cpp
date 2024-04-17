#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct InteractionsWorker : public Worker {
  // Input data
  const NumericVector date_posix;
  const NumericVector length;
  const CharacterVector from;
  const CharacterVector to;
  const double begin_date;
  const int n_subdivisions;
  
  // Output data
  std::vector<std::string> output_from;
  std::vector<std::string> output_to;
  std::vector<int> time;
  
  // Constructor
  InteractionsWorker(const NumericVector date_posix, const NumericVector length, const CharacterVector from, const CharacterVector to, const double begin_date, const int n_subdivisions)
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
            // Rprintf("as<std::string>(from[j]).c_str(): %s\n",as<std::string>(from[j]).c_str());
            // Rprintf("as<std::string>(to[j]).c_str(): %s\n",as<std::string>(to[j]).c_str());
            output_from.push_back(Rcpp::as<std::string>(from[j]));
            output_to.push_back(Rcpp::as<std::string>(to[j]));
            time.push_back(i + 1);
        }
      }
    }
  }
};

// [[Rcpp::export]]
DataFrame parallelInteractions(NumericVector date_posix, NumericVector length, CharacterVector from, CharacterVector to, double begin_date, int n_subdivisions) {
  // Create worker
  int num_workers = 2;
  int chunk_size = (n_subdivisions + num_workers - 1) / num_workers;
  InteractionsWorker worker(date_posix, length, from, to, begin_date, n_subdivisions);
  
  // Parallel execution
  parallelFor(0, n_subdivisions, worker, chunk_size);
  
  // Combine results
  return DataFrame::create(_["from"] = wrap(worker.output_from),
                            _["to"] = wrap(worker.output_to),
                            _["time"] = wrap(worker.time));
}
