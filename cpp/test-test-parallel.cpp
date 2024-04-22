
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

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
  size_t MAX_SIZE = date_posix.size();
  Rcpp::CharacterVector output_from = Rcpp::CharacterVector(MAX_SIZE);
  Rcpp::CharacterVector output_to = Rcpp::CharacterVector(MAX_SIZE);
  Rcpp::NumericVector time = Rcpp::NumericVector(MAX_SIZE);

  // Constructor
  InteractionsWorker(const Rcpp::NumericVector date_posix, const Rcpp::NumericVector length, const Rcpp::CharacterVector from, const Rcpp::CharacterVector to, const double begin_date, const int n_subdivisions)
      : date_posix(date_posix), length(length), from(from), to(to), begin_date(begin_date), n_subdivisions(n_subdivisions) {}
  // Parallel filtration to extract the interactions taking place at subdivision i
  void operator()(std::size_t begin, std::size_t end) {
    size_t index = 0; // Index to track insertion position
      for (std::size_t i = begin; i < end; ++i) {
          double t_begin = begin_date + i * 30;
          double t_end = begin_date + (i + 1) * 30;
          for (int j = 0; j < date_posix.size(); ++j) {
              double interaction_start = date_posix[j];
              double interaction_end = date_posix[j] + length[j];
              if (interaction_start <= t_begin && interaction_end >= t_end) {
                  // Insertion using the index (no memory issues????????????????)
                  output_from[index] = from[j];
                  output_to[index] = to[j];
                  time[index] = i + 1;
                  // Increment index after insertion
                  ++index;
              }
          }
      }

      // Once the worker has finished its tasks, we can resize the vectors: 
    // Determine the actual size (i.e. index!)
    size_t actual_size = index;

    // Shrink the vectors to the actual size (we initialized them as vectors of size MAX_SIZE > actual_size)
    output_from = Rcpp::CharacterVector(output_from.begin(), output_from.begin() + actual_size);
    output_to = Rcpp::CharacterVector(output_to.begin(), output_to.begin() + actual_size);
    time = Rcpp::NumericVector(time.begin(), time.begin() + actual_size);
  }

  


};

// [[Rcpp::export]]
DataFrame parallelInteractions(Rcpp::NumericVector date_posix, Rcpp::NumericVector length, Rcpp::CharacterVector from, Rcpp::CharacterVector to, double begin_date, int n_subdivisions) {
  // Create worker & call of parallelFor to the workers
  InteractionsWorker worker(date_posix, length, from, to, begin_date, n_subdivisions);
  parallelFor(0, n_subdivisions, worker);
  
  // Combine results
  return DataFrame::create(_["from"] = wrap(worker.output_from),
                            _["to"] = wrap(worker.output_to),
                            _["time"] = wrap(worker.time));
}
