#include <RcppArmadillo.h>


#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]


//' helper function for the loss computation
//'
//' @useDynLib orddid
//' @keywords internal
// [[Rcpp::export]]
arma::mat dat_block_boot(
  const arma::mat &dat,
  const arma::vec &id_cluster,
  const arma::vec &id_cluster_boot,
  const int &max_cluster_size
) {

  int J = id_cluster.n_elem;
  int n_max = max_cluster_size * J;
  arma::mat dat_out = arma::mat(n_max, 3);

  int n_used = 0;
  for (int i = 0; i < id_cluster_boot.n_elem; i++) {
    arma::uvec id_use = arma::find(id_cluster == id_cluster_boot(i));
    int n_add = id_use.n_elem;
    dat_out.rows(n_used, n_used+(n_add-1)) = dat.rows(id_use);
    n_used += n_add;
  }

  return dat_out.rows(0, n_used-1);
}
