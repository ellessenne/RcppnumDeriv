#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List genD_cpp(arma::vec const &h0,
                    Rcpp::Function const &func,
                    arma::vec const &x,
                    double const &d,
                    double const &r,
                    double const &v,
                    double const &eps,
                    double const &zero_tol) {

  arma::uword n = x.n_elem;
  arma::vec const f0 = as<arma::vec>(func(Rcpp::NumericVector(x.begin(), x.end())));

  // Define matrices
  arma::mat D(f0.n_elem, n * (n + 3) / 2, arma::fill::zeros);
  arma::mat Daprox(f0.n_elem, r, arma::fill::zeros);
  arma::mat Hdiag(f0.n_elem, n, arma::fill::zeros);
  arma::mat Haprox(f0.n_elem, r, arma::fill::zeros);

  // First derivative and hessian diagonal
  for (arma::uword i = 0; i < n; ++i) {
    arma::vec h = h0;
    for (arma::uword k = 0; k < r; ++k) {
      arma::vec idx(n);
      idx.fill(0.0);
      idx(i) = 1.0;
      arma::vec xup = x + idx % h;
      arma::vec xlo = x - idx % h;
      arma::vec const f1 = as<arma::vec>(func(Rcpp::NumericVector(xup.begin(), xup.end())));
      arma::vec const f2 = as<arma::vec>(func(Rcpp::NumericVector(xlo.begin(), xlo.end())));
      Daprox.col(k) = (f1 - f2) / (2.0 * h(i));
      Haprox.col(k) = (f1 - 2.0 * f0 + f2) / (h(i) * h(i));
      h = h / v;
    }
    for (arma::uword m = 0; m < (r - 1); ++m) {
      for (arma::uword k = 0; k < (r - m - 1); ++k) {
        Daprox.col(k) = (Daprox.col(k + 1) * std::pow(4.0, m + 1) - Daprox.col(k)) / (std::pow(4.0, m + 1) - 1.0);
        Haprox.col(k) = (Haprox.col(k + 1) * std::pow(4.0, m + 1) - Haprox.col(k)) / (std::pow(4.0, m + 1) - 1.0);
      }
    }
    D.col(i) = Daprox.col(0);
    Hdiag.col(i) = Haprox.col(0);
  }
  arma::uword u = (n - 1);

  // 2nd derivative - do lower half of hessian only
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j <= i; j++) {
      u = u + 1;
      if (i == j) {
        D.col(u) = Hdiag.col(i);
      } else {
        arma::vec h = h0;
        for (arma::uword k = 0; k < r; ++k) {
          arma::vec idx(n);
          idx.fill(0.0);
          idx(i) = 1.0;
          arma::vec jdx(n);
          jdx.fill(0.0);
          jdx(j) = 1.0;
          arma::vec xup = x + idx % h + jdx % h;
          arma::vec xlo = x - idx % h - jdx % h;
          arma::vec const f1 = as<arma::vec>(func(Rcpp::NumericVector(xup.begin(), xup.end())));
          arma::vec const f2 = as<arma::vec>(func(Rcpp::NumericVector(xlo.begin(), xlo.end())));
          Daprox.col(k) = (f1 - 2.0 * f0 + f2 - Hdiag.col(i) * h(i) * h(i) - Hdiag.col(j) * h(j) * h(j)) / (2.0 * h(i) * h(j));
          h = h / v;
        }
        for (arma::uword m = 0; m < (r - 1); ++m) {
          for (arma::uword k = 0; k < (r - m - 1); ++k) {
            Daprox.col(k) = (Daprox.col(k + 1) * std::pow(4.0, m + 1) - Daprox.col(k)) / (std::pow(4.0, m + 1) - 1);
          }
        }
        D.col(u) = Daprox.col(0);
      }
    }
  }

  // Return
  return Rcpp::List::create(Rcpp::Named("D") = D,
                            Rcpp::Named("p") = n,
                            Rcpp::Named("f0") = f0,
                            Rcpp::Named("func") = func,
                            Rcpp::Named("x") = x,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("h0") = h0,
                            Rcpp::Named("Daprox") = Daprox,
                            Rcpp::Named("Hdiag") = Hdiag,
                            Rcpp::Named("Haprox") = Haprox,
                            Rcpp::Named("eps") = eps,
                            Rcpp::Named("zero_tol") = zero_tol);

}
