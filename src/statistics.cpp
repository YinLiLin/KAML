#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "omp_set.h"
#include <R_ext/Print.h>
#include <progress.hpp>
#include "progress_bar.hpp"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

class MinimalProgressBar: public ProgressBar{
    public:
    MinimalProgressBar(const char *str)  {
        _finalized = false;
        _str = str;
    }
    ~MinimalProgressBar() {}
    void display() {}
    void update(float progress) {
        if (_finalized) return;
        int pi = (int)(progress * point_length);
        if(point[pi]){
            point[pi] = false;
            REprintf("\r");
            REprintf(_str);
            REprintf("...[finished %u%]", (int)(progress * 100));
        }
    }
    void end_display() {
    if (_finalized) return;
        REprintf("\r");
        REprintf(_str);
        REprintf("...[finished 100%]");
        REprintf("\n");
        _finalized = true;
    }
    private:
    bool _finalized;
    const char *_str;
    int point_length = 100;
    LogicalVector point = rep(true, point_length);
};

arma::mat GInv(const arma::mat A){
  
  arma::mat ginv;
  if(A.n_rows == 1){
    ginv = 1 / A;
  }else{
    arma::mat U;
    arma::vec s;
    arma::mat V;
    double tol = sqrt(datum::eps);
    
    svd(U,s,V,A);
    U = conv_to<mat>::from(conj(conv_to<cx_mat>::from(U)));
    arma::vec sMax(2); sMax.fill(0);
    sMax[1] = tol * s[0];
    arma::uvec Positive = find(s > sMax.max());
    arma::mat Up = U.cols(Positive);
    Up.each_row() %= 1/s(Positive).t();
    ginv = V.cols(Positive) * Up.t();
  }
  return ginv;
}

template <typename T>
SEXP glm_c(const arma::vec &y, const arma::mat &X, XPtr<BigMatrix> pMat, std::string barhead = "GWAS in process", const bool verbose = true, const int threads = 0){
  
  omp_setup(threads);
  
  MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

  int ind = pMat->ncol();
  int mkr = pMat->nrow();
  int q0 = X.n_cols;

  int y_len = y.n_elem;
  if(y_len != ind)
    throw Rcpp::exception("number of individuals not match.!");

  MinimalProgressBar pb(const_cast<char*>(barhead.c_str()));
  Progress progress(mkr, verbose, pb);

  arma::mat iXX = GInv(X.t() * X);
  arma::mat xy = X.t() * y;
  double yy = sum(y % y);
  arma::mat res(mkr, 3);
  arma::vec snp(ind);
  arma::mat iXXs(q0 + 1, q0 + 1);

  #pragma omp parallel for schedule(dynamic) firstprivate(snp, iXXs)
  for(int i = 0; i < mkr; i++){

    for(int ii = 0; ii < ind; ii++){
      snp[ii] = genomat[ii][i];
    }
    
    double sy = sum(snp % y);
    double ss = sum(snp % snp);
    arma::mat xs = X.t() * snp;
    arma::mat B21 = xs.t() * iXX;
    double t2 = as_scalar(B21 * xs);
    double invB22 = 1 / (ss - t2);
    arma::mat NeginvB22B21 = -1 * invB22 * B21;

    iXXs(q0, q0)=invB22;

    iXXs.submat(0, 0, q0 - 1, q0 - 1) = iXX + invB22 * B21.t() * B21;
    iXXs(q0, span(0, q0 - 1)) = NeginvB22B21;
    iXXs(span(0, q0 - 1), q0) = NeginvB22B21.t();

    // statistics
    arma::mat rhs(xy.n_rows + 1, 1);
    rhs.rows(0, xy.n_rows - 1) = xy;
    rhs(xy.n_rows, 0) = sy;
    arma::mat beta = iXXs * rhs;
    int df = ind - q0 - 1;
    double ve = (yy - as_scalar(beta.t() * rhs)) / df;

    res(i, 0) = beta[q0];
    res(i, 1) = sqrt(iXXs(q0, q0) * ve);
    res(i, 2) = 2 * R::pt(abs(res(i, 0) / res(i, 1)), df, false, false);
    progress.increment();
  }
    
  return wrap(res);
}

// [[Rcpp::export]]
SEXP glm_c(const arma::vec &y, const arma::mat &X, SEXP pBigMat, std::string barhead = "GWAS in process", const bool verbose = true, const int threads = 0){

  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()){
  case 1:
    return glm_c<char>(y, X, xpMat, barhead, verbose, threads);
  case 2:
    return glm_c<short>(y, X, xpMat, barhead, verbose, threads);
  case 4:
    return glm_c<int>(y, X, xpMat, barhead, verbose, threads);
  case 8:
    return glm_c<double>(y, X, xpMat, barhead, verbose, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template <typename T>
SEXP mlm_c(const arma::vec &y, const arma::mat &X, const arma::mat &U, const double vgs, XPtr<BigMatrix> pMat, std::string barhead = "GWAS in process", const bool verbose = true, const int threads = 0){
  
  omp_setup(threads);

  MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

  int ind = pMat->ncol();
  int mkr = pMat->nrow();
  int q0 = X.n_cols;

  int y_len = y.n_elem;
  if(y_len != ind)
    throw Rcpp::exception("number of individuals not match.!");

  MinimalProgressBar pb(const_cast<char*>(barhead.c_str()));
  Progress progress(mkr, verbose, pb);

  arma::mat Uy = U.t() * y;
  arma::mat UX = U.t() * X;
  arma::mat UXUy = UX.t() * Uy;
  arma::mat iUXUX = GInv(UX.t() * UX);
  
  arma::mat res(mkr, 3);    
  arma::vec snp(ind);
  arma::mat iXXs(q0 + 1, q0 + 1);

  #pragma omp parallel for schedule(dynamic) firstprivate(snp, iXXs)
  for(int i = 0; i < mkr; i++){

    for(int ii = 0; ii < ind; ii++){
      snp[ii] = genomat[ii][i];
    }
    arma::mat Us = U.t() * snp;
    arma::mat UXUs = UX.t() * Us;

    double UsUs = as_scalar(Us.t() * Us);
    double UsUy = as_scalar(Us.t() * Uy);
    double B22 = UsUs - as_scalar(UXUs.t() * iUXUX * UXUs);
    double invB22 = 1 / B22;
    arma::mat B21 = UXUs.t() * iUXUX;
    arma::mat NeginvB22B21 = -1 * invB22 * B21;

    iXXs(q0, q0)=invB22;
    iXXs.submat(0, 0, q0 - 1, q0 - 1) = iUXUX + invB22 * B21.t() * B21;
    iXXs(q0, span(0, q0 - 1)) = NeginvB22B21;
    iXXs(span(0, q0 - 1), q0) = NeginvB22B21.t();

    // statistics
    arma::mat rhs(UXUy.n_rows + 1, 1);
    rhs.rows(0, UXUy.n_rows - 1) = UXUy;
    rhs(UXUy.n_rows, 0) = UsUy;
    arma::mat beta = iXXs * rhs;
    int df = ind - q0 - 1;

    res(i, 0) = beta(q0, 0);
    res(i, 1) = sqrt(iXXs(q0, q0) * vgs); 
    res(i, 2) = 2 * R::pt(abs(res(i, 0) / res(i, 1)), df, false, false);
    progress.increment();
  }

  return wrap(res);
}

// [[Rcpp::export]]
SEXP mlm_c(const arma::vec &y, const arma::mat &X, const arma::mat &U, const double vgs, SEXP pBigMat, std::string barhead = "GWAS in process", const bool verbose = true, const int threads = 0){

  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()){
  case 1:
    return mlm_c<char>(y, X, U, vgs, xpMat, barhead, verbose, threads);
  case 2:
    return mlm_c<short>(y, X, U, vgs, xpMat, barhead, verbose, threads);
  case 4:
    return mlm_c<int>(y, X, U, vgs, xpMat, barhead, verbose, threads);
  case 8:
    return mlm_c<double>(y, X, U, vgs, xpMat, barhead, verbose, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template <typename T>
List BigStat(XPtr<BigMatrix> pMat, int threads = 0){

  omp_setup(threads);

  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

  int ind = pMat->ncol();
  int j, k, m = pMat->nrow();
  double p1 = 0.0;
  double scale_mean;
  arma::vec mean(m);
  arma::vec sd(m);

  #pragma omp parallel for private(p1, k)
  for (j = 0; j < m; j++){
    p1 = 0.0;
    for(k = 0; k < ind; k++){
      p1 += bigm[k][j];
    }
    mean[j] = p1 / ind;
  }

  #pragma omp parallel for private(p1, k, scale_mean)
  for (j = 0; j < m; j++){
    p1 = 0.0;
    for(k = 0; k < ind; k++){
      scale_mean = (bigm[k][j] - mean[j]);
      p1 += scale_mean * scale_mean;
    }
    double stemp = sqrt(p1 / (ind - 1));
    sd[j] = stemp == 0 ? 1 : stemp;
  }

  return List::create(Named("mean") = mean, Named("sd") = sd);
}

// [[Rcpp::export]]
List BigStat(SEXP pBigMat, int threads = 0){
  
  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()) {
  case 1:
    return BigStat<char>(xpMat, threads);
  case 2:
    return BigStat<short>(xpMat, threads);
  case 4:
    return BigStat<int>(xpMat, threads);
  case 8:
    return BigStat<double>(xpMat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template <typename T>
SEXP kin_cal_m(XPtr<BigMatrix> pMat, const Nullable<double> SUM = R_NilValue, bool scale = false, const Nullable<NumericVector> wt = R_NilValue, int threads = 0, std::string barhead = " Computing in process", bool verbose = true){

  omp_setup(threads);

  if(verbose)
    Rcout << " Computing GRM under mode: Memory" << endl;

  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

  int n = pMat->ncol();
  int m = pMat->nrow();
  int i = 0, j = 0, k = 0;

  double M;
  if(SUM.isNotNull()){
    M = as<double>(SUM);
  }else{
    M = m;
  }

  MinimalProgressBar pb(const_cast<char*>(barhead.c_str()));

  List Stat = BigStat(pMat, threads);
  NumericVector Mean =  Stat[0];
  NumericVector Sd =  Stat[1];
  
  if(!scale){
    for(int i = 0; i < m; i++){
      Sd[i] = 1;
    }
  }

  arma::mat kin(n, n);
  arma::vec coli(m);
  arma::vec colj(m);

  Progress p(n, verbose, pb);

  if(verbose)
    Rcout << " Scale the genotype matrix and compute Z'Z" << endl;

  if(wt.isNotNull()){
    NumericVector wt_ = as<NumericVector>(wt);

    #pragma omp parallel for schedule(dynamic) firstprivate(coli, colj) private(i, j, k) 
    for(i = 0; i < n; i++){
      for(k = 0; k < m; k++){
        coli[k] = sqrt(wt_[k]) * (bigm[i][k] - Mean[k]) / Sd[k];
      }
      if ( ! Progress::check_abort() ) {
        p.increment();
        for(j = i; j < n; j++){
          for(k = 0; k < m; k++){
            colj[k] = sqrt(wt_[k]) * (bigm[j][k] - Mean[k]) / Sd[k];
          }
          kin(i, j) = kin(j, i) = sum(coli % colj) / M;
        }
      }
    }
  }else{

    #pragma omp parallel for schedule(dynamic) firstprivate(coli, colj) private(i, j, k) 
    for(i = 0; i < n; i++){
      for(k = 0; k < m; k++){
        coli[k] = (bigm[i][k] - Mean[k]) / Sd[k];
      }
      if ( ! Progress::check_abort() ) {
        p.increment();
        for(j = i; j < n; j++){
          for(k = 0; k < m; k++){
            colj[k] = (bigm[j][k] - Mean[k]) / Sd[k];
          }
          kin(i, j) = kin(j, i) = sum(coli % colj) / M;
        }
      }
    }
  }
  return Rcpp::wrap(kin);
}

// [[Rcpp::export]]
SEXP kin_cal_m(SEXP pBigMat, const Nullable<double> SUM = R_NilValue, bool scale = false, const Nullable<NumericVector> wt = R_NilValue, int threads = 0, std::string barhead = " Computing in process", bool verbose = true){

  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()){
  case 1:
    return kin_cal_m<char>(xpMat, SUM, scale, wt, threads, barhead, verbose);
  case 2:
    return kin_cal_m<short>(xpMat, SUM, scale, wt, threads, barhead, verbose);
  case 4:
    return kin_cal_m<int>(xpMat, SUM, scale, wt, threads, barhead, verbose);
  case 8:
    return kin_cal_m<double>(xpMat, SUM, scale, wt, threads, barhead, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template <typename T>
SEXP kin_cal_s(XPtr<BigMatrix> pMat, const Nullable<double> SUM = R_NilValue, bool scale = false, const Nullable<NumericVector> wt = R_NilValue, int threads = 0, bool mkl = false, std::string barhead = " Computing in process", bool verbose = true){

  omp_setup(threads);

  if(verbose)
    Rcout << " Computing GRM under mode: Speed" << endl;

  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

  #ifdef _OPENMP
  #else
    if(!mkl)
      mkl = true;
  #endif
  if(threads == 1)
    mkl = true;

  int n = pMat->ncol();
  int m = pMat->nrow();
  int i = 0, j = 0;
  
  double M;
  if(SUM.isNotNull()){
    M = as<double>(SUM);
  }else{
    M = m;
  }

  MinimalProgressBar pb(const_cast<char*>(barhead.c_str()));

  List Stat = BigStat(pMat, threads);
  NumericVector Mean =  Stat[0];
  NumericVector Sd =  Stat[1];

  if(!scale){
    for(int i = 0; i < m; i++){
      Sd[i] = 1;
    }
  }
  
  arma::mat kin(n, n);
  arma::mat geno(m, n);

  if(verbose)
    Rcout << " Scale the genotype matrix" << endl;

  if(wt.isNotNull()){
    NumericVector wt_ = as<NumericVector>(wt);

    #pragma omp parallel for schedule(dynamic) private(i, j)
    for(i = 0; i < n; i++){
      for(j = 0; j < m; j++){
        geno(j, i) = sqrt(wt_[j]) * (bigm[i][j] - Mean[j]) / Sd[j];
      }
    }
  }else{

    #pragma omp parallel for schedule(dynamic) private(i, j)
    for(i = 0; i < n; i++){
      for(j = 0; j < m; j++){
        geno(j, i) = (bigm[i][j] - Mean[j]) / Sd[j];
      }
    }
  }

  if(verbose)
    Rcout << " Computing Z'Z" << endl;

  if(mkl){
    kin = geno.t() * geno / M;
  }else{

    Progress p(n, verbose, pb);
    arma::colvec coli;

    #pragma omp parallel for schedule(dynamic) private(i, j, coli)
    for(i = 0; i < n; i++){
      coli = geno.col(i);
      if ( ! Progress::check_abort() ) {
        p.increment();
        for(j = i; j < n; j++){
          kin(j, i) = kin(i, j) = sum(coli % geno.col(j)) / M;
        }
      }
    }

  }

  return Rcpp::wrap(kin);
}

// [[Rcpp::export]]
SEXP kin_cal_s(SEXP pBigMat, const Nullable<double> SUM = R_NilValue, bool scale = false, const Nullable<NumericVector> wt = R_NilValue, int threads = 0, bool mkl = false, std::string barhead = " Computing in process", bool verbose = true){

  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()){
  case 1:
    return kin_cal_s<char>(xpMat, SUM, scale, wt, threads, mkl, barhead, verbose);
  case 2:
    return kin_cal_s<short>(xpMat, SUM, scale, wt, threads, mkl, barhead, verbose);
  case 4:
    return kin_cal_s<int>(xpMat, SUM, scale, wt, threads, mkl, barhead, verbose);
  case 8:
    return kin_cal_s<double>(xpMat, SUM, scale, wt, threads, mkl, barhead, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}
