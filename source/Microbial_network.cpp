#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector rank_cpp(NumericVector X) {
  int N = X.size();
  NumericVector rank_x(N);
  for(int i = 0; i < N; ++i) {
    int r = 1,s = 1;
    for(int j = 0; j < i; ++j){
      if(X[j] < X[i]){
        r += 1;
      }
      if(X[j] == X[i]){
        s += 1;
      }
    }
    for(int j = i+1; j < N; ++j){
      if(X[j] < X[i]){
        r += 1;
      }
      if(X[j] == X[i]){
        s += 1;
      }
    }
    rank_x[i] = r + (s-1) * 0.5;
  }
  return rank_x;
}

// [[Rcpp::export]]
double pearson_cpp(NumericVector X, NumericVector Y) {
  int N = X.size();
  double sum_X = 0, sum_Y = 0, sum_XY = 0;
  double squareSum_X = 0, squareSum_Y = 0;
  for(int i = 0; i < N; ++i) {
    sum_X = sum_X + X[i];
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + X[i] * Y[i];
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }
  if(sum_X == 0 || sum_Y == 0){
    return 0;
  }
  double corr = 0, corr_up = 0, corr_down = 0;
  corr_up = (N * sum_XY - sum_X * sum_Y);
  corr_down = sqrt((N * squareSum_X - sum_X * sum_X) * (N * squareSum_Y - sum_Y * sum_Y));
  corr = corr_up/corr_down;
  return corr;
}

// [[Rcpp::export]]
double spearman_cpp(NumericVector X, NumericVector Y) {
  int N = X.size();
  X = rank_cpp(X);
  Y = rank_cpp(Y);
  double sum_X = 0, sum_Y = 0, sum_XY = 0;
  double squareSum_X = 0, squareSum_Y = 0;
  for(int i = 0; i < N; ++i) {
    sum_X = sum_X + X[i];
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + X[i] * Y[i];
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }
  if(sum_X == 0 || sum_Y == 0){
    return 0;
  }
  double corr = 0, corr_up = 0, corr_down = 0;
  corr_up = (N * sum_XY - sum_X * sum_Y);
  corr_down = sqrt((N * squareSum_X - sum_X * sum_X) * (N * squareSum_Y - sum_Y * sum_Y));
  corr = corr_up/corr_down;
  return corr;
}

// [[Rcpp::export]]
NumericMatrix pairwist_pearson_cpp(NumericMatrix X) {
  int N = X.nrow();
  NumericMatrix cor_matrix(N,N);
  for(int i = 0; i < (N-1); ++i){
    int m = i + 1;
    for(int j = m; j < N; ++j){
      NumericVector x_vac = X.row(i);
      NumericVector y_vac = X.row(j);
      double res = pearson_cpp(x_vac,y_vac);
      cor_matrix(i,j) = res;
      cor_matrix(j,i) = res;
    }
  }
  return cor_matrix;
}

// [[Rcpp::export]]
NumericMatrix pairwist_spearman_cpp(NumericMatrix X) {
  int N = X.nrow();
  NumericMatrix cor_matrix(N,N);
  for(int i = 0; i < (N-1); ++i){
    int m = i + 1;
    for(int j = m; j < N; ++j){
      NumericVector x_vac = X.row(i);
      NumericVector y_vac = X.row(j);
      double res = pearson_cpp(rank_cpp(x_vac),rank_cpp(y_vac));
      cor_matrix(i,j) = res;
      cor_matrix(j,i) = res;
    }
  }
  return cor_matrix;
}


// [[Rcpp::export]]
NumericMatrix two_pair_pearson_cpp(NumericMatrix X, NumericMatrix Y) {
  int N1 = X.nrow();
  int N2 = Y.nrow();
  NumericMatrix cor_matrix(N1,N2);
  for(int i = 0; i < N1; ++i){
    for(int j = 0; j < N2; ++j){
      NumericVector x_vac = X.row(i);
      NumericVector y_vac = Y.row(j);
      double res = pearson_cpp(x_vac,y_vac);
      cor_matrix(i,j) = res;
    }
  }
  return cor_matrix;
}

// [[Rcpp::export]]
NumericMatrix two_pair_spearman_cpp(NumericMatrix X, NumericMatrix Y) {
  int N1 = X.nrow();
  int N2 = Y.nrow();
  NumericMatrix cor_matrix(N1,N2);
  for(int i = 0; i < N1; ++i){
    for(int j = 0; j < N2; ++j){
      NumericVector x_vac = X.row(i);
      NumericVector y_vac = Y.row(j);
      double res = pearson_cpp(rank_cpp(x_vac),rank_cpp(y_vac));
      cor_matrix(i,j) = res;
    }
  }
  return cor_matrix;
}

// [[Rcpp::export]]
DataFrame edge(NumericMatrix X, CharacterVector taxa) {
  int N = X.nrow();
  CharacterVector vertice_1, vertice_2;
  NumericVector coefficient;
  for(int i = 0; i < (N-1); ++i){
    int m = i + 1;
    for(int j = m; j < N; ++j){
      if(X(i,j)*X(i,j) > 0){
        vertice_1.push_back(taxa(i));
        vertice_2.push_back(taxa(j));
        coefficient.push_back(X(i,j));
      }
    }
  }
  DataFrame df = DataFrame::create(Named("Source") = vertice_1 ,
                                   Named("Target") = vertice_2,
                                   Named("weight") = coefficient );
  return df;
}

// [[Rcpp::export]]
NumericMatrix pair_spearman_cpp(NumericMatrix X, NumericMatrix Y) {
  int Nx = X.nrow();
  int Ny = Y.nrow();
  NumericMatrix cor_matrix(Nx,Ny);
  for(int i = 0; i < (Nx-1); ++i){
    for(int j = 0; j < Ny; ++j){
      NumericVector x_vac = X.row(i);
      NumericVector y_vac = Y.row(j);
      double res = pearson_cpp(rank_cpp(x_vac),rank_cpp(y_vac));
      cor_matrix(i,j) = res;
    }
  }
  return cor_matrix;
}


// [[Rcpp::export]]
List armaPearson(arma::vec x, arma::vec y) {
  int n = x.size();
  double r = arma::as_scalar(arma::cor(x, y));
  double t = r * std::sqrt(n - 2) / std::sqrt(1 - r * r);
  double p = -2 * R::pt(std::abs(t), n - 2, false, false)*-1;
  List result;
  result["r"] = r;
  result["p"] = p;
  return result;
}

// [[Rcpp::export]]
List armaSpearman(NumericVector x, NumericVector y) {
  int n = x.size();
  x = rank_cpp(x);
  y = rank_cpp(y);
  arma::vec x_ranks(x.begin(), n, false);
  arma::vec y_ranks(y.begin(), n, false);
  double rho = arma::as_scalar(arma::cor(x_ranks, y_ranks));
  double t = rho * std::sqrt(n - 2) / std::sqrt(1 - rho * rho);
  double p = -2 * R::pt(std::abs(t), n - 2, false, false)*-1;
  List result;
  result["r"] = rho;
  result["p"] = p;
  return result;
}
// [[Rcpp::export]]
List spareZpearson(NumericVector x, NumericVector y) {
  NumericVector mating = x * y;
  LogicalVector non_zero_indices = mating != 0;
  x = x[non_zero_indices];
  y = y[non_zero_indices]; 
  if(x.size() > 2){
    double mean = Rcpp::mean(x);
    double sd = Rcpp::sd(x);
    NumericVector x_z = (x - mean) / sd;
    mean = Rcpp::mean(y);
    sd = Rcpp::sd(y);
    NumericVector y_z = (y - mean) / sd;
    List result = armaPearson(x_z, y_z);
    return result;
  } else {
    List result;
    result["r"] = 0;
    result["p"] = 1;
    return result;
  }
}


// [[Rcpp::export]]
List spareZSpearman(NumericVector x, NumericVector y) {
  NumericVector mating = x * y;
  LogicalVector non_zero_indices = mating != 0;
  x = x[non_zero_indices];
  y = y[non_zero_indices];
  if(x.size() > 2){
    double mean = Rcpp::mean(x);
    double sd = Rcpp::sd(x);
    NumericVector x_z = (x - mean) / sd;
    mean = Rcpp::mean(y);
    sd = Rcpp::sd(y);
    NumericVector y_z = (y - mean) / sd;
    List result = armaSpearman(x_z, y_z);
    return result;
  } else {
    List result;
    result["r"] = 0;
    result["p"] = 1;
    return result;
  }
}

// [[Rcpp::export]]
List sparse_z_pearson_cpp(NumericMatrix X, NumericMatrix Y) {
  int Nx = X.nrow();
  int Ny = Y.nrow();
  NumericMatrix cor_matrix(Nx,Ny);
  NumericMatrix p_matrix(Nx,Ny);
  
  for(int i = 0; i < (Nx-1); ++i){
    for(int j = 0; j < Ny; ++j){
      NumericVector x_vac = X.row(i);
      NumericVector y_vac = Y.row(j);
      
      NumericVector mating = x_vac * y_vac;
      LogicalVector non_zero_indices = mating != 0;
      x_vac = x_vac[non_zero_indices];
      y_vac = y_vac[non_zero_indices];
      int n = x_vac.size();
      if(n > 2){
        double mean = Rcpp::mean(x_vac);
        double sd = Rcpp::sd(x_vac);
        NumericVector x_z = (x_vac - mean) / sd;
        mean = Rcpp::mean(y_vac);
        sd = Rcpp::sd(y_vac);
        NumericVector y_z = (y_vac - mean) / sd;
        arma::vec xv_z(x_z.begin(), n, false);
        arma::vec yv_z(y_z.begin(), n, false);
        
        double r = arma::as_scalar(arma::cor(xv_z, yv_z));
        double t = r * std::sqrt(n - 2) / std::sqrt(1 - r * r);
        double p = -2 * R::pt(std::abs(t), n - 2, false, false)*-1;
        cor_matrix(i,j) = r;
        p_matrix(i,j) = p;
        
      } else {
        cor_matrix(i,j) = 0;
        p_matrix(i,j) = 1;
      }
    }
  }
  List result;
  result["r"] = cor_matrix;
  result["p"] = p_matrix;
  return result;
}