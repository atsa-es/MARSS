
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
SEXP eigenAB(Rcpp::NumericMatrix AA, Rcpp::NumericMatrix BB){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::Map<Eigen::MatrixXd> B = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(BB);
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapAB(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenAtB(Rcpp::NumericMatrix AA, Rcpp::NumericMatrix BB){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::Map<Eigen::MatrixXd> B = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(BB);
  Eigen::MatrixXd C = A * B.adjoint();
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapAtB(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B.adjoint();
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenAb(Rcpp::NumericMatrix AA, Rcpp::NumericVector bb){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::Map<Eigen::VectorXd> b = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(bb);
  Eigen::VectorXd C = A * b;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapAb(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> b){
  Eigen::MatrixXd C = A * b;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenaBa(Rcpp::NumericVector aa, Rcpp::NumericMatrix BB){
  Eigen::Map<Eigen::MatrixXd> B = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(BB);
  Eigen::Map<Eigen::VectorXd> a = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(aa);
  Eigen::MatrixXd C = a.adjoint() * (B * a);
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapaBa(const Eigen::Map<Eigen::VectorXd> a, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = a.adjoint() * (B *a);
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigentAB(Rcpp::NumericMatrix AA, Rcpp::NumericMatrix BB){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::Map<Eigen::MatrixXd> B = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(BB);
  Eigen::MatrixXd C = A.adjoint() * B;
  return Rcpp::wrap(C);
}
// [[Rcpp::export]]
SEXP eigenMaptAB(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A.adjoint() * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenAtBC(Rcpp::NumericMatrix AA, Rcpp::NumericMatrix BB, Rcpp::NumericMatrix CC){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::Map<Eigen::MatrixXd> B = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(BB);
  Eigen::Map<Eigen::MatrixXd> C = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(CC);
  Eigen::MatrixXd D = A * B.adjoint() * C;
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenMapAtBC(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
  Eigen::MatrixXd D = A * B.adjoint() * C;
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenABC(Rcpp::NumericMatrix AA, Rcpp::NumericMatrix BB, Rcpp::NumericMatrix CC){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::Map<Eigen::MatrixXd> B = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(BB);
  Eigen::Map<Eigen::MatrixXd> C = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(CC);
  Eigen::MatrixXd D = A * B * C;
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenMapABC(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
  Eigen::MatrixXd D = A * B * C;
  
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenABtC(Rcpp::NumericMatrix AA, Rcpp::NumericMatrix BB, Rcpp::NumericMatrix CC){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::Map<Eigen::MatrixXd> B = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(BB);
  Eigen::Map<Eigen::MatrixXd> C = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(CC);
  Eigen::MatrixXd D = A * B * C.adjoint();
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenMapABtC(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
  Eigen::MatrixXd D = A * B * C.adjoint();
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenSymm(Rcpp::NumericMatrix AA){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  Eigen::MatrixXd D = (A + A.adjoint())/2;
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenMapSymm(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::MatrixXd D = (A + A.adjoint())/2;
  
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP getDeterminant(Rcpp::NumericMatrix AA){
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(AA);
  const double D = A.determinant();
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP getMapDeterminant(const Eigen::Map<Eigen::MatrixXd> A){
  const double D = A.determinant();
  
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP dMat(Rcpp::NumericVector aa, int r, int c){
  return Rcpp::NumericMatrix(r, c, aa.begin());
}