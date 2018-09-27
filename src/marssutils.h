#ifndef _MARSSUTILS_H_
#define _MARSSUTILS_H_

#include "marssrcpp.h"


Rcpp::List degentest(const std::string& elem, Rcpp::List& MLEobj, int iter);
bool isdesign(Rcpp::NumericVector&, bool=true, SEXP=R_NilValue,
	      bool=false, bool=false);
Rcpp::LogicalVector isfixed(Rcpp::NumericVector rr, bool byrow=false);
std::map<std::string,int[3]> modeldims2map(Rcpp::List& obj);
Rcpp::List parmat(Rcpp::List& MLEobj,
		  Rcpp::CharacterVector, //elem={"B","U","Q","Z","A","R","x0","V0"},
		  const Rcpp::IntegerVector&, //tt=Rcpp::IntegerVector(1, 0),
		  Rcpp::List dims = Rcpp::List(0), std::string modelloc="marss");
Rcpp::List parmat(Rcpp::List& MLEobj,
		  Rcpp::CharacterVector
		  elem={"B","U","Q","Z","A","R","x0","V0","G","H","L"},
		  int tt=0, Rcpp::List dims=Rcpp::List(0),
		  std::string modelloc= "marss");
std::map<std::string, arma::cube> parmat_cube(Rcpp::List& MLEobj,
					      Rcpp::CharacterVector,
					      const Rcpp::IntegerVector&,
					      Rcpp::List dims = Rcpp::List(0),
					      std::string modelloc="marss");
arma::mat pcholinv(arma::mat);
Rcpp::NumericMatrix sub3D(Rcpp::NumericVector, int);
Rcpp::NumericMatrix sub3Dx(const arma::cube&, const Rcpp::List&, int);
SEXP wrapmessage(const Rcpp::CharacterVector&);

std::string packageVersion(const std::string& pkg);

// KFAS method wrappers

Rcpp::List KFS(const Rcpp::RObject& model, bool simplify);
double logLik(const Rcpp::RObject& object);
Rcpp::List SSModel3(const Rcpp::List& envlist);
#endif /* _MARSSUTILS_H_ */
