#ifndef _MARSSKEM_H_
#define _MARSSKEM_H_

#include "marssrcpp.h"

Rcpp::List marsskem(Rcpp::List& MLEobj);

// updaters
std::string Aupdate_impl(Rcpp::List& par1,
			 const std::map<std::string, arma::cube>& star,
			 Rcpp::List& kf, Rcpp::List& Ey, Rcpp::List& timevarying,
			 Rcpp::List& MLEobjiter, int TT);
std::string Bupdate_impl(Rcpp::List par1,
			 const std::map<std::string,arma::cube>& star,
			 Rcpp::List kf, const std::string& kf_x0,
			 const arma::mat& IIzV0,
			 const arma::mat& IImIIzV0, Rcpp::List& timevarying,
			 Rcpp::List MLEobjiter, int TT);
std::string Qupdate_impl(Rcpp::List& par1, Rcpp::List& kf, Rcpp::List& Ey,
			 Rcpp::List& timevarying, Rcpp::List& MLEobjiter,
			 const std::string& kf_x0, const arma::mat& IImIIzV0,
			 const arma::mat& IIzV0, int TT);
std::string Rupdate_impl(Rcpp::List& par1, Rcpp::List& kf, Rcpp::List& Ey,
			 Rcpp::List& time_varying, Rcpp::List& MLEobjiter,
			 Rcpp::List& d, int TT);
std::string Uupdate_impl(int m, Rcpp::List par1,
			 const std::map<std::string,arma::cube>& star,
			 const std::map<std::string,arma::cube>&  IIz,
			 const arma::mat& IIzV0,
			 const arma::mat&  IImIIzV0, const arma::mat& IIm,
			 Rcpp::List Ey, Rcpp::List kf, const std::string& kf_x0,
			 Rcpp::List& timvarying, Rcpp::List MLEobjiter, int TT);
Rcpp::String V0update_impl(Rcpp::List& par1, Rcpp::List& kf,
			   Rcpp::List& MLEobjiter);
Rcpp::String x0update_impl(Rcpp::List& par1,
			   std::map<std::string,arma::cube>& star,
			   Rcpp::List& kf, Rcpp::List& Ey, Rcpp::List& timevarying,
			   Rcpp::List& MLEobjiter, const std::string& kf_x0,
			   std::map<std::string, arma::cube>& IIz,
			   const arma::mat& mIIm, const arma::mat& mIImIIzV0,
			   const arma::mat& mIIzV0, int iter, int TT);


// utilities

Rcpp::List loglog_conv_test(Rcpp::List& iter_record, int iter,
			    int deltaT=9, double tol=0.5,
			    const Rcpp::CharacterVector& params_to_test=
			    {"Z","U","x0","R","Q","A","logLik"});
std::string newkf(Rcpp::List& MLEobj,
		  const std::string& elem,
		  bool controlsafe,
		  bool tagfixed,
		  Rcpp::CharacterVector& msg_kf,
		  Rcpp::CharacterVector& msg_kem,
		  Rcpp::List& kf,
		  Rcpp::List& Ey,
		  int iter);
Rcpp::List rerunkf(const std::string& elem, Rcpp::List& MLEobj, int iter);
Rcpp::String stabilitycheck(Rcpp::List& MLEobjiter,
			    const Rcpp::List& timevarying,
			    const Rcpp::NumericMatrix& XX,
			    const std::string& Xtag,
			    int thirddim, int iter);
std::string stabilitycheck_impl(Rcpp::List& MLEobjiter, bool timevarying,
				const arma::mat& XX, const std::string& Xtag,
				int thirddim, int iter);
std::string V0B_check(const Rcpp::NumericVector& par1V0,
		      const Rcpp::NumericVector& par1B);

#endif /* _MARSSKEM_H_ */
