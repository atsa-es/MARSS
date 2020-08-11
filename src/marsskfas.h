#ifndef _MARSSKFAS_
#define _MARSSKFAS_

#include "marssrcpp.h"

Rcpp::List marsskfas_impl(Rcpp::List&, bool only_logLike=false,
			  bool return_lag_one=true, bool return_kfas_model=false);

#endif /* _MARSSKFAS_ */
