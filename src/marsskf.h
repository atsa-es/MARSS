#ifndef _MARSSKF_
#define _MARSSKF_

#include "marssrcpp.h"

Rcpp::List marsskf(Rcpp::List& MLEobj);
Rcpp::List marsskf_impl(Rcpp::List& MLEobj);

#endif /* _MARSSKF_ */
