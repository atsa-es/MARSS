#include "marsskfss.h"
#include "genutils.h"
#include "marssutils.h"

using namespace arma;
using namespace Rcpp;

/**  Exported  **/

// [[Rcpp::export]]
List marsskfss(List& MLEobj)
{
    Function marsskfss_r("MARSSkfss");
    return marsskfss_r(MLEobj);
}
