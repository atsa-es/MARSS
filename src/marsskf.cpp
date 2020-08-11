#include "marsskf.h"
#include "genutils.h"
#include "marssutils.h"
#include "marsskfas.h"
#include "marsskfss.h"

using namespace arma;
using namespace Rcpp;

/**  Exported  **/

// [[Rcpp::export]]
List marsskf(List& MLEobj)
// throws invalid_argument
{
  Function marsskf_r("MARSSkf");
  return marsskf_r(MLEobj);
}

List marsskf_impl(List& MLEobj)
{
  if (MLEobj.containsElementNamed("par") && !as<List>(MLEobj["par"]).isNULL()) {
    if (as<String>(MLEobj["fun.kf"]) == "MARSSkfss") {
      return marsskfss(MLEobj);
    }
    if (as<String>(MLEobj["fun.kf"]) == "MARSSkfas") {
      return marsskfas_impl(MLEobj, false, true, false);
    }
    return List::create(
      _["ok"]=false,
      _["errors"]="kf.function does not specify a valid Kalman filter "
      "and smoother function.");
  }
  else {
    throw std::invalid_argument("MARSSkf: par element of marssMLE object "
				"is required.");
  }
}
