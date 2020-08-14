#include "marsshatyt.h"
#include "genutils.h"
#include "marssutils.h"
#include "marsskf.h"

#include <algorithm>
#include <fstream>
#include <string>

using namespace arma;
using namespace Rcpp;

List marsshatyt(Rcpp::List& MLEobjin)
{
  bool fdump = FALSE;
  if (fdump) {
    savelist(List::create(_["mhtop"]=MLEobjin), "mhtop", true);
  }

  List MLEobj = clone(MLEobjin);
  List mle;

  List modelobj = MLEobj["marss"];
  List kflist;
  if (as<List>(MLEobj["kf"]).isNULL()) {
    kflist = marsskf_impl(MLEobj);
  }
  else {
    kflist = MLEobj["kf"];
  }
  
  List model_dims = modelobj.attr("model.dims");
  int n = as<NumericVector>(model_dims["data"])[0];
  int TT = as<NumericVector>(model_dims["data"])[1];
  int m = as<NumericVector>(model_dims["x"])[0];

  NumericMatrix YM(n, TT, as<List>(modelobj["data"]).begin());
  LogicalVector YMna = is_na(YM);
  std::transform(YMna.begin(), YMna.end(), YM.begin(),
		 [](bool x) { return x ? 0 : 1; });
  
  NumericMatrix yp = as<NumericMatrix>(modelobj["data"]);

  mat y(yp.rows(), yp.cols());
  std::transform(yp.begin(), yp.end(), YMna.begin(), y.begin(),
		 [](double ypin, bool yna) {return yna ? 0 : ypin;});

  IntegerVector aargt(1);
  aargt[0] = 0;
  
  NumericMatrix v0p = as<NumericMatrix>(parmat(MLEobj, "V0", aargt)["V0"]);
  mat mv0p = as<mat>(v0p);
  vec dv0p = diagvec(mv0p);
  if (dv0p.size() != (unsigned int)m) { dv0p.resize(m); }
  mat mdv0p = diagmat(dv0p);                             // arma version of IIz$V0
  NumericMatrix IIzV0 = wrap(mdv0p);
  List IIz = List::create(_["V0"] = IIzV0);
      
  mat hatxtT = as<mat>(kflist["xtT"]);                   // was hatxt
  mat hatxtt1 = as<mat>(kflist["xtt1"]);
  NumericMatrix x0T = kflist["x0T"];

  NumericMatrix x0p = as<NumericMatrix>(parmat(MLEobj, "x0", aargt)["x0"]);
  mat Ex0 = (eye(m, m) - mdv0p) * as<mat>(x0T) + mdv0p * as<mat>(x0p);
  
  mat hatxt1 = join_rows(Ex0, hatxtT.cols(0,TT-2));

  cube hatVtT = makecube(kflist["VtT"], true);           // was hatVt
  cube hatVtt1T = makecube(kflist["Vtt1T"], true);       // was hatVtt1

  SEXP msg = R_NilValue;                                 // used in return list
  
  mat In = eye(n, n);
  bool isRdiagonal = true;
  std::vector<std::string> time_varying;
    
  vec diagR;
  
  List pari;
  std::map<std::string, cube> parim;
  
  IntegerVector iis(TT);
  for (int iiis = 0; iiis < TT; ++iiis) {
    iis[iiis] = iiis;
  }

  for (auto elem : {"R", "Z", "A"}) {
    if (as<NumericVector>(model_dims[elem])[2] == 1) {
      parim[elem] = (parmat_cube(MLEobj, elem, {0})[elem]);
      if (elem == "R") {
	diagR = diagvec(parim[elem].slice(0));
  	isRdiagonal = is_diag(parim[elem].slice(0));
      }
    }
    else {
      time_varying.push_back(elem);

      parim[elem] = (parmat_cube(MLEobj, elem, iis)[elem]);      
    }
  }
  
  mat hatytT = zeros(n, TT);
  mat hatytt1 = zeros(n, TT);
  cube hatOt = zeros(n, n, TT);
  cube hatyxt = zeros(n, m, TT);
  cube hatyxtt1 = zeros(n, m, TT);

  for (int t = 0; t < TT; ++t) {
    aargt[0] = t;
    for (auto elem : time_varying) {
      if (elem == "R") {
	diagR = diagvec(parim[elem].slice(t));
	isRdiagonal = is_diag(parim[elem].slice(t));
      }
    }
    
    if (is_true(all(YM.column(t) == 1))) {
      hatytT.col(t) = y.col(t);
      hatytt1.col(t) = y.col(t);
      hatOt.slice(t) = (hatytT.col(t)) * hatytT.col(t).t();
      hatyxt.slice(t) = (hatytT.col(t)) * hatxtT.col(t).t();
      hatyxtt1.slice(t) = (hatytT.col(t)) * hatxt1.col(t).t();
    }
    else {
      mat I2 = In;
      mat Ir = In;
      mat mYm = as<mat>(YM);
      
      I2.rows(find(mYm.col(t)==1)).zeros();
      Ir.rows(find(mYm.col(t)==0 || diagR==0)).zeros();
      
      mat Deltar {In};
      if (isRdiagonal) {Deltar -= Ir;}
      
      mat pariR = parim["R"].n_slices == 1 ? parim["R"].slice(0) :
	    parim["R"].slice(t);
      
      if (!isRdiagonal && is_true(any(YM.column(t) == 1)) &&
	        any(diagR)) {
        
      	mat mhor = Ir.rows(find(mYm.col(t)==1 && diagR!=0));
      	mat tmhor = Ir.cols(find(mYm.col(t)==1 && diagR!=0));
      	mat Rinv;
      	try {
      	  Rinv = chol(mhor * pariR * tmhor);
      	}
      	catch(std::runtime_error) {
      	  return List::create(
      	    _["ok"]=false,
      	    _["errors"]="Stopped in MARSShatyt: chol(R) error.\n"
      	    );
      	}

      	mat choltarget = trimatl(Rinv.t());
      	mat inv1 = solve(choltarget, eye(Rinv.n_rows, Rinv.n_rows));
      	Rinv = solve(trimatu(choltarget.t()), inv1);
      	Deltar = In - pariR * tmhor * Rinv * mhor;
      }
      
      mat pariA = parim["A"].n_slices == 1 ? parim["A"].slice(0) :
	    parim["A"].slice(t);
      
      mat pariZ = parim["Z"].n_slices == 1 ? parim["Z"].slice(0) :
	    parim["Z"].slice(t);
      
      hatytT.col(t) = y.col(t) - Deltar * (y.col(t) - pariZ * hatxtT.col(t) -
					   pariA);
      
      hatytt1.col(t) = pariZ * hatxtt1.col(t) + pariA;
      
      mat tmpDZ = Deltar * pariZ;
      
//      mat tDZ = (m == n) ? tmpDZ.t() : (mat(tmpDZ.begin(), m, n)).t();
      mat tDZ = (m == n) ? tmpDZ.t() : (mat(tmpDZ.begin(), n, m)).t();
      
      hatOt.slice(t) = I2 *
        	(Deltar * pariR + Deltar * pariZ * hatVtT.slice(t) * tDZ ) *
        	I2 + hatytT.col(t) * hatytT.col(t).t();
      
      hatyxt.slice(t) = hatytT.col(t) * hatxtT.col(t).t() +
      	Deltar * pariZ * hatVtT.slice(t);
      
      hatyxtt1.slice(t) = hatytT.col(t) * hatxt1.col(t).t() +
      	Deltar * pariZ * hatVtt1T.slice(t);
    }
  }
  return List::create(
    _["ytT"] = hatytT,
    _["OtT"] = hatOt,
    _["yxtT"] = hatyxt,
    _["yxt1T"] = hatyxtt1,
    _["ytt1"] = hatytt1,
    _["ok"]=true,
    _["errors"]=msg
    );
}
