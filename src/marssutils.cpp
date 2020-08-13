#include "marssutils.h"
#include "genutils.h"
#include "marsskemcheck.h"
#include "marsskf.h"
#include "marsshatyt.h"

#include <sstream>
#include <stdexcept>
#include <string>

using namespace arma;
using namespace Rcpp;

List degentest(const std::string& elem, List& MLEobj, int iter)
{
  List marssfree = as<List>(MLEobj["marss"])["free"];

  if ((isfixed(marssfree[elem]))[0])
  {
    return List::create(
      _["MLEobj"] = MLEobj,
      _["msg"] = R_NilValue,
      _["set.degen"] = false);
  }
    
  std::string constrtype = as<List>(MLEobj["constr.type"])[elem];
  if (constrtype.find("time-varying") == 0)
  {
    return List::create(
      _["MLEobj"] = MLEobj,
      _["msg"] = R_NilValue,
      _["set.degen"] = false);
  }

  if (!as<bool>(as<List>(MLEobj["control"])["allow.degen"]))
  {
    return List::create(
      _["MLEobj"] = MLEobj,
      _["msg"] = R_NilValue,
      _["set.degen"] = false);
  }
    
  if (as<int>(as<List>(MLEobj["control"])["min.degen.iter"]) >= iter)
  {
    return List::create(
      _["MLEobj"] = MLEobj,
      _["msg"] = R_NilValue,
      _["set.degen"] = false);
  }

  NumericVector nvelem = as<NumericVector>(marssfree[elem]);
  if (!isdesign(nvelem, true, R_NilValue, true, false))
  {
    return List::create(
      _["MLEobj"] = MLEobj,
      _["msg"] = R_NilValue,
      _["set.degen"] = false);
  }

  mat coefmat = as<mat>(as<List>(coef(MLEobj, "matrix", "marss"))[elem]);
  if (!is_diag(coefmat)) {
    return List::create(
      _["MLEobj"] = MLEobj,
      _["msg"] = R_NilValue,
      _["set.degen"] = false);
  }


//  diagonal, not fixed, not time-varying, allow.degen set, iter>min iter
  //  and free is a design matrix
  //  So can proceed
  //  if here then not time-varying

  colvec par_elem = abs(as<NumericVector>(as<List>(MLEobj["par"])[elem]));
  double degenlim = as<List>(MLEobj["control"])["degen.lim"];
  uvec degen_par = find(par_elem < degenlim);
    
  CharacterVector msg_degen;
  bool set_degen = false;

  List model_dims = as<List>(MLEobj["marss"]).attr("model.dims");
  int dim_elem = as<IntegerVector>(model_dims[elem])[0];

  for(unsigned int i = 0; i < degen_par.size(); ++i) {
    // best to use iterators, as we don't know which if any indices
    // are valid.
	
    List MLEobj_tmp = clone(MLEobj);
    as<NumericMatrix>(as<List>(MLEobj_tmp["par"])[elem])(i, 0) = 0;
    List marssfreetmp = as<List>(MLEobj_tmp["marss"])["free"];
    NumericVector elem_elem = marssfreetmp[elem];

    IntegerVector elemdims = elem_elem.attr("dim");

    for (int j = 0; j < elemdims[0]; ++j) {
      elem_elem[j + i * elemdims[0]] = 0;
    }

    List kemcheck = marsskemcheck(MLEobj_tmp);
	
    if (kemcheck["ok"]) {
      List new_kf = marsskf_impl(MLEobj_tmp);
      double loglike_old = MLEobj["logLik"];
      double tmp_cvg2 {1e-10};
      if (new_kf["ok"]) {
	double new_kf_loglike = as<double>(new_kf["logLik"]);
	if (is_finite(loglike_old) && is_finite(new_kf_loglike)) {
	  tmp_cvg2 = new_kf_loglike - loglike_old;
	}
	else {
	  tmp_cvg2 = datum::inf;
	}

	if (tmp_cvg2 < -std::sqrt(datum::eps)) {
	  std::ostringstream msg;
	  msg << "iter=" << iter << " Setting diagonal to 0 "
	    "blocked. logLik was lower in attempt to set "
	    "0 diagonals on " << elem << "logLike old=" <<
	    loglike_old << " new=" << new_kf_loglike << "\n";
	  msg_degen.push_back(msg.str());
	}
	else {
	  MLEobj = MLEobj_tmp;
	  MLEobj["kf"] = new_kf;

	  List Ey = marsshatyt(MLEobj);
	  MLEobj["Ey"] = Ey;
		    
	  if (as<List>(MLEobj["control"])["demean.states"]) {
	    mat kfx0T = as<mat>(new_kf["x0T"]);
	    mat kfxtT = as<mat>(new_kf["xtT"]);
	    mat cmb0t = join_rows(kfx0T, kfxtT);
	    vec xbar(cmb0t.n_rows);
	    for (unsigned int ii=0; ii < cmb0t.n_rows; ++ii) {
	      xbar[ii] = std::accumulate(cmb0t.begin_row(ii),
					 cmb0t.end_row(ii), 0);
	      xbar[ii] /= cmb0t.n_cols;
	    }
	    kfx0T.each_col() -= xbar;
	    kfxtT.each_col() -= xbar;
	    as<NumericMatrix>(as<List>(MLEobj["kf"])["xtT"]) =
	      wrap(kfxtT);
	    as<NumericMatrix>(as<List>(MLEobj["kf"])["x0T"]) =
	      wrap(kfx0T);
	  }
	  MLEobj["logLik"] = new_kf["logLik"];
	  set_degen = true;
	}
      }
      else {
	std::ostringstream msg;
	msg << "iter=" << iter << "MARSSkf returned error in attempt "
	  "to set 0 diagonals for " << elem << "\n";
	msg << "Perhaps Q and R are both going to 0?\n";
	msg_degen.push_back(msg.str());
      }
    }
  }

  return List::create(
    _["MLEobj"] = MLEobj,
    _["msg"] = wrapmessage(msg_degen),
    _["set.degen"] = set_degen
    );
}

SEXP wrapmessage(const CharacterVector& msgin)
{
  // replace zero-length message with NULL;
  // accommodates all.equals() comparisons
    
  SEXP msgreturn;
    
  if (msgin.size() == 0) {
    msgreturn = R_NilValue;
  }
  else {
    msgreturn = msgin;
  }

  return msgreturn;
}

std::string packageVersion(const std::string& pkg)
{
  Function packageVersion_r("packageVersion");
  List x = packageVersion_r(pkg);
  IntegerVector y = x[0];
  std::string s;
  if (y.length() == 3) {
    std::ostringstream oss;
    oss << y[0] << "." << y[1] << "." << y[2];
    s = oss.str();
  }
  return s;
}

// KFAS method wrappers

List KFS(const RObject& model, bool simplify)
{
  Environment KFAS("package:KFAS");
  Function KFS_r = KFAS["KFS"];

  return KFS_r(_["model"]=model, _["simplify"]=simplify);
}


double logLik(const RObject& object) {

  Function logLik_r("logLik.SSModel","KFAS");
  double res = as<double>(logLik_r(object));
  
  return res;
}

List SSModel3(const List& envlist)
{
  Function formula_r("formula");
  RObject form = as<bool>(envlist["return.lag.one"]) ?
    formula_r("yt ~ -1+SSMcustom(Z=stack.Zt, T=stack.Tt, R=stack.Rt, Q=stack.Qt, a1=stack.a1, P1=stack.P1, P1inf=stack.P1inf)") :
    formula_r("yt ~ -1+SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=a1, P1=P1, P1inf=P1inf)");
  
  Environment KFAS("package:KFAS");
  Function SSModel_r = KFAS["SSModel"];

  Environment formenv;
  Environment global_env = Environment::global_env();
  List elnames = envlist.attr("names");

  for (std::string elname : elnames) {
    formenv[elname] = wrap(envlist[elname]);
    global_env[elname] = envlist[elname];
  }

  
  form.attr(".Environment") = formenv;
  NumericVector Ht = global_env["Ht"];
  return SSModel_r(_["formula"]=form,
		   _["H"]=Ht);
}

// end KFAS

bool isdesign(NumericVector& xl, bool strict, SEXP sdim,
	      bool zero_rows_ok, bool zero_cols_ok)
{
  // throws invalid_argument

  NumericVector x = wrap(xl);
    
  RObject test = x.attr("dim");
  if (test.isNULL()) {
    throw std::invalid_argument("isdesign : function requires a 2D or"
				"3D matrix");
  }

  IntegerVector xdims = x.attr("dim");
  if (is_true(any(xdims == 0))) {
    return false;  // empty object
  }

  if ((xdims.size() == 3) && (xdims[2] != 1)) {
    throw std::invalid_argument("isdesign : if 3D, 3rd dim of matrix"
				"must be 1");
  }
  mat xx;
    
  if (sdim != R_NilValue) {
    IntegerVector dim(sdim);
    if (2 != dim.size()) {
      throw std::invalid_argument("isdesign : dim must be length 2 vector");
    }

    // Skip is.numeric() test: IntegerVector won't accept anything that
    // returns FALSE.

    if (!(is_true(all(dim > 0)))) {
      throw std::invalid_argument("isdesign : dim must be positive "
				  "whole number");
    }

    if (!(is_true(all( sapply(dim, is_wholenumber ))))) {
      throw std::invalid_argument("isdesign : dim must be positive "
				  "whole number");
    }

    if (x.size() != dim[0] * dim[1]) {
      throw std::invalid_argument("isdesign : dim is not the right size. "
				  "length(x)=dim[1]*dim[2]");
    }
    xx = unvec(x, dim);
  }
  else {
    xx = mat(x.begin(), xdims[0], xdims[1]);
  }

  // run checks on elements of x or xx, whichever is more convenient

  if (is_true(any(is_na(x)))) { return false; }
  // skip is.numeric() as input is NumericVector
  if (xx.has_nan()) { return false; }

  uvec zeroelements = find(xx == 0);
  if (!strict) {
    xx(zeroelements) = ones<vec>(zeroelements.size());
  }

  uvec nzo = find(xx == 1 || xx == 0);
  if (nzo.size() != x.size()) { return false; }

  if (!zero_cols_ok && xx.n_cols > xx.n_rows) { return false; }

  colvec rsums = sum(xx, 1);
  if (!zero_rows_ok && any(rsums != 1)) { return false; }
  nzo = find(rsums != 0 && rsums != 1);
  if (zero_rows_ok && nzo.size() > 0) { return false; }

  rowvec csums = sum(xx, 0);
  if (!zero_cols_ok && any(csums == 0)) { return false; }

  return true;
}

LogicalVector isfixed(NumericVector rr, bool byrow)
{
  // throws invalid_argument
    
  LogicalVector ret(1, false);

  RObject test = rr.attr("dim");
  if (test.isNULL()) {
    throw std::invalid_argument ("isfixed : free(D) function requires "
				 "a 2D or 3D free(D) matrix");
  }

  IntegerVector rrdims = rr.attr("dim");

  if ((2 > rrdims.size()) || (3 < rrdims.size())) {
    throw std::invalid_argument ("isfixed : function requires "
				 "a 2D or 3D free(D) matrix");
  }

  if (is_true(any(is_nan(rr))) || is_true(any(is_na(rr)))) {
    throw std::invalid_argument("isfixed : free(D) cannot have NAs or NaNs");
  }

  if (rrdims[1] == 0) {
    if (byrow) {
      ret = LogicalVector(rrdims[0], true);
    } else {
      ret[0] = true;
    }
  } else {
    if (is_true(all(rr == 0))) {
      ret[0] = true;
    } else if (byrow) {
      if (rrdims.size() == 2 || rrdims[2] == 0) {
	NumericMatrix rmat(rrdims[0], rrdims[1], rr.begin());
	ret = LogicalVector(rmat.nrow());
	for (int i = 0; i < rmat.nrow(); ++i) {
	  NumericMatrix::Row row = rmat.row(i);
	  ret[i] = is_true(all(row == 0));
	}
      } else {
	cube rcube(rr.begin(), rrdims[0], rrdims[1], rrdims[2], false);
	ret = LogicalVector(rcube.n_rows);
	for (unsigned int i = 0; i < rcube.n_rows; ++i) {
	  mat rowmat(rcube.tube(i, 0, i, rcube.n_cols - 1));
	  ret[i] = !any(vectorise(rowmat));
	}
      }
    }
  }

  return ret;
}

std::map<std::string,int[3]> modeldims2map(List& obj)
{
  std::map<std::string,int[3]> dimmap;
  List model_dims = obj.attr("model.dims");
  CharacterVector modeldimnames = model_dims.attr("names");
  return dimmap;
}

List parmat(List& MLEobj, CharacterVector elem, int tt, List dims,
	    std::string modelloc)
{
  return parmat(MLEobj, elem, IntegerVector(1,tt), dims, modelloc);
}

List parmat(List& MLEobj, CharacterVector elem, const IntegerVector& tt,
	    List dimsin, std::string modelloc)
{
  List model = MLEobj[modelloc];
  List pars = MLEobj["par"];
  List f = model["fixed"];
  List d = model["free"];
  List par_mat = List::create();
  IntegerVector t = tt;

  List dims;
  try {
    if (((int)dimsin.length()) == 0 || dimsin[0] == 0) {
      dims = model.attr("model.dims");
    }
    else {
      dims = clone(dimsin);
    }

    for (String el : elem) {
	  
      std::string tmpel = el;
      List delem;
      List felem;

      NumericVector del = d[el];
      IntegerVector del_dim = del.attr("dim");
      mat mdelem(del.begin(), del_dim[0], del_dim[1]);

      NumericVector thedims = dims[el];

      NumericVector fel = f[el];
      IntegerVector fel_dim = fel.attr("dim");

      mat mfelem(fel.begin(), fel_dim[0], fel_dim[1]);

      cube cpar_cube = (t.length() > 1) ?
	cube(thedims[0], thedims[1], t.length()) : cube(1,1,1,fill::zeros);

      if (t.length() > 1) {
	cube fillcube(thedims[0], thedims[1], t.length());
	fillcube.fill(NA_REAL);
	par_mat.push_back(fillcube, el);
      }

      mat mpar_mat;

      if (del_dim[2] == 1) {
	delem = del;
	IntegerVector ddelem = del.attr("dim");
	delem.attr("dim") = IntegerVector::create(ddelem[0], ddelem[1]);
      }

      if (fel_dim[2] == 1) {
	felem = fel;
	IntegerVector dfelem = fel.attr("dim");
	felem.attr("dim") = IntegerVector::create(dfelem[0], dfelem[1]);
      }

      vec vpar_mat;  // defn required here in order to keep it in scope
      // (hence, to keep mpar_mat from being overwritten)

      NumericVector pel;
      IntegerVector peldims;
      if (pars.containsElementNamed(el.get_cstring())) {
	pel = pars[el];
	peldims = pel.attr("dim");
      } else {
	std::cout << el.get_cstring() << " pars element not found.\n";
	peldims.push_back(0);
      }

      if (t.length() == 1) {
	if ( 1 < del_dim[2]) {
	  int di = t[0] * del_dim[0] * del_dim[1];
	  mdelem = mat(&(del[di]), del_dim[0], del_dim[1]);
	  delem = List::create(mdelem);
	}
	    
	if (1 < fel_dim[2]) {
	  int fi = t[0] * fel_dim[0] * fel_dim[1];
	  mfelem = mat(&(fel[fi]), fel_dim[0], fel_dim[1]);
	  felem = List::create(mfelem);
	}

	const vec& rpar_mat = (0 != peldims[0] && 0 != mdelem.n_cols) ?
	  mfelem + mdelem * vec(pel.begin(), pel.size()) : mfelem;
	par_mat.push_back(mat(rpar_mat.begin(), thedims[0], thedims[1]),
			  el);
      }
      else {
	cube delcube(del.begin(), del_dim[0], del_dim[1], del_dim[2], false);
	cube felcube(fel.begin(), fel_dim[0], fel_dim[1], fel_dim[2], false);

	for (IntegerVector::iterator ii = t.begin(); ii != t.end(); ii++) {

	  if ( 1 < del_dim[2]) {
	    mdelem = delcube.slice(*ii);
	    delem = List::create(mdelem);
	  }
	    
	  if (1 < fel_dim[2]) {
	    mfelem = felcube.slice(*ii);
	    felem = List::create(mfelem);
	  }

	  vpar_mat = mfelem;

	  // allow for rowless matrices (e.g., MARSSkemcheck)
	  if (0 != peldims[0] && 0 != mdelem.n_cols) {
	    vec vpel(pel.begin(), pel.size());
	    vpar_mat += mdelem * vpel;
	  }
	  mpar_mat = mat(vpar_mat.begin(), thedims[0], thedims[1], true);
	  cpar_cube.slice(*ii) = mpar_mat;
	}
	if (par_mat.containsElementNamed(el.get_cstring())) {
	  par_mat[el] = cpar_cube;
	}
	else {
	  par_mat.push_back(cpar_cube, el);
	}
      }
    }
  }
  catch(std::exception e) {
    std::cout << "parmat error caught\n";
    return List(0);
  }
  return par_mat;
}

std::map<std::string, arma::cube> parmat_cube(List& MLEobj, CharacterVector elem,
					      const IntegerVector& tt, List dimsin,
					      std::string modelloc)
{
  List model = MLEobj[modelloc];
  List pars = MLEobj["par"];
  List f = model["fixed"];
  List d = model["free"];
  std::map<std::string, arma::cube> par_mat;
  IntegerVector t = tt;

  List dims;
  try {
    if (((int)dimsin.length()) == 0 || dimsin[0] == 0) { 
      dims = model.attr("model.dims");
    }
    else {
      dims = clone(dimsin);
    }

    for (String el : elem) {
	  
      std::string tmpel = el;
      List delem;
      List felem;

      NumericVector del = d[el];
      IntegerVector del_dim = del.attr("dim");
      mat mdelem(del.begin(), del_dim[0], del_dim[1]);

      NumericVector thedims = dims[el];

      NumericVector fel = f[el];
      IntegerVector fel_dim = fel.attr("dim");

      mat mfelem(fel.begin(), fel_dim[0], fel_dim[1]);

      cube cpar_cube(thedims[0], thedims[1], t.length());

      if (1 < t.length()) {
	cube fillcube(thedims[0], thedims[1], t.length());
	fillcube.fill(NA_REAL);
	par_mat[el] = fillcube;
      }

      mat mpar_mat;

      if (del_dim[2] == 1) {
	delem = del;
	IntegerVector ddelem = del.attr("dim");
	delem.attr("dim") = IntegerVector::create(ddelem[0], ddelem[1]);
      }

      if (fel_dim[2] == 1) {
	felem = fel;
	IntegerVector dfelem = fel.attr("dim");
	felem.attr("dim") = IntegerVector::create(dfelem[0], dfelem[1]);
      }

      vec vpar_mat;  // defn required here in order to keep it in scope
      // (hence, to keep mpar_mat from being overwritten)

      NumericVector pel;
      IntegerVector peldims;
      if (pars.containsElementNamed(el.get_cstring())) {
	pel = pars[el];
	peldims = pel.attr("dim");
      } else {
	std::cout << el.get_cstring() << " pars element not found.\n";
	peldims.push_back(0);
      }

      if (t.length() == 1) {
	if ( 1 < del_dim[2]) {
	  int di = t[0] * del_dim[0] * del_dim[1];
	  mdelem = mat(&(del[di]), del_dim[0], del_dim[1]);
	  delem = List::create(mdelem);
	}
	    
	if (1 < fel_dim[2]) {
	  int fi = t[0] * fel_dim[0] * fel_dim[1];
	  mfelem = mat(&(fel[fi]), fel_dim[0], fel_dim[1]);
	  felem = List::create(mfelem);
	}

	const vec& rpar_mat = (0 != peldims[0] && 0 != mdelem.n_cols) ?
	  mfelem + mdelem * vec(pel.begin(), pel.size()) : mfelem;
	cpar_cube.slice(0) = mat(rpar_mat.begin(), thedims[0], thedims[1]);
      }
      else {
	cube delcube(del.begin(), del_dim[0], del_dim[1], del_dim[2], false);
	cube felcube(fel.begin(), fel_dim[0], fel_dim[1], fel_dim[2], false);

	for (int ii = 0; ii < t.length(); ++ii) {

	  if ( 1 < del_dim[2]) {
	    mdelem = delcube.slice(t(ii));
	    delem = List::create(mdelem);
	  }
	    
	  if (1 < fel_dim[2]) {
	    mfelem = felcube.slice(t(ii));
	    felem = List::create(mfelem);
	  }

	  vpar_mat = mfelem;

	  // allow for rowless matrices (e.g., MARSSkemcheck)
	  if (0 != peldims[0] && 0 != mdelem.n_cols) {
	    vec vpel(pel.begin(), pel.size());
	    vpar_mat += mdelem * vpel;
	  }
	  mpar_mat = mat(vpar_mat.begin(), thedims[0], thedims[1]);

	  cpar_cube.slice(ii) = mpar_mat;
	}
      }
      par_mat[el] = cpar_cube;
    }
  }
  catch(std::exception e) {
    std::cout << "parmat error caught\n";
    par_mat.clear();
    return par_mat;
  }
  return par_mat;
}


mat pcholinv(mat x)
{
  if (x.n_elem == 1) {
    return mat{1./x(0,0)};
  }
    
  mat invx;
  int dimx = x.n_rows;
  vec diagx = x.diag();
  if (!all(diagx)) {
    if (any(diagx)) {
      mat rx = trimmat(x, diagx, true, true);
      mat r = trimatu(chol(rx));
      int dimr = r.n_rows;
      mat inv1 = solve(r, eye(dimr, dimr));
      mat b = solve(trimatl(r.t()), inv1);

      mat omg = eye(dimx, dimx);
      for (int i = dimx - 1; i >= 0; --i) {
	if (diagx(i) == 0) {
	  omg.shed_row(i);
	}
      }

      invx = omg.t() * b * omg;

    } else {
      invx = zeros <mat>(dimx, dimx);
    }
  } else {
    mat r = trimatu(chol(x));
    mat inv1 = solve(r, eye(dimx, dimx));
    invx = solve(trimatl(r.t()), inv1);
  }
  return invx;
}

NumericMatrix sub3D(NumericVector x, int t = 0)
{
  IntegerVector x_dims = x.attr("dim");
  cube xcube(x.begin(), x_dims[0], x_dims[1], x_dims[2], false);
  mat xret = xcube.slice(t);
  NumericMatrix nmxret(xret.n_rows, xret.n_cols, xret.begin());
  List dimnames = x.attr("dimnames");
    
  if (dimnames.length() >= 2) {
    dimnames.erase(2, dimnames.length());
    if (dimnames.length() >= 0) {
      nmxret.attr("dimnames") = dimnames;
    }
  }

  return nmxret;
}

NumericMatrix sub3Dx(const cube& xcube, const List& dimnamesin, int t = 0)
{
  mat xret = xcube.slice(t);
  NumericMatrix nmxret(xret.n_rows, xret.n_cols, xret.begin());

  List dimnames(dimnamesin);
  if (dimnames.length() >= 2) {
    dimnames.erase(2, dimnames.length());
    if (dimnames.length() >= 0) {
      nmxret.attr("dimnames") = dimnames;
    }
  }

  return nmxret;
}

