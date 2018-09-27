#include "marsskemcheck.h"
#include "genutils.h"
#include "marssutils.h"
#include "marsshatyt.h"
#include "marsskf.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <sstream>

using namespace arma;
using namespace Rcpp;

void nonzeros_to_ones(const List& tmplist, const std::string& elem)
{
  if (tmplist.containsElementNamed("marss")) {
    if (as<List>(tmplist["marss"]).containsElementNamed("fixed") &&
	as<List>(as<List>(tmplist["marss"])["fixed"]).
	containsElementNamed(elem.c_str())){
      NumericVector fixedel =
	as<List>(as<List>(tmplist["marss"])["fixed"])[elem];
      std::transform(fixedel.begin(), fixedel.end(), fixedel.begin(),
		     [](double d){return d==0 ? 0 : 1;});
    }
    if (as<List>(tmplist["marss"]).containsElementNamed("free") &&
	as<List>(as<List>(tmplist["marss"])["free"]).
	containsElementNamed(elem.c_str())) {
      NumericVector freeel =
	as<List>(as<List>(tmplist["marss"])["free"])[elem];
      std::transform(freeel.begin(), freeel.end(), freeel.begin(),
		     [](double d){return d==0 ? 0 : 1;});
    }
  }
}

vec zero_and_fixed_diags(int pardims, int ifree, int ifixed,
			 const mat& mfree, const mat& mfixed)
{
  uvec diagrows = diag_indices(pardims);
  mat mfreerows = mfree.rows(diagrows);
	
  NumericVector nv = wrap(mfreerows);
  LogicalVector is_fixed_diags = isfixed(nv, true);

  uvec fixedres = any(mfixed.rows(diagrows) == 0, 1);
  LogicalVector zero_rows_diags = wrap(fixedres);

  LogicalVector zero_and_fixed = is_fixed_diags & zero_rows_diags;

  vec zerodiags = as<vec>(zero_and_fixed);
  return zerodiags;
}

List marsskemcheck(Rcpp::List& MLEobj)
{
  return marsskemcheck_impl(MLEobj);
}

List marsskemcheck_impl(List& MLEobjin)
{
  List MLEobj = clone(MLEobjin);
  
  List modelobj = MLEobj["marss"];
  List fixed = modelobj["fixed"];
  List free = modelobj["free"];
    
  List par_dims = modelobj.attr("model.dims");
  int m = as<NumericVector>(par_dims["x"])[0];
  int n = as<NumericVector>(par_dims["y"])[0];
  int TT = as<NumericVector>(par_dims["data"])[1];

  double pseudolim = 1e-8;
  bool ok = true;
  CharacterVector msg;

  if (1 >= TT) {
    msg.push_back("The number of time steps is <=2.\n More than 2 data points "
		  "are needed to estimate parameters.\n");
    ok = false;
  }

  // save some dimensional values for later use

  std::map<std::string,int> fixeddim2;
  std::map<std::string,int> freedim2;
    
  for (auto elem : {"A", "B", "R", "Q", "U", "Z"}) {
    fixeddim2[elem] =
      as<IntegerVector>(as<NumericVector>(fixed[elem]).attr("dim"))[2] - 1;
    freedim2[elem] =
      as<IntegerVector>(as<NumericVector>(free[elem]).attr("dim"))[2] - 1;
  }

  for (int t = 0; t < std::max(freedim2["B"], fixeddim2["B"]); ++t) {
    int ifixed = std::min(t, fixeddim2["B"]);
    mat mfreeB = makecube(as<NumericVector>(free["B"])).slice(ifixed);
    NumericMatrix nmfreeB(mfreeB.n_rows, mfreeB.n_cols, mfreeB.begin());
    if (isfixed(nmfreeB)[0]) {
      NumericMatrix zrB(0,1);
      List tmpMLEobj = List::create(_["marss"]=modelobj,
				    _["par"]=List::create(_["B"]=zrB));
      mat parB = as<mat>(parmat(tmpMLEobj,"B",t)["B"]);
      try {
	cx_vec eigval = eig_gen(parB);
	if (any(abs(eigval) > 1.)) {
	  msg.push_back("All the eigenvalues of B must be within "
			"the unit circle: !any(abs(eigval) > 1)\n");
	  ok = false;
	}
      }
      catch(std::runtime_error e) {
	std::ostringstream errmsg;
	errmsg << "Error in marsskemcheck eig_gen : " <<
	  e.what() << "\n";
	msg.push_back(errmsg.str());
	ok = false;
      }
    }
  }

  std::string el ="R";
  int Tmax = 0;
  for (auto par_test : {"R", "Z", "A"}) {
    Tmax = std::max(Tmax, std::max(fixeddim2[par_test],freedim2[par_test]));
  }

  if (!MLEobj.containsElementNamed("par") ||
      as<List>(MLEobj["par"]).isNULL()) {
  // if (as<List>(MLEobj["par"]).isNULL()) {
    list_assign(MLEobj, "par", MLEobj["start"]);
    // MLEobj["par"] = MLEobj["start"];
  }

  list_assign(MLEobj, "kf", marsskf_impl(MLEobj));
  if (!as<List>(MLEobj["kf"])["ok"]) {
    CharacterVector msgout;
    msgout.push_back(as<String>(as<List>(MLEobj["kf"])["errors"]));
    std::copy(msg.begin(), msg.end(), std::back_inserter(msgout));
    return List::create(_["ok"]=false,
			_["msg"]=msgout);
  }

  list_assign(MLEobj, "Ey", marsshatyt(MLEobj));
  
  CharacterVector msg_tmp;

  uvec diagrows = diag_indices(as<NumericVector>(par_dims[el])[0]);

  std::map<std::string, mat> IIs;
  std::map<std::string, cube> freecubes;
  std::map<std::string, List> dimnames;
    
  for (auto partest : {"Z", "A", "B", "U"}) {
    int tmpdim = as<NumericVector>(par_dims[partest])[1];
    std::string ptest = std::string(partest);
    IIs[ptest] = eye(tmpdim, tmpdim);
  }
    
  for (auto partest : {"Z", "A", "B", "U", "x0"}) {
    std::string ptest = std::string(partest);
    NumericVector free_ptest = free[partest];
    freecubes[ptest] = makecube(free[partest]);
    dimnames[ptest] = free_ptest.attr("dimnames");
  }    

  for (int t=0; t < TT; ++t) {
    int ifixed = std::min(t, fixeddim2[el]);
    int ifree = std::min(t, freedim2[el]);

    vec zerodiags =
      zero_and_fixed_diags(as<NumericVector>(par_dims[el])[0],
			   ifree, ifixed,
			   as<cube>(free[el]).slice(ifree),
			   as<cube>(fixed[el]).slice(ifixed));
    if (any(zerodiags != 0)) {
      mat II0 = diagmat(zerodiags);
      if (t <= Tmax) {
	for (auto partest : {"Z", "A"}) {
	  int ifreepar = std::min(t, freedim2[partest]);
	  mat dpart =
	    as<mat>(sub3Dx(freecubes[std::string(partest)],
			   dimnames[std::string(partest)], ifreepar));
	  mat tmpprod = kron(IIs[std::string(partest)],II0) * dpart;
	  bool par_not_fixed = any(vectorise(tmpprod) != 0);
	  if (par_not_fixed) {
	    std::ostringstream smsg;
	    smsg << "t=" << t << " For method=kem(EM), if an "
	      "element of the diagonal of R is 0, the "
	      "corresponding row of " << partest <<
	      " must be fixed.\n";
	    msg.push_back(smsg.str());
	  }
	}
      }

      mat Z = as<mat>(parmat(MLEobj, "Z", t)["Z"]);
      mat A = as<mat>(parmat(MLEobj, "A", t)["A"]);

      uvec zerodiagsindx = find(zerodiags == 1);
      mat ZR0 = Z.rows(zerodiagsindx);
      mat AR0 = A.rows(zerodiagsindx);

      mat tmpy = as<mat>(as<List>(MLEobj["Ey"])["ytT"]);
      mat yR0 = mat(tmpy.rows(zerodiagsindx)).col(t);
	    
      tmpy = as<mat>(as<List>(MLEobj["kf"])["xtT"]);
      mat yresid = yR0 - ZR0 * tmpy - AR0;

      if (any(vectorise(yresid) > pseudolim)) {
	std::ostringstream smsg;
	smsg << " Z, A, and y for the R=0 rows at time=" << t <<
	  " do not agree. For the R=0 rows, E(y) must equal "
	  "Z * E(x) + A.\n";
	msg_tmp.push_back(smsg.str());
      }
    }
  }

  if (msg_tmp.size() > 0) {
    std::copy(msg_tmp.begin(), msg_tmp.begin()+std::min<int>(9, msg_tmp.size()),
	      std::back_inserter(msg));
  }

  Tmax = std::max(0, std::max(freedim2["Q"], fixeddim2["Q"]));

  mat II0Q1;
  for (int t = 0; t <= Tmax; ++t) {
    el = "Q";
    int ifixed = std::min(t, fixeddim2[el]);
    int ifree = std::min(t, freedim2[el]);
    vec zerodiags =
      zero_and_fixed_diags(as<NumericVector>(par_dims[el])[0],
			   ifree, ifixed,
			   as<cube>(free[el]).slice(ifree),
			   as<cube>(fixed[el]).slice(ifixed));
    mat II0Q = diagmat(zerodiags);

    if (t == 0) {
      II0Q1 = II0Q;
    }
    else {
      if (any(vectorise(II0Q1 != II0Q))) {
	std::ostringstream smsg;
	smsg << "t=" << t << ": The placement of 0 variances in Q"
	  " must be time constant.\n";
	msg.push_back(smsg.str());
	ok = false;
      }
    }
  }
    
  Tmax = 0;
  for (auto partest : {"R", "Q", "U", "B", "Z"}) {
    Tmax = std::max(Tmax, std::max(fixeddim2[partest], freedim2[partest]));
  }

  std::map<std::string, mat> II0;

  int zdim1 =
    as<NumericVector>(as<IntegerVector>(fixed["Z"]).attr("dim"))[1];
  NumericVector tmpparZ(zdim1, 1.);
  tmpparZ.attr("dim") = IntegerVector {zdim1, 1};
  for (int t = 0; t <= Tmax; ++t) {
    for (auto el : {"R", "Q"}) {
      int ifixed = std::min(t, fixeddim2[el]);
      int ifree = std::min(t, freedim2[el]);
      vec zerodiags =
	zero_and_fixed_diags(as<NumericVector>(par_dims[el])[0],
			     ifree, ifixed,
			     as<cube>(free[el]).slice(ifree),
			     as<cube>(fixed[el]).slice(ifixed));
      II0[el] = diagmat(zerodiags);
    }

    for (auto el : {"B", "U"}) {
      int ifree = std::min(t, freedim2[el]);
      List tmpMLEobj = List::create(
	_["marss"]=modelobj,
	_["par"]=List::create(_["Z"]=tmpparZ)
	);

      nonzeros_to_ones(tmpMLEobj, "Z");
      mat parZ = as<mat>(parmat(tmpMLEobj, "Z", t)["Z"]);
      mat dpart =
	as<mat>(sub3Dx(freecubes[std::string(el)],
		       dimnames[std::string(el)], ifree));
      mat tmpprod = kron(IIs[std::string(el)],
			 (II0["R"] * parZ * II0["Q"])) * dpart;
      bool par_not_fixed = any(vectorise(tmpprod) != 0);
      if (par_not_fixed) {
	std::ostringstream smsg;
	smsg << "t=" << t << ": For method = kem (EM), if an element of "
	  "the diagonal of R & Q is 0, the corresponding row of " <<
	  el << "  must be fixed.\n";
	msg.push_back(smsg.str());
      }
    }
  }


  Tmax = std::max(0, std::max(fixeddim2["B"], freedim2["B"]));
  II0.clear();
  for (int t = 0; t <= Tmax; ++t) {
    el = "Q";
    int ifixed = std::min(t, fixeddim2[el]);
    int ifree = std::min(t, freedim2[el]);
    vec zerodiags =
      zero_and_fixed_diags(as<NumericVector>(par_dims[el])[0],
			   ifree, ifixed,
			   as<cube>(free[el]).slice(ifree),
			   as<cube>(fixed[el]).slice(ifixed));
    II0[el] = diagmat(zerodiags);

    el = "B";
    ifree = std::min(t, freedim2[el]);
    mat dpart =
      as<mat>(sub3Dx(freecubes[std::string(el)],
		     dimnames[std::string(el)], ifree));
    bool par_not_fixed =
      any(vectorise(kron(IIs[std::string(el)],II0["Q"])*dpart) != 0);
    if (par_not_fixed) {
      std::ostringstream smsg;
      smsg << "t=" << t << ": If an element of the diagonal of Q is 0, "
	"the corresponding row and col of " << el <<
	"  must be fixed.\n";
      msg.push_back(smsg.str());
      ok = false;
    }
  }

  Tmax = std::max(0, std::max(fixeddim2["U"], freedim2["U"]));
  II0.clear();
  el = "Q";
  int ifixed = std::min(0, fixeddim2[el]);
  int ifree = std::min(0, freedim2[el]);
    
  vec zerodiags =
    zero_and_fixed_diags(as<NumericVector>(par_dims[el])[0],
			 ifree, ifixed,
			 as<cube>(free[el]).slice(ifree),
			 as<cube>(fixed[el]).slice(ifixed));
  II0[el] = diagmat(zerodiags);
    
  mat dpart = as<mat>(sub3Dx(freecubes[std::string("x0")],
			     dimnames[std::string("x0")], 0));
  bool test_adj = any(vectorise(II0[el]*dpart) != 0);

  // skipping loop 214-218, as the loop elements do not vary

  for (int t = 1; t <= Tmax; ++t) {  // t=0 covered in initialization above
    mat dpart =
      as<mat>(sub3Dx(freecubes[std::string("x0")],
		     dimnames[std::string("x0")], 0));
    test_adj = test_adj && any(vectorise(II0[el] * dpart) != 0);
  }
    

  mat adjB_1;
  int bdim1 =
    as<NumericVector>(as<IntegerVector>(fixed["B"]).attr("dim"))[1];
  NumericVector tmpparB(bdim1, 1.);
  tmpparB.attr("dim") = IntegerVector {bdim1, 1};

  if (test_adj) {
    Tmax = 0;
    for (auto partest : {"U", "B"}) {
      Tmax = std::max(Tmax, std::max(fixeddim2[partest], freedim2[partest]));
    }
    for (int t = 0; t <= Tmax; ++t) {
      el = "B";
      int ifree = std::min(0, freedim2[el]);

      List tmpMLEobj = List::create(_["marss"]=modelobj,
				    _["par"]=List::create(_["B"]=tmpparB));

      nonzeros_to_ones(tmpMLEobj, el);
      mat adjB = as<mat>(parmat(tmpMLEobj,el,IntegerVector(1,t))[el]);
      adjB.for_each([](mat::elem_type& val) {val = (val==0) ? 0 : 1;});

      if (t == 0) {
	adjB_1 = adjB;
      }
      else {
	if (any(vectorise(adjB != adjB_1))) {
	  std::ostringstream smsg;
	  smsg << "t=" << t << ": If u^{0} or xi^{0} are estimated, "
	    "the adjacency matrix specified by B must be ",
	    "time constant.\n";
	  msg.push_back(smsg.str());
	}
      }
    }
  }
  Tmax = std::max(0, std::max(fixeddim2["Q"], freedim2["Q"]));

  std::map<std::string,mat>IIz;
  std::map<std::string,mat>IIp;
  std::map<std::string,mat>OMGz;
  std::map<std::string,mat>OMGp;
  std::map<std::string,mat>IId;
  std::map<std::string,mat>IIis;

  el = "Q";
  int diagdim = as<NumericVector>(par_dims[el])[0];
  mat tmpeye = eye(diagdim, diagdim);
  mat IId1;
  mat IIis1;
    
  for (int t = 0; t <= Tmax; ++t) {
    int ifixed = std::min(t, fixeddim2[el]);
    int ifree = std::min(t, freedim2[el]);
    vec zerodiags =
      zero_and_fixed_diags(as<NumericVector>(par_dims[el])[0],
			   ifree, ifixed,
			   as<cube>(free[el]).slice(ifree),
			   as<cube>(fixed[el]).slice(ifixed));
    IIz[el] = diagmat(zerodiags);
    IIp[el] = tmpeye - IIz[el];

    uvec IIzd1 = find(diagvec(IIz[el]) == 1);
    OMGz[el] = tmpeye.rows(IIzd1);

    uvec IIpd1 = find(diagvec(IIp[el]) == 1);
    OMGp[el] = tmpeye.rows(IIpd1);

    int onesBdim =
      as<IntegerVector>(as<NumericVector>(free["B"]).attr("dim"))[1];
    NumericMatrix onesB = wrap(ones(onesBdim, 1));
    List tmpMLEobj = List::create(_["marss"]=modelobj,
				  _["par"]=List::create(_["B"]=onesB));

    nonzeros_to_ones(tmpMLEobj, "B");
	
    mat adjmat = as<mat>(parmat(tmpMLEobj,"B",t)["B"]);
    adjmat.for_each([](mat::elem_type& val) {val = (val==0) ? 0 : 1;});

    mat adjmat_pow_m = pow(adjmat, m);
    mat Q0rows_adjmat = OMGz[el] * adjmat_pow_m;

    if (OMGp[el].n_rows != 0) { // cf use of IIpd1
      if (Q0rows_adjmat.n_rows != 0) {
	mat tmpproduct = Q0rows_adjmat * OMGp[el].t();
	mat tmpndx(tmpproduct.n_rows, 1, fill::zeros);
	for (int i = 0; i < tmpndx.n_rows; ++i) {
	  if (all(tmpproduct.row(i) == 0)) {
	    tmpndx(i, 0) = 1;
	  }
	}
	mat tmp = OMGz[el].t() * tmpndx;
	
	IId[el] = diagmat(tmp);
	IIis[el] = eye(m,m) - IId[el] - IIp[el];
      }
      else {
	IId[el] = tmpeye;
	IIis[el] = zeros(diagdim, diagdim);
      }

      if (t == 0) {
	IId1 = IId[el];
	IIis1 = IIis[el];
      }
      else {
	if (any(vectorise(IId[el] == IId1) != 1)) {
	  std::ostringstream smsg;
	  smsg << "t=" << t << ": The location of the indirectly "
	    "stochastic x's must be time constant.\n";
	  msg.push_back(smsg.str());
	  ok=false;
	}
      }
    }
  }

  return List::create(_["ok"] = ok, _["msg"] = wrapmessage(msg));
}
