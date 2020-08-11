#include "marsskem.h"
#include "genutils.h"
#include "marssutils.h"
#include "marsshatyt.h"
#include "marsskemcheck.h"
#include "marsskf.h"

#include <algorithm>
#include <exception>
#include <sstream>

using namespace arma;
using namespace Rcpp;

static Paster paster;

// [[Rcpp::export]]
List marsskem(List& MLEobj)
{
  List modelObj = MLEobj["marss"];
  std::string kf_x0 = (as<int>(modelObj["tinitx"]) == 1) ? "x10" : "x00";
  Function describe_marss_r("describe_marss","MARSS");
  List constr_type = describe_marss_r(modelObj);

  List control = MLEobj["control"];
  int control_trace = as<int>(control["trace"]);
  bool control_safe = as<bool>(control["safe"]);
  int control_silent = as<int>(control["silent"]);
  int control_maxit = as<int>(control["maxit"]);

  if (control_trace != -1) {
    List tmp = marsskemcheck(MLEobj);
    if (!tmp["ok"]) {
      std::string errmsg =
	paster("Stopped in MARSSkemcheck due to specification problem(s).\n",
	       "Errors were caught in MARSSkemcheck.\n",
	       as<std::string>(tmp["msg"]),
	       "Try using foo=MARSS(...,fit=FALSE), then summary(foo$model) ",
	       "to see what model you are trying to fit.");
      throw std::logic_error(errmsg);
    }
  }

  std::string stop_msg;
  CharacterVector msg_kem;
  CharacterVector msg_kf;
  
  NumericVector y = clone(as<NumericVector>(modelObj["data"]));
  List d = modelObj["free"];
  List f = modelObj["fixed"];
  List inits = MLEobj["start"];
  
  CharacterVector model_el = modelObj.attr("par.names");
  List model_dims = modelObj.attr("model.dims");
  CharacterVector modeldimnames = model_dims.attr("names");
    
  // copy the model_dims List to a string map for easier access,
  // avoiding as-cast to IntegerVector on each access
  std::map<std::string, int[3]> modeldims;
  for (String mdname : modeldimnames) {
    for (int i = 0; i < 3; ++i) {
      modeldims[mdname][i] = as<IntegerVector>(model_dims[mdname])[i];
    }
  }
    
  int n = modeldims["data"][0];
  int TT = modeldims["data"][1];
  int m = modeldims["x"][0];
    
  mat IIm = eye(m,m);
    
  bool stopped_with_errors = false;
  List kf = List::create();
  kf.attr("names") = CharacterVector();
  double condition_limit = 1e10;

  List MLEobjiter = clone(MLEobj);
  list_assign(MLEobjiter, "constr.type", constr_type);
  list_assign(MLEobjiter, "par", List::create());
  
  List time_varying = List::create();
  std::map<std::string, bool> fixed;
  List mlepar = MLEobjiter["par"];

  for (const String& elemi : model_el) {
    const std::string& elem = elemi.get_cstring();
    bool elemfixed = isfixed(d[elem])[0];
    list_assign(mlepar, elem,
		elemfixed ? NumericMatrix(0,1) : as<List>(MLEobj["start"])[elem]);
    fixed[elem] = elemfixed;
    time_varying.push_back((modeldims[elem][2] != 1), elem);
  }
  MLEobjiter["par"] = mlepar;


  std::map<std::string, bool> set_degen {{"Q",false},{"R",false},{"V0",false}};
  mat tmpL = as<mat>(parmat(MLEobjiter, "L", 0)["L"]);
  mat tmpV0 = as<mat>(parmat(MLEobjiter, "V0", 0)["V0"]);
  mat tmpV = tmpL * tmpV0 * tmpL.t();

  vec tmpzV0 = diagvec(tmpV);
  std::transform(tmpzV0.begin(), tmpzV0.end(), tmpzV0.begin(),
  		 [](double d){ return (d==0) ? 1 : 0; });

  mat IIzV0 = diagmat(tmpzV0);
  mat IImIIzV0 = IIm - IIzV0;
  
  y = ifelse(is_na(y), 0, y);

  List iter_record = List::create(_["par"]=NumericVector::create(),
				  _["logLik"]=NumericVector::create());

  double cvg = 1. + as<double>(control["abstol"]);
  list_assign(MLEobjiter, "logLik", wrap(NA_LOGICAL));
  NumericVector tmpcoef = coef_vec(MLEobjiter);

  list_assign(MLEobjiter, "conv.test",  List::create(
		_["convergence"]=72,
		_["messages"]=
		CharacterVector("No convergence testing performed.\n"),
		_["not.converged.params"]=tmpcoef.attr("names"),
		_["converged.params"]=R_NilValue
		));

  MLEobjiter["conv.test"] = List::create(
    _["convergence"]=72,
    _["messages"]=CharacterVector("No convergence testing performed.\n"),
    _["not.converged.params"]=tmpcoef.attr("names"),
    _["converged.params"]=R_NilValue
    );

  if (control_silent == 2) {std::cout << "EM iteration: ";}

  int iter;
  List kf_last;
  std::map<std::string,cube> IIz, star;
  double loglike_old = 0.;


  for (iter = 0; iter <= control_maxit; ++iter) {
    if (control_silent == 2) { std::cout << " "; }

    kf_last = clone(kf);
    kf = marsskf_impl(MLEobjiter);

    if (!as<bool>(kf["ok"])) {
      if (control_trace > 0) {
	msg_kf.push_back(paster("iter=",iter,as<std::string>(kf["errors"])));
      }
      else {
	for (auto err : as<std::vector<std::string>>(kf["errors"])) {
	  msg_kf.push_back(err);
	}
      }
      stop_msg = paster("Stopped at iter=", iter, " in MARSSkem() because ",
			"numerical errors were generated in the Kalman ",
			"filter.\n");
      stopped_with_errors = true;
      break;
    }

    list_assign(MLEobjiter, "kf", kf);
    list_assign(MLEobjiter, "logLik", kf, "logLik");    
    
    if (as<bool>(control["demean.states"])) {
      mat kfx0T = as<mat>(kf["x0T"]);
      mat kfxtT = as<mat>(kf["xtT"]);
      mat kfjoined = join_rows(kfx0T, kfxtT);
      vec xbar = mean(kfjoined, 1);
      as<List>(MLEobjiter["kf"])["xtT"] = kfxtT.each_col() - xbar;
      as<List>(MLEobjiter["kf"])["x0T"] = kfx0T.each_col() - xbar;
    }

    List Ey = marsshatyt(MLEobjiter);
    
    if (!as<bool>(Ey["ok"])) {
      std::vector<std::string> Eyerrors = Ey.containsElementNamed("errors") ?
	as<std::vector<std::string>>(Ey["errors"]) : std::vector<std::string>{};
      if (control_trace > 0) {
	// MARSShatyt returns at most a single "errors" value
	msg_kf.push_back(paster("iter=", iter,
				as<std::string>(Ey["errors"])));
      }
      else {
	msg_kf = Eyerrors;
      }
      stop_msg = paster("Stopped at iter=",iter," in MARSSkem() because "
			"numerical errors were generated in MARSShatyt.\n");
      stopped_with_errors = true;
      break;
    }

    list_assign(MLEobjiter, "Ey", Ey);

    double mle_loglike = as<double>(MLEobjiter["logLik"]);

    if (iter > 0 && is_finite(loglike_old) && is_finite(mle_loglike)) {
      cvg = mle_loglike - loglike_old;
    }
    if (iter > 1 && cvg < -sqrt(datum::eps)) {
      
      if (control_trace > 0) {
	msg_kem.push_back(paster("iter=", iter, " Loglike DROPPED. old=",
				 loglike_old, " new=", mle_loglike,"\n"));
      }
      else {
	msg_kem.erase(msg_kem.begin(), msg_kem.end());
	msg_kem.push_back("MARSSkem: The soln became unstable "
			  "and logLik DROPPED.\n");
      }
    }

    if (control_trace > 0) {
      NumericVector tmpcoef = coef_vec(MLEobjiter);
      iter_record["par"] = srbind(as<NumericVector>(iter_record["par"]),
				  tmpcoef);
      
      as<NumericVector>(iter_record["logLik"]). push_back(MLEobjiter["logLik"]);
      if (!as<List>(as<List>(MLEobjiter["kf"])["errors"]).isNULL()) {
	msg_kf.push_back(paster("iter=",iter,as<std::string>(kf["errors"])));
      }
      list_assign(MLEobjiter, "iter.record", iter_record);
    }
    else {
      NumericVector tmpcoef = coef_vec(MLEobjiter);
      iter_record["par"] = srbind(iter_record["par"], tmpcoef);

      NumericVector tmpll = iter_record["logLik"];
      tmpll.push_back(as<double>(MLEobjiter["logLik"]));
      iter_record["logLik"] = tmpll;
      NumericVector tmppar = iter_record["par"];

      
      int tmp_len = tmpll.length();

      if (tmp_len > as<int>(control["conv.test.deltaT"])+1) {
	NumericMatrix tmppar = iter_record["par"];
	NumericVector tmpll = iter_record["logLik"];

	int maxdelete = tmp_len - as<int>(control["conv.test.deltaT"]) - 1;
	maxdelete = std::max<int>(maxdelete, 0);
	maxdelete = std::min<int>(maxdelete, tmpll.length());
	if (maxdelete > 0) {
	  NumericMatrix reducedpar = tmppar(Range(maxdelete,tmppar.nrow()-1),
					    Range(0,tmppar.ncol()-1));
	  reducedpar.attr("dimnames") = tmppar.attr("dimnames");
	  tmppar = reducedpar;
	}
	iter_record["par"] = tmppar;
	
	tmpll.erase(0, maxdelete);
	iter_record["logLik"] = tmpll;
      }
      list_assign(MLEobjiter, "iter.record", iter_record);
    }

    if (iter >= as<int>(control["minit"]) - 1) {
      if (cvg > 0 && cvg < as<double>(control["abstol"])) {
	if (iter >= as<int>(control["min.iter.conv.test"])-1) {
	  MLEobjiter["conv.test"] =
	    loglog_conv_test(iter_record, iter,
			     as<int>(control["conv.test.deltaT"]),
			     as<double>(control["conv.test.slope.tol"]));
	  if (as<int>(as<List>(
			MLEobjiter["conv.test"])["convergence"])!= 1){
	    break;
	  }
	}
	else {
	  as<List>(MLEobjiter["conv.test"])["convergence"] = 3;
	}
      }
      else {
	as<List>(MLEobjiter["conv.test"])["convergence"] = 4;
      }
    }

    if (iter >= control_maxit) {
      iter = control_maxit - 1;
      break;
    }

    loglike_old = MLEobjiter["logLik"];

    List par1 = parmat(MLEobjiter, {"B","U","Q","Z","A","R","x0","V0"}, 0);

    if (as<bool>(control["allow.degen"])) {

      List tmp = degentest("R", MLEobjiter, iter);
      MLEobjiter = tmp["MLEobj"];

      if (tmp["msg"] != R_NilValue) {
	CharacterVector tmpmsg = as<CharacterVector>(tmp["msg"]);
	cv_append(msg_kem, tmpmsg);
      }
      if (as<bool>(tmp["set.degen"])) {
	d["R"] = as<List>(as<List>(MLEobjiter["marss"])["free"])["R"];
	f["R"] = as<List>(as<List>(MLEobjiter["marss"])["fixed"])["R"];
	kf = MLEobjiter["kf"];
	Ey = MLEobjiter["Ey"];
	fixed["R"] = isfixed(d["R"]);
	set_degen["R"] = true;
      }
    }


    mat dR;
    mat tdRdR;
    List Rcolnames;
    if (!fixed["R"]) {
	std::string resp = Rupdate_impl(par1, kf, Ey,
					time_varying, MLEobjiter, d, TT);
	
	if (!resp.empty()) {
	  stop_msg = paster("1Stopped at iter=",iter," in MARSSkem: ",resp,"\n");
	  stopped_with_errors = true;
	  break;
	}

	stop_msg = stabilitycheck_impl(MLEobjiter, time_varying["R"],
				       as<mat>(par1["R"]), "R",
				       modeldims["R"][2],
				       iter);
	if (!std::string(stop_msg).empty()) {
	  stopped_with_errors = true;
	  break;
	}

	stop_msg = newkf(MLEobjiter, "R", control_safe, fixed["R"],
			 msg_kf, msg_kem, kf, Ey, iter);
	if (!stop_msg.empty()) {
	  stopped_with_errors = true;
	  break;
	}
    }

    if (as<bool>(control["allow.degen"])) {
      List tmp = degentest("Q", MLEobjiter, iter);
      MLEobjiter = tmp["MLEobj"];
      if (tmp["msg"] != R_NilValue) {
	CharacterVector tmpmsg = as<CharacterVector>(tmp["msg"]);
	cv_append(msg_kem, tmpmsg);
      }
      if (as<bool>(tmp["set.degen"])) {
	d["Q"] = as<List>(as<List>(MLEobjiter["marss"])["free"])["Q"];
	f["Q"] = as<List>(as<List>(MLEobjiter["marss"])["fixed"])["Q"];
	kf = MLEobjiter["kf"];
	Ey = MLEobjiter["Ey"];
	fixed["Q"] = isfixed(d["Q"]);
	set_degen["Q"] = true;
      }
    }

    // Do the regular EM update
    if (!fixed["Q"]) {

      std::string resp = Qupdate_impl(par1, kf, Ey, time_varying, MLEobjiter,
				      kf_x0, IImIIzV0, IIzV0, TT);
      //Mku 344
      if (!resp.empty()) {
	stop_msg = paster("1Stopped at iter=",iter," in MARSSkem: ",resp,"\n");
	stopped_with_errors = true;
	break;
      }

      stop_msg = stabilitycheck_impl(MLEobjiter, time_varying["Q"],
				     as<mat>(par1["Q"]), "Q",
				     modeldims["Q"][2],
				     iter);
      if (!std::string(stop_msg).empty()) {
	stopped_with_errors = true;
	break;
      }

      stop_msg = newkf(MLEobjiter, "Q", control_safe, fixed["Q"],
		       msg_kf, msg_kem, kf, Ey, iter);
      if (!stop_msg.empty()) {
	stopped_with_errors = true;
	break;
      }
    }

    if (!fixed["Q"] || !fixed["R"] || !fixed["V0"] ||
	set_degen["Q"] || set_degen["R"] ||
	iter == 0) {
      std::vector<std::string> elems{"Q", "R", "V0"};
      std::string elem1;
      
      if (iter == 0) {
	// skipping IId, IIis
	for (auto elem : elems) {
	  elem1 = "L";
	  if (elem == "Q") { elem1 = "G"; }
	  if (elem == "R") { elem1 = "H"; }
	  Dimension dims(modeldims[elem][0],modeldims[elem][0],modeldims[elem][2]);
	  star[elem] = IIz[elem] = cube(modeldims[elem1][0],
					modeldims[elem1][0],
					std::max<int>(modeldims[elem][2],
						      modeldims[elem1][2]),
					fill::zeros);
	}
      }
      else {
	if (fixed["V0"]) {
	  elems.erase(std::remove(elems.begin(), elems.end(), "V0"), elems.end());
	}
	if (fixed["Q"] && !set_degen["Q"]) {
	  elems.erase(std::remove(elems.begin(), elems.end(), "Q"), elems.end());
	}
	if (fixed["R"] && !set_degen["R"]) {
	  elems.erase(std::remove(elems.begin(), elems.end(), "R"), elems.end());
	}
      }

      for (auto elem : elems) {
	int thedim = modeldims[elem][0];
	int maxT = modeldims[elem][2];

	IntegerVector iis = IdxVec::get(maxT);
	
	cube complist;
	if (maxT > 1) {
	  complist = as<cube>(parmat(MLEobjiter, elem, iis)[elem]);
	}
	
	for (int i = 0; i < maxT; ++i) {
	  mat pari = (i > 0) ? complist.slice(i) :
	    as<mat>(parmat(MLEobjiter, elem, i)[elem]);
	  if (((set_degen[elem] | iter) == 0)) {
	    vec tmppd = pari.diag();
	    std::transform(tmppd.begin(), tmppd.end(), tmppd.begin(),
			   [](double d){ return (d==0) ? 1 : 0; });
	    IIz[elem].slice(i) = diagmat(tmppd);
	    if (elem == "Q" && maxT != 0) {
	      if (any(vectorise(IIz[elem].slice(0)!=IIz[elem].slice(i)))){
		stop_msg=paster("Stopped at iter=",iter," in MARSSkem. IIz$Q "
				"(location of 0s on diagonal) must be time "
				"invariant.\n You probably want to set "
				"allow.degen=FALSE if it is true.\n");
		stopped_with_errors = true;
		break;
	      }
	    }
	  }
	  mat tmppari = pari;
	  star[elem].slice(i) = choleskyinverter(tmppari);
	}
	set_degen[elem] = false;
	if (elem == "V0" && (iter == 0 || set_degen["V0"])) {
	  IIzV0 = as<mat>(sub3Dx(IIz["V0"], List(0), 0));
	  IImIIzV0 = IIm - IIzV0;
	}
      }
      if (stopped_with_errors) break;
    }

    if (!fixed["x0"]) {

      std::string resp = x0update_impl(par1, star, kf, Ey, time_varying,
				       MLEobjiter, kf_x0, IIz, IIm,
				       IImIIzV0, IIzV0, iter, TT);
      if (!resp.empty()) {
	stop_msg=paster("1Stopped at iter=",iter," in MARSSkem: ",resp,"\n");
	stopped_with_errors = true;
	break;
      }

      stop_msg = newkf(MLEobjiter, "x0", control_safe, fixed["x0"],
		       msg_kf, msg_kem, kf, Ey, iter);
      if (!stop_msg.empty()) {
	stopped_with_errors = true;
	break;
      }
    }

    if (!fixed["V0"]) {    
      std::string resp = V0update_impl(par1, kf, MLEobjiter);

      if (!resp.empty()) {
	stop_msg=paster("1Stopped at iter=",iter," in MARSSkem: ",resp,"\n");
	stopped_with_errors = true;
	break;
      }

      resp = V0B_check(par1["V0"], par1["B"]);

      if (!stop_msg.empty()) {
	stopped_with_errors = true;
	break;
      }

      stop_msg = newkf(MLEobjiter, "V0", control_safe, fixed["V0"],
		       msg_kf, msg_kem, kf, Ey, iter);
      if (!stop_msg.empty()) {
	stopped_with_errors = true;
	break;
      }
	    
    }

    if (!fixed["A"]) {
      std::string resp = Aupdate_impl(par1, star, kf, Ey, time_varying,
				      MLEobjiter, TT);
      if (!resp.empty()) {
	stop_msg=paster("1Stopped at iter=",iter," in MARSSkem: ",resp,"\n");
	stopped_with_errors = true;
	break;
      }

      stop_msg = newkf(MLEobjiter, "A", control_safe, fixed["A"],
		       msg_kf, msg_kem, kf, Ey, iter);
      if (!stop_msg.empty()) {
	stopped_with_errors = true;
	break;
      }
    }

    if (!fixed["U"]) {
      std::string resp = Uupdate_impl(m, par1, star, IIz, IIzV0,
				      IImIIzV0, IIm, Ey, kf, kf_x0,
				      time_varying, MLEobjiter, TT);
      if (!resp.empty()) {
	stop_msg=paster("Stopped at iter=",iter," in MARSSkem: ",resp,"\n");
	stopped_with_errors = true;
	break;
      }

      stop_msg = newkf(MLEobjiter, "V0", control_safe, fixed["V0"],
		       msg_kf, msg_kem, kf, Ey, iter);
      if (!stop_msg.empty()) {
	stopped_with_errors = true;
	break;
      }
    }

    if (!fixed["B"]) {
      std::string resp = Bupdate_impl(par1, star, kf, kf_x0, IIzV0, IImIIzV0,
				      time_varying, MLEobjiter, TT);

	if (!resp.empty()) {
	stop_msg=paster("1Stopped at iter=",iter," in MARSSkem: ",resp,"\n");
	stopped_with_errors = true;
	break;
      }

      stop_msg = newkf(MLEobjiter, "B", control_safe, fixed["A"],
		       msg_kf, msg_kem, kf, Ey, iter);
      if (!stop_msg.empty()) {
	stopped_with_errors = true;
	break;
      }
    }

    // Skipping Z update : can't trigger
    
  }
  
  if (control_silent == 2) { std::cout << "\n"; }

  List MLEobj_return = clone(MLEobj);
  list_assign(MLEobj_return, "iter.record", iter_record);
  list_assign(MLEobj_return, "numIter", wrap(iter+1));

  if (stopped_with_errors) {
    if (control_silent == 2) {
      std::cout << "Stopped due to numerical instability or errors. "
	"Print $errors from output for info or set silent=FALSE.\n";
    }

    CharacterVector msg;
    msg.push_back(stop_msg);
    msg.push_back("par, kf, states, iter, loglike are the last values "
		  "before the error.\n");
    if (control_safe) {
      msg.push_back("Try control$safe=TRUE which uses a slower "
		   "but slightly more robust algorithm.\n");
    }

    if (control_trace > 0) {
      msg.push_back("Use control$trace=1 to generate a more detailed "
		   "error report. See user guide for insight.\n");
    }
    if (msg_kem.size() > 0) {
      msg.push_back("\nMARSSkem errors. Type MARSSinfo() for help.\n");
      cv_append(msg, msg_kem);
    }
    if (msg_kf.size() > 0) {
      msg.push_back("\nMARSSkf errors. Type MARSSinfo() for help.\n");
      cv_append(msg, msg_kf);
    }

    MLEobj_return["errors"] = CharacterVector(msg.begin(), msg.end());
    MLEobj_return["par"] = MLEobjiter["par"];
    MLEobj_return["kf"] = kf_last;
    list_assign(MLEobj_return, "states", kf_last["xtT"]);
    MLEobj_return["convergence"] = 52;
    MLEobj_return["logLik"] = MLEobjiter["logLik"];

    return MLEobj_return;
  }

  CharacterVector msg_conv;
    
  bool catinfo = control_silent;
  list_assign(MLEobj_return, "convergence", wrap(72));

  if (as<int>(as<List>(MLEobjiter["conv.test"])["convergence"]) == 72 ) {
    list_assign(MLEobj_return, "convergence", wrap(52));
    msg_conv = as<CharacterVector>(as<List>(MLEobjiter["conv.test"])["messages"]);
    if (catinfo) {
      std::cout << "Error! EM algorithm exited at iter=" << iter << "before minit"
	" reached.\nMinit was " << as<int>(control["minit"]) << ".\n";
    }
  }
  if (as<int>(as<List>(MLEobjiter["conv.test"])["convergence"]) < 0) {
    list_assign(MLEobj_return, "convergence", wrap(62));
    msg_conv = as<CharacterVector>(as<List>(MLEobjiter["conv.test"])["messages"]);
    if (catinfo) {
      std::cout << "Error! EM algorithm exited due to errors reported by "
	"log-log test function.\n";
    }
  }
  if (as<int>(as<List>(MLEobjiter["conv.test"])["convergence"]) == 4) {
    std::string tmp_msg = paster("Warning! Reached maxit before parameters "
				 "converged. Maxit was ", control_maxit, ".\n");
    if (iter > as<int>(control["min.iter.conv.test"])) {
      List loglog_test =
	loglog_conv_test(iter_record, iter, as<int>(control["conv.test.deltaT"]),
			 as<double>(control["conv.test.slope.tol"]));
      as<List>(MLEobjiter["conv.test"])["messages"]= loglog_test["messages"];

      int llc = as<int>(loglog_test["convergence"]);

      if (llc < 0) {
	list_assign(MLEobj_return, "convergence", wrap(63));
	if (catinfo) {
	  std::cout << " abstol not reached and log-log convergence "
	    "returned errors.\n";
	}
      }
      else if (llc == 0) {
	list_assign(MLEobj_return, "convergence", wrap(11));
	if (catinfo) {
	  std::cout << "Warning! log-log convergence only. Maxit (="
		    << control_maxit << ") reached before abstol convergence.\n";
	}
      }
      else if (llc == 1) {
	list_assign(MLEobj_return, "convergence", wrap(1));
	if (catinfo) {
	  std::cout << " neither abstol nor log-log convergence tests "
	    "were passed.\n";
	}
      }
      else if (llc > 1) {
	list_assign(MLEobj_return, "convergence", wrap(72));
	if (catinfo) {
	  std::cout << " abstol not reached and log-log convergence "
	    "returned errors.\n";
	}
      }
    } // iter > as<int>(control["min.iter.conv.test"]) 
    else {
      list_assign(MLEobj_return, "convergence", wrap(12));
      if (catinfo) {
	std::cout << " abstol not reached and no log-log test info "
	  "since maxit less than min.iter.conv.test.\n";
      }
    }
  } // as<int>(as<List>(MLEobjiter["conv.test"])["convergence"]) == 4

  int mle_iter_conv = as<int>(as<List>(MLEobjiter["conv.test"])["convergence"]);
  if (mle_iter_conv == 3) {
    list_assign(MLEobj_return, "convergence", wrap(3));
    if (catinfo) {
      std::cout << "Warning! Abstol convergence only. no info on log-log "
	"convergence.\n Maxit (=" << control_maxit << ") < min.iter.conv.test (="
		<< as<int>(control["min.iter.conv.test"])
		<< ") so not log-log test.\n";
    }
  } // mle_iter_conv == 3
  else if (mle_iter_conv == 1) {
    list_assign(MLEobj_return, "convergence", wrap(10));
    if (catinfo) {
      std::cout << "Warning! Abstol convergence only. Maxit (="
		<< control_maxit << " reached before log-log convergence.\n";
    }
  } // mle_iter_conv == 1
  else if (mle_iter_conv == 0) {
    list_assign(MLEobj_return, "convergence", wrap(0));
    if (catinfo) {
      if (iter == as<int>(control["minit"])) {
	std::cout << "Success! algorithm run for " << iter
		  << " iterations. abstol and log-log tests passed.\n";
      }
      else {
	std::cout << "Success! abstol and log-log tests passed at " << iter
		  <<" iterations.\n";
      }
      if (as<double>(control["conv.test.slope.tol"]) > 0.1) {
	std::cout << "Alert: conv.test.slope.tol is "
		  << as<double>(control["conv.test.slope.tol"]) <<
	  ".\nTest with smaller values (<0.1) to ensure convergence.\n";
      }
    }
  } // mle_iter_conv == 0

  CharacterVector msg;
  if (as<List>(MLEobjiter["conv.test"]).containsElementNamed("messages") &&
      as<List>(MLEobjiter["conv.test"])["messages"] != R_NilValue) {
    CharacterVector ttcv =
      as<CharacterVector>(as<List>(MLEobjiter["conv.test"])["messages"]);
    if (ttcv.length() > 0) {
      msg.push_back("\nConvergence warnings\n");
      cv_append(msg, ttcv);
    }
  }
  list_assign(MLEobj_return, "par", MLEobjiter, "par");
  list_assign(MLEobj_return, "states", as<List>(MLEobjiter["kf"]), "xtT");
  list_assign(MLEobj_return, "logLik", MLEobjiter, "logLik");

  if (msg_kem.size()) {
    msg.push_back("\nMARSSkem warnings. Type MARSSinfo() for help.\n");
    std::copy(msg_kem.begin(), msg_kem.end(), std::back_inserter(msg));
  }

  if (msg_kf.size()) {
    msg.push_back("\nMARSSkf warnings. Type MARSSinfo() for help.\n");
    std::copy(msg_kf.begin(), msg_kf.end(), std::back_inserter(msg));
  }

  if (msg_kem.size() > 0 || msg_kf.size() > 0) {
    if (control_trace < 1) {
      msg.push_back("\nUse control$trace=1 to generate a more detailed "
		    "error report.\n");
    }
    if (control_silent == 2) {
      std::cout << "Alert: Numerical warnings were generated. Print the "
	"$errors element of output to see the warnings.\n";
    }
  }

  MLEobj_return.attr("class") = MLEobj.attr("class");
  if (msg.length() > 0) {
    list_assign(MLEobj_return, "errors", wrap(msg));
  }
  else {
    MLEobj_return = truenull(MLEobj_return, "errors");
  }
  MLEobj_return.attr("class") = MLEobj.attr("class");
  return MLEobj_return;
}

std::string Aupdate_impl(List& par1, const std::map<std::string, arma::cube>& star,
			 List& kf, List& Ey, List& timevarying,
			 List& MLEobjiter, int TT)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {
    List model = MLEobjiter["marss"];
    List f = model["fixed"];
    List d = model["free"];

    mat numer;
    mat denom;

    NumericMatrix Z = par1["Z"];
    mat mZ = as<mat>(Z);

    mat starR = star.at("R").slice(0);

    cube cdA = makecube(d["A"]);
    cube cfA = makecube(f["A"]);

    NumericVector dsubA = d["A"];
    List dsubAdimnames = dsubA.attr("dimnames");
	
    NumericVector fsubA = f["A"];
    List fsubAdimnames = fsubA.attr("dimnames");

    mat dA = as<mat>(sub3Dx(cdA, dsubAdimnames, 0));
    mat fA = as<mat>(sub3Dx(cfA, fsubAdimnames, 0));

    numer = mat(dA.n_cols, fA.n_cols, fill::zeros);
    denom = mat(dA.n_cols, dA.n_cols, fill::zeros);

    mat mEyytT = as<mat>(Ey["ytT"]);
    mat mkfxtT = as<mat>(kf["xtT"]);

    for (int i = 0; i < TT; i++) {
      if (timevarying["Z"]) {
	Z = as<NumericMatrix>(parmat(MLEobjiter, "Z", i)["Z"]);
	mZ = as<mat>(Z);
      }

      if (timevarying["A"]) {	// i == 0 case taken care of above
	dA = as<mat>(sub3Dx(cdA, dsubAdimnames, i));
	fA = as<mat>(sub3Dx(cfA, fsubAdimnames, i));
      }

      if (timevarying["R"]) {starR = star.at("R").slice(i);}

      mat x1 = dA.t() * starR *
	(mEyytT.cols(i, i) - mZ * mkfxtT.cols(i, i) - fA);
      mat x2 = dA.t() * starR * dA;
      numer += x1;
      denom += x2;
    }

    try {
      choleskyinverter( denom );
    }
    catch(std::runtime_error) {
      return "error in Aupdate.  denom is not invertible.";
    }

    as<List>(MLEobjiter["par"])["A"] = denom * numer;
    par1["A"] = parmat(MLEobjiter, "A", 0)["A"];

    return "";
  }
  catch(std::logic_error e) {
    return paster("logic_error in Aupdate - ", e.what(), "\n");
  }
  catch(std::runtime_error e) {
    return paster("runtime_error in Aupdate - ", e.what(), "\n");
  }
  catch(std::bad_alloc e) {
    return paster("bad_alloc in Aupdate - ", e.what(), "\n");
  }
  catch(exception e) {
    return paster("error in Aupdate - ", e.what(), "\narmadillo says: ",
		  armaerror.str());
  }
}


std::string Bupdate_impl(List par1, const std::map<std::string,cube>& star,
			 List kf, const std::string& kf_x0, const mat& IIzV0,
			 const mat& IImIIzV0, List& timevarying,
			 List MLEobjiter, int TT)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {
    List model = MLEobjiter["marss"];
    List f = model["fixed"];
    List d = model["free"];

    NumericVector dsubB = d["B"];
    List dsubBdimnames = dsubB.attr("dimnames");
	
    NumericVector fsubB = f["B"];
    List fsubBdimnames = fsubB.attr("dimnames");
	
    cube cdB = makecube(d["B"]);
    cube cfB = makecube(f["B"]);
    mat dB = as<mat>(sub3Dx(cdB, dsubBdimnames, 0));
    mat fB = as<mat>(sub3Dx(cfB, fsubBdimnames, 0));

    mat U = as<mat>(par1["U"]);
    mat x0 = as<mat>(par1["x0"]);
    mat V0 = as<mat>(par1["V0"]);
    
    mat mkfxtT = as<mat>(kf["xtT"]);
    mat mkfx0T = as<mat>(kf["x0T"]);
    mat mkfV0T = as<mat>(kf["V0T"]);
    cube ckfVtT = makecube(kf["VtT"]);
    cube ckfVtt1T = makecube(kf["Vtt1T"]);
    
    mat denom;
    mat numer;

    if ("x00" == kf_x0) {
      mat hatxtm = IImIIzV0 * mkfx0T + IIzV0 * x0;
      mat hatVtm = IImIIzV0 * mkfV0T * IImIIzV0 + IIzV0 * V0 * IIzV0;
      mat hatxt = mkfxtT.col(0);
      mat Ptm = hatVtm + hatxtm * mkfx0T.t();
      mat Pttm = ckfVtt1T.slice(0) + hatxt * mkfx0T.t();
      mat kronPtmQ = kron(Ptm, star.at("Q").slice(0));
	
      denom = dB.t() * kronPtmQ * dB;
      numer = dB.t() *
	(vectorise(star.at("Q").slice(0) * (Pttm - U * mkfx0T.t())) - kronPtmQ*fB);
    } else {
      denom = mat(1,1,fill::zeros);
      numer = mat(1,1,fill::zeros);
    }

    for (int i = 1; i < TT; ++i) {
      mat hatxtm = mkfxtT.col(i-1);
      mat hatVtm = ckfVtT.slice(i-1);
      mat hatxt = mkfxtT.col(i);
      mat Ptm = hatVtm + hatxtm * hatxtm.t();
      mat Pttm = ckfVtt1T.slice(i) + hatxt * hatxtm.t();

      const mat& starQ = timevarying["Q"] ?
	star.at("Q").slice(i) : star.at("Q").slice(0);

      mat kronPtmQ = kron(Ptm, starQ);

      if (timevarying["B"]) {
	dB = as<mat>(sub3Dx(cdB, dsubBdimnames, i));
	fB = as<mat>(sub3Dx(cfB, fsubBdimnames, i));
      }

      if (timevarying["U"]) {
	U = as<mat>(parmat(MLEobjiter, "U", i)["U"]);
      }

      denom = denom + dB.t() * kronPtmQ * dB;
      numer = numer + dB.t() *
	(vectorise(starQ * (Pttm - U * hatxtm.t())) - kronPtmQ * fB);
    }

    try {
      choleskyinverter( denom );
    }
    catch(std::runtime_error) {
      return std::string("error in Bupdate.  denom is not invertible.");
    }

    as<List>(MLEobjiter["par"])["B"] = denom * numer;
    par1["B"] = parmat(MLEobjiter, "B", 0)["B"];

    return std::string("");
  }
  catch(std::logic_error e) {
    return paster("logic_error in Bupdate - ", e.what(), "\n");
  }
  catch(std::runtime_error e) {
    return paster("runtime_error in Bupdate - ", e.what(), "\n");
  }
  catch(std::bad_alloc e) {
    return paster("bad_alloc in Bupdate - ", e.what(), "\n");
  }
  catch(exception e) {
    return paster("error in Bupdate - ", e.what(), "\narmadillo says: ",
		  armaerror.str());
  }
}

std::string Qupdate_impl_orig(List& par1, List& kf, List& Ey, List& timevarying,
			 List& MLEobjiter, const std::string& kf_x0,
			 const mat& IImIIzV0, const mat& IIzV0, int TT)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);
  try {
    List model = MLEobjiter["marss"];
    List f = model["fixed"];
    List d = model["free"];

    cube cdQ = makecube(d["Q"], true);
    List dQdimnames = (as<NumericVector>(d["Q"])).attr("dimnames");
    NumericMatrix dQ = sub3Dx(cdQ, dQdimnames, 0);
    NumericMatrix B = par1["B"];
    NumericMatrix U = par1["U"];

    NumericMatrix tdQdQ;
    NumericMatrix sum1;
    mat msum1;
    mat sum1a;

    mat mdQ = as<mat>(dQ);
    mat mtdQdQ(mdQ.n_cols, mdQ.n_cols, fill::zeros);

    cube ckfVtT = makecube(kf["VtT"]);
    cube ckfVtt1T = makecube(kf["Vtt1T"]);
    mat mkfxtT = as<mat>(kf["xtT"]);
    mat mkfV0T = as<mat>(kf["V0T"]);
    mat mkfx0T = as<mat>(kf["x0T"]);
    mat mpar1x0 = as<mat>(par1["x0"]);


    mat mB = as<mat>(par1["B"]);
    
    int TT_numer;
    if (kf_x0 == "x00") {
      TT_numer = TT;
      mat X0 = IImIIzV0 * mkfx0T + IIzV0 * mpar1x0;
      mat S00 = mkfV0T + X0 * X0.t();
      mat S10 = ckfVtt1T.slice(0) + mkfxtT.col(0) * X0.t();
      mat X1 = mkfxtT.col(0);
      mat S11 = ckfVtT.slice(0) + X1 * X1.t();

      mat mB = as<mat>(B);
      mat mU = as<mat>(U);

      mat mBX0 = mB * X0;

      mat tmp1(mB);
      tmp1 = tmp1 * S10.t();
      tmp1 = S10 * mB.t();
      tmp1 = mB * S00;
      tmp1 = tmp1 * mB.t();
      tmp1 = mU * X1.t();
      tmp1 = X1 * mU.t();
      tmp1 = mU * mBX0.t();
      tmp1 = mBX0 * mU.t();
      tmp1 = mU * mU.t();
      
      sum1a = S11;
      sum1a = sum1a - mB * S10.t();
      sum1a = sum1a - S10 * mB.t();
      sum1a = sum1a + ((mB * S00) * mB.t());
      sum1a = sum1a - mU * X1.t();
      sum1a = sum1a - X1 * mU.t();
      sum1a = sum1a + mU * mBX0.t();
      sum1a = sum1a + mBX0 * mU.t();
      sum1a = sum1a+ mU * mU.t();

      sum1a = S11 - mB * S10.t() - S10 * mB.t() + ((mB * S00) * mB.t())
	- mU * X1.t() - X1 * mU.t()
	+ mU * mBX0.t() + mBX0 * mU.t() + mU * mU.t();

      sum1a = (sum1a + sum1a.t()) / 2;

      msum1 = mdQ.t() * vectorise(sum1a);

      sum1 = wrap<double>(msum1);
      mtdQdQ = mdQ.t() * mdQ;
    }
    else if (kf_x0 == "x10") {
      sum1 = 0;
      TT_numer = TT - 1;
    }

    mat minvdQ;
    IntegerVector aargt(1);
    for (int i = 1; i < TT; ++i) {
      mat X0 = mkfxtT.col(i - 1);
      mat X1 = mkfxtT.col(i);
      mat S00 = ckfVtT.slice(i - 1) + X0 * X0.t();
      mat S10 = ckfVtt1T.slice(i) + X1 * X0.t();
      mat S11 = ckfVtT.slice(i) + X1 * X1.t();

      aargt[0] = i;
	    
      if (timevarying["B"]) {
	B = as<NumericMatrix>(parmat(MLEobjiter, "B", aargt)["B"]);
      }
	    
      if (timevarying["U"]) {
	U = as<NumericMatrix>(parmat(MLEobjiter, "U", aargt)["U"]);
      }
	    
      if (timevarying["Q"]) {
	dQ = sub3Dx(cdQ, dQdimnames, i);
	mdQ = as<mat>(dQ);
	mtdQdQ = mtdQdQ + mdQ.t() * mdQ;
      }

      mat mB = as<mat>(B);
      mat mU = as<mat>(U);
      mat mBX0 = mB * X0;

      mat sum1a = S11 - mB * S10.t() - S10 * mB.t()
	+ (mB * S00) * mB.t()
	- mU * X1.t() - X1 * mU.t()
	+ mU * mBX0.t() + mBX0 * mU.t()
	+ mU * mU.t();
      sum1a = (sum1a + sum1a.t()) / 2;

      msum1 += mdQ.t() * vectorise(sum1a);
    }


    if (timevarying["Q"]) {
      try {
	minvdQ = choleskyinverter(mtdQdQ);
      }
      catch(std::runtime_error) {
	return String ("error in Qupdate. "
		       "For time-varying Q, mtdQdQ is not invertible.");
      }
    } else {
      mdQ = as<mat>(dQ);
      mtdQdQ = mdQ.t() * mdQ;
      try {
	minvdQ = choleskyinverter(mtdQdQ) / TT_numer;
      }
      catch(std::runtime_error) {
	return String ("error in Qupdate. "
		       "For non-time-varying Q, mtdQdQ is not invertible.");
      }
    }

    as<List>(MLEobjiter["par"])["Q"] = NumericMatrix(wrap(minvdQ * msum1));
    aargt[0] = 0;
    par1["Q"] = parmat(MLEobjiter, "Q", aargt)["Q"];

    return String("");
  }
  catch(std::logic_error e) {
    return paster("logic_error in Qupdate - ", e.what(), "\n");
  }
  catch(std::runtime_error e) {
    return paster("runtime_error in Qupdate - ", e.what(), "\n");
  }
  catch(std::bad_alloc e) {
    return paster("bad_alloc in Qupdate - ", e.what(), "\n");
  }
  catch(exception e) {
    return paster("error in Qupdate - ", e.what(), "\narmadillo says: ",
		  armaerror.str());
  }
}

std::string Qupdate_impl(List& par1, List& kf, List& Ey, List& timevarying,
			 List& MLEobjiter, const std::string& kf_x0,
			 const mat& IImIIzV0, const mat& IIzV0, int TT)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);
  try {
    List model = MLEobjiter["marss"];
    List f = model["fixed"];
    List d = model["free"];

    cube cdQ = makecube(d["Q"], true);
    List dQdimnames = (as<NumericVector>(d["Q"])).attr("dimnames");
    NumericMatrix dQ = sub3Dx(cdQ, dQdimnames, 0);

    NumericMatrix tdQdQ;

    mat msum1;
    mat sum1a;

    mat mdQ = as<mat>(dQ);
    mat mtdQdQ(mdQ.n_cols, mdQ.n_cols, fill::zeros);

    cube ckfVtT = makecube(kf["VtT"]);
    cube ckfVtt1T = makecube(kf["Vtt1T"]);
    mat mkfxtT = as<mat>(kf["xtT"]);
    mat mkfV0T = as<mat>(kf["V0T"]);
    mat mkfx0T = as<mat>(kf["x0T"]);
    mat mpar1x0 = as<mat>(par1["x0"]);

    mat mBpar1 = as<mat>(par1["B"]);
    mat mUpar1 = as<mat>(par1["U"]);

    int TT_numer;
    if (kf_x0 == "x00") {
      TT_numer = TT;
      mat X0 = IImIIzV0 * mkfx0T + IIzV0 * mpar1x0;
      mat S00 = mkfV0T + X0 * X0.t();
      mat S10 = ckfVtt1T.slice(0) + mkfxtT.col(0) * X0.t();
      mat X1 = mkfxtT.col(0);
      mat S11 = ckfVtT.slice(0) + X1 * X1.t();

      mat mB = mBpar1;
      mat mU = mUpar1;
      
      mat mBX0 = mB * X0;

      mat tmp1(mB);
      tmp1 = tmp1 * S10.t();
      tmp1 = S10 * mB.t();
      tmp1 = mB * S00;
      tmp1 = tmp1 * mB.t();
      tmp1 = mU * X1.t();
      tmp1 = X1 * mU.t();
      tmp1 = mU * mBX0.t();
      tmp1 = mBX0 * mU.t();
      tmp1 = mU * mU.t();
      
      sum1a = S11;
      sum1a = sum1a - mB * S10.t();
      sum1a = sum1a - S10 * mB.t();
      sum1a = sum1a + ((mB * S00) * mB.t());
      sum1a = sum1a - mU * X1.t();
      sum1a = sum1a - X1 * mU.t();
      sum1a = sum1a + mU * mBX0.t();
      sum1a = sum1a + mBX0 * mU.t();
      sum1a = sum1a+ mU * mU.t();

      sum1a = S11 - mB * S10.t() - S10 * mB.t() + ((mB * S00) * mB.t())
	- mU * X1.t() - X1 * mU.t()
	+ mU * mBX0.t() + mBX0 * mU.t() + mU * mU.t();

      sum1a = (sum1a + sum1a.t()) / 2;

      msum1 = mdQ.t() * vectorise(sum1a);
      mtdQdQ = mdQ.t() * mdQ;
    }
    else if (kf_x0 == "x10") {
      TT_numer = TT - 1;
    }

    mat minvdQ;
    IntegerVector aargt(1);
    mat mB(mBpar1);
    mat mU(mUpar1);

    IntegerVector iis(TT);
    for (int iiis = 0; iiis < TT; ++iiis) {
      iis[iiis] = iiis;
    }
    
    cube Bcube;
    if (timevarying["B"]) {
      Bcube = parmat_cube(MLEobjiter, "B", iis)["B"];
    }

    cube Ucube;
    if (timevarying["U"]) {
      Ucube = parmat_cube(MLEobjiter, "U", iis)["U"];
    }
    
    for (int i = 1; i < TT; ++i) {
      mat X0 = mkfxtT.col(i - 1);
      mat X1 = mkfxtT.col(i);
      mat S00 = ckfVtT.slice(i - 1) + X0 * X0.t();
      mat S10 = ckfVtt1T.slice(i) + X1 * X0.t();
      mat S11 = ckfVtT.slice(i) + X1 * X1.t();

      aargt[0] = i;
	    
      if (timevarying["B"]) {
	mB = Bcube.slice(i);
      }
	    
      if (timevarying["U"]) {
	mU = Ucube.slice(i);
      }
	    
      if (timevarying["Q"]) {
	dQ = sub3Dx(cdQ, dQdimnames, i);
	mdQ = as<mat>(dQ);
	mtdQdQ = mtdQdQ + mdQ.t() * mdQ;
      }

      mat mBX0 = mB * X0;

      mat sum1a = S11 - mB * S10.t() - S10 * mB.t()
	+ (mB * S00) * mB.t()
	- mU * X1.t() - X1 * mU.t()
	+ mU * mBX0.t() + mBX0 * mU.t()
	+ mU * mU.t();
      sum1a = (sum1a + sum1a.t()) / 2;

      msum1 += mdQ.t() * vectorise(sum1a);
    }


    if (timevarying["Q"]) {
      try {
	minvdQ = choleskyinverter(mtdQdQ);
      }
      catch(std::runtime_error) {
	return String ("error in Qupdate. "
		       "For time-varying Q, mtdQdQ is not invertible.");
      }
    } else {
      mdQ = as<mat>(dQ);
      mtdQdQ = mdQ.t() * mdQ;
      try {
	minvdQ = choleskyinverter(mtdQdQ) / TT_numer;
      }
      catch(std::runtime_error) {
	return String ("error in Qupdate. "
		       "For non-time-varying Q, mtdQdQ is not invertible.");
      }
    }

    as<List>(MLEobjiter["par"])["Q"] = NumericMatrix(wrap(minvdQ * msum1));
    aargt[0] = 0;
    par1["Q"] = parmat(MLEobjiter, "Q", aargt)["Q"];

    return String("");
  }
  catch(std::logic_error e) {
    return paster("logic_error in Qupdate - ", e.what(), "\n");
  }
  catch(std::runtime_error e) {
    return paster("runtime_error in Qupdate - ", e.what(), "\n");
  }
  catch(std::bad_alloc e) {
    return paster("bad_alloc in Qupdate - ", e.what(), "\n");
  }
  catch(exception e) {
    return paster("error in Qupdate - ", e.what(), "\narmadillo says: ",
		  armaerror.str());
  }
}

std::string Rupdate_impl(List& par1, List& kf, List& Ey,
		    List& time_varying, List& MLEobjiter, List& d, int TT)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {
    int dRcols = as<IntegerVector>(as<NumericVector>(d["R"]).attr("dim"))[1];
    mat sum1 = zeros(dRcols, 1);

    int Rcols = as<IntegerVector>(as<NumericVector>(d["R"]).attr("dim"))[1];
    mat tdRdR = zeros(Rcols, Rcols);

    NumericVector Z = par1["Z"];
    NumericVector A = par1["A"];

    IntegerVector zdims = Z.attr("dim");
    IntegerVector adims = A.attr("dim");
    cube zcube(zdims[0], zdims[1], TT);
    cube acube(adims[0], adims[1], TT);

    IntegerVector iis = IdxVec::get(TT);
      
    if (time_varying["Z"]) {
      zcube = as<cube>(parmat(MLEobjiter, "Z", iis)["Z"]);
    }
    if (time_varying["A"]) {
      acube = as<cube>(parmat(MLEobjiter, "A", iis)["A"]);
    }

    mat Ey_ytT = as<mat>(Ey["ytT"]);
    cube Ey_yxtT = makecube(Ey["yxtT"]);
    cube Ey_OtT = makecube(Ey["OtT"]);
    List Ey_yxtT_dn = as<List>(Ey["yxtT"]).attr("dimnames");
    List Ey_OtT_dn = as<List>(Ey["OtT"]).attr("dimnames");
    cube kf_VtT = makecube(kf["VtT"]);

    mat hatxt_mat = as<mat>(kf["xtT"]);

    mat mZ = as<mat>(Z);
    mat mA = as<mat>(A);

    mat dR;

    List Rcolnames;
    for (int i = 0; i < TT; ++i) {
      if (time_varying["Z"] && i > 0) {
	mZ = zcube.slice(i);
      }
      if (time_varying["A"] && i > 0) {
	mA = acube.slice(i);
      }

      if (time_varying["R"] || i == 0) {
	dR = as<mat>(sub3D(d["R"], i));
	tdRdR += dR.t() * dR;
	Rcolnames = as<NumericVector>(d["R"]).attr("dimnames");
      }

      mat hatyt = Ey_ytT.col(i);
      mat hatyxt = as<mat>(sub3Dx(Ey_yxtT, Ey_yxtT_dn,i));
      mat hatOt = as<mat>(sub3Dx(Ey_OtT, Ey_OtT_dn,i));
	
      mat hatxt = hatxt_mat.col(i);
      mat hatPt = kf_VtT.slice(i) + hatxt * hatxt.t();

      mat sum1a = hatOt
	- hatyxt * mZ.t() - mZ * hatyxt.t()
	- hatyt * mA.t() - mA * hatyt.t()
	+ (mZ * hatPt) * mZ.t()
	+ (mZ * hatxt) * mA.t()
	+ mA * (mZ * hatxt).t()
	+ mA * mA.t();
      sum1a = (sum1a + sum1a.t()) / 2.;
      sum1 += dR.t() * vectorise(sum1a);
    }

    mat invdR = pcholinv(tdRdR);

    if (!time_varying["R"]) {
      invdR /= TT;
    }

    as<List>(MLEobjiter["par"])["R"] = invdR * sum1;
      
    NumericVector tmpR = as<List>(MLEobjiter["par"])["R"];
    if (tdRdR.n_cols == 1) {
      tmpR.attr("dimnames") = List::create(Rcolnames[1],
    					 StringVector());
    }
    par1["R"] = parmat(MLEobjiter,"R",0)["R"];

    return("");
  }
  
  catch(std::logic_error e) {
    return paster("logic_error in Rupdate - ", e.what(), "\n");
  }
  catch(std::runtime_error e) {
    return paster("runtime_error in Rupdate - ", e.what(), "\n");
  }
  catch(std::bad_alloc e) {
    return paster("bad_alloc in Rupdate - ", e.what(), "\n");
  }
  catch(exception e) {
    return paster("error in Rupdate - ", e.what(), "\narmadillo says: ",
		  armaerror.str());
  }
}

std::string Uupdate_impl(int m, List par1,
			 const std::map<std::string,cube>& star,
			 const std::map<std::string, cube>& IIz,
			 const mat& IIzV0, const mat& IImIIzV0,
			 const mat& IIm, List Ey, List kf,
			 const std::string& kf_x0, List& timevarying,
			 List MLEobjiter, int TT)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {
    List model = MLEobjiter["marss"];
    List f = model["fixed"];
    List d = model["free"];

    cube cfU = makecube(f["U"]);
    NumericVector tmpfU = f["U"];
    List cfUdimnames = tmpfU.attr("dimnames");
    cube cdU = makecube(d["U"]);
    NumericVector tmpdU = d["U"];
    List cdUdimnames = tmpdU.attr("dimnames");
    mat mfU = as<mat>(sub3Dx(cfU, cfUdimnames, 0));
    mat mdU = as<mat>(sub3Dx(cdU, cdUdimnames, 0));

    NumericMatrix B = par1["B"];
    NumericMatrix Z = par1["Z"];
    NumericMatrix A = par1["A"];
    mat mB = as<mat>(B);
    mat mZ = as<mat>(Z);
    mat mA = as<mat>(A);

    mat mEyytT = as<mat>(Ey["ytT"]);

    mat Qinv = star.at("Q").slice(0);

    mat Rinv = star.at("R").slice(0);
	
    mat numer = mat(mdU.n_cols, 1, fill::zeros);
    mat denom = mat(mdU.n_cols, mdU.n_cols, fill::zeros);

    mat mIIzQ = IIz.at("Q").slice(0);
    vec diagQ = 1 - diagvec(mIIzQ);

    mat diagV0 = 1 - diagvec(IIzV0);

    int nQ0 = diagQ.n_elem - vec(nonzeros(diagQ)).n_elem;

    mat hatxt0;
    if (kf_x0 == "x00") {
      hatxt0 = (as<mat>(kf["x0T"]));
    } else {
      hatxt0 = mat(as<mat>(kf["xtT"])).head_cols(1);
    }

    mat E_x0 = IImIIzV0 * hatxt0 + IIzV0 * as<mat>(par1["x0"]);

    mat AdjM(mB);
    AdjM.elem(find(AdjM != 0)).ones();

    mat mkfxtT = as<mat>(kf["xtT"]);
    mat Bstar;
    mat Bstar_tm;
    mat fstar;
    mat fstar_tm;
    mat Dstar;
    mat Dstar_tm;
    mat Mt;
    mat IId_tm;
    mat IId;
    mat Delta3;
    mat Delta4;


    if (kf_x0 == "x00") {
      Bstar = mB;
      Bstar_tm = IIm;
      fstar = mfU;
      fstar_tm = mat(mfU.n_rows, mfU.n_cols, fill::zeros);
      Dstar = mdU;
      Dstar_tm = mat(mdU.n_rows, mdU.n_cols, fill::zeros);
      Mt = AdjM;
      IId_tm = IIm;
      IId = diagmat(1 - diagQ);

      if (any(diagQ) && !all(diagQ)) {
	mat tMt = trimmat(Mt, diagQ, false, true);
	vec zerorows(tMt.n_rows);
	for (unsigned int i = 0; i < tMt.n_rows; ++i) {
	  zerorows[i] = all(tMt.row(i) == 0) ? 1 : 0;
	}
	uvec therows = find(diagQ == 0);
	mat subIId(IId.submat(therows, therows));
	subIId.diag() = zerorows;
	IId.submat(therows, therows) = subIId;
      }

      Delta3 = mkfxtT.head_cols(1) - mB * E_x0 - mfU;
      Delta4 = mdU;

      numer += Delta4.t() * (Qinv * Delta3);
      denom += Delta4.t() * (Qinv * Delta4);
    } else {
      Bstar = IIm;
      fstar = mat(mfU.n_rows, mfU.n_cols, fill::zeros);
      Dstar = mat(mdU.n_rows, mdU.n_cols, fill::zeros);
      IId_tm = IIm;
      IId_tm.zeros();
      Mt = IId = IIm;
    }

    if (any(vectorise(IId) == 1)) {
      mat Delta1 = mEyytT.head_cols(1) -
	mZ * (IIm - IId) * mkfxtT.head_cols(1) -
	mZ * (IId * (Bstar * E_x0 + fstar)) - mA;
      mat Delta2 = mZ * IId * Dstar;
      numer += Delta2.t() * (Rinv * Delta1);
      denom += Delta2.t() * (Rinv * Delta2);
    }

    IntegerVector iis(TT);
    for (int iiis = 0; iiis < TT; ++iiis) {
      iis[iiis] = iiis;
    }

    cube Bcube;
    if (timevarying["B"]) {
      Bcube = parmat_cube(MLEobjiter, "B", iis)["B"];
    }
    
    cube Acube;
    if (timevarying["A"]) {
      Acube = parmat_cube(MLEobjiter, "A", iis)["A"];
    }
    
    for (int t = 1; t < TT; ++t) {
      if (timevarying["U"]) {
	mfU = as<mat>(sub3Dx(cfU, cfUdimnames, t));
	mdU = as<mat>(sub3Dx(cdU, cdUdimnames, t));
      }
      if (timevarying["B"]) {
	mB = Bcube.slice(t);
      }
      if (timevarying["A"]) {
	mA = Acube.slice(t);
      }
      if (timevarying["Z"]) {
	Z = as<NumericMatrix> (parmat(MLEobjiter, "Z", t)["Z"]);
	mZ = as<mat>(Z);
      }
      if (timevarying["R"]) { Rinv = star.at("R").slice(t); }
      if (timevarying["Q"]) { Qinv = star.at("Q").slice(t); }

      fstar_tm = fstar;
      fstar = mB * fstar + mfU;
      Dstar_tm = Dstar;
      Dstar = mB * Dstar + mdU;
      Bstar_tm = Bstar;
      Bstar = mB * Bstar;

      if (t <= m) {
	IId_tm = IId;
	IId = diagmat(1 - diagQ);
	if (any(diagQ) && !all(diagQ)) {
	  Mt = AdjM * Mt;
	  unzerorows(IId, Mt, diagQ);
	}
      }

      if (any(vectorise(IId) == 1)) {
	mat Delta1 = mEyytT.cols(t, t) -
	  mZ * ((IIm - IId) * mkfxtT.cols(t, t)) -
	  mZ * (IId * (Bstar * E_x0 + fstar)) - mA;
	mat Delta2 = mZ * IId * Dstar;
	numer += Delta2.t() * (Rinv * Delta1);
	denom += Delta2.t() * (Rinv * Delta2);
      }

      Delta3 = mkfxtT.cols(t, t) -
	mB * ((IIm - IId_tm) * mkfxtT.cols(t - 1, t - 1)) -
	mB * (IId_tm * (Bstar * E_x0 + fstar)) - mfU;
      Delta4 = mdU + mB * IId_tm * Dstar_tm;
      numer += Delta4.t() * Qinv * Delta3;
      denom += Delta4.t() * Qinv * Delta4;
    }

    try {
      choleskyinverter( denom );
    }
    catch(std::runtime_error) {
      return "error in Uupdate.  denom is not invertible.";
    }

    as<List>(MLEobjiter["par"])["U"] = denom * numer;
    par1["U"] = parmat(MLEobjiter, "U", 0)["U"];

    return "";
  }
  catch(std::logic_error e) {
    return paster("logic_error in Uupdate - ", e.what(), "\n");
  }
  catch(std::runtime_error e) {
    return paster("runtime_error in Uupdate - ", e.what(), "\n");
  }
  catch(std::bad_alloc e) {
    return paster("bad_alloc in Uupdate - ", e.what(), "\n");
  }
  catch(exception e) {
    return paster("error in Uupdate - ", e.what(), "\narmadillo says: ",
		  armaerror.str());
  }
}

String V0update_impl(List& par1, List& kf, List& MLEobjiter)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {

    List model = MLEobjiter["marss"];
    List f = model["fixed"];
    List d = model["free"];
	
    mat dV0 = as<mat>(sub3D(d["V0"], 0));
    mat kfV0T = as<mat>(kf["V0T"]);
    mat tmp = dV0.t() * (dV0.t(), vec(kfV0T));
    mat V0update = choleskyinverter(tmp);
    as<List>(MLEobjiter["par"])["V0"] = V0update;
    par1["V0"] = parmat(MLEobjiter, "V0", 0)["V0"];
	
    return String("");
  }
  catch(std::logic_error e) {
    return String(paster("logic_error in V0update - ",e.what(),"\n"));
  }
  catch(std::runtime_error e) {
    return String(paster("runtime_error in V0update - ",e.what(),"\n"));
  }
  catch(std::bad_alloc e) {
    return String(paster("bad_alloc in V0update - ",e.what(),"\n"));
  }
  catch(exception e) {
    return String(paster("error in V0update - ",e.what(),"\n",
			 "armadillo says: ", armaerror.str()));
  }
}


String x0update_impl(List& par1, std::map<std::string,cube>& star,
		     List& kf, List& Ey, List& timevarying,
		     List& MLEobjiter, const std::string& kf_x0,
		     std::map<std::string, cube>& IIz,
		     const mat& mIIm, const mat& mIImIIzV0,
		     const mat& mIIzV0, int iter, int TT)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {
    List model = MLEobjiter["marss"];
    List f = model["fixed"];
    List d = model["free"];

    List model_dims = model.attr("model.dims");
    IntegerVector model_dims_x = model_dims["x"];
    int m = model_dims_x[0];

    mat fx0 = as<mat>(sub3D(f["x0"], 0));
    mat dx0 = as<mat>(sub3D(d["x0"], 0));

    mat A = as<mat>(par1["A"]);
    mat Z = as<mat>(par1["Z"]);
    mat B = as<mat>(par1["B"]);
    mat U = as<mat>(par1["U"]);

    mat mkfx0T = as<mat>(kf["x0T"]);
    mat mkfxtT = as<mat>(kf["xtT"]);

    mat EyytT = as<mat>(Ey["ytT"]);
    mat Qinv = as<mat>(sub3Dx(star["Q"], List(0), 0));

    const mat mIIzQ = IIz["Q"].slice(0);
    vec diagQ = 1 - diagvec(mIIzQ);
    int nQ0 = uvec(find(diagQ == 0)).n_elem;

    mat Rinv = as<mat>(sub3Dx(star["R"], List(0), 0));
    mat mIIzR = IIz["R"].slice(0);
    vec diagR1 = 1 - diagvec(mIIzR);

    bool x0degenupdate = false;

    mat diagV0 = 1 - diagvec(mIIzV0);

    IntegerVector aargt(1);

    if (any(diagQ == 0)) {
      mat tmpdiag = trimmat(dx0, diagQ, true, false);
      NumericVector tmpdiagvec(tmpdiag.begin(), tmpdiag.end());
      IntegerVector tmpdims =
	IntegerVector::create(tmpdiag.n_rows, tmpdiag.n_cols);
      tmpdiagvec.attr("dims") = tmpdims;
      try {
	x0degenupdate = isfixed(tmpdiagvec);
      }
      catch(std::invalid_argument e) {
	return String(e.what());
      }
    }
    mat denom(dx0.n_cols, dx0.n_cols, fill::zeros);
    mat numer(dx0.n_cols, 1, fill::zeros);

    mat hatxt0;
    if (kf_x0 == "x00") {
      hatxt0 = mkfx0T;
    } else {
      hatxt0 = mkfxtT.col(0);
    }

    if (any(vectorise(diagV0) == 1)) {
      denom = dx0.t() * (star["V0"].slice(0) * dx0);
      numer = dx0.t() * (hatxt0 - fx0);
    }

    if (!all(vectorise(diagV0))) {
      mat AdjM = B;
      AdjM.elem(find(AdjM != 0)).ones();

      mat Bstar;
      mat Bstar_tm;
      mat Ustar;
      mat Ustar_tm;
      mat IId;
      mat IId_tm;
      mat Mt;
      mat Delta5;
      mat Delta6;
      mat Delta7;
      mat Delta8;

      if (kf_x0 == "x00") {
	Bstar = B;
	Bstar_tm = mIIm;
	Ustar = U;
	Ustar_tm = U;
	Ustar_tm.zeros();
	Mt = AdjM;
	IId_tm = mIIm;
	IId = diagmat(1 - diagQ);

	if (!all(diagQ) && any(diagQ)) {
	  mat tMt = trimmat(Mt, diagQ, false, true);
	  vec zerorows(tMt.n_rows);
	  for (unsigned int i = 0; i < tMt.n_rows; ++i) {
	    zerorows[i] = all(tMt.row(i) == 0) ? 1 : 0;
	  }
	  uvec therows = find(diagQ == 0);
	  mat subIId(IId.submat(therows, therows));
	  subIId.diag() = zerorows;
	  IId.submat(therows, therows) = subIId;
	}

	Delta7 = mkfxtT.head_cols(1) -
	  B * (mIImIIzV0 * hatxt0 + mIIzV0 * fx0) - U;
	Delta8 = B * mIIzV0 * dx0;
	numer = numer + Delta8.t() * (Qinv * Delta7);
	denom = denom + Delta8.t() * (Qinv * Delta8);
      } else {
	Bstar = mIIm;
	Ustar = U;
	Ustar.zeros();
	IId_tm = mat(1,1,fill::zeros);
	IId = mIIm;
	Mt = mIIm;
      }
	    
      if (any(vectorise(IId) == 1)) {
	Delta5 = EyytT.head_cols(1) -
	  Z * (mIIm - IId) * mkfxtT.head_cols(1) -
	  Z * IId * (Bstar * (mIImIIzV0*hatxt0+mIIzV0*fx0)+Ustar) - A;
	Delta6 = Z * IId * Bstar * mIIzV0 * dx0;
	if (!all(vectorise(diagR1))) {
	  if (any(vectorise((Delta6.t() * mIIzR) * Delta6))) {
	    return String(paster("Stopped at iter ",iter," in MARSSkem at x0 "
				 "update.\nThere are 0s on R diagonal. x0 assoc "
				 "with these must be fixed (not estimated).\n"));
	  }
	}
	numer = numer + Delta6.t() * (Rinv * Delta5);
	denom = denom + Delta6.t() * (Rinv * Delta6);
      }

      IntegerVector iis(TT);
      for (int iiis = 0; iiis < TT; ++iiis) {
	iis[iiis] = iiis;
      }
    
      cube Acube;
      if (timevarying["A"]) {
	Acube = parmat_cube(MLEobjiter, "A", iis)["A"];
      }

      cube Bcube;
      if (timevarying["B"]) {
	Bcube = parmat_cube(MLEobjiter, "B", iis)["B"];
      }

      cube Ucube;
      if (timevarying["U"]) {
	Ucube = parmat_cube(MLEobjiter, "U", iis)["U"];
      }

      for (int t = 1; t < TT; ++t) {
	aargt[0] = t;
	if (timevarying["A"]) {
	  A = Acube.slice(t);
	}
	if (timevarying["B"]) {
	  B = Bcube.slice(t);
	}
	if (timevarying["U"]) {
	  U = Ucube.slice(t);
	}
	if (timevarying["Z"]) {
	  Z = as<mat>(parmat(MLEobjiter, "Z", aargt)["Z"]);
	}
	if (timevarying["R"]) {
	  Rinv = star["R"].slice(t);
	  mIIzR = as<mat>(sub3Dx(IIz["R"],List(0), t));
	}
	if (timevarying["Q"]) {
	  Qinv = star["Q"].slice(t);
	}

	Ustar_tm = Ustar;
	Ustar = B * Ustar + U;
	Bstar_tm = Bstar;
	Bstar = B * Bstar;

	if (t <= m) {
	  IId_tm = IId;
	  Mt = AdjM * Mt;
	  IId = diagmat(1 - diagQ);
	  if (!all(diagQ) && any(diagQ)) {
	    mat tMt = trimmat(Mt, diagQ, false, true);
	    vec zerorows(tMt.n_rows);
	    for (unsigned int i = 0; i < tMt.n_rows; ++i) {
	      zerorows[i] = all(tMt.row(i) == 0) ? 1 : 0;
	    }
	    uvec therows = find(diagQ == 0);
	    mat subIId(IId.submat(therows, therows));
	    subIId.diag() = zerorows;
	    IId.submat(therows, therows) = subIId;
	  }
	}

	if (any(vectorise(IId) == 1)) {
	  Delta5 = EyytT.cols(t,t) -
	    Z * ((mIIm - IId) * mkfxtT.cols(t,t)) -
	    Z * IId *
	    (Bstar * (mIImIIzV0 * hatxt0 + mIIzV0 * fx0) + Ustar) -
	    A;
	  Delta6 = Z * IId * Bstar * (mIIzV0 * dx0);

	  if (!all(vectorise(diagR1))) {
	    if (any(vectorise((Delta6.t() * mIIzR) * Delta6))) {
	    return String(paster("Stopped at iter ",iter," in MARSSkem at x0 "
				 "update.\nThere are 0s on R diagonal. x0 assoc "
				 "with these must be fixed (not estimated).\n"));
	    }
	  }
	  numer = numer + Delta6.t() * Rinv * Delta5;
	  denom = denom + Delta6.t() * Rinv * Delta6;
	}

	if (any(vectorise(IId_tm) == 1)) {
	  Delta7 = mkfxtT.cols(t,t) -
	    B * (mIIm - IId_tm) * mkfxtT.cols(t-1,t-1) -
	    B * IId_tm *
	    (Bstar_tm * (mIImIIzV0 * hatxt0 + mIIzV0 * fx0) + Ustar) -
	    U;
	  Delta8 = B * IId_tm * Bstar_tm * mIIzV0 * dx0;
	  numer = numer + Delta8.t() * Qinv * Delta7;
	  denom = denom + Delta8.t() * Qinv * Delta8;
	}
      }
    }

    try {
      choleskyinverter( denom );
    }
    catch(std::runtime_error) {
      return String ("error in x0update.  denom is not invertible.");
    }

    as<List>(MLEobjiter["par"])["x0"] = denom * numer;

    aargt[0] = 0;
    par1["x0"] = parmat(MLEobjiter, "x0", aargt)["x0"];

    return String("");

  }
  catch(std::logic_error e) {
    return String(paster("logic_error in x0update - ",e.what(),"\n"));
  }
  catch(std::runtime_error e) {
    return String(paster("runtime_error in x0update - ",e.what(),"\n"));
  }
  catch(std::bad_alloc e) {
    return String(paster("bad_alloc in x0update - ",e.what(),"\n"));
  }
  catch(exception e) {
    return String(paster("error in x0update - ",e.what(),"\n",
			 "armadillo says: ", armaerror.str()));
  }
}


/** utilities **/

List loglog_conv_test(List& iter_record_in, int iter, int deltaT, double tol,
		      const CharacterVector& params_to_test)
{
  List iter_record = clone(iter_record_in);
  bool improper = false;
  if (iter_record.attr("names") == R_NilValue) {
    improper = true;
  }

  CharacterVector itrnames = iter_record.attr("names");
  if (!improper) {
    for (String name : CharacterVector{"par", "logLik"}) {
      if (std::find(itrnames.begin(), itrnames.end(), name)==itrnames.end()) {
	improper = true;
	break;
      }
    }
  }

  if (!improper) {
    CharacterVector itrparnames =
      (as<List>(iter_record["par"]).attr("names") == R_NilValue) ?
      CharacterVector() : as<List>(iter_record["par"]).attr("names");

    std::copy(itrparnames.begin(), itrparnames.end(), std::back_inserter(itrnames));
    
    bool found = false;
    for (String name : params_to_test) {
      if (std::find(itrnames.begin(),itrnames.end(),name)!=itrnames.end()){
	found = true;
	break;
      }
    }

    improper = !found;
  }

  NumericMatrix tmppar = iter_record["par"];
  if (!improper) {
    improper = (tmppar.attr("dim") == R_NilValue) ||
      (as<IntegerVector>(tmppar.attr("dim")).length() !=  2);
  }

  if (!improper) {
    improper = as<IntegerVector>(tmppar.attr("dim"))[0]<=1;
  }

  if (!improper) {
    improper = tmppar.attr("dimnames") == R_NilValue;
    int x = 2;
    int y = x * 8;
  }
  if (!improper) {
    List itrpardimnames = tmppar.attr("dimnames");
    if (itrpardimnames.length() < 2) {
      improper = true;
    }
    else if (itrpardimnames[1] == R_NilValue) {
      improper = true;
    }
  }
  
  if (improper) {
    return List::create(
      _["convergence"] = -1,
      _["messages"] = CharacterVector(
	"par list not a proper list (with par and logLik) "
	"or too short for conv test or has no column names.\n")
      );
  }

  NumericVector iter_record_par;  // will be a matrix

  if (std::find(params_to_test.begin(), params_to_test.end(), "logLik")
      != params_to_test.end()) {
    vec irll = as<vec>(iter_record["logLik"]);
    vec expll = exp(irll - mean(irll));
    NumericVector nvll(expll.begin(), expll.end());
    iter_record_par = scbind(iter_record["par"], nvll,
			     R_NilValue, CharacterVector("logLik"));
  }
  else {
    iter_record_par = as<NumericVector>(iter_record["par"]);
  }

  List irpar_colnames = as<List>(iter_record_par.attr("dimnames"));
  CharacterVector names_iter = irpar_colnames[1];
  std::vector<std::string> p_elems;
  p_elems.resize(names_iter.length());
  std::transform(names_iter.begin(), names_iter.end(), p_elems.begin(),
		 [](String s)->std::string {
		   std::string ss = s; return ss.substr(0,ss.find('.'));});

  std::vector<double> test_conv(names_iter.length(), 0.);
  
  int test_len2 = as<IntegerVector>(iter_record_par.attr("dim"))[0];
  int test_len1 = std::max(1, test_len2-deltaT);

  int mtl2d =std::min(test_len2-1, deltaT);
  std::vector<int> test_len;
  for (int jj = iter - mtl2d; jj <= iter; ++jj) { test_len.push_back(jj); }

  mat iter_abs_mat = abs(as<mat>(iter_record_par));
  
  for (int j = 0; j < names_iter.length(); ++j) {
    if (std::find(params_to_test.begin(), params_to_test.end(), p_elems[j].c_str())
	!= params_to_test.end()) {
      mat test_par = iter_abs_mat.submat(test_len1-1, j, test_len2-1, j);
      if (arma::any(vectorise(test_par) == 0 )) {
	test_par = test_par + ones(test_len2-test_len1+1,1);
      }
      vec v = vectorise(test_par);
      double test_loglog = (std::log(v[v.n_elem-1]) - std::log(v[0])) /
	(std::log(test_len[test_len.size()-1]) - std::log(test_len[0]));
      test_conv[j] = test_loglog;
    }
  }
  
  vec test_test_conv = test_conv;
  if (test_test_conv.has_nan()) {
    return List::create(
      _["convergence"]=-2,
      _["mesages"] = "The log-log degeneracy test produced NAs.\n");
  }

  if (arma::any(abs(test_test_conv) > tol)) {
    CharacterVector conv_iter_names;
    CharacterVector nonconv_iter_names;
    std::ostringstream msg;
    CharacterVector msglist;
    for (int j = 0; j < names_iter.length(); ++j) {
      if (std::abs(test_test_conv[j]) > tol) {
	  msg << "Warning: the  " << names_iter[j] <<
	    "  parameter value has not converged.\n";
	  msglist.push_back(msg.str());
	  msg.str("");
	  nonconv_iter_names.push_back(names_iter[j]);
	}
	else {
	  conv_iter_names.push_back(names_iter[j]);
	}
    }

    msglist.push_back(
      "Type MARSSinfo(\"convergence\") for more info on this warning.\n");
    return List::create(
      _["convergence"]=1,
      _["messages"]=msglist,
      _["not.converged.params"]=nonconv_iter_names,
      _["converged.params"]=conv_iter_names
      );
  }
  
  return List::create(
    _["convergence"] = 0,
    _["messages"] = R_NilValue,
    _["not.converged.params"] = R_NilValue,
    _["converged.params"] = names_iter
    );
}

/**
   Contains frequently-repeated testing code called following updater
   execution.  Returns an empty string for success, non-empty for failure.
*/
std::string newkf(List& MLEobj, const std::string& elem, bool controlsafe,
		  bool tagfixed,
		  CharacterVector& msg_kf,
		  CharacterVector& msg_kem,
		  List& kf, List& Ey, int iter)
{
  if (controlsafe && !tagfixed) {
    List newkf = rerunkf(elem, MLEobj, iter);
    if (!as<bool>(newkf["ok"])) {
      msg_kf.push_back(as<String>(newkf["msg.kf"]));
      return as<std::string>(newkf["stop.msg"]);
    }
    else {
      kf = newkf["kf"];
      MLEobj["kf"] = kf;
      MLEobj["logLik"] = kf["logLik"];
      Ey = marsshatyt(MLEobj);
      MLEobj["Ey"] = Ey;
      msg_kem.push_back(as<String>(newkf["msg_kem"]));
    }
  }
  return std::string();
}

List rerunkf(const std::string& elem, List& MLEobj, int iter)
{
  double cvg2 = 1e-10;;
  if (iter == 1) {
    List ctrl = MLEobj["control"];
    cvg2 = 1 + as<double>(ctrl["abstol"]);
  }

  std::string msg_kf;   // if not set at return time, replace with R_NilValue
  SEXP msg_kem = R_NilValue;
    
  double loglike_old = MLEobj["logLik"];

  List kf = marsskf_impl(MLEobj);
    
  mat kfx0T = as<mat>(kf["x0T"]);
  mat kfxtT = as<mat>(kf["xtT"]);
    
  if ((as<List>(MLEobj["control"]))["demean.states"]) {
    mat cmb0t = join_rows(kfx0T, kfxtT);
    vec xbar(cmb0t.n_rows);
    for (unsigned int i=0; i < cmb0t.n_rows; ++i) {
      xbar[i] = std::accumulate(cmb0t.begin_row(i),
				cmb0t.end_row(i), 0) / cmb0t.n_cols;
    }
    kfx0T.each_col() -= xbar;
    kfxtT.each_col() -= xbar;
  }

  if (kf["ok"]) {
    double loglike_new = kf["logLik"];
    if (1 < iter && is_finite(loglike_old) && is_finite(loglike_new)) {
      cvg2 = loglike_new - loglike_old;
    }

    if (2 < iter && cvg2 < -std::sqrt(datum::eps)) {
      if ((as<List>(MLEobj["control"]))["trace"]) {
	msg_kem = wrap(paster("iter="," LogLike DROPPED in ",elem," update, "
			      "logLik old=",loglike_old," new=",loglike_new,"\n"));
      }
      else {
	msg_kem = wrap(paster("MARSSkem: The soln became unstable and logLik "
			      "DROPPED in the ",elem," updates.\n"));
      }
    }
    return List::create(
      _["kf"] = kf,
      _["msg.kem"] = msg_kem,
      _["msg.kf"] = R_NilValue,
      _["ok"] = true
      );
  }
  else {
    return List::create(
      _["ok"] = false,
      _["msg.kf"] = paster("iter=",iter," ",elem," update",
			   as<std::string>(kf["errors"]),"\n"),
      _["stop.msg"] = paster("Stopped at iter=",iter," in MARSSkem after ",elem,
		     " update: numerical errors in ",
		     as<std::string>(MLEobj["fun.kf"]),".\n")
      );
  }
}

String stabilitycheck(List& MLEobjiter, const List& timevarying,
		      const NumericMatrix& XX,
		      const std::string& Xtag, int thirddim, int iter)
{
  bool tv = as<bool>(timevarying[Xtag]);
  mat X = as<mat>(XX);
  return stabilitycheck_impl(MLEobjiter, tv, X, Xtag, thirddim, iter);
}

std::string stabilitycheck_impl1(List& MLEobjiter, bool timevarying, const mat& XX,
				const std::string& Xtag, int thirddim, int iter)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {
    mat X = XX;
    for (int i = 0; i < thirddim; ++i) {
      if (timevarying && i > 0) {
	try {
	  NumericVector Xnv = as<NumericVector>
	    (parmat(MLEobjiter, Xtag.c_str(), i)[Xtag.c_str()]);
	  X = as<mat>(Xnv);
	}
	catch(not_compatible e) {
	  return paster("error in stabilitycheck at iter=",iter, ", timestep ",
			i+1, ": parmat() return value not "
			"compatible with NumericVector.");
	}
	catch(not_a_matrix e) {
	  return paster("error in stabilitycheck at iter=",iter, ", timestep ",
			i+1, ": parmat() return value cannot be cast as matrix.");
	}

	try {
	  vec eigval = eig_sym(X);
	  if (any(eigval < 0)) {
	    return(paster("Stopped at iter ",iter," in MARSSkem: solution became "
			  "unstable. ",Xtag," update is not positive "
			  "definite.\n"));
	  }
	}
	catch(std::runtime_error) {
	  return paster("error in stabilitycheck. eig_sym() decomposition "
			"failed for ", Xtag);
	}
	catch(std::logic_error e) {
	  return paster("error in stabilitycheck. eig_sym() decomposition "
			"failed for ",Xtag,". ",e.what()," - ");
	}
      }
    }
    return "";
  }
  catch(std::exception e) {
    std::string tmpstr = paster("error in stabilitycheck for ",Xtag,". \n");
    std::string tmpstr2 = armaerror.str().empty() ?
      "" : paster("armadillo says: ", armaerror.str());
    return String(tmpstr + tmpstr2);
  }	
}

std::string stabilitycheck_impl(List& MLEobjiter, bool timevarying, const mat& XX,
				const std::string& Xtag, int thirddim, int iter)
{
  std::ostringstream armaerror;
  set_cerr_stream(armaerror);

  try {
    IntegerVector iis(thirddim);
    for (int j = 0; j < thirddim; ++j) {
      iis[j] = j;
    }
    cube thecube(XX.n_rows, XX.n_cols, thirddim);
    try {
      if (thirddim == 1) {
	thecube.slice(0) = XX;
      }
      else {
	thecube = as<cube>(parmat(MLEobjiter, Xtag.c_str(), iis)[Xtag.c_str()]);
      }
    }
    catch(not_compatible e) {
      return paster("error in stabilitycheck at iter=",iter, 
		    ": parmat() return value not "
		    "compatible with NumericVector.");
    }
    catch(not_a_matrix e) {
      return paster("error in stabilitycheck at iter=",iter,
		    ": parmat() return value cannot be cast as matrix.");
    }
    
    for (int i = 0; i < thirddim; ++i) {
      if (timevarying && i > 0) {
	mat X = thecube.slice(i);
	try {
	  vec eigval = eig_sym(X);
	  if (any(eigval < 0)) {
	    return(paster("Stopped at iter ",iter," in MARSSkem: solution became "
			  "unstable. ",Xtag," update is not positive "
			  "definite.\n"));
	  }
	}
	catch(std::runtime_error) {
	  return paster("error in stabilitycheck. eig_sym() decomposition "
			"failed for ", Xtag);
	}
	catch(std::logic_error e) {
	  return paster("error in stabilitycheck. eig_sym() decomposition "
			"failed for ",Xtag,". ",e.what()," - ");
	}
      }
    }
    // success
    return "";
  }
  catch(std::exception e) {
    std::string tmpstr = paster("error in stabilitycheck for ",Xtag,". \n");
    std::string tmpstr2 = armaerror.str().empty() ?
      "" : paster("armadillo says: ", armaerror.str());
    return String(tmpstr + tmpstr2);
  }	
}

std::string V0B_check(const NumericVector& par1V0, const NumericVector& par1B)
{
  std::string eigsym_arg = "V0";
  try {
    mat V0 = as<mat>(par1V0);
    vec eigval = eig_sym(V0);
    if (any(eigval < 0)) {
      mat B = as<mat>(par1B);
      eigval = eig_sym(B);
      if (any(abs(eigval) > 1)) {
	return "Your B matrix is outside the unit circle. This is "
	  "likely the problem.\n";
      }
      else {
	return " ";
      }
    }
  }
  catch(not_a_matrix e) {
    return  "Error in V0B_check: parmat() return value cannot be cast to matrix.";
  }
  catch(std::runtime_error) {
    return std::string{"Error in V0B_check: eig_sym() decomposition "
	"failed for "} + eigsym_arg;
  }
  catch(std::logic_error e) {
    return std::string{"Error in V0B_check: eig_sym() decomposition "
	"failed for "} + eigsym_arg;
  }

  return "";
}
