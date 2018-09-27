#include "marsskfas.h"
#include "genutils.h"
#include "marssutils.h"

using namespace arma;
using namespace Rcpp;

List marsskfas_impl(List& MLEobj, bool only_logLik, bool return_lag_one,
		    bool return_kfas_model)
{
  List modelObj = MLEobj["marss"];
  List control = MLEobj["control"];
  bool diffuse = modelObj["diffuse"];
  List model_dims = modelObj.attr("model.dims");

  int n = as<IntegerVector>(model_dims["data"])[0];
  int TT = as<IntegerVector>(model_dims["data"])[1];
  int m = as<IntegerVector>(model_dims["x"])[0];

  int g1 = as<IntegerVector>(model_dims["Q"])[0];
  int h1 = as<IntegerVector>(model_dims["R"])[0];
  int l1 = as<IntegerVector>(model_dims["L"])[0];

  List par_1 = parmat(MLEobj);

  std::map<std::string, arma::mat> matpar_1;
  if (!RObject(par_1.attr("names")).isNULL()) {
    CharacterVector parnames = par_1.attr("names");
    for (String parname : parnames) {
      matpar_1[parname] = as<mat>(as<NumericMatrix>(par_1[parname]));
    }
  }

  mat t_B = matpar_1["B"].t();
  
  IntegerVector TTindices(TT);
  IntegerVector TTindices1(TT-1);
  for (int i = 0; i < TT; ++i) {
    TTindices[i] = i;
    if (i < TT - 1) {
      TTindices1[i] = i + 1;
    }
  }

  cube Zt;
  cube stack_Zt;
  
  if (as<IntegerVector>(model_dims["Z"])[2] == 1 &&
      as<IntegerVector>(model_dims["A"])[2] == 1) {
    mat& mAt = matpar_1["A"];
    mat& mZt = matpar_1["Z"];
    Zt = cube(mAt.n_rows, mAt.n_cols + mZt.n_cols, 1);
    Zt.subcube(0, 0, 0, mZt.n_rows-1, mZt.n_cols-1, 0) = mZt;
    Zt.subcube(0, mZt.n_cols, 0, mZt.n_rows-1, mAt.n_cols + mZt.n_cols - 1, 0) =
      mAt;
    stack_Zt = cube(n, 2*(m+1), 1, fill::zeros);
    stack_Zt.subcube(0, 0, 0, n-1, m, 0) = Zt;
  }
  else {
    Zt = cube(n, m+1, TT, fill::zeros);
    stack_Zt = cube(n, 2*(m+1), TT, fill::zeros);
    std::map<std::string, cube> pars = parmat_cube(MLEobj, {"Z", "A"}, TTindices);
    Zt.subcube(0, 0, 0, n-1, m-1, TT-1) = pars["Z"];
    Zt.subcube(0, m, 0, n-1, m, TT-1) = pars["A"];
    stack_Zt.subcube(0, 0, 0, n-1, m, TT-1) = Zt;
  }

  cube Tt;
  cube stack_Tt;

  if (as<IntegerVector>(model_dims["B"])[2] == 1 &&
      as<IntegerVector>(model_dims["U"])[2] == 1) {

    Tt = cube(m + 1, m + 1, 1, fill::zeros);
    Tt.subcube(0, 0, 0, m-1, m-1, 0) = matpar_1["B"];
    Tt.subcube(0, m, 0, m-1, m, 0) = matpar_1["U"];
    Tt(m, m, 0) = 1;

    stack_Tt = cube(2*(m+1), 2*(m+1), 1, fill::zeros);
    stack_Tt.subcube(0, 0, 0, m, m, 0) = Tt;
    mat dummy_eye(m+1, m+1, fill::eye);
    stack_Tt.subcube(m+1, 0, 0, 2*m+1, m, 0) = dummy_eye;
  }
  else {
    Tt = cube(m+1, m+1, TT, fill::zeros);
    stack_Tt = cube(2*(m+1), 2*(m+1), TT, fill::zeros);

    std::map<std::string,cube>  pars = parmat_cube(MLEobj, {"B", "U"}, TTindices);
    Tt.subcube(0, 0, 0, m-1, m-1, TT-2) =
      pars["B"].subcube(0, 0, 1, m-1, m-1, TT-1);
    Tt.subcube(0, m, 0, m-1, m, TT-2) =
      pars["U"].subcube(0, 0, 1, m-1 , 0, TT-1);
    Tt.subcube(m, m, 0, m, m, TT-1) = cube(1, 1, TT, fill::ones);
    Tt.subcube(0, 0, TT-1, m, m, TT-1) = mat(m+1, m+1, fill::ones);

    for (int i = 0; i < TT; ++i) {
      stack_Tt.subcube(m+1, 0, i, 2*m+1, m, i) = eye<mat>(m+1, m+1);
    }
    stack_Tt.subcube(0, 0, 0, m, m, TT-1) = Tt;
  }

  cube Ht;

  if (as<IntegerVector>(model_dims["R"])[2] == 1) {
    mat& mHt = matpar_1["H"];
    mat& mRt = matpar_1["R"];
    mat mtcp = (mHt * mRt) * mHt.t();
    Ht = cube(mtcp.n_rows, mtcp.n_cols, 1);
    Ht.slice(0) = mtcp;
  }
  else {
    Ht = cube(n, n, TT);
    cube tmpHcube = parmat_cube(MLEobj, {"H"}, TTindices)["H"];
    cube tmpRcube = parmat_cube(MLEobj, {"R"}, TTindices)["R"];
    for (int i = 0; i < TT; ++i) {
      mat& mHt = tmpHcube.slice(i);
      mat& mRt = tmpRcube.slice(i);
      Ht.slice(i) = (mHt * mRt) * mHt.t();
    }
  }

  cube Qt;
  cube stack_Qt;

  if (as<IntegerVector>(model_dims["Q"])[2] == 1) {
    Qt = cube(g1, g1, 1);
    Qt.slice(0) = matpar_1["Q"];
    stack_Qt = cube(2*g1, 2*g1, 1, fill::zeros);
    stack_Qt.subcube(0, 0, 0, g1-1, g1-1, 0) = Qt;
  }
  else {
    Qt = cube(g1, g1, TT, fill::zeros);
    stack_Qt = cube(2*g1, 2*g1, TT, fill::zeros);

    Qt.subcube(0, 0, 0, g1-1, g1-1, TT-2) =
      parmat_cube(MLEobj, "Q", TTindices1)["Q"];
    stack_Qt.subcube(0, 0, 0, g1-1, g1-1, TT-2) =
      Qt.subcube(0, 0, 0, g1-1, g1-1, TT-2) ;
  }

  cube Rt;
  cube stack_Rt;
  if (as<IntegerVector>(model_dims["R"])[2] == 1) {
    Rt = cube(m+1, g1, 1, fill::zeros);
    stack_Rt = cube(2*(m+1), 2*g1, 1, fill::zeros);
    
    Rt.subcube(0, 0, 0, m-1, g1-1, 0) = matpar_1["G"];
    stack_Rt.subcube(0, 0, 0, m, g1-1, 0) = Rt;
  }
  else {
    Rt = cube(m+1, g1, TT, fill::zeros);
    stack_Rt = cube(2*(m+1), 2*g1, TT, fill::zeros);

    cube tmppar=parmat_cube(MLEobj, "G", TTindices1)["G"];

    Rt.subcube(0, 0, 0, m-1, g1-1, TT-2) =
      parmat_cube(MLEobj, "G", TTindices1)["G"];
    stack_Rt.subcube(0, 0, 0, m, g1-1, TT-1) = Rt;
  }
  
  mat x00;
  mat V00;
  mat x10;
  mat V10;

  double testtinitx = modelObj["tinitx"];

  if (testtinitx == 0) {
    x00 = matpar_1["x0"];
    V00 = (matpar_1["L"] * matpar_1["V0"]) * matpar_1["L"].t();
    x10 = matpar_1["B"] * x00 + matpar_1["U"];
    V10 = matpar_1["B"] * V00 * t_B + matpar_1["Q"];
  }
  else {
    x10 = matpar_1["x0"];
    x00 = mat(10, 1, fill::zeros);
    V10 = matpar_1["V0"];
    V00 = zeros(m, m);
  }

  mat onesrow(1, x10.n_cols, fill::ones);
  mat a1 = join_cols(x10, onesrow);
  mat stack_a1 = join_cols(join_cols(join_cols(x10, onesrow),x00),onesrow);

  mat P1inf(m+1, m+1, fill::zeros);
  mat stack_P1inf(2*(m+1), 2*(m+1), fill::zeros);
  mat P1(m+1, m+1, fill::zeros);
  mat stack_P1(2*(m+1), 2*(m+1), fill::zeros);
  
  if (diffuse) {
    P1inf.submat(0, 0, m-1, m-1) = V10;
    stack_P1inf.submat(0, 0, m-1, m-1) = V10;
    stack_P1inf.submat(m+1, m+1, 2*m, 2*m) = V00;
  }
  else {
    P1.submat(0, 0, m-1, m-1) = V10;
    stack_P1.submat(0, 0, m-1, m-1) = V10;
    stack_P1.submat(m+1, m+1, 2*m, 2*m) = V00;
    stack_P1.submat(0, m+1, m-1, 2*m) = matpar_1["B"] * V00;
    stack_P1.submat(m+1, 0, 2*m, m-1) = V00 * t_B;
  }

  List kfas_model;
  List envlist = List::create(
    _["a1"]=a1,
    _["control"]=control,
    _["diffuse"]=diffuse,
    _["Ht"]=Ht,
    // _["i"]=i,
    _["m"]=m,
    _["MLEobj"]=MLEobj,
    _["model.dims"]=model_dims,
    _["modelObj"]=modelObj,
    _["n"]=n,
    _["only.logLik"]=only_logLik,
    _["P1"]=P1,
    _["P1inf"]=P1inf,
    _["par.1"]=par_1,
    _["Qt"]=Qt,
    _["return.kfas.model"]=return_kfas_model,
    _["return.lag.one"]=return_lag_one,
    _["Rt"]=Rt
    );

  envlist.push_back(stack_a1,"stack.a1");
  envlist.push_back(stack_P1,"stack.P1");
  envlist.push_back(stack_P1inf,"stack.P1inf");
  envlist.push_back(stack_Qt,"stack.Qt");
  envlist.push_back(stack_Rt,"stack.Rt");
  envlist.push_back(stack_Tt,"stack.Tt");
  envlist.push_back(stack_Zt,"stack.Zt");
  envlist.push_back(t_B,"t.B");
  envlist.push_back(Tt,"Tt");
  envlist.push_back(TT,"TT");
  envlist.push_back(V00,"V00");
  envlist.push_back(V10,"V10");
  envlist.push_back(x00,"x00");
  envlist.push_back(x10,"x10");

  NumericMatrix nmyt = clone(as<NumericMatrix>(modelObj["data"]));
  NumericVector nvyt(nmyt);
  nvyt = ifelse(is_na(nvyt), NA_REAL, nvyt);
  nmyt = transpose(nmyt);
  envlist.push_back(nmyt, "yt");
  envlist.push_back(Zt,"Zt");

  kfas_model = SSModel3(envlist);

  if (only_logLik) {
    return List::create(_["logLik"]=logLik(kfas_model));
  }

  vec diag_R = diagvec(matpar_1["R"]);

  List ks_out = KFS(kfas_model, false);

  if (packageVersion("KFAS") != "0.9.11") {
    ks_out["a"] = transpose(as<NumericMatrix>(ks_out["a"]));
    ks_out["alphahat"] = transpose(as<NumericMatrix>(ks_out["alphahat"]));
  }
  NumericMatrix ks_out_a = ks_out["a"];
  NumericMatrix ks_out_alphahat = clone(as<NumericMatrix>(ks_out["alphahat"]));

  cube ksoV(as<cube>(ks_out["V"]));
  cube cVtT = ksoV.tube(0, 0, m-1, m-1);
  NumericVector VtT = wrap(cVtT);
  
  cube ksoP(as<cube>(ks_out["P"]));
  cube cVtt1 = ksoP.subcube(0, 0, 0, m-1, m-1, TT-1);
  NumericVector Vtt1 = wrap(cVtt1);

  NumericVector Vtt1T;
  if (return_lag_one) {
    cube cVtt1T = ksoV.tube(0, m+1, m-1, 2*m);
    Vtt1T = wrap(cVtt1T);
  }

  if (any(diag_R == 0)) {
    VtT = ifelse(sapply(VtT,[](double x){return x < datum::eps;}), 0, VtT);
    Vtt1 = ifelse(sapply(Vtt1,[](double x){return x < datum::eps;}), 0, Vtt1);
    if (return_lag_one) {
      Vtt1T = ifelse(sapply(Vtt1T,[](double x){return x < datum::eps;}), 0, Vtt1T);
    }
  }

  NumericVector x10T = ks_out_alphahat.column(0);
  x10T.erase(m, ks_out_alphahat.nrow());
  IntegerVector x10Tdims {m, 1};
  x10T.attr("dim") = x10Tdims;
  
  NumericVector x0T;
  
  NumericMatrix V10T = wrap(cVtT.slice(0));
  NumericMatrix V0T;
  
  if (testtinitx == 1) {
    x00 = mat(m, 1);
    x00.fill(NA_REAL);
    V00 = mat(m, m);
    V00.fill(NA_REAL);
    x0T = x10T;
    V0T = V10T;
  }
  else {
    NumericMatrix Vtt1_1 = sub3D(Vtt1,0);
    mat mVtt1_1(Vtt1_1.begin(), Vtt1_1.nrow(), Vtt1_1.ncol(), false);
    mat Vinv = pcholinv(mVtt1_1);
    if (m != 1) {
      Vinv = (Vinv.t() + Vinv) / 2.;
    }
    mat J0 = V00 * t_B * Vinv;
    NumericVector ks_out_a_col = ks_out_a.column(0);
    int erase_end = ks_out_a_col.length();
    ks_out_a_col.erase(m, erase_end);

    vec t1 = as<vec>(x10T) - as<vec>(ks_out_a_col);
    mat t2 = J0 * t1;
    x0T = x00 + t2;
    V0T = wrap(V00 + J0 * (cVtT.slice(0) - cVtt1.slice(0)) * J0.t());
  }

  List rtn_list = List::create(
    _["xtT"] = add_default_rownames(
      ks_out_alphahat(Range(0, m-1), Range(0, ks_out_alphahat.ncol()-1))),
    _["VtT"] = VtT,
    _["Vtt1T"] = return_lag_one ? wrap(Vtt1T) : R_NilValue,
    _["x0T"] = x0T,
    _["V0T"] = V0T,
    _["x10T"] = add_default_rownames(NumericMatrix(x10T.length(), 1, x10T.begin())),
    _["V10T"] = V10T,
    _["x00T"] = x00,
    _["V00T"] = V00,
    _["Vtt"] = String("Use MARSSkfss to get Vtt"),
    _["Vtt1"] = Vtt1,
    _["J"] = String("Use MARSSkfss to get J"),
    _["J0"] = String("Use MARSSkfss to get J0"),
    _["Kt"] = String("Use MARSSkfss to get Kt"),
    _["xtt1"] = add_default_rownames(ks_out_a(Range(0, m-1), Range(0,TT-1))),
    _["xtt"] = String("Use MARSSkfss to get xtt"),
    _["Innov"] = String("Use MARSSkfss to get Innov"),
    _["Sigma"] = String("Use MARSSkfss to get Sigma"),
    _["kfas.model"] = return_kfas_model ? wrap(kfas_model) : R_NilValue,
    _["logLik"] = ks_out["logLik"]
    );

  rtn_list.push_back(true, "ok");
  rtn_list.push_back(R_NilValue, "errors");

  return rtn_list;
}

