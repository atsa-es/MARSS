#include "genutils.h"

using namespace arma;
using namespace Rcpp;

namespace dim_names {
  int push(const RObject& newnames, std::vector<std::string>& allnames,
	   int rowcol)
  {
    if (newnames != R_NilValue) {
      CharacterVector dnames;
      if (rowcol < 0) {
	dnames = as<CharacterVector>(newnames);
      }
      else {
	List lnames = as<List>(newnames);
	if (lnames[rowcol] != R_NilValue) {
	  dnames = as<CharacterVector>(lnames[rowcol]);
	}
	else {
	  return 0;
	}
      }
      for (String dname : dnames) {
	allnames.push_back(dname);
      }
    }
    return allnames.size();
  }

  int pad(std::vector<std::string>& allnames, int needed)
  {
    for (int i = allnames.size(); i < needed; ++i) {
      allnames.push_back("");
    }
    return allnames.size();
  }

  bool array_fits(const RObject& newnames, const Dimension& dims, int rowcol)
  {
    // handles 1-D arrays only (e.g., "names" attribute)
    
    if (newnames == R_NilValue) { return false; }
    
    List test = as<List>(newnames);
    
    return test.length() == dims[rowcol];
  }

  int length(const NumericVector& v, int rowcol)
  {
    if (v.attr("dim") == R_NilValue)
    {
      return 1;
    }

    IntegerVector dims = v.attr("dim");

    return (dims.size() < rowcol + 1) ? v.size() : dims[rowcol];
  }

  // Wrapper for get() which returns a CharacterVector;
  // this is empty if get() returns an R_NilValue
  CharacterVector getcv(const NumericVector& v1, const NumericVector& v2,
			const RObject& v1newnames, const RObject& v2newnames,
			const Dimension& resultdims, int rowcol, bool bindcols)
  {
    RObject ret = get(v1, v2, v1newnames, v2newnames, resultdims,
		      rowcol, bindcols);
    return (ret == R_NilValue) ? CharacterVector() : as<CharacterVector>(ret);
  }
  
  RObject get(const NumericVector& v1, const NumericVector& v2,
	      const RObject& v1newnames, const RObject& v2newnames,
	      const Dimension& resultdims, int rowcol, bool bindcols)
  {
    /*
     * rowcol = 0(1) : return row(col) dimnames
     * bindcols = true(false) : binding columns(rows)
     */

    std::vector<std::string> allnames;
  
    if (rowcol != 0 && rowcol != 1) {
      return wrap(allnames);
    }

    int needed = resultdims[rowcol];
  
    if (rowcol == (bindcols ? 1 : 0)) {  // colnames for binding cols
      // rownames for binding rows
      if (allnames.size() < (unsigned int)needed ) {
	if (push(v1.attr("dimnames"), allnames, rowcol) == 0) {
	  push(v1newnames, allnames);
	}
	pad(allnames, dim_names::length(v1, rowcol));
	int v1cnt = allnames.size();
	if (push(v2.attr("dimnames"), allnames, rowcol) == v1cnt) {
	  push(v2newnames, allnames);
	}
	pad(allnames, length(v1, rowcol) + length(v2, rowcol));
      }
    }
    else { // rownames for binding cols
      // colnames for binding rows
      if (allnames.size() < (unsigned int)needed) {
	if (push(v1.attr("dimnames"), allnames, rowcol) == 0) {
	  if (v1.attr("dim") == R_NilValue) {
	    push(v1.attr("names"), allnames);
	  }
	  if (allnames.size() == 0) {
	    if (push(v2.attr("dimnames"), allnames, rowcol) == 0) {
	      if (v2.attr("dim") == R_NilValue) {
		push(v2.attr("names"), allnames);
	      }
	    }
	  }
	}
      }
    }

    // if allnames is nothing but padding, empty it
    if (std::all_of(allnames.begin(), allnames.end(),
		    [](std::string& x) {return x == "";})) {
      allnames.clear();
    }
  
    CharacterVector cvallnames = wrap(allnames);
    if (cvallnames.size() > 0) {
      for (int i = cvallnames.size(); i < needed; ++i) {
	cvallnames.push_back(NA_STRING);
      }
      return cvallnames;
    }
    else {
      return R_NilValue;
    }
  }
}

std::map<int,Rcpp::IntegerVector> IdxVec::vecs;

void IdxVec::addvec(int len) {
  Rcpp::IntegerVector newvec(len);
  for (int i =0; i < len; ++i) { newvec[i] = i; }
  vecs[len] = newvec;
}
  
Rcpp::IntegerVector IdxVec::get(int len) {
  if (vecs.find(len) == vecs.end()) { addvec(len); }
  return vecs[len];
}


NumericMatrix add_default_rownames(NumericMatrix m)
{
  StringVector rns;
  std::ostringstream rn;
  for (int i = 1; i <= m.nrow(); ++i) {
    rn.str("");
    rn << "custom" << i;
    rns.push_back(rn.str());
  }
  m.attr("dimnames") = List::create(rns,StringVector());
  return m;
}

mat choleskyinverter(mat& target)
// throws runtime_error
{
  mat targetinverse(target.n_rows, target.n_cols, fill::zeros);
    
  vec targetdiag = target.diag();
  
  if (any(targetdiag == 0)) {
    if (any(targetdiag != 0)) {

      mat b = trimmat(target, targetdiag, true, true);
      mat choltarget = trimatl(chol(b).t());
      mat invl = solve(choltarget, eye(b.n_rows, b.n_rows));
      b = solve(trimatu(choltarget.t()), invl);

      mat omgx = eye(targetdiag.n_elem, targetdiag.n_elem);
      uvec nonzerodiag = find(targetdiag != 0);
      omgx = omgx.rows(nonzerodiag);
      targetinverse = omgx.t() * b * omgx;
    }
    else {
      targetinverse.zeros();
    }
  }
  else {
    if (target.n_elem == 1) {
      target(0, 0) = 1 / target(0, 0);
    } else {
      mat choltarget = trimatl(chol(target).t());
      mat inv1 = solve(choltarget, eye(target.n_rows, target.n_rows));
      target = solve(trimatu(choltarget.t()), inv1);
    }
    targetinverse = target;
  }
  return targetinverse;
}

void choleskyinverter_old(mat& denom)
// throws runtime_error
{
  if (denom.n_elem == 1) {
    denom(0, 0) = 1 / denom(0, 0);
  } else {
    mat choldenom = trimatl(chol(denom).t());
    mat inv1 = solve(choldenom, eye(denom.n_rows, denom.n_rows));
    denom = solve(trimatu(choldenom.t()), inv1);
  }
}

std::map<std::string, cube> cube_list_to_map(const List& star)
{
  std::map<std::string, cube> mstar;
  CharacterVector keys = star.attr("names");
  for (String key : keys) {
    mstar[key.get_cstring()] = makecube(star[key]);
  }

  return mstar;
}

List cube_map_to_list(const std::map<std::string,cube>& m)
{
  List rl = List::create();

  for (std::map<std::string,cube>::const_iterator it = m.cbegin();
       it != m.cend(); ++it) {
    cube c = it->second;
    rl[it->first] = c;
  }

  return rl;
}

CharacterVector cv_append(CharacterVector& v1, const CharacterVector& v2)
{
  for (auto v2elem : v2) {
    v1.push_back(v2elem);
  }
  // std::copy(v2.begin(), v2.end(), std::back_inserter(v1));
  return v1;
}

uvec diag_indices(int nrows)
{
  class gdiags {
    int i, m;
  public:
    gdiags(int mm) : i(-mm-1), m(mm){}
    int operator()() {return i += (m+1);}
  };

  return uvec(nrows).imbue(gdiags(nrows));
}

Dimension effective_dims(NumericVector& v, bool bindcols) {
  Dimension vdim;
  if (v.length() == 0) {
    vdim = Dimension(0, 0);
  }
  else {
    if (bindcols) {
      vdim = (v.attr("dim") == R_NilValue) ? Dimension(v.length(), 1) :
  	v.attr("dim");
    }
    else {
      vdim = (v.attr("dim") == R_NilValue) ? Dimension(1, v.length()) :
  	v.attr("dim");
    }
  }
  return vdim;
}


void ensure_dims(NumericVector& v) {
  if (v.attr("dim") == R_NilValue) {
    v.attr("dim") = IntegerVector::create(1,v.length());
  }
}

SEXP getlistelement(const List& thelist, const std::string& label)
{
  // check for any names at all
  List namelist(thelist.names());
  if (namelist.size() == 0) {
    return R_NilValue;
  }

  // check for label
  CharacterVector namevec(thelist.names());
  if (namevec.end() == std::find(namevec.begin(), namevec.end(), String(label))) {
    return R_NilValue;
  }

  return thelist[label];
}

void list_assign(List& l, const std::string& name, const RObject& value)
{
  /*
   * add name -> value to list l
   * If a name element already exists, set it to value directly.
   * If no name element exists, add the pair via push_back()
   *   in order to avoid out-of-bounds errors.
   */
  
  if (l.containsElementNamed(name.c_str())) {
    l[name] = value;
  }
  else {
    l.push_back(value, name);
  }
}

bool is_diag(const arma::mat& m)
{
  return (!any(vectorise(m - diagmat(m))));
}

bool is_wholenumber(double t)
{
  return (t - round(t) < datum::eps);
}

void list_assign(List& target, const std::string& targettag,
		 const List& source, const std::string& sourcetag,
		 SEXP defaultobj)
{
  if ((source.attr("names") == R_NilValue) || 
      (!source.containsElementNamed(sourcetag.c_str()))) {
    target[targettag] = defaultobj;
  }
  else {
    RObject val = source[sourcetag];
    list_assign(target, targettag, val);
    // target[targettag] = source[sourcetag];
  }
}

NumericVector list_to_nv(const List& l)
{
  NumericVector nv;
  for (List::const_iterator li=l.begin(); li != l.end(); ++li) {
    nv.push_back(*li);
  }

  nv.attr("names") = l.attr("names");
  nv.attr("dimnames") = l.attr("dinnames");
  nv.attr("dim") = l.attr("dim");
 
  return nv;
}

cube makecube(SEXP r, bool copy)
{
  NumericVector rr(r);
  IntegerVector rrdims = rr.attr("dim");
  if (rrdims.size() == 3) {
    return cube(rr.begin(), rrdims[0], rrdims[1], rrdims[2], copy);
  } else {
    std::cerr << "not a cube: ";
    for (int i = 0; i < rrdims.size() - 1; i++) {
      std::cerr << rrdims[i] << " x ";
    }
    std::cerr << rrdims[rrdims.size() - 1] << endl;
    return cube();
  }
}

mat makemat(SEXP r, bool copy)
{
  NumericMatrix rr(r);
  return mat(rr.begin(), rr.nrow(), rr.ncol(), copy);
}

NumericVector rcbind_impl(SEXP a1, SEXP a2, bool bindcols, SEXP rnames,
			  SEXP cnames)
{
  NumericVector v1;
  NumericVector v2;
  v2.attr("dim") = Dimension(0,0);

  // If neither argument exists, return an empty NV.
  // If there is only a single non-null argument, store it in the
  // first argument NV wrapper
    
  if (a1 == R_NilValue) {
    if (a2 == R_NilValue) {
      return v1;
    }
    else {
      v1 = as<NumericVector>(a2);
    }
  }
  else {
    v1 = as<NumericVector>(a1);
    if (a2 != R_NilValue) {
      v2 = as<NumericVector>(a2);
    }
  }

  CharacterVector colnames = (cnames==R_NilValue) ?
    CharacterVector{} : as<CharacterVector>(cnames);
  CharacterVector rownames = (rnames==R_NilValue) ?
    CharacterVector{} : as<CharacterVector>(rnames);

  Dimension v1dims = effective_dims(v1, bindcols);
  Dimension v2dims = effective_dims(v2, bindcols);
  Dimension rvdims = bindcols ?
    Dimension(std::max(v1dims[0],v2dims[0]), v1dims[1] + v2dims[1]) :
    Dimension(v1dims[0] + v2dims[0], std::max(v1dims[1],v2dims[1]));

  NumericVector rv;


  if (bindcols) {
    // If the width of the second argument doesn't match that of the first,
    // wrap its rows in order to fill the first argument's width.
    for (int j = 0; j < rvdims[1]; ++j) {
      int v10 = j * v1dims[0];
      int v20 = (j - v1dims[1]) * v2dims[0];
      for (int i = 0; i < rvdims[0]; ++i ) {
	double val = (j < v1dims[1]) ? v1[i+v10] : v2[(i%v2dims[0])+v20];
	rv.push_back(val);
      }
    }
  }
  else {
    // If the height of the second argument doesn't match that of the first,
    // wrap its columns in order to fill the first argument's height.
    
    for (int j = 0; j < rvdims[1]; ++j) {
      int v10 = j * v1dims[0];
      int v20 = j * v2dims[0];
      for (int i = 0; i < rvdims[0]; ++i) {
	double val = (i < v1dims[0]) ? v1[i+v10] :
	  v2[(i-v1dims[0]+v20) % v2.length()];
	rv.push_back(val);
      }
    }
  }

  rv.attr("dim") = rvdims;

  // Per cbind() - rbind() behavior: If the first arg has a "dimnames" 
  // attribute, propagate that to the result.  Otherwise, if the first
  // arg has a "names" attribute, use that for result's
  // colnames (dimnames[1]).
    
  CharacterVector rownamelist = dim_names::getcv(v1, v2, rownames, R_NilValue,
						 rvdims, 0, bindcols);
  CharacterVector colnamelist = dim_names::getcv(v1, v2, R_NilValue, colnames,
						 rvdims, 1, bindcols);
  rv.attr("dimnames") = List::create(rownamelist, colnamelist);
    
  return rv;
}

// bind columns, applying supplied names if available
NumericVector scbind(SEXP a1, SEXP a2, SEXP newrownames, SEXP newcolnames)
{
  return rcbind_impl(a1, a2, true, newrownames, newcolnames);
}

// bind rows, applying supplied names if available
NumericVector srbind(SEXP a1, SEXP a2, SEXP newrownames, SEXP newcolnames)
{
  return rcbind_impl(a1, a2, false, newrownames, newcolnames);
}

mat trimmat(const mat& x, const vec& diag, bool deletezerorows, bool deletezerocols)
{
  // Trims a matrix by removing rows and columns depending on diagonal values.
  // diag i       deletezerorows  deletezerocols  delete row i?  delete col i?
  //    zero          true            true             yes          yes
  // nonzero          true            true              no           no
  //    zero         false           false              no           no
  // nonzero         false           false             yes          yes
  //    zero          true           false             yes           no
  // nonzero          true           false              no          yes
  //    zero         false            true              no          yes
  // nonzero         false            true             yes           no

  mat rx(x);
  if (deletezerorows && deletezerocols) {
    for (int i = x.n_rows - 1; i >= 0; --i) {
      if (diag(i) == 0) {
	rx.shed_row(i);
	rx.shed_col(i);
      }
    }
  } else if (!deletezerorows && !deletezerocols) {
    for (int i = x.n_rows - 1; i >= 0; --i) {
      if (0 != diag(i)) {
	rx.shed_row(i);
	rx.shed_col(i);
      }
    }
  } else if (deletezerorows) {
    for (int i = x.n_rows - 1; i >= 0; --i) {
      if (diag(i) == 0) {
	rx.shed_row(i);
      } else {
	rx.shed_col(i);
      }
    }
  } else {
    for (int i = x.n_rows - 1; i >= 0; --i) {
      if (diag(i) == 0) {
	rx.shed_col(i);
      } else {
	rx.shed_row(i);
      }
    }
  }

  return rx;
}

NumericMatrix unname(NumericMatrix& x)
{
  x.attr("names") = R_NilValue;
  x.attr("dimnames") = R_NilValue;

  return x;
}

List truenull(List& l, std::string target)
{
  std::string cl = l.attr("class");
  if (l.containsElementNamed(target.c_str())) {
    int i = l.findName(target.c_str());
    l.erase(i);
    l.attr("class") = cl.c_str();
  }

  return l;
}

mat unvec(NumericVector& v, const IntegerVector& dim, bool copy)
{
  // skip the do-not-execute block
  return mat(v.begin(), dim[0], dim[1], copy);
}

void unzerorows(mat& IId, const mat& Mt, const vec& diagQ)
{
  mat tMt = trimmat(Mt, diagQ, false, false);
  vec zerorows(tMt.n_rows);
  for (unsigned int i = 0; i < tMt.n_rows; ++i) {
    zerorows[i] = all(tMt.row(i) == 0) ? 1 : 0;
  }
  uvec therows = find(diagQ == 0);
  mat subIId(IId.submat(therows, therows));
  subIId.diag() = zerorows;
  IId.submat(therows, therows) = subIId;
}

// (temporary) wrappers for R/MARSS functions

bool allequal(const List& lista, const List& listb)
{
  Function allequal_r("all.equal");
  return allequal_r(lista, listb);
}

List coef(List& object, const String& type,
	  const String& form, const String& what)
{
  object.attr("class") = "marssMLE";
  Function coef_r("coef.marssMLE","MARSS");
  return as<List>(coef_r(_["object"]=object,
			 _["type"]=type,
			 _["form"]=form,
			 _["what"]=what));
}

NumericVector coef_vec(List& object, const String& form, const String& what)
{
  object.attr("class") = "marssMLE";
  Function coef_r("coef.marssMLE","MARSS");
  NumericVector rvec = as<NumericVector>(coef_r(_["object"]=object,
						_["type"]="vector",
						_["form"]=form,
						_["what"]=what));
  rvec.attr("dim") = R_NilValue;
  
  return rvec;
}

// debugging utilities

void dumpcube(std::string name, cube thecube, bool dump = false)
{
  if (dump) {
    std::cout << name << " cube :" << thecube.n_rows << " x " <<
      thecube.n_cols;
    std::cout << " x " << thecube.n_slices << endl;
    for (unsigned int i = 0; i < thecube.n_slices; ++i) {
      std::cout << "i: " << i << " ";
      dumpmat("slice", thecube.slice(i), true);
    }
  }
}

void dumpintvec(std::string name, Rcpp::IntegerVector x, bool dump = false)
{
  std::cout << name << endl;
  for (Rcpp::IntegerVector::iterator i = x.begin(); i != x.end(); i++) {
    std::cout << *i << " ";
  }
  std::cout << endl;
}

void dumpmat(std::string name, mat themat, bool dump)
{
  if (dump) {
    dumpmat(name, themat, std::cout);
  }
}

void dumpmat(std::string name, mat themat, std::ostream& dout)
{
  dout << name << " matrix :" << themat.n_rows << " x " <<
    themat.n_cols << endl;
  for (unsigned int iii = 0; iii < themat.n_rows; iii++) {
    for (mat::row_iterator ri = themat.begin_row(iii);
	 ri != themat.end_row(iii); ri++) {
      dout << *ri << " ";
    }
    dout << endl;
  }
}


void dumpnumvec(std::string name, Rcpp::NumericVector x, bool dump = false)
{
  std::cout << name << endl;
  for (Rcpp::NumericVector::iterator i = x.begin(); i != x.end(); i++) {
    std::cout << *i << " ";
  }
  std::cout << endl;
}

void dumpobj(std::string name, Rcpp::List& MLEobj, bool dump = true)
{
  Rcpp::List model = MLEobj["marss"];
  Rcpp::List f = model["fixed"];
  Rcpp::NumericVector fel = f["A"];
  Rcpp::IntegerVector fel_dim = fel.attr("dim");

  mat mfelem(fel.begin(), fel_dim[0], fel_dim[1], false);

  cube felcube(fel.begin(), fel_dim[0], fel_dim[1], fel_dim[2], false);
  mfelem = felcube.slice(0);
  std::cout << name << mfelem(0, 0) << endl;;
}

void dumpvec(std::string name, vec thevec, bool dump = false)
{
  if (dump) {
    std::cout << name << " vector : " << thevec.n_elem << endl;
    for (vec::iterator i = thevec.begin(); i != thevec.end(); i++) {
      std::cout << *i << " ";
    }
    std::cout << endl;
  }
}

void savelist(const Rcpp::List& objects, const std::string& file, bool print)
{
  // Example : 
  // save MLEobj under the name "MLEobjname" in the file MLEobjfile
  // true ==> print a message to stdout when the file is saved
  // savelist(List::create(_["MLEobjname"]=MLEobj), "MLEobjfile", true);
  Environment e;
  CharacterVector names = objects.attr("names");
  if (names.size() > 0) {
    for (String name : names) {
      if (print) { std::cout << "Saving " << name.get_cstring() << "\n"; }
      e.assign(name, objects[name.get_cstring()]);
    }
    Function save_r("save");
    save_r(_["list"]=names, _["file"]=file, _["envir"]=e);
  }
}
