#ifndef _GENUTILS_H_
#define _GENUTILS_H_

#include "marssrcpp.h"

#include <ostream>
#include <sstream>
#include <string>

namespace dim_names{
  int push(const Rcpp::RObject& newnames, std::vector<std::string>& allnames,
	   int rowcol = -1);
  int pad(std::vector<std::string>& allnames, int needed);
  bool array_fits(const Rcpp::RObject& newnames, const Rcpp::Dimension& dims,
		  int rowcol);
  int length(const Rcpp::NumericVector& v, int rowcol);
  Rcpp::RObject get(const Rcpp::NumericVector& v1, const Rcpp::NumericVector& v2,
		    const Rcpp::RObject& v1newnames,
		    const Rcpp::RObject& v2newnames,
		    const Rcpp::Dimension& resultdims, int rowcol,
		    bool bindcols);
  Rcpp::CharacterVector getcv(const Rcpp::NumericVector& v1,
			      const Rcpp::NumericVector& v2,
			      const Rcpp::RObject& v1newnames,
			      const Rcpp::RObject& v2newnames,
			      const Rcpp::Dimension& resultdims, int rowcol,
			      bool bindcols);
}

class Paster {
private:
  std::ostringstream os;
  std::string paste_var(std::ostringstream& os) { return os.str(); }
  template<typename T, typename... Args>
  std::string paste_var(std::ostringstream& os, T val, Args... args) {
    os << val;
    return paste_var(os, args...);
  }

public :
  template<typename... Args>
  std::string operator()(Args... args) {os.str(""); return(paste_var(os, args...));}
};

class IdxVec {
private :
  static std::map<int, Rcpp::IntegerVector> vecs;
  static void addvec(int len);
  
public:
  static Rcpp::IntegerVector get(int len);
};

Rcpp::NumericMatrix add_default_rownames(Rcpp::NumericMatrix m);
arma::mat choleskyinverter(arma::mat&);
void choleskyinverter_old(arma::mat&);
std::map<std::string,arma::cube> cube_list_to_map(const Rcpp::List& star);

Rcpp::List cube_map_to_list(const std::map<std::string, arma::cube>& themap);
Rcpp::CharacterVector cv_append(Rcpp::CharacterVector& v1,
				const Rcpp::CharacterVector& v2);
arma::uvec diag_indices(int);
Rcpp::Dimension effective_dims(Rcpp::NumericVector& v, bool bindcols);
void ensure_dims(Rcpp::NumericVector& v);
SEXP getlistelement(const Rcpp::List&, const std::string&);
bool is_diag(const arma::mat&);
bool is_wholenumber(double);
void list_assign(Rcpp::List& l, const std::string& name,
		 const Rcpp::RObject& value);
arma::cube makecube(SEXP r, bool=false);
arma::mat makemat(SEXP, bool=false);
void list_assign(Rcpp::List& target, const std::string& targettag,
		 const Rcpp::List& source, const std::string& sourcetag,
		 SEXP defaultobj = R_NilValue);
Rcpp::NumericVector list_to_nv(const Rcpp::List& l);
Rcpp::NumericVector rcbind_impl(SEXP a1, SEXP a2, bool bindcols,
				SEXP rnames=R_NilValue, SEXP cnames=R_NilValue);
Rcpp::NumericVector scbind(SEXP a1, SEXP a2, SEXP newrownames=R_NilValue,
			   SEXP newcolnames=R_NilValue);
Rcpp::NumericVector srbind(SEXP a1, SEXP a2, SEXP newrownames=R_NilValue,
			   SEXP newcolnames=R_NilValue);
arma::mat trimmat(const arma::mat&, const arma::vec&, bool, bool);
Rcpp::List truenull(Rcpp::List& l, std::string target);
Rcpp::NumericMatrix unname(Rcpp::NumericMatrix&);
arma::mat unvec(Rcpp::NumericVector&, const Rcpp::IntegerVector&, bool=true);
void unzerorows(arma::mat&, const arma::mat&, const arma::vec&);

// wrappers
bool allequal(const Rcpp::List&, const Rcpp::List&);
Rcpp::List coef(Rcpp::List& object, const Rcpp::String& type="list",
		const Rcpp::String& form="", const Rcpp::String& what="par");
Rcpp::NumericVector coef_vec(Rcpp::List& object, const Rcpp::String& form="",
			     const Rcpp::String& what="par");

// debug utilities
void dumpcube(std::string, arma::cube, bool);
void dumpintvec(std::string, Rcpp::IntegerVector, bool);
void dumpmat(std::string, arma::mat, std::ostream &dout);
void dumpmat(std::string, arma::mat, bool=false);
void dumpnumvec(std::string, Rcpp::NumericVector, bool);
void dumpobj(std::string name, Rcpp::List& MLEobj, bool);
void dumpvec(std::string, arma::vec, bool);

void savelist(const Rcpp::List& objects, const std::string& file, bool print=false);
#endif /* _GENUTILS_H_ */
