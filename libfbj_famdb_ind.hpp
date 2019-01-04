// body file: libfbj_famdb_ind.cpp

#ifndef individual_class
#define individual_class

#include <map>
#include <string>
#include <vector>
#include <set>
#include <boost/logic/tribool.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using boost::logic::tribool;
using std::vector;
using std::string;
using std::pair;
class Individual ;
struct Family ;
typedef std::map<string, Individual*> IndMap;

inline std::string trb2str(boost::logic::tribool b)
{
	if		(b==true)	return "1";
	else if (b==false)	return "0";
	else				return ".";
}

inline boost::logic::tribool str2trb(std::string s)
{
	boost::to_lower(s);
	if		(s=="1" || s=="true"  || s=="yes" || s=="t" || s=="y" || s=="good" || s=="on"  || s=="positive") return true;
	else if (s=="0" || s=="false" || s=="no"  || s=="f" || s=="n" || s=="bad"  || s=="off" || s=="negative") return false;
	else return boost::logic::tribool(boost::indeterminate);
}

// I use namespace rather than a class because
// 1) there's nothing to hide; 2) one but only one instance is required
// but please use functions as much as possible than to access the data
namespace iData
{
	template <typename T>
	struct TMP { typedef pair<T*,int> _T; }; // template typedef, use as "typedef iData::TMP<double>::_T my_dbl_t;"
	typedef pair<void*,int>		degenerated_t;
	typedef pair<string*,int>	str_t;
	typedef pair<double*,int>	dbl_t;
	typedef pair<tribool*,int>	trb_t;
	//Above: used pointers so that then can be degenerated to void*
	
	// If added a var: change vNum and v_ID (famdb_ind.cpp)
	enum vNum { STR=26, NUM=43, TRB=53 }; // =MAX+1 (consider to use tuple to remove this line)
	// numbers must be unique & increase by 1 & starts from 0, because they are used by vID, vInText, vOuText.
	const str_t IID=str_t(NULL, 0); // Individual ID
	const str_t UID=str_t(NULL, 1); // Unique ID
	const str_t POP=str_t(NULL, 2); // Population
	const str_t	PID=str_t(NULL, 3);	// redundant data
	const str_t	FTH=str_t(NULL, 4);	// redundant data
	const str_t	MTH=str_t(NULL, 5);	// redundant data
	const str_t	REL=str_t(NULL, 6);	// use for output only
	const str_t	L5D=str_t(NULL, 7);	// use for output only
	const str_t	L5N=str_t(NULL, 8);	// use for output only
	const str_t	L6D=str_t(NULL, 9);	// use for output only
	const str_t	L6N=str_t(NULL,10);	// use for output only
	const str_t	L7D=str_t(NULL,11);	// use for output only
	const str_t	L7N=str_t(NULL,12);	// use for output only
	const str_t	L8D=str_t(NULL,13);	// use for output only
	const str_t	L8N=str_t(NULL,14);	// use for output only
	const str_t	L9D=str_t(NULL,15);	// use for output only
	const str_t	L9N=str_t(NULL,16);	// use for output only
	const str_t DNP=str_t(NULL,17); // downcoded PID
	const str_t DNI=str_t(NULL,18); // downcoded IID
	const str_t DNU=str_t(NULL,19); // downcoded UID
	const str_t DNF=str_t(NULL,20); // downcoded FTH
	const str_t DNM=str_t(NULL,21); // downcoded MTH
	const str_t MZT=str_t(NULL,22); // monozygosity twin
	const str_t GTP=str_t(NULL,23); // genotype
	const str_t ALT=str_t(NULL,24); // alternative name
	const str_t CMT=str_t(NULL,25); // comment
	
	const dbl_t SEX=dbl_t(NULL,26);	// sex (1/2)
	const dbl_t AFF=dbl_t(NULL,27);	// affection status (1/2)
	const dbl_t LIA=dbl_t(NULL,28);	// liability class (1+)
	const dbl_t PRB=dbl_t(NULL,29);	// proband status (0/1)
	const dbl_t AGE=dbl_t(NULL,30);	// age (>0)
	const dbl_t YOB=dbl_t(NULL,31);	// year of birth (1-9998)
	const dbl_t GTN=dbl_t(NULL,32);	// genotype counts (0+)
	const dbl_t GTE=dbl_t(NULL,33);	// genotype errors (0+)
	const dbl_t CLT=dbl_t(NULL,34);	// cluster ID (1+)
	const dbl_t INB=dbl_t(NULL,35);	// inbreeding coefficient
	const dbl_t GEN=dbl_t(NULL,36);	// generation ID
	const dbl_t MXG=dbl_t(NULL,37);	// max number of generations among descendants
	const dbl_t DES=dbl_t(NULL,38);	// tot number of descendants
	const dbl_t AL1=dbl_t(NULL,39);	// allele 1
	const dbl_t AL2=dbl_t(NULL,40);	// allele 2
	const dbl_t AFD=dbl_t(NULL,41);	// number of affected first degree relatives (nan for unaffected ind)
	const dbl_t IWT=dbl_t(NULL,42);	// individual weight for association test
	
	const trb_t INF=trb_t(NULL,43);	// informative
	const trb_t SEP=trb_t(NULL,44);	// separated
	const trb_t VST=trb_t(NULL,45);	// visited
	const trb_t B01=trb_t(NULL,46);	// spared trb, not used
	const trb_t RCG=trb_t(NULL,47);	// columns recognized
	const trb_t OTH=trb_t(NULL,48);	// columns not recognized
	const trb_t DUM=trb_t(NULL,49);	// is dummy spouse created by the program
	const trb_t SXD=trb_t(NULL,50);	// sex is deduced
	const trb_t YBI=trb_t(NULL,51);	// YoB is inferred
	const trb_t ABP=trb_t(NULL,52);	// Add By Program
	
	extern std::map<int, string>	vOuText;	// var name for output
	extern const std::string		default_vname;
	extern std::vector<string>		v_ID;
	
	template <typename T>
	inline string& vText(const pair<T*,int>& var_label) { return vOuText[var_label.second]; }

	template <typename T>
	inline int vCode(const pair<T*,int>& var_label) { return var_label.second; }

	template <typename T1, typename T2>
	inline bool is_equal(const pair<T1*,int>& label1, const pair<T2*,int>& label2) { return label1.second==label2.second; }

	template <typename T1, typename T2>
	inline bool vec_contain(const vector< pair<T1*,int> >& vec, const pair<T2*,int>& label) {
		typename vector< pair<T1*,int> >::const_iterator it;
		for (it=vec.begin(); it!=vec.end(); ++it) if (is_equal(*it, label)) return true;
//		for (each_element(vec,it)) if (is_equal(*it, label)) return true;
		return false;
	}

	inline bool is_str_t(const degenerated_t& var_label) { if (var_label.second < 0  ) return false; return var_label.second < STR; }
	inline bool is_dbl_t(const degenerated_t& var_label) { if (var_label.second < STR) return false; return var_label.second < NUM; }
	inline bool is_trb_t(const degenerated_t& var_label) { if (var_label.second < NUM) return false; return var_label.second < TRB; }
	
	string convert(const iData::degenerated_t& var_label, const string& input);
	void read_conversion_map(const string& instruction);
	void get_var_name_db(const string& filename);
	bool get_var_name(string option);		// read program option --XXX=yyyyyy
	string expand_var_names(string input);	// expand l6 => pid,iid,fa,mo,sx,af
	degenerated_t get_id_from(const string& input);
	void set_ignore_text(const iData::degenerated_t& var_label, const string& input);
	bool  to_ignore_text(const iData::degenerated_t& var_label, const string& input);
	void set_known_var(const string& input);
}

class Individual
{
public:
	struct Relationship {
		char	p_m;			// relation to proband: Paternal Maternal Both None Unknown
		int		code;			// relation to proband: code
		string	relation;		// relation to proband: description
		Relationship():p_m('U'),code(999),relation("Unknown"){}
		void clear() { p_m='U'; code=999; relation="Unknown"; }
	};
	
	Family		*					fam;	// pointer to family struct
	Individual	*					pfa;	// pointer to fater
	Individual	*					pmo;	// pointer to mother
	IndMap							spouses;// spouses
	IndMap							nextgen;// children
	std::map<string,Relationship>	rrs;	// <relates_to_who, relationship>
	Individual():fam(NULL),pfa(NULL),pmo(NULL) { }
	// type safe & fast setting data
	void		set(const iData::str_t& var_label, const string& value)	{ str[iData::vCode(var_label)]=value; }
	void		set(const iData::dbl_t& var_label, const double	 value)	{ num[iData::vCode(var_label)]=value; }
	void		set(const iData::trb_t& var_label, const tribool value) { trb[iData::vCode(var_label)]=value; }
	void		set(const iData::trb_t& var_label, const bool	 value)	{ trb[iData::vCode(var_label)]=value; }
	// type safe & fast getting data. return reference, so can change data.
	string&		get(const iData::str_t& var_label) { return str[iData::vCode(var_label)]; }
	double&		get(const iData::dbl_t& var_label) { return num[iData::vCode(var_label)]; }
	tribool&	get(const iData::trb_t& var_label) { return trb[iData::vCode(var_label)]; }
	string&		operator[](const iData::str_t& var_label) { return str[iData::vCode(var_label)]; }
	double&		operator[](const iData::dbl_t& var_label) { return num[iData::vCode(var_label)]; }
	tribool&	operator[](const iData::trb_t& var_label) { return trb[iData::vCode(var_label)]; }
	// brute-force input data, good for reading data from a file
	void	 bf_set(const iData::degenerated_t& var_label, const string& input);
	// brute-force output data, good for writing data to a file
	string	 bf_get(const iData::degenerated_t& var_label);
	void copy_data( const Individual& orig );
	void remove( iData::degenerated_t var_label );
	
private:
	Individual& operator=(const Individual& orig);	// forbidden
	Individual(const Individual& othr);				// forbidden

	// use map to store data has several advantages:
	// 1) good for sparse data; 2) easy to copy data - no need to tell what to copy; 3) natural for missing values
	// bellow -- int = iData::vCode(var_label)
	std::map<int,string>		str;	// string variables
	std::map<int,double>		num;	// numeric variables
	std::map<int,tribool>		trb;	// tribool variables, default value is false
public:
	// bellow -- int = in_file.col_num (0-based)
	std::map<int,string>		xtr;	// extra columns' data. <col, data>
};

#endif
