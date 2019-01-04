// body file: libfbj_famdb.cpp

#ifndef FamilyDB_class
#define FamilyDB_class

#include "libfbj_famdb_fam.hpp"
#include <fstream>
#include <map>
#include <vector>
#include <list>

class FamilyDB {

public:
	FamilyDB();
	~FamilyDB();
	void read_arguments(vector<string>& srce_opt);	// program.arg
	bool execute_command(vector<string>& tokens);	// one command in the form of a vector
	
private:
	typedef std::map<std::string,int> ID_Count;
	typedef void (FamilyDB::*TravFuncPtrType)(Individual *, void * UserData);
	
	// about individual
	Individual* find_ind(std::string fam_name, std::string ind_name);
	Individual* add_ind(Individual& src, bool show_2nd_input_wrn=true);		// create new ind and return its address, copy data from src. use VST
	void		del_inds_by_status(iData::trb_t status, boost::logic::tribool value);
	void		set_inds_status(Family * f, iData::trb_t status_fld, boost::logic::tribool value);
	void		set_inds_values(Family * f, iData::dbl_t values_fld, double value);
	void		set_inds_status(iData::trb_t status_fld, boost::logic::tribool value);
	void		set_inds_values(iData::dbl_t values_fld, double value);
	double		max_inds_values(Family * f, iData::dbl_t values_fld);
	double		mnp_inds_values(Family * f, iData::dbl_t values_fld);
	double		max_inds_values(iData::dbl_t values_fld);
	double		mnp_inds_values(iData::dbl_t values_fld);
	void		set_ind_ID(Individual& ind);	// check PID IID UID FTH MTH & create IID/UID. Call it bef add_ind()!
//	void		set_ind_def(Individual& ind);	// add default values. Call it before add_ind()!
	class		bad_ind_ID { public: const char *ShowReason() const { return "ID problem."; } }; // throw by set_ind_ID
	void		downcode_IDs();					// write to DNP,DNI; never change PID,IID,UID,FTH,MTH
	double		num_sibs(Individual* ind);		// number of siblings: fullsib=1 halfsib=0.5
	void		founderize(Individual* ind);	// remove the links with parents, could create floating individuals
	void		make_parents(Individual* ind, Individual* fth, Individual* mth); // all arguments cannot be NULL
	bool		is_affected(Individual* ind);	// return whether AFF==2 (AFF=0/1/2)
	bool		is_proband(Individual* ind);	// return whether PRB!=0 (PRB=0/1, but allow other numbers so that eg treat AFF as PRB)
	void		find_1st_degr_relatives(Individual* ind, IndMap& dst);	// 1st degree relatives, exclude half sibs
	void		fast_add(const std::string& parameters);
	void		_auto_add_parents(Individual* ind); // called by fast_add()
	void		_fast_add_set(vector<string>& arguments, Individual* ind); // called by fast_add()
	
	// input
	void _check_fields(std::map<int,iData::degenerated_t>& fld_4_col, std::string& d_pid);
	void read_prefile(std::string filename, bool has_title=true, const std::string& instruction=std::string());
	
	// checking common errors
	void read_family_errors();
	void check_questionable_IDs();	// Must call after read_family_errors();
	void check_mzt();
	void check_prb();
	void check_yob();
	void check_aff();
	void check_lia();
	void fill_in_missing_yob(int default_yob_for_probands);
	
	// checking loops
	void check_loop(bool to_breakup=false);
	int  _check_loop(bool to_breakup=false);
	void _find_ancestors(Individual* ind, ID_Count& ancestors);
	void _find_founders(Individual* ind, ID_Count& ancestors);
	void _find_ancestors(Individual* ind, std::map<Individual*,int>& ancestors);
	void _find_founders(Individual* ind, std::map<Individual*,int>& ancestors);
	
	// checking circles
	void check_circles(bool to_breakup=false);
	int  _check_circles(bool to_breakup=false);

	// checking self-ancestors
	void check_self_ancestor();
	void _check_self_anc(Individual* ind, ID_Count& ancestors);

	// checking separated individuals
	void check_separation();
		
	// cluster candidate INDs within N degr rel, assign cluster_ID (1+, overlapping) and print (omit 0s)
	void cal_clustering(const string& cluster_uid_f, bool cluster_affected, bool to_write);
	
	// calculate informativeness
	void cal_informativeness();
	void _is_informative(Individual * ind);

	// calculate generation number for each IND
	void cal_generation();
	void cal_generation_fast();
	
	// calculate relationship between proband and other persons
	void read_IND(const string& filename, std::map<Family*, std::vector<Individual*> >& data);
	void cal_kinship_file(const std::string& filename, bool quietly=false);
	void cal_kinship_all(bool quietly=false);
	void cal_MRCA_file(const std::string& filename);
	void cal_relationship();
	void cal_relationship_among_aff();
	void find_ancestors(std::map<std::string,std::vector<int> >* l,std::vector<int>* v,Individual* i,int g);
	void find_common_ancestors(Family* f,std::string prob,std::string subj);
	void find_nearest_common_ancestors(Family* f,std::string prob,std::string subj);
	void cal_relationship(Family* f,std::string prob,std::string subj);
	void find_founders_of(const string& filename);
	void find_connections(Individual * ind, std::set<std::string>& tca, std::set<std::string>& connectors);
	void cal_OC_file(const std::string& filename);

	// calculate number of affecte FDR
	void cal_aff_fdr();
	
	// gene drop and individual weight
	void gene_drop(bool quietly);
	void ind_weight();
	
	// parent of origin analysis
	void single_aff_par();
	
	// for traverse
	void traverse_IS(TravFuncPtrType fptr, void * UserData);			// same as input sequence [excl created / unread ones]
	void traverse_up(TravFuncPtrType fptr, void * UserData);			// ancestors visited ahead of offsprings
	void traverse_dn(TravFuncPtrType fptr, void * UserData);			// ancestors could be visited after offsprings
	void _traverse_up_subfunc(Individual* ind, TravFuncPtrType fptr, void * UserData);	// recursively called by traverse_up
	void _traverse_dn_subfunc(Individual* ind, TravFuncPtrType fptr, void * UserData);	// recursively called by traverse_dn
	void _print_line(Individual* ind, void * UserData);					// function to print a line
	void _print_header();												// function to print header

	// for output
	void set_output_fields(const std::string& instruction, bool clear_previous);// called by write_output if outFields.empty, or called by user
	void write_output(bool rows_in_input_seq=false);					// write pedigree file to program.outf
	void write_summary();												// summary of ped: 
	
	// data
	std::stringstream	add_ind_logf_;	// fam err log by add_ind(), will be corrected and cleared by read_family_errors()
	std::string			string_no_id_;	// ID string for "nobody", default is 0
	std::string			id_delimiter_;	// delimiter between PID and IID, default is ::
	std::string			given_iid_pfx;	// given_iid prefix
	bool to_remove_sep_id;				// default = no
	bool to_show_warnings;				// report warnings
	bool to_report_single;				// report warnings for single-parent problem, default=yes
	bool to_report_gender;				// report warnings for gender errors, default=yes
	bool to_report_re_inp;				// report warnings for re-input individual or pedigree, default=yes
	bool to_report_subped;				// report warnings for sub-peds (more than 1 disconnected families within 1 PID)
	bool to_report_dupUID;				// report warnings for Identical UIDs
	bool to_aggregate_con;				// aggregate connections -- do not disconnect parents if 2nd input has no parents, def=no
	bool read_uid_not_iid;				// input contains UIDs but not IIDs, def=no
	bool unaff_sib_is_inf;				// in finding informative persons, unaffected sibs of an aff are informative
	bool always_do_create;				// in fixing single_parent problems, always create a dummy parent. Will be true when reading pedfile from stdin.
	bool wr_alternativeID;				// write alternative ID instead of IID
	bool iWt_useCasesOnly;				// when calculating individual weights, use cases only
	bool iWt_SumTo1_eaPed;				// when calculating individual weights, the sum of weights in each pedigree is 1, so that I can do rare-variant association test and linkage analysis
	int	clt_degr_;						// cluster people within # degree relatives, def=3. This value affect line sorting.
	int mo_age_lb, mo_age_ub, fa_age_lb, fa_age_ub;
	std::vector<iData::degenerated_t>	outFields;	// sequence of output fields
	std::vector<iData::degenerated_t>	inpFields;	// sequence of input fields (known fields only)
	FamMap								family_db;	// pid -> Family*
	IndMap								uniqid_db;	// uid -> Individual*
	int									given_pid;	// given PIDs when input file has no PID field
	int									given_uid;	// given UIDs when input line has empty UID/IID
	int									given_iid;	// given IIDs when adding a new person
	std::map<int, std::string>			xtrColTtl;	// extra columns' title. field_number(0-based) -> title
	std::vector< std::pair< std::string,std::vector<std::string> > >	rows_input_seq;	// input sequence of families and individuals
};

#endif
