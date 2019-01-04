/*
=====limitations of informativeness check=====

This program can set the following individuals to be uninformative:
- unaffected ungenotyped founders with only 1 informative child
- unaffected ungenotyped individuals without any informative child

This will lead to a problem when
- a genotyped founder w/ 0 informative child will not be removed, creating an isolated person. Need rm_seq afterwards.
- unaffected ungenotyped nonfounders with only 1 informative child, and his/her parent is actually uninformative or not
  connecting to other informative individuals, hense he/she is actually uninformative but not exluded.
  Currently, I exclude two persons (5068-26 and 4801-17) so that it can still run correctly.

===== relationship calculation ======

algorithm taken from http://en.wikipedia.org/wiki/Granduncle 

	file<<"\n\n\n\nThis is the output file of the program <relation>.\n"
		<<"In this file, 3 fields are added:\n"
		<<"      rcode   : code of relationship\n"
		<<"      p_m     : paternal/maternal\n"
		<<"      relation: description of relationship\n\n"
		<<"The coding of rcode:                                 1st      2nd    ..  5th\n"
		<<"                                                     cousin   cousin ..  cousin\n"
		<<"                                                     n times  n times..  n times\n"
		<<"                                                     removed  removed..  removed\n"
		<<" 6: Great^4 Grand Parent  16: great^5 uncle/aunt     26: 6    36: 6  ..  66: 6\n"
		<<" 5: Great^3 Grand Parent  15: great^4 uncle/aunt     25: 5    35: 5  ..  65: 5\n"
		<<" 4: Great^2 Grand Parent  14: great^3 uncle/aunt     24: 4    34: 4  ..  64: 4\n"
		<<" 3: Great Grand Parent    13: great^2 uncle/aunt     23: 3    33: 3  ..  63: 3\n"
		<<" 2: Grand Parent          12: great uncle/aunt       22: 2    32: 2  ..  62: 2\n"
		<<" 1: parent                11: uncle/Aunt             21: 1    31: 1  ..  61: 1\n"
		<<" 0: proband               10: Siblings               20: -    30: -  ..  60: -\n"
		<<"-1: Child                -11: Niece/Nephew          -21: 1   -31: 1  .. -61: 1\n"
		<<"-2: Grand Child          -12: great niece/nephew    -22: 2   -32: 2  .. -62: 2\n"
		<<"-3: Great Grand child    -13: great^2 niece/nephew  -23: 3   -33: 3  .. -63: 3\n"
		<<"-4: Great^2 Grand child  -14: great^3 niece/nephew  -24: 4   -34: 4  .. -64: 4\n"
		<<"-5: Great^3 Grand child  -15: great^4 niece/nephew  -25: 5   -35: 5  .. -65: 5\n"
		<<"-6: Great^4 Grand child  -16: great^5 niece/nephew  -26: 6   -36: 6  .. -66: 6\n\n"
		<<"  any halfxxx is 2xx or -2xx\n"
		<<"  other blood relatives are 8\n"
		<<"  unrelated individuals are 9\n"
		<<"  unknown relationships are 9999 (should not appear)\n"
		<<"\n"
		<<"The coding of p_m \n"
		<<"  P : paternal\n"
		<<"  M : maternal\n"
		<<"  B : Both paternal and maternal\n"
		<<"  N : none (unrelated individuals)\n"
		<<"  - : proband herself or her nextgens\n"
		<<"  U : Unknown (should not appear)\n"
		<<"\n"
		<<"Be careful: \n"
		<<"  1) If there're >2 common ancestors, and the link from proband to subject have different number of meiosis,\n"
		<<"     then only the nearest one is used for relationship calculation, and underestimate the kin-coefficient.\n";

 ===== options not showed by --help ======
 --iWt-sum1        Individual weights in each pedigree sum to one.
 */

#include <sstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <queue>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_rng.hpp>
#include "libfbj_famdb.hpp"
#include "victor_par.hpp"

using namespace std;

namespace {

	const int kMAN=1;
	const int kWMN=2;
	const int kAFF=2;
	const int kUNK=0;
	
	struct str_relation // relationship between 2 persons
	{
		std::string famname,probname,subjname,relation;
		std::map<std::string,std::vector<int> > a_prob, a_subj; // g from prob, subj
		std::map<std::string,std::pair<int,int> > comanc;		// g from subj, prob
		int		min_g;	// minimum generations (number of miosis) linking the two person
		int		g_prb;	// generation from nearest common ancestor to proband
		int		g_sbj;	// generation from nearest common ancestor to subject
		int		rcode;	// relation to proband: code (see above)
		char	p_m;	// relation to proband: Paternal Maternal Both None Unknown
	}strr;
}

extern logger lns;

FamilyDB::FamilyDB()
{
	++program.nudge; // license checking
	string_no_id_="0";
	id_delimiter_="::";
	given_iid_pfx="DIABP"; // DIABP=Dummy ID Added By Pedpro. Previously "A" but easily get duplicated with a real IndID
	to_remove_sep_id=false;
	to_show_warnings=true;
	to_report_single=true;
	to_report_gender=true;
	to_report_re_inp=true;
	to_report_subped=true;
	to_report_dupUID=true;
	to_aggregate_con=false;
	read_uid_not_iid=false;
	unaff_sib_is_inf=false;
	always_do_create=false;
	wr_alternativeID=false;
	iWt_useCasesOnly=false;
	iWt_SumTo1_eaPed=false;
	mo_age_lb = 8;
	mo_age_ub = 70;
	fa_age_lb = 12;
	fa_age_ub = 90;
	given_pid = 0;
	given_uid = 0;
	given_iid = 0;
	clt_degr_ = 3;
	iData::get_var_name_db(program.prefix()+".var_name.txt");
	iData::set_ignore_text(iData::MZT,"");
	iData::read_conversion_map("sex:m/M/male/Male=1,f/F/female/Female=2;aff:affected/aff/y/Y/yes/Yes/A=2,unaffected/unaff/n/N/no/No/U=1");
}

FamilyDB::~FamilyDB()
{
}

inline void FamilyDB::founderize(Individual* ind)
{
	if (ind->pfa == NULL) exit_error("\nERROR: trying to founderize a founder");
	
	string& id=ind->get(iData::IID);
	string& fa=ind->pfa->get(iData::IID);
	string& mo=ind->pmo->get(iData::IID);
	
	ind->pfa->nextgen.erase(id);
	ind->pmo->nextgen.erase(id);
	int num_kids=0;
	for (each_element_const(ind->pfa->nextgen,it)) if (it->second->pmo->get(iData::IID)==mo) ++num_kids;
	if  (num_kids==0) { ind->pfa->spouses.erase(mo); ind->pmo->spouses.erase(fa); }
	ind->pfa=NULL;
	ind->pmo=NULL;
	ind->set(iData::FTH,string_no_id_);
	ind->set(iData::MTH,string_no_id_);
	ind->fam->nonfounders.erase(id);
	ind->fam->founders[id]=ind;	
}

inline void FamilyDB::make_parents(Individual* ind, Individual* fth, Individual* mth)
{
	string& id=ind->get(iData::IID);
	string& fa=fth->get(iData::IID);
	string& mo=mth->get(iData::IID);
	
	ind->pfa=fth;
	ind->pmo=mth;
	fth->nextgen[id]=ind;
	mth->nextgen[id]=ind;
	fth->spouses[mo]=mth;
	mth->spouses[fa]=fth;
	ind->set(iData::FTH,fa);
	ind->set(iData::MTH,mo);
	ind->fam->founders.erase(id);
	ind->fam->nonfounders[id]=ind;
}

// aft call -- iData::ABP indicates whether the person is not read but added
Individual* FamilyDB::add_ind(Individual& inds, bool show_2nd_input_wrn)
{
	string& pid		= inds.get(iData::PID);
	string& iid		= inds.get(iData::IID);
	string& father	= inds.get(iData::FTH);
	string& mother	= inds.get(iData::MTH);
	
	// check error 
	// I remove the parts to solve self-parent and same-parents, because it transforms the problem to another one -- single parent,
	// which could be common and even not a problem at all. To make it more clear that a more severe problem has happened, I would
	// like to leave it as it is, report an error and let the user to solve it manually.
	if (father==iid || mother==iid)
	{
		lns<<showe<<"fam "<<pid<<" person "<<iid<<" is self-parent, founderized and program continues."<<flush_logger;
		father=string_no_id_;
		mother=string_no_id_;
		inds.set(iData::FTH,string_no_id_);
		inds.set(iData::MTH,string_no_id_);
	}
	if ( father!=string_no_id_ &&  father==mother )
	{
		lns<<showe<<"fam "<<pid<<" person "<<iid<<" has same-parents, founderized and program continues."<<flush_logger;
		father=string_no_id_;
		mother=string_no_id_;
		inds.set(iData::FTH,string_no_id_);
		inds.set(iData::MTH,string_no_id_);
	}
	if ((father!=string_no_id_&&mother==string_no_id_)||(father==string_no_id_&&mother!=string_no_id_))
	{
		if (to_show_warnings && to_report_single)
			lns<<writew<<"fam "<<pid<<" person "<<iid<<" has a single-parent, founderized and program continues."<<flush_logger;
		add_ind_logf_<<pid<<'\t'<<iid<<'\t'<<father<<'\t'<<mother<<'\t'<<inds.get(iData::SEX)<<'\t'<<inds.get(iData::AFF)<<"\tsingle"<<'\n';
		father=string_no_id_;
		mother=string_no_id_;
		inds.set(iData::FTH,string_no_id_);
		inds.set(iData::MTH,string_no_id_);
	}
	
	//find family
	Family * f;
	if (!exist_element(family_db,pid))
	{
		f=new Family;
		if (!f) exit_lackMemory("new family");
		f->pid=pid;
		family_db[pid]=f;
	}
	else
	{
		f=family_db[pid];
	}
	
	// add father if necessary
	if (father!=string_no_id_)
	{
		if (!exist_element(f->members,father))
		{
			Individual dummy_ind;
			dummy_ind.set(iData::SEX,kMAN);
			dummy_ind.set(iData::IID,father);
			dummy_ind.set(iData::PID,pid);
			set_ind_ID(dummy_ind);
			add_ind(dummy_ind)->set(iData::ABP,true);
		}
		else
		{
			if (f->members[father]->get(iData::SEX)!=kMAN)
			{
				if (f->members[father]->get(iData::SXD)==true) 
					lns<<showe<<"fam "<<pid<<" person "<<father<<" : ambiguous-gender."<<flush_logger;
				else if (to_show_warnings && to_report_gender)
					lns<<writew<<"fam "<<pid<<" person "<<iid<<" father ("<<father<<") wrong-gender / previously unknown, corrected."<<flush_logger;
				f->members[father]->set(iData::SEX,kMAN);
			}
		}
		f->members[father]->set(iData::SXD,true);
	}

	// add mother if necessary
	if (mother!=string_no_id_)
	{
		if (!exist_element(f->members,mother))
		{
			Individual dummy_ind;
			dummy_ind.set(iData::SEX,kWMN);
			dummy_ind.set(iData::IID,mother);
			dummy_ind.set(iData::PID,pid);
			set_ind_ID(dummy_ind);
			add_ind(dummy_ind)->set(iData::ABP,true);
		}
		else
		{
			if (f->members[mother]->get(iData::SEX)!=kWMN)
			{
				if (f->members[mother]->get(iData::SXD)==true)
					lns<<showe<<"fam "<<pid<<" person "<<mother<<" : ambiguous-gender."<<flush_logger;
				else if (to_show_warnings && to_report_gender)
					lns<<writew<<"fam "<<pid<<" person "<<iid<<" mother ("<<mother<<") wrong-gender / previously unknown, corrected."<<flush_logger;
				f->members[mother]->set(iData::SEX,kWMN);
			}
		}
		f->members[mother]->set(iData::SXD,true);
	}

	// find self
	Individual *i;
	if (!exist_element(f->members,iid))
	{
		// check wrong letter case
		string lciid = boost::to_lower_copy(iid);
		if (exist_element(f->lc_iids,lciid)) { if (to_show_warnings) lns<<showw<<"Fam "<<pid<<" has >1 IIDs with different letter case of "<<lciid<<flush_logger; }
		else f->lc_iids.insert(lciid);

		// begin to add i into f
		i=new Individual;
		if (!i) exit_lackMemory("new individual");
		i->fam=f;
		f->members[iid]=i;
		i->copy_data(inds); // make_parents() needs the IID
		if (father!=string_no_id_)
			make_parents(i,f->members[father],f->members[mother]);
		else
			f->founders[iid]=i;
	}
	else
	{
		i=f->members[iid];
		if (father!=string_no_id_)
		{
			if (i->pfa)
			{
				if (i->pfa->get(iData::IID)!=father || i->pmo->get(iData::IID)!=mother)
				{
					if (i->get(iData::ABP)!=true && show_2nd_input_wrn) 
						lns<<writew<<"fam "<<pid<<" person "<<iid<<" update: change parents ("<<i->pfa->get(iData::IID)<<" & "<<i->pmo->get(iData::IID)<<") to ("<<father<<" & "<<mother<<")"<<flush_logger;
					founderize(i);
					make_parents(i,f->members[father],f->members[mother]);
				}
			}
			else
			{
				if (i->get(iData::ABP)!=true && show_2nd_input_wrn) 
					lns<<writew<<"fam "<<pid<<" person "<<iid<<" update: add parents "<<father<<" & "<<mother<<flush_logger;
				make_parents(i,f->members[father],f->members[mother]);
			}
		}
		else
		{
			if (i->pfa)
			{
				if (i->get(iData::ABP)!=true)
				{
					if (to_aggregate_con)
					{
						if (show_2nd_input_wrn)
							lns<<writew<<"fam "<<pid<<" person "<<iid<<" 2nd input has no parents, not updated"<<flush_logger;
						inds.remove(iData::FTH);
						inds.remove(iData::MTH);
						inds.remove(iData::UID);
					}
					else
					{
						if (show_2nd_input_wrn)
							lns<<writew<<"fam "<<pid<<" person "<<iid<<" update: disconnect from parents "<<i->pfa->get(iData::IID)<<" & "<<i->pmo->get(iData::IID)<<flush_logger;
						founderize(i);
					}
				}
			}
		}
		
		if (inds.get(iData::SEX)!=kUNK)
		{
			if (i->get(iData::SEX)!=kUNK && i->get(iData::SEX)!=inds.get(iData::SEX))
			{
				if (i->get(iData::SXD)==true)
				{
					if (to_show_warnings && to_report_gender)
						lns<<writew<<"fam "<<pid<<" person "<<iid<<" wrong-gender, corrected by deduction from fa/mo status."<<flush_logger;
					inds.remove(iData::SEX); // remove data to avoid copy_data() from copying it
				}
				else
				{
					if (i->get(iData::ABP)!=true && show_2nd_input_wrn) 
						lns<<writew<<"fam "<<pid<<" person "<<iid<<" update: sex different, use new."<<flush_logger;
				}
			}
		}
		else
			inds.remove(iData::SEX); // if (iData::SXD)==true, then it makes no sense to set iData::SEX = kUNK

		i->copy_data(inds);

	}
	i->set(iData::ABP,false);
	
	// set uid
	string uid_core = pid + id_delimiter_ + iid;
	if (i->get(iData::UID)==string_no_id_) i->set(iData::UID,uid_core);
	for (int j=1; ; ++j)
	{
		if (j==100) { lns<<showe<<"Create UIDs for fam "<<pid<<" ind "<<iid<<" failed!"<<fatal; }
		IndMap::iterator it = uniqid_db.find(i->get(iData::UID));
		if (it==uniqid_db.end()) break;
		string& other_pid=it->second->get(iData::PID);
		string& other_iid=it->second->get(iData::IID);
		if (pid==other_pid && iid==other_iid) break;
		if (to_show_warnings && to_report_dupUID)
			lns<<showw<<"Identical UIDs: fam "<<pid<<" ind "<<iid<<" == fam "<<other_pid<<" ind "<<other_iid<<flush_logger;
		i->set(iData::UID, uid_core+"_"+itos_format(j, 'r', 2, '0', 10));
	}
	uniqid_db[i->get(iData::UID)]=i;

	// return Individual*
	return i;
}

void FamilyDB::set_inds_status(Family * f, iData::trb_t status_fld, boost::logic::tribool value)
{
	for (each_element(f->members,iti))
		iti->second->set(status_fld,value);
}

void FamilyDB::set_inds_values(Family * f, iData::dbl_t values_fld, double value)
{
	for (each_element(f->members,iti))
		iti->second->set(values_fld,value);
}

void FamilyDB::set_inds_status(iData::trb_t status_fld, boost::logic::tribool value)
{
	for (each_element(family_db,itf))
		set_inds_status(itf->second,status_fld,value);
}

void FamilyDB::set_inds_values(iData::dbl_t values_fld, double value)
{
	for (each_element(family_db,itf))
		set_inds_values(itf->second,values_fld,value);
}

double FamilyDB::max_inds_values(Family * f, iData::dbl_t values_fld)
{
	double overall = -DBL_MAX;
	for (each_element(f->members,iti))
	{
		double r = iti->second->get(values_fld);
		if (r>overall) overall=r;
	}
	return overall;
}

double FamilyDB::mnp_inds_values(Family * f, iData::dbl_t values_fld)
{
	double overall = DBL_MAX;
	for (each_element(f->members,iti))
	{
		double r = iti->second->get(values_fld);
		if (r>0 && r<overall) overall=r;
	}
	return overall;
}

double FamilyDB::max_inds_values(iData::dbl_t values_fld)
{
	double overall = -DBL_MAX;
	for (each_element(family_db,itf))
	{
		double r = max_inds_values(itf->second,values_fld);
		if (r>overall) overall=r;
	}
	return overall;
}

double FamilyDB::mnp_inds_values(iData::dbl_t values_fld)
{
	double overall = DBL_MAX;
	for (each_element(family_db,itf))
	{
		double r = mnp_inds_values(itf->second,values_fld);
		if (r>0 && r<overall) overall=r;
	}
	return overall;
}

// Call it bef add_ind()!
void FamilyDB::set_ind_ID(Individual& ind)
{
	if (ind.get(iData::PID).empty() || ind.get(iData::PID)==".") ind.set(iData::PID,string_no_id_);
	if (ind.get(iData::FTH).empty() || ind.get(iData::FTH)==".") ind.set(iData::FTH,string_no_id_);
	if (ind.get(iData::MTH).empty() || ind.get(iData::MTH)==".") ind.set(iData::MTH,string_no_id_);
	if (ind.get(iData::IID).empty() || ind.get(iData::IID)==".") ind.set(iData::IID,string_no_id_);
	if (ind.get(iData::UID).empty() || ind.get(iData::UID)==".") ind.set(iData::UID,string_no_id_);
	if (ind.get(iData::ALT).empty() || ind.get(iData::ALT)==".") ind.set(iData::ALT,string_no_id_);
	if (ind.get(iData::PID)==string_no_id_) throw bad_ind_ID();
	if (ind.get(iData::IID)==string_no_id_ && ind.get(iData::UID)!=string_no_id_)	ind.set(iData::IID,ind.get(iData::UID));
	if (ind.get(iData::IID)!=string_no_id_ && read_uid_not_iid)						ind.set(iData::UID,ind.get(iData::IID));
	if (ind.get(iData::IID)==string_no_id_ && ind.get(iData::UID)==string_no_id_)
	{
		string g_uid = "_given_uid_"+itos(++given_uid);
		ind.set(iData::IID,g_uid);
		ind.set(iData::UID,g_uid);
	}
	// at this point, both pid and iid is fine. but uid could be no_id.
}

// return correctness
void FamilyDB::_check_fields(map<int,iData::degenerated_t>& fld_4_col, string& g_pid)
{
	set<iData::degenerated_t> known;
	for (each_element(fld_4_col,itf)) {
		if (iData::is_equal(itf->second,iData::DNP)) itf->second=iData::PID;
		if (iData::is_equal(itf->second,iData::DNI)) itf->second=iData::IID;
		if (iData::is_equal(itf->second,iData::DNU)) itf->second=iData::UID;
		if (iData::is_equal(itf->second,iData::DNF)) itf->second=iData::FTH;
		if (iData::is_equal(itf->second,iData::DNM)) itf->second=iData::MTH;
		if (exist_element(known, itf->second)) exit_error("Multiple "+iData::vText(itf->second));
		else known.insert(itf->second);
		if (iData::is_equal(itf->second,iData::L5N) || iData::is_equal(itf->second,iData::L5D) ||
			iData::is_equal(itf->second,iData::L6N) || iData::is_equal(itf->second,iData::L6D) ||
			iData::is_equal(itf->second,iData::L7N) || iData::is_equal(itf->second,iData::L7D) ||
			iData::is_equal(itf->second,iData::L8N) || iData::is_equal(itf->second,iData::L8D) ||
			iData::is_equal(itf->second,iData::L9N) || iData::is_equal(itf->second,iData::L9D) ||
			iData::is_equal(itf->second,iData::REL) || iData::is_equal(itf->second,iData::OTH) )
			exit_error("\nERROR: These fields are forbidden in an input pedigree file: rel,oth.");
	}
	if (!exist_element(known,iData::PID)) { g_pid = "_given_pid_"+itos(++given_pid); }
	if (!exist_element(known,iData::FTH)) { lns<<showe<<"No Father_ID column. Stopped reading file."<<fatal; }
	if (!exist_element(known,iData::MTH)) { lns<<showe<<"No Mother_ID column. Stopped reading file."<<fatal; }
	if (!exist_element(known,iData::IID) && !exist_element(known,iData::UID)) { lns<<showe<<"No IIDs or UIDs column. Stopped reading file."<<fatal; }
	if (!exist_element(known,iData::IID) &&  exist_element(known,iData::UID)) { read_uid_not_iid=true; }
}

void FamilyDB::read_prefile(string filename, bool has_title, const string& instruction) // <string>, true, <empty> by default
{
	if (has_title==true && !instruction.empty()) exit_error("\nERROR: Don't need input field instruction for pedfiles with a title.");
	if (has_title==false && instruction.empty()) exit_error("\nERROR: Need an input field instruction for pedfiles without a title.");
	if (is_stdin(filename)) always_do_create=true;
	
	lns<<showl<<"Read prefile "<<filename<<" ... "<<flush_logger;
	lns.sub();
	
	map<int,iData::degenerated_t> fld_4_col;
	int minCol=0;		// min number of columns required
	string g_pid;		// given PID, used when there's no PID field
	string last_pid;	// PID of last line, used for rows_input_seq
	
	if (!has_title) 
	{
		string expanded=iData::expand_var_names(instruction);
		lns<<showl<<"Input fields: "<<instruction<<flush_logger;
		lns<<showl<<"Expanded fields: "<<expanded<<flush_logger;
		for (int i=0; !expanded.empty(); ++i)
		{
			string s=extract_name(expanded);
			if (s.empty()) break;
			if (!expanded.empty()) extract_char(expanded);
			iData::degenerated_t field = iData::get_id_from(s);
			if (iData::vCode(field)==-1) { xtrColTtl[i]=s; continue; }
			fld_4_col[i]=field;
			inpFields.push_back(field);
		}
		_check_fields(fld_4_col,g_pid);
		minCol=fld_4_col.size()-1;
	}
	
	int skip_sht=0, skip_bad=0, added=0, dup_ind=0, dup_fam=0;
	for (Rows_in_File(inf,filename,0))
	{
		if (inf.RowNumber()==0 && has_title) // title line
		{
			lns<<showl<<"Among "<<inf.NumFields()<<" columns, recognized:";
			for (int i=0;i<inf.NumFields();++i)
			{
				iData::degenerated_t field = iData::get_id_from(inf[i]);
				if (iData::vCode(field)==-1) { xtrColTtl[i]=inf[i]; continue; }
				fld_4_col[i]=field;
				inpFields.push_back(field);
				lns<<" "<<iData::v_ID[iData::vCode(field)]; //lns<<" "<<iData::vText(fld_4_col[i]);
			}
			lns<<flush_logger;
			_check_fields(fld_4_col,g_pid);
			minCol=inf.NumFields()-1;
		}
		else
		{
			// read data
			if (inf.NumFields() <= minCol) { ++skip_sht; continue; }
			Individual inds;
			for (int i=0;i<inf.NumFields();++i)
			{
				if (exist_element(fld_4_col,i))
				{
					try { inds.bf_set(fld_4_col[i],inf[i]); }
					catch (boost::bad_lexical_cast &)
					{
						// if (iData::is_equal(fld_4_col[i],iData::AGE)) inds.set(iData::AGE,std::numeric_limits<double>::signaling_NaN());
						// the above work for ind in prefile, but not for those not in prefile but automatically added by add_ind()
					} // Do nothing if can't convert to double, so missing(number) better be 0. Particularly, age=0 seems reasonable but don't allow it.
				}
				else 
				{
					inds.xtr[i]=inf[i];
				}
			}
			
			// prepare IDs
			if (!g_pid.empty()) inds.set(iData::PID,g_pid);
			try 
			{
				set_ind_ID(inds);
			}
			catch (bad_ind_ID &) 
			{
				lns<<writee<<"Line "<<inf.RowNumber()+1<<" has an ID problem; skipped."<<flush_logger;
				++skip_bad;
				continue; 
			}
			
			// find duplicated individual
			string& iid=inds.get(iData::IID);
			string& pid=inds.get(iData::PID);
			bool proband=false;
			perch::read_SeqID(iid,proband);
			if (proband) inds.set(iData::PRB,1);
			Individual * old_ind = find_ind(pid,iid);
			if (old_ind!=NULL)
			{
				if (old_ind->get(iData::ABP)!=true)
				{
					if (to_show_warnings && to_report_re_inp)
						lns<<showw<<"Re-input individual "<<pid<<" - "<<iid<<flush_logger;
					++dup_ind;
				}
			}
			
			// write rows_input_seq and find duplicated family
			if (pid!=last_pid)
			{
				last_pid=pid;
				vector<string> empty_iid_vec;
				rows_input_seq.push_back(std::pair< string,vector<string> >(last_pid,empty_iid_vec));
				if (exist_element(family_db,last_pid))
				{
					if (to_show_warnings && to_report_re_inp)
						lns<<showw<<"Re-input family "<<last_pid<<flush_logger;
					++dup_fam;
				}
			}
			rows_input_seq.back().second.push_back(inds.get(iData::IID));
			
			// add
			add_ind(inds);
			++added;
		}
	}
	
	lns<<showl<<"Number of lines (excl. header): "<<added<<flush_logger;
	lns<<showl<<"Skip lines lacking some fields: "<<skip_sht<<flush_logger;
	lns<<showl<<"Skip lines whose PID=empty/0/.: "<<skip_bad<<flush_logger;
	lns<<showl<<"Multiple input of the same IND: "<<dup_ind<<flush_logger;
	lns<<showl<<"Multiple input of the same PED: "<<dup_fam<<flush_logger;
	read_family_errors();
	check_questionable_IDs();
	check_self_ancestor();
	check_loop();
	check_separation();
	check_circles();
	if (to_show_warnings && to_report_subped) cal_clustering("",false,false);
//	cal_generation();
	if (iData::vec_contain(inpFields, iData::MZT)) check_mzt();
	if (iData::vec_contain(inpFields, iData::PRB)) check_prb();
	if (iData::vec_contain(inpFields, iData::YOB)) check_yob();
	if (iData::vec_contain(inpFields, iData::AFF)) check_aff();
	if (iData::vec_contain(inpFields, iData::LIA)) check_lia();
	lns.endsub();
}

void FamilyDB::_auto_add_parents(Individual* ind)
{
	string Father_ID=given_iid_pfx+itos(++given_iid);
	string Mother_ID=given_iid_pfx+itos(++given_iid);
	Individual dad,mom;
	dad.set(iData::PID,ind->fam->pid);
	dad.set(iData::IID,Father_ID);
	dad.set(iData::SEX,kMAN);
	dad.set(iData::AFF,0);
	dad.set(iData::AGE,0);
	dad.set(iData::UID,"");
	set_ind_ID(dad);
	Individual* fth = add_ind(dad,false);
	mom.set(iData::PID,ind->fam->pid);
	mom.set(iData::IID,Mother_ID);
	mom.set(iData::SEX,kWMN);
	mom.set(iData::AFF,0);
	mom.set(iData::AGE,0);
	mom.set(iData::UID,"");
	set_ind_ID(mom);
	Individual* mth = add_ind(mom,false);
	make_parents(ind,fth,mth);
}

void FamilyDB::_fast_add_set(vector<string>& arguments, Individual* ind)
{
	string aff;	if (arguments.size()>=4) aff=arguments[3];
	string age;	if (arguments.size()>=5) age=arguments[4];
	string gtp;	if (arguments.size()>=6) gtp=arguments[5];
	string alt;	if (arguments.size()>=7) alt=arguments[6];
	string cmt;	if (arguments.size()>=8) cmt=arguments[7];

	// read aff. no_input is 0.
	int aff_code = 0;
	boost::to_lower(aff);
	if		(aff=="a" || aff=="aff"   || aff=="affected"   || aff=="2")	aff_code=2;
	else if	(aff=="u" || aff=="unaff" || aff=="unaffected" || aff=="1")	aff_code=1;
	else if (                            aff=="unknown"    || aff=="0")	aff_code=0;
	else if (read_val_ge(aff,aff_code,0)) ;
	else exit_error("Affection status "+aff+" is not allowed.");
	
	// read age. no_input or missing is nan.
	double age_num = 0;
	if (age!="." && !age.empty() && !read_val_lt(age,age_num,150.0)) exit_error("age ("+age+") wrong for "+arguments[0]+":"+arguments[1]);
	if (std::isnan(age_num) || age_num<0) age_num = 0;

	// read gtp
	boost::to_lower(gtp);
	if		(gtp=="het"||gtp=="het."||gtp=="-/+"||gtp=="+/-")	gtp="Het";
	else if (gtp=="hom"||gtp=="hom."||gtp=="+/+")				gtp="Hom";
	else if (gtp=="neg"||gtp=="neg."||gtp=="-/-")				gtp="Neg";
	else if (gtp.empty()|| gtp=="." ||gtp=="./."||gtp=="0")		gtp.clear();
	else exit_error("Genotype "+gtp+" is not allowed.");
	
	// update
	if (aff_code)		ind->set(iData::AFF,aff_code);
	if (age_num)		ind->set(iData::AGE,age_num);
	if (!gtp.empty())	ind->set(iData::GTP,gtp);
	if (!alt.empty())	ind->set(iData::ALT,alt);
	if (!cmt.empty())	ind->set(iData::CMT,cmt);
}

void FamilyDB::fast_add(const string& parameters)
{
	vector<string> arguments;
	boost::split(arguments,parameters,boost::is_any_of(","));
	if (arguments.size()<3) exit_error("insufficient arguments for --fast-add");
	string fam_name = arguments[0];
	string ind_name = arguments[1];
	string relation = arguments[2];
	
	// read fam_name and ind_name
	Individual* ind = find_ind(fam_name, ind_name);
	if (ind==NULL) exit_error("Family "+fam_name+" member "+ind_name+" not found in the database.");
	
	// determin given_iid starting point
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->members,iti))
		{
			const string& iid = iti->first;
			int num=-1;
			if (str_startsw(iid,given_iid_pfx) && read_val(iid.substr(1),num) && num>given_iid) given_iid=num;
		}
	}
	
	// read rel. must input
	boost::to_lower(relation);
	if	(relation=="set") {
		_fast_add_set(arguments,ind);
	}
	else if	(str_startsw(relation,"dad")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		_fast_add_set(arguments,ind->pfa);
	}
	else if (str_startsw(relation,"mom")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		_fast_add_set(arguments,ind->pmo);
	}
	else if (str_startsw(relation,"bro")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,ind->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kMAN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"sis")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,ind->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kWMN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"son")) {
		string Spouse_ID;
		if		(ind->spouses.empty())		Spouse_ID=given_iid_pfx+itos(++given_iid);
		else if (ind->spouses.size()==1)	Spouse_ID=ind->spouses.begin()->first;
		else exit_error("Family "+fam_name+" member "+ind_name+" has >1 spouses. I don't know how to add a son.");
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		if (ind->get(iData::SEX)==kMAN)
		{
			dummy_ind.set(iData::FTH,ind->get(iData::IID));
			dummy_ind.set(iData::MTH,Spouse_ID);
		}
		else
		{
			dummy_ind.set(iData::FTH,Spouse_ID);
			dummy_ind.set(iData::MTH,ind->get(iData::IID));
		}
		dummy_ind.set(iData::SEX,kMAN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"dau")) {
		string Spouse_ID;
		if		(ind->spouses.empty())		Spouse_ID=given_iid_pfx+itos(++given_iid);
		else if (ind->spouses.size()==1)	Spouse_ID=ind->spouses.begin()->first;
		else exit_error("Family "+fam_name+" member "+ind_name+" has >1 spouses. I don't know how to add a son.");
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		if (ind->get(iData::SEX)==kMAN)
		{
			dummy_ind.set(iData::FTH,ind->get(iData::IID));
			dummy_ind.set(iData::MTH,Spouse_ID);
		}
		else
		{
			dummy_ind.set(iData::FTH,Spouse_ID);
			dummy_ind.set(iData::MTH,ind->get(iData::IID));
		}
		dummy_ind.set(iData::SEX,kWMN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"phb")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		string Spouse_ID;
		if		(ind->pfa->spouses.size()==1)	Spouse_ID=given_iid_pfx+itos(++given_iid);
		else if (ind->pfa->spouses.size()==2) {	int i=0; string spo[2]; for (auto &it:ind->pfa->spouses) spo[i++]=it.first; if (ind->pmo->get(iData::IID)==spo[0]) Spouse_ID=spo[1]; else Spouse_ID=spo[0]; }
		else exit_error("Family "+fam_name+" member "+ind_name+"'s father has >2 spouses. I don't know how to add a paternal half-bro.");
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,Spouse_ID);
		dummy_ind.set(iData::SEX,kMAN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"phs")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		string Spouse_ID;
		if		(ind->pfa->spouses.size()==1)	Spouse_ID=given_iid_pfx+itos(++given_iid);
		else if (ind->pfa->spouses.size()==2) {	int i=0; string spo[2]; for (auto &it:ind->pfa->spouses) spo[i++]=it.first; if (ind->pmo->get(iData::IID)==spo[0]) Spouse_ID=spo[1]; else Spouse_ID=spo[0]; }
		else exit_error("Family "+fam_name+" member "+ind_name+"'s father has >2 spouses. I don't know how to add a paternal half-sis.");
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,Spouse_ID);
		dummy_ind.set(iData::SEX,kWMN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"mhb")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		string Spouse_ID;
		if		(ind->pmo->spouses.size()==1)	Spouse_ID=given_iid_pfx+itos(++given_iid);
		else if (ind->pmo->spouses.size()==2) {	int i=0; string spo[2]; for (auto &it:ind->pmo->spouses) spo[i++]=it.first; if (ind->pfa->get(iData::IID)==spo[0]) Spouse_ID=spo[1]; else Spouse_ID=spo[0]; }
		else exit_error("Family "+fam_name+" member "+ind_name+"'s mother has >2 spouses. I don't know how to add a paternal half-bro.");
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,Spouse_ID);
		dummy_ind.set(iData::MTH,ind->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kMAN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"mhs")) {
		if (ind->pfa==NULL) _auto_add_parents(ind);
		string Spouse_ID;
		if		(ind->pmo->spouses.size()==1) {	Spouse_ID=given_iid_pfx+itos(++given_iid); }
		else if (ind->pmo->spouses.size()==2) {	int i=0; string spo[2]; for (auto &it:ind->pmo->spouses) spo[i++]=it.first; if (ind->pfa->get(iData::IID)==spo[0]) Spouse_ID=spo[1]; else Spouse_ID=spo[0]; }
		else exit_error("Family "+fam_name+" member "+ind_name+"'s mother has >2 spouses. I don't know how to add a paternal half-sis.");
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,Spouse_ID);
		dummy_ind.set(iData::MTH,ind->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kWMN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"pun")) {
		if (ind->pfa==NULL)			_auto_add_parents(ind);
		if (ind->pfa->pfa==NULL)	_auto_add_parents(ind->pfa);
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pfa->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,ind->pfa->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kMAN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"pau")) {
		if (ind->pfa==NULL)			_auto_add_parents(ind);
		if (ind->pfa->pfa==NULL)	_auto_add_parents(ind->pfa);
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pfa->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,ind->pfa->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kWMN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"mun")) {
		if (ind->pmo==NULL)			_auto_add_parents(ind);
		if (ind->pmo->pfa==NULL)	_auto_add_parents(ind->pmo);
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pmo->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,ind->pmo->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kMAN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"mau")) {
		if (ind->pmo==NULL)			_auto_add_parents(ind);
		if (ind->pmo->pfa==NULL)	_auto_add_parents(ind->pmo);
		Individual dummy_ind;
		dummy_ind.set(iData::PID,ind->fam->pid);
		dummy_ind.set(iData::IID,given_iid_pfx+itos(++given_iid));
		dummy_ind.set(iData::FTH,ind->pmo->pfa->get(iData::IID));
		dummy_ind.set(iData::MTH,ind->pmo->pmo->get(iData::IID));
		dummy_ind.set(iData::SEX,kWMN);
		_fast_add_set(arguments,&dummy_ind);
		add_ind(dummy_ind,false);
	}
	else if (str_startsw(relation,"pgf")) {
		if (ind->pfa==NULL)			_auto_add_parents(ind);
		if (ind->pfa->pfa==NULL)	_auto_add_parents(ind->pfa);
		_fast_add_set(arguments,ind->pfa->pfa);
	}
	else if (str_startsw(relation,"pgm")) {
		if (ind->pfa==NULL)			_auto_add_parents(ind);
		if (ind->pfa->pfa==NULL)	_auto_add_parents(ind->pfa);
		_fast_add_set(arguments,ind->pfa->pmo);
	}
	else if (str_startsw(relation,"mgf")) {
		if (ind->pmo==NULL)			_auto_add_parents(ind);
		if (ind->pmo->pfa==NULL)	_auto_add_parents(ind->pmo);
		_fast_add_set(arguments,ind->pmo->pfa);
	}
	else if (str_startsw(relation,"mgm")) {
		if (ind->pmo==NULL)			_auto_add_parents(ind);
		if (ind->pmo->pfa==NULL)	_auto_add_parents(ind->pmo);
		_fast_add_set(arguments,ind->pmo->pmo);
	}
	else exit_error("relationship "+relation+" not allowed");
}

// requires ABP. Must call after read_family_errors();
void FamilyDB::check_questionable_IDs()
{
	lns<<showl<<"Check for questionable IIDs ... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf))
	{
		int parent_read=0;
		for (each_element(itf->second->members,iti))
		{
			Individual* ind=iti->second;
			if (ind->pfa!=NULL)
			{
				if (ind->pfa->get(iData::ABP)!=true) ++parent_read;
				if (ind->pmo->get(iData::ABP)!=true) ++parent_read;
			}
		}
		if (parent_read==0)
			; // if (to_show_warnings) lns<<showw<<"Family "<<itf->second->pid<<" IIDs never show up as Father_IDs / Mother_IDs.";
	}
	lns.endsub();
}

void FamilyDB::downcode_IDs()
{
	lns<<showl<<"Downcode PIDs/IIDs/UIDs to integers starting from 1"<<flush_logger;
	int newPID=1;
	int newUID=1;
	for (each_element(family_db,itf), ++newPID) {
		int newIID=1;
		for (each_element(itf->second->members,iti), ++newIID, ++newUID) {
			iti->second->set(iData::DNP,itos(newPID));
			iti->second->set(iData::DNI,itos(newIID));
			iti->second->set(iData::DNU,itos(newUID));
		}
	}
}

double FamilyDB::num_sibs(Individual* ind)
{
	if (ind->pfa==NULL) return 0;
	IndMap paternal;
	IndMap maternal;
	IndMap fullsibs;
	for (each_element_const(ind->pfa->nextgen,it)) paternal[it->second->get(iData::IID)]=it->second; paternal.erase(ind->get(iData::IID));
	for (each_element_const(ind->pmo->nextgen,it)) maternal[it->second->get(iData::IID)]=it->second; maternal.erase(ind->get(iData::IID));
num_sibs_find_fullsibs:
	for (each_element(paternal,it))
		if (exist_element(maternal,it->first))
		{
			fullsibs[it->first]=it->second;
			maternal.erase(it->first);
			paternal.erase(it);
			goto num_sibs_find_fullsibs;
		}
	return fullsibs.size() + (paternal.size()+maternal.size())/2.0;
}

void FamilyDB::find_1st_degr_relatives(Individual* ind, IndMap& dst)
{
	dst.clear();
	if (ind->pfa)
	{
		IndMap paternal;
		IndMap maternal;
		for (each_element_const(ind->pfa->nextgen,it)) paternal[it->second->get(iData::IID)]=it->second; paternal.erase(ind->get(iData::IID));
		for (each_element_const(ind->pmo->nextgen,it)) maternal[it->second->get(iData::IID)]=it->second; maternal.erase(ind->get(iData::IID));
		for (each_element(paternal,it))
			if (exist_element(maternal,it->first))
				dst[it->first]=it->second;
		dst[ind->pfa->get(iData::IID)]=ind->pfa;
		dst[ind->pmo->get(iData::IID)]=ind->pmo;
	}
	for (each_element_const(ind->nextgen,it)) dst[it->first]=it->second;
}

inline bool FamilyDB::is_affected(Individual* ind)
{
	return ind->get(iData::AFF) == kAFF;
}

inline bool FamilyDB::is_proband(Individual* ind)
{
	return ind->get(iData::PRB) != 0;
}

// circles are ties between 2 families (eg 2 brothers married to 2 sisters, or one man married 2 sisters, etc) in which
// the relationships (through blood / marriage) link several persons into a circle. Unlike loop that is linked only by 
// blood, circles are legal and common in the real world, but could lead to problems in genetic testing such as the 
// MelaPRO program, because the cousins have stronger genetic homogeneity than is expected. 
// This function only check for circles for up to 3 the grand-parents' generation. Between 2 nuclear families:
//
//   1-+-2       3-+-4
//     |           |
//   +-+-+     +---+---+
//   |   |     |   |   |
//  f5  m6--+--7  m8  f9
//          |
//          10
//
// there're 4 kinds of circles if these couples married, respectively: 8-5 / 3-2 / 8-2 / 6-9.
// By default, the function does not break up the circles unless it was asked for.
// return num_fam_with_circles
int FamilyDB::_check_circles(bool to_breakup)
{
	if (to_breakup) lns<<showl<<"Check for & break up circles... "<<flush_logger;
	else			lns<<showl<<"Check for circles           ... "<<flush_logger;
	lns.sub();
	typedef map< string,set<string> >	marriage_t;
	typedef set<marriage_t>				circle_t;
	int num_fam_with_circles=0;
	for (each_element_const(family_db,itf))
	{
		circle_t circles;
		for (each_element_const(itf->second->nonfounders,iti))
		{
			Individual* faptr = iti->second->pfa;
			Individual* moptr = iti->second->pmo;
			IndMap paternal; paternal[faptr->get(iData::IID)]=faptr;
			IndMap maternal; maternal[moptr->get(iData::IID)]=moptr;
			if (faptr->pfa)
			{
				paternal[faptr->pfa->get(iData::IID)]=faptr->pfa;
				paternal[faptr->pmo->get(iData::IID)]=faptr->pmo;
				for (each_element_const(faptr->pfa->nextgen,it2)) paternal[it2->second->get(iData::IID)]=it2->second;
				for (each_element_const(faptr->pmo->nextgen,it2)) paternal[it2->second->get(iData::IID)]=it2->second;
			}
			if (moptr->pfa)
			{
				maternal[moptr->pfa->get(iData::IID)]=moptr->pfa;
				maternal[moptr->pmo->get(iData::IID)]=moptr->pmo;
				for (each_element_const(moptr->pfa->nextgen,it2)) maternal[it2->second->get(iData::IID)]=it2->second;
				for (each_element_const(moptr->pmo->nextgen,it2)) maternal[it2->second->get(iData::IID)]=it2->second;
			}
			marriage_t marriages;
			for (each_element_const(paternal,it2))
				for (each_element_const(it2->second->spouses,it3))
					if (exist_element(maternal,it3->second->get(iData::IID)))
					{
						if (it2->second->get(iData::SEX)==kMAN)
							marriages[it2->second->get(iData::IID)] .insert( it3->second->get(iData::IID) );
						else
							marriages[it3->second->get(iData::IID)] .insert( it2->second->get(iData::IID) );
					}
			int marr_size=0;
			for (each_element_const(marriages,it)) marr_size+=it->second.size();
			if (marr_size>1) circles.insert(marriages);
		}
		if (circles.empty())
		{
			itf->second->circles=0;
		}
		else
		{
			itf->second->circles=1;
			++num_fam_with_circles;
			for (each_element_const(circles, it1))
			{
				lns<<writew<<"Family "<<itf->second->pid<<" circle ("; // writel => writew
				for (each_element_const(*it1, it2))
					for (each_element_const(it2->second, it3))
						lns<<" "<<it2->first<<"<=>"<<itf->second->members[*it3]->get(iData::IID);
				lns<<" )"<<flush_logger;
			}
			if (to_breakup)
			{
				double		min_sibs=DBL_MAX;	// minimum number of sibs
				set<string> min_name;			// individuals have the min_sibs
				for (each_element_const(circles, it1))
				{
					for (each_element_const(*it1, it2))
					{
						if (itf->second->members[it2->first]->pfa != NULL) // not a founder
						{
							double s = num_sibs(itf->second->members[it2->first]);
							if		(s <min_sibs) { min_sibs = s; min_name.clear();	min_name.insert(it2->first); }
							else if (s==min_sibs) {									min_name.insert(it2->first); }
						}
						for (each_element_const(it2->second,it3))
						{
							if (itf->second->members[*it3]->pfa == NULL) continue; // founder
							double s = num_sibs(itf->second->members[*it3]);
							if		(s <min_sibs) { min_sibs = s; min_name.clear(); min_name.insert(*it3); }
							else if (s==min_sibs) {									min_name.insert(*it3); }
						}
					}
				}
				// it's impossible that all candidates are founders, so no need to check for min_name.empty()
				string reason = "min(#sibs)";
				Individual* ind=NULL;
				if (min_name.size()==1) 
					ind=itf->second->members[*min_name.begin()];
				else {
					set<string> name_unaff;
					for (each_element_const(min_name, it)) 
						if (!is_affected(itf->second->members[*it])) name_unaff.insert(*it);
					if (name_unaff.empty()) name_unaff=min_name;
					else reason += " & unaff";
					if (name_unaff.size()==1) 
					{
						ind=itf->second->members[*name_unaff.begin()];
					}
					else {
						int min_affrel=INT_MAX;
						set<string> name_affrel;
						for (each_element_const(name_unaff, it))
						{
							IndMap fst_dgr;
							find_1st_degr_relatives(itf->second->members[*it],fst_dgr);
							int num_affrel=0;
							for (each_element_const(fst_dgr,it2)) if (is_affected(it2->second)) ++num_affrel;
							if (num_affrel < min_affrel) { min_affrel=num_affrel; name_affrel.clear(); name_affrel.insert(*it); }
							else if (num_affrel == min_affrel) name_affrel.insert(*it);
						}
						reason += " & min(Af_1stDgrRel)";
						if (name_affrel.size()!=1) reason += " & 1st_of("+itos(name_affrel.size())+")";
						ind=itf->second->members[*name_affrel.begin()];
					}
				}
				founderize(ind);
				lns<<" founderize "<<*min_name.begin() << " by " << reason<<flush_logger;
			}
		}
	}
	if (num_fam_with_circles) lns<<cshowl<<"Found circles in "<<num_fam_with_circles<<" families."<<flush_logger;
	lns.endsub();
	if (to_remove_sep_id) { check_separation(); del_inds_by_status(iData::SEP,true); }
	return num_fam_with_circles;
}

void FamilyDB::check_circles(bool to_breakup)
{
	if (to_breakup) while (_check_circles(to_breakup)) ;
	else _check_circles(to_breakup);
}

// add a person for clustering algorithm, called by cal_clustering() only
inline void _cl_add(Individual* ind, set<string>& all_can, queue<string>& seeds, set<Individual*>& old_rel, set<Individual*>& new_rel, int& clst_ID)
{
	if (ind==NULL) return;
	if (!exist_element(old_rel, ind)) new_rel.insert(ind);
	if (ind->get(iData::CLT)) return;
	string uid=ind->get(iData::UID);
	if (exist_element(all_can, uid)) { seeds.push(uid); all_can.erase(uid); }
	ind->set(iData::CLT,clst_ID);
}

void FamilyDB::cal_clustering(const string& cluster_uid_file, bool cluster_affected, bool to_write)
{
	int old_clt_degr_=clt_degr_;
	
	// candidate UIDs for clustering. Only these people and their relatives within N degrees are clusters.
	set<string> all_can;
	if (!cluster_uid_file.empty())
	{
		for (Rows_in_File(in, cluster_uid_file, 1))
			if (!exist_element(uniqid_db, in[0])) exit_error("UID "+in[0]+" not found in the database.");
			else all_can.insert(in[0]);
		lns<<showl<<"Do clustering of "<<all_can.size()<<" individuals in "<<cluster_uid_file<<" among "<<clt_degr_<<"-degree relatives ... "<<flush_logger;
		lns.sub();
	}
	else if (cluster_affected)
	{
		for (each_element(uniqid_db, it)) if (is_affected(it->second))			all_can.insert(it->first);
		lns<<showl<<"Do clustering of "<<all_can.size()<<" affected individuals among "<<clt_degr_<<"-degree relatives ... "<<flush_logger;
		lns.sub();
	}
	else // cluster all individuals
	{
		clt_degr_=1; // clt_degr_ doesn't matter anymore but clt_degr_=1 is faster
		for (each_element(uniqid_db, it)) all_can.insert(it->first);
		lns<<showl<<"Check for >1 sub-peds / PID ... "<<flush_logger;
		lns.sub();
	}
	
	// prepare
	set_inds_values(iData::CLT, 0);
	set_inds_status(iData::VST, false);
	if (to_write) _print_header();

	// clustering
	map<string, set<int> > cl_perFam;
	int clst_ID;
	for (clst_ID=1; !all_can.empty(); ++clst_ID)
	{
		// when UID add to seed, must: rm from all_can; set CLT, add to all_rel. Only candidates can be seeds.
		queue<string> seeds;					// candidates as seeds
		string seed_one = *all_can.begin();		// 1st candidate as 1st seed
		seeds.push(seed_one);
		all_can.erase(seed_one);
		Individual* ind = uniqid_db[seed_one];
		ind->set(iData::CLT,clst_ID);
		set<Individual*> all_rel;
		while (!seeds.empty())					// when there're seeds, work on each
		{
			string seed_ths = seeds.front();	// this seed
			seeds.pop();
			Individual* ind = uniqid_db[seed_ths];
			set<Individual*> old_rel;			// relatives of this seed, must start from this seed ONLY
			old_rel.insert(ind);
			all_rel.insert(ind);
			for (int d=0;d<clt_degr_;++d)
			{
				set<Individual*> new_rel;
				for (each_element(old_rel, it))
				{
					if ((*it)->get(iData::VST)==true) continue; else (*it)->set(iData::VST,true);
					_cl_add((*it)->pfa, all_can, seeds, old_rel, new_rel, clst_ID); // add parents
					_cl_add((*it)->pmo, all_can, seeds, old_rel, new_rel, clst_ID); // add parents
					if ((*it)->pfa) for (each_element((*it)->pfa->spouses,stp)) _cl_add(stp->second, all_can, seeds, old_rel, new_rel, clst_ID); // add step-fa/mo
					if ((*it)->pmo) for (each_element((*it)->pmo->spouses,stp)) _cl_add(stp->second, all_can, seeds, old_rel, new_rel, clst_ID); // add step-fa/mo
					if ((*it)->pfa) for (each_element((*it)->pfa->nextgen,sib)) _cl_add(sib->second, all_can, seeds, old_rel, new_rel, clst_ID); // add sibs
					if ((*it)->pmo) for (each_element((*it)->pmo->nextgen,sib)) _cl_add(sib->second, all_can, seeds, old_rel, new_rel, clst_ID); // add sibs
					for (each_element((*it)->nextgen, ch)) _cl_add(ch->second, all_can, seeds, old_rel, new_rel, clst_ID); // add childrens
					for (each_element((*it)->spouses, sp)) _cl_add(sp->second, all_can, seeds, old_rel, new_rel, clst_ID); // add spouses
				}
				if (new_rel.empty()) break;
				old_rel.insert(new_rel.begin(), new_rel.end());
				all_rel.insert(new_rel.begin(), new_rel.end());
			}
			if (clt_degr_!=1) set_inds_status(iData::VST, false); // too slow when I want to cluster everyone
		}
		
		// one cluster is made, now print it out
		if (to_write)
			for (each_element(all_rel, it)) _print_line(*it,NULL);

		// summarize
		for (each_element(all_rel, it)) cl_perFam[(*it)->fam->pid].insert((*it)->get(iData::CLT));
		
		// finishing
		set_inds_values(iData::CLT, 0);
		set_inds_status(iData::VST, false);
	}

	// update f->num_clt
	for (each_element(family_db,itf)) itf->second->num_clt=-1;
	for (each_element(cl_perFam,itc)) family_db[itc->first]->num_clt = itc->second.size();
	
	// if cluster all individuals (good for QC), output summary
	if (cluster_uid_file.empty() && !cluster_affected) {
		int num_fam_with_multi_cl=0;
		for (each_element(cl_perFam, it))
			if (it->second.size()>1) {
				lns<<writew<<"Family "<<it->first<<" has "<<it->second.size()<<" clusters."<<flush_logger; // showw=>writew, to reduce screen printing
				++num_fam_with_multi_cl;
			}
		if (num_fam_with_multi_cl) lns<<cshowl<<"Found multiple clusters in "<<num_fam_with_multi_cl<<" families."<<flush_logger;
	}
	lns.endsub();
	
	clt_degr_=old_clt_degr_;
}

/*
 iid fath moth problem corrected
 a    f    m   correct   a f m
 a    a    a   self      a . .
 a    a    m   self      a . m
 a    f    a   self      a f .
 a    .    a   self      a . .
 a    a    .   self      a . .
 a    f    f   same      a f .
 a    m    m   same      a . m
 a    o    o   same      a o .
 a    f    .   single    a f .
 a    .    m   single    a . m
*/
void FamilyDB::read_family_errors()
{
	if (add_ind_logf_.str().empty()) return;
	Individual * ind;
	Individual dummy_ind;
	string *kwnp,*uknp;
	enum problem_type {no_spouse=1,has_spouse,spouses};
	enum problem_solv {nochoice,drop,create,useold};
	problem_solv init = always_do_create ? create : nochoice;
	int  choices[4] = { init, init, init, init };
	int  tot_dummy=0;
	
	lns<<showl<<"Correct single_parent errors ... "<<flush_logger;
	lns.sub();
	for (each_line(add_ind_logf_))
	{
		string fid,iid,fth,mth,pbl;
		double sex,aff;
		add_ind_logf_ >> fid >> iid >> fth >> mth >> sex >> aff >> pbl;

		lns<<writel<<fid<<"-"<<iid<<" :";
		int type=0;
		if ((fth!=string_no_id_ && mth==string_no_id_) || (fth==string_no_id_ && mth!=string_no_id_))
		{
			lns<<" single-parent.";
			if (fth!=string_no_id_)	{	kwnp=&fth;	uknp=&mth;	}
			else					{	kwnp=&mth;	uknp=&fth;	}
			ind=find_ind(fid,*kwnp);
			if (ind!=NULL)
			{
				if ((ind->get(iData::SEX)==kWMN && fth!=string_no_id_)||(ind->get(iData::SEX)==kMAN && mth!=string_no_id_))
				{
					swap(fth, mth);		// tmp=fth;  fth=mth; mth=tmp;
					swap(kwnp, uknp);	// tmpp=kwnp; kwnp=uknp; uknp=tmpp;
					lns<<" Parent gender wrong, corrected.";
				}
				int n_spouses = 0;
				for (each_element(ind->spouses,its))
					if (its->second->get(iData::DUM)==true) *uknp=its->second->get(iData::IID);
					else									++n_spouses;
				
				if (n_spouses>1)  { type=spouses;		lns<<" Parent has >1 spouses."; }
				if (n_spouses==1) { type=has_spouse;	lns<<" Parent has  1 spouse.";  }
				if (n_spouses==0) { type=no_spouse;		lns<<" Parent has no spouse.";  }
			}
			else
			{
				type=no_spouse;
				lns<<" Parent not found.";
			}
			//now type is set, sex correct
			char c='U';
			if (choices[type]==nochoice)
			{
				lns<<showl<<">> Problem: single_parent, ";
				switch (type)
				{
					case no_spouse:
						lns<<"parent has no spouse / parent not found. ";
						lns<<showl<<"   Options: (c)reate dummy / (f)ounderize person."; 
						break;
					case spouses:
						lns<<"parent has >1 spouses.";
						lns<<showl<<"   Options: (c)reate dummy / use first (s)pouse / (f)ounderize person."; 
						break;
					case has_spouse:
						lns<<"parent has  1 spouse.";
						lns<<showl<<"   Options: (c)reate dummy / use known (s)pouse / (f)ounderize person."; 
						break;
					default:
						lns<<"unknown.";
						break;
				};
				lns<<showl<<"   Please choose and type the corresponding letter inside the parentheses";
				lns<<showl<<"   (lower case = this time only; upper case = remember your choice)";
inputchoice:
				c=lns.logger_input_char();
				if (type==no_spouse && (c=='s'||c=='S'))
				{
					lns<<showl<<"   Wrong letter inputted, please type again";
					goto inputchoice;
				}
				switch (c)
				{
					case 'c':
					case 'C': choices[type]=create; break;
					case 'f':
					case 'F': choices[type]=drop; break;
					case 's':
					case 'S': choices[type]=useold; break;
					default : lns<<showl<<"   Wrong letter inputted, please type again"; goto inputchoice;
				};
				lns<<showl<<"   You chose to ";
				switch (c)
				{
					case 'c':lns<<"create a dummy parent, for this time only.";				break;
					case 'C':lns<<"create a dummy parent, for similar problems.";			break;
					case 'f':lns<<"founderize this person, for this time only.";			break;
					case 'F':lns<<"founderize this person, for similar problems.";			break;
					case 's':lns<<"use known spouse of the parent, for this time only.";	break;
					case 'S':lns<<"use known spouse of the parent, for similar problems.";	break;
				};
				
			}
			switch (choices[type])
			{
				case nochoice:	lns<<" No choice!";							break;
				case drop:		lns<<" Founderize person "<<fid<<"-"<<iid;	break;
				case useold:	
					if (mth==string_no_id_) mth=ind->spouses.begin()->second->get(iData::IID);
					if (fth==string_no_id_) fth=ind->spouses.begin()->second->get(iData::IID);
					dummy_ind.set(iData::IID,iid);
					dummy_ind.set(iData::SEX,sex);
					dummy_ind.set(iData::PID,fid);
					dummy_ind.set(iData::FTH,fth);
					dummy_ind.set(iData::MTH,mth);
					set_ind_ID(dummy_ind);
					add_ind(dummy_ind,false);
					lns<<" Use spouse "<<ind->spouses.begin()->second->get(iData::IID);
					break;
				case create:
					if (*uknp == string_no_id_)
					{
						for (int i=1;i<=100;++i)
						{
							if (i==100) { lns<<showe<<"Create dummy failed after 100x trials."<<fatal; }
							*uknp = *kwnp + "_" + itos_format(i, 'r', 2, '0', 10);
							if (find_ind(fid,*uknp)) continue;
							else break;
						}
						++tot_dummy;
						lns<<writel<<"Created dummy "<<*uknp;
					}
					dummy_ind.set(iData::IID,iid);
					dummy_ind.set(iData::SEX,sex);
					dummy_ind.set(iData::PID,fid);
					dummy_ind.set(iData::FTH,fth);
					dummy_ind.set(iData::MTH,mth);
					set_ind_ID(dummy_ind);
					add_ind(dummy_ind,false);
					find_ind(fid, *uknp)->set(iData::DUM,true);
					break;
				default:
					lns<<" Wrong choice!";
					break;
			};
			if (c>='a' && c<='z') choices[type]=nochoice;
			continue;
		}
		lns<<" unknown type of error."<<flush_logger;
	}
	add_ind_logf_.str(std::string());
	add_ind_logf_.clear();
	lns<<showl<<tot_dummy<<" dummy individuals created."<<flush_logger;
	lns.endsub();
}

Individual* FamilyDB::find_ind(string  fam_name, string  ind_name)
{
	if (exist_element(family_db,fam_name) && exist_element(family_db[fam_name]->members,ind_name))
		return family_db[fam_name]->members[ind_name];
	else
		return NULL;
}

void FamilyDB::check_loop(bool to_breakup)
{
//	if (to_breakup)
//		while (_check_loop(to_breakup)) ;
//	else
		_check_loop(to_breakup);
}

// must be aft check_self_ancestor(), otherwise dead loop
int FamilyDB::_check_loop(bool to_breakup)
{
	if (to_breakup) lns<<showl<<"Check for & break up loops  ... "<<flush_logger;
	else			lns<<showl<<"Check for loops             ... "<<flush_logger;
	lns.sub();

	set<string> fams_with_loop;
	for (each_element(family_db,itf))
	{
	_check_loop_for_fam:
		itf->second->consang.clear();
		for (each_element(itf->second->members,iti))
		{
			FamilyDB::ID_Count ancestors;
			_find_ancestors(iti->second,ancestors); // should include self
			
			for (each_element(iti->second->spouses,its))
			{
				FamilyDB::ID_Count spouse_ancestors;
				_find_ancestors(its->second,spouse_ancestors); // should include spouse
				
				FamilyDB::ID_Count converge_anc;
				for (each_element(ancestors,ita))
					if (exist_element(spouse_ancestors,ita->first)) ++converge_anc[ita->first]; // prv: if (ita->second>1)
				
				if (!converge_anc.empty())
				{
					fams_with_loop.insert(itf->first);
					string marr_str = iti->second->get(iData::IID)+" "+its->second->get(iData::IID);
					if (exist_element(itf->second->consang,marr_str)) continue;
					itf->second->consang.insert(iti->second->get(iData::IID)+" "+its->second->get(iData::IID));
					itf->second->consang.insert(its->second->get(iData::IID)+" "+iti->second->get(iData::IID));
					lns<<writew<<"consanguineous marriage in family "<<itf->second->pid<<" : "<<marr_str<<" ancestors: "<<str_of_map_first(converge_anc,' ')<<flush_logger; // writel => writew
					
					// to break the loops, similar to br-circle except for the creation of min_name
					if (to_breakup)
					{
						double		min_sibs=DBL_MAX;	// minimum number of sibs
						set<string> min_name;			// min_name to be founderized to break the loop
						for (each_element(converge_anc,ita))
							for (each_element(itf->second->members[ita->first]->nextgen,itn))
							{
								if (exist_element(ancestors,itn->first) || exist_element(spouse_ancestors,itn->first))
								{
									double s = num_sibs(itn->second);
									if		(s <min_sibs) { min_sibs = s; min_name.clear();	min_name.insert(itn->first); }
									else if (s==min_sibs) {									min_name.insert(itn->first); }
								}
							}
						string reason = "min(#sibs)";
						Individual* ind=NULL;
						if (min_name.size()==1)
							ind=itf->second->members[*min_name.begin()];
						else {
							set<string> name_unaff;
							for (each_element_const(min_name, it))
								if (!is_affected(itf->second->members[*it])) name_unaff.insert(*it);
							if (name_unaff.empty()) name_unaff=min_name;
							else reason += " & unaff";
							if (name_unaff.size()==1)
								ind=itf->second->members[*name_unaff.begin()];
							else {
								int min_affrel=INT_MAX;
								set<string> name_affrel;
								for (each_element_const(name_unaff, it))
								{
									IndMap fst_dgr;
									find_1st_degr_relatives(itf->second->members[*it],fst_dgr);
									int num_affrel=0;
									for (each_element_const(fst_dgr,it2)) if (is_affected(it2->second)) ++num_affrel;
									if (num_affrel < min_affrel) { min_affrel=num_affrel; name_affrel.clear(); name_affrel.insert(*it); }
									else if (num_affrel == min_affrel) name_affrel.insert(*it);
								}
								reason += " & min(Af_1stDgrRel)";
								if (name_affrel.size()!=1) reason += " & 1st_of("+itos(name_affrel.size())+")";
								ind=itf->second->members[*name_affrel.begin()];
							}
						}
						founderize(ind);
						lns<<" founderize "<<*min_name.begin() << " by " << reason<<flush_logger;
						goto _check_loop_for_fam;
					}
				}
			}
		}
	}
	if (!fams_with_loop.empty()) lns<<cshowl<<"Found loops in "<<fams_with_loop.size()<<" families."<<flush_logger;
	lns.endsub();
	if (to_remove_sep_id) { check_separation(); del_inds_by_status(iData::SEP,true); }
	return fams_with_loop.size();
}

void FamilyDB::_find_ancestors(Individual* ind, FamilyDB::ID_Count& ancestors)
{
	++ancestors[ind->get(iData::IID)];
	if (ind->pfa) _find_ancestors(ind->pfa,ancestors);
	if (ind->pmo) _find_ancestors(ind->pmo,ancestors);
}

void FamilyDB::_find_founders(Individual* ind, FamilyDB::ID_Count& ancestors)
{
	if (ind->pfa)
	{
		_find_founders(ind->pfa,ancestors);
		_find_founders(ind->pmo,ancestors);
	}
	else
		++ancestors[ind->get(iData::IID)];
}

void FamilyDB::_find_ancestors(Individual* ind, std::map<Individual*,int>& ancestors)
{
	++ancestors[ind];
	if (ind->pfa) _find_ancestors(ind->pfa,ancestors);
	if (ind->pmo) _find_ancestors(ind->pmo,ancestors);
}

void FamilyDB::_find_founders(Individual* ind, std::map<Individual*,int>& ancestors)
{
	if (ind->pfa)
	{
		_find_founders(ind->pfa,ancestors);
		_find_founders(ind->pmo,ancestors);
	}
	else
		++ancestors[ind];
}

void FamilyDB::check_self_ancestor()
{
	lns<<showl<<"Check for self-ancestors    ... "<<flush_logger;
	lns.sub();
	int total_errors=0;
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->members,iti))
		{
			ID_Count ancestors;
			_check_self_anc(iti->second,ancestors);
			if (ancestors[iti->second->get(iData::IID)] > 1)
			{
				lns<<showe<<"Family "<<iti->second->fam->pid<<" person "<<iti->second->get(iData::IID)<<" is self-ancestor."<<flush_logger;
				++total_errors;
			}
		}
	}
	if (total_errors) { lns<<showl<<"Errors occured, program quitted abnormally."<<fatal; }
	lns.endsub();
}

void FamilyDB::_check_self_anc(Individual* ind, FamilyDB::ID_Count& ancestors)
{
	if (exist_element(ancestors,ind->get(iData::IID))) { ancestors[ind->get(iData::IID)]++; return; }
	++ancestors[ind->get(iData::IID)];
	if (ind->pfa) _check_self_anc(ind->pfa,ancestors);
	if (ind->pmo) _check_self_anc(ind->pmo,ancestors);
}

bool has_arg(vector<string>& tok, vector<string>::const_iterator& it, int n)
{
	int i = it - tok.begin();
	int s = tok.size();
	if ( i+n >= s ) return false;
	else return true;
}

void need_arg(vector<string>& tok, vector<string>::const_iterator& it, int n)
{
	if (!has_arg(tok, it, n)) exit_error("Insufficient arguments for "+*it);
}

void exact_arg(vector<string>& tok, vector<string>::const_iterator& it, int n)
{
	int i = it - tok.begin();
	int s = tok.size();
	if ( i+n > s-1 ) exit_error("Insufficient arguments for "+*it);
	if ( i+n < s-1 ) exit_error("No option is allowed after "+*it);
}

const string& next_arg(vector<string>& tok, vector<string>::const_iterator& it)
{
	need_arg(tok,it,1);
	++it;
	return *it;
}

//return true if executed suscessfully
bool FamilyDB::execute_command(vector<string>& tok)
{
	for (each_element_const(tok,it))
	{
		if		(*it=="--ped"   ||*it=="--read")	{ need_arg(tok,it,1); while (has_arg(tok,it,1)) { string fe=next_arg(tok,it); string fd;				  read_prefile(fe,true ,fd); } }
		else if	(*it=="--ped-nh"||*it=="--read-nh") { need_arg(tok,it,2); while (has_arg(tok,it,2)) { string fe=next_arg(tok,it); string fd=next_arg(tok,it); read_prefile(fe,false,fd); } }
		else if (*it=="--fast-add") {	exact_arg(tok,it,1); fast_add(next_arg(tok,it)); }
		else if (*it=="--id-prefix") {	exact_arg(tok,it,1); given_iid_pfx=next_arg(tok,it); }
		else if (*it=="--br-circle")	check_circles(true);
		else if (*it=="--br-loop")		check_loop(true);
		else if (*it=="--var-known")	iData::set_known_var(next_arg(tok,it));
		else if (*it=="--keep-unaff")	unaff_sib_is_inf=true;
		else if (*it=="--create-dummy")	always_do_create=true;
		else if	(*it=="--cvt")		{	iData::read_conversion_map(next_arg(tok,it)); }
		else if (*it=="--genedrop")	{	gene_drop(false); }
		else if (*it=="--ind-wt")	{	ind_weight(); }
		else if (*it=="--iWt-sum1")	{	iWt_SumTo1_eaPed=true; }
		else if (*it=="--iWt-all")	{	iWt_useCasesOnly=false; }
		else if (*it=="--iWt-cs")	{	iWt_useCasesOnly=true; }
		else if (*it=="--sap")		{	single_aff_par(); }
		else if (*it=="--rm-uninf")	{	cal_informativeness(); del_inds_by_status(iData::INF,false);
										check_separation(); del_inds_by_status(iData::SEP,true); to_remove_sep_id=true; }
		else if (*it=="--rm-sep")	{	check_separation(); del_inds_by_status(iData::SEP,true); to_remove_sep_id=true; }
		else if	(*it=="--cl-dgr")	{	clt_degr_ = boost::lexical_cast<int>(next_arg(tok,it)); }
		else if (*it=="--cl-out")	{	set_output_fields(next_arg(tok,it), true); }
		else if (*it=="--cl-aff")						{	exact_arg(tok,it,0); cal_clustering("", true, true); }
		else if	(*it=="--cl-uid")						{	exact_arg(tok,it,1); cal_clustering(next_arg(tok,it), false, true); }
		else if (*it=="--mrca-of")						{	exact_arg(tok,it,1); cal_MRCA_file(next_arg(tok,it)); }
		else if (*it=="--oc")							{	exact_arg(tok,it,1); cal_OC_file(next_arg(tok,it)); }
		else if (*it=="--kin-all")						{	exact_arg(tok,it,0); cal_kinship_all(false); }
		else if (*it=="--q-kin-all")					{	exact_arg(tok,it,0); cal_kinship_all(true); }
		else if (*it=="--kin-of")						{	exact_arg(tok,it,1); cal_kinship_file(next_arg(tok,it),false); }
		else if (*it=="--q-kin-of")						{	exact_arg(tok,it,1); cal_kinship_file(next_arg(tok,it),true); }
		else if (*it=="--var-name")						{	exact_arg(tok,it,0); program.outf<<iData::default_vname; }
		else if (*it=="--founder-of")					{	exact_arg(tok,it,1); find_founders_of(next_arg(tok,it)); }
		else if (*it=="--sum-ped")						{	exact_arg(tok,it,0); write_summary(); }
		else if (*it=="--wr"    ||*it=="--write")		{	exact_arg(tok,it,1); set_output_fields(next_arg(tok,it), true); write_output(false); }
		else if (*it=="--wr-rai"||*it=="--write-rai")	{	exact_arg(tok,it,1); set_output_fields(next_arg(tok,it), true); write_output(true);  }
		else if (*it=="--wr-alt") wr_alternativeID=true;
		else if (*it=="--yob")
		{
			try {	fill_in_missing_yob( boost::lexical_cast<int>(next_arg(tok,it)) ); }
			catch (boost::bad_lexical_cast &) { exit_error("\nERROR: Wrong parameter for --yob"); }
		}
		else if (str_startsw(*it,"--mo-age=") || str_startsw(*it,"--ma-age="))
		{
			try
			{
				string input = it->substr(9);
				mo_age_lb = extract_int(input);
				extract_char(input);
				mo_age_ub = extract_int(input);
			}
			catch (...) { exit_error("\nERROR: Wrong input for "+*it); }
		}
		else if (str_startsw(*it,"--fa-age=") || str_startsw(*it,"--pa-age="))
		{
			try
			{
				string input = it->substr(9);
				fa_age_lb = extract_int(input);
				extract_char(input);
				fa_age_ub = extract_int(input);
			}
			catch (...) { exit_error("\nERROR: Wrong input for "+*it); }
		}
		else if (str_startsw(*it,"--id-del="))
		{
			string input = it->substr(9);
			if (input.empty()) exit_error("\nERROR: Wrong input for "+*it);
			id_delimiter_ = input;
		}
		else if (*it=="-Wno-single")	{ to_report_single=false; }
		else if (*it=="-Wno-gender")	{ to_report_gender=false; }
		else if (*it=="-Wno-re-inp")	{ to_report_re_inp=false; }
		else if (*it=="-Wno-subped")	{ to_report_subped=false; }
		else if (*it=="-Wno-dupUID")	{ to_report_dupUID=false; }
		else if (*it=="-w")				{ to_show_warnings=false; }
		else if (*it=="-a")				{ to_aggregate_con=true; }
		else if (iData::get_var_name(*it))		{ }
		else return false;
	}
	return true;
}

void FamilyDB::read_arguments(vector<string>& srce_opt)
{
	vector<string> tokens;
	bool has_read_file_opt=false;
	bool has_var_name_opt=false;
	for (auto &i:srce_opt)
	{
		if (str_startsw(i,"--ped") || str_startsw(i,"--read")) has_read_file_opt=true;
		if (str_startsw(i,"--var-name")) has_var_name_opt=true;
	}
	if (!has_read_file_opt && !has_var_name_opt)
	{
		lns<<showl<<"You didn't use --ped / --ped-nh, so read a Pedigree File from stdin. --create-dummy become true."<<flush_logger;
		tokens.push_back("--ped");
		tokens.push_back(label_stdin());
		if (!execute_command(tokens)) exit_error("\nERROR: wrong program option "+str_of_container(tokens,' '));
		tokens.clear();
	}
	for (size_t argi=1; argi<srce_opt.size(); ++argi)
	{
		if	(str_startsw(srce_opt[argi],"-"))
		{
			if (!tokens.empty()) { if (!execute_command(tokens)) exit_error("\nERROR: wrong program option "+str_of_container(tokens,' ')); }
			tokens.clear();
			tokens.push_back(srce_opt[argi]);
		}
		else
		{
			tokens.push_back(srce_opt[argi]);
		}
	}
	if (!tokens.empty()) { if (!execute_command(tokens)) exit_error("\nERROR: wrong program option "+str_of_container(tokens,' ')); }
}

void FamilyDB::set_output_fields(const string& instruction, bool clear_prev)
{
	if (clear_prev) outFields.clear();
	string expanded=iData::expand_var_names(instruction);
	while (!expanded.empty())
	{
		string s=extract_name(expanded);
		if (s.empty()) break;
		if (!expanded.empty()) extract_char(expanded);
		iData::degenerated_t field = iData::get_id_from(s);
		if (iData::vCode(field)==-1) exit_error("\nERROR: Unknown field name: "+s);
		if (iData::is_equal(field,iData::RCG)) outFields.insert(outFields.end(),inpFields.begin(),inpFields.end());
		else outFields.push_back(field);
	}
}

// be careful, created people will not be visited
void FamilyDB::traverse_IS(TravFuncPtrType fptr, void * UserData)
{
	for (each_element(rows_input_seq,it1))
		for (each_element(it1->second,it2)) 
			(this->*fptr)(family_db[it1->first]->members[*it2],UserData);
}

// garantee that all ancestors are visited before their offsprings
void FamilyDB::traverse_up(TravFuncPtrType fptr, void * UserData)
{
	set<string>		ped_already_written;
	vector<Family*> ped_output_sequence;
	for (each_element(rows_input_seq,it))
		if (!exist_element(ped_already_written,it->first))
		{
			ped_output_sequence.push_back(family_db[it->first]);
			ped_already_written.insert(it->first);
		}
	
	set_inds_status(iData::VST, false);
	for (each_element(ped_output_sequence,itf))
		for (each_element((*itf)->members,iti))
			_traverse_up_subfunc(iti->second, fptr, UserData);
	set_inds_status(iData::VST, false);
}

void FamilyDB::_traverse_up_subfunc(Individual* ind, TravFuncPtrType fptr, void * UserData)
{
	if (ind->get(iData::VST)==true) return;
	ind->set(iData::VST,true); // previously at the end aft fptr

	if (ind->pfa) if (ind->pfa->get(iData::VST)==false) _traverse_up_subfunc(ind->pfa,fptr,UserData);
	if (ind->pmo) if (ind->pmo->get(iData::VST)==false) _traverse_up_subfunc(ind->pmo,fptr,UserData);
	(this->*fptr)(ind,UserData);
}

// doesn't garantee that all ancestors are printed ahead of nextgens
void FamilyDB::traverse_dn(TravFuncPtrType fptr, void * UserData)
{
	set_inds_status(iData::VST, false);
	for (each_element_const(family_db,itf))
		for (each_element_const(itf->second->founders,iti))
			_traverse_dn_subfunc(iti->second, fptr, UserData);
	set_inds_status(iData::VST, false);
}

void FamilyDB::_traverse_dn_subfunc(Individual* ind, TravFuncPtrType fptr, void * UserData)
{
	if (ind->get(iData::VST)==true) return;
	ind->set(iData::VST,true); // previously after fptr

	(this->*fptr)(ind,UserData);
	for (each_element(ind->nextgen,iti))
	{
		if (iti->second->pfa->get(iData::VST)==false) _traverse_dn_subfunc(iti->second->pfa,fptr,UserData);
		if (iti->second->pmo->get(iData::VST)==false) _traverse_dn_subfunc(iti->second->pmo,fptr,UserData);
		_traverse_dn_subfunc(iti->second,fptr,UserData);
	}
}

// ------------- functions for cal_gen, not thread-safe -------------50417

void _find_doable_founder(Individual* ind, std::map<Individual*,int>& ancestors)
{
	if (ind->pfa)
	{
		if ((ind->pmo->pfa==NULL && ind->pmo->spouses.size()==1) || ind->pmo->get(iData::GEN)) _find_doable_founder(ind->pfa,ancestors);
		if ((ind->pfa->pfa==NULL && ind->pfa->spouses.size()==1) || ind->pfa->get(iData::GEN)) _find_doable_founder(ind->pmo,ancestors);
	}
	else
		++ancestors[ind];
}

void _find_founders_g(Individual* ind, std::map<Individual*,int>& ancestors, int gen)
{
	if (ind->pfa)
	{
		_find_founders_g(ind->pfa,ancestors,gen+1);
		_find_founders_g(ind->pmo,ancestors,gen+1);
	}
	else
	{
		if (!exist_element(ancestors,ind))	ancestors[ind]=gen;
		else if (ancestors[ind]<gen)		ancestors[ind]=gen;
	}
}

struct _traverse_out_cal_gen_result_type {
	Individual * ind;
	double       cgm;
	int          gen;
	bool		 upw;
} _traverse_out_cal_gen_result;

void _traverse_out_cal_gen(Individual* ind, int gen, int flag) // 0=new_fam 1=init 2=compute 3=finish
{
	// data for each compute
	static int		num_call=0;	// number of times func is called
	static int		num_same=0;	// number of times new==old value
	static int		num_chng=0;	// number of times new!=old value
	static int		num_init=0;	// number of times old value is 0
	static int		mar_mtch=0;	// marriage that generations match
	static double	mar_indx=0;	// cross generation marriage score
	
	// data for each init
	static double			best_score = DBL_MAX;	// the best among several run for the same IND
	static int				best_input = 0;			// input  gen associated with the best_score
	static Individual *		suggestion = NULL;		// output ind associated with the best_score
	
	// data for each fam
	static set<Individual*> doable_sps;

	if (ind==NULL || gen<=0 || flag<0 || flag>3) exit_error("wrong parameter for _traverse_out_cal_gen().");
	
	if		(flag==0)
	{
		num_call=0; num_same=0; num_chng=0; num_init=0; mar_mtch=0; mar_indx=0;
		best_score = DBL_MAX;
		best_input = 0;
		suggestion = NULL;
		doable_sps.clear();
	}
	else if	(flag==1)
	{
		num_call=0; num_same=0; num_chng=0; num_init=0; mar_mtch=0; mar_indx=0;
		best_score = DBL_MAX;
		best_input = 0;
		suggestion = NULL;
	}
	else if (flag==3)
	{
		// update best gen if necessary
		if (mar_indx<=best_score)
		{
			best_score = mar_indx;
			best_input = gen;
			_traverse_out_cal_gen_result.upw = false;
		}
		else
		{
			_traverse_out_cal_gen_result.upw = true;
		}
		_traverse_out_cal_gen_result.gen = best_input;
		_traverse_out_cal_gen_result.cgm = mar_indx;
		
		// make suggestion for next ind
		suggestion=NULL;
		while (!doable_sps.empty())
		{
			map<Individual*,int> doable_fnd; // trackable founder of doable spouses
			for (each_element(doable_sps,it)) _find_doable_founder(*it,doable_fnd);
			if (doable_fnd.empty())
			{
				int			min_gen=INT_MAX;
				Individual* top_sps=NULL;
				for (each_element(doable_sps,it))
					for (each_element((*it)->spouses,its))
						if (its->second->get(iData::GEN) && its->second->get(iData::GEN)<min_gen)
						{ min_gen=its->second->get(iData::GEN); top_sps=*it; }
				map<Individual*,int> top_sps_fndr;
				_find_founders_g(top_sps,top_sps_fndr,0);
				for (each_element(top_sps_fndr,it)) _traverse_out_cal_gen(it->first,min_gen-it->second,2);
			}
			else
			{
				int	best_cnt = 0;
				for (each_element(doable_fnd,it)) if (it->second>best_cnt) { best_cnt=it->second; suggestion=it->first; }
				break;
			}
		}
		_traverse_out_cal_gen_result.ind = suggestion;
		
		num_call=0; num_same=0; num_chng=0; num_init=0; mar_mtch=0; mar_indx=0;
	}
	else if (flag==2)
	{
		++num_call;
		if (ind->get(iData::GEN)==gen)	{ ++num_same; return; }
		if (ind->get(iData::GEN)==0) ++num_init; else ++num_chng;
		
		// set self
		ind->set(iData::GEN,gen);
		ind->set(iData::VST,true);
		doable_sps.erase(ind);
		
		for (each_element(ind->spouses,its))
		{
			// set spouses
			if (its->second->pfa==NULL && its->second->spouses.size()==1 && its->second->get(iData::VST)==false)
			{
				_traverse_out_cal_gen(its->second,gen,flag);
				continue;
			}
			else if (its->second->get(iData::GEN)!=0 && its->second->get(iData::VST)==false)
			{
				if (its->second->get(iData::GEN)==gen)	mar_mtch+=1;
				else									mar_indx+=pow(100,abs(its->second->get(iData::GEN)-gen)-1);
			}
			
			// set children
			if (its->second->get(iData::GEN)!=0)
			{
				int nxt = max( int(its->second->get(iData::GEN)), gen) +1;
				for (each_element(ind->nextgen,itc))
				{
					if (ind->get(iData::SEX)==kMAN && itc->second->pmo==its->second) _traverse_out_cal_gen(itc->second,nxt,flag);
					if (ind->get(iData::SEX)==kWMN && itc->second->pfa==its->second) _traverse_out_cal_gen(itc->second,nxt,flag);
				}
			}
			else
			{
				doable_sps.insert(its->second);
			}
		}
		ind->set(iData::VST,false);
	}
	else
		exit_error("wrong parameter for _traverse_out_cal_gen().");
}

// check Gs results. It must be garantee that all INDs are checked.
void _traverse_dn_chk_gen(Individual* ind, int gen, set<string>& CgmStrng, double& CgmScore, int& ErrOrNul)
{
	if (gen==0)
	{
		lns<<writew<< "Person "<<ind->get(iData::IID)<<" G null."<<flush_logger;
		++ErrOrNul;
		return;
	}
	else if (gen>0)
	{
		// Check self
		if (ind->get(iData::GEN)!=gen)
		{
			lns<<writew<< "Person "<<ind->get(iData::IID)<<" G wrong/null."<<flush_logger;
			++ErrOrNul;
			return;
		}
		
		int g_matched=0;
		for (each_element(ind->spouses,its))
		{
			if (its->second->get(iData::GEN)!=0)
			{
				// check: founders have the same G as at least one of their spouses
				if (its->second->get(iData::GEN)==gen)
					++g_matched;
				
				// check: CgmStrng. It's normal to have some but not too many; cross many generations is abnormal.
				if (its->second->get(iData::GEN)-gen>0)
				{
					int gen_dif = its->second->get(iData::GEN)-gen;
					string marr_str = ind->get(iData::IID)+" "+its->second->get(iData::IID) + " has " + itos(gen_dif) + " generation gap";
					if (!exist_element(CgmStrng,marr_str))
					{
						CgmScore += pow(100,gen_dif-1);
						CgmStrng.insert( marr_str );
					}
				}
				
				// check: children have the G of max(parents)+1
				int nxt = max( int(its->second->get(iData::GEN)), gen) +1;
				for (each_element(ind->nextgen,itc))
				{
					if (ind->get(iData::SEX)==kMAN && itc->second->pmo==its->second) _traverse_dn_chk_gen(itc->second,nxt,CgmStrng,CgmScore,ErrOrNul);
					if (ind->get(iData::SEX)==kWMN && itc->second->pfa==its->second) _traverse_dn_chk_gen(itc->second,nxt,CgmStrng,CgmScore,ErrOrNul);
				}
			}
			else
			{
				lns<<writew<< "Person "<<its->second->get(iData::IID)<<" G null."<<flush_logger;
				++ErrOrNul;
				return;
			}
		}
//		if (ind->pfa==NULL && !ind->spouses.empty() && g_matched==0)
//		{
//			lns<<writew<< "Person "<<ind->get(iData::IID)<<" G not match with any spouse(s).";
//			++ErrOrNul;
//			return;
//		}
	}
	else
	{
		exit_error("gen should be >0 for _traverse_dn_chk_gen().");
	}
}

void _traverse_dn_cal_des(Individual* ind)
{
	int max_gen = 0; // max No. of generations
	int tot_des = 0; // total No. of descendant
	for (each_element(ind->nextgen,it))
	{
		_traverse_dn_cal_des(it->second);
		if (it->second->get(iData::MXG)>max_gen) max_gen=it->second->get(iData::MXG);
		tot_des += it->second->get(iData::DES);
	}
	ind->set(iData::MXG, max_gen + !ind->nextgen.empty());
	ind->set(iData::DES, tot_des +  ind->nextgen.size());
}

void FamilyDB::cal_generation_fast() // fast and dirty method
{
	lns<<showl<<"Calculate Gs fast-and-dirty ... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf))
		for (each_element(itf->second->founders,iti))
			_traverse_dn_cal_des(iti->second);
	lns.endsub();
}

void FamilyDB::cal_generation()
{
	lns<<showl<<"Calculate generation number ... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf))
	{
		// sort the founders by the length of lineage among descendants
		// no need to test whether spouse is doable because they will be the next
		int max_len=0;
		std::multimap<int, Individual*, std::greater<int> > fndr_map;
		for (each_element(itf->second->founders,iti))
		{
			_traverse_dn_cal_des(iti->second);
			fndr_map.insert(pair<int,Individual*>(iti->second->get(iData::MXG),iti->second));
			if (iti->second->get(iData::MXG)>max_len) max_len=iti->second->get(iData::MXG);
		}
		
		// Although robust to pedigrees w/ multiple clusters, the G of other clusters not always starts from 1
		for (;;) // for each cluster of individuals
		{
			// find the top founder w/o a G
			Individual* ind=NULL;
			for (each_element(fndr_map,ita))
				if (ita->second->get(iData::GEN)==0) { ind=ita->second; break; }
			if (ind==NULL) break;
			
			// assign G=1000 to this founder
			_traverse_out_cal_gen(ind,1000,0);
			_traverse_out_cal_gen(ind,1000,2);
			_traverse_out_cal_gen(ind,1000,3);
			
			// assign Gs following suggestions
			vector<Individual*> assign_ind;
			vector<int>			assign_gen;
			while ( _traverse_out_cal_gen_result.ind!=NULL )
			{
				ind=_traverse_out_cal_gen_result.ind;
				_traverse_out_cal_gen(ind,1,1);
				for (int gen=1000-max_len; gen<=1000+max_len; ++gen)
				{
					_traverse_out_cal_gen(ind,gen,2);
					_traverse_out_cal_gen(ind,gen,3);
					if (_traverse_out_cal_gen_result.upw) break;
				}
				_traverse_out_cal_gen(ind,_traverse_out_cal_gen_result.gen,2);
				assign_ind.push_back(ind);
				assign_gen.push_back(_traverse_out_cal_gen_result.gen);
			}
			
			/*/ fine tunning assign_gen. Removed because results are even worse or not helpful.
			set< pair<Individual*,int> > discarded;
			for (int to_cont=1; to_cont; )
			{
				to_cont=0;
				for (int i=0; i<(int)assign_gen.size(); ++i)
				{
					Individual* ind = assign_ind[i];
					int         gen = assign_gen[i];
					_traverse_out_cal_gen(ind,gen  ,1);
					_traverse_out_cal_gen(ind,gen-1,2);
					_traverse_out_cal_gen(ind,gen-1,3);
					_traverse_out_cal_gen(ind,gen+1,2);
					_traverse_out_cal_gen(ind,gen+1,3);
					_traverse_out_cal_gen(ind,gen  ,2);
					_traverse_out_cal_gen(ind,gen  ,3);
					if (_traverse_out_cal_gen_result.gen!=gen)
					{
						if  (discarded.find  (pair<Individual*,int>(ind,gen))!=discarded.end()) continue;
						else discarded.insert(pair<Individual*,int>(ind,gen));
						assign_gen[i]=_traverse_out_cal_gen_result.gen;
						_traverse_out_cal_gen(ind,assign_gen[i],2);
						to_cont=1;
						break;
					}
				}
			}//*/
		}
		
		// check Gs results. It must be garantee that all INDs are checked.
		set<string> CgmStrng;
		double		CgmScore=0;
		int			ErrOrNul=0;
		for (each_element(itf->second->founders,iti))
			if (iti->second->get(iData::GEN)!=0) _traverse_dn_chk_gen(iti->second,iti->second->get(iData::GEN),CgmStrng,CgmScore,ErrOrNul);
		if (CgmScore)	lns<<writel<<"Cross generation marriage score for fam "<<itf->first<<" : "<<CgmScore<<flush_logger;
		if (ErrOrNul)	lns<<writee<<"No. of generation number errors for fam "<<itf->first<<" : "<<ErrOrNul<<flush_logger;
		
		// downcode to [1,inf), and calculate summary statistics
		int subtract = mnp_inds_values(itf->second,iData::GEN)-1;
		for (each_element(itf->second->members,iti))
		{
			int g= iti->second->get(iData::GEN);
			if (g) iti->second->set(iData::GEN,g-subtract);
		}
	}
	lns.endsub();
}

// ------------- end of functions for cal_gen -------------

// some functions that could be used by traverse
void FamilyDB::_print_line(Individual* ind, void * UserData)
{
	static const Individual::Relationship null_rs;
	for (size_t i=0; i<outFields.size(); ++i)
	{
		if (i) program.outf<<DLMTR;
		
		if (iData::is_equal(outFields[i],iData::REL))
		{
			const Individual::Relationship& r = ind->rrs.empty() ? null_rs : (ind->rrs.begin()->second);
			program.outf<<r.code<<DLMTR<<r.p_m<<DLMTR<<r.relation;
		}
		else if (iData::is_equal(outFields[i],iData::DNF))
		{
			if (ind->pfa)	program.outf<< ind->pfa->get(iData::DNI);
			else			program.outf<< string_no_id_;
		}
		else if (iData::is_equal(outFields[i],iData::DNM))
		{
			if (ind->pfa)	program.outf<< ind->pmo->get(iData::DNI);
			else			program.outf<< string_no_id_;
		}
		else if (iData::is_equal(outFields[i],iData::OTH))
		{
			if (ind->xtr.empty())
				for (size_t i=1;i<xtrColTtl.size();++i) program.outf<<DLMTR;
			else
				print_map_second(ind->xtr,program.outf,DLMTR);
		}
		else if (iData::is_equal(outFields[i],iData::IID))
		{
			if (wr_alternativeID)
			{
				string& alt = ind->get(iData::ALT);
				if (alt==string_no_id_)	program.outf<< ind->get(iData::IID);
				else					program.outf<< alt;
			}
			else
				program.outf<< ind->bf_get(outFields[i]);
		}
		else if (iData::is_equal(outFields[i],iData::FTH))
		{
			if (wr_alternativeID && ind->pfa)
			{
				string& alt = ind->pfa->get(iData::ALT);
				// std::cerr<<"print father of "<<ind->get(iData::IID)<<" alt="<<alt<<" iid="<<ind->get(iData::FTH)<<endl;
				if (alt==string_no_id_)	program.outf<< ind->get(iData::FTH);
				else					program.outf<< alt;
			}
			else
				program.outf<< ind->bf_get(outFields[i]);
		}
		else if (iData::is_equal(outFields[i],iData::MTH))
		{
			if (wr_alternativeID && ind->pmo)
			{
				string& alt = ind->pmo->get(iData::ALT);
				// std::cerr<<"print mother of "<<ind->get(iData::IID)<<" alt="<<alt<<" iid="<<ind->get(iData::MTH)<<endl;
				if (alt==string_no_id_)	program.outf<< ind->get(iData::MTH);
				else					program.outf<< alt;
			}
			else
				program.outf<< ind->bf_get(outFields[i]);
		}
		else
		{
			program.outf<< ind->bf_get(outFields[i]);
		}
	}
	program.outf<<endl;
}

void FamilyDB::_print_header()
{
	if (outFields.empty()) exit_error("\nERROR: Output fields not defined.");
	iData::vText(iData::OTH) = str_of_map_second(xtrColTtl,DLMTR);
	lns<<"output ";
	for (size_t i=0; i<outFields.size(); ++i)
	{
		if (i) program.outf<<DLMTR;
		program.outf<<iData::vText(outFields[i]);
		lns<<' '<<iData::v_ID[iData::vCode(outFields[i])]; //lns<<' '<<iData::vText(outFields[i]);
	}
	lns<<flush_logger;
	program.outf<<endl;
}

void FamilyDB::write_summary()
{
	cal_generation_fast();
	set<string>		ped_already_written;
	vector<Family*> ped_output_sequence;
	for (each_element(rows_input_seq,it))
		if (!exist_element(ped_already_written,it->first))
		{
			ped_output_sequence.push_back(family_db[it->first]);
			ped_already_written.insert(it->first);
		}
	lns<<showl<<"Write summary stat for ped. ... "<<flush_logger;
	lns.flush_log();
	cerr<<endl;
 	
	program.outf<<"## Bits    = 2 * N(nonfounders) - N(founders)\n";
	program.outf<<"## NumAff  = number of affected individuals\n";
	program.outf<<"## NumSep  = number of separated individuals\n";
	program.outf<<"## NumClt  = number of clusters\n";
	program.outf<<"## MaxChd  = max number of children per person\n";
	program.outf<<"## MaxSib  = max number of full sibs\n";
	program.outf<<"## MaxMarr = max number of marriages\n";
	program.outf<<"## NumGen  = number of generations (fast-and-dirty)\n";
	program.outf<<"## AddInd  = number of individuals added by the program\n";
	program.outf<<"## AddDmy  = number of individuals created by the program\n";
	program.outf<<"## Consang = number of consanguineous couples\n";
	program.outf<<"## Circle  = whether circle(s) exist\n";
	size_t max_pid_len=0;
	for (each_element(family_db,itf))
		if (itf->first.size()>max_pid_len) max_pid_len=itf->first.size();
	if (max_pid_len<8)	program.outf<<"#PedID";
	else				program.outf<<"#PedigreeID";
	program.outf<<DLMTR<<"NumInd"<<DLMTR<<"Bits"<<DLMTR<<"NumAff"<<DLMTR<<"NumSep"<<DLMTR<<"NumClt"
				<<DLMTR<<"MaxChd"<<DLMTR<<"MaxSib"<<DLMTR<<"MaxMarr"<<DLMTR<<"NumGen"<<DLMTR<<"AddInd"<<DLMTR<<"AddDmy"
				<<DLMTR<<"Consang"<<DLMTR<<"Circle"<<'\n';
	for (each_element(ped_output_sequence,itf))
	{
		Family * f = *itf;
		int num_aff = 0;
		int num_sep = 0;
		int num_dum = 0;
		int num_abp = 0;
		int max_chd = 0;
		int max_sib = 0;
		int max_mar = 0;
		int tot_gen = 0;
		int bits	= 2 * f->nonfounders.size() - f->founders.size();
		for (each_element(f->members,iti))
		{
			Individual* ind=iti->second;
			num_aff += is_affected(ind);
			if (ind->get(iData::SEP)==true) ++num_sep;
			if (ind->get(iData::DUM)==true) ++num_dum;
			if (ind->get(iData::ABP)==true) ++num_abp;
			if ((int)ind->nextgen.size()>max_chd) max_chd=ind->nextgen.size();
			if ((int)ind->spouses.size()>max_mar) max_mar=ind->spouses.size();
			if (ind->get(iData::MXG)>tot_gen) tot_gen=ind->get(iData::MXG);
			for (each_element(ind->spouses,its))
			{
				int num_sib=0;
				for (each_element(ind->nextgen,itc))
				{
					if (ind->get(iData::SEX)==kMAN && itc->second->pmo==its->second) ++num_sib;
					if (ind->get(iData::SEX)==kWMN && itc->second->pfa==its->second) ++num_sib;
				}
				if (num_sib>max_sib) max_sib=num_sib;
			}
		}
		program.outf<< f->pid << DLMTR
					<< f->members.size() << DLMTR
					<< bits << DLMTR
					<< num_aff << DLMTR
					<< num_sep << DLMTR
					<< f->num_clt << DLMTR
					<< max_chd << DLMTR
					<< max_sib << DLMTR
					<< max_mar << DLMTR
					<< tot_gen << DLMTR
					<< num_abp << DLMTR
					<< num_dum << DLMTR
					<< f->consang.size()/2 << DLMTR
					<< (f->circles ? "yes" : "-")
					<< endl;
	}
}

void FamilyDB::write_output(bool rows_in_input_seq)
{
	if (outFields.empty()) set_output_fields("l5",false);
	if (iData::vec_contain(outFields,iData::DNP) ||
		iData::vec_contain(outFields,iData::DNI) ||
		iData::vec_contain(outFields,iData::DNU) ||
		iData::vec_contain(outFields,iData::L5D) ||
		iData::vec_contain(outFields,iData::L6D) ||
		iData::vec_contain(outFields,iData::L7D) ||
		iData::vec_contain(outFields,iData::L8D) ||
		iData::vec_contain(outFields,iData::L9D)) downcode_IDs();
	if (iData::vec_contain(outFields,iData::SEP)) check_separation();
	if (iData::vec_contain(outFields,iData::AFD)) cal_aff_fdr();
	if (iData::vec_contain(outFields,iData::REL)) cal_relationship();
	if (iData::vec_contain(outFields,iData::GEN)) cal_generation();
	if (iData::vec_contain(outFields,iData::INF)) cal_informativeness();
	if (iData::vec_contain(outFields,iData::CLT)) cal_clustering("", false, true); // prv check_clusters()
	else
	{
		lns<<showl<<"Write pedigree file ... "<<flush_logger;
		cerr<<endl;
		_print_header();
		if (rows_in_input_seq)	traverse_IS(&FamilyDB::_print_line,NULL);
		else					traverse_up(&FamilyDB::_print_line,NULL);
	}
}

// rm if trb[status]=value, all others remain, so may produce separated persons floating around
void FamilyDB::del_inds_by_status(iData::trb_t status, boost::logic::tribool value)
{
	lns<<showl<<"Remove individuals based on status["<<iData::vText(status)<<"]="<<value<<" ... "<<flush_logger;
	lns.sub();
	int deleted=0;
	for (each_element(family_db,itf))
	{
		Family* f=itf->second;
remove_individuals_begin:
		for (each_element(f->members,iti))
		{
			Individual* i=iti->second;
			if (i->get(status)==value)
			{
				// delete self
				uniqid_db.erase(i->get(iData::UID));
				f->members.erase(i->get(iData::IID));
				f->founders.erase(i->get(iData::IID));
				f->nonfounders.erase(i->get(iData::IID));
				if (i->pfa)
				{
					i->pfa->nextgen.erase(i->get(iData::IID));
					i->pmo->nextgen.erase(i->get(iData::IID));
				}
				// delete links with spouse
				for (each_element(i->spouses,iti2))
					iti2->second->spouses.erase(i->get(iData::IID));
				// delete links with nextgen
				if (i->get(iData::SEX)==kMAN)
				{
					for (each_element(i->nextgen,iti2))
					{
						iti2->second->pmo->nextgen.erase(iti2->first);
						iti2->second->pfa=NULL; iti2->second->set(iData::FTH, string_no_id_);
						iti2->second->pmo=NULL; iti2->second->set(iData::MTH, string_no_id_);
						f->nonfounders.erase(iti2->second->get(iData::IID));
						f->founders[iti2->second->get(iData::IID)]=iti2->second;
					}
				}
				if (i->get(iData::SEX)==kWMN)
				{
					for (each_element(i->nextgen,iti2))
					{
						iti2->second->pfa->nextgen.erase(iti2->first);
						iti2->second->pfa=NULL; iti2->second->set(iData::FTH, string_no_id_);
						iti2->second->pmo=NULL; iti2->second->set(iData::MTH, string_no_id_);
						f->nonfounders.erase(iti2->second->get(iData::IID));
						f->founders[iti2->second->get(iData::IID)]=iti2->second;
					}
				}
				// finished
				lns<<writel<<"person "<<f->pid<<id_delimiter_<<i->get(iData::IID)<<" (GTN="<<i->get(iData::GTN)<<",AFF="<<i->get(iData::AFF)<<") removed"<<flush_logger;
				delete i;
				++deleted;
				goto remove_individuals_begin;
			}
		}
	}
	lns<<cshowl<<"totally "<<deleted<<" persons removed."<<flush_logger;
	lns.endsub();
}

void FamilyDB::check_separation()
{
	int tot_num=0;
	lns<<showl<<"Check for separated persons ... ";
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->founders,iti))
		{
			Individual* i=iti->second;
			if (i->nextgen.empty()) {	i->set(iData::SEP,true); ++tot_num; }
			else						i->set(iData::SEP,false);
		}
		for (each_element(itf->second->nonfounders,iti))
			iti->second->set(iData::SEP,false);
	}
	if (tot_num) lns<<"Found "<<tot_num<<flush_logger;
}

void FamilyDB::check_prb()
{
	lns<<showl<<"Check for >1 probands/ped   ... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf))
	{
		int n_prb=0;
		for (each_element(itf->second->members,itm))
			if (is_proband(itm->second)) ++n_prb;
		if (n_prb!=1) lns<<showe<<"Family "<<itf->second->pid<<" has "<<n_prb<<" probands."<<flush_logger;
	}
	lns.endsub();
}

void FamilyDB::check_aff()
{
	lns<<showl<<"Check for affection!=0/1/2  ... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->members,itm))
		{
			int aff = itm->second->get(iData::AFF);
			if (aff!=0 && aff!=1 && aff!=2)
				lns<<showe<<"Family "<<itf->second->pid<<" person "<<itm->second->get(iData::IID)<<" aff="<<aff<<flush_logger;
		}
	}
	lns.endsub();
}

void FamilyDB::check_lia()
{
	lns<<showl<<"Check for negative liability... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf))
		for (each_element(itf->second->members,itm))
			if (itm->second->get(iData::LIA)<0)
				lns<<showe<<"Family "<<itf->second->pid<<" person "<<itm->second->get(iData::IID)<<" lia="<<itm->second->get(iData::LIA)<<flush_logger;
	lns.endsub();
}

void FamilyDB::check_yob()
{
	lns<<showl<<"Check for improbable YoB    ... "<<flush_logger;
	lns.sub();
	int min_deliver_age = 99;
	int max_deliver_age = 0;
	int min_fathers_age = 99;
	int max_fathers_age = 0;
	int num_years_a_gen = 0;
	int num_generations = 0;
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->members,itm))
			if (itm->second->get(iData::YOB)<=0 || itm->second->get(iData::YOB)>9999 || std::isnan(itm->second->get(iData::YOB))) itm->second->set(iData::YOB,9999);
		for (each_element(itf->second->members,itm))
		{
			if (itm->second->get(iData::YOB)==9999) continue;
			for (each_element(itm->second->nextgen,itng)) 
			{
				if (itng->second->get(iData::YOB)==9999) continue;
				int age = itng->second->get(iData::YOB) - itm->second->get(iData::YOB);
				if (itm->second->get(iData::SEX)==kWMN)
				{
					if (age<mo_age_lb || age>mo_age_ub)
					{
						lns<<showe<<itf->second->pid<<" - "<<itm->second->get(iData::IID)<<" improbable YoB: deliver a child at the age of "<<age<<flush_logger;
						continue;
					}
					if (age < min_deliver_age) min_deliver_age = age;
					if (age > max_deliver_age) max_deliver_age = age;
					num_years_a_gen += age;
					++num_generations;
				}
				else if (itm->second->get(iData::SEX)==kMAN) 
				{
					if (age<fa_age_lb || age>fa_age_ub)
					{
						lns<<showe<<itf->second->pid<<" - "<<itm->second->get(iData::IID)<<" improbable YoB: became a father at the age of "<<age<<flush_logger;
						continue;
					}
					if (age < min_fathers_age) min_fathers_age = age;
					if (age > max_fathers_age) max_fathers_age = age;
					num_years_a_gen += age;
					++num_generations;
				}
			}
		}
	}
	if (max_deliver_age && max_fathers_age)
	{
		lns<<showl<<"Women deliver babies at age "<<min_deliver_age<<" to "<<max_deliver_age<<flush_logger;
		lns<<showl<<"Men become fathers at age "<<min_fathers_age<<" to "<<max_fathers_age<<flush_logger;
		lns<<showl<<"One generation is about "<<num_years_a_gen*1.0/num_generations<<" years"<<flush_logger;
	}
	lns.endsub();
}

void FamilyDB::fill_in_missing_yob(int default_yob_for_probands)
{
	lns<<showl<<"Infer missing year-of-birth ... "<<flush_logger;
	lns.sub();
	// weight for son/daughter/father/mother is always 1
	int min_deliver_age = 18;
	int max_deliver_age = 45;
	int max_fathers_age = 55;
	int num_years_a_gen = 30;
	map<string,string> min_max_err;
	for (each_element(family_db,itf))
	{
		// if nobody has yob, then set default yob to all probands (even >1)
		int count=0;
		for (each_element(itf->second->members,itm))
		{
			if (itm->second->get(iData::YOB)==0) itm->second->set(iData::YOB,9999);
			if (itm->second->get(iData::YOB)!=9999) ++count;
		}
		if (count==0)
		{
			for (each_element(itf->second->members,itm))
				if (is_proband(itm->second))
				{
					itm->second->set(iData::YOB,default_yob_for_probands);
					++count;
				}
		}
		if (count==0) 
		{
			lns<<showe<<"Family "<<itf->second->pid<<" has no probands and all members have missing YoB; Infer YoB not done."<<flush_logger;
			continue;
		}
		
		// fill in missing yob
		
		double wt_input = 1; // weight for persons whose YoB is known by input
		double wt_infer = 0; // weight for persons whose YoB is known by infer
		// beginning wt_infer=0 so that they are not considered at the 1st try
		// but will have a non-zero weight after the try when no change is made
		
		for (double max_num_guess=0; ;max_num_guess=0)
		{
			// the most informative person (IID and inferred YoB), correspond to max_num_guess
			string max_ind_IID;
			int max_yob_guess=9999;
			
			for (each_element(itf->second->members,itm))
			{
				// prepare data
				Individual* ind=itm->second;
				if (ind->get(iData::YOB)!=9999) continue;
				int min_yob = 0, max_yob = 9999, avg_yob=0;
				double yob_guess = 0;
				double num_guess = 0;
				if (ind->pfa)
				{
					if (ind->pfa->get(iData::YOB)!=9999)
					{
						int new_min = ind->pfa->get(iData::YOB) + min_deliver_age;
						int new_max = ind->pfa->get(iData::YOB) + max_fathers_age;
						if (new_min > min_yob) min_yob = new_min;
						if (new_max < max_yob) max_yob = new_max;
						double weight = ind->pfa->get(iData::YBI)==true ? wt_infer : wt_input ;
						yob_guess += (ind->pfa->get(iData::YOB) + num_years_a_gen) * weight;
						num_guess += weight;
					}
					if (ind->pmo->get(iData::YOB)!=9999)
					{
						int new_min = ind->pmo->get(iData::YOB) + min_deliver_age;
						int new_max = ind->pmo->get(iData::YOB) + max_deliver_age;
						if (new_min > min_yob) min_yob = new_min;
						if (new_max < max_yob) max_yob = new_max;
						double weight = ind->pmo->get(iData::YBI)==true ? wt_infer : wt_input ;
						yob_guess += (ind->pmo->get(iData::YOB) + num_years_a_gen) * weight;
						num_guess += weight;
					}
					set<string> bro_sis; // including half-bro half-sis
					for (each_element(ind->pfa->nextgen,itng))
					{
						if (exist_element(bro_sis,itng->second->get(iData::IID))) continue;
						bro_sis.insert(itng->second->get(iData::IID));
						if (itng->second->get(iData::YOB)==9999) continue;
						double weight = itng->second->get(iData::YBI)==true ? wt_infer : wt_input ;
						yob_guess += itng->second->get(iData::YOB) * weight;
						num_guess += weight;
					}
					for (each_element(ind->pmo->nextgen,itng))
					{
						if (exist_element(bro_sis,itng->second->get(iData::IID))) continue;
						bro_sis.insert(itng->second->get(iData::IID));
						if (itng->second->get(iData::YOB)==9999) continue;
						double weight = itng->second->get(iData::YBI)==true ? wt_infer : wt_input ;
						yob_guess += itng->second->get(iData::YOB) * weight;
						num_guess += weight;
					}
				}
				for (each_element(ind->nextgen,itng)) 
				{
					if (itng->second->get(iData::YOB)!=9999)
					{
						int new_min, new_max;
						if		(ind->get(iData::SEX)==kMAN) {
							new_min = itng->second->get(iData::YOB) - max_fathers_age;
							new_max = itng->second->get(iData::YOB) - min_deliver_age;
						}
						else if (ind->get(iData::SEX)==kWMN) {
							new_min = itng->second->get(iData::YOB) - max_deliver_age;
							new_max = itng->second->get(iData::YOB) - min_deliver_age;
						}
						else { // this is unusual
							new_min = itng->second->get(iData::YOB) - (max_fathers_age + max_deliver_age)/2;
							new_max = itng->second->get(iData::YOB) - min_deliver_age;
						}
						if (new_min > min_yob) min_yob = new_min;
						if (new_max < max_yob) max_yob = new_max;
						double weight = itng->second->get(iData::YBI)==true ? wt_infer : wt_input ;
						yob_guess += (itng->second->get(iData::YOB) - num_years_a_gen) * weight;
						num_guess += weight;
					}
				}
				for (each_element(ind->spouses,itsp))
				{
					if (itsp->second->get(iData::YOB)!=9999) 
					{
						double weight = itsp->second->get(iData::YBI)==true ? wt_infer : wt_input ;
						yob_guess += itsp->second->get(iData::YOB) * weight;
						num_guess += weight;
					}
				}
				if (num_guess) avg_yob = round(yob_guess/num_guess);
				
				// decide whether YoB can be guessed
				if (min_yob == 0 || max_yob==9999 || avg_yob==0)
					continue;
				else if (min_yob > max_yob)
				{
					if (min_yob-max_yob < 5)
					{
						min_yob = max_yob = round((min_yob+max_yob)/2.0);
					}
					else
					{
						min_max_err[ind->get(iData::UID)] = itos(min_yob) + ">" + itos(max_yob) ;
						continue;
					}
				}
				
				// whether it's more informative
				if (num_guess > max_num_guess) 
				{
					max_num_guess = num_guess;
					if		(avg_yob < min_yob) max_yob_guess=min_yob;
					else if (avg_yob > max_yob) max_yob_guess=max_yob;
					else						max_yob_guess=avg_yob;
					max_ind_IID = ind->get(iData::IID);
				}
			}
			if (max_num_guess)
			{
				itf->second->members[max_ind_IID]->set(iData::YOB, max_yob_guess);
				itf->second->members[max_ind_IID]->set(iData::YBI, true);
			}
			else 
			{
				if (wt_infer) break;
				else wt_infer=.5;
			}
		}
	}
	for (each_element(min_max_err,itmme))
		lns<<showe<<"Infer YoB failed for "<<itmme->first<<" : min>max ("<<itmme->second<<")."<<flush_logger;
	lns.endsub();
}

void FamilyDB::check_mzt()
{
	lns<<showl<<"Check monozygosity twins    ... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf)) 
	{
		map<int, IndMap> twins;
		for (each_element(itf->second->members, itm))
		{
			string twin_str = itm->second->get(iData::MZT);
			if (twin_str!="") 
			{
				try 
				{
					int twin_id=boost::lexical_cast<int>(twin_str);
					if (twin_id>0) 
					{
						twins[twin_id][itm->second->get(iData::IID)]=itm->second;
					}
					else
						lns<<showe<<"fam "<<itf->second->pid<<" twin_ID "<<twin_str<<" is not greater than 0."<<flush_logger;
				}
				catch (boost::bad_lexical_cast &) { lns<<showe<<"fam "<<itf->second->pid<<" twin_ID "<<twin_str<<" is not an integer."<<flush_logger; }
			}
		}
		for (each_element(twins,itt)) {
			if (itt->second.size()==1) { lns<<showe<<"fam "<<itf->second->pid<<" twin_ID "<<itt->first<<" has  no identical sibs."<<flush_logger; continue; }
			if (itt->second.size()>2)  { lns<<showe<<"fam "<<itf->second->pid<<" twin_ID "<<itt->first<<" have >1 identical sibs."<<flush_logger; continue; }
			Individual *ind1 = itt->second.begin()->second;
			Individual *ind2 = itt->second.rbegin()->second;
			if (ind1->get(iData::SEX)!=ind2->get(iData::SEX)) lns<<showe<<"fam "<<itf->second->pid<<" twins "<<itt->first<<" are of opposite sex."<<flush_logger;
			if (ind1->get(iData::SEX)!=kMAN && ind1->get(iData::SEX)!=kWMN) lns<<showe<<"fam "<<itf->second->pid<<" ind "<<ind1->get(iData::IID)<<" has an illegal sex code for twins."<<flush_logger;
			if (ind2->get(iData::SEX)!=kMAN && ind2->get(iData::SEX)!=kWMN) lns<<showe<<"fam "<<itf->second->pid<<" ind "<<ind2->get(iData::IID)<<" has an illegal sex code for twins."<<flush_logger;
			if (ind1->get(iData::FTH)!=ind2->get(iData::FTH)) lns<<showe<<"fam "<<itf->second->pid<<" twins "<<itt->first<<" have different fathers."<<flush_logger;
			if (ind1->get(iData::MTH)!=ind2->get(iData::MTH)) lns<<showe<<"fam "<<itf->second->pid<<" twins "<<itt->first<<" have different mothers."<<flush_logger;
		}
	}
	lns.endsub();
}

void FamilyDB::_is_informative(Individual * ind)
{
	int inf_child=0;			// number of informative child
	Individual * spouse=NULL;	// spouse for the last informative child
	for (each_element(ind->nextgen,iti))
	{
		Individual * i = iti->second;
		if (boost::indeterminate(i->get(iData::INF)))	_is_informative(i);
		if (i->get(iData::INF)==true)
		{
			++inf_child;
			if		(ind->get(iData::SEX)==kMAN) spouse=i->pmo;
			else if (ind->get(iData::SEX)==kWMN) spouse=i->pfa;
			else exit_error("\nERROR: Unknown gender of individual "+ind->get(iData::IID));
		}
	}
	
	int aff_sib=0;				// number of affected sibs/half-sibs, not including self
	if (ind->pfa)
	{
		set<Individual*> aff_sib_ptr;
		for (each_element(ind->pfa->nextgen,iti)) if (is_affected(iti->second)) aff_sib_ptr.insert(iti->second);
		for (each_element(ind->pmo->nextgen,iti)) if (is_affected(iti->second)) aff_sib_ptr.insert(iti->second);
		aff_sib_ptr.erase(ind);
		aff_sib = aff_sib_ptr.size();
	}
	
	if		(unaff_sib_is_inf && aff_sib) ind->set(iData::INF,true);
	else if	(inf_child >1 || is_affected(ind) || ind->get(iData::GTN)>0) ind->set(iData::INF,true);
	else if (inf_child==0) ind->set(iData::INF,false);
	else if (inf_child==1)
	{
		ind->set(iData::INF,false);
		if (boost::indeterminate(spouse->get(iData::INF)))	_is_informative(spouse);
		if (spouse->get(iData::INF)==true)
			ind->set(iData::INF,true);
		else
		{
			if (!ind->pfa)
			{
				ind->set(iData::INF,false); // already false, but put it here to make things clear
			}
			else
			{
				ind->set(iData::INF,true);
				if (boost::indeterminate(ind->pfa->get(iData::INF))) _is_informative(ind->pfa);
				if (boost::indeterminate(ind->pmo->get(iData::INF))) _is_informative(ind->pmo);
				if (ind->pfa->get(iData::INF)==false && ind->pmo->get(iData::INF)==false)
					ind->set(iData::INF,false);
				else
				{
					ind->set(iData::INF,true);
					spouse->set(iData::INF,true);
				}
			}
		}
	}

	//set spouse=true if ind=true; otherwise, there will be a problem when a loop exist in family, such as celiac 4742
	if (ind->get(iData::INF)==true)
	{
		for (each_element(ind->fam->members,iti))
		{
			Individual * i = iti->second;
			if (exist_element(ind->nextgen,i->get(iData::IID)))
			{
				if (i->get(iData::INF)==true)
				{
					if		(ind->get(iData::SEX)==kMAN) i->pmo->set(iData::INF,true);
					else if (ind->get(iData::SEX)==kWMN) i->pfa->set(iData::INF,true);
					else exit_error("\nERROR: Unknown gender of individual "+ind->get(iData::IID));
				}
			}
		}
	}
}

void FamilyDB::cal_informativeness()
{
	lns<<showl<<"Search for uninformative Ind... "<<flush_logger;
	set_inds_status(iData::INF, boost::logic::tribool(boost::indeterminate));
	for (each_element(family_db,itf))
		for (each_element(itf->second->founders,iti))
			_is_informative(iti->second);
}

// ----------------- old version of cal_rel -----------------------

void FamilyDB::find_ancestors( map<string,vector<int> > * alist,
							   vector<int> * v,
							   Individual* ind,
							   int g)
{
	if (ind==NULL) return;
	string ind_name=ind->get(iData::IID);
	if (v!=NULL) (*alist)[ind_name]=*v;
	else (*alist)[ind_name];
	(*alist)[ind_name].push_back(g);
	if (ind->pfa) find_ancestors(alist,&(*alist)[ind_name],ind->pfa,1);
	if (ind->pmo) find_ancestors(alist,&(*alist)[ind_name],ind->pmo,2);
}

void FamilyDB::find_common_ancestors(Family* f,string prob,string subj)
{
	bool needprob=true;
	bool needsubj=true;
	if (strr.famname==f->pid)
	{
		if (strr.probname==prob) needprob=false;
		if (strr.subjname==subj) needsubj=false;
	}
	else
		strr.famname=f->pid;
	
	if (needprob)
	{
		strr.probname=prob;
		strr.a_prob.clear();
		find_ancestors(&(strr.a_prob),NULL,(f->members[prob]),0);
	}
	if (needsubj)
	{	
		strr.subjname=subj;
		strr.a_subj.clear();
		find_ancestors(&(strr.a_subj),NULL,(f->members[subj]),0);
	}
	
	if (needsubj || needprob)
	{
		strr.p_m='U';
		strr.min_g=9999;
		strr.g_prb=9999;
		strr.g_sbj=9999;
		strr.rcode=9999;
		strr.comanc.clear();
		strr.relation.clear();
		for (each_element(strr.a_prob,ita))
		{
			if (exist_element(strr.a_subj,ita->first))
			{
				int gp=ita->second.size()-1;
				int gs=strr.a_subj[ita->first].size()-1;
				int gt=gp+gs;
				strr.comanc[ita->first]=make_pair(gs,gp);
				if (gt<strr.min_g) { strr.min_g=gt; strr.g_prb=gp; strr.g_sbj=gs; }
			}
		}
	}
}

// update strr.comanc and strr.min_g
void FamilyDB::find_nearest_common_ancestors(Family* f,string prob,string subj)
{
	find_common_ancestors(f,prob,subj);
dia:
	for (each_element(strr.comanc,itca))
	{
		if ( (itca->second.first+itca->second.second)>strr.min_g)
		{
			strr.comanc.erase(itca);
			goto dia;
		}
	}
}

void FamilyDB::cal_relationship(Family* f,string prob,string subj)
{
	map<string,pair<int,int> >::iterator itca1,bkup;
	map<string,vector<int> >::iterator ita; 
	int x,y,absy,i;
	char temp[32];
	
	find_nearest_common_ancestors(f,prob,subj);
	
	if (strr.comanc.size()==0)
	{
		strr.rcode=9;
		strr.p_m='N';
		if (exist_element(f->members[prob]->spouses,subj))
			strr.relation="Spouse.";
		else
			strr.relation="Unrelated.";
		return;
	}
	if (strr.comanc.size()>2)
	{
		lns<<showl<<"fam "<<f->pid<<" ind "<<prob<<"-"<<subj<<" has more than 2 nearest common ancestors."<<flush_logger;
		strr.relation="Merge:";
	}
	if (strr.comanc.size()==2)
	{
		itca1=strr.comanc.begin();
		bkup=itca1;
		++bkup;
		if (!exist_element(f->members[itca1->first]->spouses,bkup->first))
		{
			//lns.message<<"fam "<<f->name<<" ind "<<prob<<"-"<<subj<<" has 2 nearest common ancestors without marriage.";
			//lns.show(logging::msg);
			strr.relation="Multi:";
		}
	}
	
	//set p_m
	bool pl=false,ml=false; //paternal, maternal
	for (each_element(strr.comanc,itca))
	{
		ita=strr.a_prob.find(itca->first);
		if (ita->second.size()>1)
		{
			if (ita->second[1]==1) pl=true;
			if (ita->second[1]==2) ml=true;
		}
	}
	if (pl)
	{	if (ml) strr.p_m='B';
	else	strr.p_m='P';
	}
	else
	{	if (ml) strr.p_m='M';
	else	strr.p_m='-';
	}
	
	//set relation
	for (each_element(strr.comanc,itca))
	{
		x=min ( itca->second.first , itca->second.second );
		y= itca->second.second - itca->second.first;
		absy=abs(y);
		if (x<0) exit_error("\nERROR: fam "+f->pid+" ind "+prob+"-"+subj+" x<0 !!!");
		//set rcode (careful, calculate according to the last read one)
		if (x>6 || absy>6) 
			strr.rcode=8;
		else
			strr.rcode=(x*10+absy)*(y<0?-1:1);
		if (strr.comanc.size()==1 && x!=0)
		{
			strr.rcode = strr.rcode<0 ? strr.rcode-200 : strr.rcode+200 ;
			strr.relation="Half";
		}
		switch(x)
		{
			case 0:
				if (y==0)
				{
					strr.relation+="Proband.";
					break;
				}
				if (absy>2)
					for (i=2;i<absy;++i)
						strr.relation+="Great";
				if (absy>1)
					strr.relation+="Grand";
				if (y>0)
				{
					switch (int(f->members[subj]->get(iData::SEX)))
					{
						case kMAN: strr.relation+="Father."; break;
						case kWMN: strr.relation+="Mother."; break;
						case kUNK: strr.relation+="Parent."; break;
						default:strr.relation+="Parent(sex_err)."; break;
					}
				}
				else //(y<0)
				{
					switch (int(f->members[subj]->get(iData::SEX)))
					{
						case kMAN: strr.relation+="Son."; break;
						case kWMN: strr.relation+="Daughter."; break;
						case kUNK: strr.relation+="Child."; break;
						default:strr.relation+="Child(sex_err)."; break;
					}
				}
				break;
			case 1:
				if (absy>1)
					for (i=1;i<absy;++i)
						strr.relation+="Great";
				if (y==0)
				{
					switch (int(f->members[subj]->get(iData::SEX)))
					{
						case kMAN: strr.relation+="Brother."; break;
						case kWMN: strr.relation+="Sister."; break;
						case kUNK: strr.relation+="Sibling."; break;
						default:strr.relation+="Sibling(sex_err)."; break;
					}
					break;
				}
				if (y>0)
				{
					switch (int(f->members[subj]->get(iData::SEX)))
					{
						case kMAN: strr.relation+="Uncle."; break;
						case kWMN: strr.relation+="Aunt."; break;
						case kUNK: strr.relation+="Uncle/Aunt."; break;
						default:strr.relation+="Uncle/Aunt(sex_err)."; break;
					}
				}
				else //(y<0)
				{
					switch (int(f->members[subj]->get(iData::SEX)))
					{
						case kMAN: strr.relation+="Nephew."; break;
						case kWMN: strr.relation+="Niece."; break;
						case kUNK: strr.relation+="Nephew/Niece."; break;
						default:strr.relation+="Nephew/Niece(sex_err)."; break;
					}
				}
				break;
			default:
				switch (x-1)
			{
				case 1: strr.relation+="First"; break;
				case 2: strr.relation+="Second";break;
				case 3: strr.relation+="Third"; break;
				default:strr.relation+=itoa(x-1,temp,10);
					strr.relation+="th";	break;
			};
				strr.relation+="Cousin";
				switch (absy)
			{
				case 0: strr.relation+="."; break;
				case 1: strr.relation+="OnceRemoved."; break;
				case 2: strr.relation+="TwiceRemoved."; break;
				default:strr.relation+=itoa(absy,temp,10);
					strr.relation+="TimesRemoved."; break;
			};
		};
		bkup=itca;
		++bkup;
		if (bkup!=strr.comanc.end())
			if (exist_element(f->members[itca->first]->spouses,bkup->first))
				itca=bkup;
	}
}

void FamilyDB::cal_relationship()
{
	lns<<showl<<"Calculate relationship between the probands and everyone else ... "<<flush_logger;
	lns.sub();
	for (each_element(family_db,itf))
	{
		Family* f = itf->second;
		IndMap probands;
		for (each_element_const(f->members,iti))
			if (is_proband(iti->second)) probands[iti->first]=iti->second;
		if (probands.size()<1) { lns<<showe<<"family "<<f->pid<<" has no probands, skipped."<<flush_logger; continue; }
		if (probands.size()>1) { lns<<showe<<"family "<<f->pid<<" has >1 probands, skipped."<<flush_logger; continue; }
		string prbname=probands.begin()->second->get(iData::IID);
		for (each_element(f->members,iti))
		{
			Individual* ind=iti->second;
			cal_relationship(f,prbname,ind->get(iData::IID));
			ind->rrs[prbname].code=strr.rcode;
			ind->rrs[prbname].relation=strr.relation;
			ind->rrs[prbname].p_m=strr.p_m;
		}
	}
	lns.endsub();
}
/*
void FamilyDB::cal_relationship_among_aff()
{
	lns<<showl<<"calculate relationship among affecteds ... ";
	lns.sub();
	program.outf << "FID"<<DLMTR<<"IID1"<<DLMTR<<"IID2"<<DLMTR<<"miosis"<<DLMTR<<"r_code"<<DLMTR<<"r_str\n";
	for (each_element(family_db,itf))
	{
		Family* f = itf->second;
		IndMap persons;
		for (each_element(f->members,iti))
			if (is_affected(iti->second)) persons[iti->first]=iti->second;
		for (each_element(persons,itp))
		{
			string prbname = itp->second->get(iData::IID);
			IndMap::iterator its=itp;
			for (++its; its!=persons.end(); ++its)
			{
				string sbjname = its->second->get(iData::IID);
				cal_relationship(f,prbname,sbjname);
				program.outf<< f->pid << DLMTR << prbname << DLMTR << sbjname << DLMTR
							<< strr.min_g << DLMTR << strr.rcode << DLMTR << strr.relation << endl;
			}
		}
	}
	lns.endsub();
}*/

void FamilyDB::read_IND(const string& filename, map<Family*, vector<Individual*> >& data)
{
	lns<<showl<<"Read individuals from "<<filename<<" ... "<<flush_logger;
	lns.sub();
	int nb_dropped = 0; // number of individuals dropped
	int input_type = 0; // 0=unknown 1=uid 2=pid+iid
	for (Rows_in_File(in, filename, 1))
	{
		if (in.RowNumber()==0)
		{
			if		(in[0]=="UID") input_type=1;
			else if (in.NumFields()==2 && in[0]=="PID" && in[1]=="IID") input_type=2;
			else exit_error("The first line of "+filename+" must be UID or PID IID");
			continue;
		}
		if (input_type==1)
		{
			if (exist_element(uniqid_db, in[0]))
			{
				Individual* ind=uniqid_db[in[0]];
				data[ind->fam].push_back(ind);
			}
			else ++nb_dropped;
		}
		else // must be 2
		{
			if (in.NumFields()<2) exit_error("Some lines in "+filename+" have <2 columns.");
			Individual* ind = find_ind(in[0], in[1]);
			if (ind!=NULL) data[ind->fam].push_back(ind);
			else ++nb_dropped;
		}
	}
	if (nb_dropped)
		lns<<showe<<nb_dropped<<" individuals not recognized."<<flush_logger;
	lns.endsub();
}

void FamilyDB::cal_MRCA_file(const string& filename)
{
	map<Family*, vector<Individual*> > candidates;
	read_IND(filename, candidates);
	
	lns<<showl<<"Find the most recent common ancestors ... "<<flush_logger;
	lns.sub();	
	program.outf<<"UID_1"<<DLMTR<<"UID_2"<<DLMTR<<"gen_1"<<DLMTR<<"gen_2"<<DLMTR<<"NCAs\n";
	for (each_element(candidates, itf))
	{
		for (each_element(itf->second, it1))
		{
			vector<Individual*>::iterator it2=it1;
			for (++it2; it2!=itf->second.end(); ++it2)
			{
				find_nearest_common_ancestors((*it2)->fam, (*it1)->get(iData::IID), (*it2)->get(iData::IID));
				if (!strr.comanc.empty())
				{
					program.outf << (*it1)->get(iData::UID) << DLMTR << (*it2)->get(iData::UID) << DLMTR ;
					program.outf << strr.g_prb << DLMTR << strr.g_sbj << DLMTR ;
					print_map_first(strr.comanc, program.outf, ' ', true);
				}
			}
		}
	}
	lns.endsub();
}

// find connecting individuals from individual ind to ancestors anc.
// anc themselves are not connectors, which is good for finding obligatory carriers (don't know which one of the anc is carrier)
// but ind itself is a connectors, also good for finding obligatory carriers (it's a known carrier)
void FamilyDB::find_connections(Individual * ind,
								set<string>& tca,
								set<string>& connectors)
{
	if (ind==NULL) return;
	string ind_name=ind->get(iData::IID);
	if (exist_element(tca,ind_name)) return;
	FamilyDB::ID_Count ancestors;
	_find_ancestors(ind,ancestors);
	bool not_connector=false;
	for (auto &x:tca) if (!exist_element(ancestors,x)) { not_connector=true; break; }
	if (not_connector) return;
	connectors.insert(ind_name);
	if (ind->pfa) find_connections(ind->pfa,tca,connectors);
	if (ind->pmo) find_connections(ind->pmo,tca,connectors);
}

void FamilyDB::cal_OC_file(const string& filename)
{
	cal_generation();
	map<Family*, vector<Individual*> > candidates;
	read_IND(filename, candidates);
	lns<<showl<<"Find obligatory carriers of mutations ... "<<flush_logger;
	lns.sub();
	program.outf<<"PID"<<DLMTR<<"IID"<<"\n";
	for (each_element(candidates, itf))
	{
		if (itf->second.size()<2)
		{
			if (to_show_warnings) lns<<showw<<"Fam "<<itf->first->pid<<" has less than 2 known carriers, skipped."<<flush_logger;
			continue;
		}
		set<string> common_anc;
		bool first_run=true;
		for (each_element(itf->second, it1))
		{
			FamilyDB::ID_Count ancestors;
			_find_ancestors(*it1,ancestors);
			if (first_run)
			{
				for (auto &x:ancestors) common_anc.insert(x.first);
				first_run=false;
			}
			else
			{
				set<string> tmp;
				for (auto &x:ancestors) if (exist_element(common_anc,x.first)) tmp.insert(x.first);
				common_anc=tmp;
			}
		}
		if (common_anc.empty())
		{
			if (to_show_warnings) lns<<showw<<"Fam "<<itf->first->pid<<" has no common ancestors for the specified individuals."<<flush_logger;
			continue;
		}
		if (common_anc.size()>2)
		{
			int max_gen=-1;
			for (auto &x:common_anc)
			{
				int this_gen = itf->first->members[x]->get(iData::GEN);
				if (this_gen>max_gen) max_gen=this_gen;
			}
			set<string> tmp;
			for (auto &x:common_anc)
			{
				int this_gen = itf->first->members[x]->get(iData::GEN);
				if (this_gen==max_gen) tmp.insert(x);
			}
			common_anc=tmp;
			if (common_anc.size()>2)
			{
				if (to_show_warnings) lns<<showw<<"Fam "<<itf->first->pid<<" has >2 common ancestors for the specified individuals."<<flush_logger;
				continue;
			}
		}
		if (common_anc.size()==2)
		{
			Individual* i0 = itf->first->members[*(common_anc.begin())];
			if (!exist_element(i0->spouses,*(common_anc.rbegin())))
			{
				if (to_show_warnings) lns<<showw<<"Fam "<<itf->first->pid<<" has  2 common ancestors for the specified individuals, but they are not couples."<<flush_logger;
				continue;
			}
		}
		set<string> all_carriers;
		set<string> knn_carriers;
		for (each_element(itf->second, it1))
		{
			knn_carriers.insert((*it1)->get(iData::IID));
			set<string> connectors;
			find_connections(*it1,common_anc,connectors);
			for (auto &x:connectors) all_carriers.insert(x);
		}
		for (auto &x:all_carriers)
			if (!exist_element(knn_carriers,x)) program.outf << itf->first->pid << DLMTR << x << endl;
		if (common_anc.size()==1)
			program.outf << itf->first->pid << DLMTR << *(common_anc.begin()) << endl;
	}
	lns.endsub();
}

void _cal_kinship_up(Individual* ind, map<Individual*,int>& calculated, vector< vector<double> >& kinship)
{
	if (ind->get(iData::VST)==true) return;
	if (ind->pfa) _cal_kinship_up(ind->pfa,calculated,kinship);
	if (ind->pmo) _cal_kinship_up(ind->pmo,calculated,kinship);
	
	int co_fa = calculated[ind->pfa];
	int co_mo = calculated[ind->pmo];
	int co_me = calculated.size();
	kinship[co_me][co_me] = (1+kinship[co_mo][co_fa])*.5;
	ind->set(iData::INB, kinship[co_mo][co_fa]);
	for (each_element(calculated,itc))
	{
		int co_rf = itc->second;
		double kc_ij = (kinship[co_rf][co_fa]+kinship[co_rf][co_mo])*.5;
		kinship[co_rf][co_me] = kc_ij;
		kinship[co_me][co_rf] = kc_ij;
	}
	calculated[ind]=co_me;
	ind->set(iData::VST,true);
}

void FamilyDB::cal_kinship_all(bool quietly)
{
	lns<<showl<<"Calculate kinship and inbreeding coefficient of all INDs ... "<<flush_logger;
	set_inds_status(iData::VST, false);
	if (!quietly) program.outf << "FID"<<DLMTR<<"IID1"<<DLMTR<<"IID2"<<DLMTR<<"kinship\n";
	for (each_element(family_db,itf))
	{
		// data
		Family* f = itf->second;
		map<Individual*,int> calculated; // IID -> coordinate_in_kinship_matrx
		vector< vector<double> > kinship( f->members.size(), vector<double>(f->members.size(),-1) );
		
		// calculate kinship
		for (each_element(f->founders,iti))
		{
			int co=calculated.size();
			calculated[iti->second]=co;
			kinship[co][co]=0.5;
			iti->second->set(iData::INB, 0);
			iti->second->set(iData::VST,true);
		}
		for (each_element(f->founders,it1))
		{
			int co1=calculated[it1->second];
			IndMap::iterator it2=it1;
			for (++it2; it2!=f->founders.end(); ++it2)
			{
				int co2=calculated[it2->second];
				kinship[co1][co2]=0;
				kinship[co2][co1]=0;
			}
		}
		for (each_element(f->nonfounders, itnc))
			_cal_kinship_up(itnc->second, calculated, kinship);

		// print results
		if (!quietly)
		{
			for (each_element(calculated,it1))
			{
				int co_me = it1->second;
				program.outf<< f->pid << DLMTR << it1->first->get(iData::IID) << DLMTR << it1->first->get(iData::IID) << DLMTR << kinship[co_me][co_me] << endl;
				map<Individual*,int>::iterator it2=it1;
				for (++it2; it2!=calculated.end(); ++it2)
				{
					int co_rf = it2->second;
					program.outf<< f->pid << DLMTR << it1->first->get(iData::IID) << DLMTR << it2->first->get(iData::IID) << DLMTR << kinship[co_me][co_rf] << endl;
				}
			}
		}
	}
	set_inds_status(iData::VST, false);
}

void FamilyDB::cal_kinship_file(const string& filename, bool quietly)
{
	map<Family*, vector<Individual*> > candidates;
	read_IND(filename, candidates);
	
	lns<<showl<<"Calculate kinship and inbreeding coefficient ... "<<flush_logger;
	set_inds_status(iData::VST, false);
	if (!quietly) program.outf << "FID"<<DLMTR<<"IID1"<<DLMTR<<"IID2"<<DLMTR<<"kinship\n";
	for (each_element(candidates,itf))
	{
		// data
		Family* f = itf->first;
		map<Individual*,int> calculated; // IID -> coordinate_in_kinship_matrx
		vector< vector<double> > kinship( f->members.size(), vector<double>(f->members.size(),-1) );
		
		// calculate kinship
		for (each_element(f->founders,iti))
		{
			int co=calculated.size();
			calculated[iti->second]=co;
			kinship[co][co]=0.5;
			iti->second->set(iData::INB, 0);
			iti->second->set(iData::VST,true);
		}
		for (each_element(f->founders,it1))
		{
			int co1=calculated[it1->second];
			IndMap::iterator it2=it1;
			for (++it2; it2!=f->founders.end(); ++it2)
			{
				int co2=calculated[it2->second];
				kinship[co1][co2]=0;
				kinship[co2][co1]=0;
			}
		}
		for (each_element(itf->second, itnc))
			_cal_kinship_up(*itnc, calculated, kinship);
		
		// print results
		if (!quietly)
		{
			for (each_element(itf->second,it1))
			{
				Individual* ind1 = *it1;
				int co_me = calculated[ind1];
				program.outf<< f->pid << DLMTR << ind1->get(iData::IID) << DLMTR << ind1->get(iData::IID) << DLMTR << kinship[co_me][co_me] << endl;
				vector<Individual*>::iterator it2=it1;
				for (++it2; it2!=itf->second.end(); ++it2)
				{
					Individual* ind2 = *it2;
					int co_rf = calculated[ind2];
					program.outf<< f->pid << DLMTR << ind1->get(iData::IID) << DLMTR << ind2->get(iData::IID) << DLMTR << kinship[co_me][co_rf] << endl;
				}
			}
		}
	}
	set_inds_status(iData::VST, false);
}

void FamilyDB::find_founders_of(const string& filename)
{
	map<Family*, vector<Individual*> > candidates;
	read_IND(filename, candidates);
	
	lns<<showl<<"find founders ... "<<flush_logger;
	for (each_element(candidates,itf))
	{
		for (each_element(itf->second, iti))
		{
			program.outf << (*iti)->get(iData::UID) ;
			map<Individual*,int> founders;
			_find_founders(*iti,founders);
			for (each_element(founders,ita)) program.outf << DLMTR << ita->first->get(iData::UID);
			program.outf << endl;
		}
	}
}

// allele_type=2 means risk allele, coded as negative integers. Otherwise non-risk allele, coded as positive integer.
// ind's a1 and a2 should already have genotypes (non-zero integers), and this function doesn't check whether this is the case.
int _gene_drop_get(Individual* ind, int allele_type)
{
	int a[2];
	a[0]=ind->get(iData::AL1);
	a[1]=ind->get(iData::AL2);
	if (allele_type==2)
	{
		if (a[0]<0 && a[1]<0)	return a[MyRandom(2)];
		else if (a[0]<0)		return a[0];
		else if (a[1]<0)		return a[1];
		else					exit_error("wrong genotype.");
	}
	else
	{
		if (a[0]<0 && a[1]<0)	exit_error("wrong genotype.");
		else if (a[0]<0)		return a[1];
		else if (a[1]<0)		return a[0];
		else					return a[MyRandom(2)];
	}
	exit_error("this should not happen");
	return 0;
}

// for travers from bottom up. Assume founders are populated with genotypes (non-zero integers, positive=non-risk, negative=risk allele).
// The current ind does not have genotypes yet but has allele types.
void _gene_drop(Individual* ind)
{
	if (ind->get(iData::VST)==true) return;
	if (ind->pfa) _gene_drop(ind->pfa);
	if (ind->pmo) _gene_drop(ind->pmo);
	if		(ind->get(iData::AL1)==2 && ind->get(iData::AL2)==2)
	{
		ind->set(iData::AL1,_gene_drop_get(ind->pfa,2));
		ind->set(iData::AL2,_gene_drop_get(ind->pmo,2));
	}
	else if (ind->get(iData::AL1)==2 || ind->get(iData::AL2)==2)
	{
		int fa_a[2],mo_a[2];
		fa_a[0]=ind->pfa->get(iData::AL1);
		fa_a[1]=ind->pfa->get(iData::AL2);
		mo_a[0]=ind->pmo->get(iData::AL1);
		mo_a[1]=ind->pmo->get(iData::AL2);
		bool fa_has2 = (fa_a[0]<0 || fa_a[1]<0);
		bool mo_has2 = (mo_a[0]<0 || mo_a[1]<0);
		bool fa_has1 = (fa_a[0]>0 || fa_a[1]>0);
		bool mo_has1 = (mo_a[0]>0 || mo_a[1]>0);
		if (fa_has2 && mo_has2)
		{
			if (fa_has1 && mo_has1)
			{
				if (MyRandom(2))	{ ind->set(iData::AL1,_gene_drop_get(ind->pfa,2)); ind->set(iData::AL2,_gene_drop_get(ind->pmo,1)); }
				else				{ ind->set(iData::AL1,_gene_drop_get(ind->pfa,1)); ind->set(iData::AL2,_gene_drop_get(ind->pmo,2)); }
			}
			else if (fa_has1)		{ ind->set(iData::AL1,_gene_drop_get(ind->pfa,1)); ind->set(iData::AL2,_gene_drop_get(ind->pmo,2)); }
			else if (mo_has1)		{ ind->set(iData::AL1,_gene_drop_get(ind->pfa,2)); ind->set(iData::AL2,_gene_drop_get(ind->pmo,1)); }
			else exit_error("wrong genotype.");
		}
		else if (fa_has2)			{ ind->set(iData::AL1,_gene_drop_get(ind->pfa,2)); ind->set(iData::AL2,_gene_drop_get(ind->pmo,1)); }
		else if (mo_has2)			{ ind->set(iData::AL1,_gene_drop_get(ind->pfa,1)); ind->set(iData::AL2,_gene_drop_get(ind->pmo,2)); }
		else exit_error("wrong genotype.");
	}
	else
	{
		ind->set(iData::AL1,_gene_drop_get(ind->pfa,1));
		ind->set(iData::AL2,_gene_drop_get(ind->pmo,1));
	}
	ind->set(iData::VST,true);
}

// Assume all people have alleles types but not genotypes. This function populates all members with genotypes based on allele types.
// allele_type=2 means risk allele, coded as negative integers. Otherwise non-risk allele, coded as positive integer.
// This scheme works if 1) no alleles have been read, then it generate genotypes not co-segregating with disease
//                      2) alleles were read, which were genereated by SLINK, then it generate genotypes co-segregating with the disease.
void FamilyDB::gene_drop(bool quietly)
{
	lns<<showl<<"do gene drop ... "<<flush_logger;
	set_inds_status(iData::VST, false);
	for (each_element(family_db,itf))
	{
		int maxpos=0;
		int minneg=0;
		for (each_element(itf->second->founders,iti))
		{
			if (iti->second->get(iData::AL1)==2) iti->second->set(iData::AL1,--minneg); else iti->second->set(iData::AL1,++maxpos);
			if (iti->second->get(iData::AL2)==2) iti->second->set(iData::AL2,--minneg); else iti->second->set(iData::AL2,++maxpos);
			iti->second->set(iData::VST,true);
		}
		for (each_element(itf->second->nonfounders,iti)) _gene_drop(iti->second);
	}
	set_inds_status(iData::VST, false);
}

// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2795894/
void FamilyDB::ind_weight()
{
	lns<<showl<<"calc ind wt ... "<<flush_logger;
	const int num_rounds=1000;
	set_inds_values(iData::IWT, 0);
	for (int r=0;r<num_rounds;++r)
	{
		set_inds_values(iData::AL1, 0);
		set_inds_values(iData::AL2, 0);
		gene_drop(true);
		for (each_element(family_db,itf))
		{
			map<int,int> allele_count;
			for (each_element(itf->second->members,iti))
			{
				Individual * ind = iti->second;
				if (!ind->get(iData::GTN)) continue;
				if (iWt_useCasesOnly && !is_affected(ind)) continue;
				int a1 = ind->get(iData::AL1); ++allele_count[a1];
				int a2 = ind->get(iData::AL2); ++allele_count[a2];
			}
			if (allele_count.empty()) continue; // exit_error("No member in family "+itf->second->pid+" has genotype. Please note that --ind-wt requires GTN.");
			for (each_element(itf->second->members,iti))
			{
				Individual * ind = iti->second;
				if (!ind->get(iData::GTN)) continue;
				if (iWt_useCasesOnly && !is_affected(ind)) continue;
				double this_wt = ind->get(iData::IWT);
				int a1 = ind->get(iData::AL1); this_wt+=0.5/allele_count[a1];
				int a2 = ind->get(iData::AL2); this_wt+=0.5/allele_count[a2];
				ind->set(iData::IWT,this_wt);
			}
		}
	}
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->members,iti))
		{
			Individual * ind = iti->second;
			ind->set(iData::IWT,ind->get(iData::IWT) / num_rounds);
		}
	}
	if (iWt_SumTo1_eaPed)
	{
		for (each_element(family_db,itf))
		{
			double sum_wt = 0;
			for (each_element(itf->second->members,iti))
			{
				Individual * ind = iti->second;
				sum_wt += ind->get(iData::IWT);
			}
			if (sum_wt==0) continue;
			for (each_element(itf->second->members,iti))
			{
				Individual * ind = iti->second;
				ind->set(iData::IWT,ind->get(iData::IWT) / sum_wt);
			}
		}
	}
	set_inds_values(iData::AL1, 0);
	set_inds_values(iData::AL2, 0);
}

void FamilyDB::single_aff_par()
{
	lns<<showl<<"Get discordant parents with affected children ... "<<flush_logger;
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->members,iti))
		{
			Individual * ind = iti->second;
			if (!is_affected(ind)) continue;
			if (ind->nextgen.empty()) continue;
			for (each_element(ind->spouses,its))
			{
				Individual * spo = its->second;
				if (is_affected(spo)) continue;
				int affected_children = 0;
				for (each_element(ind->nextgen,itc))
				{
					Individual * chd = itc->second;
					if (!is_affected(chd)) continue;
					if (ind->get(iData::SEX)==kMAN && chd->pfa==ind && chd->pmo==spo) ++affected_children;
					if (ind->get(iData::SEX)==kWMN && chd->pfa==spo && chd->pmo==ind) ++affected_children;
				}
				if (affected_children) ind->set(iData::B01,true);
			}
		}
	}
}

void FamilyDB::cal_aff_fdr()
{
	lns<<showl<<"Calculate the number of affected FDRs ... "<<flush_logger;
	for (each_element(family_db,itf))
	{
		for (each_element(itf->second->members,iti))
		{
			Individual * ind = iti->second;
			if (!is_affected(ind)) { ind->set(iData::AFD,std::numeric_limits<double>::signaling_NaN()); continue; }
			double affected_FDR = 0;
			if (ind->pfa!=NULL)
			{
				if (is_affected(ind->pfa)) ++affected_FDR;
				if (is_affected(ind->pmo)) ++affected_FDR;
				for (each_element(ind->pfa->nextgen,itc))
				{
					Individual * chd = itc->second;
					if (chd!=ind && is_affected(chd) && chd->pfa==ind->pfa && chd->pmo==ind->pmo) ++affected_FDR;
				}
			}
			for (each_element(ind->nextgen,itc))
			{
				Individual * chd = itc->second;
				if (is_affected(chd)) ++affected_FDR;
			}
			ind->set(iData::AFD,affected_FDR);
		}
	}
}

