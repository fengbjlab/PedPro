#include <boost/assign/list_of.hpp> // for list_of
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_log.hpp>
#include "libfbj_famdb_ind.hpp"

extern logger lns;

namespace iData
{
	using namespace std;
	
	std::vector<string> v_ID = boost::assign::list_of\
	("IID")("UID")("POP")("PID")("FTH")("MTH")("REL")("L5D")("L5N")("L6D")("L6N")("L7D")("L7N")("L8D")("L8N")("L9D")("L9N")("DNP")("DNI")("DNU")("DNF")("DNM")("MZT")("GTP")("ALT")("CMT")\
	("SEX")("AFF")("LIA")("PRB")("AGE")("YOB")("GTN")("GTE")("CLT")("INB")("GEN")("MXG")("DES")("AL1")("AL2")("AFD")("IWT")\
	("INF")("SEP")("VST")("B01")("RCG")("OTH")("DUM")("SXD")("YBI")("ABP");

	namespace
	{
		map<int, set<string> >			vInText;	// var name for input, although defaults length <=4, I allow longer names
		map<int, map<string,string> >	vTextCvt;	// var text conversion
		map<int, set<string> >			vTextIgn;	// var text to ignore
	}
	
	map<int, string>		vOuText;		// var name for output, default length <=7
	const std::string default_vname = "\
# Format of this file:Use SPACE/TAB as delimiters; multiple delimiters are treated as one.\n\
# Lines begin with # are comments. PEDPRO stops reading this file at the first blank line.\n\
# ========================================================================================\n\
# Description             Type   ID     Symbol  Synonyms(<=4 characters, case-insensitive)\n\
# ========================================================================================\n\
Pedigree_ID               STR    PID    PedID   pedi pid p_id ped fami fid f_id fam kind\n\
Individual_ID             STR    IID    IndID   indi iid i_id ind id subj pers\n\
Unique_ID                 STR    UID    UniqID  uniq uid u_id\n\
Father_ID                 STR    FTH    Father  fath fth dad fa pa PaID\n\
Mother_ID                 STR    MTH    Mother  moth mth mom mo ma MaID\n\
Downcoded_PID             STR    DNP    dPID    dpid dp\n\
Downcoded_IID             STR    DNI    dIID    diid di\n\
Downcoded_UID             STR    DNU    dUID    duid du\n\
Downcoded_FTH             STR    DNF    dFTH    dfth dfa ddad\n\
Downcoded_MTH             STR    DNM    dMTH    dmth dmo dmom\n\
Population                STR    POP    Pop     popu pop\n\
Monozygosity_Twin         STR    MZT    MZ_twin mz_t mztw mzt mz twin\n\
Genotype                  STR    GTP    Geno    geno mutn\n\
Alternative_ID            STR    ALT    AltID   alti alt\n\
Comment                   STR    CMT    Comment comm cmt\n\
Sex                       DBL    SEX    Gender  gend sex sx sex_\n\
Affection_Status          DBL    AFF    Aff     affe aff af aff_\n\
Liability                 DBL    LIA    Liab    liab lia li lia_\n\
Proband                   DBL    PRB    Proband prob prb pr prb_\n\
Age                       DBL    AGE    Age     age ag\n\
Year_Of_Birth             DBL    YOB    YoB     yob byr birt\n\
Genotype_Numer            DBL    GTN    GTN     gtn\n\
Genotype_Error            DBL    GTE    GTE     gte\n\
Cluster_Number            DBL    CLT    Cluster clus cl\n\
Inbreeding_Coefficient    DBL    INB    Inbr    inbr inb\n\
Generation_number         DBL    GEN    Gen     gen g\n\
Descendants_MaxNo.Gen     DBL    MXG    GDes    gdes\n\
Descendants_TotalNo.      DBL    DES    NDes    ndes\n\
Allele_1                  DBL    AL1    AL1     al1 a1\n\
Allele_2                  DBL    AL2    AL2     al2 a2\n\
Individual_Weight         DBL    IWT    IndWt   indw wt\n\
# ----------------------------------------------------------------------------------------\n\
# The following should not appear in a header line of a pedigree file:\n\
Affected_FirstDegreeRel   DBL    AFD    AFD     afd\n\
Informative               TRB    INF    Info    info inf\n\
Separated                 TRB    SEP    Separat sepa sepr sep isol\n\
Visited                   TRB    VST    Visited visi vst\n\
TriBool_01(not_used)      TRB    B01    TB01    tb01 b01\n\
(person)_Created_Dummy    TRB    DUM    Dummy   dumm dum\n\
Sex_Deduced               TRB    SXD    Sx_Ded  sx_d sxd\n\
Year_Of_Birth_Inferred    TRB    YBI    ybi     ybi\n\
(person)_Add_By_Program   TRB    ABP    abp     abp\n\
LinkagePed_5Col_Downcoded STR    L5D    L5D     l5d\n\
LinkagePed_6Col_Downcoded STR    L6D    L6D     l6d\n\
LinkagePed_7Col_Downcoded STR    L7D    L7D     l7d\n\
LinkagePed_8Col_Downcoded STR    L8D    L8D     l8d\n\
LinkagePed_9Col_Downcoded STR    L9D    L9D     l9d\n\
LinkagePed_5Col_Normal    STR    L5N    L5      l5\n\
LinkagePed_6Col_Normal    STR    L6N    L6      l6\n\
LinkagePed_7Col_Normal    STR    L7N    L7      l7\n\
LinkagePed_8Col_Normal    STR    L8N    L8      l8\n\
LinkagePed_9Col_Normal    STR    L9N    L9      l9\n\
Relationship              STR    REL    Relat   rela rel\n\
Columns_recognized        TRB    RCG    RECG    recg rcg\n\
Columns_not_recognized    TRB    OTH    OTHR    othr oth othe\n\
# ========================================================================================\n\
\n\
!!! Don't remove the above blank line; the program needs it to stop reading this file.\n\
Variable_description can use any character except SPACE/TAB. Case sensitive.\n\
Variable_Type is one of STR/DBL/INT/CHR/TRB/BLN. Case sensitive. (TRB=3-val-boolean)\n\
Variable_ID must be 3 characters long. Case sensitive. Better all use upper case letters.\n\
Variable_Output is for the 1st line in output. Case sensitive. Must be <8 characters.\n\
Variable_Input is for the 1st line in input. Case insensitive. Must be <=4 characters.\n\
The abbreviated Output (lower case, 4 char) should be one of the Input names.\n\
_ID _Input _Output must contain _ or alphabets or digits only, and not begins w/ digits.\n\
INF -- informative for linkage analysis \n\
SEP -- unconnected with all other family members \n\
L5x -- PID,IID,FA,MO,SEX \n\
L6x -- PID,IID,FA,MO,SEX,AFF \n\
L7x -- PID,IID,FA,MO,SEX,AFF,LIA \n\
L8x -- PID,IID,FA,MO,SEX,AFF,AL1,AL2 \n\
L9x -- PID,IID,FA,MO,SEX,AFF,LIA,AL1,AL2 \n\
";

	void set_known_var(const string& input)
	{
		string instruction = expand_var_names(input);
		map<int, set<string> > vInText2;
		vInText2[L5N.second] = vInText[L5N.second];
		vInText2[L6N.second] = vInText[L6N.second];
		vInText2[L7N.second] = vInText[L7N.second];
		vInText2[L8N.second] = vInText[L8N.second];
		vInText2[L9N.second] = vInText[L9N.second];
		vInText2[L5D.second] = vInText[L5D.second];
		vInText2[L6D.second] = vInText[L6D.second];
		vInText2[L7D.second] = vInText[L7D.second];
		vInText2[L8D.second] = vInText[L8D.second];
		vInText2[L9D.second] = vInText[L9D.second];
		vInText2[DNP.second] = vInText[DNP.second];
		vInText2[DNI.second] = vInText[DNI.second];
		vInText2[DNU.second] = vInText[DNU.second];
		vInText2[DNF.second] = vInText[DNF.second];
		vInText2[DNM.second] = vInText[DNM.second];
		vInText2[REL.second] = vInText[REL.second];
		vInText2[RCG.second] = vInText[RCG.second];
		vInText2[OTH.second] = vInText[OTH.second];
		vInText2[CLT.second] = vInText[CLT.second];
		vInText2[INB.second] = vInText[INB.second];
		vInText2[INF.second] = vInText[INF.second];
		vInText2[SEP.second] = vInText[SEP.second];
		while (!instruction.empty())
		{
			string s=extract_name(instruction);
			if (s.empty()) exit_error("Can't extract name from "+instruction);
			if (!instruction.empty()) extract_char(instruction);
			int var_code = vCode(get_id_from(s));
			vInText2[var_code] = vInText[var_code];
		}
		vInText = vInText2;
	}
	
	void set_ignore_text(const degenerated_t& var_label, const string& input) {
		vTextIgn[var_label.second].insert(input);
	}
	
	bool to_ignore_text(const degenerated_t& var_label, const string& input) {
		return vTextIgn[var_label.second].find(input) != vTextIgn[var_label.second].end() ;
	}
	
	string convert(const degenerated_t& var_label, const string& input) {
		map<int, map<string,string> >::iterator it1 = vTextCvt.find(var_label.second);
		if (it1==vTextCvt.end()) return input;
		map<string,string>::iterator it2 = it1->second.find(input);
		if (it2==it1->second.end()) return input;
		return it2->second;
	}
		
	void read_conversion_map(const string& instruction) 
	{
		for (size_t start=0,found=0;!instruction.empty() && start<instruction.size();) 
		{
			if ((found=instruction.find(':', start))==string::npos) exit_error("instruction for variable content conversion -- no ':' found");
			string var_string = instruction.substr(start, found-start);
			start = found+1 ;
			degenerated_t var = get_id_from(var_string);
			for (;;) 
			{
				if ((found=instruction.find('=', start))==string::npos) exit_error("instruction for variable content conversion -- no '=' found");
				string ori_string = instruction.substr(start, found-start);
				start = found+1 ;
				bool to_break;
				size_t found_c=instruction.find(',', start);
				size_t found_s=instruction.find(';', start);
				if		(found_c<found_s) { found=found_c;				to_break=false; }
				else if (found_c>found_s) { found=found_s;				to_break=true;  }
				else					  { found=instruction.size();	to_break=true;  }
				string new_string = instruction.substr(start, found-start);
				start = found+1 ;
				
				// add data
				typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
				boost::char_separator<char> sep("/");
				tokenizer tokens(ori_string, sep);
				for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter!=tokens.end(); ++tok_iter)
					vTextCvt[var.second][*tok_iter]=new_string;
				if (to_break) break;
			}
		}
	}
	
	string expand_var_names(string instruction)
	{
		stringstream result;
		for (int i=0; !instruction.empty(); i++)
		{
			string s=extract_name(instruction);
			if (s.empty()) exit_error("Can't extract name from "+instruction);
			if (!instruction.empty()) extract_char(instruction);
			if (i) result<<',';
			degenerated_t field = get_id_from(s);
			if		(is_equal(field,L5N)) result<<vText(PID)<<','<<vText(IID)<<','<<vText(FTH)<<','<<vText(MTH)<<','<<vText(SEX);
			else if (is_equal(field,L6N)) result<<vText(PID)<<','<<vText(IID)<<','<<vText(FTH)<<','<<vText(MTH)<<','<<vText(SEX)<<','<<vText(AFF);
			else if (is_equal(field,L7N)) result<<vText(PID)<<','<<vText(IID)<<','<<vText(FTH)<<','<<vText(MTH)<<','<<vText(SEX)<<','<<vText(AFF)<<','<<vText(LIA);
			else if (is_equal(field,L8N)) result<<vText(PID)<<','<<vText(IID)<<','<<vText(FTH)<<','<<vText(MTH)<<','<<vText(SEX)<<','<<vText(AFF)<<','<<vText(AL1)<<','<<vText(AL2);
			else if (is_equal(field,L9N)) result<<vText(PID)<<','<<vText(IID)<<','<<vText(FTH)<<','<<vText(MTH)<<','<<vText(SEX)<<','<<vText(AFF)<<','<<vText(LIA)<<','<<vText(AL1)<<','<<vText(AL2);
			else if	(is_equal(field,L5D)) result<<vText(DNP)<<','<<vText(DNI)<<','<<vText(DNF)<<','<<vText(DNM)<<','<<vText(SEX);
			else if (is_equal(field,L6D)) result<<vText(DNP)<<','<<vText(DNI)<<','<<vText(DNF)<<','<<vText(DNM)<<','<<vText(SEX)<<','<<vText(AFF);
			else if (is_equal(field,L7D)) result<<vText(DNP)<<','<<vText(DNI)<<','<<vText(DNF)<<','<<vText(DNM)<<','<<vText(SEX)<<','<<vText(AFF)<<','<<vText(LIA);
			else if (is_equal(field,L8D)) result<<vText(DNP)<<','<<vText(DNI)<<','<<vText(DNF)<<','<<vText(DNM)<<','<<vText(SEX)<<','<<vText(AFF)<<','<<vText(AL1)<<','<<vText(AL2);
			else if (is_equal(field,L9D)) result<<vText(DNP)<<','<<vText(DNI)<<','<<vText(DNF)<<','<<vText(DNM)<<','<<vText(SEX)<<','<<vText(AFF)<<','<<vText(LIA)<<','<<vText(AL1)<<','<<vText(AL2);
			else if (field.second>=0) result<<vText(field);
			else result<<s;
		}
		return result.str();
	}

	// return whether it's an expected option: --XXX=yyy, where XXX is var_name, yyy is field_name
	bool get_var_name(string option) {
		if (option.size()<7) return false;
		if (!str_startsw(option, "--")) return false;
		if (option[5]!='=') return false;
		string var_name = to_upper_copy(option.substr(2,3));
		size_t var_id;
		for (var_id=0;var_id<v_ID.size();++var_id) if (v_ID[var_id]==var_name) break;
		if (var_id==v_ID.size()) return false;
		string in_name = to_lower_copy(option.substr(6,4));
		for (each_element(vInText, it)) it->second.erase(in_name);
		vInText[var_id].insert(in_name);
		vOuText[var_id]=option.substr(6);
		return true;
	}
	
	void _get_var_name_db(std::istream& input) {
		int error1=elog.get_token("lines in has an unknown ID.");
		set<string> short_names;
		set<string> long_names;
		set<string> identifiers;
		static const set<string> type_names = boost::assign::list_of("STR")("DBL")("INT")("CHR")("TRB")("BLN").to_container(type_names);
		for (each_line(input))
		{
			// read column 1-4
			vector<string> inf(4);
			input >> inf[0] >> inf[1] >> inf[2] >> inf[3];
			
			// check column 1
			if (long_names.find(inf[0])!=long_names.end())
				exit_error("Duplicated variable description : "+inf[0]);
			else
				long_names.insert(inf[0]);
			// check column 2
			if (type_names.find(inf[1])==type_names.end())
				exit_error("Unknown variable type : "+inf[1]);
			// check column 3
			if (identifiers.find(inf[2])!=identifiers.end())
				exit_error("Duplicated variable ID : "+inf[2]);
			else
				identifiers.insert(inf[2]);
			if (inf[2].size()!=3) exit_error("Variable ID must be 3 characters long.");
			if (!is_a_valid_name(inf[2])) exit_error("Variable ID "+inf[2]+" is not valid.");
			size_t var_id;
			for (var_id=0;var_id<v_ID.size();++var_id) if (v_ID[var_id]==inf[2]) break;
			if (var_id==v_ID.size()) { elog.add(error1); continue; }
			// check column 4
			if (inf[3].size()>=8)
				exit_error("Variable Output name must be <=7 characters long.");
			if (!is_a_valid_name(inf[3])) exit_error("Variable Output name "+inf[3]+" is not valid.");
			string abbrev_name=to_lower_copy(inf[3].substr(0,4));
			vOuText[var_id]=inf[3];
			// check column 5 and up
			bool abbr_found=false;
			for (skip_whitespaces(input); !is_blank_row(input); skip_whitespaces(input))
			{
				string sn; // short name
				input >> sn;
				boost::to_lower(sn);
				if (!is_a_valid_name(sn)) exit_error("Variable short name "+sn+" is not valid.");
				if (short_names.find(sn)!=short_names.end())
					exit_error("Duplicated variable short names : "+sn);
				else
					short_names.insert(sn);
				vInText[var_id].insert(sn);
				if (abbrev_name==sn) abbr_found=true;
			}
			if (!abbr_found) exit_error("Abbreviated Output name of "+inf[3]+" is not among the short names.");
		}
//		vText(REL) = "rCode"+DLMTR+"p_m"+DLMTR+"Rlat";
		vText(L5N) = vText(PID)+DLMTR+vText(IID)+DLMTR+vText(FTH)+DLMTR+vText(MTH)+DLMTR+vText(SEX);
		vText(L6N) = vText(PID)+DLMTR+vText(IID)+DLMTR+vText(FTH)+DLMTR+vText(MTH)+DLMTR+vText(SEX)+DLMTR+vText(AFF);
		vText(L7N) = vText(PID)+DLMTR+vText(IID)+DLMTR+vText(FTH)+DLMTR+vText(MTH)+DLMTR+vText(SEX)+DLMTR+vText(AFF)+DLMTR+vText(LIA);
		vText(L8N) = vText(PID)+DLMTR+vText(IID)+DLMTR+vText(FTH)+DLMTR+vText(MTH)+DLMTR+vText(SEX)+DLMTR+vText(AFF)+DLMTR+vText(AL1)+DLMTR+vText(AL2);
		vText(L9N) = vText(PID)+DLMTR+vText(IID)+DLMTR+vText(FTH)+DLMTR+vText(MTH)+DLMTR+vText(SEX)+DLMTR+vText(AFF)+DLMTR+vText(LIA)+DLMTR+vText(AL1)+DLMTR+vText(AL2);
		vText(L5D) = vText(DNP)+DLMTR+vText(DNI)+DLMTR+vText(DNF)+DLMTR+vText(DNM)+DLMTR+vText(SEX);
		vText(L6D) = vText(DNP)+DLMTR+vText(DNI)+DLMTR+vText(DNF)+DLMTR+vText(DNM)+DLMTR+vText(SEX)+DLMTR+vText(AFF);
		vText(L7D) = vText(DNP)+DLMTR+vText(DNI)+DLMTR+vText(DNF)+DLMTR+vText(DNM)+DLMTR+vText(SEX)+DLMTR+vText(AFF)+DLMTR+vText(LIA);
		vText(L8D) = vText(DNP)+DLMTR+vText(DNI)+DLMTR+vText(DNF)+DLMTR+vText(DNM)+DLMTR+vText(SEX)+DLMTR+vText(AFF)+DLMTR+vText(AL1)+DLMTR+vText(AL2);
		vText(L9D) = vText(DNP)+DLMTR+vText(DNI)+DLMTR+vText(DNF)+DLMTR+vText(DNM)+DLMTR+vText(SEX)+DLMTR+vText(AFF)+DLMTR+vText(LIA)+DLMTR+vText(AL1)+DLMTR+vText(AL2);
	}

	void get_var_name_db(const string& filename) {
		if (!filename.empty() && FileExists(filename))
		{
			lns<<showl<<"Read variable names from "<<filename<<flush_logger;
			fstream input(filename.c_str());
			_get_var_name_db(input);
			input.close();
		}
		else
		{
			stringstream input(default_vname);
			_get_var_name_db(input);
		}
	}
	
	degenerated_t get_id_from(const string& input) {
		set<int> result;
		string s=to_lower_copy(input.substr(0,4));
		for (each_element(vInText,it))
			if (it->second.find(s)!=it->second.end()) result.insert(it->first);
		if (result.size()>1) exit_error("Duplicated variable short names : "+input);
		if (result.size()<1) return degenerated_t(NULL,-1);
		return degenerated_t(NULL,*result.begin());
	}
}

Individual& Individual::operator=(const Individual& orig) { exit_error("Individual class cannot be copied."); return *this; }
Individual::Individual(const Individual& othr) { exit_error("Individual class cannot be copied."); }

void Individual::bf_set(const iData::degenerated_t& var_label, const string& input)
{
	string value = iData::convert(var_label,input);
	if (iData::to_ignore_text(var_label,value)) return;
	if		(iData::is_str_t(var_label))   str[iData::vCode(var_label)]=value;
	else if (iData::is_dbl_t(var_label)) { double d=boost::lexical_cast<double>(value); num[iData::vCode(var_label)]=d; }
	else if (iData::is_trb_t(var_label))   trb[iData::vCode(var_label)]=str2trb(value);
	else exit_error("variable ID number out of bound");
}
// brute-force output data, good for writing data to a file
string Individual::bf_get(const iData::degenerated_t& var_label)
{
	if		(iData::is_str_t(var_label)) return str[iData::vCode(var_label)];
	else if (iData::is_dbl_t(var_label)) return boost::lexical_cast<string>(num[iData::vCode(var_label)]);
	else if (iData::is_trb_t(var_label)) return trb2str(trb[iData::vCode(var_label)]);
	else exit_error("variable ID number out of bound");
	return "This will never happen.";
}
void Individual::copy_data( const Individual& orig )
{
	if(&orig == this) return;
	for (each_element(orig.str, itstr)) this->str[itstr->first]=itstr->second;
	for (each_element(orig.num, itnum)) this->num[itnum->first]=itnum->second;
	for (each_element(orig.trb, ittrb)) this->trb[ittrb->first]=ittrb->second;
	for (each_element(orig.xtr, itstr)) this->xtr[itstr->first]=itstr->second;
}
void Individual::remove( iData::degenerated_t var_label )
{
	if		(iData::is_str_t(var_label)) str.erase(iData::vCode(var_label));
	else if (iData::is_dbl_t(var_label)) num.erase(iData::vCode(var_label));
	else if (iData::is_trb_t(var_label)) trb.erase(iData::vCode(var_label));
	else exit_error("variable ID number out of bound");
}
