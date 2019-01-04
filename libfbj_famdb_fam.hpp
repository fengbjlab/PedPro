#ifndef family_class
#define family_class

#include "libfbj_famdb_ind.hpp"

struct Family {
	// data updated when a prefile is read
	int						circles;	// has circle
	int						num_clt;	// number of clusters
	std::set<std::string>	consang;	// Consanguineous marriage: "IID1 IID2" "IID2 IID1"
	
	// data updated when a new IND is added
	std::string				pid;		// pedigree ID
	std::set<std::string>	lc_iids;	// lower case IIDs
	IndMap					members;	// family members
	IndMap					founders;	// founders
	IndMap					nonfounders;// non-founders
	
	Family():circles(-1),num_clt(-1){}
};

typedef std::map<std::string, Family *> FamMap;

#endif
