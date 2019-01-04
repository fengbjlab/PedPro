#include "libfbj_famdb.hpp"
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	default_tfile_format.forbid_option("-t");
	program.enable_option("--prefix");
	program.read_arguments(argc,argv,true,false);
	program.push_back_help(default_tfile_format.help_text());
	program.check_help_request();
	
	if (!program.quiet) lns.open(program.prefix()+".log", ios::out|ios::app);
	lns<<writel<<"Options specified:\n"<<program.commands()<<flush_logger;
	default_tfile_format.read_arguments(program.arg());
	FamilyDB fdb;
	fdb.read_arguments(program.arg());
	return 0;
}
