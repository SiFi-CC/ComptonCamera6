#include "CLog.hh"




int main(int argc, char** argv) {
  spdlog::set_level(spdlog::level::info);

    spdlog::info("Hello, Vitalii!");
	if (argc != 3) {
	    spdlog::info(
	        "type: './cm_reconstruct [INPUTFILENAME] [OUTPUTFILENAME]' to start:\n\n"
	        "where:\n\n"
	        "INPUTFILENAME - is an input file from simulations\n\n"
	        "OUTPUTFILENAME - is an output file for H matrix\n\n");
	    return 1;
  	}		

}