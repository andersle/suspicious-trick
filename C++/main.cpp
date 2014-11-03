#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "trajectory.h"

using namespace std::chrono;

int main(int argc, char** argv) {
	/* set up the instrumentation */

	high_resolution_clock::time_point start = high_resolution_clock::now();
	
	std::string filename(argv[1]);

	std::vector<Atoms::Property> properties = { Atoms::Property::ID,
		Atoms::Property::TYPE,
		Atoms::Property::X, Atoms::Property::Y, Atoms::Property::Z};
	Trajectory t(filename, properties);
	
	Atoms a = t.readFrame();
	int tsteps_processed = 0;
	int tsteps_read = 0;

	while (a.errorflag == Atoms::error::NO_ERROR) {
		tsteps_read++;
		tsteps_processed++;
		std::cout << a.timestep << std::endl;

		a = t.readFrame();
	}
	
	switch (a.errorflag) {
		case Atoms::error::NO_ERROR:
			break;
		case Atoms::error::END_OF_FILE:
			break;
		case Atoms::error::FILE_ERROR:
			std::cerr << "File error (will continue with what we"
				<< " had)" << std::endl;
			break;
		case Atoms::error::TRICLINIC_BOX:
			std::cerr << "Triclinic boxes are unsupported." <<
				std::endl;
			return 1;
			break;
		case Atoms::error::BAD_BOUNDARY:
			std::cerr << "Unsupported boundary type (not p,s,f,m)."
				<< std::endl;
			return 1;
			break;
		case Atoms::error::BAD_PROPERTY_COUNT:
			std::cerr << "The file contains " << a.num_fields <<
				" fields, but the property vector contains " <<
				properties.size() << "fields." << std::endl;
			return 1;
			break;
		case Atoms::error::FILE_CORRUPT:
			std::cerr << "The reported buffersize is not "
				"compatible with the reported number of fields"
				" per atom" << std::endl;
			return 1;
			break;
	}

	high_resolution_clock::time_point end = high_resolution_clock::now();

	duration<double> time_taken =
		duration_cast<duration<double>>(end-start);

	std::cout << "Processed " << tsteps_processed << "/" << tsteps_read << 
	       " frames in " << time_taken.count() << " seconds. " << std::endl;
	std::cout << "(" << tsteps_processed / time_taken.count() << " frames per second)" <<
		std::endl;

	return 0;
}
