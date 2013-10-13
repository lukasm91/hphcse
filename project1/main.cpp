// N-particle problem for CPU
// (c) Christian Zeman
// 08.10.2013

#include <iostream>
#include <fstream>
#include "domain.h"

int main()
{	
	// constants
	const unsigned n = 10;			// number of particles
	const double epsilon = 1;		// epsilon (Lennard-Jones)
	const double sigma = 1;			// sigma (Lennard-Jones)
	const int lb = -1;				// lower boundary
	const int ub = 1;				// upper boundary
	const unsigned t = 100;			// number of timesteps
	const double dt = 1;			// timestep size
	// file for output
	std::ofstream outfile("nbody.out", std::ios::out);
	outfile << "Timestep" << "\t" << "E_Kin" << std::endl;
	// create domain
	Domain d = Domain(n, epsilon, sigma, lb, ub);
	// timesteps
	for (unsigned i=0; i<t; ++i) {
		outfile << i << "\t" << d.calc_Ekin() << std::endl;
		d.next_timestep(dt);
	}
	outfile << t << "\t" << d.calc_Ekin() << std::endl;
	return 0;
}
