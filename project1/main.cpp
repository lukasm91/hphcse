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
	// files for output
	std::ofstream distance("distance.out", std::ios::out);
	std::ofstream outfile("nbody.out", std::ios::out);
	std::ofstream prtrack("prttrack.out", std::ios::out);
	prtrack << "i" << "\t" << "r_x" << "\t" << "r_y" << "\t" << "r_z"
			<< "\t" << "v_x" << "v_y" << "v_z" << std::endl;
	outfile << "i" << "\t" << "E_Kin" << std::endl;
	distance << "i" << "\t" << "r_x" << "\t" << "r_y" << "\t" << "r_z"
			<< "\t" << "v_x" << "v_y" << "v_z" << std::endl;
	// create domain
	Domain d = Domain(n, epsilon, sigma, lb, ub);
	// inster distances
	for (unsigned i=0; i<n; ++i) {
		for (unsigned j=0; j<n; ++j) {
			distance << d.dist[i][j] << "\t";
		}
		distance << std::endl;
	}

	// timesteps
	for (unsigned i=0; i<t; ++i) {
		outfile << i << "\t" << d.calc_Ekin() << std::endl;
		prtrack << i << "\t" << d.prt[0].x2[0] << "\t" << d.prt[0].x2[1] << "\t"
				<< d.prt[0].x2[2] << "\t" << d.prt[0].v2[0] << "\t"
				<< d.prt[0].v2[1] << "\t" << d.prt[0].v2[2] << std::endl;
		d.next_timestep(dt);
	}
	outfile << t << "\t" << d.calc_Ekin() << std::endl;
	return 0;
}
