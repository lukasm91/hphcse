// N-particle problem for CPU
// (c) Christian Zeman
// 08.10.2013

#include <iostream>
#include <fstream>
#include "domain.h"
#include "ArgumentParser.h"

void prtout(std::ofstream &out, Domain &d, unsigned i, unsigned n) {
	out << "Timestep " << i << std::endl;
	out << "-----------------------------------------------------" << std::endl;
	for (unsigned j=0; j<n; ++j) {
		out << "Particle " << j << ":" << std::endl;
		out << "\t" << "x:" << "\t" << d.prt_x[3*(2*j+1)] << "\t" << d.prt_x[3*(2*j+1)+1]
				<< "\t" << d.prt_x[3*(2*j+1)+2] << std::endl;
		out << "\t" << "v:" << "\t" << d.prt_v[3*(2*j+1)] << "\t" << d.prt_v[3*(2*j+1)+1]
				<< "\t" << d.prt_v[3*(2*j+1)+2] << std::endl;
		out << "\t" << "a:" << "\t" << d.prt_a[3*(2*j+1)] << "\t" << d.prt_a[3*(2*j+1)+1]
				<< "\t" << d.prt_a[3*(2*j+1)+2] << std::endl;
	}
	out << std::endl;
}


int main(int argc, char *argv[])
{	
	ArgumentParser parser(argc, argv);
	const unsigned t = parser("-t").asInt(10);		// number of steps, default 10
	const unsigned n = parser("-n").asInt(10);		// number of particles, default 10

	// constants
	const double m = 1;				// mass of a particle
	const double epsilon = 1;		// epsilon (Lennard-Jones)
	const double sigma = 0.1;		// sigma (Lennard-Jones)
	const double lb = -1;			// lower boundary
	const double ub = 1;			// upper boundary
	const double dt = 0.0001;		// timestep size

	// files for output
	std::ofstream outfile("nbody.out", std::ios::out);
// 	std::ofstream prtrack("prttrack.out", std::ios::out);
	outfile.precision(4);
// 	prtrack.precision(4);

	// create domain
	Domain d = Domain(n, m, dt, epsilon, sigma, lb, ub);
	outfile << std::fixed
			<< "i" << "\t" << "E_Kin" << "\t" << "E_Pot" << "\t" << "E_Tot"
			<< "\t" << "C_x" << "\t\t" << "C_y" << "\t\t" << "C_z"
			<< "\t\t" << "P_x" << "\t\t" << "P_y" << "\t\t" << "P_z"
			<< "\t\t" << "L_x" << "\t\t" << "L_y" << "\t\t" << "L_z" << std::endl;
	outfile << std::endl;

	// timesteps
	for (unsigned i=0; i<t; ++i) {
// 		prtout(prtrack, d, i, n);
		outfile << std::fixed
				<< i << "\t" << d.eKin << "\t" << d.ePot << "\t" << d.eTot
				<< "\t" << d.ctr_m[0] << "\t" << d.ctr_m[1] << "\t" << d.ctr_m[2]
				<< "\t" << d.tot_lm[0] << "\t" << d.tot_lm[1] << "\t" << d.tot_lm[2] 
				<< "\t" << d.tot_am[0] << "\t" << d.tot_am[1] << "\t" << d.tot_am[2] << std::endl;
		d.next_timestep();
	}
// 	prtout(prtrack, d, t, n);
	outfile << std::fixed
			<< t << "\t" << d.eKin << "\t" << d.ePot << "\t" << d.eTot
			<< "\t" << d.ctr_m[0] << "\t" << d.ctr_m[1] << "\t" << d.ctr_m[2]
			<< "\t" << d.tot_lm[0] << "\t" << d.tot_lm[1] << "\t" << d.tot_lm[2] 
			<< "\t" << d.tot_am[0] << "\t" << d.tot_am[1] << "\t" << d.tot_am[2] << std::endl;
	return 0;
}

