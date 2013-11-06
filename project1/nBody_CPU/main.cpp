// N-particle problem for CPU
// (c) Christian Zeman
// 08.10.2013

#include <iostream>
#include <fstream>
#include "domain.h"
#include "ArgumentParser.h"
#include "Timer.h"

// ================================================================
// output function
// ================================================================
void prtout(std::ofstream &out, Domain &d, unsigned i, unsigned n, bool hdr) {
	if (hdr) {
		out << "it" << "\t" << "prt" << "\t" << "x_1" << "\t" << "x_2" << "\t" << "x_3"
			<< "\t" << "v_1" << "\t" << "v2" << "\t" << "v3" << "\t"
			<< "a_1" << "\t" << "a_2" << "\t" << "a_3" << std::endl;
	} else {
		for (unsigned j=0; j<n; ++j) {
			out << i << "\t" << j << "\t" << d.prt_x[3*(2*j+1)] << "\t"
				<< d.prt_x[3*(2*j+1)+1]
				<< "\t" << d.prt_x[3*(2*j+1)+2]
				<< "\t" << d.prt_v[3*(2*j+1)] << "\t" << d.prt_v[3*(2*j+1)+1]
				<< "\t" << d.prt_v[3*(2*j+1)+2]
				<< "\t" << d.prt_a[3*(2*j+1)] << "\t" << d.prt_a[3*(2*j+1)+1]
				<< "\t" << d.prt_a[3*(2*j+1)+2] << std::endl;
		}
	}
}

// ================================================================
// nBody simulation
// ================================================================
void nBody(unsigned t, unsigned n, bool print)
{	
	// constants
	const double m = 1;				// mass of a particle
	const double epsilon = 1;		// epsilon (Lennard-Jones)
	const double sigma = 0.2;		// sigma (Lennard-Jones)
	const double lb = -1;			// lower boundary
	const double ub = 1;			// upper boundary
	const double dt = 0.0001;		// timestep size
	unsigned n3 = pow(n,3);			// #particle = n^3

	// create domain
	Domain d = Domain(n, m, dt, epsilon, sigma, lb, ub);
	
	// files for output
	std::ofstream outfile("nbody.out", std::ios::out);
	std::ofstream prtrack("prttrack.out", std::ios::out);
	outfile.precision(4);
	prtrack.precision(4);
	
	if (print) {
		prtout(prtrack, d, 0, n3, 1);
		outfile << std::fixed
				<< "i" << "\t" << "E_Kin" << "\t" << "E_Pot" << "\t" << "E_Tot"
				<< "\t" << "C_x" << "\t\t" << "C_y" << "\t\t" << "C_z"
				<< "\t\t" << "P_x" << "\t\t" << "P_y" << "\t\t" << "P_z"
				<< "\t\t" << "L_x" << "\t\t" << "L_y" << "\t\t" << "L_z" << std::endl;
		outfile << std::endl;
	}

	// timesteps
	for (unsigned i=0; i<t; ++i) {
		if (print) {
			prtout(prtrack, d, i, n3, 0);
			outfile << std::fixed
					<< i << "\t" << d.eKin << "\t" << d.ePot << "\t" << d.eTot
					<< "\t" << d.ctr_m[0] << "\t" << d.ctr_m[1] << "\t" << d.ctr_m[2]
					<< "\t" << d.tot_lm[0] << "\t" << d.tot_lm[1] << "\t" << d.tot_lm[2] 
					<< "\t" << d.tot_am[0] << "\t" << d.tot_am[1] << "\t" << d.tot_am[2] << std::endl;
		}
		d.next_timestep();
	}
	if (print) {
	prtout(prtrack, d, t, n3, 0);
		outfile << std::fixed
				<< t << "\t" << d.eKin << "\t" << d.ePot << "\t" << d.eTot
				<< "\t" << d.ctr_m[0] << "\t" << d.ctr_m[1] << "\t" << d.ctr_m[2]
				<< "\t" << d.tot_lm[0] << "\t" << d.tot_lm[1] << "\t" << d.tot_lm[2] 
				<< "\t" << d.tot_am[0] << "\t" << d.tot_am[1] << "\t" << d.tot_am[2] << std::endl;
	}
}

// ================================================================
// main with timer
// ================================================================

int main (int argc, char *argv[])
{
	// arguments
	ArgumentParser parser(argc, argv);
	const unsigned it = parser("-i").asInt(1);			// number of iterations, default 1
	const unsigned t = parser("-t").asInt(100);			// number of timesteps
	unsigned n = parser("-n").asInt(2);					// initial number of particles per line, default 2
	const bool print = parser("-p").asBool(1);			// print outfile?
	
	// timing
	std::ofstream times("timing.out", std::ios::out);
	times << "i" << "\t" << "n" << "\t" << "n^3" << "\t" << "time" << std::endl;
	times << "-------------------------------------------" << std::endl;
	Timer timer; // declaration of the timer

	for (unsigned i=0; i<it; ++i) {
		times << i << "\t" << n << "\t" << pow(n,3) << "\t";
		timer.start();
		nBody(t, n, print);
		times << timer.stop() << std::endl;
		n++;
	}
		
	return 0;
}


