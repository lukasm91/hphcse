// domain.h
#ifndef _DOMAIN_H
#define _DOMAIN_H

#include <ctime>
#include "particle.h"

class Domain{
private:
	const double pi;				// the number pi
	const unsigned n_pt;			// number of particles
	const double eps;				// epsilon (Lennard-Jones)
	const double sigma;				// sigma (Lennard-Jones)
	const int lb;					// lower boundary
	const int ub;					// upper boundary
	Particle *prt;		 			// particles
	double **dist;					// distance matrix
	double get_distance(Particle, Particle);	// get abs distance
	void update_dist();				// get abs distance
	void calc_acc();				// calculate acceleration
	void calc_pos(double t);		// calculate new positions
	void calc_vel(double t);		// calculate new velocities
public:
	Domain(unsigned n_pt, double eps, double sigma, int lb, int ub);	// Constructor
	~Domain();						// Destructor
	double ctr_m[3];				// center of mass
	double tot_am[3];				// total angular momentum
	double tot_lm[3];				// total linear momentum
	void next_timestep(double t);	// go one timestep further
	double calc_Ekin();				// calculate kinetic energy
	double calc_Epot();				// calculate potential energy
	double calc_Etot();				// calculate total energy
	void calc_ctr_m();				// calculate center of mass
	void calc_tot_am();				// calculate total angular momentum
	void calc_tot_lm();				// calculate total linear momentum
};

#endif
