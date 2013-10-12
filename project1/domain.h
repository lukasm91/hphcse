// domain.h
#ifndef _DOMAIN_H
#define _DOMAIN_H

#include <ctime>
#include <cstdlib>
#include "particle.h"

#define pi = 3.14159265359;

class Domain{
private:
	const int eps = 1;				// epsilon (Lennard-Jones)
	const int sigma = 1;			// sigma (Lennard-Jones)
	const int lb = -1;				// lower boundary
	const int ub = 1;				// upper boundary
	Particles *prt;		 			// particles
	double *dist;					// distance matrix
	int n_pt;						// number of particles
	double get_distance(Particle, Particle);	// get abs distance
	void update_dist();				// get abs distance
	void calc_acc();				// calculate acceleration
	void calc_pos(double t);		// calculate new positions
	void calc_vel(double t);		// calculate new velocities
public:
	Domain(int);					// Constructor
	~Domain();						// Destructor
	void next_timestep(double t);	// go one timestep further
	double calc_Ekin();				// calculate kinetic energy
	double calc_Epot();				// calcluate potential energy
};

#endif
