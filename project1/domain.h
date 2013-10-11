// domain.h
#ifndef _DOMAIN_H
#define _DOMAIN_H

#include <ctime>
#include <cstdlib>
#include <cmath>

// some constants
#define lb = -1;
#define ub = 1;
#define eps = 1;
#define pi = 3.14159265359;

// reinitialization for drand48()
srand48( static_cast<unsigned>(std::time(NULL)) );


struct Particle{
	const double mass = 1;
	// position vectors
	double x1[3], x2[3];	// position arrays
	double v1[3], v2[3];	// velocity arrays
	double a1[3], a2[3];	// acceleration arrays
};

class Domain{
private:
	Particles *prt;		 			// particles
	int n_pt;						// number of particles
	double r_vec[3];				// distance vector between 2 particles
	double f_vec[3];				// force vector between 2 particles
	void calc_acc();				// calculate acceleration
	void calc_pos(double t);		// calculate new positions
	void calc_vel(double t);		// calculate new velocities
public:
	Domain(int);					// Constructor
	~Domain();						// Destructor
	void next_timestep(double t);	// go one timestep further
	double get_distance(Particle, Particle);			// get abs distance
};

#endif
