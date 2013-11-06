// domain.h
#ifndef _DOMAIN_H
#define _DOMAIN_H

#include <ctime>
#include <cmath>
#include <cstdlib>
//#include "particle.h"

class Domain{
private:
	const double pi_;				// the number pi
	const unsigned n_;				// number of particles
	const unsigned l_;				// number of particles per line
	const double m_;				// mass of particle
	const double t_;				// timestep size
	const double e_;				// epsilon (Lennard-Jones)
	const double s_;				// sigma (Lennard-Jones)
	const double lb_;				// lower boundary
	const double ub_;				// upper boundary
	const double rc_;				// cut off radius
	const double U_ljc_;			// Lennard-Jones potential at cut off radius
	double eTotT;					// target total energy (has to remain constant)
	double get_rf();				// get rescaling factor for velocities
	void calc_acc();				// calculate acceleration
	void calc_pos();				// calculate new positions
	void calc_vel();				// calculate new velocities
	void resc_vel();				// rescale velocities
	double get_Ekin();				// calculate kinetic energy
	double get_Epot();				// calculate potential energy
	double get_Etot();				// calculate total energy
	void calc_ctr_m();				// calculate center of mass
	void calc_tot_am();				// calculate total angular momentum
	void calc_tot_lm();				// calculate total linear momentum
	void initial_grid();			// set up particles on grid
public:
	Domain(unsigned n, double m, double t, double e, double s, double lb, double ub);	// constructor
	~Domain();						// destructor
	double eTot;					// total energy
	double ePot;					// potential energy
	double eKin;					// kinetic energy
	double *prt_x;					// positions of particles (x1, y1, z1, x2, y2, z2)
	double *prt_v;					// velocities of particles (u1, v1, w1, u2, v2, w2)
	double *prt_a;					// accelerations of particles (a1, b1, c1, a2, b2, c2)
	double ctr_m[3];				// center of mass
	double tot_am[3];				// total angular momentum
	double tot_lm[3];				// total linear momentum
	void next_timestep();			// go one timestep further
};

#endif
