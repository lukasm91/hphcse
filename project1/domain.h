// domain.h
#ifndef _DOMAIN_H
#define _DOMAIN_H

#include <vector>
#include <random>

// some constants
#define lb = -1;
#define ub = 1;
#define meanv = 0;
#define sigma = 0.5;

struct Particle{
	const double mass = 1;
	// position vectors
	std::vector<double>(3) p1;  // x_i
	std::vector<double>(3) p2;  // x_i+1
	// velocity vectors
	std::vector<double>(3) v1;  // v_i
	std::vector<double>(3) v2;  // v_i+1
	// acceleration vectors
	std::vector<double>(3) a1;  // a_i
	std::vector<double>(3) a2;  // a_i+1
};

class Domain{
private:
	// generate random numbers
	std::default_random_engine generator;
	std::uniform_real_distribution<double> uniform(lb, ub);
	std::normal_distribution<double> normal(meanv, sigma);

	std::vector<Particle> particles; 	// particles
	void calc_acc();					// calculate accelerations
	void calc_pos(double);		// calculate new positions
	void calc_vel(double);		// calculate new velocities
public:
	Domain(int);					// Constructor
	void add_particle();			// add particle to domain
	void next_timestep(double);		// go one timestep further
};

#endif
