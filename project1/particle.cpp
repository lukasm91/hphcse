// particle.cpp
#include "particle.h"

// Constructor
Particle::Particle() : mass(1) { }

// get velocity of particle
double Particle::get_v() {
	return abs(sqrt(pow(v2[0],2) + pow(v2[1],2) + pow(v2[2],2)));
}
