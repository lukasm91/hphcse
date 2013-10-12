// particle.h
#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <cmath>

class Particle {
public:
	const double mass = 1;
	double x1[3], x2[3];	// position vectors
	double v1[3], v2[3];	// velocity vectors
	double a1[3], a2[3];	// acceleration vectors
	double get_v();			// calc v from v2
};

#endif
