// particle.cpp
#include "particle.h"

inline double Particle::get_v() {
	return abs(sqrt(v2[0]^2 + v2[1]^2 + v2[2]^2));
}
