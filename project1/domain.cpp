// domain.cpp
#include "domain.h"

// Constructor for domain with particles
Domain::Domain(int n){
	// add particles
	for (int i=1; i<n; i++)
		add_particle();
}

void Domain::add_particle(){
	// initialize particle
	Particle p;
	for (int i=0; i<2; i++) {
		p.p1[i] = uniform(generator);	// position uniformely distributed
		p.v1[i] = normal(generator);	// velocity normal distributed
	}
	// add particle to vector
	particle.push_back(p);
}

void Domain::next_timestep(double t){


}
