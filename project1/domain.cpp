// domain.cpp
#include "domain.h"

// Constructor
Domain::Domain(int n) : n_pt(n) {
	// create arrays for particles and random numbers
	*particles = new Particle[n];
	nn = 3*n;	// we need 3 random numbers per particle and initial condition
	int m = (nn%2 == 0) ? nn/2 : (nn+1)/2;
	double *uniform = new double[nn];
	double *normal = new double[2*m];
	// add uniformely distributed numbers to array
	for (int i=0; i<nn; ++i)
		uniform[i] = ((drand48()-0.5)+(ub+lb)/2)*abs(ub-lb);
	// add normal distributed (Box-Muller) numbers to array
	double r1;
	double r2;
	for (int i=0; i<m; ++i) {
		r1 = ((drand48()-0.5)+(ub+lb)/2)*abs(ub-lb);
		r2 = ((drand48()-0.5)+(ub+lb)/2)*abs(ub-lb);
		normal[2*i] = sqrt(-2*log(r1))*cos(2*pi*r2);
		normal[2*i+1] = sqrt(-2*log(r1))*sin(2*pi*r2);
	}
	// add particles
	for (int i=0; i<n; ++i) {
		for (int j=0; j<3; ++j) {
			particles[i].p1[j] = uniform[3*i+j];
			particles[i].v1[j] = normal[3*i+j];
		}
	}
	// calc acceleration
	calc_acc();
	// clean up
	delete[] uniform;
	delete[] normal;
}

// Destructor 
Domain::~Domain() {
	delete[] particles;
}

double get_distance(Particle PA, Particle PB) {
	double r = abs(sqrt( (PA.p1[0]-PB.p1[0])^2
						+(PA.p1[1]-PB.p1[1])^2
						+(PA.p1[2]-PB.p1[2])^2 ));
	return r;

}

void Domain::calc_acc() {
	double r_vec[3];	// distance vector between 2 particles
	double f_vec[3];	// force vector between 2 particles
	double r;			// abs distance between 2 particles
	for (unsigned i=0; i < n_pt; ++i) {
		for (unsigned j=0; j < n_pt; ++j) {
			if (i != j) {
				// calculate force and acceleration on particle
				r_vec[0] = particles[i].p1[0] - particles[j].p1[0];
				r_vec[1] = particles[i].p1[1] - particles[j].p1[1];
				r_vec[2] = particles[i].p1[2] - particles[j].p1[2];
				r = get_distance(particles[i], particles[j]);
				for (unsigned k=0; k < 3; ++k) {
					f_vec[k] =  f_vec[k]
							-24*eps*(2*(s_lj/r)^12 - (s_lj/r)^6)*r_vec[k]/r^2;
				}
			}
		}
		for (unsigned k=0; k<3; ++k)
			particles[i].a2[k] += f_vec[k];
	}
}

void Domain::calc_pos(double t) {
// TODO
}

void Domain::calc_vel(double t) {
// TODO
}

void Domain::next_timestep(double t) {
// TODO
}
