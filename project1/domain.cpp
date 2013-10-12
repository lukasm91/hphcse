// domain.cpp
#include "domain.h"

// Constructor
Domain::Domain(int n) : n_pt(n) {
	// reinitialization for drand48()
	srand48( static_cast<unsigned>(std::time(NULL)) );
	// create dynamic arrays and matrices
	*prt = new Particle[n];
	*dist = new double[n][n];
	*force = new double[n][n];
	// random numbers
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
	// add prt
	for (int i=0; i<n; ++i) {
		for (int j=0; j<3; ++j) {
			prt[i].x2[j] = uniform[3*i+j];
			prt[i].v2[j] = normal[3*i+j];
		}
		set_velocity(prt[i]);
	}
	// insert distances into distance matrix
	update_dist();
	// calc acceleration
	calc_acc();
	// clean up
	delete[] uniform;
	delete[] normal;
}

// Destructor 
Domain::~Domain() {
	delete[] prt;
}

// get distance between two particles
double Domain::get_distance(Particle PA, Particle PB) {
	double r = abs(sqrt( (PA.x2[0]-PB.x2[0])^2
						+(PA.x2[1]-PB.x2[1])^2
						+(PA.x2[2]-PB.x2[2])^2 ));
	return r;
}

// update distance matrix
void Domain::update_dist() {
	for (unsigned i=0; i<n_pt; ++i) {
		for (unsigned j=i; j<n_pt; ++j) {
			dist[i][j] = dist[j][i] = get_distance(prt[i], prt[j]);
		}
	}
}

// calculate acceleration
void Domain::calc_acc() {
	double r_vec[3];	// distance vector between 2 particles
	double f_vec[3];	// force vector between 2 particles
	double r;			// abs distance between 2 particles
	for (unsigned i=0; i<n_pt; ++i) {
		for (unsigned j=0; j<n_pt; ++j) {
			if (i!=j) {
				// calculate force and acceleration on particle
				r_vec[0] = prt[i].x1[0] - prt[j].x1[0];
				r_vec[1] = prt[i].x1[1] - prt[j].x1[1];
				r_vec[2] = prt[i].x1[2] - prt[j].x1[2];
				r = dist[i][j];
				for (unsigned k=0; k<3; ++k) {
					f_vec[k] = -24*eps*(2*(s_lj/r)^12 - (s_lj/r)^6)*r_vec[k]/r^2;
					prt[i].a2[k] += f_vec[k];
				}
			}
		}
		for (unsigned k=0; k<3; ++k)
			prt[i].a2[k] = prt[i].a2[k]/prt[i].mass;
	}
}

// calculate position
void Domain::calc_pos(double t) {
	for (unsigned i=0; i<n_pt; ++i) {
		for (unsigned k=0; k<3; ++k) {
			prt[i].x2[k] = prt[i].x1[k] + t*prt[i].v1[k] + t^2*prt[i].a1[k];
			// particles leaving the domain will enter it on the other side
			while (prt[i].x2[k]>ub)
				prt[i].x2[k] -= (ub-lb);
			while (prt[i].x2[k]<lb)
				prt[i].x2[k] += (ub-lb);
		}
	}
	// update distance matrix with new positions
	update_dist();
}

// calculate velocity
void Domain::calc_vel(double t) {
	for (unsigned i=0; i<n_pt; ++i) {
		for (unsigned k=0; k<3; ++k)
			prt[i].v2[k] = prt[i].v1[k] + t/2*(prt[i].a1[k]+prt[i].a2[k]);
	}
}

// go one timestemp further
void Domain::next_timestep(double t) {
	// new values are now old values
	for (unsigned i=0; i<n_pt; ++i) {
		for (unsigned k=0; k<3; ++k) {
			prt[i].x1[k] = prt[i].x2[k];
			prt[i].v1[k] = prt[i].v2[k];
			prt[i].a1[k] = prt[i].a2[k];
			prt[i].x2[k] = 0;
			prt[i].v2[k] = 0;
			prt[i].a2[k] = 0;
		}
	}
	// new position
	calc_pos(t);
	// new acceleration
	calc_acc();
	// new velocity
	calc_vel(t);
}

double Domain::calc_Ekin() {
	double eKin = 0;
	double v;
	for (unsigned i=0; i<n_pt; ++i) {
		v = get_v();
		eKin += prt[i].mass/2*v^2;	
	}
	return eKin;
}

double Domain::calc_Epot() {
	//TODO
}




