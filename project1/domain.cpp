// domain.cpp
#include "domain.h"

// Constructor
Domain::Domain(unsigned n, double e, double s, int l, int u) :
	pi(3.14159265359), n_pt(n), eps(e), sigma(s), lb(l), ub(u)
{
	// reinitialization for drand48()
	srand48( static_cast<unsigned>(std::time(NULL)) );
	// create dynamic arrays and matrices
	prt = new Particle[n];
	dist = new double*[n];
	for (unsigned i=0; i<n; ++i)
		dist[i] = new double[n];
	// random numbers
	unsigned nn = 3*n;	// we need 3 random numbers per particle and initial condition
	unsigned m = (nn%2 == 0) ? nn/2 : (nn+1)/2;
	double *uniform = new double[nn];
	double *normal = new double[2*m];
	// add uniformely distributed numbers to array
	for (unsigned i=0; i<nn; ++i)
		uniform[i] = ((drand48()-0.5)+(ub+lb)/2)*abs(ub-lb);
	// add normal distributed (Box-Muller) numbers to array
	double r1;
	double r2;
	for (unsigned i=0; i<m; ++i) {
//		r1 = ((drand48()-0.5)+(ub+lb)/2)*abs(ub-lb);
//		r2 = ((drand48()-0.5)+(ub+lb)/2)*abs(ub-lb);
//		normal[2*i] = sqrt(-2*log(r1))*cos(2*pi*r2);
//		normal[2*i+1] = sqrt(-2*log(r1))*sin(2*pi*r2);
		r1 = drand48();
		r2 = drand48();
		normal[2*i] = ((sqrt(-2*log(r1))*cos(2*pi*r2)-0.5)+(ub+lb)/2)*abs(ub-lb);
		normal[2*i+1] = ((sqrt(-2*log(r1))*sin(2*pi*r2)-0.5)+(ub+lb)/2)*abs(ub-lb);
	}
	// add prt
	for (unsigned i=0; i<n; ++i) {
		for (unsigned j=0; j<3; ++j) {
			prt[i].x2[j] = uniform[3*i+j];
			prt[i].v2[j] = normal[3*i+j];
		}
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
inline double Domain::get_distance(Particle PA, Particle PB) {
	double r = abs(sqrt( pow((PA.x2[0]-PB.x2[0]), 2)
						+pow((PA.x2[1]-PB.x2[1]), 2)
						+pow((PA.x2[2]-PB.x2[2]), 2)));
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
					f_vec[k] = -24*eps*(2*pow((sigma/r), 12) - pow((sigma/r), 6))*r_vec[k]/pow(r, 2);
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
			prt[i].x2[k] = prt[i].x1[k] + t*prt[i].v1[k] + pow(t, 2)*prt[i].a1[k];
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

// calculate kinetic energy
double Domain::calc_Ekin() {
	double eKin = 0;
	double v;
	for (unsigned i=0; i<n_pt; ++i) {
		v = prt[i].get_v();
		eKin += prt[i].mass/2*pow(v, 2);	
	}
	return eKin;
}

// calculate potential energy
double Domain::calc_Epot() {
	//TODO
	return 0;
}

// calculate total energy
double Domain::calc_Etot() {
	//TODO
	return 0;
}

// calculate center of mass
void Domain::calc_ctr_m() {
	double m = 0;
	// first set all values to 0
	for (unsigned k=0; k<3; ++k)
		ctr_m[k] = 0;
	// calculate...
	for (unsigned i=0; i<n_pt; ++i) {
		m += prt[i].mass;
		for (unsigned k=0; k<3; ++k)
			ctr_m[k] += prt[i].mass*prt[i].x2[k];
	}
	for (unsigned k=0; k<3; ++k)
		ctr_m[k] = ctr_m[k]/m;
}

// calculate total angular momentum
void Domain::calc_tot_am() {
	// set all values to 0
	for (unsigned k=0; k<3; ++k)
		tot_am[k] = 0;
	// calculate...
	for (unsigned i=0; i<n_pt; ++i) {
		tot_am[0] += (prt[i].x2[1]*prt[i].v2[2]-prt[i].x2[2]*prt[i].v2[1])*prt[i].mass;
		tot_am[1] += (prt[i].x2[2]*prt[i].v2[0]-prt[i].x2[0]*prt[i].v2[2])*prt[i].mass;
		tot_am[2] += (prt[i].x2[0]*prt[i].v2[1]-prt[i].x2[1]*prt[i].v2[0])*prt[i].mass;
	}
}

// calculate total linear momentum
void Domain::calc_tot_lm() {
	double tot_lm[3];
	// set all values to 0
	for (unsigned k=0; k<3; ++k)
		tot_lm[k] = 0;
	//calculate...
	for (unsigned i=0; i<n_pt; ++i) {
		for (unsigned k=0; k<3; ++k)
			tot_lm[k] += prt[i].mass*prt[i].v2[k];
	}
}
