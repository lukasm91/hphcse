// domain.cpp
#include "domain.h"
#include <iostream> // debugging

// Constructor
Domain::Domain(unsigned n, double m, double t, double e, double s, double lb, double ub) :
	pi_(3.14159265359), n_(pow(n,3)), l_(n), m_(m), t_(t), e_(e), s_(s), lb_(lb), ub_(ub), rc_(2.5*s),
	U_ljc_(4*e*(pow(0.4, 12)-pow(0.4, 6)))
{
	// reinitialization for drand48()
	srand48( static_cast<unsigned>(std::time(NULL)) );
	// create dynamic arrays and initialize them with 0
	prt_x = new double[n_*6]();
	prt_v = new double[n_*6]();
	prt_a = new double[n_*6]();
	// random numbers
	unsigned nn = 3*n_;	// we need 3 random numbers per particle
	unsigned mm = (nn%2 == 0) ? nn/2 : (nn+1)/2;
	double *normal = new double[2*mm];
	// place particles on a grid
// 	unsigned l = unsigned(pow(n_,(double)1/3)+0.5);	// particles per line
	double d = fabs(ub_-lb_)/l_;		// distance between particles
	for (unsigned i=0; i<l_; ++i) {
		for (unsigned j=0; j<l_; ++j) {
			for (unsigned k=0; k<l_; ++k) {
				prt_x[3*(2*(i*l_*l_+j*l_+k)+1)] = lb_ + d/2 + i*d;		// x2
				prt_x[3*(2*(i*l_*l_+j*l_+k)+1)+1] = lb_ + d/2 + j*d;	// y2
				prt_x[3*(2*(i*l_*l_+j*l_+k)+1)+2] = lb_ + d/2 + k*d; 	// z2
			}
		}
	}
	// add normal distributed (Box-Muller) numbers to array
	double r1;
	double r2;
	for (unsigned i=0; i<mm; ++i) {
		r1 = drand48();
		r2 = drand48();
		normal[2*i] = sqrt(-2*log(r1))*cos(2*pi_*r2);
		normal[2*i+1] = sqrt(-2*log(r1))*sin(2*pi_*r2);
	}
	// set initial velocities of particles
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned j=0; j<3; ++j) {
			prt_v[3*(2*i+1)+j] = normal[3*i+j];	// u2, v2, w2
		}
	}
	// calc acceleration
	calc_acc();
	// calc energies
	eKin = get_Ekin();
	ePot = get_Epot();
	eTot = get_Etot();
	// store total energy for velocity rescaling
	eTotT = eTot;
	// calc other values
	calc_ctr_m();
	calc_tot_am();
	calc_tot_lm();
	// clean up
	delete[] normal;
}

// Destructor 
Domain::~Domain() {
	delete[] prt_x;
	delete[] prt_v;
	delete[] prt_a;
}

// calculate acceleration
void Domain::calc_acc() {
	double r_vec[3];	// distance vector between 2 particles
	double r;			// abs distance between 2 particles
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned j=0; j<n_; ++j) {
			if (i!=j) {
				// distance vectors
				r_vec[0] = prt_x[3*(2*j+1)] - prt_x[3*(2*i+1)];
				r_vec[1] = prt_x[3*(2*j+1)+1] - prt_x[3*(2*i+1)+1];
				r_vec[2] = prt_x[3*(2*j+1)+2] - prt_x[3*(2*i+1)+2];
				// distance
				r = sqrt(pow(r_vec[0],2) + pow(r_vec[1],2) + pow(r_vec[2],2));
				// check if particle is in cut off radius
				if (r <= rc_) {
					// add force to acceleration
					for (unsigned k=0; k<3; ++k)
						prt_a[3*(2*i+1)+k] += -24*e_*(2*pow((s_/r),12) - pow((s_/r),6))*r_vec[k]/pow(r,2)/m_;
				}
			}
		}
	}
}

// calculate position
void Domain::calc_pos() {
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned k=0; k<3; ++k) {
			prt_x[3*(2*i+1)+k] = prt_x[3*2*i+k] + t_*prt_v[3*2*i+k] + pow(t_,2)*prt_a[3*2*i+k];
			// particles leaving the domain will enter it on the other side
			while (prt_x[3*(2*i+1)+k] > ub_)
				prt_x[3*(2*i+1)+k] -= (ub_-lb_);
			while (prt_x[3*(2*i+1)+k] < lb_)
				prt_x[3*(2*i+1)+k] += (ub_-lb_);
		}
	}
}

// calculate velocity
void Domain::calc_vel() {
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned k=0; k<3; ++k)
			prt_v[3*(2*i+1)+k] = prt_v[3*2*i+k] + t_/2*(prt_a[3*2*i+k]+prt_a[3*(2*i+1)+k]);
	}
}

// set recaling factor for velocities
double Domain::get_rf() {
	double eKinT = eTotT - ePot;
	return sqrt(eKinT/eKin);
}

// rescale velocities
void Domain::resc_vel() {
	double rf = get_rf();
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned k=0; k<3; ++k)
			prt_v[3*(2*i+1)+k] = rf * prt_v[3*(2*i+1)+k];
	}
}

// go one timestemp further
void Domain::next_timestep() {
	// new values are now old values
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned k=0; k<3; ++k) {
			prt_x[3*2*i+k] = prt_x[3*(2*i+1)+k];
			prt_v[3*2*i+k] = prt_v[3*(2*i+1)+k];
			prt_a[3*2*i+k] = prt_a[3*(2*i+1)+k];
			prt_x[3*(2*i+1)+k] = 0;
			prt_v[3*(2*i+1)+k] = 0;
			prt_a[3*(2*i+1)+k] = 0;
		}
	}
	// new position
	calc_pos();
	ePot = get_Epot();
	// new acceleration
	calc_acc();
	// new velocity
	calc_vel();
	resc_vel();		// rescale
	eKin = get_Ekin();
	eTot = get_Etot();
	// calc other values
	calc_ctr_m();
	calc_tot_am();
	calc_tot_lm();
}

// calculate kinetic energy
double Domain::get_Ekin() {
	double eK = 0;
	double v;
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned k=0; k<3; ++k)
			v += pow(prt_v[3*(2*i+1)+k], 2);
		v = sqrt(v);
		eK += pow(v, 2);	
	}
	eK *= m_/2;
	return eK;
}

// calculate potential energy
double Domain::get_Epot() {
	double r;
	double eP = 0;
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned j=i+1; j<n_; ++j) {
			r = sqrt(pow(prt_x[3*(2*i+1)]-prt_x[3*(2*j+1)], 2)
					+pow(prt_x[3*(2*i+1)+1]-prt_x[3*(2*j+1)+1], 2)
					+pow(prt_x[3*(2*i+1)+2]-prt_x[3*(2*j+1)+2], 2));
			eP += 4*e_*(pow(s_/r, 12)-pow(s_/r, 6)) - U_ljc_; 
		}
	}
	return eP;
}

// calculate total energy
double Domain::get_Etot() {
	return ePot + eKin;
}

// calculate center of mass
void Domain::calc_ctr_m() {
	// first set all values to 0
	for (unsigned k=0; k<3; ++k)
		ctr_m[k] = 0;
	// loop through particles
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned k=0; k<3; ++k)
			ctr_m[k] += m_*prt_x[3*(2*i+1)+k];
	}
	// divide by total mass
	for (unsigned k=0; k<3; ++k)
		ctr_m[k] = ctr_m[k]/(n_*m_);
}

// calculate total angular momentum
void Domain::calc_tot_am() {
	// set all values to 0
	for (unsigned k=0; k<3; ++k)
		tot_am[k] = 0;
	// calculate...
	for (unsigned i=0; i<n_; ++i) {
		tot_am[0] += prt_x[3*(2*i+1)+1]*prt_v[3*(2*i+1)+2]-prt_x[3*(2*i+1)+2]*prt_v[3*(2*i+1)+1];
		tot_am[1] += prt_x[3*(2*i+1)+2]*prt_v[3*(2*i+1)+0]-prt_x[3*(2*i+1)+0]*prt_v[3*(2*i+2)+1];
		tot_am[2] += prt_x[3*(2*i+1)+0]*prt_v[3*(2*i+1)+1]-prt_x[3*(2*i+1)+1]*prt_v[3*(2*i+0)+1];
	}
	for (unsigned k=0; k<3; ++k)
		tot_am[k] *= m_;
}

// calculate total linear momentum
void Domain::calc_tot_lm() {
	// set all values to 0
	for (unsigned k=0; k<3; ++k)
		tot_lm[k] = 0;
	//calculate...
	for (unsigned i=0; i<n_; ++i) {
		for (unsigned k=0; k<3; ++k)
			tot_lm[k] += prt_v[3*(2*i+1)+k];
	}
	for (unsigned k=0; k<3; ++k)
		tot_lm[k] *= m_;
}
