#include "orbital_mechanics.h"


void orbit::ephemeris(const orbitalParameters& planet_prm, const double T, double* r, double* v){ //Outputs
    /**
     * @brief Computes the ephemeris of a given planet. 
     * Returns position and speed.
     */

    double T_pastCenJ2000 = (T - 2451545.0)/36525;

    // Orbital elements
    double a = planet_prm.a0 + planet_prm.acy * T_pastCenJ2000;
    double e = planet_prm.e0 + planet_prm.ecy * T_pastCenJ2000;
    double i = planet_prm.I0 + planet_prm.Icy * T_pastCenJ2000;
    double L = planet_prm.L0 + planet_prm.Lcy * T_pastCenJ2000;
    double p = planet_prm.long_peri0 + planet_prm.long_nodecy * T_pastCenJ2000;
    double W = planet_prm.long_node0 + planet_prm.long_nodecy * T_pastCenJ2000;

    // Convert deg to radians
    i = deg2rad(i);
    L = deg2rad(L);
    p = deg2rad(p);
    W = deg2rad(W);

    double M = L - p;
    double w = p - W; 
    
    // Newton-Raphson in order to get the eccentricity
    double E = M;
    double dE = 0.0;
    unsigned int iter = 0;
    while (true){
        dE = (E - e * std::sin(E) - M)/(1 - e * std::cos(E));
        E -= dE;
        if(std::abs(dE) < 1e-6){
            break;
        }
        if(iter > 10000){
            std::cout << "Iteration limit hit on N-R eccenticity anomaly computation." << std::endl;
            break;
        }
        iter++;
    }

    // P, Q 2d coordinate system - Position
    double P = a * (std::cos(E) - e);
    double Q = a * std::sin(E) * std::sqrt(1 - std::pow(e, 2));

    // P, Q 2d coordinate system - Velocity
    double vP = - (a * std::sin(E) * planet_prm.Lcy) / (1 - e * std::cos(E));
    double vQ = (a * std::cos(E) * std::sqrt(1 - e*e) * planet_prm.Lcy) / (1 - e * std::cos(E));

    rotate_eph(w, W, i, P, Q, r); 
    rotate_eph(w, W, i, vP, vQ, v);
    
    // Save in the output 
    r[0] *= UA;
    r[1] *= UA;
    r[2] *= UA;

    v[0] *= CORRECTION; 
    v[1] *= CORRECTION; 
    v[2] *= CORRECTION; 
}


void orbit::lambert(const double *r1_in, const double *r2_in, double t, const double &mu, double *v1, double *v2){
    /**
     * @brief Lambert solver. Given the ephemeris of two planets, and time of flight T.
     * It computes the heliocentric transfer speeds. Essentialy transfer orbit in heliocentric reference frame.
     * Reference: 
     * - This algorithm froms part of the PaGMO: ESA's Open-source for massively parallel engineering optimisation.
     * @see F. Biscani, D.Izzo, C. Hong Yam: A GLOBAL OPTIMISATION TOOLBOX FOR MASSIVELY PARALLEL ENGINEERING OPTIMISATION
     */
    double	V,T,
	r2_mod = 0.0,    // R2 module
	dot_prod = 0.0, // dot product
	c,		        // non-dimensional chord
	s,		        // non dimensional semi-perimeter
	am,		        // minimum energy ellipse semi major axis
	lambda,	        //lambda parameter defined in Battin's Book
	x,x1,x2,y1,y2,x_new=0,y_new,err,alfa,beta,psi,eta,eta2,sigma1,vr1,vt1,vt2,vr2,R=0.0;
	unsigned int i;
	const double tolerance = 1e-11;
	double r1[3], r2[3], r2_vers[3];
	double ih_dum[3], ih[3], dum[3];
    double theta, a, p;

	// direction decition
	double lw_out[3];
	cross_prod(r1_in, r2_in, lw_out);
	bool lw = lw_out[2] < 0; // z pomponent of result -> left/rigth if angle > 180º

	// Increasing the tolerance does not bring any advantage as the
	// precision is usually greater anyway (due to the rectification of the tof
	// graph) except near particular cases such as parabolas in which cases a
	// lower precision allow for usual convergence.

	if (t <= 0)
	{
		throw negativeTime();
	}

	for (i = 0; i < 3; i++)
	{
		r1[i] = r1_in[i];
		r2[i] = r2_in[i];
		R += r1[i] * r1[i];
	}

	R = std::sqrt(R);
	V = std::sqrt(mu/R);
	T = R/V;

	// working with non-dimensional radii and time-of-flight
	t /= T;
	for (i = 0; i < 3; i++)  // r1 dimension is 3
	{
		r1[i] /= R;
		r2[i] /= R;
		r2_mod += r2[i] * r2[i];
	}

	// Evaluation of the relevant geometry parameters in non dimensional units
	r2_mod = std::sqrt(r2_mod);

	for (i = 0; i < 3; i++)
		dot_prod += (r1[i] * r2[i]);

	theta = std::acos(dot_prod/r2_mod);

	if (lw)
		theta= 2 * std::acos(-1.0) - theta;

	c = std::sqrt(1 + r2_mod*(r2_mod - 2.0 * std::cos(theta)));
	s = (1 + r2_mod + c)/2.0;
	am = s/2.0;
	lambda = std::sqrt(r2_mod) * std::cos(theta/2.0)/s;

	// We start finding the log(x+1) value of the solution conic:
	// NO MULTI REV --> (1 SOL)
	//	inn1=-.5233;    //first guess point
	//  inn2=.5233;     //second guess point
	x1=std::log(0.4767);
	x2=std::log(1.5233);
	y1=std::log(x2tof(-.5233, s, c, lw)) - std::log(t);
	y2=std::log(x2tof(.5233, s, c, lw)) - std::log(t);

	// Regula-falsi iterations
	err = 1;
	while ((err>tolerance) && (y1 != y2))
	{
		x_new = (x1*y2-y1*x2)/(y2-y1);
		y_new = std::log(x2tof(std::exp(x_new)-1, s, c, lw)) - std::log(t);
		x1 = x2;
		y1 = y2;
		x2 = x_new;
		y2 = y_new;
		err = std::fabs(x1 - x_new);
	}
	x = std::exp(x_new)-1;

	// The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
	// now need the conic. As for transfer angles near to pi the lagrange
	// coefficient technique goes singular (dg approaches a zero/zero that is
	// numerically bad) we here use a different technique for those cases. When
	// the transfer angle is exactly equal to pi, then the ih unit vector is not
	// determined. The remaining equations, though, are still valid.

	a = am/(1 - x * x);

	// psi evaluation
	if (x < 1)  // ellipse
	{
		beta = 2 * std::asin (std::sqrt( (s-c)/(2*a) ));
		if (lw) beta = -beta;
		alfa = 2*std::acos(x);
		psi = (alfa-beta)/2;
		eta2 = 2 * a * std::pow(std::sin(psi),2)/s;
		eta = std::sqrt(eta2);
	}
	else       // hyperbola
	{
		beta = 2*std::asinh(sqrt((c-s)/(2*a)));
		if (lw) beta = -beta;
		alfa = 2*std::acosh(x);
		psi = (alfa-beta)/2;
		eta2 = -2 * a * std::pow(std::sinh(psi),2)/s;
		eta = std::sqrt(eta2);
	}

	// parameter of the solution
	p = ( r2_mod / (am * eta2) ) * std::pow(std::sin(theta/2), 2);
	sigma1 = (1/(eta * std::sqrt(am)))* (2 * lambda * am - (lambda + x * eta));
	cross_prod(r1, r2, ih_dum);
	vers(ih_dum, ih);

	if (lw)
	{
		for (i = 0; i < 3;i++)
			ih[i] = -ih[i];
	}

	vr1 = sigma1;
	vt1 = sqrt(p);
	cross_prod(ih, r1, dum);

	for (i = 0;i < 3 ;i++)
		v1[i] = vr1 * r1[i] + vt1 * dum[i];

	vt2 = vt1 / r2_mod;
	vr2 = -vr1 + (vt1 - vt2)/tan(theta/2);
	vers(r2, r2_vers);
	cross_prod(ih, r2_vers, dum);
	for (i = 0;i < 3 ;i++)
		v2[i] = vr2 * r2[i] / r2_mod + vt2 * dum[i];

	for (i = 0;i < 3;i++)
	{
		v1[i] *= V;
		v2[i] *= V;
	}
	a *= R;
	p *= R;
}

void orbit::patched_conic(const double* Vin, const double* Vout, const double* Vplanet, const double mu, //Input
                            double& dV, double& delta, double& peri){ //Output
    /**
     * @brief Used to pacth the two Lamebrt solutions. 
     * Takes the incoming and outgoing heliocentric velocities, converts it to planetocentric velocities (relative hyperbolic)
     * Then solves the required dV, turning angle and perigee radius to join the two trajectories.
     */

    //Convert to planetocentric velocities. relative hyperbolic. 
    double VinRel[3];
    double VoutRel[3];

    minus2(Vin, Vplanet, VinRel);
    minus2(Vout, Vplanet, VoutRel);

    double ain = -mu/vec2(VinRel);
    double aout = -mu/vec2(VoutRel);
    double e;


    delta = std::acos(dot_prod(VinRel, VoutRel)/(norm(VinRel) * norm(VoutRel)));
    e = 1/std::sin((delta/2));
    peri = aout*(1-e); 
    dV = std::fabs(std::sqrt(vec2(VinRel) + ((2*mu)/(peri))) - std::sqrt(vec2(VoutRel) + ((2*mu)/(peri))));
}