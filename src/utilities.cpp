#include "utilities.h"


void minus2(const double* v1, const double* v2, double* out){
    for(unsigned int i = 0; i < 3; i++){
        out[i] = v1[i] - v2[i];
    }
}


double norm(const double* v1){
    double out = 0.0;
    for(unsigned int i = 0; i < 3; i++){
        out += v1[i] * v1[i];
    }
    return std::sqrt(out);
}


double distance(const double* v1, const double* v2){
    double out = 0.0;
    for(unsigned int i = 0; i < 3; i++){
        out += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return std::sqrt(out);
}


double dot_prod(const double* v1, const double* v2){
    double out = 0.0;
    for(unsigned int i = 0; i < 3; i++){
        out += v1[i] * v2[i]; 
    }
    return out;
}


double vec2(const double* v1){
    double out = 0.0;
    for(unsigned int i = 0; i < 3; i++){
        out += v1[i] * v1[i]; 
    }
    return out;
}


void cross_prod(const double* v1, const double* v2, double* out){
    out[0]= v1[1] * v2[2] - v1[2] * v2[1];
	out[1]= v1[2] * v2[0] - v1[0] * v2[2];
	out[2]= v1[0] * v2[1] - v1[1] * v2[0];
}


double tofabn(const double &sigma, const double &alfa, const double &beta){
    if (sigma > 0)
    {
      return (sigma * std::sqrt (sigma)* ((alfa - std::sin(alfa)) - (beta - std::sin(beta))));
    }
	else
    {
      return (-sigma * std::sqrt(-sigma)*((std::sinh(alfa) - alfa) - (std::sinh(beta) - beta)));
    }
}


void vers(const double* v1, double* v2){
    double v_mod = 0;
	unsigned int i;
    v2[2] = 3;
	for (i = 0; i < 3; i++)
    {
		v_mod += v1[i] * v1[i];
	}

	double sqrtv_mod = std::sqrt(v_mod);

	for (i = 0;i < 3;i++)
	{
		v2[i] = v1[i]/sqrtv_mod;
	}
}


double x2tof(const double &x,const double &s,const double &c, bool prograde){
    double am, a, alfa, beta;

	am = s/2;
	a = am/(1-x*x);
	if (x < 1)//ellpise
    {
		beta = 2 * std::asin(std::sqrt((s - c)/(2*a)));
		if (prograde) beta = -beta;
		alfa = 2 * std::acos(x);
    }
	else
    {
		alfa = 2 * std::acosh(x);
		beta = 2 * std::asinh(std::sqrt((s - c)/(-2 * a)));
		if (prograde) beta = -beta;
    }

	if (a > 0)
    {
      return (a * std::sqrt(a)* ( (alfa - std::sin(alfa)) - (beta - std::sin(beta)) ));
    }
	else
    {
      return (-a * std::sqrt(-a)*((std::sinh(alfa) - alfa) - (std::sinh(beta) - beta)) );
    }

}


void rotate_eph(double w, double W, double i, double P, double Q, double* vec){
    // rotate by argument of periapsis
    vec[0] = std::cos(w) * P - std::sin(w) * Q;
    vec[1] = std::sin(w) * P + std::cos(w) * Q;
    // rotate by inclination
    vec[2] = std::sin(i) * vec[1];
    vec[1] = std::cos(i) * vec[1];
    // rotate by longitude of ascending node
    double xtemp = vec[0];
    vec[0] = std::cos(W) * xtemp - std::sin(W) * vec[1];
    vec[1] = std::sin(W) * xtemp + std::cos(W) * vec[1];
}


double deg2rad(const double x){
    return x*PI/180.0;
}

double rad2deg(const double x){
    return x*180.0/PI;
}


float rand_d(){
    return (float)std::rand()/RAND_MAX;
}

int rand_rng(int min, int max){
    return min + (std::rand() % static_cast<int>(max - min + 1));
}

std::string time2bitStr(const float x) {
    int dec = int(x);
    int fract = int((x - dec) * std::pow(10, MAX_FRACTIONAL_BIT_SIZE)); // 0.46571 -> 465 if MAX_..._SIZE = 3.
    fract = std::round(fract / (std::pow(10, MAX_FRACTIONAL_BIT_SIZE) / std::pow(2,MAX_FRACTIONAL_BIT_SIZE))); // above 0.9375 -> overflow to 0. No that relevant, allowed.
    return std::bitset<MAX_BIT_SIZE>((int)x).to_string() + std::bitset<MAX_FRACTIONAL_BIT_SIZE>(fract).to_string();
}

float bitStr2Time(const std::string x) {
    float numb = std::stoi(x.substr(0, MAX_BIT_SIZE), nullptr, 2);
    float decimal = std::stoi(x.substr(MAX_BIT_SIZE, MAX_FRACTIONAL_BIT_SIZE), nullptr, 2)/std::pow(2,MAX_FRACTIONAL_BIT_SIZE); 
    return numb+decimal;
}

std::string uniformBitstrCross(const std::string s1, const std::string s2){
    if(s1.size() != s2.size()){
        throw "Invalide sizes of bit strings to uniform cross. Sizes are not equal!";
    }

    std::string res;
    for(unsigned int i = 0; i < s1.size(); i++){
        char bit = (rand_d() <= 0.5) ? s1.at(i) : s2.at(i);
        res.push_back(bit);
    }
    return res;
}


size_t callback_curl(void* ptr, size_t size, size_t nmemb, std::string* data) {
    data->append((char*)ptr, size * nmemb);
    return size * nmemb;
}



double Mean2Eccentric (const double &M, const double &e)
{

	static const double tolerance = 1e-13;
	int n_of_it = 0; // Number of iterations
	double Eccentric,Ecc_New;
	double err = 1.0;




    Eccentric = M + e* std::cos(M); // Initial guess
    while ( (err > tolerance) && (n_of_it < 100))
    {
        Ecc_New = Eccentric - (Eccentric - e* std::sin(Eccentric) -M )/(1.0 - e * std::cos(Eccentric));
        err = std::fabs(Eccentric - Ecc_New);
        Eccentric = Ecc_New;
        n_of_it++;
    }
	return Eccentric;
}



void Conversion (const double *E,
				 double *pos,
				 double *vel,
				 const double &mu)
{
	double a,e,i,omg,omp,theta;
	double b,n;
	double X_per[3],X_dotper[3];
	double R[3][3];

	a = E[0];
	e = E[1];
	i = E[2];
	omg = E[3];
	omp = E[4];
	theta = E[5];

	b = a * std::sqrt (1 - e*e);
	n = std::sqrt (mu/std::pow(a,3));

	const double sin_theta = std::sin(theta);
	const double cos_theta = std::cos(theta);

	X_per[0] = a * (cos_theta - e);
	X_per[1] = b * sin_theta;

	X_dotper[0] = -(a * n * sin_theta)/(1 - e * cos_theta);
	X_dotper[1] = (b * n * cos_theta)/(1 - e * cos_theta);

	const double cosomg = std::cos(omg);
	const double cosomp = std::cos(omp);
	const double sinomg = std::sin(omg);
	const double sinomp = std::sin(omp);
	const double cosi = std::cos(i);
	const double sini = std::sin(i);

	R[0][0] =  cosomg*cosomp-sinomg*sinomp*cosi;
	R[0][1] =  -cosomg*sinomp-sinomg*cosomp*cosi;

	R[1][0] =  sinomg*cosomp+cosomg*sinomp*cosi;
	R[1][1] =  -sinomg*sinomp+cosomg*cosomp*cosi;

	R[2][0] =  sinomp*sini;
	R[2][1] =  cosomp*sini;

	// evaluate position and velocity
	for (int i = 0;i < 3;i++)
	{
		pos[i] = 0;
		vel[i] = 0;
		for (int j = 0;j < 2;j++)
		{
				pos[i] += R[i][j] * X_per[j];
				vel[i] += R[i][j] * X_dotper[j];
		}
	}
	return;
}