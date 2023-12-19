#pragma once
#ifndef ORBITAL_MECHANICS_H
#define ORBITAL_MECHANICS_H

#define MAX_ITER 10000

#include <iostream>

#include "utilities.h"
#include "astro.h"
#include "my_exceptions.h"

namespace orbit{
    void ephemeris(const orbitalParameters& planet_prm, const double T, double* r, double* v);
    void lambert(const double *r1_in, const double *r2_in, double t, const double &mu, double *v1, double *v2);
    void patched_conic(const double* Vin, const double* Vout, const double* Vplanet, const double mu, double& dV, double& delta, double& peri);
}
#endif // ORBITAL_MECHANICS_H