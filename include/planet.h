#ifndef PLANET_H
#define PLANET_H

#include <iostream>
#include <cmath>
#include <string>

#include "astro.h"
#include "my_exceptions.h"
#include "orbital_mechanics.h"

class Planet{
public:
    std::string name;
    orbitalParameters prm;
    double mu;
    double rad;
    double mass;
    double sun_dist;
    double r_eph[3];
    double v_eph[3];

    Planet();
    Planet(int planet);
    ~Planet();
    
    void setParameters(const float param[6], const float param_cy[6], planet_props::Properties props);
    void compute_eph(float at);
    void Planet_Ephemerides_Analytical(const double &, const int &);
};

#endif //PLANET