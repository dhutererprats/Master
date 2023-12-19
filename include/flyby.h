#ifndef FLYBY_H
#define FLYBY_H

#include <vector>
#include <string>
#include <iostream>

#include "transfer.h"
#include "planet.h"

class Flyby{
public:
    Planet* planet;
    Transfer* trans1;
    Transfer* trans2;
    double v_in_rel[3];
    double v_out_rel[3];
    double delta, peri, dV;
    
    Flyby();
    ~Flyby();

    void add_planet_transfer(Transfer* t1, Transfer* t2);
    void compute_flyby();
    void print() const;
};

#endif //FLYBY_H