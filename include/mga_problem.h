#ifndef MGA_PROBLEM_H
#define MGA_PROBLEM_H

#include <vector>
#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include "nlohmann/json.hpp"

#include "utilities.h"
#include "visualize.h"
#include "planet.h"
#include "transfer.h"
#include "flyby.h"
#include "genetic.h"

class MGAProblem{
public:
    std::vector<Planet> planets;
    std::vector<Transfer> transfers;
    std::vector<Flyby> flybys;
    std::vector<float> times;
    
    MGAProblem();
    MGAProblem(const Individual ind);
    ~MGAProblem();

    void add_planet(const int planet, float at);
    void compute_ephemeris();
    void compute_transfers();
    void compute_flybys();

    void compute();
    double computeCost() const;
    bool isSolutionValid() const;

    void plot() const;
    void print() const;
};


#endif //MGA_PROBLEM_H