#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <iostream>
#include <cmath>
#include <string>

#include "astro.h"
#include "planet.h"
#include "my_exceptions.h"
#include "orbital_mechanics.h"
#include "genetic.h"
#include "mga_problem.h"

class Sequence{
public:
    std::vector<int> planets;
    std::vector<std::pair<float, float>> windows;
    float Tdep;

    float fitness; // Best of the 3rd

    Sequence(/*soemthing to define the seq*/);  //TODO:   
    ~Sequence();

    void quickEstimate(const GenOperators genOp);
    void goodEstimate(const  GenOperators genOp);

    // Used to facilitate sorting of sequence.
    bool operator< (const Sequence &other) const {
        return fitness < other.fitness;
    }

    bool operator== (const Individual &other) const {
        for(unsigned int i =0; i < this->planets.size(); i++){
            if(this->planets.at(i) != this->planets.at(i)){
                return false;
            }
        }
        return true;
    }
};

#endif //SEQUENCE