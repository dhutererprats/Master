#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <random>
#include <bitset>
#include <thread>
using namespace std;  //added feb 2nd 
//#include <bits/stdc++.h> // for std::sort //THIS NEEDS TO BE FIXED WITH IKERS HELP THE FILE I CREATED DID NOT WORK

#include "planet.h"
#include "my_exceptions.h"
#include "utilities.h"
#include "orbital_mechanics.h"
#include "numeric"

#define N_POPULATION    50000
#define GEN_LIMIT  500

#define SELECTION_ROULETTE  0
#define SELECTION_TOURNAMENT    1

#define TOURNAMENT_N    10

#define CROSS_UNIFORM   0
#define CROSS_SINGLE_GENE   1
#define CROSS_SINGLE_POINT  2
#define CROSS_DOUBLE_POINT  3
#define CROSS_PERSONALIZED  4

#define K_TIME_PENALTY 10

struct GenOperators{
    int elitism_n;
    int selectionType;
    int crossOverType;
    float crossOverProb;
    float mutationProb;  
};


class ProblemDefinition {
public:
    float departure;
    std::vector<Planet> planets;
    std::vector<std::pair<float, float>> timeWindows; //{Td(min, max), T1(min,max), T2(min, max), ... }
    float time_max = 7305; //20 years

    ProblemDefinition(const float _dep);
    ~ProblemDefinition();

    void add_planet(int _p, float min, float max);
    void set_max_time(float t);
};


class Individual {
private:
    void updateDepartureCost(double dV);
    void updateCost(const Planet& planet, double dV, double delta, double peri, double vin);
    void setChromosome(std::string chromo);
    void setGene(std::string gene, int at);

public:
    std::vector<float> flyTimes;          // Chromosome (each variable is a gene).
    ProblemDefinition* problem;         // Problem reference (planets reference to operate are in there).
    double fitness;                      // Fitness of the individual.
    double cost;                         // Total cost of the individual based on the cost function.
    double totalDV = 0.0;               // Total amount of dV.

    Individual();
    Individual(ProblemDefinition* prob);
    ~Individual();
    std::string getChromosome() const;
    std::string getGene(int at) const;
    
    void mate(const Individual& partner, int crossType); // Mate with another individual to create a new child (self parent transforms to child).
    void createMutation();

    void init();
    int getFlyTime() const;

    void evaluate();

    // Used to facilitate sorting of individuals.
    bool operator< (const Individual &other) const {
        return fitness < other.fitness;
    }

    bool operator== (const Individual &other) const {
        // All have the same flight times => They are the same.
        for(unsigned int i =0; i < this->flyTimes.size(); i++){
            if(this->flyTimes.at(i) != other.flyTimes.at(i)){
                return false;
            }
        }
        return true;
    }
};


class Population{
private:
    void sortPopulation();
    std::vector<float> fitnessEvolution;

    void evolveNewGenerationThreaded(int indx_start, int indx_end);

public:
    std::vector<Individual> population;
    std::vector<Individual> newPopulation;
    const GenOperators geParameters;
    int generationCount;

    Population(const GenOperators params, ProblemDefinition* problem);
    ~Population();

    void inception();   // let the civilization begin.
    void mateIndividuals(Individual parent1, Individual parent2);
    void plotFitnessEvolution();

    // Genetic operators.
    void selection();
    void crossOver();
    void mutate();
    void elitism();

    void evolveNewGeneration();
    void runGeneration();
};

#endif //INDIVIDUAL_H
