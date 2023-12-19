#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <random>
#include <bitset>
#include <thread>
//#include <bits/stdc++.h> // for std::sort

#include "planet.h"
#include "my_exceptions.h"
#include "utilities.h"
#include "orbital_mechanics.h"
#include "numeric"

#define SELECTION_ROULETTE  0
#define SELECTION_TOURNAMENT    1

#define TOURNAMENT_N    10

#define CROSS_UNIFORM   0
#define CROSS_SINGLE_GENE   1
#define CROSS_SINGLE_POINT  2
#define CROSS_DOUBLE_POINT  3

#define K_TIME_PENALTY 10

struct GenOperators{
    int n_population;       // Size of population
    int gen_limit;          // Limit of generations
    int elitism_n;          // Number of individuals in elit.
    int selectionType;      // Selection type
    int crossOverType;      // Crossover type
    float crossOverProb;    // Probability of doing crossover instead of reproduction. (set high like 0.9)
    float mutationProb;     // Probability of mutation
};


class ProblemDefinition {
public:
    float departure;                // Departure start window. This is the "base" date, all the times are added to this value.
    std::vector<Planet> planets;    // Sequence of planets.
    std::vector<std::pair<float, float>> timeWindows; // Time windows for each date {Td(min, max), T1(min,max), T2(min, max), ... }

    ProblemDefinition(const float _dep);    // Constructor with the departure as argument.
    ~ProblemDefinition();

    void add_planet(int _p, float min, float max);  // Adds a planet to the list with the time window as well.    
};


class Individual {
private:
    float getPenalty(const Planet& planet, double dV, double delta, double peri, double vin) const;    // Returns the penalty the cost of a flyby.
    void setChromosome(std::string chromo); // From a bitstring (chromosome) as argument, it sets the flyTimes vector. BitStr -> to -> vector of dates (solution).
    void setGene(std::string gene, int at); // Sets a single value in flyTime, i.e. one of the genes. Given the gene as bitstr, it converts it to the floating JD date.

public:
    std::vector<float> flyTimes;        // Vector of dates (solution/input of the ind) = Chromosome (each variable (time date) in the list is a gene).
    ProblemDefinition* problem;         // Problem reference (planets reference to operate are in there).
    double fitness;                     // Fitness of the individual.
    double totalDV = 0.0;               // Total amount of dV. Of that solution.

    Individual();
    Individual(ProblemDefinition* prob);    // Construct out of the problem reference.
    ~Individual();
    std::string getChromosome() const;      // Returns the individual chromosome.
    std::string getGene(int at) const;      // Returns a single gene (the one at "at" index in the time list).
    
    void mate(const Individual& partner, int crossType); // Mate with another individual to create a new child (self parent transforms to child).
    void createMutation();  // Create a mutation on the individual (using the bit flip method).

    void init();            // Randomly initizializes the vector of flyTimes.
    int getFlyTime() const; // Returns the whole duration of the interplanetary trajectory.

    void evaluate();        // Evaluates/computes the trajectory.

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
    void sortPopulation();                  // Sort the population from fitest to less fit. (fitest = less cost).
    std::vector<float> fitnessEvolution;    // Vector where the fitness evolution of the population can be saved. 

    void evolveNewGenerationThreaded(int indx_start, int indx_end); // Evolve a portion of the population, individuals between the indexes [indx_start; indx_end] in the population list. 

public:
    std::vector<Individual> population;     // Population vector: Contains all the individuals.
    std::vector<Individual> newPopulation;  // Auxiliar vector for the new generation during the evolution.
    const GenOperators geParameters;        // GA parameters
    int generationCount;                    // Count of generations runned.

    Population(const GenOperators params, ProblemDefinition* problem);
    ~Population();

    void inception();   // let the civilization begin (initis all the individuals to random time dates).

    // Genetic operators (self explanatory).
    void elitism();
    void selection();
    void crossOver();
    void mutate();

    void evolveNewGeneration();     // Evolves the current population to the next generation.
    void runGeneration();           // Runs the whole GA until the max number of generations has been reached.
};

#endif //INDIVIDUAL_H