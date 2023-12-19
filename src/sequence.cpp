#include "sequence.h"


Sequence::Sequence(/*something to represent the seq, and windows Tdep*/){
    this->fitness = 0.f;
    /*
    Set the sequence of planets //TODO:
    this->setUpSequence(...)
    */
}

Sequence::~Sequence(){

}


void Sequence::quickEstimate(const GenOperators genOp){
    /**
     * @brief Call it when you run to do the final/proper solver of the algorithm.
     * genOp should have the according number of generations and ind (LESS)
     */

    // Create the problem.
    ProblemDefinition prob = ProblemDefinition(this->Tdep);
    // Add the planets to the list, with the planet number, and the windows for it.
    for(unsigned int i = 0; i < this->planets.size(); i++){
        prob.add_planet(this->planets.at(i), this->windows.at(i).first, this->windows.at(i).second);
    }

    // Create the populaiton
    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    // Save the fitest individual. The one first and at the end of the population.
    this->fitness = population.population.at(0).fitness;
}



void Sequence::goodEstimate(const GenOperators genOp){
    /**
     * @brief Call it when you run to do the final/proper solver of the algorithm.
     * genOp should have the according number of generations and ind (MORE)
     */

    // Create the problem.
    ProblemDefinition prob = ProblemDefinition(this->Tdep);
    // Add the planets to the list, with the planet number, and the windows for it.
    for(unsigned int i = 0; i < this->planets.size(); i++){
        prob.add_planet(this->planets.at(i), this->windows.at(i).first, this->windows.at(i).second);
    }

    // Create the populaiton
    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    // Save the fitest individual. The one first and at the end of the population.
    this->fitness = population.population.at(0).fitness;

    // Also print the problem and everything
    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();  // Compute that trajectory
    mga.print();    // Print a description of it. Values of the flyby, etc.
    mga.plot();     // Plot the trajectory.
}