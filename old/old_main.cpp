#include <iostream>
#include "mga_problem.h"
#include "genetic.h"
#include <ctime>
#include <chrono>


void voyager(GenOperators genOp){
    float Tdep = 2443145.f; // 01-01-1977
    std::pair<float, float> dep_window = {0.f, 1095.f}; // NOTE: First value (min) should be 0 as it is the departure!
    std::pair<float, float> t1_window = {50.f, 2000.f};
    std::pair<float, float> t2_window = {50.f, 2000.f};

    ProblemDefinition prob = ProblemDefinition(Tdep);
    prob.add_planet(EARTH, dep_window.first, dep_window.second);
    prob.add_planet(JUPITER, t1_window.first, t1_window.second);
    prob.add_planet(SATURN, t2_window.first, t2_window.second);
    //prob.set_max_time(4.5*365.25); 

    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();
    //mga.print();
    //mga.plot();

    double inc = Tdep;
    for(const auto& ft: population.population.at(0).flyTimes){
        inc += ft;
        std::cout << ft << ", ";
    }
    std::cout << "   DV:" <<population.population.at(0).totalDV << " - fitness: "<<population.population.at(0).fitness << std::endl;

    population.plotFitnessEvolution();
}


void voyagerII(GenOperators genOp){
    float Tdep = 2443145.f; //01-01-1977
    std::pair<float, float> dep_window = {0.f, 1095.f};
    std::pair<float, float> t1_window = {50.f, 2000.f};
    std::pair<float, float> t2_window = {50.f, 2000.f};
    std::pair<float, float> t3_window = {500.f, 2500.f};
    std::pair<float, float> t4_window = {250.f, 2500.f};

    ProblemDefinition prob = ProblemDefinition(Tdep);
    prob.add_planet(EARTH, dep_window.first, dep_window.second);
    prob.add_planet(JUPITER, t1_window.first, t1_window.second);
    prob.add_planet(SATURN, t2_window.first, t2_window.second);
    prob.add_planet(URANUS, t3_window.first, t3_window.second);
    prob.add_planet(NEPTUNE, t4_window.first, t4_window.second);
    //prob.set_max_time(4.5*365.25); 

    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();
    mga.print();
    mga.plot();

    double inc = Tdep;
    for(const auto& ft: population.population.at(0).flyTimes){
        inc += ft;
        std::cout << ft << ", ";
    }
    std::cout << "   DV:" <<population.population.at(0).totalDV << " - fitness: "<<population.population.at(0).fitness << std::endl;
    
    //population.plotFitnessEvolution();
}


void gallileo(GenOperators genOp){
    float Tdep = 2447161.5f;
    std::pair<float, float> dep_window = {0.f, 883.f};
    std::pair<float, float> t1_window = {50.f, 750.f};
    std::pair<float, float> t2_window = {50.f, 750.f};
    std::pair<float, float> t3_window = {250.f, 1500.f};
    std::pair<float, float> t4_window = {250.f, 1500.f};

    ProblemDefinition prob = ProblemDefinition(Tdep);
    prob.add_planet(EARTH, dep_window.first, dep_window.second);
    prob.add_planet(VENUS, t1_window.first, t1_window.second);
    prob.add_planet(EARTH, t2_window.first, t2_window.second);
    prob.add_planet(EARTH, t3_window.first, t3_window.second);
    prob.add_planet(JUPITER, t4_window.first, t4_window.second);
    //prob.set_max_time(4.5*365.25); 

    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();
    mga.print();
    mga.plot();

    double inc = Tdep;
    for(const auto& ft: population.population.at(0).flyTimes){
        inc += ft;
        std::cout << ft << ", ";
    }
    std::cout << "   DV:" <<population.population.at(0).totalDV << " - fitness: "<<population.population.at(0).fitness << std::endl;
    
}

void test(){
    // REAL GALLILEO
    // MGAProblem mga = MGAProblem();
    // mga.add_planet(EARTH, 2447817.5); // 654
    // mga.add_planet(VENUS, 2447932.5); // 115
    // mga.add_planet(EARTH, 2448233.5); // 301
    // mga.add_planet(EARTH, 2448962.5); // 731
    // mga.add_planet(JUPITER, 2450150.5); // 1186

    // mga.compute();
    // mga.print();
    // mga.plot();

    // PAPER GALLILEO
    MGAProblem mga = MGAProblem();
    mga.add_planet(EARTH, 2443145.f + 246.0); // 687
    mga.add_planet(JUPITER, 2443145.f + 246.0 +772.5); // 97  
    mga.add_planet(SATURN, 2443145.f+ 246.0 +772.5 + 920.625); // 273
    // mga.add_planet(EARTH, 2448944.5); // 726
    // mga.add_planet(JUPITER, 2450034.5); // 1090

    mga.compute();
    mga.print();
    mga.plot();

    // // REAL VOYAGER
    // MGAProblem mga = MGAProblem();
    // mga.add_planet(EARTH, );
    // mga.add_planet(VENUS, );
    // mga.add_planet(EARTH, );
    // mga.add_planet(EARTH, );
    // mga.add_planet(JUPITER, );

    // mga.compute();
    // mga.print();
    // mga.plot();

    // // PAPER VOYAGER
    // MGAProblem mga = MGAProblem();
    // mga.add_planet(EARTH, );
    // mga.add_planet(VENUS, );
    // mga.add_planet(EARTH, );
    // mga.add_planet(EARTH, );
    // mga.add_planet(JUPITER, );

    // mga.compute();
    // mga.print();
    // mga.plot();
}

int main(int argc, char *argv[]){
    std::srand(std::time(nullptr));
    
    GenOperators genOp;
    genOp.elitism_n = 10;
    genOp.selectionType = SELECTION_TOURNAMENT;
    genOp.crossOverType = CROSS_DOUBLE_POINT;
    genOp.crossOverProb = 0.9;
    genOp.mutationProb = 0.2;

    for(unsigned int i = 0; i < 1; i++){
        auto start = std::chrono::high_resolution_clock::now();
        voyagerII(genOp);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Elapsed time:  " << elapsed.count() << " s" << std::endl;
    }
}
