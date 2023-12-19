#include <iostream>
#include "mga_problem.h"
#include "genetic.h"
#include <ctime>
#include <chrono>
#include <cmath>


// Convert Bayesian Julian Date to Gregorian calendar date
//std::array<int,3> bayesianJulianToDate(double bjd, int& year, int& month, int& day)
int dJulianToDate(double bjd)
{
    int jd = static_cast<int>(std::round(bjd)); // Convert to integer Julian date

    int a = jd + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b / 4);

    int d = (4 * c + 3) / 1461;
    int e = c - (1461 * d / 4);
    int m = (5 * e + 2) / 153;

    int day = e - (153 * m + 2) / 5 + 1;
    int month = m + 3 - 12 * (m / 10);
    int year = 100 * b + d - 4800 + (m / 10);

    //std::array<int,3> date;
    int date[3]={0,0,0};
    std::cout << "date inside bayJultoDate function "<< date[0] << std::endl;
    //std::cout << "- Fitness of the solution: " << population.population.at(0).fitness << std::endl;
    date[0]=day;
    std::cout << "day inside bayJultoDate function "<< date[0] << std::endl;
    date[1]=month;
    date[2]=year;

    return (date[0],date[1],date[2]);
}

/*
int dJulianToDate(double bjd)
{
    int jd = static_cast<int>(std::round(bjd)); // Convert to integer Julian date

    int a = jd + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b / 4);

    int d = (4 * c + 3) / 1461;
    int e = c - (1461 * d / 4);
    int m = (5 * e + 2) / 153;

    int day = e - (153 * m + 2) / 5 + 1;
    int month = m + 3 - 12 * (m / 10);
    int year = 100 * b + d - 4800 + (m / 10);

    //std::array<int,3> date;
    int date[3]={0,0,0};
    std::cout << "date inside bayJultoDate function "<< date[0] << std::endl;
    //std::cout << "- Fitness of the solution: " << population.population.at(0).fitness << std::endl;
    date[0]=day;
    date[1]=month;
    date[2]=year;

    return date[0],date[1],date[2];
}
int dJulianToDate(double bjd)
{
    int jd = static_cast<int>(std::round(bjd)); // Convert to integer Julian date

    int a = jd + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b / 4);

    int d = (4 * c + 3) / 1461;
    int e = c - (1461 * d / 4);
    int m = (5 * e + 2) / 153;

    int day = e - (153 * m + 2) / 5 + 1;
    int month = m + 3 - 12 * (m / 10);
    int year = 100 * b + d - 4800 + (m / 10);

    //std::array<int,3> date;
    int date[3]={0,0,0};
    std::cout << "date inside bayJultoDate function "<< date[0] << std::endl;
    //std::cout << "- Fitness of the solution: " << population.population.at(0).fitness << std::endl;
    date[0]=day;
    date[1]=month;
    date[2]=year;

    return date[0],date[1],date[2];
}
*/
//#include <iostream>
//using namespace std;

void julianToGregorian(int jd, int& year, int& month, int& day) {
    int i, j, k, l, n;

    l = jd + 68569;
    n = (4 * l) / 146097;
    l = l - (146097 * n + 3) / 4;
    i = (4000 * (l + 1)) / 1461001;
    l = l - (1461 * i) / 4 + 31;
    j = (80 * l) / 2447;
    day = l - (2447 * j) / 80;
    l = j / 11;
    month = j + 2 - (12 * l);
    year = 100 * (n - 49) + i + l;
}


void voyager(GenOperators genOp){
    /**
     * @brief Example for voyager I mission.
     * genOp: Parameters for the genetic algorithms.
     * The function will optimize the the voyager I mission, provide the dates, display the result and trajectory.
     */

    // Start of the departure windoo.
    float Tdep = 2443145.5f; // 01-01-1977

    // Events windows (range of valid days).
    std::pair<float, float> dep_window = {0.f, 1095.f}; // Departure date: Earth.
    std::pair<float, float> t1_window = {50.f, 2000.f}; // Flyby date: Jupiter.
    std::pair<float, float> t2_window = {50.f, 2000.f}; // Arrival date: Saturn. 

    // Define the problem: Planets in the sequence and windows for each.
    ProblemDefinition prob = ProblemDefinition(Tdep);
    prob.add_planet(EARTH, dep_window.first, dep_window.second);
    prob.add_planet(JUPITER, t1_window.first, t1_window.second);
    prob.add_planet(SATURN, t2_window.first, t2_window.second); 

    // Create the population, intantiate it, and run the algorithm.
    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    // After the genetic algorithm has finished.
    // We create a mga problem to better analyze the optimal trajectory found by the GA.
    // We use the inputs of the BEST individual at the end (the optimal one) and compute what is its trajectory.
    // Individual at begining of population is the fittest one. 
    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();  // Compute that trajectory
    mga.print();    // Print a description of it. Values of the flyby, etc.
    mga.plot();     // Plot the trajectory.


    // Print also the raw solution (values of t0, t1, t2) find by the algorithm
    std::cout << "- Raw fly times: " << std::endl;
    for(const auto& ft: population.population.at(0).flyTimes){
        std::cout << ft << std::endl;
    }
    std::cout << "- Fitness of the solution: " << population.population.at(0).fitness << std::endl;
    std::cout << "- Delta-v cost of the solution: " << population.population.at(0).totalDV << std::endl;
}


void voyagerII(GenOperators genOp){
    /**
     * @brief Example for voyager II mission.
     * genOp: Parameters for the genetic algorithms.
     * The function will optimize the the voyager II mission, provide the dates, display the result and trajectory.
     */

    // Look at voyager I comments for explantions.
    // For this mission, the search space is higher, therefore more individuals and generations should be used.

    float Tdep = 2443145.5f; //01-01-1977
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

    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();
    mga.print();
    mga.plot();


    std::cout << "- Raw fly times: " << std::endl;
    for(const auto& ft: population.population.at(0).flyTimes){
        std::cout << ft << std::endl;
    }
    std::cout << "- Fitness of the solution: " << population.population.at(0).fitness << std::endl;
    std::cout << "- Delta-v cost of the solution: " << population.population.at(0).totalDV << std::endl;
}

//VOYAGER 3 TEST...
void voyager3(GenOperators genOp){
    /**
     * @brief Example for voyager 3 mission that is a copy of voyager II but leaves on my birthday in 2000-feb-6th
     * genOp: Parameters for the genetic algorithms.
     * The function will optimize the the voyager II mission, provide the dates, display the result and trajectory.
     */

    // Look at voyager I comments for explantions.
    // For this mission, the search space is higher, therefore more individuals and generations should be used.

    //float Tdep = 2451580.5f; //06-02-2000  converted using https://www.aavso.org/cgi-bin/cal2jd.pl
    float Tdep = 2443145.5f; //01-01-1977
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

    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();
    mga.print();
    mga.plot();


    std::cout << "- Raw fly times: " << std::endl;
    int counter=0;
    for(const auto& ft: population.population.at(0).flyTimes){
        std::cout<<"Counter is "<<counter<<std::endl;
        if(counter == 0)
        {
            std::array<int,3> gregdate;
            gregdate[0], gregdate[1], gregdate[2] = dJulianToDate(Tdep+ft);
            /*gregdate[0]=06;
            gregdate[1]=02;
            gregdate[2]=2000;*/
            //date=bayesianJulianToDate(Tdep+ft);
            std::cout << "Departure day: " << gregdate[0] <<"/"<<gregdate[1]<<"/"<<gregdate[2]<< std::endl;
        }
        std::cout << ft << std::endl;
        counter+=1;
    }
    std::cout << "- Fitness of the solution: " << population.population.at(0).fitness << std::endl;
    std::cout << "- Delta-v cost of the solution: " << population.population.at(0).totalDV << std::endl;
}


void gallileo(GenOperators genOp){
    /**
     * @brief Example for Galileo mission.
     * genOp: Parameters for the genetic algorithms.
     * The function will optimize the the Galileo mission, provide the dates, display the result and trajectory.
     */

    // Look at voyager I comments for explantions.
    // For this mission, the search space is higher, therefore more individuals and generations should be used.

    float Tdep = 2447161.5f;
    std::pair<float, float> dep_window = {0.f, 750.f};
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

    Population population = Population(genOp, &prob);
    population.inception();
    population.runGeneration();

    MGAProblem mga = MGAProblem(population.population.at(0));
    mga.compute();
    mga.print();
    mga.plot();

    std::cout << "- Raw fly times: " << std::endl;
    for(const auto& ft: population.population.at(0).flyTimes){
        std::cout << ft << std::endl;
    }
    std::cout << "- Fitness of the solution: " << population.population.at(0).fitness << std::endl;
    std::cout << "- Delta-v cost of the solution: " << population.population.at(0).totalDV << std::endl;
}



int main(int argc, char *argv[]){

    // Random seed updated here to obtain random values between different runs in all the random functions.
    std::srand(std::time(nullptr));
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();
    
    

    // Define the GA parameters.
    GenOperators genOp;
    genOp.n_population = 50000;         // Size of population.
    genOp.gen_limit = 50;               // Number of generations.
    genOp.elitism_n = 10;                        // Number of individuals in the elit group
    genOp.selectionType = SELECTION_TOURNAMENT; // Method of selection
    genOp.crossOverType = CROSS_DOUBLE_POINT;   // Method of cross over
    genOp.crossOverProb = 0.9;                  // Probability of cross over instead of reproduction
    genOp.mutationProb = 0.2;                   // Probability of mutation.
    
    voyager3(genOp);

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout<<"Time measured: %.3f seconds: "<< elapsed.count() * 1e-9<< std::endl;
}
