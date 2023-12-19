#include "genetic.h"

ProblemDefinition::ProblemDefinition(const float _dep) : departure(_dep){
}

ProblemDefinition::~ProblemDefinition(){
}

void ProblemDefinition::add_planet(int _p, float min, float max){
    /**
     * @brief Adds a planet to the list, and its time windows to.
     * Planet is added to its list. The time dates are added to its list.
     * We will consider: this->planets.at(i) has time dates at the same index, thus at: this->timeWindows.at(i).
     */
    this->planets.push_back(Planet(_p));
    
    // Sanity check.
    if(min >= max){
        throw timeRange();
    }

    // Add the time windows 
    this->timeWindows.push_back(std::make_pair(min, max));
}


/*
====================================================================================================================
                                                I N D I V I D U A L                                               ==
====================================================================================================================
*/

Individual::Individual(ProblemDefinition* prob){
    this->problem = prob;
    this->fitness = 0.f;
    this->totalDV = 0.f;
}

Individual::~Individual(){

}

std::string Individual::getChromosome() const{
    /**
     * @brief Get full chromosome in a bitstring
     * Flytimes are converted to bit string and concatenates to create the full chromosome.
     */
    std::string chromo;
    for(const auto& t: this->flyTimes){
        chromo.append(time2bitStr(t));
    }
    return chromo;
}

std::string Individual::getGene(int at) const{
    /**
     * @brief Get a given gene in bitstring format;
     */
    // Sanity check.
    if(at >= this->flyTimes.size() || at < 0){
        throw "Gene position out of range";
    }
    return time2bitStr(this->flyTimes.at(at));
}

void Individual::setChromosome(std::string chromo){
    /**
     * @brief Given a full chromosome, converts it to flight times and update it.
     * Calls setGene to set each time at flyTimes.
     */
    for(unsigned int i = 0; i < this->flyTimes.size(); i++){
        this->setGene(chromo.substr(i * (MAX_BIT_SIZE + MAX_FRACTIONAL_BIT_SIZE),(MAX_BIT_SIZE + MAX_FRACTIONAL_BIT_SIZE)), i);
    }   
}

void Individual::setGene(std::string gene, int at){
    /**
     * @brief Given a gene and position convert to fly time and update.
     * Upates the value at the corresponding index at the flyTimes vector.
     */
    // Sanity check.
    if(at >= this->flyTimes.size() || at < 0){
        throw "Gene position out of range";
    }

    float new_time = bitStr2Time(gene); // Gene to float (JD date)

    // Because the resolution of the gene allows for dates that may be outside the time window of that gene, we must consider this.
    // If the gene as a float is outside the window, we take limit of it.
    // Above the window = we take the upper limit.
    // Below the window = we take the lower limit.
    new_time = (new_time >= this->problem->timeWindows.at(at).first) ? new_time: this->problem->timeWindows.at(at).first;
    new_time = (new_time <= this->problem->timeWindows.at(at).second) ? new_time : this->problem->timeWindows.at(at).second;
    
    // Update the new time in the vector at the corresponding index.
    this->flyTimes.at(at) = new_time;
}


void Individual::mate(const Individual& partner, int crossType){
    /**
     * @brief Mating function that creates a new individual. 
     * The new individual is actually the current one, value on him are updated. Parent -> Transforms to "child".
     */
    if(crossType == CROSS_UNIFORM){
        // Take the full chromosomes. Pick randmly the bit for the child from one of the parents.
        std::string p1_ch = this->getChromosome();
        std::string p2_ch = partner.getChromosome();
        std::string child_ch = uniformBitstrCross(p1_ch, p2_ch);

        this->setChromosome(child_ch);
    }
    else if(crossType == CROSS_SINGLE_GENE){
        // Select a random Gene. Do CROSS_UNIFORM only on that gene.
        const int r_gene = rand() % this->flyTimes.size();          // Random chosen gene.
        std::string p1_g = time2bitStr(this->flyTimes.at(r_gene));  // Get gene of parent 1 as bitstr
        std::string p2_g = time2bitStr(partner.flyTimes.at(r_gene));// Get gene of parent 2 as bitstr
        std::string child_g = uniformBitstrCross(p1_g, p2_g);       // Cross the genes using the uniform technique

        this->setGene(child_g, r_gene); // Set the new gene on parent 1 (current this).
    }
    else if(crossType == CROSS_SINGLE_POINT){
        // From the full chromosome. At a given random point, interchange the parent 1 bits by the partner's bits.
        std::string p1_chromo = this->getChromosome();      // Chromosome of parent 1 (current this)
        std::string p2_chromo = partner.getChromosome();    // Chromosome of parent 2.
        
        const int cut_at = std::rand() % p1_chromo.size();  // Random cut point on the chromosomes.
        std::string child;
        child.append(p1_chromo.substr(0, cut_at)); // Start to cut_at.      
        child.append(p2_chromo.substr(cut_at)); // from cut_at pos to end.

        this->setChromosome(child); // Set the current new chromosome
    }
    else if(crossType == CROSS_DOUBLE_POINT){
        // From the full chromosome. At a given random point 1, interchange the parent 1 bits by the partner's bits. Do the same a second time at point 2.
        std::string p1_chromo = this->getChromosome();      // Chromosome of parent 1 (current this)
        std::string p2_chromo = partner.getChromosome();    // Chromosome of parent 2.
        
        const int cut_at_1 = rand_rng(1, p1_chromo.size() - 2); // From above 1, to at least leave 1 for the second cut.
        const int cut_at_2 = rand_rng(cut_at_1 + 1,  p1_chromo.size() - 1); // From at least one bit after cut one to end leave at least one for this cut.
        std::string child;
        child.append(p1_chromo.substr(0, cut_at_1));                    // Append the first part
        child.append(p2_chromo.substr(cut_at_1, cut_at_2 - cut_at_1));  // Append the second part
        child.append(p1_chromo.substr(cut_at_2));                       // Append the last part

        this->setChromosome(child);     // Set new chromosome.
    }
    else{
        throw "Unkown crossover type!";
    }
    
}

void Individual::createMutation(){
    /**
     * @brief Generates mutation in the individual.
     * Single bit flip technique.
     */
    std::string chromo = this->getChromosome();     // get the chromosome
    int mutate_at = rand_rng(0, chromo.size() - 1); // Chose a random flip from the chromosome
    chromo.at(mutate_at) = (chromo.at(mutate_at) == '0')? '1' : '0';    // Flip that bit

    this->setChromosome(chromo);    // Set again the new chromosome
}


void Individual::init(){
    /**
     * @brief Used to initizilize the individual to a random vector of times (within the windows).
     * Essentialy it is like the birth of the individuals.
     */
    for(const auto& window: this->problem->timeWindows){
        float t = window.first + rand() % (int)(window.second - window.first - 1);  // Random value inside the window for that date.
        t += rand_d();  // Add decimals to inrease resolution.
        this->flyTimes.push_back(t);
    }
}


void Individual::evaluate(){
    /**
     * @brief Solves the mga problem. Updates the cost.
     * It is done in a direct way basis, no info about the trajectory is saved more than the dates, dv and cost/fitness.
     * It is a direct version of compute from MGAProblem class. 
     * This function computes the trajectory coded by the dates of the individual.
     */
    // Reset from old execution.
    this->fitness = 0.f;
    this->totalDV = 0.f;

    float total_penalty = 0.f;
    float alpha=1.0; //weight of total delta V
    float beta=1-alpha; // weight of arrival time 
    float sf= 0.5f; // scaling factor 
    //float normalizing_dV=10000; //this should be mission dependent and coded as one of the mission's normalizing constants in the main.cpp and not here (long term)
    float normalizing_dV=7500; //Galileo mission normalizing constant
    for(const auto& window: this->problem->timeWindows){
        sf += (window.second -  window.first)/2;  // Use the mean value of the time window
    }

    sf = 1.0/sf;
    sf = normalizing_dV*sf;    //now the scaling factor, sf, is as described in the iPad notes "Meeting with Sebastian"


    // Per planets.
    // We take three planets each time P1 - P2 - P3.
    // We compute the transfers P1 - P2 and P2 - P3
    // We compute the flyby at P2.
    // This is done each time: Ex we have P1-P2-P3-P4-P5. We will have (P1-P2-P3) then (P2-P3-P4) then (P3-P4-P5) -> Three flybys.
    for(unsigned int i = 0; i < this->problem->planets.size() - 2; i++){
        double r1[3], v1[3], r2[3], v2[3], r3[3], v3[3];    // Variables to save the ephemeris.
        // Remenber that times are accumulative. We need real time, thus we need to compute them.
        double T1 = this->problem->departure + std::accumulate(this->flyTimes.begin(), this->flyTimes.begin() + i + 1, 0.0); // Times are in format [Departure T, +ΔT1, +ΔT2,..]
        double T2 = this->problem->departure + std::accumulate(this->flyTimes.begin(), this->flyTimes.begin() + (i + 1) + 1, 0.0);
        double T3 = this->problem->departure + std::accumulate(this->flyTimes.begin(), this->flyTimes.begin() + (i + 2) + 1, 0.0);
        // Ephemeris of the planets
        orbit::ephemeris(this->problem->planets.at(i).prm, T1, r1, v1);
        orbit::ephemeris(this->problem->planets.at(i + 1).prm, T2, r2, v2);
        orbit::ephemeris(this->problem->planets.at(i + 2).prm, T3, r3, v3);

        double v_dep1[3], v_arr2[3], v_dep2[3], v_arr3[3];  // Variables to save the transfers.
        // Transfers (lambert) computation.
        orbit::lambert(r1, r2, (T2 - T1) * DAY2SEC, MU_SUN, v_dep1, v_arr2);
        orbit::lambert(r2, r3, (T3 - T2) * DAY2SEC, MU_SUN, v_dep2, v_arr3);

        double dV, delta, peri; // Variables to save the flyby
        // Flyby computation
        orbit::patched_conic(v_arr2, v_dep2, v2, this->problem->planets.at(i + 1).mu, dV, delta, peri);

        // If departure planet: we have to add the departure delta-v
        if(i == 0){
            // Departure is speed must be considered at infinity, accouting for the speed of earth, automaticlly given to the prove.
            double v_dep_at_infi[3];
            minus2(v_dep1, v1, v_dep_at_infi);
            this->totalDV += norm(v_dep_at_infi); // Departure dV.
        }

        // Cost of flyby update
        total_penalty += this->getPenalty(this->problem->planets.at(i + 1), dV, delta, peri, norm(v_arr2)); // Flyby dV cost
        this->totalDV += dV;                                                                                // Flyby dV.
    }


    int length_flyTimes=sizeof(this->flyTimes) / sizeof(int);
    int timesincedeparturewindow = 0;

     // for loop runs from 0 to size - 1
    for(int i = 0; i < length_flyTimes; i++)
    {
    timesincedeparturewindow = timesincedeparturewindow + this->flyTimes[i];
    }

    // Update: the fitness is the cost. Objective function is the total amount of dV    ###DHP: Iker's old function
    //this->fitness = this->totalDV + total_penalty;

    // Update: the fitness is the cost. Objective function is the total amount of dV    ###DHP add the arrival date as a term here of a linear equation
    this->fitness = alpha*(this->totalDV + total_penalty)+beta*sf*(timesincedeparturewindow);
    // the new function including time shall look something like:
    //this->fitness = alpha*(this->totalDV + total_penalty)+beta*Arrival_Date;  ### however in order to obtain the arrival date one needs to be at the arrival planet (the end)! 
    //###hence the cost function for the previous steps needs to be evaluated differently (as it is done above) compared to the final evaluation
    //###so trying to eliminate/safe some computational time by starting to eliminate individuals before reaching the final planet because they are already a certain % above the current fittest indivudal is not an option
    //###however this also may be a good thing in terms of keeping individuals that help stay away from local minimums --> check "Meeting with Sebastian" IPAD notes 

}


float Individual::getPenalty(const Planet& planet, double dV, double delta, double peri, double vin) const{
    /**
     * @brief Returns the penalty cost.
     * - If the flyby is not feasable (too small peri for instance) it increases the cost a lot (to avoid solution evolving).
     */
    float penalty = 0.f;

    // Low perigee radius HARD penalty. ###this might need to be adjusted once time is added to the cost
    if(peri < 1.1 * planet.rad){
        penalty += 1e10;
    }

    // Low speed flyby penalty. Using the formula from Jacob A. Englander, Bruce A. Conway, and Trevor Williams. 
    // Automated mission planning via evolutionary algorithms
    double rsoi = std::pow((planet.mass / MASS_SUN), 2/5) * planet.sun_dist;
    double E = ((vin * vin) / 2) - planet.mu/rsoi; // Penalize low velocities flybys (can lead to spacecraft planet capture).
    if(E < 0){
        penalty += 1/vin;
    }
    return penalty;
}


/*
====================================================================================================================
                                                P O P U L A T I O N                                               ==
====================================================================================================================
*/
 
Population::Population(const GenOperators params, ProblemDefinition* problem) : geParameters(params), 
                                                                                population(params.n_population, Individual(problem)), 
                                                                                newPopulation(params.n_population, Individual(problem)){

}

Population::~Population(){

}

void Population::inception(){
    /**
     * @brief Initializes each individual randomly.
     */
    for(auto& ind: this->population){
        ind.init();     // Init random the individual
        ind.evaluate(); // Evaluate is trajectory and fitness. So we can start the evolution.
    }
}


void Population::sortPopulation(){
    /**
     * @brief Sorts the population vector from fitest to less fit.
     */
    std::sort(this->population.begin(), this->population.end());
}

// ----- Genetic operators. -----
void Population::selection(){
    /**
     * @brief Selects the individuals that will undergo reporduction/crossover.
     */
    if(this->geParameters.selectionType == SELECTION_ROULETTE){
        float sum_adjusted_ft = 0.f;
        float roulette[this->geParameters.n_population];
        
        for(unsigned int i = 0; i < this->geParameters.n_population; i++){
            sum_adjusted_ft += 1/(1 + this->population.at(i).fitness);
        }
        // Roulette wheel
        for(unsigned int j = 0; j < this->geParameters.n_population; j++){
            roulette[j] = (1/(1 + this->population.at(j).fitness)) / (sum_adjusted_ft);
        }
        // Do the selection
        float mean = 0;
        for(unsigned int k = 0; k < this->geParameters.n_population - this->geParameters.elitism_n; k++){
            float r = rand_d();
            float curr = roulette[0];
            int idx = 0;
            while (r >= curr){
                idx++;
                curr += roulette[idx]; 
            }
            idx = (idx <= this->geParameters.n_population - 1)? idx: this->geParameters.n_population - 1; //Sanity-check (the >= for floating points is not excat enough! In the end they value the same.)
            this->newPopulation.push_back(this->population.at(idx));
        }        
    }
    else if(this->geParameters.selectionType == SELECTION_TOURNAMENT){
        // Tournament method
        float mean = 0.f;
        for(unsigned int i = 0; i < this->geParameters.n_population - this->geParameters.elitism_n; i++){

            int rnd_idx[TOURNAMENT_N];
            for(int r = 0; r < TOURNAMENT_N; r++){
                rnd_idx[r] = rand_rng(0, this->geParameters.n_population - 1);
            }
            
            float best_fit = 1e20;
            int fitest_idx = -1;
            for(const auto& idx: rnd_idx){
                if(this->population.at(idx).fitness < best_fit){
                    fitest_idx = idx;
                    best_fit = this->population.at(idx).fitness;
                }
            }
            if(fitest_idx == -1){ // sanity check
                throw "Something went wrong during the tournament selection. No fittest than -1.f indivduales where found!";
            }
            mean += fitest_idx;
            this->newPopulation.push_back(this->population.at(fitest_idx)); // selected!
        }
    }
    else{
        throw "Unknow selection method";
    }

}

void Population::crossOver(){
    /**
     * @brief Cross over two parents
     */

    // Sanity check
    if(this->newPopulation.size() != this->geParameters.n_population){
        throw "The population expected to be crossed is missing individuales. Size is different from N_POPULATION";
    }
    
    // From elitism size upwards, there are the selected individuals.
    for(unsigned int i = this->geParameters.elitism_n; i < this->geParameters.n_population - 1; i = i+2){
        if(rand_d() <= this->geParameters.crossOverProb){
            // CrossOver (Each cross creates a chill out from the parent - 2 parents -> 2 cross -> 2 childs).
            Individual temp_p1 = this->newPopulation.at(i); // Direct update of child in parent. Should temp parent 1 (first to be changed).

            this->newPopulation.at(i).mate(this->newPopulation.at(i+1), this->geParameters.crossOverType);
            this->newPopulation.at(i+1).mate(temp_p1, this->geParameters.crossOverType);
        }
        else {
            // "Reproduction": Do nothing, include in new generation as they are.
        }
    }
}

void Population::mutate(){
    /**
     * @brief Mutate the population.
     */
    for(unsigned int i = this->geParameters.elitism_n - 1; i < this->geParameters.n_population -  this->geParameters.elitism_n; i++){
        if(rand_d() <= this->geParameters.mutationProb){ // if probability of mutation
            this->newPopulation.at(i).createMutation();
        }
    }
}

void Population::elitism(){
    /**
     * @brief Select the #n most fit individuals and insert them into the next generation.
     */

    // Sanity check
    if(this->newPopulation.size() != 0){
        throw "The new population vector must be cleaned before elitism operator";
    }

    for(unsigned int i = 0; i < this->geParameters.elitism_n; i++){
        this->newPopulation.push_back(this->population.at(i));
    }
}

// ***** Evolution ******
void Population::evolveNewGeneration(){
    /**
     * @brief Evaluates each individual.
     * In a multhreaded maneur.
     */
    int int_size = this->population.size() / 10;
    std::thread t1(&Population::evolveNewGenerationThreaded, this,  int_size * 0, int_size * 1);
    std::thread t2(&Population::evolveNewGenerationThreaded, this, int_size * 1, int_size * 2);
    std::thread t3(&Population::evolveNewGenerationThreaded, this,  int_size * 2, int_size * 3);
    std::thread t4(&Population::evolveNewGenerationThreaded, this, int_size * 3, int_size * 4);
    std::thread t5(&Population::evolveNewGenerationThreaded, this,  int_size * 4, int_size * 1);
    std::thread t6(&Population::evolveNewGenerationThreaded, this, int_size * 5, int_size * 6);
    std::thread t7(&Population::evolveNewGenerationThreaded, this,  int_size * 6, int_size * 7);
    std::thread t8(&Population::evolveNewGenerationThreaded, this, int_size * 7, int_size * 8);
    std::thread t9(&Population::evolveNewGenerationThreaded, this,  int_size * 8, int_size * 9);
    std::thread t10(&Population::evolveNewGenerationThreaded, this, int_size * 9, int_size * 10);
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
    t9.join();
    t10.join();
    // The join() waits for the threads to end. Its a blocking call.
}

void Population::evolveNewGenerationThreaded(int indx_start, int indx_end){
    /**
     * @brief Base evolution for multithreaded evolution
     */
    for(unsigned int i = indx_start; i < indx_end; i++){
        this->population.at(i).evaluate();  // Evaluate them
    }
}

void Population::runGeneration(){
    /**
     * @brief Runs the full genetic algorithm
     */

    int g = 0;
    this->newPopulation.clear();

    while (g < this->geParameters.gen_limit){
        this->sortPopulation();

        // Record evolution of the fitness or whatever if you need it. 
        //this->fitnessEvolution.push_back(this->population.at(0).fitness);

        this->elitism();
        this->selection();
        this->crossOver();
        this->mutate();

        this->population = this->newPopulation;
        this->newPopulation.clear();
        
        //Evolve
        this->evolveNewGeneration(); 
        g++;
    }

}
